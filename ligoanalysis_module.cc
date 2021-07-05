////////////////////////////////////////////////////////////////////////
/// \brief The ligoanalysis module looks at GW coincidences.
///
/// Given a time window, specified as an absolute time and a delta, it
/// writes out a livetime histogram and an ntuple of supernova event
/// candidates.
///
/// \author M. Strait
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "DAQDataFormats/RawEvent.h"
#include "DAQDataFormats/RawTriggerMask.h"
#include "DAQDataFormats/RawDataBlock.h"
#include "RawData/FlatDAQData.h"
#include "RawData/RawTrigger.h"

#include "Geometry/Geometry.h"
#include "GeometryObjects/CellGeo.h"
#include "StandardRecord/SREnums.h"
#include "RecoBase/Track.h"
#include "MCCheater/BackTracker.h"

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

// For getting the event count when the file is opened
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/IO/Root/RootFileBlock.h"

#include <string>
#include <vector>
#include <algorithm>
#include <signal.h>

#include "progress.cpp"
#include "func/timeutil.h"


static const int TDC_PER_US = 64;
static const int US_PER_MICROSLICE = 50; // I hope this is always true

// Rough effective light speed in fiber
//
// PhotonTransport uses n=1.59 -> 18.9cm/ns. But for rays at just the
// critical angle (see doc-2665), it is effectively 16.9cm/ns. So the
// mean is 17.9cm/ns. PhotonTransport adds a time from a histogram to
// that gotten from the index. Does that model the angular effect?
static const float lightspeed = 17.9; // cm/ns
static const float invlightspeed = 1/lightspeed; // ns/cm

// Copied from GeometryBase so I can get pigtails with less overhead
static const double pigtailtimes[32] = {
   34.5738/lightspeed,   38.4379/lightspeed,
   42.3020/lightspeed,   46.1660/lightspeed,
   50.0301/lightspeed,   53.8942/lightspeed,
   57.7582/lightspeed,   61.6223/lightspeed,
   64.7504/lightspeed,   68.6144/lightspeed,
   72.4785/lightspeed,   76.3426/lightspeed,
   80.2067/lightspeed,   84.0707/lightspeed,
   87.9348/lightspeed,   91.0790/lightspeed,
   95.3301/lightspeed,   99.1941/lightspeed,
  103.058 /lightspeed,  106.922 /lightspeed,
  110.786 /lightspeed,  114.650 /lightspeed,
  118.514 /lightspeed,  122.379 /lightspeed,
  125.507 /lightspeed,  129.371 /lightspeed,
  133.235 /lightspeed,  137.099 /lightspeed,
  140.963 /lightspeed,  144.827 /lightspeed,
  148.691 /lightspeed,  150.751 /lightspeed
};

// z extent of a plane (slightly wrong at block boundaries, but for
// my purposes, it doesn't matter).  In cm.
static const double plnz = 6.6479;

// w extent of a cell (on average). It's the extrusion width plus the
// glue width (tiny), divided by the number of cells in an extrusion. cm.
static const double cellw = (63.455+0.048)/16.;

// Set from the FCL parameters
static double gwevent_unix_double_time = 0;
static long long WindowSize = 1000;

// Unsliced hits above these ADC levels are treated like slices instead
// of being allowed to be parts of supernova-like candidate clusters.
static const int16_t fd_high_adc = 1500;
static const int16_t nd_high_adc = 2500;

// At least at the ND, want to allow non-adjacent cells and planes. This
// makes it easier to accept neutrons, and also helps with positrons.
//
// What's wrong with attaching hits too far away from each other? (1)
// Lots background making huge file sizes that will just have to be cut
// later (2) Spoiling signal hits by attaching them to background hits.
// At the FD, the second motivates requiring hits to be adjacent.
//
// Values overwritten by fcl parameter.
static int MaxPlaneDist = 10;
static int MaxCellDist = 16;

// Generous box for regular Michel rejection. The first is mostly useful
// for doing a simple cut to reject most of the background. The second
// is for tricky cases that end up being important when all the easy stuff
// has already been rejected.
const int trk_pln_buf = 5, trk_cel_buf = 9;
const int big_trk_pln_buf = 40, big_trk_cel_buf = 67;

// Size of wedge forward of a track in which we'll look for stealth
// Michels. The first number is the distance in "planes", where a
// "plane" is just cm times the plane extent, so not really planes.
// The second number is the half-angle of the wedge in degrees. It's
// really a triangle because it is done in 2D. The number of planes is
// deliberately quite large, because I keep finding cases of stealth
// Michels implausibly far away in real data (I've seen one 30 planes
// and some cells away).
const double trkproj_pln_buf = 60, trkproj_ang_buf = 15;
const double trkproj_cm_buf = trkproj_pln_buf * plnz;
const double trkproj_cell_buf = trkproj_pln_buf * cellw;
const double cos_trkproj_ang_buf = cos(trkproj_ang_buf*M_PI/180);

// The "shape" variables measure how close in time a hit or cluster is
// to a slice, taking into account the closest distance between any hit
// in the slice and the candidate supernova hit. This defines how far,
// in units of plane widths, counts as being close to a slice for the
// "shape" variables.
const float shape_pln_buf = 24.;

// Box to put around each slice for determining if we are somewhere near it.
const int far_slc_pln_buf  = 48, far_slc_cel_buf  = 96;

static const int nplanefd = 896;
static const int ncellfd  = 384;

static bool chmaskv[nplanefd][ncellfd];
static float noiseratesv[nplanefd][ncellfd];
static float tposoverc[nplanefd][ncellfd];

static int gDet = caf::kUNKNOWN;

// Hit information, distilled from CellHit, plus more info
struct mhit{
  unsigned int hitid; // index into the noise slice hit array
  int32_t tdc; // coarse timing, in tdc ticks (1/64 us)
  float tns; // fine timing, in nanoseconds
  float tposoverc; // transverse position divided by speed of light
  float pe; // photoelectrons
  int16_t adc;
  int16_t plane;
  bool isx;
  int cell;
  bool used; // has this hit been used in a cluster yet?
  int truepdg; // For MC, what the truth is
  float trueE; // for MC, what the true initial particle kinetic energy is
  int16_t trueplane;
  int16_t truecellx;
  int16_t truecelly;

  // The minimum time to the next/previous nearby track end
  double totrkend_s;
  double sincetrkend_s;

  // The minimum time to the next/previous not-so-nearby track end
  double tofartrkend_s;
  double sincefartrkend_s;

  // The minimum time to the next/previous projected-forward track wedge
  double totrkproj_s;
  double sincetrkproj_s;

  // Same, for a region defined as a distance to any hit
  double toshapeslc_s;
  double sinceshapeslc_s;

  // The minimum time to the next/previous slice that's
  // overlapping this hit in space, with a big buffer
  double tofarslc_s;
  double sincefarslc_s;

  // Same, but special for trackless slices
  double totlslc_s;
  double sincetlslc_s;

  // The minimum time to the next slice anywhere
  double toanyslice_s;
  double sinceanyslice_s;

  // Roughly how noisy the channel with this hit is, where 1.0 is the
  // median noisiness. Really it is the channel's hit rate relative to
  // the median channel's hit rate, which of course is higher for some
  // channels because they have more real hits. Also for most input
  // files, it is the mean of that number over several subruns, which is
  // a slightly strange object, but probably meaningful enough.
  float noiselevel;
};

typedef std::vector<mhit *> sncluster;

struct mslice{
  // Number of tracks in this slice.  I only plan to use this to check
  // if it is zero, in which case the slice is going to be treated with
  // extra suspicion, but may as well store the number.
  int ntrack;

  float mintns, maxtns;

  // Extents
  int16_t minplane, maxplane;
  int16_t mincellx, maxcellx, mincelly, maxcelly;

  // Extents with buffers, which might save time. This "far" isn't the Far
  // Detector, but the size of the box around the slice.
  int16_t bminplane_far, bmaxplane_far;
  int16_t bmincellx_far, bmaxcellx_far, bmincelly_far, bmaxcelly_far;

  // Time in integer seconds and nanoseconds since Jan 1, 1970
  uint32_t time_s;
  uint32_t time_ns;

  // Pairs of plane, cell. (Heavy, but not too heavy.)
  vector< std::pair<int16_t, int> > xhits, yhits;
};

// Minimal track information, distilled from rb::Track, for use
// in correlating later activity to track ends
struct mtrack{
  // The index of the slice that this track is in.
  int slice;

  // The plane the track ended. If this was an X plane, the X cell the
  // track ended, and the estimated nearest Y cell. If it was a Y plane,
  // vice versa.
  int16_t endplane, endcellx, endcelly;

  // The direction unit vectors in each view.
  float yview_dirz, yview_diry;
  float xview_dirz, xview_dirx;

  float tns;
};

class ligoanalysis : public art::EDAnalyzer {
  public:
  explicit ligoanalysis(fhicl::ParameterSet const& pset);
  virtual ~ligoanalysis() { }; // compiles, but does not run, without this
  void analyze(const art::Event& evt);

  void beginSubRun(const art::SubRun& subrun);

  // Used to get the number of events in the file
  void respondToOpenInputFile(art::FileBlock const &fb);

  /// \brief If true, only report livetime withot examining data.
  bool fBlind;

  /// \brief The user-supplied time of the gravitational wave burst, or
  /// whatever time you want to center the search on.
  ///
  /// Expressed in RFC-3339 format, always in UTC, always with a Z at
  /// the end (to emphasize that it is UTC).
  std::string fGWEventTime;

  /// \brief The user-supplied length of the search window in seconds.
  ///
  /// We search for half this length on either side of the given time.
  float fWindowSize;

  /// \brief Largest number of planes away a hit can be to be clustered.
  ///
  /// If set to 1, planes must be contiguous.
  int fMaxPlaneDist;

  /// \brief Largest number of cells away a hit can be to be clustered.
  ///
  /// If set to 1, cells must be contiguous.  Note that this applies across
  /// planes, and cells are staggered, so the definition is a little odd.
  int fMaxCellDist;
};

static TH1D * livetimehist = NULL;

// Given the livetime of the event in seconds and the TDC count, return
// true if this object should be accepted:
//
// If the livetime is under 0.005 seconds, then this is not part of a
// multi-sub-trigger trigger, so there is no risk of overlapping data
// (not generally true, but true for the files we are choosing to read
// and the trigger bits we're accepting!), so return true.
//
// If the livetime is 0.005 seconds (the maximum possible), then return
// true if the TDC is non-negative (at a time equal to or greater than
// the time requested for this sub-trigger) and less than 5000*64, so as
// to exclude the parts of the first and last microslice that we'll use
// in the adjacent subtriggers.
//
// This works in the usual case in which the trigger requested 5000us
// and got 5050us because the trigger time was in the middle of a
// microslice. It also works in the case that the trigger time was on
// the microslice boundary and we got 5000us readouts without overlaps.
// In this case, there are no negative TDC values.
static bool uniquedata_tdc(const double livetime, const int32_t tdc)
{
  if(livetime < 0.005) return true;
  return tdc >= 0 && tdc < 5000 * TDC_PER_US;
}

static int trigger(
  const art::Handle< std::vector<rawdata::RawTrigger> > & rawtrigger)
{
  if(rawtrigger->empty()) return -1;
  return (*rawtrigger)[0].fTriggerMask_TriggerType;
}

static bool goodtriggertype(const int trigger)
{
  // See DAQDataFormats/cxx/include/TriggerDefines.h, and note the
  // off-by-one error between those defines and what appears in files

  // Ignore SNEWS fast beat and LIGO fast beat, because they will
  // overlap with real SNEWS/LIGO triggers, and slow beat SNEWS/LIGO
  // and be in the same files.
  if(trigger+1 == daqdataformats::TRIG_ID_SNEWS_BEAT_FAST) return false;
  if(trigger+1 == daqdataformats::TRIG_ID_LIGO_BEAT_FAST) return false;
  return true;
}

static bool longtriggertype(const int trigger)
{
  if(trigger+1 == daqdataformats::TRIG_ID_SNEWS_BEAT_SLOW) return true;
  if(trigger+1 == daqdataformats::TRIG_ID_LIGO_TRIGGER) return true;
  return false;
}

static bool is_complete_event(
  const art::Handle< std::vector<rawdata::FlatDAQData> > & flatdaq)
{
  daqdataformats::RawEvent raw;
  if(flatdaq->empty()) return false;

  raw.readData((*flatdaq)[0].getRawBufferPointer());
  if(raw.getDataBlockNumber() == 0) return false;

  return !raw.getHeader()->isEventIncomplete();
}

// Get the length of the event in TDC ticks, and "delta_tdc", the time
// between the trigger time and the time the event starts. You can
// subtract this off of the time that the offline gives each hit to
// get the time since the beginning of the readout, and with the event
// length, the time until the end of the readout.
//
// delta_tdc is a signed 64 bit integer, even though it should always be
// a small positive number, just in case. Ditto for the event length.
//
// Returns whether this information was successfully extracted.
static bool delta_and_length(int64_t & event_length_tdc,
  int64_t & delta_tdc,
  const art::Handle< std::vector<rawdata::FlatDAQData> > & flatdaq,
  const art::Handle< std::vector<rawdata::RawTrigger> > & rawtrigger)
{
  daqdataformats::RawEvent raw;
  if(flatdaq->empty()) return false;

  raw.readData((*flatdaq)[0].getRawBufferPointer());
  if(raw.getDataBlockNumber() == 0) return false;

  raw.setFloatingDataBlock(0);
  daqdataformats::RawDataBlock& datablock = *raw.getFloatingDataBlock();

  uint64_t event_start_time = 0xffffffffffffffff;
  uint64_t event_end_time   = 0x0000000000000000;

  for(unsigned int di = 0; di < raw.getDataBlockNumber(); di++){
    raw.setFloatingDataBlock(di);
    datablock = (*raw.getFloatingDataBlock());

    if(datablock.getHeader()->getMarker() ==
         daqdataformats::datablockheader::SummaryBlock_Marker ||
       !datablock.getHeader()->checkMarker()) continue;

    for(unsigned int mi = 0; mi < datablock.getNumMicroBlocks(); mi++){
      datablock.setFloatingMicroBlock(mi);
      daqdataformats::RawMicroBlock * ub = datablock.getFloatingMicroBlock();

      // The time is always in the second and third words of the
      // microslice, which follows two words of microblock header, so
      // just get it. Justin says you can also get it from getTime(),
      // but this already works and I'm not touching it.
      const uint32_t t_marker_low  = ((uint32_t *)(ub->getBuffer()))[3];
      const uint32_t t_marker_high = ((uint32_t *)(ub->getBuffer()))[4];

      uint64_t time_marker = t_marker_low;
      time_marker |= (uint64_t)t_marker_high << 32;
      if(time_marker < event_start_time) event_start_time = time_marker;
      if(time_marker > event_end_time  ) event_end_time   = time_marker;
    }
  }

  delta_tdc = (int64_t)((*rawtrigger)[0].fTriggerTimingMarker_TimeStart
                        - event_start_time);

  // Assume that microblocks are always 50us. I hope that's true for all
  // relevant data.
  event_length_tdc = ((int64_t)(event_end_time - event_start_time))
                     + US_PER_MICROSLICE*TDC_PER_US;
  return true; // ok
}

/***************************** Tree stuff ****************************/

static TTree * sntree = NULL;

// Structure of the output supernova analysis tree
struct sninfo_t{
  // Fraction of the hits that are caused by Monte Carlo particles.
  // Typically 1 when the cluster is caused by a Monte Carlo event and
  // 0 when it is caused by real detector activity, presumed to be
  // unrelated to any supernova. Of course, it will always be zero if
  // you are analyzing pure real data without Monte Carlo. It can be
  // between 0 and 1 if a Monte Carlo hit was accidentally next to a
  // background hit and got clustered with it.
  //
  // I *think* occasionally a hit can be a combination of real data
  // and Monte Carlo which were very close in time, in which case it
  // will either be counted as one or the other. I *think* any amount
  // of MC will cause it to be counted as an MC hit. In any case, the
  // numerator of 'truefrac' is always an integer.
  //
  // Occasionally a hit can be caused by an MC particle, but without
  // enough information stored in the intermediate file for BackTracker
  // to know that. This will cause 'truefrac' to be lower than it should
  // be. I believe this to affect << 1% of hits. I studied the case of
  // hits caused *only* by a positron annihilation gamma, whose true
  // tracks are by default not saved. Brems might also be relevant.
  float truefrac;

  // The PDG ID of the particle that contributed the plurality of the
  // PE to this cluster.  This might not be the majority of the energy
  // because of differing attenuation, but is a reasonable proxy for
  // the majority of the energy.  Anyway, in most cases one particle
  // contributes all of the energy.
  //
  // The main cases of interest are:
  //
  //  11: Electron from elastic scattering (rare)
  // -11: Positron from inverse beta decay
  //  22: Gamma from neutron capture from inverse beta decay, or
  //      Gamma from neutral current excitation of oxygen (except these
  //      aren't simulated yet)
  //
  // Note that you never see a hit from a neutron, because the
  // simulation considers the gammas that result from the neutron
  // capture to be responsible for the energy deposition. This is
  // arbitrary. What actually ionizes the scintillator is mainly the
  // Compton electrons from gamma scattering (neutrino -> neutron ->
  // gamma -> electron). Currently there's no way to distinguish between
  // gammas resulting from different processes, but that would be nice.
  int truepdg;

  // For Monte Carlo, the true initial kinetic energy, in MeV, of the
  // particle making the plurality of the cluster. This isn't the energy
  // of the neutrino, but if it is a positron, it is close to the energy
  // of the neutrino.
  float trueE;

  // For positrons, the true position of the neutrino interaction. For
  // gammas from neutron capture, the true position the neutron was
  // captured. The interaction is in either an X cell or a Y cell, and
  // the other cell number is the one closest to the interaction in the
  // other view.
  int16_t trueplane;
  int16_t truecellx;
  int16_t truecelly;

  // Number of hits in the cluster.  Only hits above a certain ADC
  // are accepted, so there may be other nearby hits that don't get
  // included because they are presumed noise.  "A certain ADC" could
  // be zero, in which case all hits are accepted.
  int nhit;

  // The TDC (64ths of microseconds since the trigger time) of the
  // earliest hit in the cluster. This is mainly useful for finding the
  // cluster in an event display where the TDC is available, but not the
  // absolute time.
  int32_t tdc;

  // Number of TDC from the earliest hit in the cluster to the beginning
  // and end of the trigger. This is useful primarily for the FD 10Hz
  // trigger where we don't have contiguous readout and need to know
  // where our ignorance starts.
  int32_t tdc_tobeginning;
  int32_t tdc_toend;

  // Mean time of the hits in this cluster in seconds and nanoseconds
  // since Jan 1, 1970.
  uint32_t time_s;
  uint32_t time_ns;

  // Total number of photoelectrons in this cluster in x cells and y
  // cells, respectively. The sum is a rough proxy for energy, but
  // doesn't take into account attenuation. When we have hits in both
  // views, we can calculate energy better than this.
  float pex;
  float pey;

  // Photoelectrons, corrected for attenuation, to get the number we
  // would have gotten if the fibers were perfectly clear. Only
  // filled if the cluster has both x and y hits. (At a later stage
  // of processing, can be filled for matched prompt/delayed pairs if
  // between them, they have both x and y hits.)
  float unattpe;

  // Lowest ADC of any hit in this cluster.  This is a measure of how
  // likely part of the cluster is to be APD noise.
  int16_t minhitadc;
  // And the highest ADC
  int16_t maxhitadc;

  // Number of nanoseconds from the first hit to the last hit. Large
  // times may indicate background. WARNING: We've seen that the Monte
  // Carlo is overly optimistic about how tight the detector timing is.
  // I have added an ad hoc correction which makes it look more like the
  // data, but it probably isn't super reliable.
  float timeext_ns;

  // Largest gap in time between two hits in the cluster.  Same as
  // timeext_ns if there are 2 hits.
  float maxtimegap_ns;

  // Smallest gap in time between two hits.  Same as timeext_ns
  // if there are 2 hits.  If there are more than 2 hits, may
  // tell you that 2 of the hits are related, and a third is
  // background.
  float mintimegap_ns;

  // First and last planes in the cluster
  int minplane;
  int maxplane;

  // Largest gap in planes. If there are fewer than three hits, this
  // is just maxplane-minplane-1, or zero if maxplane == minplane. If
  // there are three or more hits, this adds some information.
  int planegap;

  // Just maxplane - minplane + 1, for convenience.
  int planeextent;

  // First and last cells in the X view. Cell zero is on the east side
  // of the detector. Currently always defined, since we require a hit
  // in each view to form a cluster. If we relax this requirement, and
  // there is no X hit, mincellx and maxcellx will each be -1.
  int mincellx;
  int maxcellx;

  // Same idea as planegap, planeextent.  If there are no x hits, each
  // of these will be -1.
  int cellxgap, cellxextent;

  // First and last cells in the Y view. Cell zero is on the bottom. See
  // commentary on the X versions of these variables.
  int mincelly;
  int maxcelly;

  // See corresponding x variables.
  int cellygap, cellyextent;

  // The time in seconds until and since the last track end near this
  // cluster.
  double x_totrkend_s, y_totrkend_s;
  double x_sincetrkend_s, y_sincetrkend_s;

  // The time in seconds until and since the last track end not-so-near this
  // cluster.
  double x_tofartrkend_s, y_tofartrkend_s;
  double x_sincefartrkend_s, y_sincefartrkend_s;

  // The time in seconds until and since the last time a wedge projected
  // a long way forward from a track end included this cluster.
  double x_totrkproj_s, y_totrkproj_s;
  double x_sincetrkproj_s, y_sincetrkproj_s;

  double x_toshapeslc_s, y_toshapeslc_s;
  double x_sinceshapeslc_s, y_sinceshapeslc_s;

  // Same, but only for slices with no tracks in them, which probably
  // actually do contain muons that are nearly aligned with the planes.
  // This uses a smaller spatial buffer.
  double x_totlslc_s, y_totlslc_s;
  double x_sincetlslc_s, y_sincetlslc_s;

  // The time until and since the last slice that overlaps this cluster
  // in space, plus a big buffer of several cells and planes.  If this
  // hit is during such a slice, both are zero.  If there is no such
  // slice for one or the other case, that one is set to 1e9.
  double x_tofarslc_s, y_tofarslc_s;
  double x_sincefarslc_s, y_sincefarslc_s;

  // The time until and since the last slice anywhere. If this cluster is
  // during such a slice, both are zero. If there is no such slice for
  // one or the other case, that one is set to 1e9.
  double x_toanyslice_s, y_toanyslice_s;
  double x_sinceanyslice_s, y_sinceanyslice_s;

  // The NOvA run number
  unsigned int run;

  // The NOvA subrun number
  unsigned int subrun;

  // The event number, where in this case "event" means a 5ms data
  // segment (or whatever length of time was read out for each "event"
  // for this data).
  unsigned int event;

  // The hit number inside the event for one of the hits of this cluster.
  // It's the index of the hit in the noise slice, because CellHit::ID()
  // seems to always return zero.
  unsigned int hitid;

  // The maximum noiselevel of any hit in this cluster.  See definition
  // of noiselevel in struct mhit.
  float maxnoise;
} sninfo;

static void init_supernova()
{
  art::ServiceHandle<art::TFileService> t;
  sntree = t->make<TTree>("sn", "");

  // These reduce the nonsense of typos that goes with TTree:Branch()
  #define Q(x) #x
  #define BRN(x, t) sntree->Branch(Q(x), &sninfo.x, Q(x/t))
  BRN(nhit,             I);
  BRN(truefrac,         F);
  BRN(truepdg,          I);
  BRN(trueE,            F);
  BRN(trueplane,        S);
  BRN(truecellx,        S);
  BRN(truecelly,        S);
  BRN(tdc,              I);
  BRN(tdc_tobeginning,  I);
  BRN(tdc_toend,        I);
  BRN(time_s,           i);
  BRN(time_ns,          i);

  BRN(x_totrkend_s,       D);
  BRN(y_totrkend_s,       D);
  BRN(x_sincetrkend_s,    D);
  BRN(y_sincetrkend_s,    D);
  BRN(x_tofartrkend_s,    D);
  BRN(y_tofartrkend_s,    D);
  BRN(x_sincefartrkend_s, D);
  BRN(y_sincefartrkend_s, D);
  BRN(x_totrkproj_s,      D);
  BRN(y_totrkproj_s,      D);
  BRN(x_sincetrkproj_s,   D);
  BRN(y_sincetrkproj_s,   D);
  BRN(x_toshapeslc_s,     D);
  BRN(y_toshapeslc_s,     D);
  BRN(x_sinceshapeslc_s,  D);
  BRN(y_sinceshapeslc_s,  D);
  BRN(x_totlslc_s,        D);
  BRN(y_totlslc_s,        D);
  BRN(x_sincetlslc_s,     D);
  BRN(y_sincetlslc_s,     D);
  BRN(x_tofarslc_s,       D);
  BRN(y_tofarslc_s,       D);
  BRN(x_sincefarslc_s,    D);
  BRN(y_sincefarslc_s,    D);
  BRN(x_toanyslice_s,     D);
  BRN(y_toanyslice_s,     D);
  BRN(x_sinceanyslice_s,  D);
  BRN(y_sinceanyslice_s,  D);

  BRN(pex,              F);
  BRN(pey,              F);
  BRN(unattpe,          F);
  BRN(minhitadc,        S);
  BRN(maxhitadc,        S);
  BRN(timeext_ns,       F);
  BRN(maxtimegap_ns,    F);
  BRN(mintimegap_ns,    F);
  BRN(minplane,         I);
  BRN(maxplane,         I);
  BRN(planegap,         I);
  BRN(planeextent,      I);
  BRN(mincellx,         I);
  BRN(maxcellx,         I);
  BRN(cellxgap,         I);
  BRN(cellxextent,      I);
  BRN(mincelly,         I);
  BRN(maxcelly,         I);
  BRN(cellygap,         I);
  BRN(cellyextent,      I);
  BRN(run,              i);
  BRN(subrun,           i);
  BRN(event,            i);
  BRN(hitid,            i);
  BRN(maxnoise,         F);
}

/*************************** End tree stuff **************************/

// If we can't get the real number, this is roughly right for some FD files
static int64_t eventsinfile = 237;

void ligoanalysis::respondToOpenInputFile(const art::FileBlock &fb)
{
  // Get the number of events as soon as the file opens. Really
  // fragile. We get the number of entries in *some* tree,
  // which at the moment turns out to be the right one, but
  // EDAnalyzer::respondToOpenInputFile is totally undocumented as far
  // as I can see. But if it's wrong, it only affects the progress
  // indicator.
  //
  // If the job has more than one file, we don't know that until the
  // second one triggers this function. This will also just make the
  // progress indicator wrong.
  //
  // If the user gave -n and/or --nskip to limit the number of events,
  // we don't pick those up either. Chris Backhouse says "you can
  // introspect what the fcl configuration was for other modules in
  // the path (see CAFMaker for an example of trying to dig out genie
  // settings) so maybe you can get at the InputSource config". But I
  // don't think it's worth the effort.
  auto const* rfb = dynamic_cast<art::RootFileBlock const*>(&fb);

  if(rfb == NULL)
    printf("Can't get event count. Raw input? Assuming %lu.\n", eventsinfile);
  else
    eventsinfile = rfb->tree()->GetEntries();
}

void ligoanalysis::beginSubRun(const art::SubRun& subrun)
{
  if(fBlind) return;

  const int run = subrun.run(), sr = subrun.subRun();

  TFile * noisy = NULL;

  // Since when does TFile::TFile() throw an exception when it fails?
  // How obnoxious.  Catch its stupid exception, throw it away, and
  // handle errors the same way as always worked before.
  try{
    noisy = new TFile(
      Form("noisychannels_r%08d_s%02d.root", run, sr), "read");
  }catch(...){}

  if(noisy == NULL || noisy->IsZombie()){
    try{
      noisy = new TFile("noisychannels.root", "read");
    }catch(...){}
  }

  if(noisy == NULL || noisy->IsZombie()){
    fprintf(stderr, "Could not open noisychannels_r%08d_s%02d.root\n"
                    "or noisychannels.root. "
                    "Need this for low-energy analysis.\n"
                    "Try nova -c noisychannelsjob.fcl inputfile.root -T "
                    "noisychannels_r%08d_s%02d.root\n",
                    run, sr, run, sr);
    _exit(1);
  }

  TH2C * chmask = dynamic_cast<TH2C*>(noisy->Get("noisychannels/chmask"));
  if(chmask == NULL){
    fprintf(stderr, "Could not get chmask histogram from "
            "noisychannels_r%08d_s%02d.root", run, sr);
    _exit(1);
  }

  TH2D * noiserates = dynamic_cast<TH2D*>(noisy->Get("noisychannels/rates"));
  if(noiserates == NULL){
    fprintf(stderr, "Could not get rates histogram from "
            "noisychannels_r%08d_s%02d.root", run, sr);
    _exit(1);
  }

  art::ServiceHandle<geo::Geometry> geo;

  for(int i = 0; i < nplanefd; i++)
   for(int j = 0; j < ncellfd; j++){
     chmaskv[i][j] = chmask->GetBinContent(i, j);
     noiseratesv[i][j] = noiserates->GetBinContent(i, j);

     // Really lazy way of handling the ND.  Try all possible
     // channels, and ignore it when it doesn't work.
     try{
       tposoverc[i][j] = geo->CellTpos(i, j) * invlightspeed;
     }catch(...){}
   }
}

ligoanalysis::ligoanalysis(const fhicl::ParameterSet & pset) : EDAnalyzer(pset),
  fBlind(pset.get<bool>("Blind")),
  fGWEventTime(pset.get<std::string>("GWEventTime")),
  fWindowSize(pset.get<unsigned long long>("WindowSize")),
  fMaxPlaneDist(pset.get<int>("MaxPlaneDist")),
  fMaxCellDist(pset.get<int>("MaxCellDist"))
{
  // expose to static functions
  MaxPlaneDist = fMaxPlaneDist;
  MaxCellDist = fMaxCellDist;
  WindowSize = fWindowSize;
  gwevent_unix_double_time = rfc3339_to_unix_double(fGWEventTime);

  art::ServiceHandle<art::TFileService> t;

  livetimehist = t->make<TH1D>("blindlive", "",
    WindowSize, -WindowSize/2, WindowSize/2);

  if(!fBlind) init_supernova();
}

// Get the FlatDAQData, either from "minbias", in the case of supernova MC
// overlays, or "daq", for everything else.
static void getflatdaq(
  art::Handle< std::vector<rawdata::FlatDAQData> > & flatdaq,
  const art::Event & evt)
{
  evt.getByLabel("minbias", flatdaq);
  if(flatdaq.failedToGet())
    evt.getByLabel("daq", flatdaq);
}

static void getrawtrigger(
  art::Handle< std::vector<rawdata::RawTrigger> > & trg,
  const art::Event & evt)
{
  evt.getByLabel("minbias", trg);
  if(trg.failedToGet())
    evt.getByLabel("daq", trg);
}

static void getrawdigits(
  art::Handle< std::vector<rawdata::RawDigit> > & digits,
  const art::Event & evt)
{
  evt.getByLabel("minbias", digits);
  if(digits.failedToGet())
    evt.getByLabel("daq", digits);
}

// Return the livetime in this event in seconds, as it is relevant for
// raw hits (same as for anything else? Maybe not if we don't trust
// tracks close to the time-edges of events).
static double rawlivetime(const art::Event & evt, const bool veryraw = false)
{
  art::Handle< std::vector<rawdata::FlatDAQData> > flatdaq;
  getflatdaq(flatdaq, evt);
  if(flatdaq.failedToGet()){
    static bool first = true;
    if(first)
      puts("No FlatDAQ. Probably unoverlayed MC. Returning zero live time");
    first = false;
    return 0;
  }

  art::Handle< std::vector<rawdata::RawTrigger> > rawtrigger;
  getrawtrigger(rawtrigger, evt);
  if(flatdaq.failedToGet()){
    fprintf(stderr, "Unexpectedly failed to get flatdaq, returning -1\n");
    return -1;
  }

  if(rawtrigger->empty()) return -1;

  int64_t event_length_tdc = 0, delta_tdc = 0;

  if(!delta_and_length(event_length_tdc, delta_tdc, flatdaq, rawtrigger))
    return -1;

  const double wholeevent = event_length_tdc / TDC_PER_US * 1e-6;

  // Special case: If this is a 5ms subtrigger of a long trigger,
  // we are going to ignore the overlapping portions for the usual
  // case in which we get 0.00505 seconds, so report a livetime of
  // only the unignored portion.
  if(!veryraw && wholeevent > 0.005) return 0.005;

  return wholeevent;
}

// Return the bin number for this event, i.e. the number of seconds from the
// beginning of the window, plus 1.
//
// This does not handle an event overlapping a bin boundary particularly
// well, but does consistently assign it to the earlier bin.
static int timebin(const art::Event & evt, const bool verbose = false)
{
  const double evt_time = art_time_to_unix_double(evt.time().value());
  const double delta = evt_time - gwevent_unix_double_time;

#if 0
  if(verbose) printf("Accepted delta = %f seconds\n", delta);
#endif

  return floor(delta) + WindowSize/2
         + 1; // stupid ROOT 1-based numbering!
}

// Same thing as CellHit::operator== does, except this can be
// inlined, and I skip the check that the ADCs are the same, because
// how could there be two hits in the same cell at the same time?
// (I mean, except for overlaid MC)
static bool celleq(const art::Ptr<rb::CellHit> & a,
                   const art::Ptr<rb::CellHit> & b)
{
  return a->Cell() == b->Cell() &&
         a->Plane() == b->Plane() &&
         a->TDC() == b->TDC();
}

// return the index in the regular slice array that the given track is in.
static int which_slice_is_this_track_in(
  const rb::Track & t,
  const art::Handle< std::vector<rb::Cluster> > & slice)
{
  if(t.NCell() == 0) return -1;
  const art::Ptr<rb::CellHit> ahit =
    t.Cell(0); // some random hit on the track

  // Skip slice 0 since it is the noise slice
  const unsigned int nslice = slice->size();
  for(unsigned int i = 1; i < nslice; i++){
    const rb::Cluster & slc = (*slice)[i];

    if(ahit->View() == geo::kX){
      const unsigned int ncell = slc.NXCell();
      for(unsigned int j = 0; j < ncell; j++){
        const art::Ptr<rb::CellHit> & shit = slc.XCell(j);
        if(celleq(ahit, shit)) return i;
      }
    }
    if(ahit->View() == geo::kY){
      const unsigned int ncell = slc.NYCell();
      for(unsigned int j = 0; j < ncell; j++){
        const art::Ptr<rb::CellHit> & shit = slc.YCell(j);
        if(celleq(ahit, shit)) return i;
      }
    }
  }
  return -1;
}

struct hit { int16_t plane, cell; };

static bool comparecell(const hit & a, const hit & b)
{
  return a.cell < b.cell;
}

static bool compareplane(const hit & a, const hit & b)
{
  return a.plane < b.plane;
}

static bool comparemtracktime(const mtrack & a, const mtrack & b)
{
  return a.tns < b.tns;
}

// Returns a useful list of (window)track information, with all tracks
// going downward, and cell numbers for both x and y.
static std::vector<mtrack> make_trackinfo_list(const art::Event & evt,
  const art::Handle< std::vector<rb::Cluster> > & slice)
{
  art::ServiceHandle<geo::Geometry> geo;

  art::Handle< std::vector<rb::Track> > tracks;
  evt.getByLabel("windowtrack", tracks);
  if(tracks->empty() && gDet == caf::kFARDET)
    fprintf(stderr, "Unexpected empty windowtrack vector\n");

  // Don't assume input tracks are stored in time order. Process them
  // all, sort, and write out in time order.
  std::vector<mtrack> mtracks;

  for(const auto & trk : *tracks){
    if(!trk.Is3D()) continue;

    // I don't know of any documentation that explains what order
    // hits are stored, and I can't find a pattern that convinces
    // me I can trust the order for anything.  Sort them myself.

    std::vector<hit> xhits, yhits;

    for(unsigned int c = 0; c < trk.NYCell(); c++){
      hit h;
      h.cell  = trk.YCell(c)->Cell();
      h.plane = trk.YCell(c)->Plane();
      yhits.push_back(h);
    }
    for(unsigned int c = 0; c < trk.NXCell(); c++){
      hit h;
      h.cell  = trk.XCell(c)->Cell();
      h.plane = trk.XCell(c)->Plane();
      xhits.push_back(h);
    }

    // Efficiency note:
    // Don't have to do full sorts here, but it's not slow, so...

    // Sort Y cells into time order, assuming a regular downwards cosmic
    std::sort(yhits.begin(), yhits.end(), comparecell);
    std::reverse(yhits.begin(), yhits.end());

    // Ok, this is the bottom of the track in the Y view
    int16_t minycell = yhits[yhits.size()-1].cell;
    const int16_t lastyplane = yhits[yhits.size()-1].plane;

    // Figure out if the track end is in X or Y

    // Sort X cells into rough time order according to Y information
    std::sort(xhits.begin(), xhits.end(), compareplane);
    const bool increasingz = yhits[0].plane < yhits[yhits.size()-1].plane;
    if(!increasingz) std::reverse(xhits.begin(), xhits.end());

    const int16_t lastxplane = xhits[xhits.size()-1].plane;

    // In time order, is X cell number increasing?
    const bool increasingx = xhits[0].cell < xhits[xhits.size()-1].cell;

    int16_t endxcell = (*(increasingx?
      std::max_element(xhits.begin(), xhits.end(), comparecell):
      std::min_element(xhits.begin(), xhits.end(), comparecell))).cell;

    const int16_t endplane = increasingz? std::max(lastyplane, lastxplane)
                                         : std::min(lastyplane, lastxplane);

    const bool endisx = endplane == lastxplane;

    // Now, one of endxcell and minycell needs to be improved using
    // reconstructed track information.

    // Naively expect the end of the track to be lower than the
    // beginning.  Is it true?
    const bool forwardstrack = trk.Stop().Y() < trk.Start().Y();

    const double endx = forwardstrack? trk.Stop().X(): trk.Start().X();
    const double endy = forwardstrack? trk.Stop().Y(): trk.Start().Y();
    const double endz = forwardstrack? trk.Stop().Z(): trk.Start().Z();

    // Get a better cell position for the view that doesn't have the last
    // hit(s).  But only do this if the reco agrees with me that the muon
    // stopped in the last plane with hits
    int recoplane = 0, recocell = 0;
    try{
      geo->getPlaneAndCellID(endx, endy, endz-plnz, recoplane, recocell);

      recoplane++; // because I shifted by one to find the orthogonal cell

      if(recoplane == endplane){
        if(endisx) minycell = recocell;
        else       endxcell = recocell;
      }

    }
    catch(...){
      // Failed to get a better cell.  Just use the previous best estimate
      // in this rare case.  (Don't throw out the track, we need it.)
    }

    struct mtrack thisone;

    thisone.slice = which_slice_is_this_track_in(trk, slice);

    thisone.endcelly = minycell;
    thisone.endcellx = endxcell;
    thisone.endplane = endplane;

    if(forwardstrack){
      // Convert the unit vector in 3D to one in 2D since I want to work
      // in 2D for projecting tracks forward to look for stealth Michels.
      const double normxview = 1/hypot(trk.StopDir().X(), trk.StopDir().Z());
      if(isnan(normxview)){
        thisone.xview_dirz = 0;
        thisone.xview_dirx = 0;
      }
      else{
        thisone.xview_dirz = trk.StopDir().Z()*normxview;
        thisone.xview_dirx = trk.StopDir().X()*normxview;
      }

      const double normyview = 1/hypot(trk.StopDir().Y(), trk.StopDir().Z());
      if(isnan(normyview)){
        thisone.yview_dirz = 0;
        thisone.yview_diry = 0;
      }
      else{
        thisone.yview_dirz = trk.StopDir().Z()*normyview;
        thisone.yview_diry = trk.StopDir().Y()*normyview;
      }
    }
    else{
      const double normxview = 1/hypot(trk.Dir().X(), trk.Dir().Z());
      if(isnan(normxview)){
        thisone.xview_dirz = 0;
        thisone.xview_dirx = 0;
      }
      else{
        thisone.xview_dirz = -trk.Dir().Z()*normxview;
        thisone.xview_dirx = -trk.Dir().X()*normxview;
      }

      const double normyview = 1/hypot(trk.Dir().Y(), trk.Dir().Z());
      if(isnan(normyview)){
        thisone.yview_dirz = 0;
        thisone.yview_diry = 0;
      }
      else{
        thisone.yview_dirz = -trk.Dir().Z()*normyview;
        thisone.yview_diry = -trk.Dir().Y()*normyview;
      }
    }

    thisone.tns = trk.MeanTNS();

    mtracks.push_back(thisone);
  }

  std::sort(mtracks.begin(), mtracks.end(), comparemtracktime);

  return mtracks;
}

// Builds list of distilled slice information
static std::vector<mslice> make_sliceinfo_list(const art::Event & evt,
  const art::Handle< std::vector<rb::Cluster> > & slice,
  const vector<mtrack> & trkinfo)
{
  std::vector<mslice> sliceinfo;

  // Start at 1 to skip "noise" slice. Slice indices will all be off by
  // one, but makes it easier not to accidentally use the noise slice.
  for(unsigned int i = 1; i < slice->size(); i++){
    mslice slc;
    memset(&slc, 0, sizeof(mslice));

    if((*slice)[i].NCell() == 0) continue;

    for(const auto & hit : (*slice)[i].XCells())
      slc.xhits.push_back(std::pair<int16_t, int>(hit->Plane(), hit->Cell()));
    for(const auto & hit : (*slice)[i].YCells())
      slc.yhits.push_back(std::pair<int16_t, int>(hit->Plane(), hit->Cell()));

    for(const mtrack & t : trkinfo) if(t.slice == (int)i) slc.ntrack++;

    slc.mintns = (*slice)[i].MinTNS();
    slc.maxtns = (*slice)[i].MaxTNS();

    slc.time_s = art_time_plus_some_ns(evt.time().value(),
                      (*slice)[i].MeanTNS()).first;
    slc.time_ns = art_time_plus_some_ns(evt.time().value(),
                      (*slice)[i].MeanTNS()).second;

    slc.minplane = (*slice)[i].MinPlane();
    slc.maxplane = (*slice)[i].MaxPlane();

    slc.bminplane_far = (*slice)[i].MinPlane() - far_slc_pln_buf;
    slc.bmaxplane_far = (*slice)[i].MaxPlane() + far_slc_pln_buf;

    // If the slice has no hits in a view, we will exclude only plane 0,
    // cell 0 in that view
    if((*slice)[i].NCell(geo::kX) > 0){
      slc.mincellx = (*slice)[i].MinCell(geo::kX);
      slc.maxcellx = (*slice)[i].MaxCell(geo::kX);
      slc.bmincellx_far = (*slice)[i].MinCell(geo::kX) - far_slc_cel_buf;
      slc.bmaxcellx_far = (*slice)[i].MaxCell(geo::kX) + far_slc_cel_buf;
    }
    if((*slice)[i].NCell(geo::kY) > 0){
      slc.mincelly = (*slice)[i].MinCell(geo::kY);
      slc.maxcelly = (*slice)[i].MaxCell(geo::kY);
      slc.bmincelly_far = (*slice)[i].MinCell(geo::kY) - far_slc_cel_buf;
      slc.bmaxcelly_far = (*slice)[i].MaxCell(geo::kY) + far_slc_cel_buf;
    }

    // Take all slices, even if they are outside the 5ms window for
    // long trigger sub-triggers, because this is a list of things to
    // *exclude*.
    sliceinfo.push_back(slc);
  }

  // Make pseudo-slices for all unsliced hits >= high_adc.  These, near
  // signal-like clusters, are hints that they are background.
  const int16_t high_adc = gDet == caf::kNEARDET? nd_high_adc: fd_high_adc;
  for(unsigned int i = 0; i < (*slice)[0].NCell(); i++){
    const art::Ptr<rb::CellHit> & hit = (*slice)[0].Cell(i);

    if(hit->ADC() < high_adc) continue;

    const int cell  = hit->Cell();
    const int plane = hit->Plane();

    // Don't treat noisy channels as slices.
    // We're operating entirely with bin numbers and not bin ranges,
    // so this is not an off-by-one error.
    if(chmaskv[plane][cell]) continue;

    mslice slc;
    memset(&slc, 0, sizeof(mslice));
    slc.mintns = slc.maxtns = hit->TNS();

    slc.time_s =art_time_plus_some_ns(evt.time().value(), slc.mintns).first;
    slc.time_ns=art_time_plus_some_ns(evt.time().value(), slc.mintns).second;

    slc.bminplane_far = plane - far_slc_pln_buf;
    slc.bmaxplane_far = plane + far_slc_pln_buf;

    // In the view the slice isn't in, we will only exclude plane 0, cell 0.
    if(hit->View() == geo::kX){
      slc.bmincellx_far = cell - far_slc_cel_buf;
      slc.bmaxcellx_far = cell + far_slc_cel_buf;
    }
    else{ // Y-view
      slc.bmincelly_far = cell - far_slc_cel_buf;
      slc.bmaxcelly_far = cell + far_slc_cel_buf;
    }

    sliceinfo.push_back(slc);
  }

  return sliceinfo;
}

static void trueposition(int16_t & trueplane, int16_t & truecellx,
                         int16_t & truecelly,
                         const TLorentzVector & pos)
{
  art::ServiceHandle<geo::Geometry> g;
  int pln = -1, cl = -1;

  const double x = pos.X(), y = pos.Y(), z = pos.Z();

  // Find the cell nearest to the true starting position of this particle
  // (or whatever position we were passed).  Since "cell" means "inside
  // scintillator", we have to hunt around nearby for it.  Look 2cm in
  // 0.2mm steps in all.  This still doesn't work perfectly, but it has
  // 99.9% success, and most failures are around the detector edges, where
  // maybe they aren't really failures (if the particle interaction was
  // really outside the detector and something scattered in).
  for(int i = 0; i < 2; i++){
    const double Z = z - i*plnz;
    try{ g->getPlaneAndCellID(g->CellId(x, y, Z, 1.0, 0.0, 0.0, 2), pln, cl);
    }catch(...){
    try{ g->getPlaneAndCellID(g->CellId(x, y, Z, 0.0, 1.0, 0.0, 2), pln, cl);
    }catch(...){
    try{ g->getPlaneAndCellID(g->CellId(x, y, Z, 0.0, 0.0, 1.0, 2), pln, cl);
    }catch(...){
    try{ g->getPlaneAndCellID(g->CellId(x, y, Z, 0.7, 0.7, 0.0, 2), pln, cl);
    }catch(...){
    try{ g->getPlaneAndCellID(g->CellId(x, y, Z, 0.7, 0.0, 0.7, 2), pln, cl);
    }catch(...){
    try{ g->getPlaneAndCellID(g->CellId(x, y, Z, 0.0, 0.7, 0.7, 2), pln, cl);
    }catch(...){
    try{ g->getPlaneAndCellID(g->CellId(x, y, Z,-0.7, 0.7, 0.0, 2), pln, cl);
    }catch(...){
    try{ g->getPlaneAndCellID(g->CellId(x, y, Z,-0.7, 0.0, 0.7, 2), pln, cl);
    }catch(...){
    try{ g->getPlaneAndCellID(g->CellId(x, y, Z, 0.0,-0.7, 0.7, 2), pln, cl);
    }catch(...){
    try{ g->getPlaneAndCellID(g->CellId(x, y, Z, 0.6, 0.6, 0.6, 2), pln, cl);
    }catch(...){
    try{ g->getPlaneAndCellID(g->CellId(x, y, Z,-0.6, 0.6, 0.6, 2), pln, cl);
    }catch(...){
    try{ g->getPlaneAndCellID(g->CellId(x, y, Z, 0.6,-0.6, 0.6, 2), pln, cl);
    }catch(...){
    try{ g->getPlaneAndCellID(g->CellId(x, y, Z,-0.6,-0.6, 0.6, 2), pln, cl);
    }catch(...){
      #if 0
      printf("No %d cell for %.1f %.1f %.1f\n", i, x, y, Z);
      #endif
      pln = cl = -1;
    } } } } } } } } } } } } }

    if(i == 0) trueplane = pln;

    if(pln % 2 == 0) truecelly = cl;
    else truecellx = cl;
  }

  // I don't know why, but I get truecelly = 32690 in a small fraction
  // of cases. Protect against this. No apparent trouble with truecellx,
  // but protect it too.
  if(truecellx > ncellfd) truecellx = -1;
  if(truecelly > ncellfd) truecelly = -1;
}

static void tosince(double & oldsince, double & oldto,
                    const double since, const double to)
{
  if(since <= 0 && to <= 0){ oldto = oldsince = 0; return; }
  if(since > 0 && since < oldsince) oldsince = since;
  if(to > 0 && to < oldto) oldto = to;
}

static bool inside(const int low, const int high, const int x)
{
  return x > low && x < high;
}

static bool hitinbox(const int lowplane, const int highplane,
                     const int lowcellx, const int highcellx,
                     const int lowcelly, const int highcelly,
                     const int plane, const int cell,
                     const int view)
{
  if(!inside(lowplane, highplane, plane)) return false;

  if(view == geo::kX) return inside(lowcellx, highcellx, cell);
  else                return inside(lowcelly, highcelly, cell);
}

static bool compare_plane(const mhit & a, const mhit & b)
{
  return a.plane < b.plane;
}

// Helper function for supernova(). Selects hits that are candidates to
// be put into supernova clusters.
static std::vector<mhit> hits_for_supernova(const rb::Cluster & slice,
                                            const double livetime)
{
  std::vector<mhit> mhits;

  const int16_t high_adc = gDet == caf::kNEARDET? nd_high_adc: fd_high_adc;

  for(unsigned int i = 0; i < slice.NCell(); i++){
    const art::Ptr<rb::CellHit> & hit = slice.Cell(i);

    if(!uniquedata_tdc(livetime, hit->TDC())) continue;

    if(hit->ADC() >= high_adc) continue;

    const int cell  = hit->Cell();
    const int plane = hit->Plane();

    // Don't use noisy channels. We're operating entirely with bin
    // numbers and not bin ranges, so this is not an off-by-one error.
    if(chmaskv[plane][cell]) continue;

    // Ok, it's a good hit

    // Set only the things that are needed to determine if it goes
    // into a cluster.
    mhit h;
    memset(&h, 0, sizeof(mhit));
    h.hitid = i; // Because CellHit::ID() seems to always be zero
    h.plane = plane;
    h.cell  = cell;
    h.isx   = hit->View() == geo::kX;

    // CellHit::TNS() does *not* have pigtail corrections, it is just
    // the raw time that the electronics registered the hit. Correct
    // the time so it is as though the electronics were right at the
    // end of the module for all cells.
    h.tns = hit->TNS() - pigtailtimes[h.cell%32];

    // Smear out MC timing as per my study in doc-45041.
    //
    // Use a tolerable approximation (smear between 17.1ns and 43.4ns,
    // depending on position) that takes into account that the
    // data looks worse compared to the MC for the dim side of FD
    // modules.
    //
    // TODO: I haven't studied the ND at all, but just assume it is
    // the same as the FD as a function of distance to readout.
    if(hit->IsMC()){
      art::ServiceHandle<cheat::BackTracker> bt;
      std::vector<sim::FLSHit> flshits = bt->HitToFLSHit(hit);

      double true_w = 0;
      if(!flshits.empty()){
        const double loc[3] = { flshits[0].GetXAverage(),
                                flshits[0].GetYAverage(),
                                flshits[0].GetZAverage() };
        double wor[3] = { 0 };
        art::ServiceHandle<geo::Geometry> geo;
        geo->Plane(flshits[0].GetPlaneID())
           ->Cell (flshits[0].GetCellID())
           ->LocalToWorld(loc, wor);
        true_w = h.isx? wor[1]: wor[0];
      }

      const double extrulen = gDet == caf::kFARDET?1549.4:399.28;

      const double true_d = true_w + extrulen/2;

      static TRandom3 randfortiming;
      h.tns += randfortiming.Gaus() * (17.1 + true_d*0.017);
    }

    h.adc   = hit->ADC();
    h.tposoverc = tposoverc[plane][cell];

    mhits.push_back(h);
  }
  std::sort(mhits.begin(), mhits.end(), compare_plane);
  return mhits;
}

static bool hitinshape(const mslice & slc, const int16_t plane,
                       const int cell, const int view)
{
  const auto & hits = view == geo::kX? slc.xhits: slc.yhits;
  for(const auto & slchit : hits){
    if(abs(plane - slchit.first) > shape_pln_buf) continue;
    if(fabs((cell -  slchit.second)*cellw/plnz) > shape_pln_buf) continue;
    const float dist = hypot(plane - slchit.first,
                             (cell -  slchit.second)*cellw/plnz);
    if(dist < shape_pln_buf) return true;
  }

  return false;
}

// Fill in the details of this hit, mostly (but not exclusively) its
// proximity to cosmics.  Return false if the hit's characteristics
// should make the cluster fail preselection.  In this case, we probably
// haven't filled in everything about the hit.
static bool fill_in_hit(mhit & h,
  const rb::Cluster & noiseslice, const std::vector<mslice> & sliceinfo,
  const std::vector<mtrack> & trackinfo)
{
  const art::Ptr<rb::CellHit> & hit = noiseslice.Cell(h.hitid);

  const unsigned int nslice = sliceinfo.size(), ntrack = trackinfo.size();

  const int view  = hit->View();

  h.noiselevel = noiseratesv[h.plane][h.cell];

  h.sincetrkend_s = 1e9;
  h.totrkend_s = 1e9;
  h.sincefartrkend_s = 1e9;
  h.tofartrkend_s = 1e9;
  h.sincetrkproj_s = 1e9;
  h.totrkproj_s = 1e9;
  h.sinceshapeslc_s = 1e9;
  h.toshapeslc_s = 1e9;
  h.sincetlslc_s = 1e9;
  h.totlslc_s = 1e9;
  h.sincefarslc_s = 1e9;
  h.tofarslc_s = 1e9;
  h.sinceanyslice_s = 1e9;
  h.toanyslice_s = 1e9;

  for(unsigned int j = 0; j < nslice; j++){
    const mslice & slc = sliceinfo[j];
    const double timesince_s = (h.tns - slc.maxtns)*1e-9;
    const double timeto_s    = (slc.mintns - h.tns)*1e-9;

    tosince(h.sinceanyslice_s, h.toanyslice_s, timesince_s, timeto_s);

    if(hitinshape(slc, h.plane, h.cell, view))
      tosince(h.sinceshapeslc_s, h.toshapeslc_s, timesince_s, timeto_s);

    // Preselections
    if(h.toshapeslc_s    <= 2e-6) return false;
    if(h.sinceshapeslc_s <= 6e-6) return false;

    // Now with "far" spatial restriction
    if(hitinbox(slc.bminplane_far, slc.bmaxplane_far,
                slc.bmincellx_far, slc.bmaxcellx_far,
                slc.bmincelly_far, slc.bmaxcelly_far,
                h.plane, h.cell, view))
      tosince(h.sincefarslc_s, h.tofarslc_s, timesince_s, timeto_s);

    // Now cast special suspicion on trackless slices.  Treat them as though
    // they are made entirely of track ends.
    if(slc.ntrack == 0 &&
       hitinbox(slc.minplane - trk_pln_buf, slc.maxplane + trk_pln_buf,
                slc.mincellx - trk_cel_buf, slc.maxcellx + trk_cel_buf,
                slc.mincelly - trk_cel_buf, slc.maxcelly + trk_cel_buf,
                h.plane, h.cell, view))
      tosince(h.sincetlslc_s, h.totlslc_s, timesince_s, timeto_s);
  }

  // Construct measures of how close we are to track ends
  for(unsigned int j = 0; j < ntrack; j++){
    const mtrack & trk = trackinfo[j];
    const double since_s = (h.tns - trk.tns)*1e-9;
    const double to_s    = -since_s;

    // Simple boxes around the track end
    if(hitinbox(trk.endplane - trk_pln_buf, trk.endplane + trk_pln_buf,
                trk.endcellx - trk_cel_buf, trk.endcellx + trk_cel_buf,
                trk.endcelly - trk_cel_buf, trk.endcelly + trk_cel_buf,
                h.plane, h.cell, view))
      tosince(h.sincetrkend_s, h.totrkend_s, since_s, to_s);

    // Preselections
    if(h.totrkend_s    <=  2e-6) return false;
    if(h.sincetrkend_s <= 20e-6) return false;

    if(hitinbox(trk.endplane - big_trk_pln_buf,
                trk.endplane + big_trk_pln_buf,
                trk.endcellx - big_trk_cel_buf,
                trk.endcellx + big_trk_cel_buf,
                trk.endcelly - big_trk_cel_buf,
                trk.endcelly + big_trk_cel_buf,
                h.plane, h.cell, view))
      tosince(h.sincefartrkend_s, h.tofartrkend_s, since_s, to_s);

    // forward cone (well, really just a wedge (not a biggs))
    const double trkend2hit_z_pln = h.plane - trk.endplane;
    if(fabs(trkend2hit_z_pln) < trkproj_pln_buf){
      const double trkend2hit_w_cell =
        h.cell - (view == geo::kX?trk.endcellx
                               :trk.endcelly);

      // Avoid expensive hypot() if we're definitely too far away.
      if(fabs(trkend2hit_w_cell) < trkproj_cell_buf){

        const double trkend2hit_z = trkend2hit_z_pln * plnz;
        const double trkend2hit_w = trkend2hit_w_cell * cellw;
        const double dist_cm = hypot(trkend2hit_w, trkend2hit_z);

        if(dist_cm < trkproj_cm_buf){
          const double dotprod = view == geo::kX?
             (trk.xview_dirz*trkend2hit_z +
              trk.xview_dirx*trkend2hit_w)
            :(trk.yview_dirz*trkend2hit_z +
              trk.yview_diry*trkend2hit_w);

          // save lots of time by comparing cosines and not running acos()
          const double cosdtheta = dotprod/dist_cm;

          if(cosdtheta > cos_trkproj_ang_buf)
            tosince(h.sincetrkproj_s, h.totrkproj_s, since_s, to_s);
        }
      }
    }
  }

  h.tdc   = hit->TDC();
  h.pe    = hit->PE();
  h.adc   = hit->ADC();

  if(hit->IsMC()){
    art::ServiceHandle<cheat::BackTracker> bt;
    const std::vector<cheat::TrackIDE> & trackIDEs=bt->HitToTrackIDE(*hit);

    // Very occasionally, this is empty, maybe when "uninteresting"
    // particles are pruned, but turned out to give the only hits?
    if(!trackIDEs.empty()){
      const cheat::TrackIDE & trackIDE = trackIDEs.at(0);
      const sim::Particle* particle = bt->TrackIDToParticle(trackIDE.trackID);
      h.truepdg = particle->PdgCode();

      // Warning: sim::Particle::T() is the *time*, not the kinetic
      // energy like it should be!
      h.trueE = particle->E() - particle->Mass();

      trueposition(h.trueplane, h.truecellx, h.truecelly,
                   particle->Position());
    }else{
      h.truepdg = 81; // reserved for MC internal use, says PDG
      h.trueE = 0;
      h.trueplane = -1;
      h.truecellx = -1;
      h.truecelly = -1;
    }
  }
  else{
    h.truepdg = 0;
    h.trueE = 0;
    h.trueplane = -1;
    h.truecellx = -1;
    h.truecelly = -1;
  }

  return true;
}

// Return true iff this hit does cluster with the existing cluster.
static bool does_cluster(const sncluster & clu, const mhit & h)
{
  // 1. Check that the hit is in time with the first hit of the cluster.
  // Maybe it should have to be in time with the mean time of the
  // cluster so far or something. I think checking against the first hit
  // is good enough.

  // Used to be 250, but then I discovered that I never wanted clusters
  // with a greater total extent than 100ns.  But wait, with a better
  // analysis, I can accept more.  150ns seems close to optimal for now
  const float timewindow = 150; // ns

  if(clu[0]->isx == h.isx){
    // These hits are in the same view, so the times can be compared directly
    if(fabs(clu[0]->tns - h.tns) > timewindow) return false;
  }
  else{
    // These hits are in different views, so for each, use the other's
    // transverse position to correct the time.
    const float
      time_1st_corr = clu[0]->tns +       h.tposoverc,
      time_new_corr =       h.tns + clu[0]->tposoverc;

    if(fabs(time_new_corr - time_1st_corr) > timewindow) return false;
  }

  // 2. Check that if this hit is in the same view as another hit
  // in the cluster, that they are nearby.  Otherwise, fail it.  This
  // introduces an order dependence, and in principle we could pick
  // a bad hit first, and then fail a subsequent good hit.  We could
  // take hits tentatively waiting for a better one to come along.  But
  // let's indefinitely shelve that idea for the day when the simple
  // approach doesn't seem sufficient, which it probably is.

  bool in_same_view_as_another_hit = false;
  bool close_to_another_hit_in_w = false;

  for(unsigned int i = 0; i < clu.size(); i++)
    if(h.isx == clu[i]->isx){
      in_same_view_as_another_hit = true;
      break;
    }

  // TODO: When hits are in the same view in different planes, the
  // notion of cell distance is complicated. In the idealized version of
  // the detectors, each plane in a view is staggered half a cell width
  // from the next. So probably we should count in half-cells. On the
  // other hand, in the real detector, the plane-to-plane alignment is
  // only good to about half a cell width, so it's not clear how much it
  // matters.

  for(unsigned int i = 0; i < clu.size(); i++)
    if(h.isx == clu[i]->isx &&
       abs(h.cell - clu[i]->cell) <= MaxCellDist){
      close_to_another_hit_in_w = true;
      break;
    }

  if(in_same_view_as_another_hit && !close_to_another_hit_in_w) return false;

  return true;
}

// For the given supernova cluster, return the energy of the particle that
// contributed the most PE-weighted hits.  Ignore hits with no truth.
static float plurality_of_E(const sncluster & c)
{
  if(c.empty()) return 0;

  std::map<float, float> m;
  for(const mhit * h : c) if(h->trueE != 0) m[h->trueE] += h->pe;

  float most = 0;
  float best = 0;

  for(const auto & x : m)
    if(x.second > most){
      most = x.second;
      best = x.first;
    }

  return best*1000; // convert to MeV
}

// For the given supernova cluster, return the true plane of the particle
// that contributed the most PE-weighted hits.
static int16_t plurality_of_trueplane(const sncluster & c)
{
  if(c.empty()) return -1;

  std::map<int16_t, float> m;
  for(const mhit * h : c)
    if(h->trueplane != -1)
      m[h->trueplane] += h->pe;

  float most = 0;
  int16_t best = -1;

  for(const auto & x : m)
    if(x.second > most){
      most = x.second;
      best = x.first;
    }

  return best;
}

// For the given supernova cluster, return the true x cell of the particle
// that contributed the most PE-weighted hits.
static int16_t plurality_of_truecellx(const sncluster & c)
{
  if(c.empty()) return -1;

  std::map<int16_t, float> m;
  for(const mhit * h : c)
    if(h->truecellx != -1)
      m[h->truecellx] += h->pe;

  float most = 0;
  int16_t best = -1;

  for(const auto & x : m)
    if(x.second > most){
      most = x.second;
      best = x.first;
    }

  return best;
}

// For the given supernova cluster, return the true x cell of the particle
// that contributed the most PE-weighted hits.
static int16_t plurality_of_truecelly(const sncluster & c)
{
  if(c.empty()) return -1;

  std::map<int16_t, float> m;
  for(const mhit * h : c)
    if(h->truecelly != -1)
      m[h->truecelly] += h->pe;

  float most = 0;
  int16_t best = -1;

  for(const auto & x : m)
    if(x.second > most){
      most = x.second;
      best = x.first;
    }

  return best;
}

// For the given supernova cluster, return the true PDG id that contributed
// the most PE-weighted hits.  Ignore hits with no truth information.
static int plurality_of_truth(const sncluster & c)
{
  if(c.empty()) return 0;

  std::map<int, float> m;
  for(const mhit * h : c) if(h->truepdg != 0) m[h->truepdg] += h->pe;

  float most = 0;
  int best = 0;

  for(const auto & x : m)
    if(x.second > most){
      most = x.second;
      best = x.first;
    }

  return best;
}

// Return fraction of hits in this supernova cluster with truth information
static float fractrue(const sncluster & c)
{
  unsigned int ntrue = 0;
  for(const mhit * h : c) if(h->truepdg != 0) ntrue++;
  return (float)ntrue/c.size();
}

static double mean_tns(const sncluster & c)
{
  double sum = 0;
  for(const mhit * h : c) sum += h->tns;
  return sum / c.size();
}

// Return <pe in x cells, pe in y cells>
static std::pair<float, float> sum_pe(const sncluster & c)
{
  std::pair<float, float> sum(0, 0);
  for(const mhit * h : c) (h->isx?sum.first:sum.second) += h->pe;
  return sum;
}

static int min_plane(const sncluster & c)
{
  int ans = 9999; // There's always a plane, so result will never be 9999
  for(const mhit * h : c) if(h->plane < ans) ans = h->plane;
  return ans;
}

static int max_plane(const sncluster & c)
{
  int ans = -1; // There's always a plane, so result will never be -1
  for(const mhit * h : c) if(h->plane > ans) ans = h->plane;
  return ans;
}

static int plane_gap(const sncluster & c)
{
  std::vector<int> planes;
  for(const auto h : c) planes.push_back(h->plane);

  if(planes.size() == 1) return 0;

  std::sort(planes.begin(), planes.end());

  int maxgap = 0;
  for(unsigned int i = 0; i < planes.size()-1; i++){
    const int gap = planes[i+1] - planes[i] - 1;
    if(gap > maxgap)  maxgap = gap;
  }

  return maxgap;
}

static int cell_gap(const sncluster & c, const bool x)
{
  std::vector<int> cells;
  for(const auto h : c) if((h->isx ^ !x)) cells.push_back(h->cell);

  if(cells.size() == 0) return -1;
  if(cells.size() == 1) return 0;

  std::sort(cells.begin(), cells.end());

  int maxgap = 0;
  for(unsigned int i = 0; i < cells.size()-1; i++){
    const int gap = cells[i+1] - cells[i] - 1;
    if(gap > maxgap)  maxgap = gap;
  }

  return maxgap;
}

static int min_cell(const sncluster & c, const bool x)
{
  int ans = 9999;
  for(const auto h : c) if((h->isx ^ !x) && h->cell < ans) ans = h->cell;

  // It's more convenient if the invalid value is always -1
  if(ans == 9999) ans = -1;
  return ans;
}

static int max_cell(const sncluster & c, const bool x)
{
  int ans = -1;
  for(const auto h : c) if((h->isx ^ !x) && h->cell > ans) ans = h->cell;
  return ans;
}

static double to_trkend(const sncluster & c, const bool x)
{
  // Previously used the mean instead of the smallest. That's dumb,
  // because one of the hits could be right on top of a trkend, but if
  // the other one is too far away in space to see that trkend, it will
  // dilute away that information.
  double least = 1e9;
  for(const auto h : c)
    if((h->isx ^ !x) && h->totrkend_s < least)
      least = h->totrkend_s;
  return least;
}

static double since_trkend(const sncluster & c, const bool x)
{
  double least = 1e9;
  for(const auto h : c)
    if((h->isx ^ !x) && h->sincetrkend_s < least)
      least = h->sincetrkend_s;
  return least;
}

static double to_fartrkend(const sncluster & c, const bool x)
{
  double least = 1e9;
  for(const auto h : c)
    if((h->isx ^ !x) && h->tofartrkend_s < least)
      least = h->tofartrkend_s;
  return least;
}

static double since_fartrkend(const sncluster & c, const bool x)
{
  double least = 1e9;
  for(const auto h : c)
    if((h->isx ^ !x) && h->sincefartrkend_s < least)
      least = h->sincefartrkend_s;
  return least;
}

static double to_trkproj(const sncluster & c, const bool x)
{
  double least = 1e9;
  for(const auto h : c)
    if((h->isx ^ !x) && h->totrkproj_s < least)
      least = h->totrkproj_s;
  return least;
}

static double since_trkproj(const sncluster & c, const bool x)
{
  double least = 1e9;
  for(const auto h : c)
    if((h->isx ^ !x) && h->sincetrkproj_s < least)
      least = h->sincetrkproj_s;
  return least;
}

static double to_shapeslc(const sncluster & c, const bool x)
{
  double least = 1e9;
  for(const auto h : c)
    if((h->isx ^ !x) && h->toshapeslc_s < least)
      least = h->toshapeslc_s;
  return least;
}

static double since_shapeslc(const sncluster & c, const bool x)
{
  double least = 1e9;
  for(const auto h : c)
    if((h->isx ^ !x) && h->sinceshapeslc_s < least)
      least = h->sinceshapeslc_s;
  return least;
}

static double to_tlslc(const sncluster & c, const bool x)
{
  double least = 1e9;
  for(const auto h : c)
    if((h->isx ^ !x) && h->totlslc_s < least)
      least = h->totlslc_s;
  return least;
}

static double since_tlslc(const sncluster & c, const bool x)
{
  double least = 1e9;
  for(const auto h : c)
    if((h->isx ^ !x) && h->sincetlslc_s < least)
      least = h->sincetlslc_s;
  return least;
}

static double to_farslc(const sncluster & c, const bool x)
{
  double least = 1e9;
  for(const auto h : c)
    if((h->isx ^ !x) && h->tofarslc_s < least)
      least = h->tofarslc_s;
  return least;
}

static double since_farslc(const sncluster & c, const bool x)
{
  double least = 1e9;
  for(const auto h : c)
    if((h->isx ^ !x) && h->sincefarslc_s < least)
      least = h->sincefarslc_s;
  return least;
}

static double to_anyslice(const sncluster & c, const bool x)
{
  double least = 1e9;
  for(const auto h : c)
    if((h->isx ^ !x) && h->toanyslice_s < least)
      least = h->toanyslice_s;
  return least;
}

static double since_anyslice(const sncluster & c, const bool x)
{
  double least = 1e9;
  for(const auto h : c)
    if((h->isx ^ !x) && h->sinceanyslice_s < least)
      least = h->sinceanyslice_s;
  return least;
}

static int32_t first_tdc(const sncluster & c)
{
  int32_t ans = INT_MAX;
  for(const auto h : c) if(h->tdc < ans) ans = h->tdc;
  return ans;
}

static int16_t min_hit_adc(const sncluster & c)
{
  int16_t ans = SHRT_MAX;
  for(const auto h : c) if(h->adc < ans) ans = h->adc;
  return ans;
}

static int16_t max_hit_adc(const sncluster & c)
{
  int16_t ans = 0;
  for(const auto h : c) if(h->adc > ans) ans = h->adc;
  return ans;
}

static float max_noise(const sncluster & c)
{
  float ans = 0;
  for(const auto h : c) if(h->noiselevel > ans) ans = h->noiselevel;
  return ans;
}

static float time_ext_ns(float &mingap, float & maxgap,
                         const sncluster & c)
{
  double mintime = FLT_MAX, maxtime = FLT_MIN;
  if(c.size() <= 1) return 0;
  std::vector<float> reltimes;
  for(const auto & h : c){
    // Calculate all times as relative to the first hit.  Absolute time and
    // overall sign don't matter because we're finding the extent.
    float deltat = 0;
    if(c[0]->plane%2 == h->plane%2)
      deltat = c[0]->tns - h->tns;
    else{
      const float
        time_1st_corr = c[0]->tns +    h->tposoverc,
        time_oth_corr =    h->tns + c[0]->tposoverc;
      deltat = time_1st_corr - time_oth_corr;
    }

    if(deltat < mintime) mintime = deltat;
    if(deltat > maxtime) maxtime = deltat;
    reltimes.push_back(deltat);
  }

  std::sort(reltimes.begin(), reltimes.end());
  maxgap = 0, mingap = FLT_MAX;
  for(unsigned int i = 0; i < reltimes.size()-1; i++){
    const float gap = reltimes[i+1] - reltimes[i];
    if(gap > maxgap) maxgap = gap;
    if(gap < mingap) mingap = gap;
  }

  return maxtime - mintime;
}

// Should probably cache answers, since there are only 2*384 of them,
// except this isn't a hotspot, so shouldn't.
static double atten(const double w, const bool isx)
{
  // My light level tune from early 2020.  Good enough.
  const double viewfactor = isx? 0.6374: 0.5546;
  const double B = 0.4493, a1 = 203.6, a2 = 755.4;

  const double celllength = gDet == caf::kFARDET?1650.:399.28;

  const double att = viewfactor*(
       B *exp(-(celllength-w)/a1) +
    (1-B)*exp(-(celllength-w)/a2)

    +  B *exp(-(2*celllength-(celllength-w))/a1) +
    (1-B)*exp(-(2*celllength-(celllength-w))/a2)
  );

  return 1/att;
}

// If this event has XY position, calibrate it, i.e. write something useful
// into 'unattpe'.  This uses a (fairly) simple attenuation curve to get
// something proportional to energy.
static void calibrate(sninfo_t & ev)
{
  if(ev.maxcellx == -1 || ev.maxcelly == -1){
    ev.unattpe = -1;
    return;
  }

  const int16_t cellx  = (ev.maxcellx + ev.mincellx)/2;
  const int16_t celly  = (ev.maxcelly + ev.mincelly)/2;

  // This is close enough for our purposes (deviations on order of
  // a few mm because outer module walls are larger, etc.)
  static const double extruwidth_cm = 63.455;
  static const double two2onegluewidth_cm = 0.48;
  static const unsigned int cellspermodule = 32;
  static const double cellwidth = (2*extruwidth_cm + two2onegluewidth_cm)/
                             cellspermodule;

  // These are distances from the readout end, not the center
  const double x = cellx * cellwidth;
  const double y = celly * cellwidth;

  ev.unattpe = ev.pex * atten(y, true)
             + ev.pey * atten(x, false);
}

static void savecluster(const art::Event & evt, const sncluster & c)
{
  memset(&sninfo, 0, sizeof(sninfo));

  if(c.empty()){
    fprintf(stderr, "savecluster: got empty cluster, skipping\n");
    return;
  }

  // Already got this before, so (maybe) don't need to check it
  art::Handle< std::vector<rawdata::RawTrigger> > rawtrigger;
  getrawtrigger(rawtrigger, evt);

  art::Handle< std::vector<rawdata::FlatDAQData> > flatdaq;
  getflatdaq(flatdaq, evt);
  if(flatdaq.failedToGet()){
    fprintf(stderr, "Couldn't get flatdaq in savecluster()");
    exit(1);
  }

  int64_t event_length_tdc, delta_tdc;
  delta_and_length(event_length_tdc, delta_tdc, flatdaq, rawtrigger);

  sninfo.pex = sum_pe(c).first;
  sninfo.pey = sum_pe(c).second;

  sninfo.mincellx = min_cell(c, true);
  sninfo.maxcellx = max_cell(c, true);
  sninfo.mincelly = min_cell(c, false);
  sninfo.maxcelly = max_cell(c, false);

  calibrate(sninfo); // set unattpe

  sninfo.minhitadc  = min_hit_adc(c);

  // Preselection
  {
    if(sninfo.unattpe < 0){
      const double sumpe = sninfo.pex + sninfo.pey;
      if(sumpe <= 50 || sumpe >= 250) return;
    }
    else{
      if(sninfo.unattpe <= 500 || sninfo.unattpe >= 2500) return;
    }

    if(sninfo.minhitadc <= 40) return;
  }

  sninfo.nhit = c.size();
  sninfo.truefrac = fractrue(c);
  sninfo.truepdg = plurality_of_truth(c);
  sninfo.trueE = plurality_of_E(c);
  sninfo.trueplane = plurality_of_trueplane(c);
  sninfo.truecellx = plurality_of_truecellx(c);
  sninfo.truecelly = plurality_of_truecelly(c);
  sninfo.tdc = first_tdc(c);
  sninfo.tdc_tobeginning = sninfo.tdc + delta_tdc;
  sninfo.tdc_toend       = event_length_tdc - (sninfo.tdc + delta_tdc);
  sninfo.time_s = art_time_plus_some_ns(evt.time().value(),
                    mean_tns(c)).first;
  sninfo.time_ns = art_time_plus_some_ns(evt.time().value(),
                    mean_tns(c)).second;
  sninfo.maxhitadc  = max_hit_adc(c);
  sninfo.timeext_ns = time_ext_ns(sninfo.mintimegap_ns,
                                  sninfo.maxtimegap_ns, c);
  sninfo.minplane = min_plane(c);
  sninfo.maxplane = max_plane(c);
  sninfo.planegap = plane_gap(c);
  sninfo.planeextent = sninfo.maxplane - sninfo.minplane + 1;
  sninfo.cellxgap = cell_gap(c, true);
  sninfo.cellxextent = sninfo.maxcellx == -1? -1:
                       sninfo.maxcellx - sninfo.mincellx + 1;
  sninfo.cellygap = cell_gap(c, false);
  sninfo.cellyextent = sninfo.maxcelly == -1? -1:
                       sninfo.maxcelly - sninfo.mincelly + 1;

  sninfo.x_sincetrkend_s = since_trkend(c, true);
  sninfo.y_sincetrkend_s = since_trkend(c, false);
  sninfo.   x_totrkend_s =    to_trkend(c, true);
  sninfo.   y_totrkend_s =    to_trkend(c, false);

  sninfo.x_sincefartrkend_s = since_fartrkend(c, true);
  sninfo.y_sincefartrkend_s = since_fartrkend(c, false);
  sninfo.   x_tofartrkend_s =    to_fartrkend(c, true);
  sninfo.   y_tofartrkend_s =    to_fartrkend(c, false);

  sninfo.x_sincetrkproj_s = since_trkproj(c, true);
  sninfo.y_sincetrkproj_s = since_trkproj(c, false);
  sninfo.   x_totrkproj_s =    to_trkproj(c, true);
  sninfo.   y_totrkproj_s =    to_trkproj(c, false);

  sninfo.x_sinceshapeslc_s = since_shapeslc(c, true);
  sninfo.y_sinceshapeslc_s = since_shapeslc(c, false);
  sninfo.   x_toshapeslc_s = to_shapeslc(c, true);
  sninfo.   y_toshapeslc_s = to_shapeslc(c, false);

  sninfo.x_sincetlslc_s = since_tlslc(c, true);
  sninfo.y_sincetlslc_s = since_tlslc(c, false);
  sninfo.   x_totlslc_s =    to_tlslc(c, true);
  sninfo.   y_totlslc_s =    to_tlslc(c, false);

  sninfo.x_sincefarslc_s = since_farslc(c, true);
  sninfo.y_sincefarslc_s = since_farslc(c, false);
  sninfo.   x_tofarslc_s =    to_farslc(c, true);
  sninfo.   y_tofarslc_s =    to_farslc(c, false);

  sninfo.x_sinceanyslice_s = since_anyslice(c, true);
  sninfo.y_sinceanyslice_s = since_anyslice(c, false);
  sninfo.   x_toanyslice_s =    to_anyslice(c, true);
  sninfo.   y_toanyslice_s =    to_anyslice(c, false);

  sninfo.run = evt.run();
  sninfo.subrun = evt.subRun();
  sninfo.event = evt.event();
  sninfo.hitid = c[0]->hitid;

  sninfo.maxnoise = max_noise(c);

  sntree->Fill();
}

static bool comparebytime(const sncluster & a, const sncluster & b)
{
  return mean_tns(a) < mean_tns(b);
}

static mhit hitwithplane(const unsigned int plane)
{
  mhit ans;
  ans.plane = plane;
  return ans;
}

// Search for supernova-like events.
static void supernova(const art::Event & evt,
                      const std::vector<mslice> & sliceinfo,
                      const std::vector<mtrack> & trackinfo)
{
  art::Handle< std::vector<rb::Cluster> > slice;
  evt.getByLabel("slicer", slice);

  if(slice->empty()){
    fprintf(stderr, "Unexpected empty slice vector, skipping\n");
    return;
  }

  const double livetime = rawlivetime(evt);

  // Find hits which we'll accept for possible membership in pairs.
  std::vector<mhit> mhits = hits_for_supernova((*slice)[0], livetime);

  std::vector<sncluster> snclusters;

  // Starting with every eligible hit, form a cluster (possibly of size one)
  const unsigned int nmhits = mhits.size(); // really seems to help speed
  for(unsigned int i = 0; i < nmhits; i++){
    if(mhits[i].used) continue;

    sncluster clu;
    mhits[i].used = true;
    clu.push_back(&mhits[i]);

    bool done = false;
    do{
      done = true;

      // mhits is sorted by plane number. Find the first other hit that
      // has a plane number close enough that we'd use it. (Must find
      // something, because this is satisfied by the first hit added to
      // the cluster.)
      const unsigned int startat =
        std::lower_bound(mhits.begin(), mhits.end(),
                         hitwithplane(min_plane(clu) - MaxPlaneDist),
                         compare_plane)
        - mhits.begin();

      const unsigned int endat =
        std::lower_bound(mhits.begin(), mhits.end(),
                         hitwithplane(max_plane(clu) + MaxPlaneDist + 1),
                         compare_plane)
        - mhits.begin();


      for(unsigned int j = startat; j < endat; j++){
        if(mhits[j].used) continue;

        if(!does_cluster(clu, mhits[j])) continue;

        mhits[j].used = true;
        clu.push_back(&mhits[j]);
        done = false;
      }
    }while(!done);

    if(clu.size() < 2 || clu.size() > 7) continue;

    // Some preselection here to save time calculating slice
    // distance variables.  This is a massive time savings.
    if(gDet == caf::kFARDET){
      if(clu.size() == 2){
        if(clu[0]->adc + clu[1]->adc <= 160) continue;
        if(std::min(clu[0]->adc, clu[1]->adc) <= 50) continue;
      }

      bool passfid = true;

      for(auto h : clu)
        if(h->plane <= 4 || h->plane >= 876)
          passfid = false;

      if(!passfid) continue;

      bool hasx = false, hasy = false;

      for(auto h : clu){
        if(h->isx) hasx = true;
        else       hasy = true;
      }

      // Fail events that are near x edge regardless of 2D or 3D
      if(hasx)
        for(auto h : clu)
          if( h->isx && (h->cell <=  8 || h->cell >= 383 - 8))
            passfid = false;

      // Fail events near the top or bottom for y-only 2D events.
      if(!hasx && hasy)
        for(auto h : clu)
          if(!h->isx && (h->cell <= 15 || h->cell >= 350))
            passfid = false;

      if(!passfid) continue;
    }

    // Now that we're going to keep it, do the hard work
    bool passpresel = true;
    for(unsigned int j = 0; j < clu.size(); j++){
      if(!fill_in_hit(*clu[j], (*slice)[0], sliceinfo, trackinfo)){
        passpresel = false;
        break;
      }
    }

    if(!passpresel) continue;

    snclusters.push_back(clu);
  }

  std::sort(snclusters.begin(), snclusters.end(), comparebytime);

  // While we throw most of these events out in preselection later,
  // writing them out is an insignificant minority of the time usage.
  for(const auto & c : snclusters) savecluster(evt, c);
}

static void count_livetime(const art::Event & evt)
{
  const double livetime = rawlivetime(evt);
#if 0
  const double veryrawlivetime = rawlivetime(evt, true);

  art::Handle< std::vector<rawdata::RawTrigger> > rawtrigger;
  getrawtrigger(rawtrigger, evt);

  printf("Trigger type %d: %g seconds, %g accepted\n",
         trigger(rawtrigger), veryrawlivetime, livetime);
#endif
  int bin = timebin(evt, true);

  // Ensure out-of-range falls into the overflow bins rather than being lost
  if(bin < 0) bin = 0;
  if(bin > livetimehist->GetNbinsX()+1) bin = livetimehist->GetNbinsX()+1;

  // Use SetBinContent instead of Fill(x, weight) to skip one bin lookup
  livetimehist->SetBinContent(bin, livetimehist->GetBinContent(bin)
                                   + livetime);
}

void ligoanalysis::analyze(const art::Event & evt)
{
  {
    static unsigned int n = 0;
    if(n == 0) printf("Processing first event\n");
    // Start at 1 because 1st event takes forever to initialize everything
    if(n == 1) initprogressindicator(eventsinfile-1, 3);
    if(n > 0) progressindicator(n - 1);
    n++;
  }

  // Must be called on every event to prevent SIGPIPE loops.
  signal(SIGPIPE, SIG_DFL);

  art::Handle< std::vector<rawdata::RawTrigger> > rawtrigger;
  getrawtrigger(rawtrigger, evt);
  if(!goodtriggertype(trigger(rawtrigger))) return;

  const bool is_long_trigger = longtriggertype(trigger(rawtrigger));

  {
    art::Handle< std::vector<rawdata::FlatDAQData> > flatdaq;
    getflatdaq(flatdaq, evt);
    if(!flatdaq.failedToGet() && !is_complete_event(flatdaq)){
      printf("WARNING: Incomplete event, but assuming can trust livetime\n");
    }

    if(is_long_trigger){
      int64_t event_length_tdc, delta_tdc;
      delta_and_length(event_length_tdc, delta_tdc, flatdaq, rawtrigger);

      // Check for empty 50us events in long readouts
      if(event_length_tdc/TDC_PER_US == 50){
        art::Handle< std::vector< rawdata::RawDigit > > rd;
        getrawdigits(rd, evt);
        if(rd->empty()){
          printf("Rejecting trigger: 50us event with zero hits\n");
          return;
        }
      }

      // Check for truly incomplete readouts (not to be confused with events
      // marked incomplete).  Could use these, probably, but would need to
      // modify code to handle overlaps to determine when the overlaps are in
      // these cases, and the complexity doesn't seem worth it.
      const int len = event_length_tdc/TDC_PER_US;
      if(len != 5050 && len != 5000){
        printf("Rejecting trigger of length %dus != 5050 or 5000\n", len);
        return;
      }
    }
  }

  count_livetime(evt);

  if(fBlind) return;

  art::Handle< std::vector<rb::Cluster> > slice;

  art::ServiceHandle<geo::Geometry> geo;
  gDet = geo->DetId();

  evt.getByLabel("slicer", slice);
  // Event probably filtered out in reco file
  if(slice.failedToGet()) return;

  if(slice->empty()){
    fprintf(stderr, "Unexpected event with zero slices!\n");
    return;
  }

  const std::vector<mtrack> trackinfo = make_trackinfo_list(evt, slice);
  const std::vector<mslice> sliceinfo = make_sliceinfo_list(evt, slice,
                                                            trackinfo);

  // For holding tracks and slices from the previous trigger so (in the
  // rare case that triggers are 5 and not 5.05ms and therefore have no
  // overlaps) we can flag Michels, etc at very beginning of the trigger.
  static std::vector<mtrack> prev_trackinfo;
  static std::vector<mslice> prev_sliceinfo;

  // Will pass in the concatenation of both sets.  Duplicates
  // are not a problem, because we always use the closest one.
  std::vector<mtrack> trackinfo_wprev = prev_trackinfo;
  trackinfo_wprev.insert(trackinfo_wprev.end(),
                         trackinfo.begin(), trackinfo.end());

  std::vector<mslice> sliceinfo_wprev = prev_sliceinfo;
  sliceinfo_wprev.insert(sliceinfo_wprev.end(),
                         sliceinfo.begin(), sliceinfo.end());

  supernova(evt, sliceinfo_wprev, trackinfo_wprev);

  // Translate times for the next event, 5ms later. This ONLY makes
  // sense for long readouts, and is wrong if a trigger is dropped. If
  // it's not a long trigger, don't fill these vectors. The effect will
  // be to admit more background because we don't know about cosmics
  // as long ago, but usually only the previous ~25us matters. And, as
  // always we measure the background with the data, so it comes out in
  // the wash.
  if(is_long_trigger){
    prev_trackinfo = trackinfo;
    prev_sliceinfo = sliceinfo;

    for(unsigned int i = 0; i < prev_trackinfo.size(); i++)
      prev_trackinfo[i].tns -= 5e6; // 5 ms in ns
    for(unsigned int i = 0; i < prev_sliceinfo.size(); i++){
      prev_sliceinfo[i].mintns -= 5e6;
      prev_sliceinfo[i].maxtns -= 5e6;
    }
  }
}

DEFINE_ART_MODULE(ligoanalysis)
