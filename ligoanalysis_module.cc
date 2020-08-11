////////////////////////////////////////////////////////////////////////
/// \brief The ligoanalysis module looks at GW coincidences.
///
/// Given a time window, specified as an absolute time and a delta, it
/// searches for a set of interesting event types so that you can see
/// if there is a spike coincident with a gravitational wave event, or,
/// more broadly, any sort of external event you might think of. The
/// output is a set of histograms.
///
/// \author M. Strait
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"

#include "Geometry/Geometry.h"

#include "DAQDataFormats/RawTriggerTime.h"
#include "DAQDataFormats/RawEvent.h"
#include "DAQDataFormats/RawTrigger.h"
#include "DAQDataFormats/RawTriggerMask.h"
#include "DAQDataFormats/RawDataBlock.h"
#include "StandardRecord/SREnums.h"

#include "RawData/FlatDAQData.h"
#include "RawData/RawTrigger.h"

#include "RecoBase/Track.h"

#include "MCCheater/BackTracker.h"

#include "CMap/service/DetectorService.h"
#include "CelestialLocator/CelestialLocator.h"
#include "NovaTimingUtilities/TimingUtilities.h"

#include "TH1.h"
#include "TTree.h"
#include "TRandom3.h"

#include <string>
#include <vector>
#include <set>
#include <algorithm>

#include <signal.h>

#include "progress.cpp"

#include "func/timeutil.h"
#include "func/DaqChannelMask.h"


// "`-._,-'"`-._,-'"`-._,-' BEGIN sky map stuff "`-._,-'"`-._,-'"`-._,-'
#include "healpix_base.h"
#include "healpix_map.h"
#include "healpix_map_fitsio.h"

#include "alm.h" // Alm<T>
#include "alm_healpix_tools.h" // alm2map(), etc.
#include "alm_powspec_tools.h" // smoothWithGauss()

// Define this if you want to print extra information
//#define LOUD

// The sky map from LIGO/Virgo, if available and necessary, i.e. if we
// are analyzing events with pointing.  Two copies, the first smeared
// with our pointing resolution and the other smeared also with a larger
// angle to catch low energy numuCC-like events.
static const unsigned int npointres = 2;
static Healpix_Map<float> * healpix_skymap[npointres] = { NULL };

// The critical probability density values in the sky maps above which
// we are in the 90% confidence level region.  Set in beginRun().
static double skymap_crit_val[npointres] = { 0 };
// "`-._,-'"`-._,-'"`-._,-'  END sky map stuff  "`-._,-'"`-._,-'"`-._,-'


static const int TDC_PER_US = 64;
static const int US_PER_MICROSLICE = 50; // I hope this is always true

// Rough effective light speed in fiber
//
// Note: PhotonTransport uses an index of 1.59 -> 18.9cm/ns. But for
// rays at just the critical angle (see doc-2665), it is effectively
// 16.9cm/ns. So the mean speed is 17.9cm/ns. PhotonTransport adds a
// time from a histogram to the time gotten from the index. Does that
// model the angular effect?
static const float lightspeed = 17.; // cm/ns
static const float invlightspeed = 1/lightspeed; // ns/cm

// z extent of a plane (slightly wrong at block boundaries, but for
// my purposes, it doesn't matter)
static const double plnz = 6.6479;

// Set from the FCL parameters
static double gwevent_unix_double_time = 0;
static double needbgevent_unix_double_time = 0;
static long long window_size_s = 1000;

// Relaxed for the supernova analysis
// Were: 85, 600, 107, 2500
//
// At least for the FD, I want to apply the cuts later, so don't
// cut here.
static const int16_t fd_high_adc =10000;
static const int16_t nd_high_adc = 2500;

// Smallest ADC to save a single hit as a cluster. Without a cut here,
//the output is extremely large.
//
// About 2/3 of the FD signal with any hits has only one hit. 70% of
// that is above 65 ADC.
//
// Not clear what the best value is here
static const int16_t MINSINGLETONADC = 130;

// Generous box for regular Michel rejection. For stealth Michels,
// we'll rely on the huge slice box, so there's no need for this
// to 100% efficient (indeed, it's impossible for it to be).
const int trk_pln_buf = 5, trk_cel_buf = 9;

// Box to put around each slice for determining if we are near it.
// Not optimized.  In ND, probably just throw out whole detector.
const int slc_pln_buf = 24;
const int slc_cel_buf = 48;


static int gDet = caf::kUNKNOWN;

// Types of analysis, dependent on which trigger we're looking at
enum analysis_class_t { NDactivity, DDenergy, MinBiasFD, MinBiasND,
                        MichelFD, SNonlyFD, SNonlyND,
                        Blind, MAX_ANALYSIS_CLASS
} static analysis_class = MAX_ANALYSIS_CLASS;

// Hit information, distilled from CellHit, plus more info
struct mhit{
  unsigned int hitid; // index into the noise slice hit array
  float tns; // fine timing, in nanoseconds
  float tposoverc; // transverse position divided by speed of light
  float pe; // photoelectrons
  int16_t adc;
  int16_t plane;
  bool isx;
  int cell;
  bool used; // has this hit been used in a cluster yet?
  bool paired; // has this hit been used in a cluster of size 2+ yet?
  int truepdg; // For MC, what the truth is
  float trueE; // for MC, what the true initial particle kinetic energy is
  int16_t trueplane;
  int16_t truecellx;
  int16_t truecelly;

  // The minimum time to the next nearby track end
  double totrkend_s;
  // The time since the previous nearby track end
  double sincetrkend_s;

  // The minimum time to the next slice overlapping this hit in space.
  double toslice_s;
  // The time since the previous slice overlapping this hit in space.
  double sinceslice_s;

  // The minimum time to the next slice anywhere
  double toanyslice_s;
  // The time since the previous slice anywhere
  double sinceanyslice_s;

  // The time since the previous big shower *anywhere*
  double sincelastbigshower_s;

  // Whether DaqChannelMask has marked this channel as noisy
  bool noisy;
};

struct sncluster{
  std::vector<mhit *> hits;
  // I thought I'd store more in here, but maybe not and it should be a
  // typedef?
};

// Minimal slice information, distilled from rb::Cluster
struct mslice{
  float mintns, maxtns;

  // Extents with buffers
  int16_t bminplane, bmaxplane;
  int16_t bmincellx, bmaxcellx, bmincelly, bmaxcelly;

  float sliceduration;

  // Used as a proxy for "big shower"
  bool longslice;

  // From rb::Cluster::TotalADC()
  double totaladc;

  // From rb::Cluster::MeanX(), etc.
  double meanx, meany, meanz;

  // Time in integer seconds and nanoseconds since Jan 1, 1970
  uint32_t time_s;
  uint32_t time_ns;

  // Index into a reduced set of slices where if two overlap in time,
  // they are considered one.  This index is not contiguous.
  int mergeslice;
};

// Minimal track information, distilled from rb::Track, for use
// in correlating later activity to track ends
struct mtrack{
  // The plane the track ended. If this was an X plane, the X cell the
  // track ended, and the estimated nearest Y cell. If it was a Y plane,
  // vice versa.
  int16_t endplane, endcellx, endcelly;

  // Whether the last plane is X.  This is redundant with 'endplane',
  // but I can never remember whether X planes are even or odd, and this
  // massively reduces the potential for error.
  bool lastisx;

  // True if the rb::Track end point doesn't seem to match the collection of
  // hits.
  bool endconfused;

  float tns;
} mtrackMCMLXXXI;

class ligoanalysis : public art::EDProducer {
  public:
  explicit ligoanalysis(fhicl::ParameterSet const& pset);
  virtual ~ligoanalysis() { }; // compiles, but does not run, without this
  void produce(art::Event& evt);

  void beginJob();

  void beginRun(art::Run& run);

  /// \brief User-supplied type of trigger being examined.
  ///
  /// Which histograms we make depends on this.
  analysis_class_t fAnalysisClass;

  /// \brief The user-supplied time of the gravitational wave burst, or
  /// whatever time you want to center the search on.
  ///
  /// Expressed in RFC-3339 format, always in UTC, always with a Z at
  /// the end (to emphasize that it is UTC).
  std::string fGWEventTime;

  // If we want to measure the skymap-dependent background for a gravitational
  // wave at fNeedBGEventTime, set this to the time of the GW, and set
  // fGWEventTime to the time we're using as a background sample.  Otherwise,
  // leave this as the empty string.
  std::string fNeedBGEventTime;

  /// The file name of the LIGO/Virgo skymap, or the empty string if this
  /// analysis does not care about pointing or no map is available.
  std::string fSkyMap;

  /// \brief The user-supplied length of the search window in seconds.
  ///
  /// We search for half this length on either side of the given time.
  float fWindowSize;

  /// \brief Whether to cut ND events with multiple slices
  ///
  /// This is an effort to remove NuMI events that sneak past other filters.
  bool fCutNDmultislices;
};

/**********************************************************************/
/*                             Histograms                             */
/**********************************************************************/

struct ligohist{
  // Things per second, for whatever it is we're selecting
  TH1D * sig = NULL;

  // Livetime per second.  This means the amount of time read out by
  // triggers in the analysis window, not the amount of time a trigger
  // was enabled in the DAQ or anything similar.
  TH1D * live = NULL;

  // Base name for the histograms.  Name for the livetime histogram will be
  // this with "live" appended.
  std::string name;

  // Answers the question "Are livetime histograms meaningful for this?"
  //
  // They aren't meaningful if the sort of event we're looking at causes
  // the trigger to fire. For instance, it doesn't mean anything to look
  // at the livetime provided by the Upmu trigger when looking at a
  // count of upward going muons. However, it may be meaningful to look
  // at the livetime of Upmu if looking at low-energy events that happen
  // to get read out because of Upmu.
  //
  // If false, we won't write out the second histogram. In some cases,
  // whether or not they are meaningful depends on the trigger that is
  // using this histogram. We'll enable this if it is at least sometimes
  // meaningful. (Or would it be better to make different histograms and
  // avoid ever having irrelevant output?)
  bool dolive;

  ligohist(const std::string & name_, const bool dolive_)
  {
    name = name_;
    dolive = dolive_;
  }

  ligohist() {}
};

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
static bool uniquedata_tdc(const double livetime, const int tdc)
{
  if(livetime < 0.005) return true;
  return tdc >= 0 && tdc < 5000 * TDC_PER_US;
}

// Same as uniquedata_tns but for a floating point time. This is for
// deciding whether to keep tracks and clusters, and is necessarily
// sloppier than hit timing. But probably slices and tracks are
// reconstructed exactly the same way most of the time regardless of
// adjacent microslices, so we should accept them exactly once. It's
// probably possible to end up accepting the same track twice (or zero
// times), though, in unlucky cases.
__attribute__((unused)) static bool uniquedata_tns(const double livetime,
                                                   const double tns)
{
  if(livetime < 0.005) return true;
  return tns >= 0 && tns < 5e6;
}

// Construct the histograms and hook them up with TFileService to be saved
static void init_lh(ligohist & lh)
{
  if(lh.sig != NULL || lh.live != NULL){
    fprintf(stderr, "%s already initialized.\n", lh.name.c_str());
    exit(1);
  }

  art::ServiceHandle<art::TFileService> t;

  lh.sig = t->make<TH1D>(lh.name.c_str(), "",
    window_size_s, -window_size_s/2, window_size_s/2);

  if(lh.dolive)
    lh.live = t->make<TH1D>(Form("%slive", lh.name.c_str()), "",
      window_size_s, -window_size_s/2, window_size_s/2);
}

static void init_lh_name_live(ligohist & lh, const std::string & name,
                              const bool dolive)
{
  lh.name = name;
  lh.dolive = dolive;
  init_lh(lh);
}

/**********************************************************************/

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

/*
  Get the length of the event in TDC ticks, typically 550*64, and
  "delta_tdc", the time between the trigger time and the time the event
  starts. You can subtract this off of the time that the offline gives
  each hit to get the time since the beginning of the readout, and with
  the event length, the time until the end of the readout.

  delta_tdc is a signed 64 bit integer, even though it should always be
  a small positive number, just in case. Ditto for the event length.

  Returns whether this information was successfully extracted.
*/
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

/*********************************************************************/
/************************* Begin histograms **************************/
/*********************************************************************/

// Count of triggers, with no examination of the data within
static ligohist lh_rawtrigger("rawtrigger", true);

// Number of hits that are not in any Slicer4D slice, i.e. they are in the
// Slicer4D noise slice.  And away from slices.
static ligohist lh_unsliced_hits("unslicedhits", true);

// My experience with my neutron capture analysis in the ND is that
// things become well-behaved at about this level.
static const double bighit_threshold_nd = 35; // PE
static const double bighit_threshold_fd = 30;
static       double bighit_threshold    =  0; // set when we know det

// Number of unsliced hits that are over the above PE threshold
static ligohist lh_unsliced_big_hits("unslicedbighits", true);

// Number of unsliced hits that are over some PE threshold and are paired
// with a hit in an adjacent plane.  (Not an adjacent cell, because then
// we select lots of pairs from noisy modules.  We could build a noise map
// to fix that, but it would require two passes through the data.)
static ligohist lh_supernovalike("supernovalike", true);

// Count of slices with nothing around the edges, regardless of what sorts of
// objects are inside.
static ligohist lh_contained_slices("contained_slices", true);

// Number of triggers above two cuts for DDEnergy, the first pair
// in raw ADC, the second ADC per unit time.
static ligohist lh_ddenergy_locut("energy_low_cut", false);
static ligohist lh_ddenergy_hicut("energy_high_cut", false);
static ligohist lh_ddenergy_vhicut("energy_vhigh_cut", false);
static ligohist lh_ddenergy_lopertime("energy_low_cut_pertime", false);
static ligohist lh_ddenergy_hipertime("energy_high_cut_pertime", false);
static ligohist lh_ddenergy_vhipertime("energy_vhigh_cut_pertime", false);

/********************** Histogram with pointing ***********************/

// Raw number of tracks
static ligohist lh_tracks("tracks", true);
static ligohist lh_tracks_point[npointres];

// Number of tracks with at least one contained endpoint, with some
// additional sanity checks.
static ligohist lh_halfcontained_tracks("halfcontained_tracks", true);
static ligohist lh_halfcontained_tracks_point[npointres];

// Number of tracks with at two contained endpoints, with some additional
// sanity checks.
static ligohist lh_fullycontained_tracks("fullycontained_tracks", true);
static ligohist lh_fullycontained_tracks_point[npointres];

// Number of tracks that pass the Upmu analysis
static ligohist lh_upmu_tracks("upmu_tracks", true);
static ligohist lh_upmu_tracks_point[2];

// Just report livetime
static ligohist lh_blind("blind", true);

/*********************************************************************/
/**************************  End histograms **************************/
/*********************************************************************/

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

  // Number of nanoseconds from the first hit to the last hit. Large
  // times may indicate background. WARNING: We've seen that the Monte
  // Carlo is overly optimistic about how tight the detector timing is.
  // I have added an ad hoc correction which makes it look more like the
  // data, but it probably isn't super reliable.
  float timeext_ns;

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
  double totrkend_s;
  double sincetrkend_s;

  // The time until and since the last slice that overlaps this cluster
  // in space, plus a buffer of several cells and planes.  If this
  // hit is during such a slice, both are zero.  If there is no such
  // slice for one or the other case, that one is set to 1e9.
  double toslice_s;
  double sinceslice_s;

  // The time until and since the last slice anywhere. If this cluster is
  // during such a slice, both are zero. If there is no such slice for
  // one or the other case, that one is set to 1e9.
  double toanyslice_s;
  double sinceanyslice_s;

  // Time in seconds since the last "big" cosmic ray shower. This is
  // of interest because big showers can dump many (thousands) of
  // neutrons into the detector, which subsequently capture over the
  // following 100s of microseconds. It's not clear which, if any of the
  // following are relevant to NOvA, but they can also produce (1) B-12
  // and N-12, high energy beta decay isotopes with lifetimes in the
  // 10s of milliseconds (2) cosmogenic isotopes such as Li-9 that have
  // decays that mimic supernova IBDs and have lifetimes in the 100s
  // of milliseconds (3) other isotopes with lower energy simple beta
  // decays, but even longer lifetimes.
  //
  // A "big" shower is defined as one with a length over 2 microseconds.
  // This works, but is probably not so much a physical time duration as
  // it is using the fact that large energy depositions cause "flasher"
  // hits over the next few microseconds, a detector artifact.  It may also
  // be connecting the beginning of the wave of neutron captures or muons
  // truly separated in time by up to ~1-2 microseconds.
  //
  // If there was no previous big shower, this is set to the number of
  // seconds between the cluster and the Unix Epoch (~2e9).
  //
  // In my short tests so far, it looks like this variable does not
  // help identify background events that have not already been cut by
  // removing hits spatially close to big showers for 200 microseconds.
  // We could play with the definitions of "big shower" and the size of
  // this hard cut.
  double sincebigshower_s;

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

  // Fraction of the hits in this cluster that are from noisy channels
  float noisefrac;

} sninfo;

/*********************************************************************/

static void init_mev_stuff()
{
  init_lh(lh_unsliced_hits);
  init_lh(lh_unsliced_big_hits);
  init_lh(lh_supernovalike);

  art::ServiceHandle<art::TFileService> t;
  sntree = t->make<TTree>("sn", "");

  // These reduce the nonsense of typos that goes with TTree:Branch()
  #define Q(x) #x
  #define BRN(x, t) sntree    ->Branch(Q(x), &sninfo.x, Q(x/t))
  BRN(nhit,             I);
  BRN(truefrac,         F);
  BRN(truepdg,          I);
  BRN(trueE,            F);
  BRN(trueplane,        S);
  BRN(truecellx,        S);
  BRN(truecelly,        S);
  BRN(time_s,           i);
  BRN(time_ns,          i);
  BRN(totrkend_s,       D);
  BRN(sincetrkend_s,    D);
  BRN(toslice_s,        D);
  BRN(sinceslice_s,     D);
  BRN(toanyslice_s,     D);
  BRN(sinceanyslice_s,  D);
  BRN(sincebigshower_s, D);
  BRN(pex,              F);
  BRN(pey,              F);
  BRN(unattpe,          F);
  BRN(minhitadc,        S);
  BRN(timeext_ns,       F);
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
  BRN(noisefrac,        F);
}

static void init_blind_hist()
{
  init_lh(lh_blind);
}

static void init_track_and_contained_hists()
{
  init_lh(lh_tracks);
  init_lh(lh_halfcontained_tracks);
  init_lh(lh_fullycontained_tracks);
  init_lh(lh_contained_slices);

  for(unsigned int q = 0; q < npointres; q++){
    init_lh_name_live(lh_tracks_point[q],
                      Form("tracks_point_%d", q), true);
    init_lh_name_live(lh_halfcontained_tracks_point[q],
                      Form("halfcontained_tracks_point_%d", q), true);
    init_lh_name_live(lh_fullycontained_tracks_point[q],
                      Form("fullycontained_tracks_point_%d", q), true);
  }
}

// Probabilities for a point in the sky map. One is the raw value and
// the second is scaled by area, i.e. multiplied by sin(theta)
struct ac_raw{
  double raw;
  double area_corrected;
};

static bool compare_ac_raw(const ac_raw & a, const ac_raw & b)
{
  return a.raw > b.raw;
}

// Takes a skymap and smears it with a Gaussian to take into account our
// detector resolution.
//
// Mostly copied from the example code in the Healpix package in
// smoothing_cxx_module.cc
//
// I'd prefer to return a smeared map from the const& argument, but I
// can't immediately see how to manage that.
static void smear_skymap(Healpix_Map<float> * map, const double degrees)
{
  // Some power of two between 32 and 1024, by grepping.
  // Is this a good value?
  const int nlmax = 512;

  // 3 or 5 in all examples, more often 3.
  const int num_iter = 3;

  const tsize nmod = map->replaceUndefWith0();
  if(nmod!=0)
    printf("smear_skymap() WARNING: replaced %lu undefined map pixels "
           "with a value of 0\n", nmod);

  const float avg = map->average();
  map->Add(-avg);

  // I *think* the LIGO/Virgo skymaps are not weighted, so this is
  // easier than using get_right_weights(), which requires a paramfile.
  arr<double> weight;
  weight.alloc(2*map->Nside());
  weight.fill(1);

  Alm<xcomplex<float> > alm(nlmax, nlmax);
  if(map->Scheme() == NEST) map->swap_scheme();

  map2alm_iter(*map, alm, num_iter, weight);

  // What really is the best resolution to use? We don't have a solid
  // number, and obviously it is a function of energy, too.
  const double fwhm = degrees  * 2.355 * M_PI/180.;

  smoothWithGauss(alm, fwhm);
  alm2map(alm, *map);

  map->Add(avg);
}

static double find_critical_value(const int q)
{
  double sumprob = 0;

  // We need to integrate the probability and find the critical value
  // above which we are in the 90% CL region. The convenient way to
  // retrieve probabilities from the map is interpolated values at
  // unevenly spaced points. Scale these by their effective area to do
  // the integration, but set the critical value using unscaled values.
  std::vector<ac_raw> vals;

  int ni = 1000, nj = 1000;
  for(int i = 1; i < ni; i++){
    for(int j = 0; j < nj; j++){
      const double theta = i*M_PI/ni,
                   phi   = j*2.*M_PI/nj;
      const float val = healpix_skymap[q]->interpolated_value(
        pointing(theta,//dec: except this is 0 to pi, and dec is pi/2 to -pi/2
                 phi));//ra: as normal: 0h, 24h = 2pi
      sumprob += val * sin(theta);
      ac_raw new_ac_raw;
      new_ac_raw.raw = val;
      new_ac_raw.area_corrected = val*sin(theta);
      vals.push_back(new_ac_raw);
    }
  }

  sort(vals.begin(), vals.end(), compare_ac_raw);

  const float CL = 0.9; // 90% confidence level

  double crit = 0;

  float acc = 0;
  for(unsigned int i = 0; i < vals.size(); i++){
    acc += vals[i].area_corrected/sumprob;
    if(acc > CL){
      crit = vals[i].raw;
      break;
    }
  }

  // Print the map to the screen just so we know something is happening
  for(int which = 0; which < 2; which++){
    printf("Sky map %.0f%% region, in %s:\n",
           CL*100, which == 0?"ra/dec":"zen/azi");
    const int down = 40;
    for(int i = 0; i < down; i++){
      const double theta = i*M_PI/down;
      const int maxacross = 2*down;
      const int across = 2*int(down*sin(theta) + 0.5);

      for(int j = 0; j < (maxacross - across)/2; j++)
        printf(" ");

      for(int j = across-1; j >= 0; j--){
        const double phi = j*2*M_PI/across;

        // If which == 0, then theta and phi are ra and dec.  Otherwise, they
        // are the zen and azi, and we need to find the ra and dec.
        double ra = 0, dec = 0;
        if(which == 1){
          art::ServiceHandle<locator::CelestialLocator> celloc;
          celloc->GetRaDec(theta,phi,
                           (time_t)(needbgevent_unix_double_time?
                                    needbgevent_unix_double_time:
                                        gwevent_unix_double_time),
                           ra,dec);
        }

        const float val = healpix_skymap[q]->interpolated_value(
          which == 0?pointing(theta, phi)
                    :pointing(dec+M_PI_2, ra));
        printf("%c", val > crit?'X':'-');
      }
      printf("\n");
    }
  }

  return crit;
}

void ligoanalysis::beginJob()
{
}

void ligoanalysis::beginRun(art::Run& run)
{
  if(fAnalysisClass == DDenergy || fAnalysisClass == Blind ||
     fAnalysisClass == MichelFD ||
     fAnalysisClass == SNonlyFD || fAnalysisClass == SNonlyND) return;

  art::ServiceHandle<ds::DetectorService> ds;
  art::ServiceHandle<locator::CelestialLocator> celloc;
  celloc->SetDetector(ds->DetId(run));

  if(fSkyMap == "") return;

  for(unsigned int q = 0; q < npointres; q++){
    healpix_skymap[q] = new Healpix_Map<float>;
    printf("Opening FITS file %s\n", fSkyMap.c_str());
    try       { read_Healpix_map_from_fits(fSkyMap, *healpix_skymap[q]); }
    catch(...){ exit(1);                                                 }
  }

  printf("Smearing skymap\n");
  /* doc-16860 */
  smear_skymap(healpix_skymap[0], 1.3);

  /* doc-26828 */
  smear_skymap(healpix_skymap[1], acos(0.96) * 180/M_PI /* ~16 degrees */);

  for(unsigned int q = 0; q < npointres; q++)
    skymap_crit_val[q] = find_critical_value(q);
}

ligoanalysis::ligoanalysis(fhicl::ParameterSet const& pset) : EDProducer(),
  fGWEventTime(pset.get<std::string>("GWEventTime")),
  fNeedBGEventTime(pset.get<std::string>("NeedBGEventTime")),
  fSkyMap(pset.get<std::string>("SkyMap")),
  fWindowSize(pset.get<unsigned long long>("WindowSize")),
  fCutNDmultislices(pset.get<bool>("CutNDmultislices"))
{
  const std::string analysis_class_string(
    pset.get<std::string>("AnalysisClass"));

  if     (analysis_class_string == "NDactivity") fAnalysisClass = NDactivity;
  else if(analysis_class_string == "DDenergy")   fAnalysisClass = DDenergy;
  else if(analysis_class_string == "MinBiasFD")  fAnalysisClass = MinBiasFD;
  else if(analysis_class_string == "MinBiasND")  fAnalysisClass = MinBiasND;
  else if(analysis_class_string == "MichelFD")   fAnalysisClass = MichelFD;
  else if(analysis_class_string == "SNonlyFD")   fAnalysisClass = SNonlyFD;
  else if(analysis_class_string == "SNonlyND")   fAnalysisClass = SNonlyND;
  else if(analysis_class_string == "Blind")      fAnalysisClass = Blind;
  else{
    fprintf(stderr, "Unknown AnalysisClass \"%s\" in job fcl. See list "
            "in ligoanalysis.fcl.\n", analysis_class_string.c_str());
    exit(1);
  }

  analysis_class = fAnalysisClass; // expose to static functions

  gwevent_unix_double_time = rfc3339_to_unix_double(fGWEventTime);
  needbgevent_unix_double_time = fNeedBGEventTime==""?0
                                 :rfc3339_to_unix_double(fNeedBGEventTime);
  window_size_s = fWindowSize;

  switch(fAnalysisClass){
    case NDactivity:
      init_track_and_contained_hists();
      break;
    case DDenergy:
      init_lh(lh_rawtrigger);
      init_lh(lh_ddenergy_locut);
      init_lh(lh_ddenergy_hicut);
      init_lh(lh_ddenergy_vhicut);
      init_lh(lh_ddenergy_lopertime);
      init_lh(lh_ddenergy_hipertime);
      init_lh(lh_ddenergy_vhipertime);
      break;
    case MinBiasFD:
      init_track_and_contained_hists();
      init_lh(lh_upmu_tracks);
      for(unsigned int q = 0; q < npointres; q++)
        init_lh_name_live(lh_upmu_tracks_point[q],
                          Form("upmu_tracks_point_%d",q), false);
      // Fall-through
    case SNonlyFD:
    case MichelFD:
      init_lh(lh_rawtrigger);
      init_mev_stuff();
      break;
    case MinBiasND:
      init_track_and_contained_hists();
      // Fall-through
    case SNonlyND:
      init_lh(lh_rawtrigger);
      init_mev_stuff();
      break;
    case Blind:
      init_blind_hist();
      break;
    default:
      printf("No case for type %d\n", fAnalysisClass);
  }
}

/**********************************************************************/
/*                          The meat follows                          */
/**********************************************************************/

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

// Add 'sig' to the signal and 'live' to the livetime in bin number 'bin'.
static void THplusequals(ligohist & lh, int bin, const double sig,
                         const double live)
{
  // Ensure out-of-range falls into the overflow bins rather than being lost
  if(bin < 0) bin = 0;
  if(bin > lh.sig->GetNbinsX()+1) bin = lh.sig->GetNbinsX()+1;

  // Use SetBinContent instead of Fill(x, weight) to avoid having to look up
  // the bin number twice.
  lh.sig->SetBinContent(bin, lh.sig->GetBinContent(bin) + sig);
  if(lh.dolive)
    lh.live->SetBinContent(bin, lh.live->GetBinContent(bin) + live);
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

  if(verbose) printf("Accepted delta = %f seconds\n", delta);

  return floor(delta) + window_size_s/2
         + 1; // stupid ROOT 1-based numbering!
}

// Returns true if a point is "contained" for purposes of deciding if a track
// or physics slice is contained.
static bool contained(const TVector3 & v)
{
  if(gDet == caf::kNEARDET)
    return fabs(v.X()) < 150 && fabs(v.Y()) < 150 &&
           v.Z() > 40 && v.Z() < 1225;
  if(gDet == caf::kFARDET)
    return fabs(v.X()) < 650 &&
      v.Y() < 500 && v.Y() > -650 &&
      v.Z() > 75 && v.Z() < 5900;
  fprintf(stderr, "Unknown detector %d\n", gDet);
  exit(1);
}

// Returns true if the track enters and exits.  Or if both of its ends
// are pretty close to the edge.  But not if it is just steep.
static bool un_contained_track(const rb::Track & t)
{
  return !contained(t.Start()) && !contained(t.Stop());
}

// Check that the reconstruction (probably BreakPointFitter) doesn't think
// either end of the track points nearly along a plane AND check that if we
// just draw a line from one end of the track to the other, that that doesn't
// either.  This second part is to make sure that some tiny kink at the start
// or end can't fool us into thinking it is well-contained.
static bool good_track_direction(const rb::Track & t)
{
  static const double tan_track_cut = atan(15 * M_PI/180);

  const double tot_dx = t.Start().X() - t.Stop().X();
  const double tot_dy = t.Start().Y() - t.Stop().Y();
  const double tot_dz = t.Start().Z() - t.Stop().Z();

  const double rec_dx1 = t.    Dir().X();
  const double rec_dx2 = t.StopDir().X();
  const double rec_dy1 = t.    Dir().Y();
  const double rec_dy2 = t.StopDir().Y();
  const double rec_dz1 = t.    Dir().Z();
  const double rec_dz2 = t.StopDir().Z();

  return fabs(tot_dz /tot_dx ) > tan_track_cut &&
         fabs(tot_dz /tot_dy ) > tan_track_cut &&
         fabs(rec_dz1/rec_dx1) > tan_track_cut &&
         fabs(rec_dz2/rec_dx2) > tan_track_cut &&
         fabs(rec_dz1/rec_dy1) > tan_track_cut &&
         fabs(rec_dz2/rec_dy2) > tan_track_cut;
}

// Fewest planes we'll accept a track having, or in the case of
// fully-contained tracks, the fewest *contiguous* planes.
static const int min_plane_extent = 10;

// Returns true if the track starts AND stops inside the detector,
// and we're pretty sure that this isn't artefactual.
static bool fully_contained_track(const rb::Track & t)
{
  // Note how here we're looking at all the hits' positions, not
  // just the line whose endpoints are Start() and Stop().  For steep
  // tracks, this reveals their uncontainedness.
  return contained(t.MaxXYZ()) && contained(t.MinXYZ()) &&
          good_track_direction(t) &&

          // Don't accept a track just because of a stray hit that gets
          // swept into it.
          t.MostContiguousPlanes(geo::kXorY) >= min_plane_extent;
}

// Returns true if either end of the track is uncontained or might be
static bool half_uncontained_track(const rb::Track & t)
{
  return !contained(t.Start()) || !contained(t.Stop()) ||
          !good_track_direction(t) ||
          t.ExtentPlane() < min_plane_extent;
}

// Returns true if the track either starts or stops in the detector, or both,
// and we're pretty sure that this isn't artefactual.
static bool half_contained_track(const rb::Track & t)
{
  return (contained(t.Start()) || contained(t.Stop())) &&
          good_track_direction(t) &&
          t.ExtentPlane() >= min_plane_extent;
}

// return the index in the slice array that the given track is in.  But use
// my merged slice index if slices overlap.
static int which_slice_is_this_track_in(
  const rb::Track & t,
  const art::Handle< std::vector<rb::Cluster> > & slice,
  const std::vector<mslice> & sliceinfo)
{
  // I'm sure this is not the best way, but I have had it with trying to
  // figure out what the best way is, and I'm just going to do it *some*
  // way.
  if(t.NCell() == 0) return -1;
  const art::Ptr<rb::CellHit> ahit =
    t.Cell(0); // some random hit on the track

  // Skip slice 0 since it is the noise slice
  for(unsigned int i = 1; i < slice->size(); i++){
    const rb::Cluster & slc = (*slice)[i];
    for(unsigned int j = 0; j < slc.NCell(); j++){
      const art::Ptr<rb::CellHit> shit = slc.Cell(j);
      if(*ahit == *shit) return sliceinfo[i].mergeslice;
    }
  }
  return -1;
}

/************************************************/
/* Here come the functions that fill histograms */
/************************************************/

// Put put raw triggers into a histogram
static void count_triggers(const art::Event & evt)
{
  const double rawtime = rawlivetime(evt);

  // rawtime doesn't make sense for counting, e.g. DDEnergy triggers,
  // but it is a useful diagonistic for seeing if we're within a long
  // SNEWS/LIGO trigger.
  THplusequals(lh_rawtrigger, timebin(evt), 1, rawtime);
}

static void count_ddenergy(const art::Event & evt)
{
  art::Handle< std::vector<rawdata::RawDigit> > rawhits;
  getrawdigits(rawhits, evt);

  int64_t sumadc = 0;

  for(unsigned int i = 0; i < rawhits->size(); i++)
    sumadc += (*rawhits)[i].ADC();

  const double rawtime = rawlivetime(evt);

  printf("ADC: %ld %12.0f\n", sumadc, sumadc/rawtime);

  if(sumadc/rawtime > 5e10)
    THplusequals(lh_ddenergy_lopertime, timebin(evt), 1, rawtime);
  if(sumadc/rawtime > 5e11)
    THplusequals(lh_ddenergy_hipertime, timebin(evt), 1, rawtime);
  if(sumadc/rawtime > 5e12)
    THplusequals(lh_ddenergy_vhipertime, timebin(evt), 1, rawtime);

  if(sumadc > 5e6)
    THplusequals(lh_ddenergy_locut, timebin(evt), 1, rawtime);
  if(sumadc > 5e7){
    THplusequals(lh_ddenergy_hicut, timebin(evt), 1, rawtime);
    printf("Event passing high cut: event number %d\n", evt.event());
  }
  if(sumadc > 5e8){
    THplusequals(lh_ddenergy_vhicut, timebin(evt), 1, rawtime);
    printf("Event passing very high cut: event number %d\n", evt.event());
  }
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
static std::vector<mtrack> save_mtracks(const art::Event & evt)
{
  art::ServiceHandle<geo::Geometry> geo;

  art::Handle< std::vector<rb::Track> > tracks;
  evt.getByLabel("windowtrack", tracks);
  if(tracks->empty())
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
    }
    catch(...){
      // Assume this means the end is outside of the detector,
      // so we don't want it
      continue;
    }

    recoplane++; // because I shifted by one to find the orthogonal cell

    if(recoplane == endplane){
      if(endisx) minycell = recocell;
      else       endxcell = recocell;
    }

    // exiters
    if(endplane == 0 || endplane == 895) continue;
    if(minycell < 10) continue;
    if(endxcell < 5 || endxcell > 378) continue;

    struct mtrack thisone;

    thisone.endconfused = recoplane != endplane;
    thisone.endcelly = minycell;
    thisone.endcellx = endxcell;
    thisone.endplane = endplane;
    thisone.lastisx = endisx;

    thisone.tns = trk.MeanTNS();

    mtracks.push_back(thisone);
  }

  std::sort(mtracks.begin(), mtracks.end(), comparemtracktime);

  return mtracks;
}

// Builds list of distilled slice information
static std::vector<mslice> make_sliceinfo_list(const art::Event & evt,
  const art::Handle< std::vector<rb::Cluster> > & slice)
{
  std::vector<mslice> sliceinfo;

  // Start at 1 to skip "noise" slice. Slice indicies will all be off by
  // one, but makes it easier not to accidentally use the noise slice.
  for(unsigned int i = 1; i < slice->size(); i++){
    mslice slc;
    memset(&slc, 0, sizeof(mslice));

    if((*slice)[i].NCell() == 0) continue;

    slc.mintns = (*slice)[i].MinTNS();
    slc.maxtns = (*slice)[i].MaxTNS();

    const double LONGSLICEDURATION = 2000;
    slc.sliceduration = slc.maxtns - slc.mintns;
    slc.longslice = slc.sliceduration > LONGSLICEDURATION;

    slc.totaladc = (*slice)[i].TotalADC();

    slc.time_s = art_time_plus_some_ns(evt.time().value(),
                      (*slice)[i].MeanTNS()).first;
    slc.time_ns = art_time_plus_some_ns(evt.time().value(),
                      (*slice)[i].MeanTNS()).second;

    slc.meanx = (*slice)[i].MeanX();
    slc.meany = (*slice)[i].MeanY();
    slc.meanz = (*slice)[i].MeanZ();

    slc.bminplane = (*slice)[i].MinPlane() - slc_pln_buf;
    slc.bmaxplane = (*slice)[i].MaxPlane() + slc_pln_buf;
    if((*slice)[i].NCell(geo::kX)){
      slc.bmincellx = (*slice)[i].MinCell(geo::kX) - slc_cel_buf;
      slc.bmaxcellx = (*slice)[i].MaxCell(geo::kX) + slc_cel_buf;
    }
    if((*slice)[i].NCell(geo::kY)){
      slc.bmincelly = (*slice)[i].MinCell(geo::kY) - slc_cel_buf;
      slc.bmaxcelly = (*slice)[i].MaxCell(geo::kY) + slc_cel_buf;
    }

    if(i == 1 || i == 2)
      slc.mergeslice = i;
    else if(slc.mintns < sliceinfo[sliceinfo.size()-1].maxtns)
      slc.mergeslice = sliceinfo[sliceinfo.size()-1].mergeslice;
    else
      slc.mergeslice = i;

    // Take all slices, even if they are outside the 5ms window for
    // long trigger sub-triggers, because this is a list of things to
    // *exclude*.
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
  // of cases. Protect agains this. No apparent trouble with truecellx,
  // but protect it too.
  if(truecellx > 384) truecellx = -1;
  if(truecelly > 384) truecelly = -1;
}

// Helper function for count_mev(). Selects hits that
// are candidates to be put into hit pairs. Apply a low ADC cut if
// 'adc_cut'. Always apply a high ADC cut. This is intended to be true
// if we are searching for supernova-like events and false if we are
// searching for unmodeled bursts.
static std::vector<mhit> select_hits_for_mev_search(
  const rb::Cluster & noiseslice, const std::vector<mslice> & sliceinfo,
  const std::vector<mtrack> & trackinfo,
  const double livetime, const bool adc_cut, unsigned long long evttime)
{
  art::ServiceHandle<geo::Geometry> geo;

  // profiling indicates that it is helpful to save these.
  const unsigned int nslice = sliceinfo.size(), ntrack = trackinfo.size();

  std::vector<mhit> mhits;

  // Units are hits/event: cold threshold, hot threshold.
  // A normal channel has about 120Hz.
  // Andrey's more sophisticated treatment in
  // novaddt/OnlineCalibration/HotMapMaker_module.cc uses thresholds
  // of 1kHz sustained over 50ms (hot) and >90% of 50ms windows with
  // no hits (cold). We aren't looking at enough data to find cold
  // channels, so disable that. We're going to sample for 50ms,
  // so to avoid masking off channels that are up to 10x the mean
  // noisiness (~50Hz) 99% of the time, we need to set the threshold
  // at 1200Hz -- a 1200Hz channel has a mean of 60 hits in 50ms.  Assuming
  // Poisson, it has a 99% chance of < 80 hits in 50ms, so that's 8
  // hits per event.
  static sn::DaqChannelMask chmask(-1, 8.0);
  {
    static unsigned int evcollected = 0;
    if(evcollected < 10){
      for(unsigned int i = 0; i < noiseslice.NCell(); i++)
        chmask.AddHit(*(noiseslice.Cell(i)));
      evcollected++;
    }
    if(evcollected == 10){
      chmask.IncrementDuration(10 * 0.005);
      chmask.CalculateRates();
    }
  }

  const int16_t high_adc = gDet == caf::kNEARDET? nd_high_adc: fd_high_adc;

  for(unsigned int i = 0; i < noiseslice.NCell(); i++){
    const art::Ptr<rb::CellHit> & hit = noiseslice.Cell(i);

    if(!uniquedata_tdc(livetime, hit->TDC())) continue;

    if(hit->ADC() >= high_adc) continue;

    const int cell  = hit->Cell();
    const int plane = hit->Plane();

    // Cut only hit location, but don't do that for the supernova analysis,
    // because we'll figure out the best such cuts downstream.
    if(analysis_class != SNonlyFD && analysis_class != MichelFD &&
       analysis_class != SNonlyND){
      if(gDet == caf::kNEARDET){
        const int ndnplaneedge = 4;

        // Exclude whole muon catcher. Can't reasonably have a supernova
        // event hit planes on each side of a steel plane, and while
        // they might hit adjacent scintillator planes in the muon
        // catcher, I don't want to deal with the complications.
        if(plane <= ndnplaneedge) continue;
        if(plane >= 192-ndnplaneedge) continue;

        // At top of detector, stricter cut based on observed backgrounds
        if(hit->View() == geo::kY && cell >= 96 - 20) continue;

        // Drop outermost cells. This excludes most muon tracks that
        // just barely enter the detector, but don't get reconstructed.
        const int ndncelledge_x = 4;
        if(hit->View() == geo::kX &&
           (cell <= ndncelledge_x || cell >= 96 - ndncelledge_x)) continue;
      }
      else{ // far
        const int fdnplaneedge = 2;
        if(plane <= fdnplaneedge) continue;
        if(plane >= 895-fdnplaneedge) continue;

        if(hit->View() == geo::kY && cell >= 383 - 50)
          continue;

        const int fdncelledge_x = 10;
        if(hit->View() == geo::kX &&
           (cell <= fdncelledge_x || cell >= 383 - fdncelledge_x)) continue;
      }
    }

    static TRandom3 randfortiming;

    // Smear out MC timing as per my study shown in doc-45041
    const float tns = hit->TNS()
      + (hit->IsMC()?randfortiming.Gaus()*23.:0);

    mhit h;
    h.sincelastbigshower_s = 1e9;
    h.sincetrkend_s = 1e9;
    h.totrkend_s = 1e9;
    h.sinceslice_s = 1e9;
    h.toslice_s = 1e9;
    h.sinceanyslice_s = 1e9;
    h.toanyslice_s = 1e9;

    // How long since the last "big shower" (a.k.a. long slice) anywhere.
    // This looks back at least ~500us since we're looking at all slices
    // from this trigger and the previous one.
    //
    // Construct measures of how close we are to slices. Start at 0 because
    // we have not included the "noise" slice in sliceinfo.
    for(unsigned int j = 0; j < nslice; j++){
      if(!sliceinfo[j].longslice) continue;

      const double timesince_s = (tns - sliceinfo[j].maxtns)*1e-9;

      if(timesince_s > 0 && timesince_s < h.sincelastbigshower_s)
        h.sincelastbigshower_s = timesince_s;
    }

    for(unsigned int j = 0; j < nslice; j++){
      const double timesince_s = (tns - sliceinfo[j].maxtns)*1e-9;
      const double timeto_s    = (sliceinfo[j].mintns - tns)*1e-9;

      if(timesince_s > 0 && timesince_s < h.sinceanyslice_s)
        h.sinceanyslice_s = timesince_s;

      if(timeto_s > 0 && timeto_s < h.toanyslice_s)
        h.toanyslice_s = timeto_s;

      // It's in the middle of the time range of the slice
      if(timesince_s <= 0 && timeto_s <= 0)
        h.toanyslice_s = h.sinceanyslice_s = 0;

      // Now with a spatial restriction
      if(plane > sliceinfo[j].bminplane && plane < sliceinfo[j].bmaxplane){
        if(
           (hit->View() == geo::kX &&
            (cell > sliceinfo[j].bmincellx &&
             cell < sliceinfo[j].bmaxcellx))
           ||
           (hit->View() != geo::kX &&
            (cell > sliceinfo[j].bmincelly &&
             cell < sliceinfo[j].bmaxcelly))){

          if(timesince_s > 0 && timesince_s < h.sinceslice_s)
            h.sinceslice_s = timesince_s;

          if(timeto_s > 0 && timeto_s < h.toslice_s)
            h.toslice_s = timeto_s;

          if(timesince_s <= 0 && timeto_s <= 0){
            h.toslice_s = h.sinceslice_s = 0;
            // Safe to break because we've also set *anyslice above to 0.
            break;
          }
        }
      }
    }

    // Construct measures of how close we are to track ends
    for(unsigned int j = 0; j < ntrack; j++){
      const double since_s = (tns - trackinfo[j].tns)*1e-9;
      const double to_s    = (trackinfo[j].tns - tns)*1e-9;

      if(plane > trackinfo[j].endplane - trk_pln_buf &&
         plane < trackinfo[j].endplane + trk_pln_buf){
        if(
           (hit->View() == geo::kX &&
            (cell > trackinfo[j].endcellx - trk_cel_buf &&
             cell < trackinfo[j].endcellx + trk_cel_buf))
           ||
           (hit->View() != geo::kX &&
            (cell > trackinfo[j].endcelly - trk_cel_buf &&
             cell < trackinfo[j].endcelly + trk_cel_buf))){

          if(since_s > 0 && since_s < h.sincetrkend_s)
            h.sincetrkend_s = since_s;

          if(to_s > 0 && to_s < h.totrkend_s)
            h.totrkend_s = to_s;

          if(since_s <= 0 && to_s <= 0){
            h.totrkend_s = h.sincetrkend_s = 0;
            break;
          }
        }
      }
    }

    h.hitid = i; // Because CellHit::ID() seems to always be zero
    h.noisy = chmask.RatesCalculated() && chmask.ChannelIsMasked(*hit);
    h.used  = false;
    h.paired= false;
    h.tns   = tns;
    h.plane = plane;
    h.isx   = hit->View() == geo::kX;
    h.cell  = cell;
    h.tposoverc = geo->CellTpos(plane, cell) * invlightspeed;
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
      }
    }
    else{
      h.truepdg = 0;
      h.trueE = 0;
      h.trueplane = -1;
      h.truecellx = -1;
      h.truecelly = -1;
    }

    mhits.push_back(h);
  }

  return mhits;
}

// https://thedailywtf.com/articles/What_Is_Truth_0x3f_
enum boOl { falSe, trUe, past_plane_of_interest };

// Should we allow non-adjacent planes to get neutron hit clusters?
// Yes, it seems so. It looks like 2.6% of neutrons do the golden
// thing of having hits in adjacent planes, but another 3.8% have hits
// in two (or more) non-adjacent planes. Unfortunately, they're spread
// out over many planes, with 1.1% two planes separated, 0.4% three
// planes, 0.7% four planes, 0.3% five planes, etc. (I think even is
// favored over odd because then both are either on the bright side or
// the dark side instead of being random.)
//
// What's the risk of accepting hits too far away from each other? (1)
// Lots background making huge file sizes that will just have to be
// cut later (2) Spoiling signal hits by attaching them to background
// hits. Well, tests indicate that it about triples the file size. Since
// it also doubles the potential signal, that seems tolerable.
// We'll just have to see about (2).
static const int MAXPLANESAWAY = 10;

// Return trUe if this hit does cluster with the existing cluster,
// falSe if it does not, but further hits should be checked,
// and past_plane_of_interest if it does not and no more hits
// should be checked because the hits are sorted by plane number
// and this hit is too far downstream already.
static boOl does_cluster(const sncluster & clu, const mhit & h)
{
  // 1. Check if this hit is close enough in planes to an existing hit
  // in the cluster. If not, fail it. If there's another hit between
  // which would bridge them, that's fine because we will come back and
  // check this hit again.
  bool another_hit_few_enough_planes_away = false;
  int leastpast = 1000;
  for(unsigned int i = 0; i < clu.hits.size(); i++){
    if(abs(h.plane - clu.hits[i]->plane) <= MAXPLANESAWAY)
      another_hit_few_enough_planes_away = true;
    const int pastby = h.plane - clu.hits[i]->plane;
    if(pastby < leastpast) leastpast = pastby;
  }

  if(leastpast > MAXPLANESAWAY) return past_plane_of_interest;

  if(!another_hit_few_enough_planes_away) return falSe;

  // 2. Check that the hit is in time with the first hit of the cluster.
  // Maybe it should have to be in time with the mean time of the
  // cluster so far or something. I think checking against the first hit
  // is good enough.

  // Intentionally bigger than optimum value suggested by MC (150ns FD,
  // 10ns(!) ND) because we know that the data has a bigger time spread
  // from, for instance, the very well known slice duration discrepancy
  // (see, e.g., doc-19053).
  const float timewindow = 250; // ns

  if(clu.hits[0]->plane%2 == h.plane%2){
    // These hits are in the same view, so the times can be compared directly
    if(fabs(clu.hits[0]->tns - h.tns) > timewindow) return falSe;
  }
  else{
    // These hits are in different views, so for each, use the other's
    // transverse position to correct the time;
    const float
      time_1st_corr = clu.hits[0]->tns +            h.tposoverc,
      time_new_corr =            h.tns + clu.hits[0]->tposoverc;

    if(fabs(time_new_corr - time_1st_corr) > timewindow) return falSe;
  }

  // 3. Check that if this hit is in the same view as another hit
  // in the cluster, that they are nearby.  Otherwise, fail it.  This
  // introduces an order dependence, and in principle we could pick
  // a bad hit first, and then fail a subsequent good hit.  We could
  // take hits tentatively waiting for a better one to come along.  But
  // let's indefinitely shelve that idea for the day when the simple
  // approach doesn't seem sufficient, which it probably is.

  // Should be fairly large to admit neutron clusters, which are made
  // out of gammas with mean free path ~25cm
  const int MAXCELLDIST = 16;

  bool in_same_view_as_another_hit = false;
  bool close_to_another_hit_in_w = false;

  for(unsigned int i = 0; i < clu.hits.size(); i++)
    if(h.plane%2 == clu.hits[i]->plane%2)
      in_same_view_as_another_hit = true;

  for(unsigned int i = 0; i < clu.hits.size(); i++)
    if(h.plane%2 == clu.hits[i]->plane%2 &&
       abs(h.cell - clu.hits[i]->cell) <= MAXCELLDIST)
      close_to_another_hit_in_w = true;

  if(in_same_view_as_another_hit && !close_to_another_hit_in_w) return falSe;

  return trUe;
}

// For the given supernova cluster, return the energy of the particle that
// contributed the most pe-weighted hits.  Ignore hits with no truth
// information.
static float plurality_of_E(const sncluster & c)
{
  if(c.hits.empty()) return 0;

  std::map<float, float> m;
  for(const mhit * h : c.hits) if(h->trueE != 0) m[h->trueE] += h->pe;

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
// that contributed the most pe-weighted hits.
static int16_t plurality_of_trueplane(const sncluster & c)
{
  if(c.hits.empty()) return -1;

  std::map<int16_t, float> m;
  for(const mhit * h : c.hits)
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
// that contributed the most pe-weighted hits.
static int16_t plurality_of_truecellx(const sncluster & c)
{
  if(c.hits.empty()) return -1;

  std::map<int16_t, float> m;
  for(const mhit * h : c.hits)
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
// that contributed the most pe-weighted hits.
static int16_t plurality_of_truecelly(const sncluster & c)
{
  if(c.hits.empty()) return -1;

  std::map<int16_t, float> m;
  for(const mhit * h : c.hits)
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

static float noise_frac(const sncluster & c)
{
  unsigned int nn = 0;
  for(const mhit * h : c.hits) nn += h->noisy;
  return float(nn)/c.hits.size();
}

// For the given supernova cluster, return the true PDG id that contributed
// the most pe-weighted hits.  Ignore hits with no truth information.
static int plurality_of_truth(const sncluster & c)
{
  if(c.hits.empty()) return 0;

  std::map<int, float> m;
  for(const mhit * h : c.hits) if(h->truepdg != 0) m[h->truepdg] += h->pe;

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
  for(const mhit * h : c.hits) if(h->truepdg != 0) ntrue++;
  return (float)ntrue/c.hits.size();
}

static double mean_tns(const sncluster & c)
{
  double sum = 0;
  for(const mhit * h : c.hits) sum += h->tns;
  return sum / c.hits.size();
}

// Return <pe in x cells, pe in y cells>
static std::pair<float, float> sum_pe(const sncluster & c)
{
  std::pair<float, float> sum(0, 0);
  for(const mhit * h : c.hits) (h->isx?sum.first:sum.second) += h->pe;
  return sum;
}

static int min_plane(const sncluster & c)
{
  int ans = 9999; // There's always a plane, so result will never be 9999
  for(const mhit * h : c.hits) if(h->plane < ans) ans = h->plane;
  return ans;
}

static int max_plane(const sncluster & c)
{
  int ans = -1; // There's always a plane, so result will never be -1
  for(const mhit * h : c.hits) if(h->plane > ans) ans = h->plane;
  return ans;
}

static int plane_gap(const sncluster & c)
{
  std::vector<int> planes;
  for(const auto h : c.hits) planes.push_back(h->plane);

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
  for(const auto h : c.hits) if((h->isx ^ !x)) cells.push_back(h->cell);

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
  for(const auto h : c.hits) if((h->isx ^ !x) && h->cell < ans) ans = h->cell;

  // It's more convenient if the invalid value is always -1
  if(ans == 9999) ans = -1;
  return ans;
}

static int max_cell(const sncluster & c, const bool x)
{
  int ans = -1;
  for(const auto h : c.hits) if((h->isx ^ !x) && h->cell > ans) ans = h->cell;
  return ans;
}

static double to_trkend(const sncluster & c)
{
  // Previously used the mean instead of the smallest. That's dumb,
  // because one of the hits could be right on top of a trkend, but if
  // the other one is too far away in space to see that trkend, it will
  // dilute away that information.
  double least = 1e9;
  for(const auto h : c.hits)
    if(h->totrkend_s < least)
      least = h->totrkend_s;
  return least;
}

static double since_trkend(const sncluster & c)
{
  double least = 1e9;
  for(const auto h : c.hits)
    if(h->sincetrkend_s < least)
      least = h->sincetrkend_s;
  return least;
}

static double to_slice(const sncluster & c)
{
  double least = 1e9;
  for(const auto h : c.hits)
    if(h->toslice_s < least)
      least = h->toslice_s;
  return least;
}

static double since_slice(const sncluster & c)
{
  double least = 1e9;
  for(const auto h : c.hits)
    if(h->sinceslice_s < least)
      least = h->sinceslice_s;
  return least;
}

static double to_anyslice(const sncluster & c)
{
  double least = 1e9;
  for(const auto h : c.hits)
    if(h->toanyslice_s < least)
      least = h->toanyslice_s;
  return least;
}

static double since_anyslice(const sncluster & c)
{
  double least = 1e9;
  for(const auto h : c.hits)
    if(h->sinceanyslice_s < least)
      least = h->sinceanyslice_s;
  return least;
}

static double since_last_big_shower(const sncluster & c)
{
  double least = 1e9;
  for(const auto h : c.hits)
    if(h->sincelastbigshower_s < least)
      least = h->sincelastbigshower_s;
  return least;
}

static int16_t min_hit_adc(const sncluster & c)
{
  int16_t ans = SHRT_MAX;
  for(const auto h : c.hits) if(h->adc < ans) ans = h->adc;
  return ans;
}

static float time_ext_ns(const sncluster & c)
{
  double mintime = FLT_MAX, maxtime = FLT_MIN;
  if(c.hits.size() <= 1) return 0;
  for(const auto & h : c.hits){
    // Calculate all times as relative to the first hit.  Absolute time and
    // overall sign don't matter because we're finding the extent.
    float deltat = 0;
    if(c.hits[0]->plane%2 == h->plane%2)
      deltat = c.hits[0]->tns - h->tns;
    else{
      const float
        time_1st_corr = c.hits[0]->tns +         h->tposoverc,
        time_oth_corr =         h->tns + c.hits[0]->tposoverc;
      deltat = time_1st_corr - time_oth_corr;
    }

    if(deltat < mintime) mintime = deltat;
    if(deltat > maxtime) maxtime = deltat;
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

  // XXX FD only -- we may well want this for the ND as well
  const double celllength = 1650;

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

  if(c.hits.empty()){
    fprintf(stderr, "savecluster: got empty cluster, skipping\n");
    return;
  }

  sninfo.nhit = c.hits.size();
  sninfo.truefrac = fractrue(c);
  sninfo.truepdg = plurality_of_truth(c);
  sninfo.trueE = plurality_of_E(c);
  sninfo.trueplane = plurality_of_trueplane(c);
  sninfo.truecellx = plurality_of_truecellx(c);
  sninfo.truecelly = plurality_of_truecelly(c);
  sninfo.time_s = art_time_plus_some_ns(evt.time().value(),
                    mean_tns(c)).first;
  sninfo.time_ns = art_time_plus_some_ns(evt.time().value(),
                    mean_tns(c)).second;
  sninfo.pex = sum_pe(c).first;
  sninfo.pey = sum_pe(c).second;
  sninfo.minhitadc  = min_hit_adc(c);
  sninfo.timeext_ns = time_ext_ns(c);
  sninfo.minplane = min_plane(c);
  sninfo.maxplane = max_plane(c);
  sninfo.planegap = plane_gap(c);
  sninfo.planeextent = sninfo.maxplane - sninfo.minplane + 1;
  sninfo.mincellx = min_cell(c, true);
  sninfo.maxcellx = max_cell(c, true);
  sninfo.cellxgap = cell_gap(c, true);
  sninfo.cellxextent = sninfo.maxcellx == -1? -1:
                       sninfo.maxcellx - sninfo.mincellx + 1;
  sninfo.mincelly = min_cell(c, false);
  sninfo.maxcelly = max_cell(c, false);
  sninfo.cellygap = cell_gap(c, false);
  sninfo.cellyextent = sninfo.maxcelly == -1? -1:
                       sninfo.maxcelly - sninfo.mincelly + 1;

  sninfo.sincetrkend_s = since_trkend(c);
  sninfo.   totrkend_s =    to_trkend(c);

  sninfo.toslice_s = to_slice(c);
  sninfo.sinceslice_s = since_slice(c);

  sninfo.toanyslice_s = to_anyslice(c);
  sninfo.sinceanyslice_s = since_anyslice(c);

  sninfo.sincebigshower_s = since_last_big_shower(c);

  sninfo.run = evt.run();
  sninfo.subrun = evt.subRun();
  sninfo.event = evt.event();
  sninfo.hitid = c.hits[0]->hitid;
  sninfo.noisefrac = noise_frac(c);

  calibrate(sninfo); // set unattpe

  sntree->Fill();
}

static bool compare_plane(const mhit & a, const mhit & b)
{
  return a.plane < b.plane;
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

// Search for supernova-like events if 'supernovalike'.  Otherwise
// search for the sum of supernova like events and unpaired hits, once
// with a PE cut and once without.
static void count_mev(const art::Event & evt, const bool supernovalike,
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
  // Return value is not const because we modify ::used below.
  std::vector<mhit> mhits =
    select_hits_for_mev_search((*slice)[0], sliceinfo, trackinfo, livetime,
      supernovalike, evt.time().value());

  std::sort(mhits.begin(), mhits.end(), compare_plane);

  unsigned int hitclusters = 0;

  std::vector<sncluster> snclusters;

  // Starting with every eligible hit, form a cluster (possibly of size one)
  const unsigned int nmhits = mhits.size(); // really seems to help speed
  for(unsigned int i = 0; i < nmhits; i++){
    if(mhits[i].used) continue;

    sncluster clu;
    mhits[i].used = true;
    clu.hits.push_back(&mhits[i]);

    bool done = false;
    do{
      done = true;

      // mhits is sorted by plane number. Find the first other hit that
      // has a plane number close enough that we'd use it. (Must find
      // something, because this is satisfied by the first hit added to
      // the cluster.)
      const unsigned int startat =
        std::lower_bound(mhits.begin(), mhits.end(),
                         hitwithplane(min_plane(clu) - MAXPLANESAWAY),
                         compare_plane)
        - mhits.begin();

      for(unsigned int j = startat; j < nmhits; j++){
        if(mhits[j].used) continue;

        const boOl res = does_cluster(clu, mhits[j]);
        if(res == past_plane_of_interest) break;
        if(res == falSe) continue;

        mhits[i].paired = mhits[j].paired = mhits[j].used = true;
        clu.hits.push_back(&mhits[j]);
        done = false;
      }
    }while(!done);

    if((clu.hits.size() >= 2 && clu.hits.size() <= 7) ||
       (clu.hits.size() == 1 && clu.hits[0]->adc >= MINSINGLETONADC)){
      snclusters.push_back(clu);
      hitclusters++;
    }
  }

  std::sort(snclusters.begin(), snclusters.end(), comparebytime);

  for(const auto & c : snclusters) savecluster(evt, c);

  unsigned int unpairedbighits = 0, unpairedsmallhits = 0;
  for(unsigned int i = 0; i < nmhits; i++){
    if(mhits[i].paired) continue;
    (mhits[i].pe > bighit_threshold? unpairedbighits: unpairedsmallhits)++;
  }

  if(supernovalike){
    THplusequals(lh_supernovalike, timebin(evt), hitclusters, livetime);
  }else{
    const unsigned int big = hitclusters + unpairedbighits;
    const unsigned int all = big + unpairedsmallhits;
    THplusequals(lh_unsliced_big_hits, timebin(evt), big, livetime);
    THplusequals(lh_unsliced_hits,     timebin(evt), all, livetime);
  }
}

// Passes back the {ra, dec} of a track direction given an event and that
// direction.
static void track_ra_dec(double & ra, double & dec,
                         const art::Event & evt, const TVector3 & dir)
{
  art::ServiceHandle<locator::CelestialLocator> celloc;

  long int unixtime = 0;
  if(needbgevent_unix_double_time == 0){
    unixtime = int(art_time_to_unix_double(evt.time().value()));
  }
  else{
    // If we're measuring background, set the time to what it would be
    // for the real event at the same point in the search window.
    const double evt_time = art_time_to_unix_double(evt.time().value());
    const double offset = evt_time - gwevent_unix_double_time;
    const double timetouse = needbgevent_unix_double_time + offset;
    unixtime = int(timetouse);
  }

  celloc->GetTrackRaDec(dir, unixtime, ra, dec);
}

// Returns whether the given ra and dec are inside the 90% region of map
// number 'mapi'. If 'allow_backwards' accept tracks that point in the
// reverse direction, since the orientation of the track is arbitrary.
// Otherwise, take the track orientation literally.
static bool points(const double ra, const double dec, const int mapi,
                   const bool allow_backwards)
{
  // not sure if normalizing is necessary, but can't hurt
  pointing forward ( dec+M_PI_2, ra);
  forward.normalize();
  pointing backward(-dec+M_PI_2, ra+M_PI);
  backward.normalize();

  if(allow_backwards)
    return skymap_crit_val[mapi] < std::max(
    healpix_skymap[mapi]->interpolated_value(forward),
    healpix_skymap[mapi]->interpolated_value(backward));
  else
    return skymap_crit_val[mapi] <
    healpix_skymap[mapi]->interpolated_value(forward);
}

static void count_upmu(const art::Event & evt)
{
  art::Handle< std::vector<rb::Track> > upmu;

  evt.getByLabel("upmuanalysis", upmu);
  if(upmu.failedToGet()){
    fprintf(stderr, "No UpMu product to read\n");
    return;
  }

  const double livetime = rawlivetime(evt);

  unsigned int upmucount = 0;
  for(unsigned int i = 0; i < upmu->size(); i++)
    if(uniquedata_tns(livetime, (*upmu)[i].MeanTNS()))
      upmucount++;

  printf("Up-mu tracks: %d\n", upmucount);

  THplusequals(lh_upmu_tracks, timebin(evt), upmucount, livetime);

  int npoint[npointres] = { 0 };
  for(unsigned int i = 0; i < upmu->size(); i++){
    if(!uniquedata_tns(livetime, (*upmu)[i].MeanTNS())) continue;

    if(healpix_skymap[0] == NULL) continue;

    // upmu tracks can still point the wrong way even though they have
    // been selected to be upwards-going, because they are just copied
    // from the input track array.
    const TVector3 trackdir = (*upmu)[i].Dir();
    TVector3 correctdir = trackdir;
    if(trackdir.Z() < 0){
      correctdir.SetX(-trackdir.X());
      correctdir.SetY(-trackdir.Y());
      correctdir.SetZ(-trackdir.Z());
    }

    double ra, dec;
    track_ra_dec(ra, dec, evt, correctdir);
    for(unsigned int q = 0; q < npointres; q++){
      printf("Up-mu track pointing at region %d with ra,dec = %f %f\n",
             q, ra, dec);
      npoint[q] += points(ra, dec, q, false);
    }
  }

  for(unsigned int q = 0; q < npointres; q++){
    if(npoint[q])printf("%d up-mu tracks point at region %d\n", npoint[q], q);
    THplusequals(lh_upmu_tracks_point[q], timebin(evt),
                 npoint[q], rawlivetime(evt));
  }
}

// Counts tracks with various cuts and contained slices and adds the
// results to the output histograms
static void
count_tracks_containedslices(const art::Event & evt,
                             const std::vector<mslice> & sliceinfo)
{
  art::Handle< std::vector<rb::Cluster> > slice;
  evt.getByLabel("slicer", slice);

  art::Handle< std::vector<rb::Track> > tracks;
  evt.getByLabel("breakpoint", tracks);
  if(tracks.failedToGet()){
    fprintf(stderr, "Unexpected lack of breakpoint product!\n");
    return;
  }

  const int timbin = timebin(evt);
  const double livetime = rawlivetime(evt);

  // Count tracks with either a contained start or end, but not if
  // they are very steep in x or y, since those probably aren't really
  // contained and may not even be complete. Also don't count more than
  // one per slice, nor count slices with uncontained tracks. This
  // protects against counting a brem as a contained track.

  // Find out what slice each track is in and make a list of slices with
  // tracks that aren't fully contained
  std::vector<int> trk_slices(tracks->size(), -1);
  std::set<int> slices_with_tracks,
                slices_with_uc_tracks,
                slices_with_huc_tracks,
                contained_slices,
                contained_shower_slices;
  for(unsigned int i = 0; i < slice->size(); i++){
    if(!(*slice)[i].NCell(geo::kX) || !(*slice)[i].NCell(geo::kY))
      continue;

    if(!uniquedata_tns(livetime, (*slice)[i].MeanTNS())) continue;

    if(!contained((*slice)[i].MinXYZ()) ||
       !contained((*slice)[i].MaxXYZ())) continue;

    contained_slices.insert(i);

    if((*slice)[i].NCell() < 10) continue;

    const int planesx = ((*slice)[i].ExtentPlane(geo::kX)+1)/2;
    const int planesy = ((*slice)[i].ExtentPlane(geo::kY)+1)/2;

    if(planesx < 5 || planesy < 5) continue;

    const int cellsx = (*slice)[i].ExtentCell(geo::kX);
    const int cellsy = (*slice)[i].ExtentCell(geo::kY);

    // Must be somewhat extended in z (not a straight-down track), and
    // whatever box is defined by the cell and plane extent must neither
    // be too sparse (random hits plus a straight-down track) nor too full
    // (FEB flashing).  A block of flashers plus a random hit is not
    // really excludable with this cut.
    //
    // These numbers are reasonable, but not very rigorous.
    const double min_boxupancy = 0.02,
                 max_boxupancy = 0.10;
    const double boxupancy_x = (*slice)[i].NXCell()/double(cellsx*planesx),
                 boxupancy_y = (*slice)[i].NYCell()/double(cellsy*planesy);

    if(boxupancy_x < max_boxupancy && boxupancy_y < max_boxupancy &&
       boxupancy_x > min_boxupancy && boxupancy_y > min_boxupancy)
      contained_shower_slices.insert(i);
  }

  // After next loop, these are true if the corresponding track in 'tracks'
  // points back to the sky region from LIGO/Virgo
  std::vector<bool> track_point[npointres];

  for(unsigned int i = 0; i < tracks->size(); i++){
    if(!uniquedata_tns(livetime, (*tracks)[i].MeanTNS())){
      for(unsigned int q = 0; q < npointres; q++)
        track_point[q].push_back(false);
      continue;
    }

    trk_slices[i] =
      which_slice_is_this_track_in((*tracks)[i], slice, sliceinfo);

    if(un_contained_track((*tracks)[i]))
      slices_with_uc_tracks.insert(trk_slices[i]);
    if(half_uncontained_track((*tracks)[i]))
      slices_with_huc_tracks.insert(trk_slices[i]);

    bool track_points_to_event[npointres] = { false };

    if(healpix_skymap[0] != NULL){
      double ra, dec;
      track_ra_dec(ra, dec, evt, (*tracks)[i].Dir());
      for(unsigned int q = 0; q < npointres; q++)
        track_points_to_event[q] = points(ra, dec, q, true);
    }
    for(unsigned int q = 0; q < npointres; q++)
      track_point[q].push_back(track_points_to_event[q]);
  }

  // Find any tracks that are half-contained/fully-contained
  // and don't share a slice with any tracks
  // that are, like, totally uncontained, man.
  std::set<int> slices_with_hc_tracks, slices_with_fc_tracks;

  // Same, but the tracks must point towards the LIGO/Virgo event
  std::set<int> slices_with_tracks_point[npointres],
                slices_with_hc_tracks_point[npointres],
                slices_with_fc_tracks_point[npointres];

  for(unsigned int i = 0; i < tracks->size(); i++){
    if(!uniquedata_tns(livetime, (*tracks)[i].MeanTNS())) continue;

    // Exclude 2D tracks
    if((*tracks)[i].Stop().X() == 0 || (*tracks)[i].Stop().Y() == 0) continue;

    slices_with_tracks.insert(trk_slices[i]);
    for(unsigned int q = 0; q < npointres; q++)
      if(track_point[q][i])
        slices_with_tracks_point[q].insert(trk_slices[i]);

    // To be called a slice with half-contained tracks, it must not have any
    // tracks that both enter and exit
    if(!slices_with_uc_tracks.count(trk_slices[i]) &&
       half_contained_track((*tracks)[i])){
      slices_with_hc_tracks.insert(trk_slices[i]);
      for(unsigned int q = 0; q < npointres; q++)
        if(track_point[q][i])
          slices_with_hc_tracks_point[q].insert(trk_slices[i]);
    }

    // To be called a slice with fully contained tracks, it must not have any
    // track that either enters or exits, and the slice itself must be
    // contained.  So, in other words, no hits around the edges, and no tracks
    // reconstructed to be near the edges even in the absence of hits.
    if(!slices_with_huc_tracks.count(trk_slices[i]) &&
       fully_contained_track((*tracks)[i]) &&
       contained_slices.count(trk_slices[i])){
      slices_with_fc_tracks.insert(trk_slices[i]);
      for(unsigned int q = 0; q < npointres; q++)
        if(track_point[q][i])
          slices_with_fc_tracks_point[q].insert(trk_slices[i]);
    }
  }

#ifdef LOUD
  printf("Slices with tracks: %2lu (",
         slices_with_tracks.size());
  for(unsigned int q = 0; q < npointres; q++) printf("%lu%s",
    slices_with_tracks_point[q].size(), q==npointres-1?" pointing)\n":", ");

  printf("Slices with half-contained tracks: %2lu (",
         slices_with_hc_tracks.size());
  for(unsigned int q = 0; q < npointres; q++) printf("%lu%s",
    slices_with_hc_tracks_point[q].size(),q==npointres-1?" pointing)\n":", ");

  printf("Slices with fully-contained tracks: %2lu (",
         slices_with_fc_tracks.size());
  for(unsigned int q = 0; q < npointres; q++) printf("%lu%s",
    slices_with_fc_tracks_point[q].size(),q==npointres-1?" pointing)\n":", ");


  for(std::set<int>::iterator i = slices_with_fc_tracks.begin();
      i != slices_with_fc_tracks.end(); i++)
    printf("  slice %3d\n", *i);

  printf("Contained GeV physics-like slices: %lu\n",
         contained_shower_slices.size());
  for(std::set<int>::iterator i = contained_shower_slices.begin();
      i != contained_shower_slices.end(); i++)
    printf("  slice %3d\n", *i);
#endif

  THplusequals(
    lh_tracks, timbin, slices_with_tracks.size(), livetime);
  THplusequals(
    lh_halfcontained_tracks, timbin, slices_with_hc_tracks.size(), livetime);
  THplusequals(
    lh_fullycontained_tracks, timbin, slices_with_fc_tracks.size(), livetime);
  THplusequals(
    lh_contained_slices, timbin, contained_shower_slices.size(), livetime);

  for(unsigned int q = 0; q < npointres; q++){
    THplusequals(lh_tracks_point[q], timbin,
                 slices_with_tracks_point[q].size(), livetime);
    THplusequals(lh_halfcontained_tracks_point[q], timbin,
                 slices_with_hc_tracks_point[q].size(), livetime);
    THplusequals(lh_fullycontained_tracks_point[q], timbin,
                 slices_with_fc_tracks_point[q].size(), livetime);
  }
  fflush(stdout);
}

static void count_all_mev(art::Event & evt,
                          const std::vector<mslice> & sliceinfo,
                          const std::vector<mtrack> & trackinfo)
{
  count_mev(evt, true , sliceinfo, trackinfo); // supernova-like
  count_mev(evt, false, sliceinfo, trackinfo); // unmodeled low energy
}

static void count_livetime(const art::Event & evt)
{
  const double livetime = rawlivetime(evt);
  const double veryrawlivetime = rawlivetime(evt, true);
  printf("Livetime %g seconds, %g used\n", veryrawlivetime, livetime);
  THplusequals(lh_blind, timebin(evt, true), 0, livetime);
}

static bool more_than_one_physics_slice(art::Event & evt)
{
  art::Handle< std::vector<rb::Cluster> > slice;
  evt.getByLabel("slicer", slice);
  // Event probably filtered out in reco file
  if(slice.failedToGet()) return false;
  return slice->size() > 2;
}

void ligoanalysis::produce(art::Event & evt)
{
  {
    static unsigned int n = 0;
    // Start at 1 because the first event takes forever as everything
    // is initialized.
    if(n == 1){
      // art provides no way of knowing how much events we will process,
      // and its own progress indicator is of limited use.
      #if 1
      printf("For progress indicator, assuming a 2000 event long readout\n");
      initprogressindicator(2000-1, 3);
      #else
      printf("For progress indicator, assuming a 1e6 event MC\n");
      initprogressindicator(1000000-1, 6);
      #endif
    }
    else if(n > 1){
      progressindicator(n - 1);
    }
    n++;
  }

  // Must be called on every event to prevent SIGPIPE loops.
  signal(SIGPIPE, SIG_DFL);

  {
    art::Handle< std::vector<rawdata::RawTrigger> > rawtrigger;
    getrawtrigger(rawtrigger, evt);
    if(!goodtriggertype(trigger(rawtrigger))) return;

    art::Handle< std::vector<rawdata::FlatDAQData> > flatdaq;
    getflatdaq(flatdaq, evt);
    if(!flatdaq.failedToGet() && !is_complete_event(flatdaq)){
      printf("WARNING: Incomplete event, but assuming can trust livetime\n");
    }

    if(longtriggertype(trigger(rawtrigger))){
      int64_t event_length_tdc, delta_tdc;
      delta_and_length(event_length_tdc, delta_tdc, flatdaq, rawtrigger);

      // Check for empty 50us events in long readouts
      if(event_length_tdc/TDC_PER_US == 50){
        art::Handle< std::vector< rawdata::RawDigit > > rd;
        getrawdigits(rd, evt);
        if(rd->empty()){
          printf("Rejecting LIGO_TRIGGER 50us event with zero hits\n");
          return;
        }
      }

      // Check for truly incomplete readouts (not to be confused with events
      // marked incomplete).  Could use these, probably, but would need to
      // modify code to handle overlaps to determine when the overlaps are in
      // these cases, and the complexity doesn't seem worth it.
      const int len = event_length_tdc/TDC_PER_US;
      if(len != 5050 && len != 5000){
        printf("Rejecting LIGO_TRIGGER length %dus != 5050 or 5000\n", len);
        return;
      }
    }
  }

  art::ServiceHandle<geo::Geometry> geo;
  gDet = geo->DetId();

  bighit_threshold = gDet == caf::kNEARDET?
    bighit_threshold_nd : bighit_threshold_fd;

  art::Handle< std::vector<rb::Cluster> > slice;
  if(fAnalysisClass != Blind){
    // Reject any ND trigger with multiple physics slices. This is a
    // crude way of getting rid of NuMI events, previously used by the
    // seasonal multi-mu analysis because RemoveBeamSpills doesn't
    // really work.
    if(fCutNDmultislices &&
       (fAnalysisClass == NDactivity || fAnalysisClass == MinBiasND) &&
       more_than_one_physics_slice(evt)){
      printf("Cut event with multiple slices\n");
      return;
    }

    if(fAnalysisClass != DDenergy){
      evt.getByLabel("slicer", slice);
      // Event probably filtered out in reco file
      if(slice.failedToGet()) return;

      if(slice->empty()){
        fprintf(stderr, "Unexpected event with zero slices!\n");
        return;
      }
    }
  }

  const std::vector<mtrack> trackinfo = 
  (fAnalysisClass == SNonlyFD || fAnalysisClass == MichelFD)?
    save_mtracks(evt): std::vector<mtrack>();

  // Make list of all times and locations of "physics" slices. For the
  // MeV search, we will exclude other hits near them in time to drop
  // Michels, FEB flashers & neutrons. For the track search, we will
  // merge slices that overlap in time to make slices-with-track counts
  // more Poissonian.
  const std::vector<mslice> sliceinfo =
    (fAnalysisClass == Blind || fAnalysisClass == DDenergy)?
      std::vector<mslice>(): make_sliceinfo_list(evt, slice);

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

  switch(fAnalysisClass){
    case NDactivity:
      count_tracks_containedslices(evt, sliceinfo);
      break;
    case DDenergy:
      count_triggers(evt);
      count_ddenergy(evt);
      break;
    case MinBiasFD:
      count_tracks_containedslices(evt, sliceinfo);
      count_upmu(evt);
      count_all_mev(evt, sliceinfo_wprev, trackinfo_wprev);
      break;
    case SNonlyFD:
    case MichelFD:
      count_triggers(evt);
      count_mev(evt, true /* SN-like */, sliceinfo_wprev, trackinfo_wprev);
      break;
    case MinBiasND:
      count_tracks_containedslices(evt, sliceinfo);
      count_all_mev(evt, sliceinfo_wprev, trackinfo_wprev);
      break;
    case SNonlyND:
      count_triggers(evt);
      count_mev(evt, true /* SN-like */, sliceinfo_wprev, trackinfo_wprev);
      break;
    case Blind:
      count_livetime(evt);
      break;
    default:
      printf("No case for type %d\n", fAnalysisClass);
  }

  prev_trackinfo = trackinfo;
  prev_sliceinfo = sliceinfo;

  // Translate times for the next event, 5ms later. This ONLY makes
  // sense for long readouts, and is wrong if a trigger is dropped.
  for(unsigned int i = 0; i < prev_trackinfo.size(); i++)
    prev_trackinfo[i].tns -= 5e6; // 5 ms in ns
  for(unsigned int i = 0; i < prev_sliceinfo.size(); i++){
    prev_sliceinfo[i].mintns -= 5e6;
    prev_sliceinfo[i].maxtns -= 5e6;
  }
}

DEFINE_ART_MODULE(ligoanalysis)
