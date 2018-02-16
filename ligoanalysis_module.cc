////////////////////////////////////////////////////////////////////////
/// \brief   This module is named ligoanalysis and looks at GW coincidences
/// \author  M. Strait
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

#include "TH1.h"

#include <string>
using std::string;
#include <vector>
#include <set>
#include <algorithm>

#include <signal.h>

#include "func/timeutil.h"

static const int TDC_PER_US = 64;
static const int US_PER_MICROSLICE = 50; // I hope this is always true

// Set from the FCL parameters
static double gwevent_unix_double_time = 0;
static long long window_size_s = 1000;

static int gDet = caf::kUNKNOWN;

// Types of analysis, dependent on which trigger we're looking at
enum analysis_class_t { NDactivity, LiveTime, UpMu, DDenergy,
                        JustTrigger, NDMeV, MAX_ANALYSIS_CLASS };

class ligoanalysis : public art::EDProducer {
  public:
  explicit ligoanalysis(fhicl::ParameterSet const& pset);
  virtual ~ligoanalysis() { }; // compiles, but does not run, without this
  void produce(art::Event& evt);

  /// \brief User-supplied type of trigger being examined.
  ///
  /// Which histograms we make depends on this.
  analysis_class_t fAnalysisClass;

  /// \brief The user-supplied time of the gravitational wave burst, or
  /// whatever time you want to center the search on.
  ///
  /// Expressed in RFC-3339 format, always in UTC, always with a Z at
  /// the end (to emphasize that it is UTC).
  string fGWEventTime;

  /// \brief The user-supplied length of the search window in seconds.
  ///
  /// We search for half this length on either side of the given time.
  float fWindowSize;
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
  string name;

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

  ligohist(const string & name_, const bool dolive_)
  {
    name = name_;
    dolive = dolive_;
  }
};

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

/**********************************************************************/

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

// Count of triggers, with no examination of the data within
static ligohist lh_rawtrigger("rawtrigger", false);

// Number of hits, with no filtering of any sort
static ligohist lh_rawhits("rawhits", true);

// Number of hits that are not in any Slicer4D slice, i.e. they are in the
// Slicer4D noise slice.
static ligohist lh_unslice4ddhits("unslice4ddhits", true);

// My experience with my neutron capture analysis in the ND is that
// things become well-behaved at about this level.
const double bighit_threshold = 35; // PE

// At this point (in the ND) the hit is almost definitely physics.  It's about
// 2/3 of a MIP, and is past the end of the steeply falling noise curve.
const double biggerhit_threshold = 100; // PE

// Number of unsliced hits that are over the above PE threshold
static ligohist lh_unsliced_big_hits("unslicedbighits", true);

// Number of unsliced hits that are over some PE threshold and are paired
// with a hit in an adjacent plane.  (Not an adjacent cell, because then
// we select lots of pairs from noisy modules.  We could build a noise map
// to fix that, but it would require two passes through the data.)
static ligohist lh_unsliced_hit_pairs("unslicedhitpairs", true);

// Raw number of tracks
static ligohist lh_tracks("tracks", true);

// Number of tracks with at least one contained endpoint, with some additional
// sanity checks.
static ligohist lh_halfcontained_tracks("halfcontained_tracks", true);

// Number of tracks with at two contained endpoints, with some additional
// sanity checks.
static ligohist lh_fullycontained_tracks("fullycontained_tracks", true);

// Count of slices with nothing around the edges, regardless of what sorts of
// objects are inside.
static ligohist lh_contained_slices("contained_slices", true);

// Number of tracks that pass the Upmu analysis
static ligohist lh_upmu_tracks("upmu_tracks", false);

// Number of triggers above two cuts for DDEnergy, the first pair
// in raw ADC, the second ADC per unit time.
static ligohist lh_ddenergy_locut("energy_low_cut", false);
static ligohist lh_ddenergy_hicut("energy_high_cut", false);
static ligohist lh_ddenergy_lopertime("energy_low_cut_pertime", false);
static ligohist lh_ddenergy_hipertime("energy_high_cut_pertime", false);

static void init_mev_hists()
{
  init_lh(lh_rawhits);
  init_lh(lh_unslice4ddhits);
  init_lh(lh_unsliced_big_hits);
  init_lh(lh_unsliced_hit_pairs);
}

static void init_track_and_contained_hists()
{
  init_lh(lh_tracks);
  init_lh(lh_halfcontained_tracks);
  init_lh(lh_fullycontained_tracks);
  init_lh(lh_contained_slices);
}

ligoanalysis::ligoanalysis(fhicl::ParameterSet const& pset) : EDProducer(),
  fGWEventTime(pset.get<string>("GWEventTime")),
  fWindowSize(pset.get<unsigned long long>("WindowSize"))
{
  const string analysis_class_string(pset.get<string>("AnalysisClass"));

  if     (analysis_class_string == "NDactivity") fAnalysisClass = NDactivity;
  else if(analysis_class_string == "LiveTime")   fAnalysisClass = LiveTime;
  else if(analysis_class_string == "UpMu")       fAnalysisClass = UpMu;
  else if(analysis_class_string == "DDenergy")   fAnalysisClass = DDenergy;
  else if(analysis_class_string == "JustTrigger")fAnalysisClass = JustTrigger;
  else if(analysis_class_string == "NDMeV")      fAnalysisClass = NDMeV;
  else{
    fprintf(stderr, "Unknown AnalysisClass \"%s\" in job fcl. See list "
            "in ligoanalysis.fcl.\n", analysis_class_string.c_str());
    exit(1);
  }

  gwevent_unix_double_time = rfc3339_to_unix_double(fGWEventTime);
  window_size_s = fWindowSize;

  switch(fAnalysisClass){
    case NDactivity:
      init_lh(lh_rawtrigger);
      init_track_and_contained_hists();
      break;
    case JustTrigger:
      init_lh(lh_rawtrigger);
      break;
    case UpMu:
      init_lh(lh_upmu_tracks);
      break;
    case LiveTime:
      init_mev_hists();
      init_track_and_contained_hists();
      break;
    case DDenergy:
      init_lh(lh_rawtrigger);
      init_lh(lh_ddenergy_locut);
      init_lh(lh_ddenergy_hicut);
      init_lh(lh_ddenergy_lopertime);
      init_lh(lh_ddenergy_hipertime);
      break;
    case NDMeV:
      init_mev_hists();
      break;
    default:
      printf("No case for type %d\n", fAnalysisClass);
  }
}

/**********************************************************************/
/*                          The meat follows                          */
/**********************************************************************/

// Return the livetime in this event in seconds, as it is relevant for
// raw hits (same as for anything else? Maybe not if we don't trust
// tracks close to the time-edges of events).
static double rawlivetime(const art::Event & evt)
{
  art::Handle< std::vector<rawdata::FlatDAQData> > flatdaq;
  evt.getByLabel("daq", flatdaq);

  art::Handle< std::vector<rawdata::RawTrigger> > rawtrigger;
  evt.getByLabel("daq", rawtrigger);

  if(rawtrigger->empty()) return -1;

  int64_t event_length_tdc = 0, delta_tdc = 0;

  if(!delta_and_length(event_length_tdc, delta_tdc, flatdaq, rawtrigger))
    return -1;

  return event_length_tdc / TDC_PER_US * 1e-6;
}

// Add 'sig' to the signal and 'live' to the livetime in bin number 'bin'.
static void THplusequals(ligohist & lh, const int bin, const double sig,
                         const double live)
{
  // Use SetBinContent instead of Fill(x, weight) to avoid having to look up
  // the bin number twice.
  lh.sig ->SetBinContent(bin, lh.sig ->GetBinContent(bin) + sig);
  if(lh.dolive)
    lh.live->SetBinContent(bin, lh.live->GetBinContent(bin) + live);
}

// Return the bin number for this event, i.e. the number of seconds from the
// beginning of the window, plus 1.
static int timebin(const art::Event & evt)
{
  const double evt_time = art_time_to_unix_double(evt.time().value());
  return floor(evt_time - gwevent_unix_double_time) + window_size_s/2
         + 1; // stupid ROOT 1-based numbering!
}

static bool contained(const TVector3 & v)
{
  if(gDet == caf::kNEARDET)
    return fabs(v.X()) < 150 && fabs(v.Y()) < 150 && v.Z() > 40 && v.Z() < 1225;
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

static const int min_plane_extent = 10;

static const double tan_track_cut = atan(15 * M_PI/180);

// Check that the reconstruction (probably BreakPointFitter) doesn't think
// either end of the track points nearly along a plane AND check that if we
// just draw a line from one end of the track to the other, that that doesn't
// either.  This second part is to make sure that some tiny kink at the start
// or end can't fool us into thinking it is well-contained.
static bool good_track_direction(const rb::Track & t)
{
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

// return the index in the slice array that the given track is in
static int which_slice_is_this_track_in(
  const rb::Track & t,
  const art::Handle< std::vector<rb::Cluster> > & slice)
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
      if(*ahit == *shit) return i;
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
  THplusequals(lh_rawtrigger, timebin(evt), 1, 1);
}

static void count_ddenergy(const art::Event & evt)
{
  art::Handle< std::vector<rawdata::RawDigit> > rawhits;
  evt.getByLabel("daq", rawhits);

  int64_t sumadc = 0;

  for(unsigned int i = 0; i < rawhits->size(); i++)
    sumadc += (*rawhits)[i].ADC();

  const double rawtime = rawlivetime(evt);

  printf("ADC: %ld %12.0f\n", sumadc, sumadc/rawtime);

  if(sumadc/rawtime > 5e10)
    THplusequals(lh_ddenergy_lopertime, timebin(evt), 1, rawtime);
  if(sumadc/rawtime > 5e11)
    THplusequals(lh_ddenergy_hipertime, timebin(evt), 1, rawtime);

  if(sumadc >  5000000)
    THplusequals(lh_ddenergy_locut, timebin(evt), 1, rawtime);
  if(sumadc > 50000000)
    THplusequals(lh_ddenergy_hicut, timebin(evt), 1, rawtime);
}

// Count the number of raw hits in the event and fill the appropriate histograms
static void count_hits(const art::Event & evt)
{
  art::Handle< std::vector<rawdata::RawDigit> > rawhits;
  evt.getByLabel("daq", rawhits);

  printf("Hits in this event: %lu\n", rawhits->size());

  THplusequals(lh_rawhits, timebin(evt), rawhits->size(), rawlivetime(evt));
}

struct mhit{ float tns; uint16_t plane; bool bigger; };

static bool mhit_by_time(const mhit & a, const mhit & b)
{
  return a.tns < b.tns;
}

static void count_unsliced_hit_pairs(const art::Event & evt)
{
  art::Handle< std::vector<rb::Cluster> > slice;
  evt.getByLabel("slicer", slice);

  if(slice->empty()){
    printf("Unexpected event with zero slices!\n");
    return;
  }

  std::vector<mhit> mhits;

  // Make a list of all the times of all the hits in "physics" slices.
  // We're going to exclude other hits near them in time to drop Michels
  // and maybe also neutrons.
  std::set<float> slicetimes;
  for(unsigned int i = 1; i < slice->size(); i++)
    for(unsigned int j = 0; j < (*slice)[i].NCell(); j++)
       slicetimes.insert((*slice)[i].Cell(j)->TNS());

  for(unsigned int i = 0; i < (*slice)[0].NCell(); i++){
    if((*slice)[0].Cell(i)->PE() <= bighit_threshold) continue;

    const int cell  = (*slice)[0].Cell(i)->Cell();
    const int plane = (*slice)[0].Cell(i)->Plane();

    if(gDet == caf::kNEARDET){
      if(plane == 0) continue;

      // Exclude last regular plane and the whole muon catcher. Can't
      // reasonably have a supernova-type event hit planes on each
      // side of a steel plane, and while they might hit adjacent
      // scintillator planes in the muon catcher, I don't want to deal
      // with all the additional complications there.
      if(plane > 190) continue;

      // Drop outermost three cells. This excludes most muon tracks that
      // just barely enter the detector, but don't get reconstructed.
      if(cell <= 2 || cell >= 95) continue;

      // At top of detector, stricter cut based on observed backgrounds
      if((*slice)[0].Cell(i)->View() == geo::kY && cell >= 80) continue;
    }

    const float tns = (*slice)[0].Cell(i)->TNS();

    if(!slicetimes.empty()){
      std::set<float>::iterator nextslicetime = slicetimes.upper_bound(tns);

      if(nextslicetime != slicetimes.begin()){
        nextslicetime--;
        const float time_since_slc = tns - (*nextslicetime);
        const float time_since_slc_cut = 50e3; // 1 neutron, many muon lifetimes

        // Just drop entire detector in this case. That's fine for the
        // ND, but probably not for the FD, so come back here for that.
        if(time_since_slc < time_since_slc_cut) continue;
      }
    }

    mhit h;
    h.tns   = tns;
    h.plane = plane;
    h.bigger = (*slice)[0].Cell(i)->PE() > biggerhit_threshold;

    mhits.push_back(h);
  }

  std::sort(mhits.begin(), mhits.end(), mhit_by_time);

  // Same cut as I use for clustering Michels and neutron capture hits
  const float timewindow = 500.0; // ns

  unsigned int hitpairs = 0;

  for(unsigned int i = 0; i < mhits.size(); i++){
    for(unsigned int j = i+1; j < mhits.size(); j++){
      // Require at least one to be pretty definite physics
      if(!mhits[j].bigger && ! mhits[i].bigger) break;

      if(mhits[j].tns - mhits[i].tns > timewindow) break;
      if(abs(mhits[j].plane - mhits[i].plane) != 1) continue;

      hitpairs++;

      // Do not allow these hits to be used again.  (The first one is automatic
      // since we won't ever look at it again.) Is this hacky?  Yes.
      mhits.erase(mhits.begin()+j);
      break;
    }
  }

  printf("Big hit pairs in the noise slice: %u\n", hitpairs);

  THplusequals(lh_unsliced_hit_pairs, timebin(evt), hitpairs,
               rawlivetime(evt));
}

static void count_unsliced_big_hits(const art::Event & evt)
{
  art::Handle< std::vector<rb::Cluster> > slice;
  evt.getByLabel("slicer", slice);

  if(slice->empty()){
    printf("Unexpected event with zero slices!\n");
    return;
  }

  unsigned int bighits = 0;

  for(unsigned int i = 0; i < (*slice)[0].NCell(); i++)
    bighits += (*slice)[0].Cell(i)->PE() > bighit_threshold;

  printf("Big hits in this event in the noise slice: %u\n", bighits);

  THplusequals(lh_unsliced_big_hits, timebin(evt), bighits,
               rawlivetime(evt));
}

// Counts hits in the noise slices and adds the results to the output histograms
static void count_unslice4dd_hits(const art::Event & evt)
{
  art::Handle< std::vector<rb::Cluster> > slice;
  evt.getByLabel("slicer", slice);

  if(slice->empty()){
    printf("Unexpected event with zero slices!\n");
    return;
  }

  printf("Hits in this event in the Slicer4D noise slice: %d\n",
         (*slice)[0].NCell());

  THplusequals(lh_unslice4ddhits, timebin(evt), (*slice)[0].NCell(),
              rawlivetime(evt));
}

static void count_upmu(const art::Event & evt)
{
  art::Handle< std::vector<rb::Track> > upmu;
  if(!evt.getByLabel("upmuanalysis", upmu)){
    printf("No UpMu product to read\n");
    _exit(1);
  }

  THplusequals(lh_upmu_tracks, timebin(evt),
               upmu->size(), rawlivetime(evt));
}

// Counts tracks with various cuts and contained slices and adds the
// results to the output histograms
static void count_tracks_containedslices(const art::Event & evt)
{
  art::Handle< std::vector<rb::Cluster> > slice;
  evt.getByLabel("slicer", slice);

  art::Handle< std::vector<rb::Track> > tracks;
  evt.getByLabel("breakpoint", tracks);

  printf("Tracks in this event: %lu\n", tracks->size());
  THplusequals(lh_tracks, timebin(evt), tracks->size(), rawlivetime(evt));

  // Count tracks with either a contained start or end, but not if they are
  // very steep in x or y, since those probably aren't really contained and may
  // not even be complete. Also don't count more than one per slice, nor count
  // slices with uncontained tracks.  This protects against counting a brem as
  // a contained track.

  // Find out what slice each track is in and make a list of slices with tracks
  // that aren't fully contained
  std::vector<int> trk_slices(tracks->size(), -1);
  std::set<int> slices_with_uc_tracks, slices_with_huc_tracks;
  std::set<int> contained_slices;
  std::set<int> contained_shower_slices;
  for(unsigned int i = 0; i < slice->size(); i++){
    if(!(*slice)[i].NCell(geo::kX) || !(*slice)[i].NCell(geo::kY))
      continue;

    if(!contained((*slice)[i].MinXYZ()) ||
       !contained((*slice)[i].MaxXYZ())) continue;

    contained_slices.insert(i);

    if((*slice)[i].NCell() < 10) continue;

    const int planesx = (*slice)[i].ExtentPlane(geo::kX);
    const int planesy = (*slice)[i].ExtentPlane(geo::kY);
    const int cellsx = (*slice)[i].ExtentCell(geo::kX);
    const int cellsy = (*slice)[i].ExtentCell(geo::kY);

    // Must be somewhat extended in z (not a straight-down track), and
    // whatever box is defined by the cell and plane extent must neither
    // be too sparse (random hits plus a straight-down track) nor too full
    // (FEB flashing).  A block of flashers plus a random hit is not
    // really excludable with this cut.  Note that the maximum is ~0.5
    // because the count of planes is in both views.
    const double min_boxupancy = 0.1,
                 // max_boxupancy = 0.45; // seems to work
                 max_boxupancy = 1; // shouldn't work

    if(planesx >= 9 && planesy >= 9 && // always odd
       (*slice)[i].NXCell()/double(cellsx*planesx) < max_boxupancy &&
       (*slice)[i].NYCell()/double(cellsy*planesy) < max_boxupancy &&
       (*slice)[i].NXCell()/double(cellsx*planesx) > min_boxupancy &&
       (*slice)[i].NYCell()/double(cellsy*planesy) > min_boxupancy)
      contained_shower_slices.insert(i);
  }

  for(unsigned int i = 0; i < tracks->size(); i++){
    trk_slices[i] = which_slice_is_this_track_in((*tracks)[i], slice);
    if(un_contained_track((*tracks)[i]))
      slices_with_uc_tracks.insert(trk_slices[i]);
    if(half_uncontained_track((*tracks)[i]))
      slices_with_huc_tracks.insert(trk_slices[i]);
  }

  // Find any tracks that are half-contained/fully-contained and do not share a
  // slice with any tracks that are, like, totally uncontained, man
  std::set<int> slices_with_hc_tracks, slices_with_fc_tracks;
  for(unsigned int i = 0; i < tracks->size(); i++){
    // Exclude 2D tracks
    if((*tracks)[i].Stop().X() == 0 || (*tracks)[i].Stop().Y() == 0) continue;

    // To be called a slice with half-contained tracks, it must not have any
    // tracks that both enter and exit
    if(!slices_with_uc_tracks.count(trk_slices[i]) &&
       half_contained_track((*tracks)[i]))
      slices_with_hc_tracks.insert(trk_slices[i]);

    // To be called a slice with fully contained tracks, it must not have any
    // track that either enters or exits, and the slice itself must be
    // contained.  So, in other words, no hits around the edges, and no tracks
    // reconstructed to be near the edges even in the absence of hits.
    if(!slices_with_huc_tracks.count(trk_slices[i]) &&
       fully_contained_track((*tracks)[i]) &&
       contained_slices.count(trk_slices[i])){
      slices_with_fc_tracks.insert(trk_slices[i]);
      printf("start cosx %6.3f cosy %6.3f\n",
             (*tracks)[i].Dir().X(), (*tracks)[i].Dir().Y());
      printf("stop  cosx %6.3f cosy %6.3f\n",
             (*tracks)[i].StopDir().X(), (*tracks)[i].StopDir().Y());
    }
  }

  printf("Slices with half-contained tracks: %lu\n",
         slices_with_hc_tracks.size());
  printf("Slices with fully-contained tracks: %lu\n",
         slices_with_fc_tracks.size());
  for(std::set<int>::iterator i = slices_with_fc_tracks.begin();
      i != slices_with_fc_tracks.end(); i++)
    printf("  %3d\n", *i);
  printf("Contained GeV physics-like slices: %lu\n",
         contained_shower_slices.size());
  for(std::set<int>::iterator i = contained_shower_slices.begin();
      i != contained_shower_slices.end(); i++)
    printf("  %3d\n", *i);

  THplusequals(lh_halfcontained_tracks,  timebin(evt),
               slices_with_hc_tracks.size(), rawlivetime(evt));
  THplusequals(lh_fullycontained_tracks, timebin(evt),
               slices_with_fc_tracks.size(), rawlivetime(evt));
  THplusequals(lh_contained_slices, timebin(evt),
               contained_shower_slices.size(), rawlivetime(evt));
}

static void count_mev(art::Event & evt)
{
  count_hits(evt);
  count_unslice4dd_hits(evt);
  count_unsliced_big_hits(evt);
  count_unsliced_hit_pairs(evt);
}

// TODO Somehow deal with overlapping triggers?  I found a case where a NuMI
// spill caused three overlapping ddactivity1 triggers, each offset slightly.
// And I noticed it because RemoveBeamSpills fails to remove it.  Ideally, I'd
// drop all the hits that were seen before and only reconstruct the new ones,
// but that is non-trivial.
void ligoanalysis::produce(art::Event & evt)
{
  // Must be called on every event to prevent SIGPIPE loops.
  signal(SIGPIPE, SIG_DFL);

  art::ServiceHandle<geo::Geometry> geo;
  gDet = geo->DetId();

  switch(fAnalysisClass){
    case NDactivity:
      count_triggers(evt);
      count_tracks_containedslices(evt);
      break;
    case JustTrigger:
      count_triggers(evt);
      break;
    case UpMu:
      count_upmu(evt);
      break;
    case LiveTime:
      count_mev(evt);
      count_tracks_containedslices(evt);
      break;
    case DDenergy:
      count_triggers(evt);
      count_ddenergy(evt);
      break;
    case NDMeV:
      count_mev(evt);
      break;
    default:
      printf("No case for type %d\n", fAnalysisClass);
  }
}

DEFINE_ART_MODULE(ligoanalysis);
