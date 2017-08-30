////////////////////////////////////////////////////////////////////////
/// \brief   This module is named ligo and looks at LIGO coincidences
/// \author  M. Strait
/// \date
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"

#include "Geometry/Geometry.h"

#include <DAQDataFormats/RawTriggerTime.h>
#include "DAQDataFormats/RawEvent.h"
#include "DAQDataFormats/RawTrigger.h"
#include "DAQDataFormats/RawTriggerMask.h"
#include "DAQDataFormats/RawDataBlock.h"
#include "DAQDataFormats/RawMicroBlock.h"
#include "DAQDataFormats/RawMicroSlice.h"
#include "DAQDataFormats/RawNanoSlice.h"
#include "DAQChannelMap/DAQChannelMap.h"
#include "StandardRecord/SREnums.h"

#include "RawData/RawSumDropMB.h"
#include "RawData/FlatDAQData.h"
#include "RawData/RawTrigger.h"

#include "NovaTimingUtilities/TimingUtilities.h"

#include "RecoBase/Track.h"

#include "TH1.h"

#include <string>
#include <vector>
#include <set>
using std::string;
using std::vector;
using std::set;
#include <algorithm>

#include <signal.h>

const int TDC_PER_US = 64;
const int US_PER_MICROSLICE = 50; // I hope this is always true

namespace ligo {

class ligo : public art::EDProducer {
  public:
  explicit ligo(fhicl::ParameterSet const& pset);
  virtual ~ligo();
  void produce(art::Event& evt);
  void endJob();

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


  /// True if we see the beginning of the window, i.e. a transition from
  /// being before the window to being in the window.
  bool risingEdge = false;

  /// True if we see the end of the window, i.e. a transition from being
  /// in the window to being after it.
  bool fallingEdge = false;

  /// True if the first event we saw was in the window.
  bool startedHigh = false;

  /// True if we ever see any events in the window.
  bool sawTheWindow = false;

  /// True if we encounter an event in the window after already seeing
  /// a transition from being in the window to being out of the window.
  /// Could conceivably happen if leap seconds are involved, or if the
  /// user constructs some sort of odd file.
  bool tooManyEdges = false;
};

// Set from the FCL parameters
double gwevent_unix_double_time = 0;
long long window_size_s = 1000;

int gDet = caf::kUNKNOWN;

/**********************************************************************/
/*                             Histograms                             */
/**********************************************************************/

struct ligohist{
  // Things per second, for whatever it is we're selecting
  TH1D * sig = NULL;

  // Livetime per second (XXX same for all histograms? I think not necessarily
  // for DDTs)
  TH1D * live = NULL;
};

void init_lh(ligohist & lh, const char * const name)
{
  art::ServiceHandle<art::TFileService> t;
  lh.sig  = t->make<TH1D>(name,                 "",
    window_size_s, -window_size_s/2, window_size_s/2);
  lh.live = t->make<TH1D>(Form("%slive", name), "",
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

// Convert an art time, which is a 64 bit number where the upper 32
// bits are the number of seconds since the UNIX Epoch and the lower 32
// bits are the number of nanoseconds to be added to that, to a double
// which is the number of seconds since the UNIX Epoch.
//
// Note that some precision is lost in this process since a double
// only holds about 16 decimal digits. The granularity at any relevant
// time stamp is 2**-22 seconds = 238 nanoseconds, which does not matter
// for these purposes.
double art_time_to_unix_double(const unsigned long long at)
{
  return (at >> 32) + (at & 0xffffffffULL)*1e-9;
}

// Take a time string like 2000-01-01T00:00:00.123Z and return the time
// in Unix time with fractional seconds. That is, a floating point
// number where the integer part is the same as what you get from "date
// +%s", the number of seconds since the UNIX Epoch, ignoring leap
// seconds, and there's also a fractional part.
double rfc3339_to_unix_double(const string & stime)
{
  tm tm_time;
  memset(&tm_time, 0, sizeof(tm));
  const unsigned int rfc3339length = sizeof("2000-01-01T00:00:00") - 1;

  const char * const timehelpmessage = "Must look like "
    "2000-01-01T00:00:00.123Z where the Z denotes UTC time "
    "and you have to give it in UTC time.  The fractional part "
    "of the second is optional.";

  if(stime.size() < rfc3339length+1 || stime[stime.size()-1] != 'Z'){
    fprintf(stderr, "Malformed time string for LIGO. %s\n", timehelpmessage);
    exit(1);
  }
  string dateforstrptime = stime.substr(0, rfc3339length);
  strptime(dateforstrptime.c_str(),
           "%Y-%m-%dT%H:%M:%S", &tm_time);

  setenv("TZ", "", 1); // Make sure we are interpreting the time as UTC
  tzset();
  int32_t unix_s = mktime(&tm_time);

  const string fractional_second_s =
    stime.substr(rfc3339length, stime.size() - 1 - rfc3339length);

  double unix_fraction = 0;
  if(fractional_second_s.size() > 1)
    sscanf(fractional_second_s.c_str(), "%lf", &unix_fraction);

  if(unix_fraction < 0 || unix_fraction >= 1){
    fprintf(stderr, "Your time string, \"%s\", gave fractional seconds outside"
            "the range [0-1). I guess it is in the wrong format. %s\n",
            stime.c_str(), timehelpmessage);
    exit(1);
  }

  return unix_s + unix_fraction;
}

void ligo::endJob()
{
  if(tooManyEdges)                     printf("Data was out of time order! :-0\n");
  else if( risingEdge &&  fallingEdge) printf("Saw a whole window :-)\n");
  else if( risingEdge && !fallingEdge) printf("Saw only beginning of window :-\\\n");
  else if(!risingEdge &&  fallingEdge) printf("Saw only end of the window :-\\\n");
  else if(sawTheWindow)                printf("Saw only middle of window >:-\\\n");
  else                                 printf("Didn't see any data in window :-(\n");
}

// Number of hits, with no filtering of any sort
ligohist lh_rawhits;

// Number of hits that are not in any Slicer4D slice, i.e. they are in the
// Slicer4D noise slice.
ligohist lh_unslice4ddhits;

// Raw number of tracks
ligohist lh_tracks;

// Number of tracks with at least one contained endpoint, with some additional
// sanity checks.
ligohist lh_halfcontained_tracks;

// Number of tracks with at two contained endpoints, with some additional
// sanity checks.
ligohist lh_fullycontained_tracks;

ligo::ligo(fhicl::ParameterSet const& pset) : EDProducer(),
  fGWEventTime(pset.get<string>("GWEventTime")),
  fWindowSize(pset.get<unsigned long long>("WindowSize"))
{
  gwevent_unix_double_time = rfc3339_to_unix_double(fGWEventTime);
  window_size_s = fWindowSize;

  init_lh(lh_rawhits,          "rawhits");
  init_lh(lh_unslice4ddhits,   "unslice4ddhits");
  init_lh(lh_tracks,           "tracks");
  init_lh(lh_halfcontained_tracks, "halfcontained_tracks");
  init_lh(lh_fullycontained_tracks, "fullycontained_tracks");
}

ligo::~ligo() { }

/**********************************************************************/
/*                          The meat follows                          */
/**********************************************************************/

// Returns true if the event time is within the window defined by the
// user.
//
// I wondered if evt.time() correctly returns true Unix time or if there
// was perhaps a leap second problem of one flavor or another in there.
// Since absolute timing doesn't matter for anything in NOvA except
// SNEWS triggering and LIGO follow-ups, this is certainly a risk.
//
// Here is code to test that by examining the raw trigger. It
// shows that evt.time() is always within a microsecond or so of
// fTriggerTimingMarker_TimeStart converted into Unix time by code that
// explicitly does a leap second correction:
//
////////////////////////////////////////////////////////////////////////
//  art::Handle< std::vector<rawdata::RawTrigger> > rawtrigger;
//  evt.getByLabel("daq", rawtrigger);
//
//  // Consulting example code in
//  // Online/NovaGlobalTrigger/src/test/GTReceiver.cc
//  // But I have a rawdata::RawTrigger.  I need a daqdataformats::RawTrigger
//  // Or I just need the triggertime directly.  Anyway, I can't directly
//  // do this:
//  //
//  // daqdataformats::RawTriggerTime* triggerTime =
//  //   (*rawtrigger)[0].getTriggerTime();
//  //
//  // daqdataformats:RawTrigger is full of "clever" preprocesor macros
//  // that make it impossible to read or understand. Sure, I'm guilty of
//  // doing this from time to time, but usually not in code I intend a
//  // whole collaboration to use.
//
//  timeval tv; // seconds and microseconds
//
//  // Args in the opposite of the conventional order. tv gets set.
//  novadaq::timeutils::convertNovaTimeToUnixTime(
//    // Actually want fTriggerTimingMarker_ExtractionStart, but that is
//    // always zero, which is probably a bug. TimeStart is the time the
//    // trigger *asked* for rather than the one it *got*, but that's
//    // probably close enough.
//   (*rawtrigger)[0].fTriggerTimingMarker_TimeStart, tv);
//
//  const double evt_time_from_header = tv.tv_sec + tv.tv_usec * 1e-6;
//  printf("DEBUG: %6u %16f %16f %16f : %16f %16f : %.9f\n", evt.id().event(),
//    evt_time, evt_time_from_header, gw_time,
//    evt_time - gw_time, evt_time_from_header - gw_time,
//    evt_time - evt_time_from_header);
////////////////////////////////////////////////////////////////////////
//
// Now, Unix time is not a great time stamp because a leap second and
// the second prior to a leap second are represented by the same value.
// It's like using Daylight Savings Time and not having any marker for
// whether or not it is in effect. This means that if there is a leap
// second in our window it is going to be 1s too long and/or have 2
// seconds of events piled up in 1s, depending on how the rest of this
// module is implemented. If there were to be a negative leap second
// (never happened as of 2017), there'd be the opposite problem.
//
// We could solve this problem by using the time from the header and
// converting it ourselves *without* leap seconds, so it is the (simple
// monotonic) number of seconds since the NOvA epoch. But we'd have to
// correct the user supplied time then, and different headaches ensue.
// Probably better to check by hand that we're not near a leap second,
// which is pretty unlikely since they come only every year or two, and
// so far gravitational wave events only come a few per year. There's
// only a ~2e-5 chance of it being a problem for a given event with a
// 1000s window. I don't know. Watch out for events at the end of June
// or December (or, to some extent, the end of any month).
bool inwindow(const art::Event & evt)
{
  const double evt_time = art_time_to_unix_double(evt.time().value());
  printf("DEBUG: %6u %16f %16f %16f\n", evt.id().event(),
    evt_time, gwevent_unix_double_time, evt_time - gwevent_unix_double_time);
  return fabsl(evt_time - gwevent_unix_double_time) < window_size_s/2.;
}


// Return the livetime in this event in seconds, as it is relevant for
// raw hits (same as for anything else? Maybe not if we don't trust
// tracks close to the time-edges of events).
double rawlivetime(const art::Event & evt)
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
void THplusequals(ligohist & lh, const int bin, const double sig, const double live)
{
  // Use SetBinContent instead of Fill(x, weight) to avoid having to look up
  // the bin number twice.
  lh.sig ->SetBinContent(bin, lh.sig ->GetBinContent(bin) + sig);
  lh.live->SetBinContent(bin, lh.live->GetBinContent(bin) + live);
}

// Return the bin number for this event, i.e. the number of seconds from the
// beginning of the window, plus 1.
int timebin(const art::Event & evt)
{
  const double evt_time = art_time_to_unix_double(evt.time().value());
  return floor(evt_time - gwevent_unix_double_time) + window_size_s/2
         + 1; // stupid ROOT 1-based numbering!
}

bool contained(const TVector3 & v)
{
  if(gDet == caf::kNEARDET)
    return fabs(v.X()) < 165 && fabs(v.Y()) < 165 && v.Z() > 25 && v.Z() < 1225;
  if(gDet == caf::kFARDET)
    return fabs(v.X()) < 700 && v.Y() < 500 && v.Y() > -700 && v.Z() > 25 && v.Z() < 5950;
  fprintf(stderr, "Unknown detector %d\n", gDet);
  exit(1);
}

// Returns true if the track enters and exits.  Or if both of its ends
// are pretty close to the edge.  But not if it is just steep.
bool un_contained_track(const rb::Track & t)
{
  return !contained(t.Start()) && !contained(t.Stop());
}

const int min_plane_extent = 10;

// Returns true if the track starts AND stops inside the detector, 
// and we're pretty sure that this isn't artifactual.
bool fully_contained_track(const rb::Track & t)
{
  return contained(t.Start()) && contained(t.Stop()) &&
          fabs(t.Dir().Y()) < 0.95 &&
          fabs(t.Dir().X()) < 0.95 &&
          t.ExtentPlane() >= min_plane_extent;
}

// Returns true if either end of the track is uncontained or might be
bool half_uncontained_track(const rb::Track & t)
{
  return !contained(t.Start()) || !contained(t.Stop()) ||
          fabs(t.Dir().Y()) > 0.95 ||
          fabs(t.Dir().X()) > 0.95 ||
          t.ExtentPlane() < min_plane_extent;
}

// Returns true if the track either starts or stops in the detector, or both, 
// and we're pretty sure that this isn't artifactual.
bool half_contained_track(const rb::Track & t)
{
  return (contained(t.Start()) || contained(t.Stop())) &&
          fabs(t.Dir().Y()) < 0.95 &&
          fabs(t.Dir().X()) < 0.95 &&
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

// Count the number of raw hits in the event and fill the appropriate histograms
void count_hits(const art::Event & evt)
{
  art::Handle< std::vector<rawdata::RawDigit> > rawhits;
  evt.getByLabel("daq", rawhits);

  printf("Hits in this event: %lu\n", rawhits->size());

  THplusequals(lh_rawhits, timebin(evt), rawhits->size(), rawlivetime(evt));
}

void count_unslice4dd_hits(const art::Event & evt)
{
  art::Handle< std::vector<rb::Cluster> > slice;
  evt.getByLabel("slicer", slice);

  if(slice->empty()){
    printf("Unexpected event with zero slices!\n");
    return;
  }

  printf("Hits in this event in the Slicer4D noise slice: %d\n", (*slice)[0].NCell());

  THplusequals(lh_unslice4ddhits, timebin(evt), (*slice)[0].NCell(), rawlivetime(evt));
}

void count_tracks(const art::Event & evt)
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
  vector<int> trk_slices(tracks->size(), -1);
  std::set<int> slices_with_uc_tracks, slices_with_huc_tracks;
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
    // track that either enters or exits.
    if(!slices_with_huc_tracks.count(trk_slices[i]) &&
       fully_contained_track((*tracks)[i]))
      slices_with_fc_tracks.insert(trk_slices[i]);
  }

  printf("Slices with half-contained tracks: %lu\n", slices_with_hc_tracks.size());
  printf("Slices with fully-contained tracks: %lu\n", slices_with_fc_tracks.size());
  for(set<int>::iterator i = slices_with_fc_tracks.begin();
      i != slices_with_fc_tracks.end(); i++)
    printf("  %3d\n", *i);

  THplusequals(lh_halfcontained_tracks,  timebin(evt),
               slices_with_hc_tracks.size(), rawlivetime(evt));
  THplusequals(lh_fullycontained_tracks, timebin(evt),
               slices_with_fc_tracks.size(), rawlivetime(evt));
}

void ligo::produce(art::Event & evt)
{
  // Must be called on every event to prevent SIGPIPE loops.
  signal(SIGPIPE, SIG_DFL);

  art::ServiceHandle<geo::Geometry> geo;
  gDet = geo->DetId();

  const bool is_in_window = inwindow(evt);

  // Keep track of whether we've seen the beginning and end of the window.
  {
    static bool firstevent = true;
    if(is_in_window)                                 sawTheWindow = true;
    if(firstevent && is_in_window)                   startedHigh  = true;
    if(!startedHigh && is_in_window && !risingEdge)  risingEdge   = true;
    if(is_in_window && fallingEdge)                  tooManyEdges = true;
    if((startedHigh || risingEdge) && !is_in_window) fallingEdge  = true;
    firstevent = false;
  }

  if(is_in_window){
    count_hits(evt);
    count_unslice4dd_hits(evt);
    count_tracks(evt);
  }

  printf("In window: %d\n", is_in_window);
}

DEFINE_ART_MODULE(ligo);

} // end namespace ligo
//////////////////////////////////////////////////////////////////////////
