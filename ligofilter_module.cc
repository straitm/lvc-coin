////////////////////////////////////////////////////////////////////////
/// \brief   This module filters out a time window for the GW analysis
/// \author  M. Strait
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include <string>

#include "func/timeutil.h"

// Set from the FCL parameters
static double gwevent_unix_double_time = 0;
static long long window_size_s = 1000;

class ligofilter : public art::EDFilter {
  public:
  explicit ligofilter(const fhicl::ParameterSet & pset);
  virtual ~ligofilter() { }; // compiles, but does not run, without this
  bool filter(art::Event& evt);
  void endJob();

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
  /// Observed in some DDenergy files. Could also conceivably happen if
  /// leap seconds are involved, or if the user constructs some sort of
  /// odd file.
  bool tooManyEdges = false;
};

ligofilter::ligofilter(fhicl::ParameterSet const& pset) : EDFilter(),
  fGWEventTime(pset.get<std::string>("GWEventTime")),
  fWindowSize(pset.get<unsigned long long>("WindowSize"))
{
  gwevent_unix_double_time = rfc3339_to_unix_double(fGWEventTime);
  window_size_s = fWindowSize;
}

void ligofilter::endJob()
{
  if(tooManyEdges)                     printf("Data out of time order :-0\n");
  else if( risingEdge &&  fallingEdge) printf("Saw whole window :-)\n");
  else if( risingEdge && !fallingEdge) printf("Saw beginning of window :-\\\n");
  else if(!risingEdge &&  fallingEdge) printf("Saw end of window :-/\n");
  else if(sawTheWindow)                printf("Saw middle of window >:-|\\n");
  else                                 printf("Saw no data in window :-(\n");
}

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
bool ligofilter::filter(art::Event & evt)
{
  const double evt_time = art_time_to_unix_double(evt.time().value());
#if 1
  // This is very useful for understanding if we got the whole window
  // in files like DDEnergy where the events are only loosely ordered.
  printf("DEBUG: %6u %16f %16f %16f\n", evt.id().event(),
    evt_time, gwevent_unix_double_time, evt_time - gwevent_unix_double_time);
#endif
  const bool inwindow =
    fabsl(evt_time - gwevent_unix_double_time) < window_size_s/2.;

  // Keep track of whether we've seen the beginning and end of the window.
  static bool firstevent = true;
  if(inwindow)                                 sawTheWindow = true;
  if(firstevent && inwindow)                   startedHigh  = true;
  if(!startedHigh && inwindow && !risingEdge)  risingEdge   = true;
  if(inwindow && fallingEdge)                  tooManyEdges = true;
  if((startedHigh || risingEdge) && !inwindow) fallingEdge  = true;
  firstevent = false;

  if(inwindow) printf("In GW coincidence window\n");
  return inwindow;
}
DEFINE_ART_MODULE(ligofilter)
