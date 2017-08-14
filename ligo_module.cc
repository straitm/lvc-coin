////////////////////////////////////////////////////////////////////////
/// \brief   This module is named ligo and looks at LIGO coincidences
/// \author  M. Strait
/// \date
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include <DAQDataFormats/RawTriggerTime.h>
#include "DAQDataFormats/RawEvent.h"
#include "DAQDataFormats/RawTrigger.h"
#include "DAQDataFormats/RawTriggerMask.h"
#include "DAQDataFormats/RawDataBlock.h"
#include "DAQDataFormats/RawMicroBlock.h"
#include "DAQDataFormats/RawMicroSlice.h"
#include "DAQDataFormats/RawNanoSlice.h"
#include "DAQChannelMap/DAQChannelMap.h"
#include "RawData/RawSumDropMB.h"

#include "RawData/FlatDAQData.h"
#include "RawData/RawTrigger.h"

#include "NovaTimingUtilities/TimingUtilities.h"


#include "RecoBase/Track.h"

#include <string>
using std::string;
#include <algorithm>

#include <signal.h>

namespace ligo {

class ligo : public art::EDAnalyzer {
  public:
  explicit ligo(fhicl::ParameterSet const& pset);
  virtual ~ligo();
  void analyze(const art::Event& evt);

  string fGWEventTime;
  float fWindowSize;
};

unsigned long long gwevent_unix_double_time = 0;
unsigned long long window_size_s = 1000;

// Convert an art time, which is a 64 bit number where the upper 32
// bits are the number of seconds since the UNIX Epoch and the lower 32
// bits are the number of nanoseconds to be added to that, to a double
// which is the number of seconds since the UNIX Epoch.
//
// Note that some precision is lost in this process since a double
// only holds about 16 decimal digits. The granularity at any relevant
// timestamp is 2**-22 seconds = 238 nanoseconds, which does not matter
// for these purposes.
double art_time_to_unix_double(const unsigned long long at)
{
  return (at >> 32) + (at & 0xffffffffULL)*1e-9;
}

// Convert a nova timestamp which represents Unix time with integer
// and fractional seconds, into a double that also represents Unix time.
// as with art_time_to_unix_double(), this incurs an unimportant
// loss of precision.
//
// For the name of the argument to this function cf. 
// http://www.catb.org/esr/time-programming/#_broken_down_time
// "In what is probably unintentional humor, various manual pages and
// standards documents refer to it as 'broken-down time'."
double nova_unix_time_to_double(const timeval & nova_broken_down_time)
{
  return nova_broken_down_time.tv_sec +
         nova_broken_down_time.tv_usec * 1e-6;
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

ligo::ligo(fhicl::ParameterSet const& pset) : EDAnalyzer(pset),
  fGWEventTime(pset.get<string>("GWEventTime")),
  fWindowSize(pset.get<float>("WindowSize"))
{
  gwevent_unix_double_time = rfc3339_to_unix_double(fGWEventTime);
  window_size_s = fWindowSize;
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
//  const double evt_time_from_header = nova_unix_time_to_double(tv);
// 
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
// 1000s window. I don't know.
bool inwindow(const art::Event & evt)
{
  const double evt_time = art_time_to_unix_double(evt.time().value());
  printf("DEBUG: %6u %16f %16f %16f\n", evt.id().event(),
    evt_time, gwevent_unix_double_time, evt_time - gwevent_unix_double_time);
  return fabsl(evt_time - gwevent_unix_double_time) < window_size_s/2.;
}


void ligo::analyze(const art::Event & evt)
{
  // Must be called on every event to prevent SIGPIPE loops.
  signal(SIGPIPE, SIG_DFL);

  {
    static int NOvA = printf("ntuple: inwindow\n");
    NOvA = NOvA;
  }

  printf("ntuple: %d", inwindow(evt));

  printf("\n");
}

DEFINE_ART_MODULE(ligo);

} // end namespace ligo
//////////////////////////////////////////////////////////////////////////
