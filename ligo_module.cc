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

  }; // class ligo
}

static unsigned long long gwevent_time_us = 0;
static unsigned long long window_size_s = 1000;

namespace ligo{

// Convert an art time, which is a 64 bit number where the upper 32
// bits are the number of seconds since the UNIX Epoch and the lower 32
// bits are the number of microseconds to be added to that, to a double
// which is the number of seconds since the UNIX Epoch. Note that some
// precision is lost in this process since a double only holds about 16
// decimal digits.
static double art_time_to_unix_double(const unsigned long long at)
{
  return (at >> 32) + (at & 0xffffffffULL)*1e-6;
}

/* XXX wait. Leap seconds. GPS. TIA.  UTC.  Oh no. */

// Take a time string like 2000-01-01T00:00:00.123Z and return the time
// in art format. That is, a 64 bit number where the upper 32 bits are
// the number of seconds since the UNIX Epoch and the lower 32 bits are
// the number of microseconds to be added to that.
static unsigned long long rfc3339_to_art_time(const string & stime)
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
  printf("dateforstrptime = %s\n", dateforstrptime.c_str());
  strptime(dateforstrptime.c_str(),
           "%Y-%m-%dT%H:%M:%S", &tm_time);
  printf("Year: %d\n", tm_time.tm_year);
  printf("Month: %d\n", tm_time.tm_mon);
  printf("Day: %d\n", tm_time.tm_mday);
  printf("Hour: %d\n", tm_time.tm_hour);
  printf("Minute: %d\n", tm_time.tm_min);
  int32_t utc_s = mktime(&tm_time); // XXX local time :-(

  utc_s += 6 * 3600; // XXX :-(

  // tzset() followed by correction by timezone and somehow fix DST?

  double f_us = 0;
  const string fractional_second_s =
    stime.substr(rfc3339length, stime.size() - 1 - rfc3339length);

  if(fractional_second_s.size() > 1)
    sscanf(fractional_second_s.c_str(), "%lf", &f_us);

  printf("DEBUG: seconds %d, fractional seconds: %f\n", utc_s, f_us);

  if(f_us < 0 || f_us > 1){
    fprintf(stderr, "Your time string, \"%s\", gave a fractional number of "
            "seconds outside the range [0-1]. I guess it in the wrong format. %s\n",
            stime.c_str(), timehelpmessage);
    exit(1);
  }

  return (((unsigned long long)utc_s) << 32) + (unsigned long long)(f_us * 1e6);
}

ligo::ligo(fhicl::ParameterSet const& pset) : EDAnalyzer(pset),
  fGWEventTime(pset.get<string>("GWEventTime")),
  fWindowSize(pset.get<float>("WindowSize"))
{
  gwevent_time_us = rfc3339_to_art_time(fGWEventTime);
  window_size_s = fWindowSize;
}

ligo::~ligo() { }

/**********************************************************************/
/*                          The meat follows                          */
/**********************************************************************/

// Returns true if the event time is within the window defined by the user
static bool inwindow(const art::Event & evt)
{
  const double evt_time = art_time_to_unix_double(evt.time().value());
  const double gw_time  = art_time_to_unix_double(gwevent_time_us);
  printf("DEBUG: %16f %16f %16f\n", evt_time, gw_time, evt_time - gw_time);
  return fabsl(evt_time - gw_time) < window_size_s/2.;
}


void ligo::analyze(const art::Event & evt)
{
  // I'm so sorry that I have to do this.  And, my goodness, doing
  // it in the constructor isn't sufficient.
  signal(SIGPIPE, SIG_DFL);

  {
    static int NOvA = printf(
      "ntuple: inwindow\n");
    NOvA = NOvA;
  }

  printf("ntuple: %d", inwindow(evt));

  printf("\n");
}

DEFINE_ART_MODULE(ligo);

} // end namespace ligo
//////////////////////////////////////////////////////////////////////////
