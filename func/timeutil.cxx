#include <string>
#include <string.h>
#include <time.h>

double rfc3339_to_unix_double(const std::string & stime)
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
  std::string dateforstrptime = stime.substr(0, rfc3339length);
  strptime(dateforstrptime.c_str(),
           "%Y-%m-%dT%H:%M:%S", &tm_time);

  setenv("TZ", "", 1); // Make sure we are interpreting the time as UTC
  tzset();
  int32_t unix_s = mktime(&tm_time);

  const std::string fractional_second_s =
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

// Use with caution! The result has less precision than the input.  On the way
// in, there is exactly nanosecond precision.  On the way out, for all dates
// relevant to NOvA the granularity is ~238ns. That's true through the Unix
// Time Apocalypse in 2038.
double art_time_to_unix_double(const unsigned long long at)
{
  return (at >> 32) + (at & 0xffffffffULL)*1e-9;
}

double delta_art_time(const unsigned long long a, const unsigned long long b)
{
  uint32_t a_s = a >> 32;
  uint32_t b_s = b >> 32;
  uint32_t a_ns = a;
  uint32_t b_ns = b;

  // Use 64 bit in case the difference is bigger than (1 << 31)
  int64_t delta_s  = (int64_t)a_s - (int64_t)b_s;

  // These aren't going to overflow because they only go up to 999,999,999
  int32_t delta_ns = (int32_t)a_ns - (int32_t)b_ns;

  return delta_s + delta_ns*1e-9;
}
