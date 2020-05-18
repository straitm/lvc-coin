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

static std::pair<uint32_t, uint32_t>
art_time_minus_some_ns(const unsigned long long art, const double somens)
{
  std::pair<uint32_t, uint32_t> ans;
  ans.first = art >> 32;
  ans.second = art;

  const uint64_t total_ns_to_subtract = somens + 0.5;
  
  const uint32_t s_to_subtract  = total_ns_to_subtract / 1000000000;
  const uint32_t ns_to_subtract = total_ns_to_subtract % 1000000000;

  if(ns_to_subtract > ans.second){
    ans.second += 1000000000;
    ans.first--;
  }

  ans.second -= ns_to_subtract;
  ans.first -= s_to_subtract;

  return ans;
}

// Adds exactly up to rounding to the nearest nanosecond
std::pair<uint32_t, uint32_t>
art_time_plus_some_ns(const unsigned long long art, const double somens)
{
  if(somens < 0) return art_time_minus_some_ns(art, -somens);

  std::pair<uint32_t, uint32_t> ans;
  ans.first = art >> 32;
  ans.second = art;

  const uint64_t total_ns_to_add = somens + 0.5;

  const uint32_t s_to_add  = total_ns_to_add / 1000000000;
  const uint32_t ns_to_add = total_ns_to_add % 1000000000;

  ans.second += ns_to_add;
  if(ans.second > 1000000000){
    ans.second -= 1000000000;
    ans.first++;
  }
  ans.first += s_to_add;

  return ans;
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
