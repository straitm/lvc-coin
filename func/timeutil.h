// Take a time string like 2000-01-01T00:00:00.123Z and return the time
// in Unix time with fractional seconds. That is, a floating point
// number where the integer part is the same as what you get from "date
// +%s", the number of seconds since the UNIX Epoch, ignoring leap
// seconds, and there's also a fractional part.
double rfc3339_to_unix_double(const std::string & stime);

// Convert an art time, which is a 64 bit number where the upper 32
// bits are the number of seconds since the UNIX Epoch and the lower 32
// bits are the number of nanoseconds to be added to that, to a double
// which is the number of seconds since the UNIX Epoch.
//
// Note that some precision is lost in this process since a double
// only holds about 16 decimal digits. The granularity at any relevant
// time stamp is 2**-22 seconds = 238 nanoseconds, which does not matter
// for my purposes.
double art_time_to_unix_double(const unsigned long long at);
