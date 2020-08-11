using namespace std;

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <algorithm>

// Size of character buffers that we're going to snprintf into. This is
// way more than is actually needed, so we're not going to check whether
// snprintf output was truncated. In the extremely unlikely case that it
// is, the only consequence is that the user doesn't get to see the end
// of the string.
static const int CHARMAX = 256;

// Given a double greater than 1, round to 2 digits or less
static int sigfigs(const double in)
{
  if(in >= 2.147483648e9){
    fprintf(stderr, "I'm not going to be able to store %f in an int!\n", in);
    return 2147483647;
  }

  int n = int(in+0.5);
  int divided = 0;
  if(n > 100){
    int lastdig = n%10;
    while((n /= 10) > 100){
      lastdig = n%10;
      divided++;
    }
    if(lastdig >= 5) n++;
    for(int i = 0; i <= divided; i++) n *= 10;
  }
  return n;
}

static void formatestimate(char * answer, const int etasec)
{
  if(etasec < 55){ // 1-55 etasec
    snprintf(answer, CHARMAX, "%ds", etasec);
  }
  else if(etasec < 570){ // 1 minute to 9m30s
    snprintf(answer, CHARMAX, "%dm%ds", (etasec+5)/60,
            (etasec+5)%60/10*10);
  }
  else if(etasec < 3570){ // 10 minutes to 59 minutes
    snprintf(answer, CHARMAX, "%dm", (etasec+30)/60);
  }
  else if(etasec < 34200){ // 1 hour to 9h50m
    snprintf(answer, CHARMAX, "%dh%dm",
             (etasec+300)/3600, (etasec+300)%3600/600*10);
  }
  else if(etasec < 84600){ // 10 hours to 23 hours
    snprintf(answer, CHARMAX, "%dh", (etasec+1800)/3600);
  }
  // This one is weird. What is two sig figs when you're between 1 and
  // 10 days? Well, a tenth of a day is about 2 hours, so let's use
  // that: 1 day to 9d22h
  else if(etasec < 856800){ 
    snprintf(answer, CHARMAX, "%dd%dh", 
             (etasec+3600)/86400, (etasec+3600)%86400/7200*2);
  }
  else{
    snprintf(answer, CHARMAX, "%dd", sigfigs(double(etasec)/86400));
  }
  // Easy enough to extend to weeks, etc if you like. If my programs
  // want to run that long, I either rewrite them, use a bigger cluster,
  // or change my goals.
}

static int etasec(const double ince, const double tote, 
                  const double frac)
{ 
  return int(round(ince == -1? tote: 
             frac < 0.5? (tote+ince)/2:
             frac*ince+(1-frac)*tote) + 0.5);
}

// Returns a string describing the estimated time left. Returns a
// pointer to a string representing the time to be printed. Caller must
// free the string when done with it.
static char * eta(const double ince, const double tote, 
                  const double frac)
{
  char * answer = (char*)malloc(CHARMAX);

  formatestimate(answer, etasec(ince, tote, frac)); 

  return answer;
}

// Returns a string describing the estimated total running time. Returns
// a pointer to a string representing the time to be printed. Caller
// must free the string when done with it. If the total time improves by
// 10% or more, sets status to 1. If it gets worse by 10% or more, sets
// status to -1. Otherwise, sets it to 0.
static char * etotal(int & status, const double tottime, 
                     const double ince, const double tote, 
                     const double frac)
{
  char * answer = (char*)malloc(CHARMAX);

  static int previous;
  static bool firsttime = true; 
  int current = etasec(ince, tote, frac) + int(tottime);

  formatestimate(answer, current); 

  if(firsttime || fabs(previous - current) < 0.1*current) status = 0;
  else if(previous < current) status = -1;
  else status = 1;

  previous = current;
  firsttime = false;

  return answer;
}

// Returns a pointer to a string representing the time to be printed.
// Caller must free the string when done with it.
static char * disptime(const double ttime)
{
  char * buf = (char*)malloc(CHARMAX);

  // Could upgrade to 64-bit, but seriously, it would only be useful
  // for programs running longer than 78 years, so who cares?
  if(ttime >= 2.147483648e9){
    fprintf(stderr, "I'm not going to be able to store %f in an int!\n",
            ttime);
    snprintf(buf, CHARMAX, "more than 78 years");
    return buf;
  }

  int t = int(ttime);
  bool reqzero = false, showseconds = true;
  int printed = 0; // keep track of how many characters have been used

  if(t >= 86400){
    printed += snprintf(buf, CHARMAX, "%dd ",  t/86400); 
    t %= 86400;
    reqzero = true;
    showseconds = false;
  }

  if(t >= 3600 || reqzero){
    printed += snprintf(buf+printed, CHARMAX-printed, 
                        "%0*dh%s", reqzero?2:1, t/3600, reqzero?"":" ");
    t %= 3600;
    reqzero = true;
  }

  if(t >= 60 || reqzero){
    printed += snprintf(buf+printed, CHARMAX-printed, 
                        "%0*dm", reqzero?2:1, t/60);
    t %= 60;
    reqzero = true;
  }

  if(showseconds)
    snprintf(buf+printed, CHARMAX-printed, "%0*ds", reqzero?2:1, t);

  return buf;
} 


// Not meant to have any generality. Just a helper function for
// generateprintpoints.
static uint64_t iexp10(const int ep)
{
  switch(ep){
    case 2: return 100;
    case 3: return 1000;
    case 4: return 10000;
    case 5: return 100000;
    case 6: return 1000000;
    case 7: return 10000000;
    case 8: return 100000000;
    case 9: return 1000000000;
    default: 
      fprintf(stderr, "Bad exponent in iexp10()\n");
      return 1;
  }
}

/* Given the total number of events and the most digits to print in the
reports, generate the events on which progress should be reported. */
vector<uint64_t> generateprintpoints(const uint64_t total,
                                     const int maxe)
{
  vector<uint64_t> ppoints;
  // makes 10% - 90% print
  // Look at the parentheses around (total/10).  Those are required
  // to make this work for numbers bigger than 0xffffffff/10
  for(uint64_t i = 1; i <= 9; i++) ppoints.push_back(i * (total/10));

  // Makes 1%-9% and 91%-99%, 0.1%-0.9% and 99.1%-99.9%, etc.
  for(int ep = 2; ep <= maxe; ep++){
    for(uint64_t i = 1; i <= 9; i++){
      ppoints.push_back(i * (total/iexp10(ep)));
      ppoints.push_back(total - i * (total/iexp10(ep)));
    }
  }

  // cat ppoints | sort | uniq
  sort(ppoints.begin(), ppoints.end());
  ppoints.resize(unique(ppoints.begin(), ppoints.end()) - ppoints.begin());

  return ppoints;
}

// Each time, new is set to the current time. old is set to the current
// time the first time, then subsequently is set to new at the bottom
static int firsttime, oldtime;
static vector<uint64_t> ppoints; // the values of sofar to print
static uint64_t nextprint = 0xffffffff;
static double lastfrac;
static uint64_t total;

static void printprogress(const uint64_t sofar)
{
  // we're never going to find this one or any one before it again,
  // so remove them to save time in future searches.
  if(ppoints.size()){
    ppoints.erase(ppoints.begin());
    nextprint = ppoints[0];
  }

  double frac = double(sofar)/total;

  int newtime = time(NULL);
  double tottime = newtime - firsttime;
  double inctime = newtime - oldtime;

  // If the total time or incremental time is
  // too short, don't print or save timestamp
  if(tottime < 3 || inctime < 1) return;

  double tote = tottime/frac - tottime;
  double ince = frac-lastfrac > 0? (1-frac)*inctime/(frac-lastfrac): -1;

  int ep;
  if     (frac < 0.0000000099 || frac > 0.99999999) ep = 8;
  else if(frac < 0.000000099  || frac > 0.9999999)  ep = 7;
  else if(frac < 0.00000099   || frac > 0.999999)   ep = 6;
  else if(frac < 0.0000099    || frac > 0.99999)    ep = 5;
  else if(frac < 0.000099     || frac > 0.9999)     ep = 4;
  else if(frac < 0.00099      || frac > 0.999)      ep = 3;
  else if(frac < 0.0099       || frac > 0.99)       ep = 2;
  else if(frac < 0.099        || frac > 0.9)        ep = 1;
  else                                              ep = 0;

  char * thedisptime = disptime(tottime);
  char * theeta = eta(ince, tote, frac);
  int status;
  char * theetot = etotal(status, tottime, ince, tote, frac);

  // If your background wasn't black, this makes it black; You'll have
  // to say 'reset' afterwards if you don't like black.
  // Colors: 37=white, 31=red, 32=green
  printf("****"
         " %8.*f%% "
         "Elapsed: %-10s  ETA: %-6s  Est total: %c[%s;%s;40m%-6s%c[0;37;40m "
         "****\n",
         ep-1 > 0? ep-1: 0, 
         round(pow(10, ep+1)*frac)/pow(10, ep-1),
         thedisptime, 
         theeta,
         0x1b,
         status == 0?"0":"1",
         status == 0?"37":status == 1?"32":"31",
         theetot, 
         0x1b);

  fflush(stdout);

  free(thedisptime);
  free(theeta);

  oldtime = newtime;
  lastfrac = frac;
}

void initprogressindicator(const uint64_t totin, const int maxe)
{
  if     (maxe > 9) fprintf(stderr, "maxe may not be > 9. Using 9\n");
  else if(maxe < 1) fprintf(stderr, "maxe may not be < 1. Using 1\n");
   
  total = totin;
  ppoints = generateprintpoints(total, maxe>9? 9: maxe<1? 1: maxe);
  nextprint = ppoints[0];
  firsttime = time(NULL);
  oldtime = firsttime;
  lastfrac = 0;
}

// If this does not print, all it does is compares two integers and then
// returns. The exception to this is if it is not printing because not
// enough time has elapsed since the previous print. This is about as
// rare as the prints themselves, so not a performance concern.
//
// inline makes a serious impact (factor of 10 improvement), at least if
// using -O2 and if this cpp file is included directly rather than via
// the header. (If using the header, this file is compiled separately
// and so *can't* be inlined.)
//
// Without inline, with or without including this file directly, it is
// a factor of 2 speed improvement to have printprogress as a separate
// function, presumably because we don't invalidate as much of the cache
// or pipeline or something by having a big dangling function body that
// turns out to be unused.
#ifdef PROGRESS_INDICATOR_HEADER_USED
  void progressindicator(const uint64_t sofar)
#else
  inline void progressindicator(const uint64_t sofar)
#endif
{
  if(sofar == nextprint) printprogress(sofar);
}
