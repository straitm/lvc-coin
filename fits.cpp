#include <stdio.h>
#include <math.h>
#include "healpix_base.h"
#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include <vector>
#include <algorithm>
using std::vector;
using std::pair;

// g++ -o test test.cpp -I$HOME/Healpix_3.31/src/cxx/optimized_gcc/include
// -L$HOME/Healpix_3.31/src/cxx/optimized_gcc/lib -L$HOME/lib -lhealpix_cxx
// -lcxxsupport -lsharp -lfftpack -lc_utils -lcfitsio -lcurl &&
// LD_LIBRARY_PATH=$HOME/Healpix_3.31/lib ./test

bool comp(const pair<double, double> & a, const pair<double, double> & b)
{
  return a.second > b.second;
}

int main()
{
  Healpix_Map<float> map;

  // LOCALIZATION AND BROADBAND FOLLOW-UP OF GW150914
  // "Since GW150914 is a CBC event, we
  // consider the LALInference map to be the most accurate, au-
  // thoritative, and final localization for this event."
  read_Healpix_map_from_fits("LALInference_skymap.fits", map);

  double sumprob = 0;

  vector< pair<double, double> > vals; // area corrected, raw

  int ni = 1000, nj = 1000;
  for(int i = 1; i < ni; i++){
    for(int j = 0; j < nj; j++){
      const double theta = i*M_PI/ni,
                   phi   = j*2.*M_PI/nj;
      const float val = map.interpolated_value(pointing(theta, // dec - except this is 0-pi, and dec is pi/2 to -pi/2
                                                        phi    // ra - straight up 0h = 0, 24h = 2pi
                                               ));
      sumprob += val * sin(theta);
      vals.push_back(pair<double, double>(val*sin(theta), val));
    }
  }

  sort(vals.begin(), vals.end(), comp); 

  float acc = 0;
  float critval = 0;
  for(int i = 0; i < vals.size(); i++){
    acc += vals[i].first/sumprob;
    if(acc > 0.9 /* CL */){
      critval = vals[i].second;
      break;
    }
  }

  for(int i = 0; i < 80; i++){
    for(int j = 159; j >= 0; j--){
      const double theta = i*M_PI/80, phi = j*2*M_PI/160;
      const float val = map.interpolated_value(pointing(theta, phi));
      printf("%c", val > critval?'X':'-');
    }
    printf("\n");
  }

  return 0;
}
