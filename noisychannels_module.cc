////////////////////////////////////////////////////////////////////////
/// \brief The noisychannels module looks for noisy channels
///
/// \author M. Strait
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "DAQDataFormats/RawTriggerTime.h"
#include "DAQDataFormats/RawEvent.h"
#include "DAQDataFormats/RawTrigger.h"
#include "DAQDataFormats/RawTriggerMask.h"
#include "DAQDataFormats/RawDataBlock.h"
#include "StandardRecord/SREnums.h"

#include "RawData/FlatDAQData.h"
#include "ChannelInfo/BadChanList.h"
#include "RawData/RawDigit.h"
#include "RawData/RawTrigger.h"

#include "RecoBase/Track.h"

#include "TH2.h"
#include "TTree.h"

#include <vector>
#include <algorithm>

#include <signal.h>

#include "progress.cpp"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "Calibrator/Calibrator.h"
#include "ChannelInfo/BadChanList.h"
#include "RawData/RawDigit.h"
#include "RecoBase/CellHit.h"
#include "RunHistory/service/RunHistoryService.h"

// Always make it FD size and sometimes don't use the whole thing
static const unsigned int nplane = 32*28, ncell = 12*32;

static uint64_t counts[nplane][ncell] = { 0 };
static uint64_t hicounts[nplane][ncell] = { 0 };
static TH2C * chmask = NULL;
static TH2D * hrates = NULL, * hhirates = NULL;

class noisychannels : public art::EDAnalyzer {
  public:
  explicit noisychannels(fhicl::ParameterSet const& pset);
  virtual ~noisychannels() { }; // compiles, but does not run, without this
  void analyze(const art::Event& evt);

  void beginJob();
  void endJob();
};

noisychannels::noisychannels(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset)
{
}

void noisychannels::beginJob()
{
  art::ServiceHandle<art::TFileService> t;

  chmask = t->make<TH2C>("chmask",  "", nplane, 0, nplane, ncell, 0, ncell);
  hrates = t->make<TH2D>("rates",   "", nplane, 0, nplane, ncell, 0, ncell);
  hhirates=t->make<TH2D>("hirates", "", nplane, 0, nplane, ncell, 0, ncell);
}

void noisychannels::endJob()
{
  // Count active channels so this works regardless of detector,
  // and even works for partial detectors. 
  uint64_t activechannels = 0;

  vector<uint64_t> countssort, hicountssort;
  
  for(unsigned int i = 0; i < nplane; i++){
    for(unsigned int j = 0; j < ncell; j++){
      if(counts[i][j] > 0){
        activechannels++;
        countssort.push_back(counts[i][j]);
      }
      if(hicounts[i][j] > 0)
        hicountssort.push_back(hicounts[i][j]);
    }
  }

  printf("%lu channel%s seem%s active\n", activechannels,
         activechannels == 1?"":"s",
         activechannels == 1?"s":"");

  if(activechannels == 0) return;

  std::sort(countssort.begin(), countssort.end());
  std::sort(hicountssort.begin(), hicountssort.end());

  const double median = countssort[countssort.size()/2];
  const double himedian = hicountssort.size()?
    hicountssort[hicountssort.size()/2]:0;

  const double threshold = 10*median,
             hithreshold = 10*himedian;

  for(unsigned int i = 0; i < nplane; i++){
    for(unsigned int j = 0; j < ncell; j++){
      // We're operating entirely with bin numbers and not bin ranges,
      // so this is not an off-by-one error.
      if(counts[i][j] > threshold) chmask->SetBinContent(i, j, 1);
      if(hicounts[i][j] > hithreshold) chmask->SetBinContent(i, j, 1);
      hrates  ->SetBinContent(i, j, counts  [i][j]/median);
      hhirates->SetBinContent(i, j, hicounts[i][j]/himedian);
    }
  }
}

void noisychannels::analyze(const art::Event & evt)
{
  {
    static unsigned int n = 0;
    // Start at 1 because the first event takes forever as everything
    // is initialized.
    if(n == 1){
      // art provides no way of knowing how much events we will process,
      // and its own progress indicator is of limited use.
      printf("For progress indicator, assuming a 2000 event long readout\n");
      initprogressindicator(2000-1, 3);
    }
    else if(n > 1){
      progressindicator(n - 1);
    }
    n++;
  }

  // Must be called on every event to prevent SIGPIPE loops.
  signal(SIGPIPE, SIG_DFL);

  art::Handle< std::vector<rawdata::RawDigit> > hits;
  evt.getByLabel("daq", hits);

  art::ServiceHandle<cmap::CMap> map;

  for(const auto & hit : *hits){
    const int p = map->GetPlane(&hit),
              c = map->GetCell(&hit);
    counts[p][c]++;
    if(hit.ADC() > 100) hicounts[p][c]++;
  }
}

DEFINE_ART_MODULE(noisychannels)
