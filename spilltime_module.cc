////////////////////////////////////////////////////////////////////////
/// \brief Print the time of each NuMI spill
///
/// \author M. Strait
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "RawData/RawTrigger.h"
#include "DAQDataFormats/TriggerDefines.h"

#include <signal.h>

class spilltime : public art::EDAnalyzer {
  public:
  explicit spilltime(fhicl::ParameterSet const& pset);
  virtual ~spilltime() { }; // compiles, but does not run, without this
  void analyze(const art::Event& evt);

  void beginJob() {};
};

spilltime::spilltime(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset)
{
}

static bool goodtriggertype(const int trigger)
{
  // See DAQDataFormats/cxx/include/TriggerDefines.h, and note the
  // off-by-one error between those defines and what appears in files
  return trigger+1 == daqdataformats::TRIG_ID_BEAM_NUMI;
}

static void getrawtrigger(
  art::Handle< std::vector<rawdata::RawTrigger> > & trg,
  const art::Event & evt)
{
  evt.getByLabel("minbias", trg);
  if(trg.failedToGet()) evt.getByLabel("daq", trg);
}

static int trigger(
  const art::Handle< std::vector<rawdata::RawTrigger> > & rawtrigger)
{
  if(rawtrigger->empty()) return -1;
  return (*rawtrigger)[0].fTriggerMask_TriggerType;
}

void spilltime::analyze(const art::Event & evt)
{
  // Must be called on every event to prevent SIGPIPE loops.
  signal(SIGPIPE, SIG_DFL);

  art::Handle< std::vector<rawdata::RawTrigger> > rawtrigger;
  getrawtrigger(rawtrigger, evt);
  if(rawtrigger->empty()) return;

  if(!goodtriggertype(trigger(rawtrigger))) return;

  // doc-34154, good to ~0.05us
  const double trigger2spill_s =  217.85e-6;
  const double tdc_per_second = 64e6;
  const uint64_t trigger2spills_tdc = uint64_t(trigger2spill_s * tdc_per_second);

  // This is a monotonically increasing timestamp with no worries about leap
  // seconds.  The only leap second trouble is that since AD systems don't handle
  // leap seconds well, we might be taking triggers that are 1s (or other amounts)
  // off of the physical spills near leap seconds.
  const uint64_t triggertime = (*rawtrigger)[0].fTriggerTimingMarker_TimeStart;
  
  const uint64_t spilltime = triggertime + trigger2spills_tdc;

  printf("Spilltime: %lu\n", spilltime);

}

DEFINE_ART_MODULE(spilltime)
