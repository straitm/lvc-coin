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

#include "func/timeutil.h"

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
  if(!goodtriggertype(trigger(rawtrigger))) return;

  printf("Spilltime: %f\n", art_time_to_unix_double(evt.time().value()));

}

DEFINE_ART_MODULE(spilltime)
