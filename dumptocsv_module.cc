////////////////////////////////////////////////////////////////////////
/// \brief Dump supernova-trigger-like data to a CSV file
///
/// \author M. Strait
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "Geometry/Geometry.h"
#include "DAQDataFormats/RawEvent.h"
#include "DAQDataFormats/RawTriggerMask.h"
#include "DAQDataFormats/RawDataBlock.h"
#include "RawData/FlatDAQData.h"
#include "RawData/RawTrigger.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/CellHit.h"
#include "StandardRecord/SREnums.h"
#include "MCCheater/BackTracker.h"
#include "GeometryObjects/CellGeo.h"

// For smearing MC hit timing
#include "TRandom3.h"

#include <vector>

static const int TDC_PER_US = 64;
static const int US_PER_MICROSLICE = 50;

class dumptocsv : public art::EDAnalyzer {
  public:
  explicit dumptocsv(fhicl::ParameterSet const& pset);
  virtual ~dumptocsv() { }; // compiles, but does not run, without this
  void analyze(const art::Event& evt);
};

dumptocsv::dumptocsv(const fhicl::ParameterSet & pset) : EDAnalyzer(pset)
{
}

static int trigger(
  const art::Handle< std::vector<rawdata::RawTrigger> > & rawtrigger)
{
  if(rawtrigger->empty()) return -1;
  return (*rawtrigger)[0].fTriggerMask_TriggerType;
}

static bool goodtriggertype(const int trigger)
{
  // See DAQDataFormats/cxx/include/TriggerDefines.h, and note the
  // off-by-one error between those defines and what appears in files

  // Ignore SNEWS fast beat and LIGO fast beat, because they will
  // overlap with real SNEWS/LIGO triggers, and slow beat SNEWS/LIGO
  // and be in the same files.
  if(trigger+1 == daqdataformats::TRIG_ID_SNEWS_BEAT_FAST) return false;
  if(trigger+1 == daqdataformats::TRIG_ID_LIGO_BEAT_FAST) return false;
  return true;
}

static bool longtriggertype(const int trigger)
{
  if(trigger+1 == daqdataformats::TRIG_ID_SNEWS_BEAT_SLOW) return true;
  if(trigger+1 == daqdataformats::TRIG_ID_LIGO_TRIGGER) return true;
  return false;
}

// Get the length of the event in TDC ticks, and "delta_tdc", the time
// between the trigger time and the time the event starts. You can
// subtract this off of the time that the offline gives each hit to
// get the time since the beginning of the readout, and with the event
// length, the time until the end of the readout.
//
// delta_tdc is a signed 64 bit integer, even though it should always be
// a small positive number, just in case. Ditto for the event length.
//
// Returns whether this information was successfully extracted.
static bool delta_and_length(int64_t & event_length_tdc,
  int64_t & delta_tdc,
  const art::Handle< std::vector<rawdata::FlatDAQData> > & flatdaq,
  const art::Handle< std::vector<rawdata::RawTrigger> > & rawtrigger)
{
  daqdataformats::RawEvent raw;
  if(flatdaq->empty()) return false;

  raw.readData((*flatdaq)[0].getRawBufferPointer());
  if(raw.getDataBlockNumber() == 0) return false;

  raw.setFloatingDataBlock(0);
  daqdataformats::RawDataBlock& datablock = *raw.getFloatingDataBlock();

  uint64_t event_start_time = 0xffffffffffffffff;
  uint64_t event_end_time   = 0x0000000000000000;

  for(unsigned int di = 0; di < raw.getDataBlockNumber(); di++){
    raw.setFloatingDataBlock(di);
    datablock = (*raw.getFloatingDataBlock());

    if(datablock.getHeader()->getMarker() ==
         daqdataformats::datablockheader::SummaryBlock_Marker ||
       !datablock.getHeader()->checkMarker()) continue;

    for(unsigned int mi = 0; mi < datablock.getNumMicroBlocks(); mi++){
      datablock.setFloatingMicroBlock(mi);
      daqdataformats::RawMicroBlock * ub = datablock.getFloatingMicroBlock();

      // The time is always in the second and third words of the
      // microslice, which follows two words of microblock header, so
      // just get it. Justin says you can also get it from getTime(),
      // but this already works and I'm not touching it.
      const uint32_t t_marker_low  = ((uint32_t *)(ub->getBuffer()))[3];
      const uint32_t t_marker_high = ((uint32_t *)(ub->getBuffer()))[4];

      uint64_t time_marker = t_marker_low;
      time_marker |= (uint64_t)t_marker_high << 32;
      if(time_marker < event_start_time) event_start_time = time_marker;
      if(time_marker > event_end_time  ) event_end_time   = time_marker;
    }
  }

  delta_tdc = (int64_t)((*rawtrigger)[0].fTriggerTimingMarker_TimeStart
                        - event_start_time);

  // Assume that microblocks are always 50us; true for all relevant data.
  event_length_tdc = ((int64_t)(event_end_time - event_start_time))
                     + US_PER_MICROSLICE*TDC_PER_US;
  return true; // ok
}

static void getrawtrigger(
  art::Handle< std::vector<rawdata::RawTrigger> > & trg,
  const art::Event & evt)
{
  evt.getByLabel("daq", trg);
  if(trg.failedToGet())
    evt.getByLabel("minbias", trg);
}

static void getrawdigits(
  art::Handle< std::vector<rawdata::RawDigit> > & digits,
  const art::Event & evt)
{
  evt.getByLabel("daq", digits);
  if(digits.failedToGet())
    evt.getByLabel("minbias", digits);
}


// Get the FlatDAQData, either from "minbias", in the case of supernova MC
// overlays, or "daq", for everything else.
static void getflatdaq(
  art::Handle< std::vector<rawdata::FlatDAQData> > & flatdaq,
  const art::Event & evt)
{
  evt.getByLabel("daq", flatdaq);
  if(flatdaq.failedToGet())
    evt.getByLabel("minbias", flatdaq);
}


void dumptocsv::analyze(const art::Event & evt)
{
  // Must be called on every event to prevent SIGPIPE loops.
  signal(SIGPIPE, SIG_DFL);

  art::Handle< std::vector<rawdata::RawTrigger> > rawtrigger;
  getrawtrigger(rawtrigger, evt);
  if(!goodtriggertype(trigger(rawtrigger))) return;

  {
    art::Handle< std::vector<rawdata::FlatDAQData> > flatdaq;
    getflatdaq(flatdaq, evt);

    if(longtriggertype(trigger(rawtrigger))){
      int64_t event_length_tdc = 0, delta_tdc = 0;
      delta_and_length(event_length_tdc, delta_tdc, flatdaq, rawtrigger);

      // Check for empty 50us events in long readouts
      if(event_length_tdc/TDC_PER_US == 50){
        art::Handle< std::vector< rawdata::RawDigit > > rd;
        getrawdigits(rd, evt);
        if(rd->empty()){
          printf("Rejecting trigger: 50us event with zero hits\n");
          return;
        }
      }

      // Check for truly incomplete readouts (not to be confused with events
      // marked incomplete).  Could use these, probably, but would need to
      // modify code to handle overlaps to determine when the overlaps are in
      // these cases, and the complexity doesn't seem worth it.
      const int len = event_length_tdc/TDC_PER_US;
      if(len != 5050 && len != 5000){
        printf("Rejecting trigger of length %dus != 5050 or 5000\n", len);
        return;
      }
    }
  }


  printf("Event %d\n", evt.id().event());

  art::Handle< std::vector< rb::CellHit > > cellhits;
  evt.getByLabel("calhit", cellhits);
  if(cellhits.failedToGet()){
    fprintf(stderr, "Failed to get CellHits.  You need to run calhit\n");
    return;
  }

  art::ServiceHandle<geo::Geometry> geo;
  int whichDet = geo->DetId();

  for(auto & hit : *cellhits){
    float tns = hit.TNS();

    // Smear out MC timing as per my study in doc-45041.
    //
    // Use a tolerable approximation (smear between 17.1ns and 43.4ns,
    // depending on position) that takes into account that the data
    // looks worse compared to the MC for the dim side of FD modules.
    //
    // TODO: I haven't studied the ND at all, but just assume it is the
    // same as the FD as a function of distance to readout.
    if(hit.IsMC()){
      art::ServiceHandle<cheat::BackTracker> bt;
      std::vector<sim::FLSHit> flshits = bt->HitToFLSHit(hit);

      double true_w = 0;
      if(!flshits.empty()){
        const double loc[3] = { flshits[0].GetXAverage(),
                                flshits[0].GetYAverage(),
                                flshits[0].GetZAverage() };
        double wor[3] = { 0 };
        art::ServiceHandle<geo::Geometry> geo;
        geo->Plane(flshits[0].GetPlaneID())
           ->Cell (flshits[0].GetCellID())
           ->LocalToWorld(loc, wor);
        true_w = hit.View() == geo::kX? wor[1]: wor[0];
      }

      const double extrulen = whichDet == caf::kFARDET?1549.4:399.28;

      const double true_d = true_w + extrulen/2;

      static TRandom3 randfortiming;
      tns += randfortiming.Gaus() * (17.1 + true_d*0.017);
    }

    printf("%d %d %.0f %d %d\n", hit.Plane(), hit.Cell(), tns, hit.ADC(), hit.IsMC());
  }
}

DEFINE_ART_MODULE(dumptocsv)
