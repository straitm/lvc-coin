////////////////////////////////////////////////////////////////////////
/// \brief   A new version of RemoveBeamSpills that takes a list of
//           NuMI event times as input, and does not depend on the
//           mysterious database that RemoveBeamSpills does.
/// \author  M. Strait
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include "DAQDataFormats/RawEvent.h"
#include "DAQDataFormats/RawDataBlock.h"
#include "RawData/FlatDAQData.h"
#include "RawData/RawTrigger.h"

#include <string>
#include <algorithm>
#include <fstream>
#include <signal.h>

const int TDC_PER_US = 64;

class EliminateBeamSpills : public art::EDFilter {
  public:
  explicit EliminateBeamSpills(const fhicl::ParameterSet & pset);
  virtual ~EliminateBeamSpills() { }; // compiles, but does not run, without this
  bool filter(art::Event& evt);
  void endJob();
  void beginJob();

  private:

  // If the trigger time of an event is before the beginning spill by
  // less than this amount, in seconds, it will be removed by the filter.
  //
  // You probably want to set this to at least the length of your events
  // 50us to 5005us, depending on the trigger plus 50us since the trigger
  // time could be anywhere in the first 50us.  Additional buffer is probably
  // not needed since no activity should preceed the spill.
  double bufferbefore;

  // If the trigger time of an event is after the *beginning* of a
  // spill (not the end, because it's not trivial to detminine the
  // length of each spill) by less than this amount, in seconds, it
  // will be removed by the filter.
  //
  // You probably want at least 50us since the trigger time of your
  // event could be anywhere in the first 50us of readout, and for many
  // purposes you probably want more. For instance, another 5-10us is
  // needed for FEB flashers to die down, 10us for Michels to decay
  // away, several 100us for most neutron captures to end, and several
  // milliseconds for a long tail of neutrons diffusing out of the rock
  // and drifting through the air.
  double bufferafter;

  // Name of the file containing the list of spill times.  This file should
  // contain three columns.  The first is ignored, the second is the NOvA time
  // of the spill (64ths of microseconds since the NOvA epoch).
  std::string spillfile;

  std::vector<uint64_t> spilltimes;
};

EliminateBeamSpills::EliminateBeamSpills(fhicl::ParameterSet const& pset) : EDFilter(),
  bufferbefore(pset.get<double>("bufferbefore")),
  bufferafter(pset.get<double>("bufferafter")),
  spillfile(pset.get<std::string>("spillfile"))
{
}

void EliminateBeamSpills::endJob()
{
}

void EliminateBeamSpills::beginJob()
{
  std::ifstream in(spillfile.c_str());
  if(!in.is_open()){
    fprintf(stderr, "Could not open spill timestamp file %s\n", spillfile.c_str());
    exit(1);
  }

  std::string dum;
  uint64_t t;

  while(in >> dum >> t) spilltimes.push_back(t);
  std::sort(spilltimes.begin(), spilltimes.end());
}

static bool start_and_end(uint64_t & event_start_time,
  uint64_t & event_end_time,
  const art::Handle< std::vector<rawdata::FlatDAQData> > & flatdaq,
  const art::Handle< std::vector<rawdata::RawTrigger> > & rawtrigger)
{
  daqdataformats::RawEvent raw;
  if(flatdaq->empty()) return false;

  raw.readData((*flatdaq)[0].getRawBufferPointer());
  if(raw.getDataBlockNumber() == 0) return false;

  raw.setFloatingDataBlock(0);
  daqdataformats::RawDataBlock& datablock = *raw.getFloatingDataBlock();

  event_start_time = 0xffffffffffffffff;
  event_end_time   = 0x0000000000000000;

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

  // Assume that microblocks are always 50us. I hope that's true for all
  // relevant data.
  event_end_time += 50 * TDC_PER_US;
  return true; // ok
}

static void getrawtrigger(
  art::Handle< std::vector<rawdata::RawTrigger> > & trg,
  const art::Event & evt)
{
  evt.getByLabel("minbias", trg);
  if(trg.failedToGet()) evt.getByLabel("daq", trg);
}


// Returns true if the event is to be kept, false if it is too close to a spill
bool EliminateBeamSpills::filter(art::Event & evt)
{
  // Must be called on every event to prevent SIGPIPE loops.
  signal(SIGPIPE, SIG_DFL);

  art::Handle< std::vector<rawdata::FlatDAQData> > flatdaq;
  evt.getByLabel("daq", flatdaq);

  art::Handle< std::vector<rawdata::RawTrigger> > rawtrigger;
  getrawtrigger(rawtrigger, evt);
  if(rawtrigger->empty()) return false;

  uint64_t start, end;

  if(!start_and_end(start, end, flatdaq, rawtrigger)){
    fprintf(stderr, "Could not get event start and end\n");
    return false;
  }

  // Take all spills to be 10us long.  Could be wrong depending on what AD
  // is doing, but shouldn't be off by more than a few us.
  const double spill_len_us = 10;
  const uint64_t spill_len_tdc = spill_len_us * TDC_PER_US;

  const uint64_t startbuf = start - (uint64_t)(bufferbefore * TDC_PER_US * 1e6);
  const uint64_t endbuf   = end   + (uint64_t)(bufferafter * TDC_PER_US * 1e6);

  // Find the next spill after the end of the event and its buffer
  auto nextspill = upper_bound(spilltimes.begin(), spilltimes.end(), endbuf);

  // If all spills are after the event, then we're good.  (Although in this case,
  // the user perhaps should have provided a longer list of spills.)
  if(nextspill == spilltimes.begin()) return true;

  auto previousspill = nextspill-1;

  return *previousspill + spill_len_tdc <= startbuf;
}
DEFINE_ART_MODULE(EliminateBeamSpills)
