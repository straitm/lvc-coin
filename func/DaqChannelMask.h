////////////////////////////////////////////////////////////////////////
// \file:   DaqChannelMask.h
// \brief   Book-keeping for channels that need to be masked from an event
// \author  Justin Vasel <justin.vasel@gmail.com>
// \date    2019-12-05
////////////////////////////////////////////////////////////////////////

#ifndef DAQCHANNELMASK_H
#define DAQCHANNELMASK_H

#include <cstdint>

#include "RawData/RawDigit.h"
#include "RecoBase/CellHit.h"

namespace sn {
  class DaqChannelMask {
  public:
    DaqChannelMask(float coldThreshold, float hotThreshold);
    ~DaqChannelMask();

    void AddHit(rawdata::RawDigit d);
    void AddHit(rb::CellHit h);

    void IncrementDuration(float i);
    float CurrentDuration();

    void CalculateRates();
    bool RatesCalculated();

    float Rate(rawdata::RawDigit d);
    float Rate(rb::CellHit h);

    bool ChannelIsMasked(rawdata::RawDigit d);
    bool ChannelIsMasked(rb::CellHit h);

    size_t MaskSize();
    size_t RatesSize();

    std::map<unsigned int, float> GetRates();

    void Print();

  private:
    std::set<unsigned int> fChannelMask;
    std::map<unsigned int, float> fChannelRates;

    bool fRatesCalculated = false;
    float fIntegratedDuration = 0;

    int fNumColdChannels = 0;
    int fNumHotChannels = 0;

    float fColdThreshold;
    float fHotThreshold;

    float fHitsAdded = 0;
  };
}

#endif
