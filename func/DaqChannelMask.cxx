////////////////////////////////////////////////////////////////////////
// \file:   DaqChannelMask.cxx
// \brief   Book-keeping for channels that need to be masked from an event
// \author  Justin Vasel <justin.vasel@gmail.com>
// \date    2019-12-05
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iterator>
#include <set>

#include "RawData/RawDigit.h"
#include "RecoBase/CellHit.h"
#include "Supernova/BackgroundRejection/DaqChannelMask.h"


// ............................................................................
sn::DaqChannelMask::DaqChannelMask(float coldThreshold, float hotThreshold):
fColdThreshold(coldThreshold),
fHotThreshold(hotThreshold)
{

}

// ............................................................................
sn::DaqChannelMask::~DaqChannelMask() {}

// ............................................................................
void sn::DaqChannelMask::AddHit(rawdata::RawDigit d)
{
  if (this->fRatesCalculated) {
    std::cout << "Channel rates have already been computed, so the mask is now immutable." << std::endl;
  }

  std::map<unsigned int, float>::iterator it = this->fChannelRates.find(d.DaqChannel());
  if (it != this->fChannelRates.end()) {
    it->second += 1.0;
  } else {
    this->fChannelRates.insert(std::pair<unsigned int, float>(d.DaqChannel(), 1.0));
  }

  ++this->fHitsAdded;

  return;
}

// ............................................................................
void sn::DaqChannelMask::AddHit(rb::CellHit h)
{
  if (this->fRatesCalculated) {
    std::cout << "Channel rates have already been computed, so the mask is now immutable." << std::endl;
  }

  std::map<unsigned int, float>::iterator it = this->fChannelRates.find(h.DaqChannel());
  if (it != this->fChannelRates.end()) {
    it->second += 1.0;
  } else {
    this->fChannelRates.insert(std::pair<unsigned int, float>(h.DaqChannel(), 1.0));
  }

  return;
}

// ............................................................................
void sn::DaqChannelMask::IncrementDuration(float i)
{
  this->fIntegratedDuration += i;
}

// ............................................................................
float sn::DaqChannelMask::CurrentDuration()
{
  return this->fIntegratedDuration;
}

// ............................................................................
void sn::DaqChannelMask::CalculateRates()
{
  for (std::map<unsigned int, float>::iterator channel = this->fChannelRates.begin(); channel != this->fChannelRates.end(); ++channel) {
    // compute rate
    channel->second /= (this->fIntegratedDuration / 5e-3);

    // determine if channel needs to be masked
    if (channel->second < this->fColdThreshold) {
      this->fChannelMask.insert(channel->first);
      ++this->fNumColdChannels;
    }
    if (channel->second > this->fHotThreshold) {
      this->fChannelMask.insert(channel->first);
      ++this->fNumHotChannels;
    }
  }

  this->fRatesCalculated = true;

  return;
}

// ............................................................................
bool sn::DaqChannelMask::RatesCalculated()
{
  return this->fRatesCalculated == true;
}


// ............................................................................
float sn::DaqChannelMask::Rate(rawdata::RawDigit d)
{
  return this->fChannelRates[d.DaqChannel()];
}


// ............................................................................
float sn::DaqChannelMask::Rate(rb::CellHit h)
{
  return this->fChannelRates[h.DaqChannel()];
}


// ............................................................................
bool sn::DaqChannelMask::ChannelIsMasked(rawdata::RawDigit d)
{
  std::set<unsigned int>::iterator loc = this->fChannelMask.find(d.DaqChannel());
  return loc != this->fChannelMask.end();
}


// ............................................................................
bool sn::DaqChannelMask::ChannelIsMasked(rb::CellHit h)
{
  std::set<unsigned int>::iterator loc = this->fChannelMask.find(h.DaqChannel());
  return loc != this->fChannelMask.end();
}


// ............................................................................
size_t sn::DaqChannelMask::MaskSize()
{
  return this->fChannelMask.size();
}


// ............................................................................
size_t sn::DaqChannelMask::RatesSize()
{
  return this->fChannelRates.size();
}


// ............................................................................
std::map<unsigned int, float> sn::DaqChannelMask::GetRates()
{
  return this->fChannelRates;
}


// ............................................................................
void sn::DaqChannelMask::Print()
{
  std::cout << "COLD/HOT CHANNEL MASK DETAILS" << std::endl;
  std::cout << "-----------------------------" << std::endl;
  std::cout << "Cold Channels: " << this->fNumColdChannels << std::endl;
  std::cout << "Hot Channels:  " << this->fNumHotChannels  << std::endl;
  std::cout << "-----------------------------" << std::endl;
}
