//////////////////////////////////////////////////////////////////
//
//   TSBSSimADC

#include "TSBSSimADC.h"
#include "THaEvData.h"
#include "THaSlotData.h"
#include "TMath.h"

#include <unistd.h>
#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>  // for memset
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <cassert>

using namespace std;

//#define DEBUG
//#define WITH_DEBUG

namespace Decoder {

  Module::TypeIter_t TSBSSimADC::fgThisType =
    DoRegister( ModuleType( "Decoder::TSBSSimADC" , 50250 ));

  TSBSSimADC::TSBSSimADC()
  {
    fadc_data.resize(NADCCHAN);
  }

  TSBSSimADC::TSBSSimADC(Int_t crate, Int_t slot)
  {
    fadc_data.resize(NADCCHAN);
    IsInit = kFALSE;
    Init();
  }

  TSBSSimADC::~TSBSSimADC() {
#if defined DEBUG && defined WITH_DEBUG
    // delete fDebugFile; fDebugFile = 0;
#endif
  }

  /*
  Bool_t TSBSSimADC::HasCapability(Decoder::EModuleType type) {
    return ( type == kSampleADC || type == kPulseIntegral || type == kPulseTime
        || type == kPulsePeak || type == kPulsePedestal || type == kCoarseTime || type == kFineTime);
  } */

  // Clear all data vectors
  void TSBSSimADC::ClearDataVectors() {
    // Clear all data objects
    assert(fadc_data.size() == NADCCHAN);  // Initialization error in constructor
    for (uint32_t i = 0; i < NADCCHAN; i++) {
      fadc_data[i].clear();
    }
  }

  void TSBSSimADC::Clear( const Option_t* opt) {
    // Clear event-by-event data
    ClearDataVectors();
  }

  void TSBSSimADC::Init() {
    Clear();
    IsInit = kTRUE;
    fName = "SBSSimFADC250 JLab Flash ADC Module";
  }

  void TSBSSimADC::CheckDecoderStatus() const {
    cout << "FADC250 Decoder has been called" << endl;
  }

  Int_t TSBSSimADC::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer,
      const UInt_t *pstop) {
    Clear();
    int chan = 0, type = 0, num_samples = 0;
    UInt_t raw_buff;
    bool printed = false;
    while(evbuffer < pstop) {
      // First get channel number
      chan = *evbuffer++;
      type = *evbuffer++;
      if(type == 0) { // Samples mode
        num_samples = *evbuffer++;
        for(int i = 0; i < num_samples; i++) {
          raw_buff = *evbuffer++;
          fadc_data[chan].samples.push_back(raw_buff);
          sldat->loadData("adc",chan,raw_buff,raw_buff);
        }
      } else if (type==1) { // integral of adc
        num_samples = *evbuffer++;
        for(int i = 0; i < num_samples; i++) {
          raw_buff = *evbuffer++;
          std::cerr << " [" << chan << ", " << raw_buff << "]";
          printed = true;
          fadc_data[chan].integrals.push_back(raw_buff);
          sldat->loadData("adc",chan,raw_buff,raw_buff);
        }
      }
    }
    if(printed)
      std::cerr << std::endl;
   return 0;
  }

  Int_t TSBSSimADC::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer,
      Int_t pos, Int_t len) {
    return LoadSlot(sldat,evbuffer+pos,evbuffer+pos+len);
    //return TSBSSimADC::LoadSlot(sldat,evbuffer,len);
  }

}

ClassImp(Decoder::TSBSSimADC)
