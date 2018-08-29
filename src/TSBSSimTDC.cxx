//////////////////////////////////////////////////////////////////
//
//   TSBSSimTDC

#include "TSBSSimTDC.h"
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

  Module::TypeIter_t TSBSSimTDC::fgThisType =
    DoRegister( ModuleType( "Decoder::TSBSSimTDC" , 53204 ));

  TSBSSimTDC::TSBSSimTDC()
  {
    tdc_data.resize(NADCCHAN);
  }

  TSBSSimTDC::TSBSSimTDC(Int_t crate, Int_t slot)
  {
    tdc_data.resize(NADCCHAN);
    IsInit = kFALSE;
    Init();
  }

  TSBSSimTDC::~TSBSSimTDC() {
#if defined DEBUG && defined WITH_DEBUG
    // delete fDebugFile; fDebugFile = 0;
#endif
  }

  /*
  Bool_t TSBSSimTDC::HasCapability(Decoder::EModuleType type) {
    return ( type == kSampleADC || type == kPulseIntegral || type == kPulseTime
        || type == kPulsePeak || type == kPulsePedestal || type == kCoarseTime || type == kFineTime);
  } */

  // Clear all data vectors
  void TSBSSimTDC::ClearDataVectors() {
    // Clear all data objects
    assert(tdc_data.size() == NADCCHAN);  // Initialization error in constructor
    for (uint32_t i = 0; i < NADCCHAN; i++) {
      tdc_data[i].lead_time = 0;
      tdc_data[i].trail_time = 0;
    }
  }

  void TSBSSimTDC::Clear( const Option_t* opt) {
    // Clear event-by-event data
    ClearDataVectors();
  }

  void TSBSSimTDC::Init() {
    Clear();
    IsInit = kTRUE;
    fName = "SBSSimFADC250 JLab Flash ADC Module";
  }

  void TSBSSimTDC::CheckDecoderStatus() const {
    cout << "FADC250 Decoder has been called" << endl;
  }

  Int_t TSBSSimTDC::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer,
      const UInt_t *pstop) {
    Clear();
    int chan = 0, type = 0;//, num_samples = 0; //not needed so far
    UInt_t raw_buff;
    bool printed = false;
    while(evbuffer < pstop) {
      // First get channel number
      chan = *evbuffer++;
      type = *evbuffer++;
      if(type == 0) { // Leading Edge
        raw_buff = *evbuffer++;
        tdc_data[chan].lead_time = raw_buff;
        sldat->loadData("tdc",chan,raw_buff,raw_buff);
      } else if (type==1) { // Trailing edge
        raw_buff = *evbuffer++;
        tdc_data[chan].trail_time = raw_buff;
        sldat->loadData("tdc",chan,raw_buff,raw_buff);
      }
    }
    if(printed)
      std::cerr << std::endl;
   return 0;
  }

  Int_t TSBSSimTDC::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer,
      Int_t pos, Int_t len) {
    return LoadSlot(sldat,evbuffer+pos,evbuffer+pos+len);
    //return TSBSSimTDC::LoadSlot(sldat,evbuffer,len);
  }

}

ClassImp(Decoder::TSBSSimTDC)
