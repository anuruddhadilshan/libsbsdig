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
      fadc_data[i].integral = 0;
      fadc_data[i].samples.clear();
    }
  }

  void TSBSSimADC::Clear( const Option_t* opt) {
    // Clear event-by-event data
    ClearDataVectors();
  }

  void TSBSSimADC::Init() {
    Clear();
    IsInit = kTRUE;
    fName = "SBSSimADC (Simple JLab Flash ADC Simulated Module)";
  }

  void TSBSSimADC::CheckDecoderStatus() const {
  }

  Int_t TSBSSimADC::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer,
      const UInt_t *pstop) {
    Clear();
    unsigned int nwords = 0;
    unsigned short chan = 0, type = 0;
    UInt_t raw_buff;
    bool printed = false;
    while(evbuffer < pstop) {
      // First, decode the header
      TSBSSimDataEncoder::DecodeHeader(*evbuffer++,type,chan,nwords);
      TSBSSimDataEncoder *enc = TSBSSimDataEncoder::GetEncoder(type);
      if(enc && nwords > 0) {
        if(enc->IsFADC()) { // FADC with samples
          SimEncoder::fadc_data tmp_fadc_data;
          enc->DecodeFADC(tmp_fadc_data,evbuffer,nwords);
          evbuffer += nwords; // skip ahead the total number of words read
          for(size_t i = 0; i < tmp_fadc_data.samples.size(); i++) {
            raw_buff = tmp_fadc_data.samples[i];
            fadc_data[chan].samples.push_back(tmp_fadc_data.samples[i]);
            sldat->loadData("adc",chan,raw_buff,raw_buff);
            //std::cout << " " << raw_buff;
            //printed = true;
          }
        } /*else if (type==1) { // integral of adc
            num_samples = *evbuffer++;
            for(int i = 0; i < num_samples; i++) {
            raw_buff = *evbuffer++;
        //std::cerr << " [" << chan << ", " << raw_buff << "]";
        //printed = true;
        fadc_data[chan].integrals.push_back(raw_buff);
        sldat->loadData("adc",chan,raw_buff,raw_buff);
        }
        } */
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
