//////////////////////////////////////////////////////////////////
//
//   TSBSSimTDC

#include "TSBSSimTDC.h"
#include "TSBSSimDataEncoder.h"
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
      tdc_data[i].lead_time.clear();
      tdc_data[i].trail_time.clear();
    }
  }

  void TSBSSimTDC::Clear( const Option_t* opt) {
    // Clear event-by-event data
    ClearDataVectors();
  }

  void TSBSSimTDC::Init() {
    Clear();
    IsInit = kTRUE;
    fName = "SBSSimTDC Generic TDC Module";
  }

  void TSBSSimTDC::CheckDecoderStatus() const {
  }

  Int_t TSBSSimTDC::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer,
      const UInt_t *pstop) {
    Clear();
    unsigned int nwords = 0;
    unsigned short chan = 0, type = 0;
    UInt_t raw_buff;
    SimEncoder::tdc_data tmp_tdc_data;
    while(evbuffer < pstop) {
      // First, decode the header
      chan = type = nwords = 0;
      TSBSSimDataEncoder::DecodeHeader(*evbuffer++,type,chan,nwords);
      TSBSSimDataEncoder *enc = TSBSSimDataEncoder::GetEncoder(type);
      evbuffer += nwords; // Skip ahead in the buffer
      if(!enc) {
        std::cerr << "Could not find TDC decoder of type: " << type
          << std::endl;
      } else {
        if(!enc->IsTDC()) {
          std::cerr << "Encoder " << enc->GetName() << " of type " << type
            << " is not a TDC!" << std::endl;
        } else if ( nwords > 0 ) {
          enc->DecodeTDC(tmp_tdc_data,evbuffer,nwords);
          std::cerr << "Got TDC encoder for type: " << type
            << ", name: " << enc->GetName() << std::endl;
          for(size_t i = 0; i < tmp_tdc_data.time.size(); i++ ) {
            raw_buff = tmp_tdc_data.getTime(i);
            if(tmp_tdc_data.getEdge(i)) { // Trail
              tdc_data[chan].lead_time.push_back(raw_buff);
            } else { // Lead
              tdc_data[chan].trail_time.push_back(raw_buff);
            }
            // TODO: Figure out what to do with the edge information
            // I'd imagine we need to distinguish it somehow!
            sldat->loadData("tdc",chan,raw_buff,raw_buff);
          }
          tmp_tdc_data.time.clear(); // Clear it to prepare for next read
        }
      }
    }
   return 0;
  }

  Int_t TSBSSimTDC::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer,
      Int_t pos, Int_t len) {
    return LoadSlot(sldat,evbuffer+pos,evbuffer+pos+len);
    //return TSBSSimTDC::LoadSlot(sldat,evbuffer,len);
  }

}

ClassImp(Decoder::TSBSSimTDC)
