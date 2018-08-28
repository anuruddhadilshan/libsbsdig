//*-- Author :    Ole Hansen (ole@jlab.org)    9-Dec-2011

/////////////////////////////////////////////////////////////////////
//
//   TSBSSimDecoder
//
//   Decoder for SoLID simulation data
//
//   Interprets event buffer from input as TSBSSimEvent objects
//   (containing digitized simulation data) and unpacks them into
//   crateslot arrays for low-level decoding by detectors.
//
/////////////////////////////////////////////////////////////////////

#include "TSBSSimDecoder.h"
#include "THaCrateMap.h"
#include "THaBenchmark.h"
#include "VarDef.h"
#include "TSBSDBManager.h"
#include "THaSlotData.h"

#include "TError.h"
#include "TSystem.h"
#include "TMath.h"
#include "TDatabasePDG.h"
#include "TRandom.h"
#include "THaVarList.h"
#include "TSBSSimAuxi.h"

#include <SBSSimFadc250Module.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <utility>
#include <stdexcept>

using namespace std;
using namespace Podd;

//EFuchey: 2016/12/10: it is necessary to declare the TSBSDBManager as a static instance here 
// (and not make it a member) because it is used by functions whic are defined as "static inline".
static TSBSDBManager* fManager = TSBSDBManager::GetInstance();
static const Int_t kPrimaryType = 1, kPrimarySource = 0;
// Projection types must match the definitions in TreeSearch
enum EProjType { kUPlane = 0, kVPlane =1, kXPlane = 2, kYPlane = 3};

typedef vector<int>::size_type vsiz_t;

//-----------------------------------------------------------------------------
TSBSSimDecoder::TSBSSimDecoder() : fCheckedForEnabledDetectors(false)
{
  // Constructor
  DefineVariables();

  gSystem->Load("libEG.so");  // for TDatabasePDG
}

//-----------------------------------------------------------------------------
TSBSSimDecoder::~TSBSSimDecoder() {

  DefineVariables( THaAnalysisObject::kDelete );
  
}

//-----------------------------------------------------------------------------
Int_t TSBSSimDecoder::DefineVariables( THaAnalysisObject::EMode mode )
{
  // Define global variables for the MC quantities. Extends the base
  // class method.

  const char* const here = "TSBSSimDecoder::DefineVariables";
  
  if( mode == THaAnalysisObject::kDefine && fIsSetup )
    return THaAnalysisObject::kOK;
  
  SimDecoder::DefineVariables( mode );
  
  cout << "Read TSBSSimDecoder variables " << endl;
  
  RVarDef vars[] = {
    { 0 }
  };

  return THaAnalysisObject::
    DefineVarsFromList( vars, THaAnalysisObject::kRVarDef,
			mode, "", this, MC_PREFIX, here );
}

//-----------------------------------------------------------------------------
void TSBSSimDecoder::Clear( Option_t* opt )
{
  // Clear track and plane data

  SimDecoder::Clear(opt);   // clears fMCCherHits, fMCCherClus
  
  fPMTMap.clear(); 
}

//-----------------------------------------------------------------------------
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
int TSBSSimDecoder::LoadEvent(const UInt_t* evbuffer )
#else
int TSBSSimDecoder::LoadEvent(const Int_t* evbuffer )
#endif
{
  // Wrapper around DoLoadEvent so we can conveniently stop the benchmark
  // counter in case of errors

  int ret = DoLoadEvent( evbuffer );

  if( fDoBench ) fBench->Stop("physics_decode");

  return ret;
}

/*
//-----------------------------------------------------------------------------
static inline
void PMTtoROC( Int_t h_chan,
	       Int_t& crate, Int_t& slot, Int_t& chan )
{
  // Convert location parameters (row, col, chan) of the given PMT
  // to hardware channel (crate,slot,chan)
  // The (crate,slot,chan) assignment must match the detmap definition in
  // the database!  See TreeSearch/dbconvert.cxx
  // In the case of GRINCH/RICH: 
  // crate = GTP; slot = VETROC; chan = PMT. (NINOs are "transparent", in a similar way to the MPDs)
  
  //div_t d = div( h_chan, fManager->GetChanPerSlot() );
  div_t d = div( h_chan, 1 );
  slot = d.quot;
  chan = d.rem;

  d = div( slot, 1 );
  crate = d.quot;
  slot  = d.rem;
}

//-----------------------------------------------------------------------------
static inline
Int_t MakeROCKey( Int_t crate, Int_t slot, Int_t chan )
{
  return chan;// +
  //fManager->GetChanPerSlot()*( slot + fManager->GetSlotPerCrate()*crate );
}

//-----------------------------------------------------------------------------
Int_t TSBSSimDecoder::PMTfromROC( Int_t crate, Int_t slot, Int_t chan ) const
{
  // Return index of digitized strip correspomding to hardware channel
  // (crate,slot,chan)

  if( fPMTMap.empty() )
    return -1;

  PMTMap_t::const_iterator found = fPMTMap.find( MakeROCKey(crate,slot,chan) );
  if( found == fPMTMap.end() )
    return -1;

  return found->second;
}
*/

//-----------------------------------------------------------------------------
static inline Int_t NumberOfSetBits( UInt_t v )
{
  // Count number of bits set in 32-bit integer. From
  // http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel

  v = v - ((v >> 1) & 0x55555555);
  v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
  return (((v + (v >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

//-----------------------------------------------------------------------------
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
Int_t TSBSSimDecoder::DoLoadEvent(const UInt_t* evbuffer )
#else
Int_t TSBSSimDecoder::DoLoadEvent(const Int_t* evbuffer )
#endif
{
  // Fill crateslot structures with Monte Carlo event data in 'evbuffer'

  static const char* const here = "TSBSSimDecoder::LoadEvent";

#if ANALYZER_VERSION_CODE < ANALYZER_VERSION(1,6,0)
  Bool_t fNeedInit = fgNeedInit;
#endif
  assert( fMap || fNeedInit );

  // Local copy of evbuffer pointer, used in GetMCHitInfo
  buffer = evbuffer;

  // Cast the evbuffer pointer back to exactly the event type that is present
  // in the input file (in TSBSSimFile). The pointer-to-unsigned integer is
  // needed compatibility with the standard decoder.
  const TSBSSimEvent* simEvent = reinterpret_cast<const TSBSSimEvent*>(buffer);
  
  Int_t ret = HED_OK;
  if (first_decode || fNeedInit) {
    fMap->print();
    if( (ret = init_cmap()) != HED_OK )
      return ret;
#if ANALYZER_VERSION_CODE < ANALYZER_VERSION(1,6,0)
    if( (ret = init_slotdata(fMap)) != HED_OK)
#else
    if( (ret = init_slotdata()) != HED_OK)
#endif
      return ret;
    first_decode = false;
  }

  if( fDoBench ) fBench->Begin("clearEvent");
  Clear();
  for( int i=0; i<fNSlotClear; i++ )
    crateslot[fSlotClear[i]]->clearEvent();
  if( fDoBench ) fBench->Stop("clearEvent");

  // FIXME: needed?
  evscaler = 0;
  event_length = 0;
  
  event_type = 1;
  event_num = simEvent->fEvtID;
  recent_event = event_num;

  // Event weight
  fWeight = simEvent->fWeight;

  //
  if( fDoBench ) fBench->Begin("physics_decode");
  
  
  Bool_t newclus;
  Int_t crate, slot, chan,lchan;
  
  // We must check at least once which detectors are enabled
  // before we try to load up data for that detector
  if(!fCheckedForEnabledDetectors)
    CheckForEnabledDetectors();

  std::vector<std::map<Decoder::THaSlotData*, std::vector<UInt_t> > > detmaps;
  detmaps.resize(fDetNames.size());

  // looks kinda dumb done this way, but it avoids unnecessary loop on events.
  //std::map<Decoder::THaSlotData*, std::vector<UInt_t> > grinchmap;
  //std::map<Decoder::THaSlotData*, std::vector<UInt_t> > bbpsmap;
  //std::map<Decoder::THaSlotData*, std::vector<UInt_t> > hodomap;
  //std::map<Decoder::THaSlotData*, std::vector<UInt_t> > bbshmap;
  //std::map<Decoder::THaSlotData*, std::vector<UInt_t> > cdetmap;
  //std::map<Decoder::THaSlotData*, std::vector<UInt_t> > hcalmap;
  
  std::cerr << "\n\n\n\n\nStart Processing event: " << event_num << std::endl;
  for(std::vector<TSBSSimEvent::DetectorData>::const_iterator it =
      simEvent->fDetectorData.begin(); it != simEvent->fDetectorData.end();
      ++it )
  {
  /*
    std::cout << "Found detid: " << (*it).fDetID << std::endl;
    if((*it).fDetID == HCAL_UNIQUE_DETID) { // HCal
      TDetInfo detInfo = fManager->GetDetInfo("hcal");
      bool done_with_prev_data = true;
      int buff;
      for( int j = 0; j < (*it).fData.size(); j++ ) {
        int mod =  (*it).fChannel;
        int lchan = mod + (*it).fData[j]*detInfo.NChan();
        TDigChannelInfo chinfo = detInfo.FindLogicalChannelSlot(lchan);
        crate = chinfo.crate;
        slot = chinfo.slot;
        chan = chinfo.ch;
        Decoder::THaSlotData *sldat = crateslot[idx(crate,slot)];
        if(sldat) { // meaning the module is available
          std::vector<UInt_t> *myev = &(hcalmap[sldat]);
          myev->push_back(chan);
          int dat_type = (*it).fData[++j];
          std::cout << "type: ";
          if(dat_type == 0 || dat_type == 1) { // ADC data
            myev->push_back(dat_type); // Push back the type
            std::cout << "adc(" << dat_type << ") ";
          } else {
            std::cout << "tdc ";
          }
          std::cout << "[ ";
          // Push back number of items in data
          int nhits = (*it).fData[++j];
          myev->push_back(nhits);
          for(int k = 0; k < nhits && j < (*it).fData.size(); k++) {
            std::cout << " " << (*it).fData[j+1];
            myev->push_back((*it).fData[++j]);
          }
          std::cout << " ]" << std::endl;
        }
        //if((*it).fData[0] == 1) {
        //  std::cerr << "M: " << mod << ", C: " << crate << ", S: " << slot
        //    << ", C: " << chan << ", I: " << (*it).fData[2] << std::endl;
        //}
      }
    }
    */
    for(size_t d = 0; d < fDetNames.size(); d++) {
      LoadDetector(detmaps[d], fDetNames[d], (*it), fDetIDs[d]);
    }
    //LoadDetector(hcalmap, "hcal", (*it), HCAL_UNIQUE_DETID);
    //LoadDetector(cdetmap, "cdet", (*it), CDET_UNIQUE_DETID);
    //LoadDetector(bbshmap, "sh", (*it), BBSH_UNIQUE_DETID);
    //LoadDetector(hodomap, "hodo", (*it), HODO_UNIQUE_DETID);
    //LoadDetector(bbpsmap, "ps", (*it), BBPS_UNIQUE_DETID);
    //LoadDetector(grinchmap, "grinch", (*it), GRINCH_UNIQUE_DETID);
    
    // what if we were just coding the stuff above in a function ?
    // what would this function need ? name (or CPS/SPC) +detID of det, and map ???? 
    // go for it ?
    // OK, faisons l'exercise bete de copier en adaptant pour e.g. CDet brouillon
    /*
    if((*it).fDetID == CDET_UNIQUE_DETID && (*it).fData.size() > 0) { // 
      int mod =  (*it).fChannel;
      //This should be *general* and work for *every* subsystem
      chan = mod%CPS_cdet;
      slot = ((mod-chan)/CPS_cdet)%SPC_cdet;//+first_slot
      crate = (mod-slot*CPS_cdet-chan)/SPC_cdet;//+first_crate
      
      Decoder::THaSlotData *sldat = crateslot[idx(crate,slot)];
      if(sldat) { // meaning the module is available
	std::vector<UInt_t> *myev = &(map[sldat]);
	myev->push_back(chan);
	for(size_t k = 0; k < (*it).fData.size(); k++) {
	  myev->push_back((*it).fData[k]);
	}
      }
      if((*it).fData[0] == 1) {
	std::cerr << "M: " << mod << ", C: " << crate << ", S: " << slot
		  << ", C: " << chan << ", I: " << (*it).fData[2] << std::endl;
      }
    }
    */
    
  }

  // Now call LoadSlot for the different detectors
  for(size_t d = 0; d < fDetNames.size(); d++) {
    //std::cout << "About to load data for: " << fDetNames[d] << " with ID: " << fDetIDs[d] << std::endl;
    for( std::map<Decoder::THaSlotData*, std::vector<UInt_t> >::iterator it =
        detmaps[d].begin(); it != detmaps[d].end(); ++it) {
      //std::cout << ">>Loading data for: " << fDetNames[d] << " with ID: " << fDetIDs[d] << std::endl;
      std::cout << "Module: " << it->first->GetModule()
        << " c: " << it->first->getCrate()
        << ", s: " << it->first->getSlot() << std::endl;
      it->first->GetModule()->LoadSlot(it->first,
          it->second.data(),0,it->second.size() );
    }
  }

  /*
  
  // CDET
  for( std::map<Decoder::THaSlotData*, std::vector<UInt_t> >::iterator it =
	 cdetmap.begin(); it != hcalmap.end(); ++it) {
    it->first->GetModule()->LoadSlot(it->first,
        it->second.data(),0,it->second.size() );
  }
  // BBSH
  for( std::map<Decoder::THaSlotData*, std::vector<UInt_t> >::iterator it =
	 bbshmap.begin(); it != bbshmap.end(); ++it) {
    it->first->GetModule()->LoadSlot(it->first,
        it->second.data(),0,it->second.size() );
  }
  // HODO
  for( std::map<Decoder::THaSlotData*, std::vector<UInt_t> >::iterator it =
	 hodomap.begin(); it != hodomap.end(); ++it) {
    it->first->GetModule()->LoadSlot(it->first,
        it->second.data(),0,it->second.size() );
  }
  // BBPS
  for( std::map<Decoder::THaSlotData*, std::vector<UInt_t> >::iterator it =
	 bbpsmap.begin(); it != bbpsmap.end(); ++it) {
    it->first->GetModule()->LoadSlot(it->first,
        it->second.data(),0,it->second.size() );
  }
  // GRINCH
  for( std::map<Decoder::THaSlotData*, std::vector<UInt_t> >::iterator it =
	 grinchmap.begin(); it != grinchmap.end(); ++it) {
    it->first->GetModule()->LoadSlot(it->first,
        it->second.data(),0,it->second.size() );
  }
  */
  
  std::cerr << "End Processing event:   " << event_num << std::endl;
  return HED_OK;
}

Int_t TSBSSimDecoder::RetrieveDetMapParam(const char* detname, 
					  int& chanperslot, int& slotpercrate, 
					  int& firstcrate, int& firstslot)
{
  // chanperslot = ((TDetInfo &)fManager->GetDetInfo("hcal")).ChanPerSlot();
  // slotpercrate = ((TDetInfo &)fManager->GetDetInfo("hcal")).SlotPerCrate();
  // firstslot = ((TDetInfo &)fManager->GetDetInfo("hcal")).FirstSlot();
  // firstcrate = ((TDetInfo &)fManager->GetDetInfo("hcal")).FirstCrate();
  TDetInfo detinfo = fManager->GetDetInfo(detname);
  chanperslot = detinfo.ChanPerSlot();
  slotpercrate = detinfo.SlotPerCrate();
  firstslot = detinfo.FirstSlot();
  firstcrate = detinfo.FirstCrate();
}

Int_t TSBSSimDecoder::LoadDetector( std::map<Decoder::THaSlotData*, std::vector<UInt_t> > &map,
      const char *detname, TSBSSimEvent::DetectorData detdata, const int detid)
{
  TDetInfo detinfo = fManager->GetDetInfo(detname);

//Int_t TSBSSimDecoder::LoadDetector( std::map<Decoder::THaSlotData*, std::vector<UInt_t> > map,
//				    TSBSSimEvent::DetectorData detdata, 
//				    const int detid, 
//				    const int chanperslot, const int slotpercrate, 
//				    const int firstcrate, const int firstslot)
//{
  Int_t crate, slot, chan;
  TDetInfo detInfo = fManager->GetDetInfo(detname);

  if(detdata.fDetID == detid && detdata.fData.size() > 0) { // Data to process
    int mod =  detdata.fChannel;
    //This should be *general* and work for *every* subsystem
    // Loop over all raw data in this event
    for( int j = 0; j < detdata.fData.size(); j++ ) {
      // Identify the "logical" channel number for this event
      // based on the first integer in the raw data
      int lchan = mod + detdata.fData[j]*detInfo.NChan();
      // Get information about this logical channel from TDetInfo
      TDigChannelInfo chinfo = detInfo.FindLogicalChannelSlot(lchan);
      crate = chinfo.crate;
      slot = chinfo.slot;
      chan = chinfo.ch;
      Decoder::THaSlotData *sldat = 0;
      if( crate >= 0 || slot >=  0 || chan >= 0) {
        sldat = crateslot[idx(crate,slot)];
      }
      // Now get the corresponding THaSlotData based on crate and slot
      // and load it with the data
      if(sldat) { // If module available, we are free to store data in it
        std::vector<UInt_t> *myev = &(map[sldat]);
        myev->push_back(chan);
        for(size_t k = 0; k < detdata.fData.size(); k++) {
          myev->push_back(detdata.fData[k]);
          j++;
        }
      } else {
        std::cerr << "Yikes!! No data for " << detname << " in c: "
          << crate << " s: " << slot << " c: " << chan << std::endl;
      }
    }
  }

  return HED_OK;
}


void TSBSSimDecoder::CheckForEnabledDetectors()
{
  fDetNames.clear();
  fDetIDs.clear();
  CheckForDetector("hcal",HCAL_UNIQUE_DETID);
  CheckForDetector("cdet", CDET_UNIQUE_DETID);
  CheckForDetector("sh", BBSH_UNIQUE_DETID);
  CheckForDetector("hodo", HODO_UNIQUE_DETID);
  CheckForDetector("ps", BBPS_UNIQUE_DETID);
  CheckForDetector("grinch", GRINCH_UNIQUE_DETID);
  fCheckedForEnabledDetectors = true;
}

void TSBSSimDecoder::CheckForDetector(const char *detname, short id)
{
  if(fManager->IsDetInfoAvailable(detname)) {
    fDetNames.push_back(detname);
    fDetIDs.push_back(id);
  }
}
