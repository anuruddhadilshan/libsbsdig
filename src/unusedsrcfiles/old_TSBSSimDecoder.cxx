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
TSBSSimDecoder::TSBSSimDecoder()
{
  // Constructor

  fMCCherHits = new TClonesArray( "TSBSSimPMTHit",      2000 );
  fMCCherClus = new TClonesArray( "TSBSSimCherCluster", 500 );
  
  DefineVariables();

  gSystem->Load("libEG.so");  // for TDatabasePDG
}

//-----------------------------------------------------------------------------
TSBSSimDecoder::~TSBSSimDecoder() {

  DefineVariables( THaAnalysisObject::kDelete );
  
  delete fMCCherHits;
  delete fMCCherClus;
  
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
    // PMT hits
    { "cer.pmt.id",      "MC Cher PMT ID",             "fMCCherHits.TSBSSimPMTHit.fID"          },
    { "cer.pmt.src",     "MC Cher PMT Source",         "fMCCherHits.TSBSSimPMTHit.fSource"      },
    { "cer.pmt.type",    "MC Cher PMT Type",           "fMCCherHits.TSBSSimPMTHit.fType"        },
    { "cer.pmt.MCID",    "MC Cher PMT MC track ID",    "fMCCherHits.TSBSSimPMTHit.fMCtrackID"   },
    { "cer.pmt.MCPID",   "MC Cher PMT MC track PID",   "fMCCherHits.TSBSSimPMTHit.fMCtrackPID"  },
    { "cer.pmt.origVol", "MC Cher PMT origin volume",  "fMCCherHits.TSBSSimPMTHit.fOrigVolFlag" },
    { "cer.pmt.x",       "MC Cher PMT X (m)",          "fMCCherHits.TSBSSimPMTHit.fXPMT"        },
    { "cer.pmt.y",       "MC Cher PMT Y (m)",          "fMCCherHits.TSBSSimPMTHit.fYPMT"        },
    { "cer.pmt.Npe",     "MC Cher PMT N pe",           "fMCCherHits.TSBSSimPMTHit.fNpe"         },
    { "cer.pmt.t0",      "MC Cher PMT t0 (ns)",        "fMCCherHits.TSBSSimPMTHit.fTime"        },
    { "cer.pmt.tr",      "MC Cher PMT rise time (ns)", "fMCCherHits.TSBSSimPMTHit.fTDCtime1"    },
    { "cer.pmt.tf",      "MC Cher PMT fall time (ns)", "fMCCherHits.TSBSSimPMTHit.fTDCtime2"    },
    { "cer.pmt.detID",   "MC Cher PMT detector ID",    "fMCCherHits.TSBSSimPMTHit.fDetID"       },
    { "cer.pmt.chan",    "MC Cher PMT chan",           "fMCCherHits.TSBSSimPMTHit.fChannel"     },
    { "cer.pmt.row",     "MC Cher PMT row",            "fMCCherHits.TSBSSimPMTHit.fPMTrow"      },
    { "cer.pmt.col",     "MC Cher PMT col",            "fMCCherHits.TSBSSimPMTHit.fPMTcol"      },
    
    // PMT clusters
    { "cer.cl.size",    "MC Cher clus size",                "fMCCherClus.TSBSSimCherCluster.fSize"            },
    { "cer.cl.X",       "MC Cher clus avg X (m)",           "fMCCherClus.TSBSSimCherCluster.fXcenter"         },
    { "cer.cl.Y",       "MC Cher clus avg Y (m)",           "fMCCherClus.TSBSSimCherCluster.fYcenter"         },
    { "cer.cl.X_w",     "MC Cher clus avg X (w Npe) (m)",   "fMCCherClus.TSBSSimCherCluster.fXcenter_w"       },
    { "cer.cl.Y_w",     "MC Cher clus avg X (w Npe) (m)",   "fMCCherClus.TSBSSimCherCluster.fYcenter_w"       },
    { "cer.cl.Npe",     "MC Cher clus Npe total",           "fMCCherClus.TSBSSimCherCluster.fNpe"             },
    { "cer.cl.tr_avg",  "MC Cher clus avg rise time (ns)",  "fMCCherClus.TSBSSimCherCluster.fMeanRisingTime"  },
    { "cer.cl.tr_rms",  "MC Cher clus rise time rms (ns)",  "fMCCherClus.TSBSSimCherCluster.fRisingTimeRMS"   },
    { "cer.cl.tf_avg",  "MC Cher clus avg fall time (ns)",  "fMCCherClus.TSBSSimCherCluster.fMeanFallingTime" },
    { "cer.cl.tf_rms",  "MC Cher clus fall time rms (ns)",  "fMCCherClus.TSBSSimCherCluster.fFallingTimeRMS"  },
    { "cer.cl.MCID",    "MC Cher clus track G4ID",          "fMCCherClus.TSBSSimCherCluster.fMCtrackID"       },
    { "cer.cl.MCPID",   "MC Cher clus track G4PID",         "fMCCherClus.TSBSSimCherCluster.fMCtrackPID"      },
    /*
      if something is to be added, can be found easily in libsbsgem::TSBSSimDecoder
    */
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
  
  fMCCherHits->Clear(opt);
  fMCCherClus->Clear(opt);
  
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
  
  div_t d = div( h_chan, fManager->GetChanPerSlot() );
  slot = d.quot;
  chan = d.rem;

  d = div( slot, fManager->GetSlotPerCrate() );
  crate = d.quot;
  slot  = d.rem;
}

//-----------------------------------------------------------------------------
static inline
Int_t MakeROCKey( Int_t crate, Int_t slot, Int_t chan )
{
  return chan +
    fManager->GetChanPerSlot()*( slot + fManager->GetSlotPerCrate()*crate );
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

/*
//-----------------------------------------------------------------------------
MCHitInfo TSBSSimDecoder::GetMCHitInfo( Int_t crate, Int_t slot, Int_t chan ) const
{
  // Get MC truth info for the given hardware channel
  
  const char* const here = "TSBSSimDecoder::GetMCHitInfo";

  // Int_t istrip = StripFromROC( crate, slot, chan );
  // assert( istrip >= 0 );  // else logic error in caller or bad fPMTMap
  
  // assert( buffer );       // Must still have the event buffer
  // const TSBSSimEvent* simEvent = reinterpret_cast<const TSBSSimEvent*>(buffer);

  // assert( static_cast<vsiz_t>(istrip) < simEvent->fGEMStrips.size() );
  // const TSBSSimEvent::DigiGEMStrip& strip = simEvent->fGEMStrips[istrip];
  // assert( strip.fProj >= 0 && strip.fProj < fManager->GetNReadOut() );

  MCHitInfo mc;
  
  // if(strip.fProj==0 && strip.fPlane==4 && strip.fTime1>50.0)
  //   printf("%f \n", strip.fTime1);
  
  if (TESTBIT(strip.fSigType, kInducedStrip) && !TESTBIT(strip.fSigType, kPrimaryStrip) &&
      !TESTBIT(strip.fSigType, kSecondaryStrip) ){
    mc.fMCTrack = 0;
    mc.fMCPos = fManager->GetPosFromSectorStrip(strip.fProj, strip.fPlane, strip.fSector, strip.fChan);
    mc.fMCTime = strip.fTime1;
    //printf("strip num: %d, proj: %d, type: %d, time: %f\n", strip.fChan, strip.fProj, strip.fSigType, strip.fTime1);
    //cout << "strip = " << strip.fChan << ", time = " << mc.fMCTime << ", pos = " <<  mc.fMCPos << endl;
    return mc;
  }
  
  Double_t nOverlapSignal = 0.;
  for( Int_t i = 0; i<strip.fClusters.GetSize(); ++i ) {
    Int_t iclust = strip.fClusters[i] - 1;  // yeah, array index = clusterID - 1
    assert( iclust >= 0 && static_cast<vsiz_t>(iclust) < simEvent->fGEMClust.size() );
    const TSBSSimEvent::GEMCluster& c = simEvent->fGEMClust[iclust];
    assert( c.fID == iclust+1 );
    assert( strip.fPlane == c.fPlane && strip.fSector == c.fSector );
    
    Int_t signalID = -1;
    for (unsigned int ii = 0; ii<fSignalInfo.size(); ii++){
      if (c.fType == fSignalInfo.at(ii).tid && c.fPID == fSignalInfo.at(ii).pid)
	signalID = ii;
    }
    // cout << "Plane " << strip.fPlane << ", proj (x: 0, y: 1) " << strip.fProj 
    //  	 << ": pos[proj] = "  << c.fXProj[strip.fProj] << endl;
    if( signalID >= 0 && c.fSource == kPrimarySource ) {
      if( mc.fMCTrack > 0 ) {
        //this means that there two signal hits overlapping
        //for now I keep the fMCTrack to the first one, by average the fMCPos nad fMCTime
        //Weizhi Xiong
        //assert(manager->GetNSigParticle() > 1); //otherwise should not happen
        
        mc.fMCPos += c.fXProj[strip.fProj];
        mc.fMCTime += c.fTime; 
      }else{
        // Strip contains a contribution from a primary particle hit :)
        mc.fMCTrack = fSignalInfo.at(signalID).tid; 
        mc.fMCPos   = c.fXProj[strip.fProj];
        mc.fMCTime  = c.fTime;
      }
      nOverlapSignal++;
    } else {
      ++mc.fContam;
      if( mc.fMCTrack == 0 ) {
	mc.fMCPos += c.fXProj[strip.fProj];
	//printf("strip num: %d, proj: %d, type: %d, time: %f, clust: %d, time: %f \n", strip.fChan, strip.fProj, strip.fSigType, strip.fTime1, iclust, c.fTime);
      }
    }
  }
  assert( strip.fClusters.GetSize() == 0 || mc.fMCTrack > 0 || mc.fContam > 0 );
  
  if( mc.fMCTrack == 0 ) {
    if( mc.fContam > 1 ) {
      // If only background hits, report the mean position of all those hits
      mc.fMCPos /= static_cast<Double_t>(mc.fContam);
    }
    mc.fMCTime = strip.fTime1;
  }else{
    mc.fMCPos /= nOverlapSignal;
    mc.fMCTime /= nOverlapSignal;
  }
  
  
  for( Int_t i = 0; i<strip.fClusters.GetSize(); ++i ) {
    Int_t iclust = strip.fClusters[i] - 1;  // yeah, array index = clusterID - 1
    assert( iclust >= 0 && static_cast<vsiz_t>(iclust) < simEvent->fGEMClust.size() );
    const TSBSSimEvent::GEMCluster& c = simEvent->fGEMClust[iclust];
    assert( c.fID == iclust+1 );
    assert( strip.fPlane == c.fPlane && strip.fSector == c.fSector );
    if( c.fType == kPrimaryType && c.fSource == kPrimarySource ) {
      if( mc.fMCTrack > 0 ) {
	Warning( Here(here), "Event %d: Multiple hits of primary particle "
		 "in plane %d\nShould never happen. Call expert.",
		 simEvent->fEvtID, strip.fPlane );
	continue;
      }
      // Strip contains a contribution from a primary particle hit :)
      mc.fMCTrack = 1;    // currently only ever one primary particle per event
      mc.fMCPos   = c.fXProj[strip.fProj];
      mc.fMCTime  = c.fTime;
    } else {
      ++mc.fContam;
      if( mc.fMCTrack == 0 ) {
	mc.fMCPos += c.fXProj[strip.fProj];
      }
    }
  }
  assert( strip.fClusters.GetSize() == 0 || mc.fMCTrack > 0 || mc.fContam > 0 );

  if( mc.fMCTrack == 0 ) {
    if( mc.fContam > 1 ) {
      // If only background hits, report the mean position of all those hits
      mc.fMCPos /= static_cast<Double_t>(mc.fContam);
    }
    mc.fMCTime = strip.fTime1;
  }

  // cout << "TSBSSimDecoder: GetMCHitInfo: strip ADC sample: crate " << crate << " slot " << slot << " chan " << chan<< endl;
  // for( int k=0; k<strip.fNsamp; k++ ) {
  //   cout << strip.fADC[k] << " ";
  // }
  // cout << endl;
  
  
  return mc;
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
  Int_t crate, slot, chan;

  std::map<Decoder::THaSlotData*, std::vector<UInt_t> > hcalmap;
  std::cerr << "\n\n\n\n\nStart Processing event: " << event_num << std::endl;
  for(std::vector<TSBSSimEvent::DetectorData>::const_iterator it =
      simEvent->fDetectorData.begin(); it != simEvent->fDetectorData.end();
      ++it )
  {
    if((*it).fDetID == 2 && (*it).fData.size() > 0) { // HCal
      //crate = (*it).fCrate;
      //slot  = (*it).fSlot;
      //chan  = (*it).fChannel;
      int mod =  (*it).fChannel;
      if(mod < 192) {
        crate = 10;
        slot = 4+(mod/16);
      } else {
        crate = 11;
        slot = 4+(mod-192)/16;
      }
      chan = mod%16;
      Decoder::THaSlotData *sldat = crateslot[idx(crate,slot)];
      if(sldat) { // meaning the module is available
        std::vector<UInt_t> *myev = &(hcalmap[sldat]);
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
  }

  // Now that all data is prepared, loop through the hcal map and load
  // the appropriate slots
  std::cerr << "HCal hits: " << hcalmap.size() << std::endl;
  for( std::map<Decoder::THaSlotData*, std::vector<UInt_t> >::iterator it =
      hcalmap.begin(); it != hcalmap.end(); ++it) {
    it->first->GetModule()->LoadSlot(it->first,
        it->second.data(),0,it->second.size() );
        //myevbuff.resize((*it).fData.size());
        //myevbuff[0] = chan; // Only value that needs to be changed
        //for(size_t k = 1; k < myevbuff.size(); k++) {
        //  myevbuff[k] = (*it).fData[k];
        //}
        //if(myevbuff.size() > 1) { // i.e. have actual data
        //  std::cerr << "Det data, mod: " << mod << ", Crate: " << crate
        //    << ", slot=" << slot << ", chan=" << chan
        //    << ", detID: " << (*it).fDetID
        //    << ", size: " << myevbuff.size() << std::endl;
        //  sldat->GetModule()->LoadSlot(sldat,myevbuff.data(),0,myevbuff.size());
        //}
  }

  // Decode the digitized PMT data.  Populate crateslot array.
  for( vector<TSBSSimEvent::PMTHit>::size_type i = 0;
       i < simEvent->fPMTHits.size() && false; ++i ) {
    const TSBSSimEvent::PMTHit& h = simEvent->fPMTHits[i];
    
    PMTtoROC(h.fChannel, crate, slot, chan);
    
    /*
    cout << "h.fChannel, crate, slot, chan " << h.fChannel << " " << crate << " " << slot << " " << chan << endl;
    cout << "index, h.fTDC1, 2 " << idx(crate,slot) << " " << h.fTDCtime[0] << " " << h.fTDCtime[1] << endl;
    UInt_t test = h.fTDC[0];
    cout << "TDC word 0 = " << test << " -2^31 = ";
    test+=pow(2, 31);
    cout << test << "; =>  +2^31 = ";
    test+=pow(2, 31);
    cout << test << endl;
    */
    
    // NB: the VETROC words are subtracted 2^31 to convert them from Uint to Int... 
    // This *should be temporary* !!!
    if( crateslot[idx(crate,slot)]->loadData("tdc",chan, h.fTDC[0]-pow(2, 31), h.fTDC[0]-pow(2, 31)) == SD_ERR )
      return HED_ERR;
    
    if( crateslot[idx(crate,slot)]->loadData("tdc",chan, h.fTDC[1]-pow(2, 31), h.fTDC[1]-pow(2, 31)) == SD_ERR )
      return HED_ERR;
    
    /*
    cout << "crateslot crate, slot, NumChan " << crateslot[idx(crate,slot)]->getCrate() << " " 
	 << crateslot[idx(crate,slot)]->getSlot() << " " 
	 << crateslot[idx(crate,slot)]->getNumChan() << endl;
    for(int j_ = 0; j_<510; j_++)
      {
	if(crateslot[idx(crate,slot)]->getNumHits(j_)>0)
	  cout << "number of hits for channel " << j_ << ": " << crateslot[idx(crate,slot)]->getNumHits(j_) << endl;
	for(int j__ = 0; j__<crateslot[idx(crate,slot)]->getNumHits(j_); j__++)
	  cout << crateslot[idx(crate,slot)]->getData(j_, j__) << endl;
      }
    */
    
    newclus = true;
    
    new( (*fMCCherHits)[GetNPMThits()] ) TSBSSimPMTHit(h);
    
    for(int j = 0; j<GetNPMTclus(); j++){
      
      TSBSSimCherCluster* c = GetPMTclus(j);
      
      if(h.fMCtrackID==c->fMCtrackID){
	//Add a hit to the cluster
	newclus = false;
	c->fHitList.push_back(h.fChannel);
	//if( !fHitList ) fHitList = new TList;
	//fHitList->AddLast( theHit );
	c->fSize+=1;
	Int_t listnewsize = c->fSize;
	c->fXcenter = (c->fXcenter*((Double_t)(listnewsize-1))+h.fXPMT)/((Double_t)listnewsize);
	c->fYcenter = (c->fYcenter*((Double_t)(listnewsize-1))+h.fYPMT)/((Double_t)listnewsize);
	
	c->fXcenter_w = c->fXcenter_w*c->fNpe;
	c->fYcenter_w = c->fYcenter_w*c->fNpe;
	c->fNpe+= h.fNpe;
	c->fXcenter_w+= h.fNpe*h.fXPMT;
	c->fYcenter_w+= h.fNpe*h.fYPMT;
	c->fXcenter_w = c->fXcenter_w/c->fNpe;
	c->fYcenter_w = c->fYcenter_w/c->fNpe;
	
	c->fMeanRisingTime = (c->fMeanRisingTime*((Double_t)(listnewsize-1))+h.fTDCtime[0])/((Double_t)listnewsize);
	c->fMeanFallingTime = (c->fMeanFallingTime*((Double_t)(listnewsize-1))+h.fTDCtime[1])/((Double_t)listnewsize);
	c->fRisingTimeRMS = sqrt((pow(c->fRisingTimeRMS, 2)*((Double_t)(listnewsize-1))+ pow(h.fTDCtime[0], 2))/
			      ((Double_t)listnewsize));
	c->fFallingTimeRMS = sqrt((pow(c->fFallingTimeRMS, 2)*((Double_t)(listnewsize-1))+ pow(h.fTDCtime[1], 2))/
				  ((Double_t)listnewsize));
      }
      
    }
    
    if(newclus){
      TSBSSimCherCluster* c = new( (*fMCCherClus)[GetNPMTclus()] ) TSBSSimCherCluster();
      
      c->fHitList.push_back(h.fChannel);
      //if( !fHitList ) fHitList = new TList;
      //fHitList->AddLast( theHit );
      c->fMCtrackID = h.fMCtrackID;
      c->fMCtrackPID = h.fMCtrackPID;
      
      c->fSize+=1;
      Int_t listnewsize = c->fSize;
      c->fXcenter = h.fXPMT;
      c->fYcenter = h.fYPMT;
      
      c->fXcenter_w = c->fXcenter_w*c->fNpe;
      c->fYcenter_w = c->fYcenter_w*c->fNpe;
      c->fNpe = h.fNpe;
      c->fXcenter_w = h.fNpe*h.fXPMT;
      c->fYcenter_w = h.fNpe*h.fYPMT;
      
      c->fMeanRisingTime = h.fTDCtime[0];
      c->fMeanFallingTime = h.fTDCtime[1];
      c->fRisingTimeRMS = pow(h.fTDCtime[0], 2);
      c->fFallingTimeRMS = pow(h.fTDCtime[1], 2);	
    }
  }
  
  /*
  Bool_t newclus = true;
  //fill the 
  for(int i_ = 0; i_<fEvent->fMCClusterHitID.size(); i_++){
    if(tscd.GetParticleType(ih)==fEvent->fMCClusterHitID[i_].first){
      fEvent->fMCClusterHitID[i_].second.push_back(ih);
      newclus = false;
      break;
    }
  }
  if(newclus){
    Int_t trkID = tscd.GetParticleType(ih);
    std::vector<Short_t> ClusterPMTlist;
    ClusterPMTlist.push_back(ih);
    fEvent->fMCClusterHitID.push_back(make_pair(trkID, ClusterPMTlist));
  }
  */
  
  /*
  // Decode the digitized strip data.  Populate crateslot array.
  for( vector<TSBSSimEvent::DigiGEMStrip>::size_type i = 0;
       i < simEvent->fGEMStrips.size(); i++) {
    //cout << "i " << i << endl;
    const TSBSSimEvent::DigiGEMStrip& s = simEvent->fGEMStrips[i];
    Int_t crate, slot, chan;
    //cout << "striptoroc: " << endl;
    StripToROC( s.fPlane, s.fSector, s.fProj, s.fChan, crate, slot, chan );
    //cout << "crate = " << crate << ", slot = " << slot << ", chan " << chan << endl;
    //cout << "samples: " << endl;
    for( Int_t k = 0; k < s.fNsamp; k++ ) { 
      Int_t raw = s.fADC[k];
      //cout << raw << " ";
      if( crateslot[idx(crate,slot)]->loadData("adc",chan,raw,raw) == SD_ERR )
	return HED_ERR;
    }
    //cout << endl;
    //cout << "stripmap : " << endl;
    // Build map from ROC address to strip index. This is needed to extract
    // the MC truth info later in the tracking detector decoder via GetMCChanInfo.
#ifndef NDEBUG
    pair<PMTMap_t::const_iterator,bool> ins =
#endif
      fPMTMap.insert( make_pair( MakeROCKey(crate,slot,chan), i ) );
    // cout << "ROC key inserted in strip map " << endl;
    // cout << "ins.second ? " << ins.second << endl;
    assert( ins.second );
  }
  
  // Create lists of two types of tracks:
  // 1) Physics tracks, as generated at the target
  // 2) "Back tracks": hits in any GEM plane from the primary particle

  // Physics tracks. We need to copy them here so we can export them as global
  // variables.
  TClonesArray* tracks = simEvent->fMCTracks;
  assert( tracks );
  for( Int_t i = 0; i < tracks->GetLast()+1; i++ ) {
    TSBSSimTrack* trk = static_cast<TSBSSimTrack*>(tracks->UncheckedAt(i));
   new( (*fMCTracks)[i] ) TSBSSimTrack(*trk);
  }
  assert( GetNMCTracks() > 0 );

  // MC hit data ("clusters") and "back tracks"
  Int_t best_primary = -1, best_primary_plane = fManager->GetNChamber(), primary_sector = -1;
  UInt_t primary_hitbits = 0, ufail = 0, vfail = 0;
  for( vector<TSBSSimEvent::GEMCluster>::size_type i = 0;
       i < simEvent->fGEMClust.size(); ++i ) {
    const TSBSSimEvent::GEMCluster& c = simEvent->fGEMClust[i];

    if( c.fPlane < 0 || c.fPlane >= fManager->GetNChamber() ) {
      Error( here, "Illegal plane number = %d in cluster. "
	     "Should never happen. Call expert.", c.fPlane );
      simEvent->Print("clust");
      return HED_FATAL;
    }

    // Save hits in the GEMs
    new( (*fMCHits)[GetNMCHits()] ) TSBSSimGEMHit(c);

    // Extra bookkeeping for primary tracks, used for making back tracks below
    if( c.fType == kPrimaryType && c.fSource == kPrimarySource ) {
      // Record the primary track's points for access via the SimDecoder interface.
      // Record one point per projection so that we can study residuals.
      Int_t itrack = 1;
      primary_sector = c.fSector;
      MCTrackPoint* upt = // kUPlane changed to kXPlane: necessary to match TreeSearch EProjType
	new( (*fMCPoints)[GetNMCPoints()] ) MCTrackPoint( itrack,
							  c.fPlane, kXPlane,
							  c.fMCpos, c.fP );
      upt->fMCTime = c.fTime;
      MCTrackPoint* vpt =// kVPlane changed to kYPlane: necessary to match TreeSearch EProjType
	new( (*fMCPoints)[GetNMCPoints()] ) MCTrackPoint( itrack,
							  c.fPlane, kYPlane,
							  c.fMCpos, c.fP );
      vpt->fMCTime = c.fTime;
      
      // //debug...
      // cout << "TSBSSimDecoder.cxx: Print MC points " << endl;
      // cout << "kXplane ? " << kXPlane << endl;
      // upt->Print("");
      // cout << "kVYlane ? " << kYPlane << endl;
      // vpt->Print("");
      
      // Keep bitpattern of planes crossed by this primary
      SETBIT(primary_hitbits,c.fPlane);

      //cout << "Plane number " << c.fPlane << ", primary hitbits " << primary_hitbits << endl; 
      
      // Save index of the primary particle hit closest to plane 0
      if( c.fPlane < best_primary_plane ) {
	best_primary = i;
	best_primary_plane = c.fPlane;
      }
      // Determine digitization hit inefficiency: Check if this MC hit
      // activated GEM strips in both readout planes
      if( c.fSize[0] == 0 ) {
	SETBIT(ufail, c.fPlane);
	CLRBIT(upt->fStatus, MCTrackPoint::kDigitized);
      } else {
	SETBIT(upt->fStatus, MCTrackPoint::kDigitized);
      }
      if( c.fSize[1] == 0 ) {
	SETBIT(vfail, c.fPlane);
	CLRBIT(vpt->fStatus, MCTrackPoint::kDigitized);
      } else {
	SETBIT(vpt->fStatus, MCTrackPoint::kDigitized);
      }
    }
  }
  
  // Sort fMCPoints by type (u,v) and plane number, then calculate plane-to-plane
  // differences. The following assumes that all points are from the same track
  // (ensured above). If that is no longer so one day, fMCPoints will need to
  // be sorted by track number as well, and the algo below needs to be changed.
  fMCPoints->Sort();
  Double_t mass = 0;
  TSBSSimTrack* trk = static_cast<TSBSSimTrack*>(fMCTracks->UncheckedAt(0));
  assert(trk);
  if( TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(trk->fPID) )
    mass = particle->Mass();
  else
    Warning( "LoadEvent", "No enrty in PDG database for PID = %d", trk->fPID );

  MCTrackPoint* prev_pt = 0;
  for( Int_t i = 0; i < GetNMCPoints(); ++i ) {
    MCTrackPoint* pt = static_cast<MCTrackPoint*>( fMCPoints->UncheckedAt(i) );
    assert(pt);
    if( prev_pt && prev_pt->fType == pt->fType ) {
      assert( pt->fMCTrack == prev_pt->fMCTrack );
      if( prev_pt->fPlane+1 == pt->fPlane ) {
	pt->fDeltaE = TMath::Sqrt(prev_pt->fMCP.Mag2() + mass*mass) -
	  TMath::Sqrt(pt->fMCP.Mag2() + mass*mass);
	pt->fDeflect = prev_pt->fMCP.Angle(pt->fMCP);
	pt->fToF = pt->fMCTime - prev_pt->fMCTime;
      }
    }
    prev_pt = pt;
  }

  // Keep statistics in the MC track
  trk->fNHits = 2*NumberOfSetBits(primary_hitbits);
  trk->fHitBits = primary_hitbits;

  // "Back tracks"
  // Record the apparent track from the primary particle
  // of the signal data here, i.e. type == 1 and source == 0.
  // There is only ever one primary particle per event.
  if( best_primary >= 0 ) {

    Int_t nback = GetNBackTracks();
    assert( nback == 0 );

    TSBSSimBackTrack* btr = new( (*fBackTracks)[nback] )
      TSBSSimBackTrack(simEvent->fGEMClust[best_primary]);

    //cout << "Backtrack primary hitbits " << primary_hitbits << endl;
    
    btr->SetHitBits(primary_hitbits);
    btr->SetUfailBits(ufail);
    btr->SetVfailBits(vfail);

    // Use the back track to emulate calorimeter hits.
    // Assumptions:
    // - Only tracks crossing all fManager->GetNChamber() GEMs (points in all planes)
    //   make a calorimeter hit. This is a crude model for the trigger.
    // - The track propagates without deflection from the last GEM plane
    //   to the front of the emulated calorimeter.
    // - The measured calorimeter position is independent of the incident
    //   track angle.
    if( fManager->DoCalo() && trk->fNHits == 2*fManager->GetNChamber() ) {
      // Retrieve last MC track point
      assert( GetNMCPoints() == 2*fManager->GetNChamber() );
      MCTrackPoint* pt =
	static_cast<MCTrackPoint*>( fMCPoints->UncheckedAt(2*fManager->GetNChamber()-1) );
      assert( pt );
      const TVector3& pos = pt->fMCPoint;
      TVector3 dir = pt->fMCP.Unit();
      if( fManager->GetCaloZ() <= pos.Z() ) {
	Error( here, "Calorimeter z = %lf less than z of last GEM plane = %lf. "
	       "Set correct value with SetCaloZ() or turn off calo emulation.",
	       fManager->GetCaloZ(), pos.Z() );
	return HED_FATAL;
      }
      if( TMath::Abs(dir.Z()) < 1e-6 ) {
	Error( here, "Illegal primary track direction (%lf,%lf,%lf). "
	       "Should never happen. Call expert.", dir.X(), dir.Y(), dir.Z() );
	return HED_ERR;
      }
      dir *= 1.0/dir.Z();  // Make dir a transport vector
      TVector3 hitpos = pos + (fManager->GetCaloZ()-pos.Z()) * dir;

      // Smear the position with the given resolution
      // Assumes z-axis normal to calorimeter plane. Otherwise we would have to
      // get the plane's fXax and fYax
      TVector3 res( gRandom->Gaus(0.0, fManager->GetCaloRes()),
		    gRandom->Gaus(0.0, fManager->GetCaloRes()), 0.0 );
      hitpos += res;

      // Encode the raw hit data for the dummy GEM planes.
      // The actual coordinate transformation to u or v takes place in each
      // plane's Decode() where all the required geometry information is at hand.
      // This bypasses any type of digitization. That should be fine for dummy
      // planes where we want to inject known lab hit positions into the tracking.
      //
      // Because of the way the detector map is layed out at the moment,
      // we place the calorimeter in fake sector 31, so the data are in two slots
      // (for u and v coordinates, respectively) in the ROC immediately
      // following the GEM trackers for sector 30. In each slot, channels
      // 0-29 correspond to the sector of the MC track sector (should always be 0
      // if mapping sectors. Each "hit" corresponds to one measured position.
      // Currently, there is only ever one hit per channel since there is only
      // one MC track. The hit's raw data are hitpos.X(), the data, hitpos.Y(),
      // each a Float_t value interpreted as Int_t.
      assert( primary_sector == 0 );

      union FloatIntUnion {
	Float_t f;
	Int_t   i;
      } datx, daty;
      datx.f = static_cast<Float_t>(hitpos.X());
      daty.f = static_cast<Float_t>(hitpos.Y());

      Int_t crate, slot, chan;
      //StripToROC( 0, fManager->GetNSector(), kUPlane, primary_sector, crate, slot, chan );
      StripToROC( 0, fManager->GetNSector(), kXPlane, primary_sector, crate, slot, chan );
       if( crateslot[idx(crate,slot)]->loadData("adc",chan,datx.i,daty.i)
	  == SD_ERR )
	return HED_ERR;
      //StripToROC( 0, fManager->GetNSector(), kVPlane, primary_sector, crate, slot, chan );
      StripToROC( 0, fManager->GetNSector(), kYPlane, primary_sector, crate, slot, chan );
      if( crateslot[idx(crate,slot)]->loadData("adc",chan,datx.i,daty.i)
	  == SD_ERR )
	return HED_ERR;
    }
  }

  // DEBUG:
  //cout << "SimDecoder: nTracks = " << GetNMCTracks() << endl;
  //fMCTracks.Print();
  */

  std::cerr << "End Processing event:   " << event_num << std::endl;
  return HED_OK;
}

// --------------- //
// Hit functions  //
// ------------- //

//-----------------------------------------------------------------------------
TSBSSimPMTHit::TSBSSimPMTHit( )
  : fID(0), 
    fSource(0), fType(0), 
    fMCtrackID(0), fMCtrackPID(0), 
    fOrigVolFlag(0), 
    fXPMT(0), fYPMT(0), 
    fNpe(0), fTime(0), 
    fTDCtime1(0), fTDCtime2(0),
    fDetID(0), fChannel(0), 
    fPMTrow(0), fPMTcol(0), 
    fTDC1(0), fTDC2(0)
{

}

//-----------------------------------------------------------------------------
TSBSSimPMTHit::TSBSSimPMTHit( const TSBSSimEvent::PMTHit& h )
  : fID(h.fID), 
    fSource(h.fSource), fType(h.fType), 
    fMCtrackID(h.fMCtrackID), fMCtrackPID(h.fMCtrackPID), 
    fOrigVolFlag(h.fOrigVolFlag), 
    fXPMT(h.fXPMT), fYPMT(h.fYPMT), 
    fNpe(h.fNpe), fTime(h.fTime), 
    fTDCtime1(h.fTDCtime[0]), fTDCtime2(h.fTDCtime[1]),
    fDetID(h.fDetID), fChannel(h.fChannel), 
    fPMTrow(h.fPMTrow), fPMTcol(h.fPMTcol), 
    fTDC1(h.fTDC[0]), fTDC2(h.fTDC[1])
{

}

//-----------------------------------------------------------------------------
void TSBSSimPMTHit::Print( const Option_t* ) const
{
  cout << "hit = " << fID 
       << ", type = "    << fType 
       << ", src = "    << fSource 
       << endl;
  cout << "PMT row = "   << fPMTrow 
       << ", PMT col = "   << fPMTcol 
       << ", X_PMT = "   << fXPMT 
       << " m, Y_PMT = "   << fYPMT 
       << " m, number of photoelectrons = " << fNpe 
       << ", time = " << fTime 
       << " ns; TDC rising time = " << fTDCtime1
       << " ns, TDC falling time = " << fTDCtime2 
       << endl;
  cout << " Channel = " << fChannel 
       << ", TDC 'words' for rising time: " << fTDC1
       << "; for falling time: " << fTDC2 
       << endl;
}

// ------------------- //
// Cluster functions  //
// ----------------- //

//-----------------------------------------------------------------------------
TSBSSimCherCluster::TSBSSimCherCluster()
  : fSize(0),
    fXcenter(0), fYcenter(0),
    fXcenter_w(0), fYcenter_w(0),
    fNpe(0), 
    fMeanRisingTime(0), fMeanFallingTime(0),
    fRisingTimeRMS(0), fFallingTimeRMS(0),
    fMCtrackID(0), fMCtrackPID(0)
{
  fHitList.clear();
}

//-----------------------------------------------------------------------------
void TSBSSimCherCluster::Print( const Option_t* ) const
{
  cout << "Cluster size: " << fSize << "." << endl;
  cout << " hit list: " << endl;
  for(int i = 0; i<fSize; i++){
    cout << fHitList.at(i) << "; ";
  }
  cout << endl;
  
  cout << "Cluster center: X = " << fXcenter 
       << " m, Y = " << fYcenter 
       << " m." << endl;
  cout << "(weighted with Npe: X = " << fXcenter_w 
       << " m, Y = " << fYcenter_w 
       << " m.)" << endl;
  cout << "Cluster mean rising time " << fMeanRisingTime
       << " ns (RMS " << fRisingTimeRMS 
       << " ns)" << endl;
  cout << "Cluster mean falling time " << fMeanFallingTime
       << " ns (RMS " << fFallingTimeRMS 
       << " ns)" << endl;
  cout << "total number of photoelectrons " << fNpe 
       << "MC track PID " << fMCtrackPID 
       << endl;
}

