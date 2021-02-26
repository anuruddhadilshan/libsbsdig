#include "TSBSSimDigitizer.h"
#include "TSBSGeant4File.h"
#include "TSBSSimDetector.h"
#include "TSBSSimEvent.h"
#include "TSBSSimAuxi.h"
#include "TSBSDBManager.h"
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TChainElement.h>
#include "TGEMSBSDBManager.h"

TSBSSimDigitizer::TSBSSimDigitizer(const char* outputfilename) :
  fMarkInterval(0), fVerbose(1)
{
  cout << "Initialize TSBSSimDigitzer " << endl;
  //clear everything - in case.
  fSources.clear();
  fSourceChainMap.clear();
  fSourceWeightMap.clear();
  
  fManager = TSBSDBManager::GetInstance();
  fRN = TRndmManager::GetInstance();
  
  cout << "instantiated DB and RN managers" << endl;
  
  fEvent = new TSBSSimEvent();//problem here
  cout << "instantiated SimEvent" << endl;
  fOutFile = new TFile(outputfilename,"RECREATE");
  fOutTree = new TTree("digtree","");
  //fOutTree->Branch("fEvents",&fEvent);
  //fOutTree->Branch("event",&fEvent);
  fOutTree->Branch("RunID",&fEvent->fRunID);
  fOutTree->Branch("EvtID",&fEvent->fEvtID);
  fOutTree->Branch("Weight",&fEvent->fWeight);
  fOutTree->Branch("NSignal",&fEvent->fNSignal);
  
  cout << "declared event info for the output tree" << endl;
  
  TGEMSBSDBManager* GEMDBManager;
  std::string fullgemname;
  
  const std::vector<TDetInfo> AllDetInfo = fManager->GetAllDetInfo();
  for(uint i = 0; i<AllDetInfo.size(); i++){
    //SimDetData_Channel
    TDetInfo DetInfo_i = AllDetInfo.at(i);
    std::string fulldetname = DetInfo_i.DetFullName();
    det_type dettype = DetInfo_i.DetType();
    cout << fulldetname.c_str() << endl;

    switch(DetInfo_i.DetType()){
    case(kGEM):
      //MC track info
      fOutTree->Branch(Form("%s.mctrack.ntracks", fulldetname.c_str()),&fEvent->fMCTrackOutData[fulldetname.c_str()].fNMCTracks);
      fOutTree->Branch(Form("%s.mctrack.source", fulldetname.c_str()),&fEvent->fMCTrackOutData[fulldetname.c_str()].fMCTrackSource);
      fOutTree->Branch(Form("%s.mctrack.trid", fulldetname.c_str()),&fEvent->fMCTrackOutData[fulldetname.c_str()].fMCTrackTRID);
      fOutTree->Branch(Form("%s.mctrack.pid", fulldetname.c_str()),&fEvent->fMCTrackOutData[fulldetname.c_str()].fMCTrackPID);
      fOutTree->Branch(Form("%s.mctrack.x", fulldetname.c_str()),&fEvent->fMCTrackOutData[fulldetname.c_str()].fMCTrackX);
      fOutTree->Branch(Form("%s.mctrack.y", fulldetname.c_str()),&fEvent->fMCTrackOutData[fulldetname.c_str()].fMCTrackY);
      fOutTree->Branch(Form("%s.mctrack.t", fulldetname.c_str()),&fEvent->fMCTrackOutData[fulldetname.c_str()].fMCTrackT);
      fOutTree->Branch(Form("%s.mctrack.p", fulldetname.c_str()),&fEvent->fMCTrackOutData[fulldetname.c_str()].fMCTrackP);
      fOutTree->Branch(Form("%s.mctrack.dx", fulldetname.c_str()),&fEvent->fMCTrackOutData[fulldetname.c_str()].fMCTrackdX);
      fOutTree->Branch(Form("%s.mctrack.dy", fulldetname.c_str()),&fEvent->fMCTrackOutData[fulldetname.c_str()].fMCTrackdY);
      fOutTree->Branch(Form("%s.mctrack.xv", fulldetname.c_str()),&fEvent->fMCTrackOutData[fulldetname.c_str()].fMCTrackXv);
      fOutTree->Branch(Form("%s.mctrack.yv", fulldetname.c_str()),&fEvent->fMCTrackOutData[fulldetname.c_str()].fMCTrackYv);	
      fOutTree->Branch(Form("%s.mctrack.zv", fulldetname.c_str()),&fEvent->fMCTrackOutData[fulldetname.c_str()].fMCTrackZv);	
      fOutTree->Branch(Form("%s.mctrack.pxv", fulldetname.c_str()),&fEvent->fMCTrackOutData[fulldetname.c_str()].fMCTrackPXv);
      fOutTree->Branch(Form("%s.mctrack.pyv", fulldetname.c_str()),&fEvent->fMCTrackOutData[fulldetname.c_str()].fMCTrackPYv);	
      fOutTree->Branch(Form("%s.mctrack.pzv", fulldetname.c_str()),&fEvent->fMCTrackOutData[fulldetname.c_str()].fMCTrackPZv);
      fOutTree->Branch(Form("%s.mctrack.weight", fulldetname.c_str()),&fEvent->fMCTrackOutData[fulldetname.c_str()].fMCTrackWeight);
      
      //MC hit info
      fOutTree->Branch(Form("%s.simhit.nhits", fulldetname.c_str()),&fEvent->fSimGEMHitMCOutData[fulldetname.c_str()].fNSimHits);
      fOutTree->Branch(Form("%s.simhit.src", fulldetname.c_str()),&fEvent->fSimGEMHitMCOutData[fulldetname.c_str()].fSimSource);
      fOutTree->Branch(Form("%s.simhit.trid", fulldetname.c_str()),&fEvent->fSimGEMHitMCOutData[fulldetname.c_str()].fSimTRID);
      fOutTree->Branch(Form("%s.simhit.pid", fulldetname.c_str()),&fEvent->fSimGEMHitMCOutData[fulldetname.c_str()].fSimPID);
      fOutTree->Branch(Form("%s.simhit.plane", fulldetname.c_str()),&fEvent->fSimGEMHitMCOutData[fulldetname.c_str()].fPlane);
      fOutTree->Branch(Form("%s.simhit.module", fulldetname.c_str()),&fEvent->fSimGEMHitMCOutData[fulldetname.c_str()].fModule);
      fOutTree->Branch(Form("%s.simhit.edep", fulldetname.c_str()),&fEvent->fSimGEMHitMCOutData[fulldetname.c_str()].fSimEdep);
      fOutTree->Branch(Form("%s.simhit.time", fulldetname.c_str()),&fEvent->fSimGEMHitMCOutData[fulldetname.c_str()].fSimTime);
      fOutTree->Branch(Form("%s.simhit.xpos", fulldetname.c_str()),&fEvent->fSimGEMHitMCOutData[fulldetname.c_str()].fXpos);
      fOutTree->Branch(Form("%s.simhit.ypos", fulldetname.c_str()),&fEvent->fSimGEMHitMCOutData[fulldetname.c_str()].fYpos);
      fOutTree->Branch(Form("%s.simhit.px", fulldetname.c_str()),&fEvent->fSimGEMHitMCOutData[fulldetname.c_str()].fPX);
      fOutTree->Branch(Form("%s.simhit.py", fulldetname.c_str()),&fEvent->fSimGEMHitMCOutData[fulldetname.c_str()].fPY);
      fOutTree->Branch(Form("%s.simhit.pz", fulldetname.c_str()),&fEvent->fSimGEMHitMCOutData[fulldetname.c_str()].fPZ);
      fOutTree->Branch(Form("%s.simhit.sizex", fulldetname.c_str()),&fEvent->fSimGEMHitMCOutData[fulldetname.c_str()].fSizeX);
      fOutTree->Branch(Form("%s.simhit.sizey", fulldetname.c_str()),&fEvent->fSimGEMHitMCOutData[fulldetname.c_str()].fSizeY);
      fOutTree->Branch(Form("%s.simhit.startx", fulldetname.c_str()),&fEvent->fSimGEMHitMCOutData[fulldetname.c_str()].fStartX);
      fOutTree->Branch(Form("%s.simhit.starty", fulldetname.c_str()),&fEvent->fSimGEMHitMCOutData[fulldetname.c_str()].fStartY);
      
      //digitized info
      GEMDBManager = DetInfo_i.GetGEMDB();
      for(int ipl = 0; ipl<GEMDBManager->GetNGEMPlane(); ipl++){
	//for(int imod = 0; imod<GEMDBManager->GetNModule(ipl); imod++){
	for(int ipr = 0; ipr<GEMDBManager->GetNReadOut(); ipr++){
	  // fullgemname = Form("%s.p%d.m%d.%s", 
	  // 		     fulldetname.c_str(), 
	  // 		     ipl+1, imod+1, kProj_str[ipr].c_str());
	  fullgemname = Form("%s.%d.%s", 
			     fulldetname.c_str(), 
			     ipl+1, kProj_str[ipr].c_str());
	  if(fDebug>=3)cout << fullgemname.c_str() << endl;
	  fOutTree->Branch(Form("%s.hit.nhits", fullgemname.c_str()),&fEvent->fSimDigSampOutData[fullgemname.c_str()].fNHits);
	  fOutTree->Branch(Form("%s.hit.chan", fullgemname.c_str()),&fEvent->fSimDigSampOutData[fullgemname.c_str()].fChannel);
	  //  !!! "dataword" is the number of words here !!!
	  fOutTree->Branch(Form("%s.hit.dataword", fullgemname.c_str()),&fEvent->fSimDigSampOutData[fullgemname.c_str()].fDataWord);
	  fOutTree->Branch(Form("%s.hit.adc", fullgemname.c_str()),&fEvent->fSimDigSampOutData[fullgemname.c_str()].fADC);
	  fOutTree->Branch(Form("%s.hit.nsamps", fullgemname.c_str()),&fEvent->fSimDigSampOutData[fullgemname.c_str()].fNsamps);
	  fOutTree->Branch(Form("%s.hit.samps_adc", fullgemname.c_str()),&fEvent->fSimDigSampOutData[fullgemname.c_str()].fADC_samps);
	  fOutTree->Branch(Form("%s.hit.samps_datawords", fullgemname.c_str()),&fEvent->fSimDigSampOutData[fullgemname.c_str()].fDataWord_samps);
	  
	}
	//}
      }
      //}
      break;
    case(kHCal):
      //MC track projection
      fOutTree->Branch(Form("%s.trackmchit.nhits", fulldetname.c_str()),&fEvent->fTrackMCHitOutData[fulldetname.c_str()].fNTrackMCHits);
      fOutTree->Branch(Form("%s.trackmchit.source", fulldetname.c_str()),&fEvent->fTrackMCHitOutData[fulldetname.c_str()].fTrackMCSource);
      fOutTree->Branch(Form("%s.trackmchit.trid", fulldetname.c_str()),&fEvent->fTrackMCHitOutData[fulldetname.c_str()].fTrackMCTRID);
      fOutTree->Branch(Form("%s.trackmchit.pid", fulldetname.c_str()),&fEvent->fTrackMCHitOutData[fulldetname.c_str()].fTrackMCPID);
      fOutTree->Branch(Form("%s.trackmchit.xhit", fulldetname.c_str()),&fEvent->fTrackMCHitOutData[fulldetname.c_str()].fTrackMCXhit);
      fOutTree->Branch(Form("%s.trackmchit.yhit", fulldetname.c_str()),&fEvent->fTrackMCHitOutData[fulldetname.c_str()].fTrackMCYhit);
      fOutTree->Branch(Form("%s.trackmchit.thit", fulldetname.c_str()),&fEvent->fTrackMCHitOutData[fulldetname.c_str()].fTrackMCThit);
      fOutTree->Branch(Form("%s.trackmchit.e", fulldetname.c_str()),&fEvent->fTrackMCHitOutData[fulldetname.c_str()].fTrackMCE);
      fOutTree->Branch(Form("%s.trackmchit.weight", fulldetname.c_str()),&fEvent->fTrackMCHitOutData[fulldetname.c_str()].fTrackMCWeight);
      
      //MC info
      fOutTree->Branch(Form("%s.simhit.nhits", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fNSimHits);
      fOutTree->Branch(Form("%s.simhit.src", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fSimSource);
      fOutTree->Branch(Form("%s.simhit.trid", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fSimTRID);
      fOutTree->Branch(Form("%s.simhit.pid", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fSimPID);
      fOutTree->Branch(Form("%s.simhit.chan", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fSimChannel);
      fOutTree->Branch(Form("%s.simhit.edep", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fSimEdep);
      fOutTree->Branch(Form("%s.simhit.npe", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fSimNpe);
      fOutTree->Branch(Form("%s.simhit.time", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fSimTime);
      fOutTree->Branch(Form("%s.simhit.t_lead", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fSimLeadTime);
      fOutTree->Branch(Form("%s.simhit.t_trail", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fSimTrailTime);

      //digitized info
      fOutTree->Branch(Form("%s.hit.nhits", fulldetname.c_str()),&fEvent->fSimDigSampOutData[fulldetname.c_str()].fNHits);
      fOutTree->Branch(Form("%s.hit.chan", fulldetname.c_str()),&fEvent->fSimDigSampOutData[fulldetname.c_str()].fChannel);
      // "dataword" could be the number of words or the TDC word
      fOutTree->Branch(Form("%s.hit.dataword", fulldetname.c_str()),&fEvent->fSimDigSampOutData[fulldetname.c_str()].fDataWord);
      fOutTree->Branch(Form("%s.hit.adc", fulldetname.c_str()),&fEvent->fSimDigSampOutData[fulldetname.c_str()].fADC);
      fOutTree->Branch(Form("%s.hit.nsamps", fulldetname.c_str()),&fEvent->fSimDigSampOutData[fulldetname.c_str()].fNsamps);
      fOutTree->Branch(Form("%s.hit.samps_adc", fulldetname.c_str()),&fEvent->fSimDigSampOutData[fulldetname.c_str()].fADC_samps);
      fOutTree->Branch(Form("%s.hit.samps_datawords", fulldetname.c_str()),&fEvent->fSimDigSampOutData[fulldetname.c_str()].fDataWord_samps);
      if(DetInfo_i.DigInfo().TDCBits()>0){
	fOutTree->Branch(Form("%s.hit.tdc_l", fulldetname.c_str()),&fEvent->fSimDigSampOutData[fulldetname.c_str()].fTDC_L);
	fOutTree->Branch(Form("%s.hit.tdc_t", fulldetname.c_str()),&fEvent->fSimDigSampOutData[fulldetname.c_str()].fTDC_T);
      }
      break;
    default:
      //Track MC hits
      fOutTree->Branch(Form("%s.trackmchit.nhits", fulldetname.c_str()),&fEvent->fTrackMCHitOutData[fulldetname.c_str()].fNTrackMCHits);
      fOutTree->Branch(Form("%s.trackmchit.source", fulldetname.c_str()),&fEvent->fTrackMCHitOutData[fulldetname.c_str()].fTrackMCSource);
      fOutTree->Branch(Form("%s.trackmchit.trid", fulldetname.c_str()),&fEvent->fTrackMCHitOutData[fulldetname.c_str()].fTrackMCTRID);
      fOutTree->Branch(Form("%s.trackmchit.pid", fulldetname.c_str()),&fEvent->fTrackMCHitOutData[fulldetname.c_str()].fTrackMCPID);
      fOutTree->Branch(Form("%s.trackmchit.xhit", fulldetname.c_str()),&fEvent->fTrackMCHitOutData[fulldetname.c_str()].fTrackMCXhit);
      fOutTree->Branch(Form("%s.trackmchit.yhit", fulldetname.c_str()),&fEvent->fTrackMCHitOutData[fulldetname.c_str()].fTrackMCYhit);
      fOutTree->Branch(Form("%s.trackmchit.thit", fulldetname.c_str()),&fEvent->fTrackMCHitOutData[fulldetname.c_str()].fTrackMCThit);
      fOutTree->Branch(Form("%s.trackmchit.e", fulldetname.c_str()),&fEvent->fTrackMCHitOutData[fulldetname.c_str()].fTrackMCE);
      fOutTree->Branch(Form("%s.trackmchit.weight", fulldetname.c_str()),&fEvent->fTrackMCHitOutData[fulldetname.c_str()].fTrackMCWeight);
      
      //MC info
      fOutTree->Branch(Form("%s.simhit.nhits", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fNSimHits);
      fOutTree->Branch(Form("%s.simhit.src", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fSimSource);
      fOutTree->Branch(Form("%s.simhit.trid", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fSimTRID);
      fOutTree->Branch(Form("%s.simhit.pid", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fSimPID);
      fOutTree->Branch(Form("%s.simhit.chan", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fSimChannel);
      if(dettype!=kCher)fOutTree->Branch(Form("%s.simhit.edep", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fSimEdep);
      fOutTree->Branch(Form("%s.simhit.npe", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fSimNpe);
      fOutTree->Branch(Form("%s.simhit.time", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fSimTime);
      if(DetInfo_i.DigInfo().TDCBits()>0){
	fOutTree->Branch(Form("%s.simhit.t_lead", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fSimLeadTime);
	fOutTree->Branch(Form("%s.simhit.t_trail", fulldetname.c_str()),&fEvent->fSimHitMCOutData[fulldetname.c_str()].fSimTrailTime);
      }
      
      //digitized info
      fOutTree->Branch(Form("%s.hit.nhits", fulldetname.c_str()),&fEvent->fSimDigOutData[fulldetname.c_str()].fNHits);
      fOutTree->Branch(Form("%s.hit.chan", fulldetname.c_str()),&fEvent->fSimDigOutData[fulldetname.c_str()].fChannel);
      fOutTree->Branch(Form("%s.hit.dataword", fulldetname.c_str()),&fEvent->fSimDigOutData[fulldetname.c_str()].fDataWord);
      if(DetInfo_i.DigInfo().ADCBits()>0)
	fOutTree->Branch(Form("%s.hit.adc", fulldetname.c_str()),&fEvent->fSimDigOutData[fulldetname.c_str()].fADC);
      if(DetInfo_i.DigInfo().TDCBits()>0){
	fOutTree->Branch(Form("%s.hit.tdc_l", fulldetname.c_str()),&fEvent->fSimDigOutData[fulldetname.c_str()].fTDC_L);
	fOutTree->Branch(Form("%s.hit.tdc_t", fulldetname.c_str()),&fEvent->fSimDigOutData[fulldetname.c_str()].fTDC_T);
      }
      break;
    }
    
  }
}

TSBSSimDigitizer::~TSBSSimDigitizer()
{
  fDetectors.clear();
  fSources.clear();
  fSourceChainMap.clear();
  fSourceWeightMap.clear();
}

int TSBSSimDigitizer::AddFileToEvent(TSBSGeant4File *f, double tjitter)
{
  if(fDebug>=3)cout << "Add in full file " << f->GetName() << endl; 
  
  if(!f)
    return 0;
  
  // Now loop through the file and digitize entries
  //bool has_data;
  while( f->ReadNextEvent(fDebug) ) {
    // Loop through all detectors and have them parse data vector
     if(fDebug>=3)cout << "f->GetDataVector().size() " << f->GetDataVector().size() << endl;
    for(size_t det = 0; det < fDetectors.size(); det++) {
      if(fDebug>=3){
	cout << "load event for det " << fDetectors[det]->GetName() << endl;
      }
      double t0 = tjitter+fRN->Uniform(-fManager->GetBkgdSpreadTimeWindowHW(), fManager->GetBkgdSpreadTimeWindowHW());
      fDetectors[det]->SetTimeZero(t0);
      fDetectors[det]->LoadAccumulateData(f->GetDataVector());
      fDetectors[det]->LoadMCTrackData(f->GetGenDataVector(), *fEvent);
      if(fDebug>=3)cout << "Done loading data in " << fDetectors[det]->GetName() << endl;
    }
  }
  
  return 1;
}

int TSBSSimDigitizer::Process(ULong_t max_events)
{
  if(fSourceChainMap.empty())return 0;
  
  int source = 0;
  ULong_t nevent = 0;
  int nfile = 0;
  Int_t nevt_b;
  //int ngood = 0;
  bool has_data;  
  double t0 = 0;
  double tjitter = 0;

  TObjArray *SigFileList = fSourceChainMap[0]->GetListOfFiles();
  TIter next_sig(SigFileList);
  TChainElement *chEl_sig = 0;

  std::vector<Int_t> Sources;
  std::vector<TObjArray*> BkgdFileLists;
  std::vector<TIter> BkgdIters;
    
  for(std::set<Int_t>::iterator it = fSources.begin(); it!=fSources.end(); ++it){
    source = *it;
    if(source==0)continue;//signal - we're already treating it
    if(fDebug>=2)cout << "source number " << source << endl;
    Sources.push_back(source);
    BkgdFileLists.push_back(fSourceChainMap[source]->GetListOfFiles());
    BkgdIters.push_back(TIter(BkgdFileLists.back()));
  }
  TChainElement *chEl_bkgd = 0;

  // Set a reasonable print interval of every 10% if fMarkInterval is 0
  if(fMarkInterval==0) {
    fMarkInterval = max_events*0.1;
  }
  
  while( (chEl_sig=(TChainElement*)next_sig()) && 
	 nevent<max_events ){
    if(fDebug>=2)cout << "Opening file " << chEl_sig->GetTitle() << endl;
    TSBSGeant4File* f = new TSBSGeant4File(chEl_sig->GetTitle());
    if(f->IsZombie())continue;
    f->Open();
    if(fDebug>=4)cout << "Set source 0 for file " << chEl_sig->GetTitle() << endl;
    f->SetSource(0);
    //G4SBSRunData *rd;
    
    while( nevent<max_events && f->ReadNextEvent(fDebug)) {
      if(fVerbose && fMarkInterval > 0 && nevent%fMarkInterval==0)
        std::cout << "Event: " << nevent << " / " << max_events << std::endl;
      if(fDebug>=3)cout << "clear event " << endl;
      if(fDebug>=1)cout << "Process event " << nevent << endl;
      fEvent->Clear();
      has_data = false;
      // Accumulate data here...
      
      //while( f->ReadNextEvent(fDebug) ){
      //if(!f->ReadNextEvent(fDebug))continue;
      if(f->GetDataVector().size()==0)continue;
      
      //tjitter is the trigger jitter, which I assume to be global for all detectors, and one per event, that needs to be propagated for all background that is superimposed to this event
      
      tjitter = fRN->Gaus(0.0, fManager->GetTriggerJitter());
      for(size_t det = 0; det < fDetectors.size(); det++) {
	if(fDebug>=2){
	  cout << "load event " << f->GetEvNum() << " for file " << f->GetName() 
	       << " det " << fDetectors[det]->GetName() << endl;
	}
	//"LoadEventData" for signal - we want everything cleaned up for signal
	if(fDebug>=3)cout << "f->GetDataVector().size() " << f->GetDataVector().size() << endl;
	fDetectors[det]->SetTimeZero(tjitter);
	fDetectors[det]->LoadEventData(f->GetDataVector());
	fDetectors[det]->LoadMCTrackData(f->GetGenDataVector(), *fEvent);
	if(fDebug>=3)cout << " event loaded for  " << fDetectors[det]->GetName() << endl;
      }
      //now loop on other chains to add background
      
      //std::set<Int_t>::iterator it = fSources.begin();
      //for(it; it!=fSources.end(); ++it){
      //  source = *it;
      for(uint i = 0; i<Sources.size(); i++){
	source = Sources[i];
	if(source==0)continue;//signal - we're already treating it
	if(fDebug>=2)cout << "source number " << source << endl;
	//TObjArray *BkgdFileList = fSourceChainMap[source]->GetListOfFiles();
	//TIter next_bkgd(BkgdFileLists[i]);
	
	//if the weight has a negative value, by convention we add a number of files, not a number of events
	if(fSourceWeightMap[source]<0){
	  nfile = 0;
	  //while( (chEl_bkgd=(TChainElement*)next_bkgd()) && 
	  while( (chEl_bkgd=(TChainElement*)BkgdIters[i]()) && 
		 nfile<abs(fSourceWeightMap[source]) ){
	    if(fDebug>=2)cout << "bkgd file counter " << nfile << "; Opening file " << chEl_bkgd->GetTitle() << endl;
	    TSBSGeant4File* f_b = new TSBSGeant4File(chEl_bkgd->GetTitle());
	    if(f_b->IsZombie())continue;
	    f_b->SetSource(source);
	    f_b->Open();

	    AddFileToEvent(f_b, tjitter);
	    //chEl_bkgd
	    //fSourceChainMap[source]->RecursiveRemove(f_b);
	    f_b->~TSBSGeant4File();
	    nfile++;
	  }//end loop on background files
	}else{
	  nevt_b = 0;
	  //while( (chEl_bkgd=(TChainElement*)next_bkgd()) && 
	  while( (chEl_bkgd=(TChainElement*)BkgdIters[i]()) && 
		 nevt_b<fSourceWeightMap[source] ){
	    if(fDebug>=3)cout << "Opening file " << chEl_bkgd->GetTitle() << endl;
	    //for that we will need to do it smarter
	    TSBSGeant4File* f_b = new TSBSGeant4File(chEl_bkgd->GetTitle());
	    if(f_b->IsZombie())continue;
	    f_b->SetSource(source);
	    if(!f_b->IsOpen())f_b->Open();
	    
	    if(fDebug>=3)cout << "f->GetDataVector().size() " << f_b->GetDataVector().size() << endl;
	    for(size_t det = 0; det < fDetectors.size(); det++) {
	      if(fDebug>=3){
		cout << "load event " << f_b->GetEvNum() << " for file " << f_b->GetName() 
		     << " det " << fDetectors[det]->GetName() << endl;
	      }
	      //"LoadAccumulateData" for any other stuff we want to superimpose to signal
	      //double 
	      t0 = tjitter+fRN->Uniform(-fManager->GetBkgdSpreadTimeWindowHW(), fManager->GetBkgdSpreadTimeWindowHW());
	      fDetectors[det]->SetTimeZero(t0);
	      fDetectors[det]->LoadAccumulateData(f->GetDataVector());
	      fDetectors[det]->LoadMCTrackData(f->GetGenDataVector(), *fEvent);
	      if(fDebug>=3)cout << "Done loading data in " << fDetectors[det]->GetName() << endl;
	    }
	    //AddFileToEvent(f_b);
	    //fSourceChainMap[source]->RecursiveRemove(f_b);
	  }//end loop on background events
	}
      }//end loop on sources
      
      for(size_t det = 0; det < fDetectors.size(); det++) {
	if(fDebug>=3)cout << "digitize det " << fDetectors[det]->GetName() << endl;
	fDetectors[det]->Digitize(*fEvent);
      }
      fEvent->fEvtID = nevent;
      
      // Fill in tree if any of the detectors have data to fill
      for(size_t det = 0; det < fDetectors.size() && !has_data; det++) {
	has_data = fDetectors[det]->HasData();
      }
      if(has_data) {
	// Write to the tree
	if(fDebug>=2)std::cout << "Have data for event: " << nevent << std::endl;
	fOutTree->Fill();
      }
      nevent++;
    }//end loop on events
    
  }//end loop on signal files
  
  //if(fDebug>=1)
  cout << "Done processing all events, write output file " << endl;
  fOutFile->Write();
  // Close files
  
  if(fDebug>=2)cout << "close output file " << endl;
  fOutFile->Close();
  
  return 0;
}


void TSBSSimDigitizer::AddDetector(TSBSSimDetector* detector)
{
  fDetectors.push_back(detector);
}

void TSBSSimDigitizer::AddInputFile(const char* filename, Int_t source, Int_t weight)
{
  fSources.insert(source);
  //if Chain does not exist yet, create it.
  if(fSourceChainMap[source]==0)fSourceChainMap[source] = new TChain("T");
  fSourceChainMap[source]->Add(filename);
  fSourceWeightMap[source] = weight;
}



ClassImp(TSBSSimDigitizer)
