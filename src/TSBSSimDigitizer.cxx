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

TSBSSimDigitizer::TSBSSimDigitizer(const char* outputfilename)
{
  cout << "Initialize TSBSSimDigitzer " << endl;
  //clear everything - in case.
  fSources.clear();
  fSourceChainMap.clear();
  fSourceWeightMap.clear();
  
  fManager = TSBSDBManager::GetInstance();
  fRN = TRndmManager::GetInstance();
  
  cout << "instantiated DB and RN managers" << endl;
  
  //fSigChain = new TChain("T");
  //fBkgdChain = new TChain("T");
  
  fEvent = new TSBSSimEvent();//problem here
  cout << "instantiated SimEvent" << endl;
  fOutFile = new TFile(outputfilename,"RECREATE");
  //fOutTree = new TTree("TSBSDigi","");
  fOutTree = new TTree("digtree","");
  //fOutTree->Branch("fEvents",&fEvent);
  //fOutTree->Branch("event",&fEvent);
  fOutTree->Branch("RunID",&fEvent->fRunID);
  fOutTree->Branch("EvtID",&fEvent->fEvtID);
  fOutTree->Branch("Weight",&fEvent->fWeight);
  fOutTree->Branch("NSignal",&fEvent->fNSignal);
  
  cout << "declared event info for the output tree" << endl;
  
  /*
  fOutTree->Branch("SimDetData_size",&fEvent->NSimDetData);
  fOutTree->Branch("SimDetData_DetID",&fEvent->SimDetID);
  fOutTree->Branch("SimDetData_Channel",&fEvent->SimDetChannel);
  fOutTree->Branch("SimDetData_DataType",&fEvent->SimDetDataType);
  fOutTree->Branch("SimDetData_Ndata",&fEvent->SimDetNData);
  fOutTree->Branch("SimDetData_Data",&fEvent->SimDetData);
  fOutTree->Branch("DetData_size",&fEvent->NDetData);
  fOutTree->Branch("DetData_DetID",&fEvent->DetID);
  fOutTree->Branch("DetData_Channel",&fEvent->DetChannel);
  fOutTree->Branch("DetData_Ndata",&fEvent->DetNData);
  fOutTree->Branch("DetData_Data",&fEvent->DetData);
  */
  
  const std::vector<TDetInfo> AllDetInfo = fManager->GetAllDetInfo();
  for(uint i = 0; i<AllDetInfo.size(); i++){
    //SimDetData_Channel
    TDetInfo DetInfo_i = AllDetInfo.at(i);
    std::string fulldetname = DetInfo_i.DetFullName();
    det_type dettype = DetInfo_i.DetType();
    cout << fulldetname.c_str() << endl;

    /*
    fOutTree->Branch(Form("NSimData_%s", fulldetname.c_str()),&fEvent->NSimDetData[fulldetname.c_str()]);
    fOutTree->Branch(Form("SimData_%s_Chan", fulldetname.c_str()),&fEvent->SimDetChannel[fulldetname.c_str()]);
    fOutTree->Branch(Form("SimData_%s_Type", fulldetname.c_str()),&fEvent->SimDetDataType[fulldetname.c_str()]);
    fOutTree->Branch(Form("SimData_%s_Ndata", fulldetname.c_str()),&fEvent->SimDetNData[fulldetname.c_str()]);
    fOutTree->Branch(Form("SimData_%s_Data", fulldetname.c_str()),&fEvent->SimDetData[fulldetname.c_str()]);
    
    fOutTree->Branch(Form("NData_%s", fulldetname.c_str()),&fEvent->NDetData[fulldetname.c_str()]);
    fOutTree->Branch(Form("Data_%s_Chan", fulldetname.c_str()),&fEvent->DetChannel[fulldetname.c_str()]);
    fOutTree->Branch(Form("Data_%s_Ndata", fulldetname.c_str()),&fEvent->DetNData[fulldetname.c_str()]);
    fOutTree->Branch(Form("Data_%s_Data", fulldetname.c_str()),&fEvent->DetData[fulldetname.c_str()]);
    */

    fOutTree->Branch(Form("%s_Nsimhits", fulldetname.c_str()),&fEvent->NSimDetHits[fulldetname.c_str()]);
    fOutTree->Branch(Form("%s_simhit_chan", fulldetname.c_str()),&fEvent->SimDetChannel[fulldetname.c_str()]);
    if(dettype!=kCher)fOutTree->Branch(Form("%s_simhit_Edep", fulldetname.c_str()),&fEvent->SimDetEdep[fulldetname.c_str()]);
    fOutTree->Branch(Form("%s_simhit_npe", fulldetname.c_str()),&fEvent->SimDetNpe[fulldetname.c_str()]);
    fOutTree->Branch(Form("%s_simhit_time", fulldetname.c_str()),&fEvent->SimDetTime[fulldetname.c_str()]);
    if(DetInfo_i.DigInfo().TDCBits()>0){//if(dettype!=kECal){
      fOutTree->Branch(Form("%s_simhit_t_lead", fulldetname.c_str()),&fEvent->SimDetLeadTime[fulldetname.c_str()]);
      fOutTree->Branch(Form("%s_simhit_t_trail", fulldetname.c_str()),&fEvent->SimDetTrailTime[fulldetname.c_str()]);
    }
    
    fOutTree->Branch(Form("%s_Nhits", fulldetname.c_str()),&fEvent->NDetData[fulldetname.c_str()]);
    fOutTree->Branch(Form("%s_hit_chan", fulldetname.c_str()),&fEvent->DetChannel[fulldetname.c_str()]);
    fOutTree->Branch(Form("%s_hit_dataword", fulldetname.c_str()),&fEvent->DetDataWord[fulldetname.c_str()]);
    // if(DetInfo_i.DigInfo().ADCBits()>0)
    //   fOutTree->Branch(Form("%s_adc", fulldetname.c_str()),&fEvent->DetADC[fulldetname.c_str()]);
    // if(DetInfo_i.DigInfo().TDCBits()>0)
    //   fOutTree->Branch(Form("%s_tdc", fulldetname.c_str()),&fEvent->DetTDC[fulldetname.c_str()]);
  }
  
  //fOutTree->Branch("SimDetectorData",&fEvent->fSimDetectorData);
  //fOutTree->Branch("DetectorData",&fEvent->fSimDetectorData);
}

TSBSSimDigitizer::~TSBSSimDigitizer()
{
  fDetectors.clear();
  fSources.clear();
  fSourceChainMap.clear();
  fSourceWeightMap.clear();
}

int TSBSSimDigitizer::AddFileToEvent(TSBSGeant4File *f)
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
      double t0 = fRN->Uniform(-fManager->GetBkgdSpreadTimeWindowHW(), 
			       fManager->GetBkgdSpreadTimeWindowHW());
      fDetectors[det]->SetTimeZero(t0);
      fDetectors[det]->LoadAccumulateData(f->GetDataVector());
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
  
  while( (chEl_sig=(TChainElement*)next_sig()) && 
	 nevent<max_events ){
    if(fDebug>=2)cout << "Opening file " << chEl_sig->GetTitle() << endl;
    TSBSGeant4File* f = new TSBSGeant4File(chEl_sig->GetTitle());
    if(f->IsZombie())continue;
    f->Open();
    if(fDebug>=4)cout << "Set source 0 for file " << chEl_sig->GetTitle() << endl;
    f->SetSource(0);
    //G4SBSRunData *rd;
    
    while( nevent<max_events ) {
      if(nevent%100==0)cout << nevent << " / " << max_events << endl;
      if(fDebug>=3)cout << "clear event " << endl;
      if(fDebug>=1)cout << "Process event " << nevent << endl;
      fEvent->Clear();
      has_data = false;
      // Accumulate data here...
      //f->
      
      // if(!f->ReadNextEvent(fDebug)){
      // 	if(nevt_b==0)nfile = -1;
      // 	//cout << "youhoo" << endl;
      // 	break;
      // }
      //while( f->ReadNextEvent(fDebug) ){
      if(!f->ReadNextEvent(fDebug))continue;
 
      for(size_t det = 0; det < fDetectors.size(); det++) {
	if(fDebug>=2){
	  cout << "load event " << f->GetEvNum() << " for file " << f->GetName() 
	       << " det " << fDetectors[det]->GetName() << endl;
	}
	//"LoadEventData" for signal - we want everything cleaned up for signal
	if(fDebug>=3)cout << "f->GetDataVector().size() " << f->GetDataVector().size() << endl;
	fDetectors[det]->SetTimeZero(0.);
	fDetectors[det]->LoadEventData(f->GetDataVector());
	/*
	//"LoadAccumulateData" for any other stuff we want to superimpose to signal
	if(fDebug>=3)cout << "f->GetDataVector().size() " << f->GetDataVector().size() << endl;
	double t0 = fRN->Uniform(-fManager->GetBkgdSpreadTimeWindowHW(), 
				 fManager->GetBkgdSpreadTimeWindowHW());
	fDetectors[det]->SetTimeZero(t0);
	fDetectors[det]->LoadAccumulateData(f->GetDataVector());
	*/
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

	    AddFileToEvent(f_b);
	    //chEl_bkgd
	    //fSourceChainMap[source]->RecursiveRemove(f_b);
	    f_b->~TSBSGeant4File();
	    nfile++;
	  }
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
	      double t0 = fRN->Uniform(-fManager->GetBkgdSpreadTimeWindowHW(), 
				       fManager->GetBkgdSpreadTimeWindowHW());
	      fDetectors[det]->SetTimeZero(t0);
	      fDetectors[det]->LoadAccumulateData(f->GetDataVector());
	      if(fDebug>=3)cout << "Done loading data in " << fDetectors[det]->GetName() << endl;
	    }
	    //AddFileToEvent(f_b);
	    //fSourceChainMap[source]->RecursiveRemove(f_b);
	  }
	}
	/*
	  TObjArray *BkgdFileList = fSourceChainMap[source]->GetListOfFiles();
	  TIter next_bkgd(BkgdFileList);
	  
	  TChainElement *chEl_bkgd = 0;
	*/
      }
      
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
