#include "TSBSSimDigitizer.h"
#include "TSBSGeant4File.h"
#include "TSBSSimDetector.h"
#include "TSBSSimEvent.h"
#include "TSBSSimAuxi.h"
#include "TSBSDBManager.h"
#include <TTree.h>
#include <TFile.h>

TSBSSimDigitizer::TSBSSimDigitizer(const char* outputfilename)
{
  if(fDebug>=1)cout << "Initialize TSBSSimDigitzer " << endl;
  fManager = TSBSDBManager::GetInstance();
  
  fEvent = new TSBSSimEvent();
  fOutFile = new TFile(outputfilename,"RECREATE");
  //fOutTree = new TTree("TSBSDigi","");
  fOutTree = new TTree("digtree","");
  //fOutTree->Branch("fEvents",&fEvent);
  fOutTree->Branch("event",&fEvent);
  fOutTree->Branch("RunID",&fEvent->fRunID);
  fOutTree->Branch("EvtID",&fEvent->fEvtID);
  fOutTree->Branch("Weight",&fEvent->fWeight);
  fOutTree->Branch("NSignal",&fEvent->fNSignal);

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
    fOutTree->Branch(Form("SimData_%s_Chan", fulldetname.c_str()),&fEvent->SimDetChannel[fulldetname.c_str()]);
    fOutTree->Branch(Form("SimData_%s_Type", fulldetname.c_str()),&fEvent->SimDetDataType[fulldetname.c_str()]);
    fOutTree->Branch(Form("SimData_%s_Ndata", fulldetname.c_str()),&fEvent->SimDetNData[fulldetname.c_str()]);
    fOutTree->Branch(Form("SimData_%s_Data", fulldetname.c_str()),&fEvent->SimDetData[fulldetname.c_str()]);
    
    fOutTree->Branch(Form("Data_%s_Chan", fulldetname.c_str()),&fEvent->DetChannel[fulldetname.c_str()]);
    fOutTree->Branch(Form("Data_%s_Ndata", fulldetname.c_str()),&fEvent->DetNData[fulldetname.c_str()]);
    fOutTree->Branch(Form("Data_%s_Data", fulldetname.c_str()),&fEvent->DetData[fulldetname.c_str()]);
  }
  
  //const std::vector<TSpectroInfo> = fManager->GetAllSpectroInfo();
  //std::vector<TDetInfo> = fManager->GetAllDetInfo();
  
  //const TSpectroInfo SpecInfo = fManager->GetAllSpectroInfo();
  //for()
  //fOutTree->Branch("SimDetectorData",&fEvent->fSimDetectorData);
  //fOutTree->Branch("DetectorData",&fEvent->fSimDetectorData);
}

TSBSSimDigitizer::~TSBSSimDigitizer()
{
  fDetectors.clear();
  fG4FileStack_.clear();
}

int TSBSSimDigitizer::Process(TSBSGeant4File *f, int max_events)
{
  if(!f)
    return 0;

  // Can we open the file?
  f->SetSource(0);
  int res = f->Open();
  if( res != 1) {
    std::cerr << "Failed to open g4sbs rootfile. Failed with error code: "
      << res << std::endl;
    return 0;
  }

  if ( max_events <= 0 || max_events > f->GetEntries() )
    max_events = f->GetEntries();

  // Now loop through the file and digitize entries
  //int d_flag_readevent = 0;
  int nevent = 0;
  //int ngood = 0;
  bool has_data;
  while( f->ReadNextEvent(fDebug) && nevent<max_events ) {
    if(fDebug>=3)cout << "clear event " << endl;
    fEvent->Clear();
    has_data = false;
    // Tell detectors a new event has started
    for(size_t det = 0; det < fDetectors.size(); det++) {
      if(fDebug>=3)cout << "event start det " << fDetectors[det]->GetName() << endl;
      fDetectors[det]->EventStart();
    }
   
    TSBSSimDetector::SetEventNum(nevent); ///< Needed by GEMs for some reason
    // Loop through all detectors and have them parse data vector
    for(size_t det = 0; det < fDetectors.size(); det++) {
      if(fDebug>=3){
	cout << "load event for det " << fDetectors[det]->GetName() << endl;
	cout << "f->GetDataVector().size() " << f->GetDataVector().size() << endl;
      }
      fDetectors[det]->LoadEventData(f->GetDataVector());
    }
    // Now digitize all detectors
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
      fOutTree->Fill();
      if(fDebug>=1)std::cout << "Have data for event: " << nevent << std::endl;
    }
    // Tell detectors the event ended
    for(size_t det = 0; det < fDetectors.size(); det++) {
      if(fDebug>=3)cout << "event end det " << fDetectors[det]->GetName() << endl;
      fDetectors[det]->EventEnd();
    }


    nevent++;
  }
  
  if(fDebug>=1)cout << "write output file " << endl;
  fOutFile->Write();
  // Close files
  
  //cout << "close output file " << endl;
  //fOutFile->Close();
  
  //cout << "close geeant 4 file " << endl;
  //f->Close();
  
  if(fDebug>=1)cout << "Done closing files " << endl;

  return 0;
}

int TSBSSimDigitizer::Process(int max_events)
{
  cout << "Warning:  TSBSSimDigitizer::Process(int) is not functional yet." << endl 
       << "Please use int TSBSSimDigitizer::Process(TSBSGeant4File, int)" << endl;
  if(fG4FileStack_.size()==0)
    return 0;
  
  //Loop on files:
  //TSBSGeant4File* f; = fG4FileStack.at(0);
  //TSBSGeant4File* f_b;
  // Can we open the files?
  /*
  for(int i_f = 0; i_f<fG4FileStack.size(); i_f++){
    for(int j_f = 0; j_f<fG4FileStack.at(i_f).size(); j_f++){
      int res = fG4FileStack_.at(i_f)->Open();
      if( res != 1) {
	std::cerr << "Failed to open g4sbs rootfile " << std::endl
		  << fG4FileStack_.at(i_f)->GetFileName() << std::endl 
		  << "Failed with error code: " << res << std::endl;
	return 0;
      }
    }
  }
  */
  //go through the file stack and open all of them...
  // determine which is primary:
  Double_t PrimWeight = 0;
  for(size_t i_f = 0; i_f<fG4FileStack_.size(); i_f++){
    int res = fG4FileStack_.at(i_f)->Open();
    if( res != 1) {
      std::cerr << "Failed to open g4sbs rootfile " << std::endl
		<< fG4FileStack_.at(i_f)->GetFileName() << std::endl 
		<< "Failed with error code: " << res << std::endl;
      return 0;
    }
    if(fG4FileStack_.at(i_f)->GetSource()==0)PrimWeight = fG4FileWeights.at(i_f);
  }
  
  if ( max_events <= 0 || max_events > fG4FileStack_.at(0)->GetEntries() )
    max_events = fG4FileStack_.at(0)->GetEntries();
  
  // Now loop through the file and digitize entries
  //int d_flag_readevent = 0;
  int nevent = 0;
  int i_f = 0;
  UInt_t nevt_b;
  //int ngood = 0;
  bool has_data;
  
  while( nevent<max_events ) {
    if(fDebug>=3)cout << "clear event " << endl;
    fEvent->Clear();
    has_data = false;
    i_f = 0;
    // Accumulate data here...
    //for(size_t i_f = 0; i_f<fG4FileStack_.size(); i_f++){
    for(std::vector<TSBSGeant4File*>::const_iterator it = fG4FileStack_.begin(); it!=fG4FileStack_.end(); ++it){
      TSBSGeant4File* f = (*it);
      cout << " i_f "  << i_f << endl;
      nevt_b = 0;
      //f_b = fG4FileStack.at(i_f);
      //if(fG4FileWeights.at(i_f)>=0){}
      
      // while( f->ReadNextEvent(fDebug) && 
      // 	     nevt_b<(fG4FileWeights.at(i_f)/PrimWeight) ) {//Keep adding as many events as indicated by the weight
      while( f->ReadNextEvent(fDebug) ) {
	if(fG4FileWeights.at(i_f)>=0 && 
	   nevt_b>=(fG4FileWeights.at(i_f)/PrimWeight) )break;
	// Loop through all detectors and have them parse data vector
	for(size_t det = 0; det < fDetectors.size(); det++) {
	  if(fDebug>=3){
	    cout << "load event " << nevt_b << " for file " << f->GetFileName() 
		 << " det " << fDetectors[det]->GetName() << endl;
	  }
	  if(f->GetSource()==0){
	    //"LoadEventData" for signal - we want teverything cleanded up for signal
	    if(fDebug>=3)cout << "f->GetDataVector().size() " << f->GetDataVector().size() << endl;
	    fDetectors[det]->LoadEventData(f->GetDataVector());
	  }else{
	    //"LoadAccumulateData" for any other stuff we want to superimposeto signal
	    fDetectors[det]->LoadAccumulateData(f->GetDataVector());
	  }
	}
	nevt_b++;
      }
      i_f++;
    }
    
    // Now digitize all detectors
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
      fOutTree->Fill();
      if(fDebug>=1)std::cout << "Have data for event: " << nevent << std::endl;
    }
    
    nevent++;
  }
  
  if(fDebug>=1)cout << "write output file " << endl;
  fOutFile->Write();
  // Close files
  
  //cout << "close output file " << endl;
  //fOutFile->Close();
  
  //cout << "close geeant 4 file " << endl;
  //f->Close();
  
  if(fDebug>=1)cout << "Done closing files " << endl;

  return 0;
}


void TSBSSimDigitizer::AddDetector(TSBSSimDetector* detector)
{
  fDetectors.push_back(detector);
}

void TSBSSimDigitizer::AddInputFile(TSBSGeant4File* file, UInt_t weight)
{
  //fG4FileStack.push_back(std::make_pair<file, weight>);
  // Files and weights are added at the same time, and cannot be accessed from the outside... 
  // No need to inforce them to be bound together by a quirky object
  fG4FileStack_.push_back(file);
  fG4FileWeights.push_back(weight);
  
}

ClassImp(TSBSSimDigitizer)
