#include "TSBSGeant4File.h"
#include "TSBSSimData.h"
#include "g4sbs_types.h"
#include "fstream"

#ifndef __CINT__

// Set following variables to 1 (and recompile) t get some useful printouts
#ifndef D_FLAG
#define D_FLAG 1 //0: nothing; 1: warning; 2: debug;
#endif

//List of detector unique IDs: 
// by (proposed) convention: DetUniqueID = DetType*10+DetID
// DetType of type det_type defined in g4sbs_types: kHCal(0), kECal(1), kCher(2), kScint(3), kGEM(4);
#define HCAL_UNIQUE_DETID 0
#define BBPS_UNIQUE_DETID 10
#define BBSH_UNIQUE_DETID 11
#define ECAL_UNIQUE_DETID 12
#define GRINCH_UNIQUE_DETID 20
#define RICH_UNIQUE_DETID 21
#define HODO_UNIQUE_DETID 30
#define CDET_UNIQUE_DETID 31
#define BBGEM_UNIQUE_DETID 40
#define SBSGEM_UNIQUE_DETID 41
#define FT_UNIQUE_DETID 42
#define FPP1_UNIQUE_DETID 43
#define FPP2_UNIQUE_DETID 44
//#define CHER_HIT_ID 0

TSBSGeant4File::TSBSGeant4File() : fFile(0), fSource(0) {
  fFilename[0] = '\0';
}

TSBSGeant4File::TSBSGeant4File(const char *f) : fFile(0), fSource(0) {
  //TSBSGeant4File::TSBSGeant4File(const char *f) : fFile(0), fSource(0) {
  SetFilename(f);
  fManager = TSBSDBManager::GetInstance();
  
  //Filling the table that will be used to calculate the low energy electron range in the gas. 
#if D_FLAG>1
  cout << "Initialization completed" << endl;
#endif
}

TSBSGeant4File::~TSBSGeant4File() {
  Clear();
  delete fFile;
}

void TSBSGeant4File::SetFilename( const char *f ){
  if( !f ) return;
  strcpy( fFilename, f );
}

Int_t TSBSGeant4File::Open(){
    // Return 0 on fail, 1 on success
    if( fFilename[0] == '\0' ){ return 0; }

    delete fFile;
    fFile = new TFile(fFilename);
    
    if( !fFile->IsOpen() ){ 
      fprintf(stderr, "%s: File could not be made\n",__PRETTY_FUNCTION__);
      return 0; 
    }
    
    TChain* C1 = (TChain*)fFile->Get("T");//Get the tree from the file
    
    /*
#if D_FLAG>1 
    cout << "Detector option " << fManager->Getg4sbsDetType() << endl;
#endif
    
    // Detector flag. See the printout content below.
    if(fManager->Getg4sbsDetType()<1 && fManager->Getg4sbsDetType()>2){
      cout << "Invalid detector option: Set correct option in db_generalinfo.dat" << endl;
      cout << "(remider: 1 - GRINCH; 2 - RICH)" << endl;
      return 0;
    }
    
    fTree = new g4sbs_tree(C1, fManager->Getg4sbsExpType());
    // g4sbs_tree declare all variables, branches, etc... 
    // to read, event by event, the varaibles stored in the tree. 
    // See comments in g4sbs_tree for more details...
    */
    
    fEvNum = -1;
 
    return 1;
}

Int_t TSBSGeant4File::Close(){
    // Return 0 on fail, 1 on success
    Int_t ret = 1;
    
    if( !fFile->IsOpen() ){ return 0; }
    
    fFile->Close();
    
    delete fFile; fFile = 0;
    return ret;
}

Int_t TSBSGeant4File::ReadNextEvent(int d_flag){
  // Return 1 on success
  
  // Channel not open
  if( !fFile->IsOpen() ){ 
    fprintf(stderr, "%s %s line %d Channel not open\n",
	    __FILE__,__PRETTY_FUNCTION__,__LINE__ );
    return 0; 
  }
  
  Clear();
    
  int n_hits = 0;//total number of hits at the end of the event
  int n_gen = 0;//total number of tracks at the end of the event
  // bool newtrk, dupli;// These variables help avoid store many times the same MC track info
  bool res = false;
  
  fEvNum++;

  //cout << "Read Next Event: Evt " << fEvNum << endl;
  
  res = fTree->GetEntry(fEvNum);
  //Test that the next entry exist
  if( !res ){
    // Don't need to print this out.  Not really an error
    if(d_flag>1){
      fprintf(stderr, "%s %s line %d: Channel read return is false...  probably end of file\n",
	      __FILE__, __FUNCTION__, __LINE__ );
    } //DEBUG
    return 0;
  }
  
  double weight = fTree->ev_solang*fTree->ev_sigma; 
    
  int det_id;//2017/02/09: now corresponds to fManager->Getg4sbsDetType()
  
  int PMTID;
  double XPMT;
  double YPMT;
  int Npe;
  double t;
  double trms;
  int type;
  int PID_MCtrack;
  double pz_MCtrack;
  TVector3 Mom_MCtrack;
  TVector3 Vtx_MCtrack;
  TVector3 Pos_det;
  int origvolflag;
  double nph = 0;
  
  double hit_data_temp[19];
  double gen_data_temp[9];
  
  // Process Scint data
  if(fTree->Earm_BBHodoScint.nhits){
    for(int i = 0; i<fTree->Earm_BBHodoScint.nhits; i++){
      g4sbshitdata *hodoscinthit = new g4sbshitdata(HODO_UNIQUE_DETID, 3);
      hodoscinthit->SetData(0, fTree->Earm_BBHodoScint.cell->at(i));
      hodoscinthit->SetData(1, 1);
      hodoscinthit->SetData(2, fTree->Earm_BBHodoScint.sumedep->at(i));
      fg4sbsHitData.push_back(hodoscinthit);
      
      for(int j = 0; j<2; j++){//j = 0: close beam PMT, j = 1: far beam PMT
	nph = 1.0e7*fTree->Earm_BBHodoScint.sumedep->at(i)*0.113187*exp(-(0.3+pow(-1, j)*fTree->Earm_BBHodoScint.xhit->at(i))/1.03533);
	g4sbshitdata *hodopmthit = new g4sbshitdata(HODO_UNIQUE_DETID, 3);
	hodopmthit->SetData(0, fTree->Earm_BBHodoScint.cell->at(i)*2+j);
	hodopmthit->SetData(1, 0);
	hodopmthit->SetData(2, nph*0.24);
	fg4sbsHitData.push_back(hodopmthit);
      }
    }
  }


  // Now process HCAL data
  if(fTree->hcalpart.E) {
    for(size_t k = 0; k < fTree->hcalpart.E->size(); k++) {
      if(fTree->hcalpart.detected->at(k)) {
        g4sbshitdata *hcalpmthit = new g4sbshitdata(HCAL_UNIQUE_DETID,3);
        hcalpmthit->SetData(0,fTree->hcalpart.part_PMT->at(k));
        hcalpmthit->SetData(1,0); // data type 0 == optical photon
        hcalpmthit->SetData(2,fTree->hcalpart.t->at(k));
        fg4sbsHitData.push_back(hcalpmthit);
      }
    }
  }
  // For now, get the adc signal from the total energy deposited on the
  // scintillators. This can be changed later, I suppose...
  if(fTree->hcalscint.sumedep) {
    for(size_t k = 0; k < fTree->hcalscint.sumedep->size(); k++) {
      g4sbshitdata *hcalscinthit = new g4sbshitdata(HCAL_UNIQUE_DETID,3);
      hcalscinthit->SetData(0,fTree->hcalscint.cell->at(k));
      hcalscinthit->SetData(1,1); // data type 1 == sumedep
      hcalscinthit->SetData(2,fTree->hcalscint.sumedep->at(k));
      fg4sbsHitData.push_back(hcalscinthit);
    }
  }

  return 1;
}


void TSBSGeant4File::Clear(){
  // Clear out hit and generated data

#if D_FLAG>1
  fprintf(stderr, "%s %s line %d: Deleting hits\n",
	  __FILE__, __FUNCTION__, __LINE__);
#endif //DEBUG

  unsigned int i;
  for( i = 0; i < fg4sbsHitData.size(); i++ ){
    delete fg4sbsHitData[i];
  }

  for( i = 0; i < fg4sbsGenData.size(); i++ ){
    delete fg4sbsGenData[i];
  }

  fg4sbsHitData.clear();
  fg4sbsGenData.clear();

#if D_FLAG>1
  fprintf(stderr, "%s %s line %d: Hits deleted\n",
	  __FILE__, __FUNCTION__, __LINE__);
#endif //DEBUG

  return;
}

#endif//__CINT__
