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
  
  double hit_data_temp[19];
  double gen_data_temp[9];

  // Electron Arm
  
  // Process GRINCH data
  if(fTree->Earm_GRINCH.nhits){
    for(int i = 0; i<fTree->Earm_GRINCH.nhits; i++){
      g4sbshitdata *grinchhit = new g4sbshitdata(GRINCH_UNIQUE_DETID, 4);
      grinchhit->SetData(0, fTree->Earm_GRINCH.PMT->at(i));
      grinchhit->SetData(1, 1);
      grinchhit->SetData(2, fTree->Earm_GRINCH.Time_avg->at(i));
      grinchhit->SetData(3, fTree->Earm_GRINCH.NumPhotoelectrons->at(i));
      fg4sbsHitData.push_back(grinchhit);
    }
  }
    
  // Process hodoscope data
  if(fTree->Earm_BBHodoScint.nhits){
    for(int i = 0; i<fTree->Earm_BBHodoScint.nhits; i++){
      g4sbshitdata *hodoscinthit = new g4sbshitdata(HODO_UNIQUE_DETID, 3);
      hodoscinthit->SetData(0, fTree->Earm_BBHodoScint.cell->at(i));
      hodoscinthit->SetData(1, 1);
      hodoscinthit->SetData(2, fTree->Earm_BBHodoScint.sumedep->at(i));
      //hodoscinthit->SetData(3, fTree->Earm_BBHodoScint.tavg->at(i));
      fg4sbsHitData.push_back(hodoscinthit);
      
      for(int j = 0; j<2; j++){//j = 0: close beam PMT, j = 1: far beam PMT
	Npe = 1.0e7*fTree->Earm_BBHodoScint.sumedep->at(i)*0.113187*exp(-(0.3+pow(-1, j)*fTree->Earm_BBHodoScint.xhit->at(i))/1.03533)* 0.24;
	t = fTree->Earm_BBHodoScint.tavg->at(i)+(0.55+pow(-1, j)*fTree->Earm_BBHodoScint.xhit->at(i))/0.15;
	g4sbshitdata *hodopmthit = new g4sbshitdata(HODO_UNIQUE_DETID, 4);
	hodopmthit->SetData(0, fTree->Earm_BBHodoScint.cell->at(i)*2+j);
	hodopmthit->SetData(1, 0);
	hodopmthit->SetData(2, t);
	hodopmthit->SetData(3, Npe*0.24);
	fg4sbsHitData.push_back(hodopmthit);
      }
    }
  }

  // Process BB PS data
  if(fTree->Earm_BBPSTF1.nhits){
    for(int i = 0; i<fTree->Earm_BBPSTF1.nhits; i++){
      g4sbshitdata *bbpstf1hit = new g4sbshitdata(BBPS_UNIQUE_DETID, 3);
      bbpstf1hit->SetData(0, fTree->Earm_BBPSTF1.cell->at(i));
      bbpstf1hit->SetData(1, 1);
      bbpstf1hit->SetData(2, fTree->Earm_BBPSTF1.sumedep->at(i));
      //bbpstf1hit->SetData(3, fTree->Earm_BBPSTF1.tavg->at(i));
      fg4sbsHitData.push_back(bbpstf1hit);
    }
  }
  
  if(fTree->Earm_BBPS.nhits){
    for(int i = 0; i<fTree->Earm_BBPS.nhits; i++){
      g4sbshitdata *bbpshit = new g4sbshitdata(BBPS_UNIQUE_DETID, 4);
      bbpshit->SetData(0, fTree->Earm_BBPS.PMT->at(i));
      bbpshit->SetData(1, 0);
      bbpshit->SetData(2, fTree->Earm_BBPS.Time_avg->at(i));
      bbpshit->SetData(3, fTree->Earm_BBPS.NumPhotoelectrons->at(i));
      fg4sbsHitData.push_back(bbpshit);
    }
  }
  
  // Process BB SH data
  if(fTree->Earm_BBSHTF1.nhits){
    for(int i = 0; i<fTree->Earm_BBSHTF1.nhits; i++){
      g4sbshitdata *bbshtf1hit = new g4sbshitdata(BBSH_UNIQUE_DETID, 3);
      bbshtf1hit->SetData(0, fTree->Earm_BBSHTF1.cell->at(i));
      bbshtf1hit->SetData(1, 1);
      bbshtf1hit->SetData(2, fTree->Earm_BBSHTF1.sumedep->at(i));
      //bbshtf1hit->SetData(3, fTree->Earm_BBSHTF1.tavg->at(i));
      fg4sbsHitData.push_back(bbshtf1hit);
    }
  }
  
  if(fTree->Earm_BBSH.nhits){
    for(int i = 0; i<fTree->Earm_BBSH.nhits; i++){
      g4sbshitdata *bbshhit = new g4sbshitdata(BBSH_UNIQUE_DETID, 4);
      bbshhit->SetData(0, fTree->Earm_BBSH.PMT->at(i));
      bbshhit->SetData(1, 0);
      bbshhit->SetData(2, fTree->Earm_BBSH.Time_avg->at(i));
      bbshhit->SetData(3, fTree->Earm_BBSH.NumPhotoelectrons->at(i));
      fg4sbsHitData.push_back(bbshhit);
    }
  }
  
  
  // Hadron Arm
  // Process CDet data
  if(fTree->Harm_CDET_Scint.nhits){
    for(int i = 0; i<fTree->Harm_CDET_Scint.nhits; i++){
      g4sbshitdata *cdetscinthit = new g4sbshitdata(CDET_UNIQUE_DETID, 3);
      cdetscinthit->SetData(0, fTree->Harm_CDET_Scint.cell->at(i));
      cdetscinthit->SetData(1, 1);
      cdetscinthit->SetData(2, fTree->Harm_CDET_Scint.sumedep->at(i));
      //cdetscinthit->SetData(3, fTree->Harm_CDET_Scint.tavg->at(i));
      fg4sbsHitData.push_back(cdetscinthit);
    }
  }
  // For the time being, use the g4sbs npe estimation for CDET, divided by 5.
  // (the p.e. yield for cosmics form g4sbs is about 5 times more the one measured)
  if(fTree->Harm_CDET.nhits){
    for(int i = 0; i<fTree->Harm_CDET.nhits; i++){
      g4sbshitdata *cdetpmthit = new g4sbshitdata(CDET_UNIQUE_DETID, 4);
      cdetpmthit->SetData(0, fTree->Harm_CDET.PMT->at(i));
      cdetpmthit->SetData(1, 0);
      cdetpmthit->SetData(2, fTree->Harm_CDET.Time_avg->at(i));
      cdetpmthit->SetData(3, round(fTree->Harm_CDET.NumPhotoelectrons->at(i)/5.0) );
      fg4sbsHitData.push_back(cdetpmthit);
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
