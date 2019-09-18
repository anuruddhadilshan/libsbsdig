#include "TSBSGeant4File.h"
#include "TSBSSimAuxi.h"
#include "TSBSSimData.h"
#include "g4sbs_types.h"
#include "fstream"

#ifndef __CINT__

// Set following variables to 1 (and recompile) t get some useful printouts
#ifndef D_FLAG
#define D_FLAG 0 //0: nothing; 1: warning; 2: debug;
#endif

TSBSGeant4File::TSBSGeant4File() : fFile(0), fSource(0) {
  fFilename[0] = '\0';
}

TSBSGeant4File::TSBSGeant4File(const char *f) : fFile(0), fSource(0) {
  //TSBSGeant4File::TSBSGeant4File(const char *f) : fFile(0), fSource(0) {
  SetFilename(f);
  fManager = TSBSDBManager::GetInstance();
  fRN = TRndmManager::GetInstance();
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

#if D_FLAG>1 
    cout << "Been there, done that" << endl;
#endif
    
    fTree = new g4sbs_tree(C1, fManager->GetExpType());
    // g4sbs_tree declare all variables, branches, etc... 
    // to read, event by event, the varaibles stored in the tree. 
    // See comments in g4sbs_tree for more details...
    
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
    
  bool res = false;
  if(d_flag>1){
    printf("Reading event %d\n", fEvNum);
  }
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

  //Useful variables for the processing
  int Npe;
  double t;
  
  /*
  // and here a bunch of variables that had a similar use, but lost it.
  // int n_hits = 0;//total number of hits at the end of the event
  // int n_gen = 0;//total number of tracks at the end of the event
  // bool newtrk, dupli;// These variables help avoid store many times the same MC track info
  // double weight = fTree->ev_solang*fTree->ev_sigma; 
  
  //int det_id;//2017/02/09: now corresponds to fManager->Getg4sbsDetType()
  
  // int PMTID;
  // double XPMT;
  // double YPMT;
  // double trms;
  // int type;
  // int PID_MCtrack;
  // double pz_MCtrack;
  // TVector3 Mom_MCtrack;
  // TVector3 Vtx_MCtrack;
  // TVector3 Pos_det;
  // int origvolflag;
  
  // double hit_data_temp[19];
  // double gen_data_temp[9];
  */

  // Electron Arm
  
  // Process GRINCH data
  
  if(fTree->Earm_GRINCH.nhits){
    for(int i = 0; i<fTree->Earm_GRINCH.nhits; i++){
      g4sbshitdata *grinchhit = new g4sbshitdata(GRINCH_UNIQUE_DETID, 5);
      grinchhit->SetData(0, fSource);
      grinchhit->SetData(1, int(fTree->Earm_GRINCH.PMT->at(i)/5));
      grinchhit->SetData(2, 0);
      grinchhit->SetData(3, fTree->Earm_GRINCH.Time_avg->at(i));
      grinchhit->SetData(4, fTree->Earm_GRINCH.NumPhotoelectrons->at(i));
      fg4sbsHitData.push_back(grinchhit);
    }
  }
    
  // Process hodoscope data
  if(fTree->Earm_BBHodoScint.nhits){
    for(int i = 0; i<fTree->Earm_BBHodoScint.nhits; i++){
      for(int j = 0; j<2; j++){//j = 0: close beam PMT, j = 1: far beam PMT
	g4sbshitdata *hodoscinthit = new g4sbshitdata(HODO_UNIQUE_DETID, 5);
	hodoscinthit->SetData(0, fSource);
	hodoscinthit->SetData(1, fTree->Earm_BBHodoScint.cell->at(i)*2+j);
	hodoscinthit->SetData(2, 1);
	hodoscinthit->SetData(3, fTree->Earm_BBHodoScint.tavg->at(i));
	hodoscinthit->SetData(4, fTree->Earm_BBHodoScint.sumedep->at(i));
	fg4sbsHitData.push_back(hodoscinthit);
	
	Npe = fRN->Poisson(1.0e7*fTree->Earm_BBHodoScint.sumedep->at(i)*0.113187*exp(-(0.3+pow(-1, j)*fTree->Earm_BBHodoScint.xhit->at(i))/1.03533)* 0.24);
	t = fTree->Earm_BBHodoScint.tavg->at(i)+(0.55+pow(-1, j)*fTree->Earm_BBHodoScint.xhit->at(i))/0.15;
	g4sbshitdata *hodopmthit = new g4sbshitdata(HODO_UNIQUE_DETID, 5);
	hodopmthit->SetData(0, fSource);
	hodopmthit->SetData(1, fTree->Earm_BBHodoScint.cell->at(i)*2+j);
	hodopmthit->SetData(2, 0);
	hodopmthit->SetData(3, t);
	hodopmthit->SetData(4, Npe);
	fg4sbsHitData.push_back(hodopmthit);
      }
    }
  }

  // Process BB PS data
  if(fTree->Earm_BBPSTF1.nhits){
    for(int i = 0; i<fTree->Earm_BBPSTF1.nhits; i++){
      g4sbshitdata *bbpstf1hit = new g4sbshitdata(BBPS_UNIQUE_DETID, 5);
      bbpstf1hit->SetData(0, fSource);
      bbpstf1hit->SetData(1, fTree->Earm_BBPSTF1.cell->at(i));
      bbpstf1hit->SetData(2, 1);
      bbpstf1hit->SetData(3, fTree->Earm_BBPSTF1.tavg->at(i));
      bbpstf1hit->SetData(4, fTree->Earm_BBPSTF1.sumedep->at(i));
      fg4sbsHitData.push_back(bbpstf1hit);
      
      Npe = fRN->Poisson(454.0*fTree->Earm_BBPSTF1.sumedep->at(i));
      t = fRN->Gaus(3.2-5.805*fTree->Earm_BBPSTF1.zhit->at(i)-17.77*pow(fTree->Earm_BBPSTF1.zhit->at(i), 2), 0.5);
      g4sbshitdata *bbpshit = new g4sbshitdata(BBPS_UNIQUE_DETID, 5);
      bbpshit->SetData(0, fSource);
      bbpshit->SetData(1, fTree->Earm_BBPSTF1.cell->at(i));
      bbpshit->SetData(2, 0);
      bbpshit->SetData(3, t);
      bbpshit->SetData(4, Npe);
      fg4sbsHitData.push_back(bbpshit);
    }
  }
  /*
  if(fTree->Earm_BBPS.nhits){
    for(int i = 0; i<fTree->Earm_BBPS.nhits; i++){
      g4sbshitdata *bbpshit = new g4sbshitdata(BBPS_UNIQUE_DETID, 5);
      bbpshit->SetData(0, fSource);
      bbpshit->SetData(1, fTree->Earm_BBPS.PMT->at(i));
      bbpshit->SetData(2, 0);
      bbpshit->SetData(3, fTree->Earm_BBPS.Time_avg->at(i));
      bbpshit->SetData(4, fTree->Earm_BBPS.NumPhotoelectrons->at(i));
      fg4sbsHitData.push_back(bbpshit);
    }
  }
  */
  
  // Process BB SH data
  if(fTree->Earm_BBSHTF1.nhits){
    for(int i = 0; i<fTree->Earm_BBSHTF1.nhits; i++){
      g4sbshitdata *bbshtf1hit = new g4sbshitdata(BBSH_UNIQUE_DETID, 5);
      bbshtf1hit->SetData(0, fSource);
      bbshtf1hit->SetData(1, fTree->Earm_BBSHTF1.cell->at(i));
      bbshtf1hit->SetData(2, 1);
      bbshtf1hit->SetData(3, fTree->Earm_BBSHTF1.tavg->at(i));
      bbshtf1hit->SetData(4, fTree->Earm_BBSHTF1.sumedep->at(i));
      fg4sbsHitData.push_back(bbshtf1hit);
      
      Npe = fRN->Poisson(932.0*fTree->Earm_BBSHTF1.sumedep->at(i));
      t = fRN->Gaus(2.216-8.601*fTree->Earm_BBSHTF1.zhit->at(i)-7.469*pow(fTree->Earm_BBSHTF1.zhit->at(i), 2), 0.8);
      g4sbshitdata *bbshhit = new g4sbshitdata(BBSH_UNIQUE_DETID, 5);
      bbshhit->SetData(0, fSource);
      bbshhit->SetData(1, fTree->Earm_BBSHTF1.cell->at(i));
      bbshhit->SetData(2, 0);
      bbshhit->SetData(3, t);
      bbshhit->SetData(4, Npe);
      fg4sbsHitData.push_back(bbshhit);
     }
  }
  /*
  if(fTree->Earm_BBSH.nhits){
    for(int i = 0; i<fTree->Earm_BBSH.nhits; i++){
      g4sbshitdata *bbshhit = new g4sbshitdata(BBSH_UNIQUE_DETID, 5);
      bbshhit->SetData(0, fSource);
      bbshhit->SetData(1, fTree->Earm_BBSH.PMT->at(i));
      bbshhit->SetData(2, 0);
      bbshhit->SetData(3, fTree->Earm_BBSH.Time_avg->at(i));
      bbshhit->SetData(4, fTree->Earm_BBSH.NumPhotoelectrons->at(i));
      fg4sbsHitData.push_back(bbshhit);
    }
  }
  */
  
  // Hadron Arm
  // Process CDet data
  if(fTree->Harm_CDET_Scint.nhits){
    for(int i = 0; i<fTree->Harm_CDET_Scint.nhits; i++){
      g4sbshitdata *cdetscinthit = new g4sbshitdata(CDET_UNIQUE_DETID, 5);
      cdetscinthit->SetData(0, fSource);
      cdetscinthit->SetData(1, fTree->Harm_CDET_Scint.cell->at(i));
      cdetscinthit->SetData(2, 1);
      cdetscinthit->SetData(3, fTree->Harm_CDET_Scint.tavg->at(i));
      cdetscinthit->SetData(4, fTree->Harm_CDET_Scint.sumedep->at(i));
      fg4sbsHitData.push_back(cdetscinthit);
      
      Npe = fRN->Poisson( fTree->Harm_CDET_Scint.sumedep->at(i)*5.634e3 );
      t = fTree->Harm_CDET_Scint.tavg->at(i)+5.75+fTree->Harm_CDET_Scint.xhit->at(i)/0.16;
      g4sbshitdata *cdetpmthit = new g4sbshitdata(CDET_UNIQUE_DETID, 5);
      cdetpmthit->SetData(0, fSource);
      cdetpmthit->SetData(1, fTree->Harm_CDET_Scint.cell->at(i));
      cdetpmthit->SetData(2, 0);
      cdetpmthit->SetData(3, t);
      cdetpmthit->SetData(4, Npe);
      fg4sbsHitData.push_back(cdetpmthit);
    }
  }

  // Now process the GEM data
  for(std::vector<gem_branch>::iterator it = fTree->GEMs.begin();
       it != fTree->GEMs.end(); it++) {

    TSBSGeant4::GEMData_t &t = it->tree;
    if(t.plane) {
      for(int k = 0; k < t.nhits; k++) {
        // Don't bother with events that deposited no energy
        if(t.edep->at(k)>0) {
          g4sbshitdata *gemhit = new g4sbshitdata(it->id, 35);
          // For now, just copy the whole tree, we'll trim it later I guess
          // (with slight modifications as prescribed in libsbsgem)
          gemhit->SetData(0,fSource);
          gemhit->SetData(1,t.plane->at(k));
          gemhit->SetData(2,t.strip->at(k));
          gemhit->SetData(3,t.x->at(k));
          gemhit->SetData(4,t.y->at(k));
          gemhit->SetData(5,t.z->at(k));
          gemhit->SetData(6,t.polx->at(k));
          gemhit->SetData(7,t.poly->at(k));
          gemhit->SetData(8,t.polz->at(k));
          gemhit->SetData(9,t.t->at(k));
          gemhit->SetData(10,t.trms->at(k));
          gemhit->SetData(11,t.tmin->at(k));
          gemhit->SetData(12,t.tmax->at(k));
          gemhit->SetData(13,t.tx->at(k));
          gemhit->SetData(14,t.ty->at(k));
          gemhit->SetData(15,t.txp->at(k));
          gemhit->SetData(16,t.typ->at(k));
          gemhit->SetData(17,t.xg->at(k));
          gemhit->SetData(18,t.yg->at(k));
          gemhit->SetData(19,t.zg->at(k));
          gemhit->SetData(20,t.trid->at(k));
          gemhit->SetData(21,t.mid->at(k)+1);
          gemhit->SetData(22,t.pid->at(k));
          gemhit->SetData(23,t.vx->at(k));
          gemhit->SetData(24,t.vy->at(k));
          gemhit->SetData(25,t.vz->at(k));
          gemhit->SetData(26,t.p->at(k));
          gemhit->SetData(27,t.edep->at(k)); // convert to MeV?
          gemhit->SetData(28,t.beta->at(k));
          gemhit->SetData(29,t.xin->at(k));
          gemhit->SetData(30,t.yin->at(k));
          gemhit->SetData(31,t.zin->at(k));
          gemhit->SetData(32,t.xout->at(k));
          gemhit->SetData(33,t.yout->at(k));
          gemhit->SetData(34,t.zout->at(k));
          fg4sbsHitData.push_back(gemhit);
        }
      }
    }
  }
  /*
  // For the time being, use the g4sbs npe estimation for CDET, divided by 5.
  // (the p.e. yield for cosmics form g4sbs is about 5 times more the one measured)
  if(fTree->Harm_CDET.nhits){
    for(int i = 0; i<fTree->Harm_CDET.nhits; i++){
      g4sbshitdata *cdetpmthit = new g4sbshitdata(CDET_UNIQUE_DETID, 5);
      cdetpmthit->SetData(0, fSource);
      cdetpmthit->SetData(1, fTree->Harm_CDET.PMT->at(i));
      cdetpmthit->SetData(2, 0);
      cdetpmthit->SetData(3, fTree->Harm_CDET.Time_avg->at(i));
      cdetpmthit->SetData(4, round(fTree->Harm_CDET.NumPhotoelectrons->at(i)/5.0) );
      fg4sbsHitData.push_back(cdetpmthit);
    }
  }
  */
  
  // Now process HCAL data
  if(fTree->hcalpart.E) {
    for(size_t k = 0; k < fTree->hcalpart.E->size(); k++) {
      if(fTree->hcalpart.detected->at(k)) {
        g4sbshitdata *hcalpmthit = new g4sbshitdata(HCAL_UNIQUE_DETID,4);
        hcalpmthit->SetData(0,fSource);
        hcalpmthit->SetData(1,fTree->hcalpart.part_PMT->at(k));
        hcalpmthit->SetData(2,0); // data type 0 == optical photon
        hcalpmthit->SetData(3,fTree->hcalpart.t->at(k));
        fg4sbsHitData.push_back(hcalpmthit);
      }
    }
  }
  // For now, get the adc signal from the total energy deposited on the
  // scintillators. This can be changed later, I suppose...
  if(fTree->hcalscint.sumedep) {
    for(size_t k = 0; k < fTree->hcalscint.sumedep->size(); k++) {
      g4sbshitdata *hcalscinthit = new g4sbshitdata(HCAL_UNIQUE_DETID,4);
      hcalscinthit->SetData(0,fSource);
      hcalscinthit->SetData(1,fTree->hcalscint.cell->at(k));
      hcalscinthit->SetData(2,1); // data type 1 == sumedep
      hcalscinthit->SetData(3,fTree->hcalscint.sumedep->at(k));
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
