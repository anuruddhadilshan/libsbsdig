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

TSBSGeant4File::TSBSGeant4File() : //fFile(0), 
  TFile::TFile(""), fTree(0), fSource(0) {
  //fFilename[0] = '\0';
}

TSBSGeant4File::TSBSGeant4File(const char *f) : //fFile(0), 
  TFile::TFile(f), fTree(0), fSource(0) {
  //TSBSGeant4File::TSBSGeant4File(const char *f) : fFile(0), fSource(0) {
  //SetFilename(f);
  //TFile::TFile(f)
  fManager = TSBSDBManager::GetInstance();
  fRN = TRndmManager::GetInstance();
  //Filling the table that will be used to calculate the low energy electron range in the gas. 
#if D_FLAG>1
  cout << "Initialization completed" << endl;
#endif
}

TSBSGeant4File::~TSBSGeant4File() {
  Clear();
  Delete();
  //delete fFile;
}

/* // 2019/10/18: TSBSGeant4File now inherits of TFile (EF)
void TSBSGeant4File::SetFilename( const char *f ){
  if( !f ) return;
  strcpy( fFilename, f );
}
*/

Int_t TSBSGeant4File::Open(){
    // Return 0 on fail, 1 on success
  if( GetName() == '\0' ){ return 0; }

  //delete fFile;
  //fFile = new TFile(fFilename);
    
  //if( !fFile->IsOpen() ){ 
  if( !IsOpen() ){ 
    fprintf(stderr, "%s: File could not be made\n",__PRETTY_FUNCTION__);
    return 0; 
  }
  
  TChain* C1 = (TChain*)fFile->Get("T");//Get the tree from the file
  
  fTree = new g4sbs_tree(C1, fManager->GetExpType());
  // g4sbs_tree declare all variables, branches, etc... 
  // to read, event by event, the varaibles stored in the tree. 
  // See comments in g4sbs_tree for more details...

  fThetaSBS = fdHCal = 0;
  fzGrinch = fzHodo = fzPS = fzSH = 0;
  
  if(fManager->IsDetInfoAvailable("hcal")){
    TSpectroInfo sbsinfo = fManager->GetSpectroInfo("sbs");
    fThetaSBS = sbsinfo.MCAngle();
    TDetInfo hcalinfo = fManager->GetDetInfo("hcal");
    fdHCal = hcalinfo.GeoInfo(0).ZPos();
  }
  if(fManager->IsDetInfoAvailable("grinch")){
    TDetInfo grinchinfo = fManager->GetDetInfo("grinch");
    fzGrinch = grinchinfo.GeoInfo(0).ZPos();
  }
  if(fManager->IsDetInfoAvailable("hodo")){
    TDetInfo hodoinfo = fManager->GetDetInfo("hodo");
    fzHodo = hodoinfo.GeoInfo(0).ZPos();
  }
  if(fManager->IsDetInfoAvailable("ps")){
    TDetInfo psinfo = fManager->GetDetInfo("ps");
    fzPS = psinfo.GeoInfo(0).ZPos();
  }
  if(fManager->IsDetInfoAvailable("sh")){
    TDetInfo shinfo = fManager->GetDetInfo("sh");
    fzSH = shinfo.GeoInfo(0).ZPos();
  }
  fEvNum = -1;
  
#if D_FLAG>1 
  cout << "Just opened file " << GetName() << endl;
#endif
  
  return 1;
}

/* // 2019/10/18: TSBSGeant4File now inherits of TFile (EF)
Int_t TSBSGeant4File::Close(){
  // Return 0 on fail, 1 on success
  Int_t ret = 1;
  
  if( !fFile->IsOpen() ){ return 0; }
  
  fFile->Close();
  
  delete fFile; fFile = 0;
  return ret;
}
*/

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
  if(d_flag>=3){
    printf("Reading event %d\n", fEvNum+1);
  }
  fEvNum++;

  res = fTree->GetEntry(fEvNum);

  //Test that the next entry exist
  if( !res ){
    // Don't need to print this out.  Not really an error
    if(d_flag>=2){
      fprintf(stderr, "%s %s line %d: Channel read return is false...  probably end of file\n",
	      __FILE__, __FUNCTION__, __LINE__ );
    } //DEBUG
    return 0;
  }

  //Useful variables for the processing
  int Npe;
  double t;
  
  double x_ref = fdHCal*sin(fThetaSBS);
  double z_ref = fdHCal*cos(fThetaSBS);
  
  double z_hit, Npe_Edep_ratio, sigma_tgen;

  double z_det, pz, E, beta;
  
  // Electron Arm
  
  // Process GRINCH data
  
  if(d_flag>=3)printf("Source  = %d\n", fSource);
  
  if(d_flag>=3)printf("Unfolding MC info \n");
  
  //if(fManager->IsDetInfoAvailable("hcal")){
  if(fSource==0 && fTree->hcalscint.sumedep){
    Double_t alpha = atan2(fTree->ev_npx,fTree->ev_npz);
    Double_t tanbeta = fTree->ev_npy/sqrt(fTree->ev_npx*fTree->ev_npx+fTree->ev_npz*fTree->ev_npz);
    Double_t v = fTree->ev_vz + fTree->ev_vx*fTree->ev_npx/fTree->ev_npz;
    Short_t PID = 2112;
    if(fTree->ev_nucl==1)PID = 2122;
    //SIDIS ?
    if(fManager->GetExpType()==kSIDIS){
      switch(TMath::Abs(fTree->ev_hadr)){
      case(0):
	PID = 111;
	break;
      case(1):
	PID = 211;
	break;
      case(2):
	PID = 311;
	break;
      default:
	PID = 2122;
	break;
      }
      if(fTree->ev_hadr<0)PID = -PID;
    }
    E = sqrt(fTree->ev_np*fTree->ev_np) + M_p(PID)*M_p(PID);
    beta = fTree->ev_np/E;
    g4sbsgendata *hcalgenhit = new g4sbsgendata(HCAL_UNIQUE_DETID, 8);
    hcalgenhit->SetData(0, fSource); 
    hcalgenhit->SetData(1, 2);
    hcalgenhit->SetData(2, PID);
    //Assumes NO bending (which is *wrong* for anything except neutrons)
    hcalgenhit->SetData(3, -(fdHCal-v*cos(fThetaSBS))*tanbeta/cos(alpha-fThetaSBS)-fTree->ev_vy); 
    hcalgenhit->SetData(4, (fdHCal-v*cos(fThetaSBS))*tan(alpha-fThetaSBS)-v*sin(fThetaSBS)); 
    hcalgenhit->SetData(5, fdHCal*fTree->ev_np/fTree->ev_npz/(beta*0.299792458) );
    hcalgenhit->SetData(6, E);
    hcalgenhit->SetData(7, 1.0);// TODO: weight
    fg4sbsGenData.push_back(hcalgenhit);
  }
  //}
  
  //was redoing that loop over and over again... what a waste... need to condense
  for(std::vector<gem_branch>::iterator it = fTree->GEMs.begin();
      it != fTree->GEMs.end(); it++) {
    if(d_flag>=3)printf("GEM tree %ld\n", std::distance(fTree->GEMs.begin(), it));
    
    TSBSGeant4::TrackerData_t &t = it->Track_tree;
    if(t.ntracks) {
      if(d_flag>=3)printf("%d tracks \n", t.ntracks);
      for(int k = 0; k < t.ntracks; k++) {
	//if(t.NumPlanes->at(k)<5)continue;
	
	pz = sqrt( t.P->at(k)*t.P->at(k) / ( t.Xp->at(k)*t.Xp->at(k) + t.Yp->at(k)*t.Yp->at(k) + 1.0) );
	E = sqrt( t.P->at(k)*t.P->at(k) + M_p(t.PID->at(k))*M_p(t.PID->at(k)) );
	beta = t.P->at(k)/E;
	
	g4sbsgendata *bbgemgentrack = new g4sbsgendata(BBGEM_UNIQUE_DETID, 16);
	bbgemgentrack->SetData(0, fSource);
	bbgemgentrack->SetData(1, t.TID->at(k));
	bbgemgentrack->SetData(2, t.PID->at(k));
	bbgemgentrack->SetData(3, t.X->at(k));
	bbgemgentrack->SetData(4, t.Y->at(k));
	bbgemgentrack->SetData(5, t.T->at(k));
	bbgemgentrack->SetData(6, t.P->at(k));
	bbgemgentrack->SetData(7, t.Xp->at(k));
	bbgemgentrack->SetData(8, t.Yp->at(k));
	if(fSource==0 && t.TID->at(k)==1){
	  bbgemgentrack->SetData(9, fTree->ev_vx);
	  bbgemgentrack->SetData(10, fTree->ev_vy);
	  bbgemgentrack->SetData(11, fTree->ev_vz);
	  bbgemgentrack->SetData(12, fTree->ev_epx);
	  bbgemgentrack->SetData(13, fTree->ev_epy);
	  bbgemgentrack->SetData(14, fTree->ev_epz);
	}else{
	  for(int j = 9; j<15; j++){
	    bbgemgentrack->SetData(j, 0.);
	  }
	}
	bbgemgentrack->SetData(15, 1.0);// TODO: weight
	fg4sbsGenData.push_back(bbgemgentrack);

	
	if(fTree->Earm_GRINCH.nhits){
	  z_det = fzGrinch-0.80;
	  g4sbsgendata *grinchgenhit = new g4sbsgendata(GRINCH_UNIQUE_DETID, 8);
	  grinchgenhit->SetData(0, fSource); 
	  grinchgenhit->SetData(1, t.TID->at(k)); 
	  grinchgenhit->SetData(2, t.PID->at(k)); 
	  grinchgenhit->SetData(3, t.X->at(k)+t.Xp->at(k)*z_det); 
	  grinchgenhit->SetData(4, t.Y->at(k)+t.Yp->at(k)*z_det); 
	  grinchgenhit->SetData(5, t.T->at(k)+z_det*t.P->at(k)/pz/(beta*0.299792458) );
	  grinchgenhit->SetData(6, E);
	  grinchgenhit->SetData(7, 1.0);// TODO: weight
	  fg4sbsGenData.push_back(grinchgenhit);
	}
	
	if(fTree->Earm_BBHodoScint.nhits){  
	  z_det = fzHodo-0.80;
	  g4sbsgendata *hodogenhit = new g4sbsgendata(HODO_UNIQUE_DETID, 8);
	  hodogenhit->SetData(0, fSource); 
	  hodogenhit->SetData(1, t.TID->at(k)); 
	  hodogenhit->SetData(2, t.PID->at(k)); 
	  hodogenhit->SetData(3, t.X->at(k)+t.Xp->at(k)*z_det); 
	  hodogenhit->SetData(4, t.Y->at(k)+t.Yp->at(k)*z_det); 
	  hodogenhit->SetData(5, t.T->at(k)+z_det*t.P->at(k)/pz/(beta*0.299792458) );
	  hodogenhit->SetData(6, E);	  
	  hodogenhit->SetData(7, 1.0);// TODO: weight
	  fg4sbsGenData.push_back(hodogenhit);
	}
	
	if(fTree->Earm_BBPSTF1.nhits){
	  z_det = fzPS-0.80;
	  g4sbsgendata *bbpsgenhit = new g4sbsgendata(BBPS_UNIQUE_DETID, 8);
	  bbpsgenhit->SetData(0, fSource); 
	  bbpsgenhit->SetData(1, t.TID->at(k)); 
	  bbpsgenhit->SetData(2, t.PID->at(k)); 
	  bbpsgenhit->SetData(3, t.X->at(k)+t.Xp->at(k)*z_det); 
	  bbpsgenhit->SetData(4, t.Y->at(k)+t.Yp->at(k)*z_det); 
	  bbpsgenhit->SetData(5, t.T->at(k)+z_det*t.P->at(k)/pz/(beta*0.299792458) );
	  bbpsgenhit->SetData(6, E);	  
	  bbpsgenhit->SetData(7, 1.0);// TODO: weight
	  fg4sbsGenData.push_back(bbpsgenhit);
	}
	
	if(fTree->Earm_BBSHTF1.nhits){
	  double z_det = fzSH-0.80;
	  g4sbsgendata *bbshgenhit = new g4sbsgendata(BBSH_UNIQUE_DETID, 8);
	  bbshgenhit->SetData(0, fSource); 
	  bbshgenhit->SetData(1, t.TID->at(k)); 
	  bbshgenhit->SetData(2, t.PID->at(k)); 
	  bbshgenhit->SetData(3, t.X->at(k)+t.Xp->at(k)*z_det); 
	  bbshgenhit->SetData(4, t.Y->at(k)+t.Yp->at(k)*z_det); 
	  bbshgenhit->SetData(5, t.T->at(k)+z_det*t.P->at(k)/pz/(beta*0.299792458) );
	  bbshgenhit->SetData(6, E);	  
	  bbshgenhit->SetData(7, 1.0);// TODO: weight
	  fg4sbsGenData.push_back(bbshgenhit);
	}
      }// end loop on tracks
    }
  }//end loop on GEM tree; 
  
  if(fTree->Earm_GRINCH.nhits){
    if(d_flag>=3)printf("Nhits in GRINCH = %d\n", fTree->Earm_GRINCH.nhits);
    for(int i = 0; i<fTree->Earm_GRINCH.nhits; i++){
      g4sbshitdata *grinchhit = new g4sbshitdata(GRINCH_UNIQUE_DETID, 4);
      grinchhit->SetData(0, fSource);
      grinchhit->SetData(1, int(fTree->Earm_GRINCH.PMT->at(i)/5)-1);
      // -1 is to range from 0 to 509 instead of 1 to 510 :/
      grinchhit->SetData(2, fTree->Earm_GRINCH.Time_avg->at(i));
      grinchhit->SetData(3, fTree->Earm_GRINCH.NumPhotoelectrons->at(i));
      fg4sbsHitData.push_back(grinchhit);
    }
    if(d_flag>=3)printf("Accumulated data = %lu\n", fg4sbsHitData.size());
  }
    
  // Process hodoscope data
  if(fTree->Earm_BBHodoScint.nhits){
    if(d_flag>=3)printf("Nhits in BBhodo = %d\n", fTree->Earm_BBHodoScint.nhits);
    for(int i = 0; i<fTree->Earm_BBHodoScint.nhits; i++){
      for(int j = 0; j<2; j++){//j = 0: close beam PMT, j = 1: far beam PMT
	// Evaluation of number of photoelectrons and time from energy deposit documented at:
	// https://hallaweb.jlab.org/dvcslog/SBS/170711_172759/BB_hodoscope_restudy_update_20170711.pdf
	// TODO: put that stuff in DB...
	Npe = fRN->Poisson(1.0e7*fTree->Earm_BBHodoScint.sumedep->at(i)*0.113187*exp(-(0.3+pow(-1, j)*fTree->Earm_BBHodoScint.xhit->at(i))/1.03533)* 0.24);
	t = fTree->Earm_BBHodoScint.tavg->at(i)+(0.55+pow(-1, j)*fTree->Earm_BBHodoScint.xhit->at(i))/0.15;
	g4sbshitdata *hodopmthit = new g4sbshitdata(HODO_UNIQUE_DETID, 5);
	hodopmthit->SetData(0, fSource);
	hodopmthit->SetData(1, fTree->Earm_BBHodoScint.cell->at(i)*2+j);
	hodopmthit->SetData(2, t);
	hodopmthit->SetData(3, Npe);
	hodopmthit->SetData(4, fTree->Earm_BBHodoScint.sumedep->at(i));
	fg4sbsHitData.push_back(hodopmthit);
      }
    }
    if(d_flag>=3)printf("Accumulated data = %lu\n", fg4sbsHitData.size());
  }

  // Process BB PS data
  if(fTree->Earm_BBPSTF1.nhits){
    if(d_flag>=3)printf("Nhits in BBPSTF1 = %d\n", fTree->Earm_BBPSTF1.nhits);
    for(int i = 0; i<fTree->Earm_BBPSTF1.nhits; i++){
      // Evaluation of number of photoelectrons and time from energy deposit documented at:
      // 
      // TODO: put that stuff in DB...
      //Hard cutoff at 0.3MeV
      if(fTree->Earm_BBPSTF1.sumedep->at(i)<3.e-4)continue;
      //check probability to generate p.e. yield
      bool genpeyield = true;
      if(fTree->Earm_BBPSTF1.sumedep->at(i)<1.e-2)
	genpeyield = fRN->Uniform(0, 1)<=1-exp(0.29-950.*fTree->Earm_BBPSTF1.sumedep->at(i));
      //if we're go, let's generate
      if(genpeyield){
	//Used to be 454.: just wrong
	Npe = fRN->Poisson(1500.0*fTree->Earm_BBPSTF1.sumedep->at(i));
	t = fTree->Earm_BBPSTF1.tavg->at(i)+fRN->Gaus(3.2-5.805*fTree->Earm_BBPSTF1.zhit->at(i)-17.77*pow(fTree->Earm_BBPSTF1.zhit->at(i), 2), 0.5);
	g4sbshitdata *bbpshit = new g4sbshitdata(BBPS_UNIQUE_DETID, 5);
	bbpshit->SetData(0, fSource);
	bbpshit->SetData(1, fTree->Earm_BBPSTF1.cell->at(i));
	bbpshit->SetData(2, t);
	bbpshit->SetData(3, Npe);
	bbpshit->SetData(4, fTree->Earm_BBPSTF1.sumedep->at(i));
	fg4sbsHitData.push_back(bbpshit);
      }
    }
    if(d_flag>=3)printf("Accumulated data = %lu\n", fg4sbsHitData.size());
  }
  
  // Process BB SH data
  if(fTree->Earm_BBSHTF1.nhits){
    if(d_flag>=3)printf("Nhits in BBSHTF1 = %d\n", fTree->Earm_BBSHTF1.nhits);
    for(int i = 0; i<fTree->Earm_BBSHTF1.nhits; i++){
      // Evaluation of number of photoelectrons and time from energy deposit documented at:
      // 
      // TODO: put that stuff in DB...
      //Used to be 932.: just wrong
       if(fTree->Earm_BBSHTF1.sumedep->at(i)<3.e-4)continue;
      //check probability to generate p.e. yield
      bool genpeyield = true;
      if(fTree->Earm_BBSHTF1.sumedep->at(i)<1.e-2)genpeyield = fRN->Uniform(0, 1)<=1-exp(0.29-950.*fTree->Earm_BBSHTF1.sumedep->at(i));
      //if we're go, let's generate
      if(genpeyield){
	Npe = fRN->Poisson(1800.0*fTree->Earm_BBSHTF1.sumedep->at(i));
	t = fTree->Earm_BBSHTF1.tavg->at(i)+fRN->Gaus(2.216-8.601*fTree->Earm_BBSHTF1.zhit->at(i)-7.469*pow(fTree->Earm_BBSHTF1.zhit->at(i), 2), 0.8);
	g4sbshitdata *bbshhit = new g4sbshitdata(BBSH_UNIQUE_DETID, 5);
	bbshhit->SetData(0, fSource);
	bbshhit->SetData(1, fTree->Earm_BBSHTF1.cell->at(i));
	//bbshhit->SetData(2, 0);
	bbshhit->SetData(2, t);
	bbshhit->SetData(3, Npe);
	bbshhit->SetData(4, fTree->Earm_BBSHTF1.sumedep->at(i));
	fg4sbsHitData.push_back(bbshhit);
      }
     }
    if(d_flag>=3)printf("Accumulated data = %lu\n" , fg4sbsHitData.size());
  }
  
  // Hadron Arm
  // Process CDet data
  if(fTree->Harm_CDET_Scint.nhits){
    if(d_flag>=3)printf("Nhits in CDet = %d\n", fTree->Harm_CDET_Scint.nhits);
    for(int i = 0; i<fTree->Harm_CDET_Scint.nhits; i++){
      // Evaluation of number of photoelectrons and time from energy deposit documented at:
      // 
      // TODO: put that stuff in DB...
      Npe = fRN->Poisson( fTree->Harm_CDET_Scint.sumedep->at(i)*5.634e3 );
      t = fTree->Harm_CDET_Scint.tavg->at(i)+5.75+fTree->Harm_CDET_Scint.xhit->at(i)/0.16;
      g4sbshitdata *cdetpmthit = new g4sbshitdata(CDET_UNIQUE_DETID, 5);
      cdetpmthit->SetData(0, fSource);
      cdetpmthit->SetData(1, fTree->Harm_CDET_Scint.cell->at(i));
      //cdetpmthit->SetData(2, 0);
      cdetpmthit->SetData(2, t);
      cdetpmthit->SetData(3, Npe);
      cdetpmthit->SetData(4, fTree->Harm_CDET_Scint.sumedep->at(i));
      fg4sbsHitData.push_back(cdetpmthit);

    }
    if(d_flag>=3)printf("Accumulated data = %lu\n" , fg4sbsHitData.size());
  }

  // GEp Electron Arm
  // Process CDet data
  if(fTree->Earm_CDET_Scint.nhits){
    if(d_flag>=3)printf("Nhits in CDet = %d\n", fTree->Earm_CDET_Scint.nhits);
    for(int i = 0; i<fTree->Earm_CDET_Scint.nhits; i++){
      // Evaluation of number of photoelectrons and time from energy deposit documented at:
      // 
      // TODO: put that stuff in DB...
      Npe = fRN->Poisson( fTree->Earm_CDET_Scint.sumedep->at(i)*5.634e3 );
      t = fTree->Harm_CDET_Scint.tavg->at(i)+5.75+fTree->Earm_CDET_Scint.xhit->at(i)/0.16;
      g4sbshitdata *cdetpmthit = new g4sbshitdata(CDET_UNIQUE_DETID, 5);
      cdetpmthit->SetData(0, fSource);
      cdetpmthit->SetData(1, fTree->Earm_CDET_Scint.cell->at(i));
      //cdetpmthit->SetData(2, 0);
      cdetpmthit->SetData(2, t);
      cdetpmthit->SetData(3, Npe);
      cdetpmthit->SetData(4, fTree->Earm_CDET_Scint.sumedep->at(i));
      fg4sbsHitData.push_back(cdetpmthit);
    }
    if(d_flag>=3)printf("Accumulated data = %lu\n" , fg4sbsHitData.size());
  }
  
  // GEN RP detectors: Active analyzer 
  if(fTree->Harm_ActAnScint.nhits){
    if(d_flag>=3)printf("Nhits in GEn RP active analyzer = %d\n", fTree->Harm_ActAnScint.nhits);
    for(int i = 0; i<fTree->Harm_ActAnScint.nhits; i++){
      //TODO: update with accurate paramters
      Npe = fRN->Poisson( fTree->Harm_ActAnScint.sumedep->at(i)*5.634e3 ); //
      t = fTree->Harm_ActAnScint.tavg->at(i)+5.75+fTree->Harm_ActAnScint.xhit->at(i)/0.16; //
      g4sbshitdata *actanapmthit = new g4sbshitdata(ACTIVEANA_UNIQUE_DETID, 5);
      actanapmthit->SetData(0, fSource);
      actanapmthit->SetData(1, fTree->Harm_ActAnScint.cell->at(i));
      actanapmthit->SetData(2, t);
      actanapmthit->SetData(3, Npe);
      actanapmthit->SetData(4, fTree->Harm_ActAnScint.sumedep->at(i));
      fg4sbsHitData.push_back(actanapmthit);
    }
    if(d_flag>=3)printf("Accumulated data = %lu\n" , fg4sbsHitData.size());
  }
  
  
  
  // Now process the GEM data
  if(d_flag>=3)printf("about to digitize GEMs: %lu tree(s)\n", fTree->GEMs.size());
  for(std::vector<gem_branch>::iterator it = fTree->GEMs.begin();
       it != fTree->GEMs.end(); it++) {
    if(d_flag>=3)printf("GEM tree %ld\n", std::distance(fTree->GEMs.begin(), it));
    
    TSBSGeant4::GEMData_t &t = it->tree;
    if(t.plane) {
      if(d_flag>=3)printf("%d hits \n", t.nhits);
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
	  
	  /*
	  if(d_flag>=3){
	    TSBSGeant4::TrackerData_t &tt = it->Track_tree;
	    if(tt.ntracks) {
	      if(d_flag>=3)printf("%d tracks \n", tt.ntracks);
	      for(int kk = 0; kk < tt.ntracks; kk++) {
		if(t.trid->at(k)==tt.TID->at(kk)){
		  printf("plane %d, z = %f\n", t.plane->at(k), t.z->at(k)); 
		  double dx = t.tx->at(k)-tt.X->at(kk);
		  double dy = t.ty->at(k)-tt.Y->at(kk);
		  printf("x: %f-%f = %f\n", t.tx->at(k), tt.X->at(kk), dx);
		  printf("y: %f-%f = %f\n", t.ty->at(k), tt.Y->at(kk), dy);
		  double dz1 = dx/t.txp->at(k);
		  double dz2 = dy/t.typ->at(k);
		  double dz3 = dx/tt.Xp->at(kk);
		  double dz4 = dy/tt.Yp->at(kk);
		  printf("dx/dz: %f (hit) => dz = %f ; %f (track) => dz = %f \n",
			 t.txp->at(k), dz1, tt.Xp->at(kk), dz3);
		  printf("dy/dz: %f (hit) => dz = %f ; %f (track) => dz = %f \n",
			 t.typ->at(k), dz2, tt.Yp->at(kk), dz4);
		  printf("z0 = %f, %f, %f, %f\n", 
			 t.z->at(k)-dz1, t.z->at(k)-dz2, 
			 t.z->at(k)-dz3, t.z->at(k)-dz4);
		}
	      }
	    }//
	  }//end if d_flag
	  */
        }
      }
    }
    if(d_flag>=3)printf("Accumulated data = %lu\n" , fg4sbsHitData.size());
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
  if(d_flag>=3)printf("about to digitize HCal\n");
  // Now process HCAL data
  /*
  if(fTree->hcalpart.E) {
    if(d_flag>=3)printf("Nhits in HCal = %lu\n", fTree->hcalpart.E->size());
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
    if(d_flag>=3)printf("Accumulated data = %lu\n" , fg4sbsHitData.size());
  }
  */
  // For now, get the adc signal from the total energy deposited on the
  // scintillators. This can be changed later, I suppose...
  if(fTree->hcalscint.sumedep) {
    if(d_flag>=3)printf("Nhits in HCalScint = %lu\n", fTree->hcalscint.sumedep->size());
    for(size_t k = 0; k < fTree->hcalscint.sumedep->size(); k++) {
      g4sbshitdata *hcalscinthit = new g4sbshitdata(HCAL_UNIQUE_DETID,4);
      hcalscinthit->SetData(0,fSource);
      hcalscinthit->SetData(1,fTree->hcalscint.cell->at(k));
      hcalscinthit->SetData(2,1); // data type 1 == sumedep
      hcalscinthit->SetData(3,fTree->hcalscint.sumedep->at(k));
      fg4sbsHitData.push_back(hcalscinthit);

      z_hit = -(fTree->hcalscint.xhitg->at(k)-x_ref)*sin(fThetaSBS)+(fTree->hcalscint.zhitg->at(k)-z_ref)*cos(fThetaSBS);
      
      // Evaluation of number of photoelectrons from energy deposit documented at:
      // https://sbs.jlab.org/DocDB/0000/000043/002/HCal_Digi_EdepOnly_2.pdf
      // TODO: put that stuff in DB...
      Npe_Edep_ratio = 5.242+11.39*z_hit+10.41*pow(z_hit, 2);
      Npe = fRN->Poisson(Npe_Edep_ratio*fTree->hcalscint.sumedep->at(k)*1.0e3);
      t = fRN->Gaus(fTree->hcalscint.tavg->at(k)+10.11, 1.912);
      
      sigma_tgen = 0.4244+11380/pow(Npe+153.4, 2);
      for(int l = 0; l<Npe; l++){
	//Generate here,...
	g4sbshitdata *hcalpmthit = new g4sbshitdata(HCAL_UNIQUE_DETID,4);
	hcalpmthit->SetData(0,fSource);
	hcalpmthit->SetData(1,fTree->hcalscint.cell->at(k));
	hcalpmthit->SetData(2,0); // data type 0 == optical photon
	hcalpmthit->SetData(3,fRN->Landau(t, sigma_tgen));
	fg4sbsHitData.push_back(hcalpmthit);
      }
    }
    if(d_flag>=3)printf("Accumulated data = %lu\n", fg4sbsHitData.size());
  }
  
  if(d_flag>=5)for(uint i = 0; i<fg4sbsHitData.size(); i++){
      printf("Det ID %d : chan %1.0f \n", fg4sbsHitData.at(i)->GetDetUniqueID(), fg4sbsHitData.at(i)->GetData(1));
    }
  
  if(d_flag>=3){
    printf("Just read event %d\n", fEvNum);
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

double TSBSGeant4File::M_p(int pid)
{
  switch(TMath::Abs(pid)){
  case (11):
    return 0.000511;
    break;
  case (13):
    return 0.105658;
    break;
  case (111):
    return 0.134977;
    break;
  case (211):
    return 0.139570;
    break;
  case (321):
    return 0.493677;
    break;
  case (2212):
    return 0.938272;
    break;
  default:
    return 0.0;
    break;
  }
  
}


#endif//__CINT__
