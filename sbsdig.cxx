//includes: standard
#include <iostream>
#include <fstream>
#include <string>
#include <map>

//includes: root
#include <TROOT.h>
#include "TString.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TCut.h"
#include "TEventList.h"
#include "TMath.h"
#include "TRandom3.h"

//includes: specific
#include "G4SBSRunData.hh"
#include "g4sbs_types.h"
#include "gmn_tree.h"
#include "SBSDigAuxi.h"
#include "SBSDigPMTDet.h"
#include "SBSDigGEMDet.h"
#include "SBSDigGEMSimDig.h"
#include "SBSDigBkgdGen.h"

#ifdef __APPLE__
#include "unistd.h"
#endif

#define qe 1.602e-19
#define spe_unit 1.0e-9

#define NPlanes_BBGEM 32 // modules...
#define NChan_BBPS 52
#define NChan_BBSH 189
#define NChan_BBHODO 180 
#define NChan_GRINCH 510 
#define NChan_HCAL 288

#define ROimpedance 50.0 //Ohm
#define TriggerJitter 3.0 //ns
#define ADCbits 12 
#define gatewidth_PMT 100 //ns
#define gatewidth_GEM 400 //ns

#define FADC_sampsize 4.0 //ns

//DB???
#define sigmapulse_BBPSSH 3.0 //ns / 1.2 of 
#define sigmapulse_BBHODO 1.6 //ns
#define sigmapulse_GRINCH 3.75 //ns

#define gain_BBPS 2.e6
#define ped_BBPS 600.0 // ADC channel
#define pedsigma_BBPS 3.0 // ADC channel
#define trigoffset_BBPS 18.2 //ns
#define ADCconv_BBPS 50 //fC/ch

#define gain_BBSH 7.5e5
#define ped_BBSH 500.0 // ADC channel
#define pedsigma_BBSH 4.5 // ADC channel
#define trigoffset_BBSH 18.5 //ns
#define ADCconv_BBSH 50 //fC/ch

#define gain_GRINCH 7.0e6
#define ped_GRINCH 0.0 
#define pedsigma_GRINCH 0.0
#define trigoffset_GRINCH 15.3 //ns
#define threshold_GRINCH 3.e-3 //V
#define ADCconv_GRINCH 100 //fC/ch
#define TDCconv_GRINCH 1.0 //ns/channel
#define TDCbits_GRINCH 16 //ns/channel

#define gain_BBHODO 1.0e5
#define ped_BBHODO 0.0 
#define pedsigma_BBHODO 0.0
#define trigoffset_BBHODO 18.6 //ns
#define threshold_BBHODO 3.e-3 //V
#define ADCconv_BBHODO 100 //fC/ch
#define TDCconv_BBHODO 0.1 //ns/channel
#define TDCbits_BBHODO 19 //ns/channel

#define gain_HCAL 1.0e6
#define ped_HCAL 0.0 
#define pedsigma_HCAL 0.0
#define trigoffset_HCAL 81.0 //ns
#define threshold_HCAL 3.e-3 //V
#define ADCconv_HCAL 1.0 //fC/ch //??
#define TDCconv_HCAL 0.12 //ns/channel
#define TDCbits_HCAL 16 //ns/channel




using namespace std;
//____________________________________________________
int main(int argc, char** argv){
  
  // Step 0: read out arguments
  string inputsigfile, inputbkgdfile = "";//sources of files
  ULong64_t Nentries = -1;//number of events to process
  //UShort_t Nbkgd = 0;//number of background files to add to each event
  double LumiFrac = 0;
  
  if(argc<2){
    cout << "*** Not enough arguments! ***" << endl
	 << " Arguments: list_of_sig_input_files (str, mandatory); " << endl
	 << "           nb_of_sig_evts_to_process (int, def=-1); " << endl
	 << "          list_of_bkgd_input_files (str, def=''); " << endl
	 << "         nb_of_bkgd_files_to_add_to_sig_evt (int, def=0); " << endl;
    return(-1);
  }
  
  inputsigfile = argv[1];
  cout << " Signal input files from: " << inputsigfile << endl;
  if(argc>2)Nentries = atoi(argv[2]);
  cout << " Number of (signal) events to process = " << Nentries << endl;
  if(argc>4){
    inputbkgdfile = argv[3];
    cout << " Background hsotgrams from: " << inputbkgdfile << endl;
    LumiFrac = max(0., atof(argv[4]));
    cout << " Fraction of background to superimpose to signal = " << LumiFrac << endl;
  }
  
  TFile* f_bkgd;
  SBSDigBkgdGen* BkgdGenerator;
  if(LumiFrac>0){
    f_bkgd = TFile::Open(inputbkgdfile.c_str());
    if(f_bkgd->IsZombie()){
      LumiFrac = 0;
    }else{
      BkgdGenerator = new SBSDigBkgdGen(f_bkgd);
    }
  }
  
  // ------------------- // dev notes // ------------------- //
  // The loop on the input signal and background chains 
  // is going to happen here in the main I guess.
  //
  // First, we want to extend the input tree (for signal only!!!)
  // I guess in order to avoid adding extra layers of code, 
  // the tree extension might have to be coded in the custom tree class
  
  // Step 1: read input files build the input chains
  TRandom3* R = new TRandom3(0);
  TString currentline;
  
  //double nstrips_bbgem[NPlanes_BBGEM] = {3840, 3072, 3840, 3072, 3840, 3072, 3840, 3072, 3840, 6144};
  double nstrips_bbgem[NPlanes_BBGEM] = {1280, 1024, 1280, 1024, 1280, 1024, 
					 1280, 1024, 1280, 1024, 1280, 1024, 
					 1280, 1024, 1280, 1024, 1280, 1024, 
					 1280, 1024, 1280, 1024, 1280, 1024, 
					 1280, 1536, 1280, 1536, 
					 1280, 1536, 1280, 1536};
  double offset_bbgem[NPlanes_BBGEM] = {-0.512, 0., 0., 0., 0.512, 0., 
					-0.512, 0., 0., 0., 0.512, 0., 
					-0.512, 0., 0., 0., 0.512, 0., 
					-0.512, 0., 0., 0., 0.512, 0., 
					-0.768, 0., -0.256, 0., 
					 0.256, 0.,  0.768, 0.};
  double angle_bbgem[NPlanes_BBGEM] = {0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 
				       0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 
				       0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 
				       0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 
				       0.0, 90.0, 0.0, 90.0, 
				       0.0, 90.0, 0.0, 90.0};
  for(int i = 0; i<NPlanes_BBGEM; i++){
    angle_bbgem[i]*= TMath::DegToRad();
    //cout << nstrips_bbgem[i] << " ";
  }//cout << endl;
  
  double triggeroffset[NPlanes_BBGEM/2] = {121., 121., 121., 121.5, 121.5, 121.5,  
					   122., 122., 122., 122.5, 122.5, 122.5,  
					   126., 126., 126., 126.};
  
  double commonmode_array[1] = {1500.};
    
  //double NpeChargeConv
  
  //Declaring detectors
  SBSDigGEMDet* bbgem = new SBSDigGEMDet(BBGEM_UNIQUE_DETID, NPlanes_BBGEM, nstrips_bbgem, offset_bbgem, angle_bbgem, 6, 240.);
  SBSDigGEMSimDig* gemdig = new SBSDigGEMSimDig(NPlanes_BBGEM, triggeroffset, 240., 1, commonmode_array);
  
  SBSDigPMTDet* bbps = new SBSDigPMTDet(BBPS_UNIQUE_DETID, NChan_BBPS, gain_BBPS*qe, sigmapulse_BBPSSH, gatewidth_PMT);
  SBSDigPMTDet* bbsh = new SBSDigPMTDet(BBSH_UNIQUE_DETID, NChan_BBSH, gain_BBSH*qe, sigmapulse_BBPSSH, gatewidth_PMT);
  SBSDigPMTDet* grinch = new SBSDigPMTDet(GRINCH_UNIQUE_DETID, NChan_GRINCH, gain_GRINCH*qe, sigmapulse_GRINCH, gatewidth_PMT);
  SBSDigPMTDet* bbhodo = new SBSDigPMTDet(HODO_UNIQUE_DETID, NChan_BBHODO, gain_BBHODO*qe, sigmapulse_BBHODO, gatewidth_PMT);
  SBSDigPMTDet* hcal = new SBSDigPMTDet(HCAL_UNIQUE_DETID, NChan_HCAL);
  
  bbps->fGain = gain_BBPS;
  bbps->fPedestal = ped_BBPS;
  bbps->fPedSigma = pedsigma_BBPS;
  bbps->fTrigOffset = trigoffset_BBPS;
  bbps->fGateWidth = gatewidth_PMT;
  bbps->fADCconv = ADCconv_BBPS;
  bbps->fADCbits = ADCbits;
  
  bbsh->fGain = gain_BBSH;
  bbsh->fPedestal = ped_BBSH;
  bbsh->fPedSigma = pedsigma_BBSH;
  bbsh->fTrigOffset = trigoffset_BBSH;
  bbsh->fGateWidth = gatewidth_PMT;
  bbsh->fADCconv = ADCconv_BBSH;
  bbsh->fADCbits = ADCbits;
  
  grinch->fGain = gain_GRINCH;
  grinch->fPedestal = ped_GRINCH;
  grinch->fPedSigma = pedsigma_GRINCH;
  grinch->fTrigOffset = trigoffset_GRINCH;
  grinch->fThreshold = threshold_GRINCH*spe_unit/ROimpedance;
  grinch->fGateWidth = gatewidth_PMT;
  grinch->fADCconv = ADCconv_GRINCH;
  grinch->fADCbits = ADCbits;
  grinch->fTDCconv = TDCconv_GRINCH;
  grinch->fTDCbits = TDCbits_GRINCH; 
  
  bbhodo->fGain = gain_BBHODO;
  bbhodo->fPedestal = ped_BBHODO;
  bbhodo->fPedSigma = pedsigma_BBHODO;
  bbhodo->fTrigOffset = trigoffset_BBHODO;
  bbhodo->fThreshold = threshold_BBHODO*spe_unit/ROimpedance;
  bbhodo->fGateWidth = gatewidth_PMT;
  bbhodo->fADCconv = ADCconv_BBHODO;
  bbhodo->fADCbits = ADCbits;
  bbhodo->fTDCconv = TDCconv_BBHODO;
  bbhodo->fTDCbits = TDCbits_BBHODO; 
  
  hcal->fGain = gain_HCAL;
  hcal->fPedestal = ped_HCAL;
  hcal->fPedSigma = pedsigma_HCAL;
  hcal->fTrigOffset = trigoffset_HCAL;
  hcal->fThreshold = threshold_HCAL*spe_unit/ROimpedance;
  hcal->fGateWidth = gatewidth_PMT-20.;
  hcal->fADCconv = ADCconv_HCAL;
  hcal->fADCbits = ADCbits;
  hcal->fTDCconv = TDCconv_HCAL;
  hcal->fTDCbits = TDCbits_HCAL; 
  hcal->SetSamples(FADC_sampsize);
  
  std::vector<SBSDigPMTDet*> PMTdetectors;
  std::vector<int> detmap;
  std::vector<SBSDigGEMDet*> GEMdetectors;
  std::vector<int> gemmap;
  
  //ordered by increasing uinque id
  PMTdetectors.push_back(hcal);  detmap.push_back(HCAL_UNIQUE_DETID);
  PMTdetectors.push_back(bbps);  detmap.push_back(BBPS_UNIQUE_DETID);
  PMTdetectors.push_back(bbsh);  detmap.push_back(BBSH_UNIQUE_DETID);
  PMTdetectors.push_back(grinch);  detmap.push_back(GRINCH_UNIQUE_DETID);
  PMTdetectors.push_back(bbhodo);  detmap.push_back(HODO_UNIQUE_DETID);

  GEMdetectors.push_back(bbgem); gemmap.push_back(BBGEM_UNIQUE_DETID);
  /*  
  std::map<int, SBSDigPMTDet*> PMTdetectors;
  PMTdetectors[HCAL_UNIQUE_DETID] = hcal;
  PMTdetectors[HODO_UNIQUE_DETID] = bbhodo;
  PMTdetectors[BBPS_UNIQUE_DETID] = bbps;
  PMTdetectors[BBSH_UNIQUE_DETID] = bbsh;
  PMTdetectors[GRINCH_UNIQUE_DETID] = grinch;
  std::map<int, SBSDigGEMDet*> GEMdetectors;
  GEMdetectors[BBGEM_UNIQUE_DETID] = bbgem;
  */
  
  // build signal chain
  ifstream sig_inputfile(inputsigfile);
  TChain *C_s = new TChain("T");
  while( currentline.ReadLine(sig_inputfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C_s->Add(currentline.Data());
    }
  }
  TObjArray *fileElements_s=C_s->GetListOfFiles();
  TIter next_s(fileElements_s);
  TChainElement *chEl_s=0;
  
  /* need to change this... 
  // build background chain
  ifstream beam_inputfile(inputbkgdfile);
  TChain *C_b = new TChain("T");
  //TCut global_cut = "";
  //TEventList *eblist = new TEventList("eblist");
  if(Nbkgd!=0){
    while( currentline.ReadLine(beam_inputfile) && !currentline.BeginsWith("endlist") ){
      if( !currentline.BeginsWith("#") ){
	C_b->Add(currentline.Data());
	//global_cut += currentline.Data();
	//cout << currentline.Data() << endl;
      }
    }
    //C_b->Draw(">>eblist",global_cut);
  }
  //TObjArray *fileElements_b=C_b->GetListOfFiles();
  //TIter next_b(fileElements_b);
  //TChainElement *chEl_b=0;
  */
  
  G4SBSRunData* run_data;
  
  double Theta_SBS, D_HCal;
  
  gmn_tree *T_s;//, *T_b;

  ULong64_t Nev_fs;//, Nev_fb;
  ULong64_t ev_s;//, ev_b;
  
  ULong64_t NEventsTotal = 0;
  //UShort_t nbkgd = 0;
  //int treenum = 0;
  //int oldtreenum = 0;
  
  int i_fs = 0;
  bool has_data;
  
  double timeZero;
  
  //T_b = new gmn_tree(C_b);
  //ev_b = 0;
  
  while (( chEl_s=(TChainElement*)next_s() )) {
    if(NEventsTotal>=Nentries){
      break;
    }
    TFile f_s(chEl_s->GetTitle(), "UPDATE");
    run_data = (G4SBSRunData*)f_s.Get("run_data");
    Theta_SBS = run_data->fSBStheta;
    D_HCal = run_data->fHCALdist;
    //TFile fs_c(Form("digitized/simdigtest_%d.root", i_fs), "UPDATE");
    //f_s.Cp(Form("digitized/simdigtest_%d.root", i_fs));
    //if(fs_c.IsOpen())cout << "copy of file is open" << endl;
    //cout << fs_c->ReOpen("UPDATE") << endl;
    //C_s = (TChain*)fs_c.Get("T");
    C_s = (TChain*)f_s.Get("T");
    T_s = new gmn_tree(C_s);
    
    // Expend tree here! (again, for signal only!!!)
    T_s->AddDigBranches();
    
    Nev_fs = C_s->GetEntries();
    
    for(ev_s = 0; ev_s<Nev_fs; ev_s++, NEventsTotal++){
      if(NEventsTotal>=Nentries)break;
      if(NEventsTotal%1000==0)
	cout << NEventsTotal << "/" << Nentries << endl;
      
      timeZero = R->Gaus(0.0, TriggerJitter);
      
      bbgem->Clear();
      //for(int i = 0; i<NPlanes_BBGEM; i++){
      //cout << bbgem->GEMPlanes[i].GetNStrips() << " ";
      //}cout << endl;
      bbps->Clear();
      bbsh->Clear();
      grinch->Clear();
      bbhodo->Clear();
      hcal->Clear(true);
      
      has_data = false;
      
      T_s->ClearDigBranches();
      T_s->GetEntry(ev_s);
      
      // unfold the thing then... but where???
      has_data = UnfoldData(T_s, Theta_SBS, D_HCal, R, PMTdetectors, detmap, GEMdetectors, gemmap, timeZero, 0);
      if(!has_data)continue;
      
      
      if(LumiFrac>0){
	BkgdGenerator->GenerateBkgd(R, PMTdetectors, detmap, GEMdetectors, gemmap, LumiFrac);
      }
      /*
      // loop here for background
      if(Nbkgd>0){
	nbkgd = 0;
	while( T_b->GetEntry(ev_b++) ){
	//while( T_b->GetEntry(eblist->GetEntry(ev_b++)) ){
	  treenum = C_b->GetTreeNumber();
	  if(treenum!=oldtreenum){
	    oldtreenum = treenum;
	    nbkgd++;
	    if(nbkgd>=Nbkgd)break;
	  }
	  timeZero = R->Uniform( -gatewidth_GEM-50., gatewidth_GEM/2.-50. );
	  
	  UnfoldData(T_b, Theta_SBS, D_HCal, R, PMTdetectors, detmap, GEMdetectors, gemmap, timeZero, 1);
	  //if(treenum)
	}
	
	while (( chEl_b=(TChainElement*)next_b() )) {
	  if(nbkgd>=Nbkgd)break;
	  cout << chEl_b->GetTitle() << endl;
	  TFile f_b(chEl_b->GetTitle());
	  C_b = (TChain*)f_b.Get("T");
	  T_b = new gmn_tree(C_b);
	  Nev_fb = C_b->GetEntries();
	  for(ev_b = 0; ev_b<Nev_fb; ev_b++){
	    T_b->GetEntry(ev_b);
	    UnfoldData(T_b, Theta_SBS, D_HCal, R);
	  }// end loop on background events
	  nbkgd++;
	}// end loop on background files
      }//end if Nbkgd>0
      */
      
      bbps->Digitize(T_s,R);
      bbsh->Digitize(T_s,R);
      grinch->Digitize(T_s,R);
      bbhodo->Digitize(T_s,R);
      hcal->Digitize(T_s,R);
      gemdig->Digitize(bbgem, R);
      //cout << " hou hou " << bbgem->GEMPlanes[4].GetADCSum(400) << endl;
      // for(int j = 570; j<580; j++){
      // 	cout << " * " << j << "   ";
      // 	for(int b = 0; b<6; b++){
      // 	  cout << " " << bbgem->GEMPlanes[30].GetADC(j, b);
      // 	}cout << endl;
      // }
      gemdig->CheckOut(bbgem, R, T_s);
      
      //How come this function is so taxing in time??? 
      // Answer in the function... hope we've found a workaround
      //FillDigTree(T_s, PMTdetectors, GEMdetectors);

      T_s->FillDigBranches();
      //T_s->fChain->Fill();
    }// end loop on signal events 
    
    T_s->fChain->Write("", TObject::kOverwrite);
    //fs_c.Write();
    //fs_c.Close();
    f_s.Write();
    f_s.Close();
    i_fs++;
  }// end loop on signal files
  
  
  exit(0);
}
