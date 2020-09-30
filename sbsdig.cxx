//includes: standard
#include <iostream>
#include <fstream>
#include <string>
#include <map>

//includes: root
#include <TROOT.h>
#include "TString.h"
#include "TObjString.h"
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
#include "g4sbs_tree.h"
#include "SBSDigAuxi.h"
#include "SBSDigPMTDet.h"
#include "SBSDigGEMDet.h"
#include "SBSDigGEMSimDig.h"
#include "SBSDigBkgdGen.h"

#ifdef __APPLE__
#include "unistd.h"
#endif

/*
//Defining here the parameters for the new detectors.
//TODO: write a list of parameters that are not "frozen" (e.g. gain, pedestal parameters, etc...) and switch them into databases...
#define NPlanes_BBGEM 32 // modules...
#define NChan_BBPS 52
#define NChan_BBSH 189
#define NChan_BBHODO 180 
#define NChan_GRINCH 510 
#define NChan_HCAL 288

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
*/



using namespace std;
//____________________________________________________
int main(int argc, char** argv){
  
  // Step 0: read out arguments
  string db_file, inputsigfile, inputbkgdfile = "";//sources of files
  ULong64_t Nentries = -1;//number of events to process
  //UShort_t Nbkgd = 0;//number of background files to add to each event
  double LumiFrac = 0;
      
  if(argc<3){
    cout << "*** Not enough arguments! ***" << endl
	 << " Arguments: database (mandatory); " << endl
	 << "           list_of_sig_input_files (str, mandatory); " << endl
	 << "          nb_of_sig_evts_to_process (int, def=-1); " << endl
	 << "         list_of_bkgd_input_files (str, def=''); " << endl
	 << "        nb_of_bkgd_files_to_add_to_sig_evt (int, def=0); " << endl;
    return(-1);
  }
  
  db_file = argv[1];
  cout << " database file " << db_file << endl;
  inputsigfile = argv[2];
  cout << " Signal input files from: " << inputsigfile << endl;
  if(argc>3)Nentries = atoi(argv[3]);
  cout << " Number of (signal) events to process = " << Nentries << endl;
  if(argc>5){
    inputbkgdfile = argv[4];
    cout << " Background histgrams from: " << inputbkgdfile << endl;
    LumiFrac = max(0., atof(argv[5]));
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
  //f_bkgd->Close();
  
  // ------------------- // dev notes // ------------------- //
  // The loop on the input signal and background chains 
  // is going to happen here in the main I guess.
  //
  // First, we want to extend the input tree (for signal only!!!)
  // I guess in order to avoid adding extra layers of code, 
  // the tree extension might have to be coded in the custom tree class
  
  /*
  double nstrips_bbgem[NPlanes_BBGEM] = {3840, 3072, 3840, 3072, 3840, 3072, 3840, 3072, 3840, 6144};

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
  
  double ZsupThr_bbgem = 240.;
  
  double commonmode_array[1] = {1500.};
  */
  
  std::vector<SBSDigPMTDet*> PMTdetectors;
  std::vector<int> detmap;
  std::vector<SBSDigGEMDet*> GEMdetectors;
  std::vector<SBSDigGEMSimDig*> GEMsimDig;
  std::vector<int> gemdetmap;
  
  // Variable parameters. 
  // Can be configured with the database, but are provided with defaults.
  Int_t Rseed = 0;
  Double_t TriggerJitter = 3.0;
  
  std::vector<TString> detectors_list;
  
  Int_t NChan_BBPS = 52;
  Double_t gatewidth_BBPS = 100.;
  Double_t gain_BBPS = 2.e6;
  Double_t ped_BBPS = 600.;//
  Double_t pedsigma_BBPS = 3.;//
  Double_t trigoffset_BBPS = 18.2;//
  Double_t ADCconv_BBPS = 50.;
  Int_t ADCbits_BBPS = 12;
  Double_t sigmapulse_BBPS = 3.0;
    
  Int_t NChan_BBSH = 189;
  Double_t gatewidth_BBSH = 100.;
  Double_t gain_BBSH = 7.5e5;
  Double_t ped_BBSH = 500.;
  Double_t pedsigma_BBSH = 4.5;
  Double_t trigoffset_BBSH = 18.5;
  Double_t ADCconv_BBSH = 50.;
  Int_t ADCbits_BBSH = 12;
  Double_t sigmapulse_BBSH = 3.0;
    
  Int_t NChan_GRINCH = 510;
  Double_t gatewidth_GRINCH = 100.;
  Double_t gain_GRINCH = 7.e6;
  Double_t ped_GRINCH = 0.;
  Double_t pedsigma_GRINCH = 0.;
  Double_t trigoffset_GRINCH = 15.3;
  Double_t threshold_GRINCH = 3.e-3;
  Double_t ADCconv_GRINCH = 100;
  Int_t ADCbits_GRINCH = 12;
  Double_t TDCconv_GRINCH = 1.;
  Int_t TDCbits_GRINCH = 16;
  Double_t sigmapulse_GRINCH = 3.75;
 
  Int_t NChan_BBHODO = 180;
  Double_t gatewidth_BBHODO = 100.;
  Double_t gain_BBHODO = 1.e5;
  Double_t ped_BBHODO = 0.;
  Double_t pedsigma_BBHODO = 0.;
  Double_t trigoffset_BBHODO = 18.6;
  Double_t threshold_BBHODO = 3.e3;
  Double_t ADCconv_BBHODO = 100.;
  Int_t ADCbits_BBHODO = 12;
  Double_t TDCconv_BBHODO = 0.1;
  Int_t TDCbits_BBHODO = 19;
  Double_t sigmapulse_BBHODO = 1.6;

  Int_t NChan_HCAL = 288;
  Double_t gatewidth_HCAL = 80;
  Double_t gain_HCAL = 1.e6;
  Double_t ped_HCAL = 0.;
  Double_t pedsigma_HCAL = 0.;
  Double_t trigoffset_HCAL = 81.;
  Double_t threshold_HCAL = 3.e-3;
  Double_t ADCconv_HCAL = 1.;
  Double_t TDCconv_HCAL = 0.12;
  Int_t TDCbits_HCAL = 16;
  Int_t FADC_ADCbits = 12;
  Double_t FADC_sampsize = 4.0;
 
  Int_t NPlanes_BBGEM = 32;// 32 // number of planes/modules/readout
  Double_t gatewidth_BBGEM = 400.;
  Int_t* nstrips_bbgem;
  Double_t* offset_bbgem;
  Double_t* strip_angle_bbgem;
  Double_t* triggeroffset_bbgem;
  Double_t ZsupThr_bbgem = 240.;
  Double_t* commonmode_array_bbgem;
  
  // ** How to add a new subsystem **
  // Add param for new detectors there...
  Int_t NChan_POLSCINT_BS = 180;
  Double_t gatewidth_POLSCINT_BS = 100.;
  Double_t gain_POLSCINT_BS = 1.e5;
  Double_t ped_POLSCINT_BS = 0.;
  Double_t pedsigma_POLSCINT_BS = 0.;
  Double_t trigoffset_POLSCINT_BS = 18.6;
  Double_t threshold_POLSCINT_BS = 3.e3;
  Double_t ADCconv_POLSCINT_BS = 100.;
  Int_t ADCbits_POLSCINT_BS = 12;
  Double_t TDCconv_POLSCINT_BS = 0.1;
  Int_t TDCbits_POLSCINT_BS = 19;
  Double_t sigmapulse_POLSCINT_BS = 1.6;
  
  //-----------------------------
  //  Read database
  //-----------------------------
  cout << "read database: " << db_file.c_str() << endl;
  ifstream in_db(db_file.c_str());
  if(!in_db.is_open()){
    cout << "database " << db_file.c_str() << " does not exist!!!" << endl;
    exit(-1);
  }
  
  TString currentline;
  while( currentline.ReadLine(in_db) && !currentline.BeginsWith("endconfig")){
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      int ntokens = tokens->GetEntries();
      
      if( ntokens >= 2 ){
	TString skey = ( (TObjString*) (*tokens)[0] )->GetString();
	
	if(skey=="Rseed"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Rseed = stemp.Atoi();
	}
	
	if(skey=="TriggerJitter"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TriggerJitter = stemp.Atof();
	}
	
	if(skey=="detectors_list"){
	  for(int k = 1; k<ntokens; k++){
	    TString sdet = ( (TObjString*) (*tokens)[k] )->GetString();
	    detectors_list.push_back(sdet);
	  }
	}
	
	//BBPS
	if(skey=="NChan_BBPS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_BBPS = stemp.Atoi();
	}
	
	if(skey=="gatewidth_BBPS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_BBPS = stemp.Atof();
	}
	
	if(skey=="gain_BBPS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_BBPS = stemp.Atof();
	}
	
	if(skey=="ped_BBPS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_BBPS = stemp.Atof();
	}
	
	if(skey=="pedsigma_BBPS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_BBPS = stemp.Atof();
	}
	
	if(skey=="trigoffset_BBPS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_BBPS = stemp.Atof();
	}
	
	if(skey=="ADCconv_BBPS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_BBPS = stemp.Atof();
	}	

	if(skey=="ADCbits_BBPS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_BBPS = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_BBPS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_BBPS = stemp.Atof();
	}
	
	//BBSH
	if(skey=="NChan_BBSH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_BBSH = stemp.Atoi();
	}
	
	if(skey=="gatewidth_BBSH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_BBSH = stemp.Atof();
	}
	
	if(skey=="gain_BBSH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_BBSH = stemp.Atof();
	}
	
	if(skey=="ped_BBSH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_BBSH = stemp.Atof();
	}
	
	if(skey=="pedsigma_BBSH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_BBSH = stemp.Atof();
	}
	
	if(skey=="trigoffset_BBSH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_BBSH = stemp.Atof();
	}
	
	if(skey=="ADCconv_BBSH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_BBSH = stemp.Atof();
	}	

	if(skey=="ADCbits_BBSH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_BBSH = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_BBSH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_BBSH = stemp.Atof();
	}

	//GRINCH
	if(skey=="NChan_GRINCH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_GRINCH = stemp.Atoi();
	}
	
	if(skey=="gatewidth_GRINCH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_GRINCH = stemp.Atof();
	}
	
	if(skey=="gain_GRINCH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_GRINCH = stemp.Atof();
	}
	
	if(skey=="ped_GRINCH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_GRINCH = stemp.Atof();
	}
	
	if(skey=="pedsigma_GRINCH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_GRINCH = stemp.Atof();
	}
	
	if(skey=="trigoffset_GRINCH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_GRINCH = stemp.Atof();
	}
	
	if(skey=="threshold_GRINCH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_GRINCH = stemp.Atof();
	}
	
	if(skey=="ADCconv_GRINCH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_GRINCH = stemp.Atof();
	}	

	if(skey=="ADCbits_GRINCH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_GRINCH = stemp.Atoi();
	}	
	
	if(skey=="TDCconv_GRINCH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_GRINCH = stemp.Atof();
	}	
	
	if(skey=="TDCbits_GRINCH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_GRINCH = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_GRINCH"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_GRINCH = stemp.Atof();
	}
	
	//BBHODO
	if(skey=="NChan_BBHODO"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_BBHODO = stemp.Atoi();
	}
	
	if(skey=="gatewidth_BBHODO"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_BBHODO = stemp.Atof();
	}
	
	if(skey=="gain_BBHODO"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_BBHODO = stemp.Atof();
	}
	
	if(skey=="ped_BBHODO"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_BBHODO = stemp.Atof();
	}
	
	if(skey=="pedsigma_BBHODO"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_BBHODO = stemp.Atof();
	}
	
	if(skey=="trigoffset_BBHODO"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_BBHODO = stemp.Atof();
	}
	
	if(skey=="threshold_BBHODO"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_BBHODO = stemp.Atof();
	}
	
	if(skey=="ADCconv_BBHODO"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_BBHODO = stemp.Atof();
	}	

	if(skey=="ADCbits_BBHODO"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_BBHODO = stemp.Atoi();
	}	
	
	if(skey=="TDCconv_BBHODO"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_BBHODO = stemp.Atof();
	}	
	
	if(skey=="TDCbits_BBHODO"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_BBHODO = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_BBHODO"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_BBHODO = stemp.Atof();
	}
	
	//HCal
	if(skey=="NChan_HCAL"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_HCAL = stemp.Atoi();
	}
	
	if(skey=="gatewidth_HCAL"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_HCAL = stemp.Atof();
	}
	
	if(skey=="gain_HCAL"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_HCAL = stemp.Atof();
	}
	
	if(skey=="ped_HCAL"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_HCAL = stemp.Atof();
	}
	
	if(skey=="pedsigma_HCAL"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_HCAL = stemp.Atof();
	}
	
	if(skey=="trigoffset_HCAL"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_HCAL = stemp.Atof();
	}
	
	if(skey=="threshold_HCAL"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_HCAL = stemp.Atof();
	}
	
	if(skey=="ADCconv_HCAL"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_HCAL = stemp.Atof();
	}	
	
	if(skey=="TDCconv_HCAL"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_HCAL = stemp.Atof();
	}	
	
	if(skey=="TDCbits_HCAL"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_HCAL = stemp.Atoi();
	}	
	
	if(skey=="FADC_ADCbits"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  FADC_ADCbits = stemp.Atoi();
	}
	
	if(skey=="FADC_sampsize"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  FADC_sampsize = stemp.Atof();
	}
	
	// ** How to add a new subsystem **
	// Add reading of param from other detectors there...
	//GEn-RP Hodoscopes
	if(skey=="NChan_POLSCINT_BS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_POLSCINT_BS = stemp.Atoi();
	}
	
	if(skey=="gatewidth_POLSCINT_BS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_POLSCINT_BS = stemp.Atof();
	}
	
	if(skey=="gain_POLSCINT_BS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_POLSCINT_BS = stemp.Atof();
	}
	
	if(skey=="ped_POLSCINT_BS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_POLSCINT_BS = stemp.Atof();
	}
	
	if(skey=="pedsigma_POLSCINT_BS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_POLSCINT_BS = stemp.Atof();
	}
	
	if(skey=="trigoffset_POLSCINT_BS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_POLSCINT_BS = stemp.Atof();
	}
	
	if(skey=="threshold_POLSCINT_BS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_POLSCINT_BS = stemp.Atof();
	}
	
	if(skey=="ADCconv_POLSCINT_BS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_POLSCINT_BS = stemp.Atof();
	}	

	if(skey=="ADCbits_POLSCINT_BS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_POLSCINT_BS = stemp.Atoi();
	}	
	
	if(skey=="TDCconv_POLSCINT_BS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_POLSCINT_BS = stemp.Atof();
	}	
	
	if(skey=="TDCbits_POLSCINT_BS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_POLSCINT_BS = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_POLSCINT_BS"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_POLSCINT_BS = stemp.Atof();
	}
	
	
	//GEMs
	if(skey=="NPlanes_BBGEM"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NPlanes_BBGEM = stemp.Atoi();
	  
	  nstrips_bbgem = new Int_t[NPlanes_BBGEM];
	  offset_bbgem = new Double_t[NPlanes_BBGEM];
	  strip_angle_bbgem = new Double_t[NPlanes_BBGEM];
	  triggeroffset_bbgem = new Double_t[NPlanes_BBGEM];
	  triggeroffset_bbgem = new Double_t[NPlanes_BBGEM];
	}
	
	if(skey=="gatewidth_BBGEM"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_BBGEM = stemp.Atof();
	}
	
	if(skey=="ZsupThr_bbgem"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ZsupThr_bbgem = stemp.Atof();
	}
	
	if(skey=="nstrips_bbgem"){
	  if(ntokens==NPlanes_BBGEM+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      nstrips_bbgem[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for nstrips_bbgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_BBGEM << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="offset_bbgem"){
	  if(ntokens==NPlanes_BBGEM+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      offset_bbgem[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for offset_bbgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_BBGEM << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="strip_angle_bbgem"){
	  if(ntokens==NPlanes_BBGEM+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      strip_angle_bbgem[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for strip_angle_bbgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_BBGEM << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="triggeroffset_bbgem"){
	  if(ntokens==NPlanes_BBGEM/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      triggeroffset_bbgem[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for strip_angle_bbgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_BBGEM << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="commonmode_array_bbgem"){
	  commonmode_array_bbgem = new Double_t[ntokens-1];
	  for(int k = 1; k<ntokens; k++){
	    TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	    commonmode_array_bbgem[k-1] = stemp.Atof();
	  }
	}
	
      }//end if( ntokens >= 2 )
    }//end if( !currentline.BeginsWith("#"))
  }//end while
  
  //-----------------------------
  //  Declare detectors
  //-----------------------------
  cout << " declaring detectors " << endl;
  for(int k = 0; k<detectors_list.size(); k++){
    cout << "detector: " << detectors_list[k].Data() << "... " << endl;
    if(detectors_list[k] == "bbgem"){
      SBSDigGEMDet* bbgem = new SBSDigGEMDet(BBGEM_UNIQUE_DETID, NPlanes_BBGEM, nstrips_bbgem, offset_bbgem, strip_angle_bbgem, 6, ZsupThr_bbgem);
      SBSDigGEMSimDig* gemdig = new SBSDigGEMSimDig(NPlanes_BBGEM, triggeroffset_bbgem, ZsupThr_bbgem, 1, commonmode_array_bbgem);
      
      GEMdetectors.push_back(bbgem);
      gemdetmap.push_back(BBGEM_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "bbps"){
      SBSDigPMTDet* bbps = new SBSDigPMTDet(BBPS_UNIQUE_DETID, NChan_BBPS, gain_BBPS*qe, sigmapulse_BBPS, gatewidth_BBPS);

      bbps->fGain = gain_BBPS;
      bbps->fPedestal = ped_BBPS;
      bbps->fPedSigma = pedsigma_BBPS;
      bbps->fTrigOffset = trigoffset_BBPS;
      bbps->fGateWidth = gatewidth_BBPS;
      bbps->fADCconv = ADCconv_BBPS;
      bbps->fADCbits = ADCbits_BBPS;
      
      PMTdetectors.push_back(bbps);
      detmap.push_back(BBPS_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "bbsh"){
      SBSDigPMTDet* bbsh = new SBSDigPMTDet(BBSH_UNIQUE_DETID, NChan_BBSH, gain_BBSH*qe, sigmapulse_BBPS, gatewidth_BBSH);
      
      bbsh->fGain = gain_BBSH;
      bbsh->fPedestal = ped_BBSH;
      bbsh->fPedSigma = pedsigma_BBSH;
      bbsh->fTrigOffset = trigoffset_BBSH;
      bbsh->fGateWidth = gatewidth_BBSH;
      bbsh->fADCconv = ADCconv_BBSH;
      bbsh->fADCbits = ADCbits_BBSH;    
      
      PMTdetectors.push_back(bbsh);
      detmap.push_back(BBSH_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "grinch"){
      SBSDigPMTDet* grinch = new SBSDigPMTDet(GRINCH_UNIQUE_DETID, NChan_GRINCH, gain_GRINCH*qe, sigmapulse_GRINCH, gatewidth_GRINCH);
  
      grinch->fGain = gain_GRINCH;
      grinch->fPedestal = ped_GRINCH;
      grinch->fPedSigma = pedsigma_GRINCH;
      grinch->fTrigOffset = trigoffset_GRINCH;
      grinch->fThreshold = threshold_GRINCH*spe_unit/ROimpedance;
      grinch->fGateWidth = gatewidth_GRINCH;
      grinch->fADCconv = ADCconv_GRINCH;
      grinch->fADCbits = ADCbits_GRINCH;
      grinch->fTDCconv = TDCconv_GRINCH;
      grinch->fTDCbits = TDCbits_GRINCH;
      
      PMTdetectors.push_back(grinch);
      detmap.push_back(GRINCH_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "bbhodo"){
      SBSDigPMTDet* bbhodo = new SBSDigPMTDet(HODO_UNIQUE_DETID, NChan_BBHODO, gain_BBHODO*qe, sigmapulse_BBHODO, gatewidth_BBHODO);
      
      bbhodo->fGain = gain_BBHODO;
      bbhodo->fPedestal = ped_BBHODO;
      bbhodo->fPedSigma = pedsigma_BBHODO;
      bbhodo->fTrigOffset = trigoffset_BBHODO;
      bbhodo->fThreshold = threshold_BBHODO*spe_unit/ROimpedance;
      bbhodo->fGateWidth = gatewidth_BBHODO;
      bbhodo->fADCconv = ADCconv_BBHODO;
      bbhodo->fADCbits = ADCbits_BBHODO;
      bbhodo->fTDCconv = TDCconv_BBHODO;
      bbhodo->fTDCbits = TDCbits_BBHODO; 
      
      PMTdetectors.push_back(bbhodo);
      detmap.push_back(HODO_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "hcal"){
      SBSDigPMTDet* hcal = new SBSDigPMTDet(HCAL_UNIQUE_DETID, NChan_HCAL);
      
      hcal->fGain = gain_HCAL;
      hcal->fPedestal = ped_HCAL;
      hcal->fPedSigma = pedsigma_HCAL;
      hcal->fTrigOffset = trigoffset_HCAL;
      hcal->fThreshold = threshold_HCAL*spe_unit/ROimpedance;
      hcal->fGateWidth = gatewidth_HCAL;
      hcal->fADCconv = ADCconv_HCAL;
      hcal->fADCbits = FADC_ADCbits;
      hcal->fTDCconv = TDCconv_HCAL;
      hcal->fTDCbits = TDCbits_HCAL; 
      hcal->SetSamples(FADC_sampsize);
      
      //ordered by increasing uinque id
      PMTdetectors.push_back(hcal);
      detmap.push_back(HCAL_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    // ** How to add a new subsystem **
    // Add the new detector here!
    if(detectors_list[k] == "prpolscint_bs"){
      SBSDigPMTDet* polscint_bs = new SBSDigPMTDet(PRPOLBS_SCINT_UNIQUE_DETID, NChan_POLSCINT_BS, gain_POLSCINT_BS*qe, sigmapulse_POLSCINT_BS, gatewidth_POLSCINT_BS);
      
      polscint_bs->fGain = gain_POLSCINT_BS;
      polscint_bs->fPedestal = ped_POLSCINT_BS;
      polscint_bs->fPedSigma = pedsigma_POLSCINT_BS;
      polscint_bs->fTrigOffset = trigoffset_POLSCINT_BS;
      polscint_bs->fThreshold = threshold_POLSCINT_BS*spe_unit/ROimpedance;
      polscint_bs->fGateWidth = gatewidth_POLSCINT_BS;
      polscint_bs->fADCconv = ADCconv_POLSCINT_BS;
      polscint_bs->fADCbits = ADCbits_POLSCINT_BS;
      polscint_bs->fTDCconv = TDCconv_POLSCINT_BS;
      polscint_bs->fTDCbits = TDCbits_POLSCINT_BS; 
      
      PMTdetectors.push_back(polscint_bs);
      detmap.push_back(PRPOLBS_SCINT_UNIQUE_DETID);
      cout << " set up! " << endl;
    } 
  }
  
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
  
  TRandom3* R = new TRandom3(Rseed);
  
  // Step 1: read input files build the input chains
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
  
  //gmn_tree *T_s_;//, *T_b;
  g4sbs_tree *T_s;
  
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
    if(f_s.IsZombie())cout << "File " << chEl_s->GetTitle() << " cannot be found. Please check the path of your file." << endl; 
    run_data = (G4SBSRunData*)f_s.Get("run_data");
    Theta_SBS = run_data->fSBStheta;
    D_HCal = run_data->fHCALdist;
    //TFile fs_c(Form("digitized/simdigtest_%d.root", i_fs), "UPDATE");
    //f_s.Cp(Form("digitized/simdigtest_%d.root", i_fs));
    //if(fs_c.IsOpen())cout << "copy of file is open" << endl;
    //cout << fs_c->ReOpen("UPDATE") << endl;
    //C_s = (TChain*)fs_c.Get("T");
    C_s = (TChain*)f_s.Get("T");
    //T_s = new gmn_tree(C_s);
    T_s = new g4sbs_tree(C_s, detectors_list);
    
    // Expend tree here! (again, for signal only!!!)
    //T_s->AddDigBranches();
    
    Nev_fs = C_s->GetEntries();
    
    for(ev_s = 0; ev_s<Nev_fs; ev_s++, NEventsTotal++){
      if(NEventsTotal>=Nentries)break;
      if(NEventsTotal%1000==0)
	cout << NEventsTotal << "/" << Nentries << endl;
      
      timeZero = R->Gaus(0.0, TriggerJitter);
      
      for(int k = 0; k<PMTdetectors.size(); k++){
	if(detmap[k]==HCAL_UNIQUE_DETID){
	  PMTdetectors[k]->Clear(true);
	}else{
	  PMTdetectors[k]->Clear();
	}
      }
      for(int k = 0; k<GEMdetectors.size(); k++){
	GEMdetectors[k]->Clear();
      }
      /*
      bbgem->Clear();
      //for(int i = 0; i<NPlanes_BBGEM; i++){
      //cout << bbgem->GEMPlanes[i].GetNStrips() << " ";
      //}cout << endl;
      bbps->Clear();
      bbsh->Clear();
      grinch->Clear();
      bbhodo->Clear();
      hcal->Clear(true);
      */
      
      has_data = false;
      
      T_s->ClearDigBranches();
      T_s->GetEntry(ev_s);
      
      // unfold the thing then... but where???
      has_data = UnfoldData(T_s, Theta_SBS, D_HCal, R, PMTdetectors, detmap, GEMdetectors, gemdetmap, timeZero, 0);
      if(!has_data)continue;
      
      
      if(LumiFrac>0){
	BkgdGenerator->GenerateBkgd(R, PMTdetectors, detmap, GEMdetectors, gemdetmap, LumiFrac);
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
	  
	  UnfoldData(T_b, Theta_SBS, D_HCal, R, PMTdetectors, detmap, GEMdetectors, gemdetmap, timeZero, 1);
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
      
      for(int k = 0; k<PMTdetectors.size(); k++){
	PMTdetectors[k]->Digitize(T_s,R);
	
	// bbps->Digitize(T_s,R);
	// bbsh->Digitize(T_s,R);
	// grinch->Digitize(T_s,R);
	// bbhodo->Digitize(T_s,R);
	// hcal->Digitize(T_s,R);
      }
      
      for(int k = 0; k<GEMdetectors.size(); k++){
	GEMsimDig[k]->Digitize(GEMdetectors[k], R);
	GEMsimDig[k]->CheckOut(GEMdetectors[k], gemdetmap[k], R, T_s);
      }
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
