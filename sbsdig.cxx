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
#define NPlanes_bbgem 32 // modules...
#define NChan_bbps 52
#define NChan_bbsh 189
#define NChan_bbhodo 180 
#define NChan_grinch 510 
#define NChan_hcal 288

#define TriggerJitter 3.0 //ns
#define ADCbits 12 
#define gatewidth_PMT 100 //ns
#define gatewidth_GEM 400 //ns

#define FADC_sampsize 4.0 //ns

//DB???
#define sigmapulse_bbpsSH 3.0 //ns / 1.2 of 
#define sigmapulse_bbhodo 1.6 //ns
#define sigmapulse_grinch 3.75 //ns

#define gain_bbps 2.e6
#define ped_bbps 600.0 // ADC channel
#define pedsigma_bbps 3.0 // ADC channel
#define trigoffset_bbps 18.2 //ns
#define ADCconv_bbps 50 //fC/ch

#define gain_bbsh 7.5e5
#define ped_bbsh 500.0 // ADC channel
#define pedsigma_bbsh 4.5 // ADC channel
#define trigoffset_bbsh 18.5 //ns
#define ADCconv_bbsh 50 //fC/ch

#define gain_grinch 7.0e6
#define ped_grinch 0.0 
#define pedsigma_grinch 0.0
#define trigoffset_grinch 15.3 //ns
#define threshold_grinch 3.e-3 //V
#define ADCconv_grinch 100 //fC/ch
#define TDCconv_grinch 1.0 //ns/channel
#define TDCbits_grinch 16 //ns/channel

#define gain_bbhodo 1.0e5
#define ped_bbhodo 0.0 
#define pedsigma_bbhodo 0.0
#define trigoffset_bbhodo 18.6 //ns
#define threshold_bbhodo 3.e-3 //V
#define ADCconv_bbhodo 100 //fC/ch
#define TDCconv_bbhodo 0.1 //ns/channel
#define TDCbits_bbhodo 19 //ns/channel

#define gain_hcal 1.0e6
#define ped_hcal 0.0 
#define pedsigma_hcal 0.0
#define trigoffset_hcal 81.0 //ns
#define threshold_hcal 3.e-3 //V
#define ADCconv_hcal 1.0 //fC/ch //??
#define TDCconv_hcal 0.12 //ns/channel
#define TDCbits_hcal 16 //ns/channel
*/

using namespace std;
//____________________________________________________
int main(int argc, char** argv){
  
  // Step 0: read out arguments
  string db_file, inputsigfile, inputbkgdfile = "";//sources of files
  ULong64_t Nentries = -1;//number of events to process
  //UShort_t Nbkgd = 0;//number of background files to add to each event
  double BkgdTimeWindow = 0, LumiFrac = 0;
      
  if(argc<3 && argc>4){
    cout << "*** Inadequate number of arguments! ***" << endl
	 << " Arguments: database (mandatory); " << endl
	 << "           list_of_sig_input_files (str, mandatory); " << endl
	 << "          nb_of_sig_evts_to_process (int, def=-1); " << endl;
      //<< "         bkgd_histo_input_file (str, def=''); " << endl
      // << "        bkgd_lumi_frac (double, def=0); " << endl;
    return(-1);
  }
  
  db_file = argv[1];
  cout << " database file " << db_file << endl;
  inputsigfile = argv[2];
  cout << " Signal input files from: " << inputsigfile << endl;
  if(argc>3)Nentries = atoi(argv[3]);
  cout << " Number of (signal) events to process = " << Nentries << endl;
  /*
  if(argc>5){
    inputbkgdfile = argv[4];
    cout << " Background histgrams from: " << inputbkgdfile << endl;
    LumiFrac = max(0., atof(argv[5]));
    cout << " Fraction of background to superimpose to signal = " << LumiFrac << endl;
  }
  */
  
  // ------------------- // dev notes // ------------------- //
  // First, we want to extend the input tree (for signal only!!!)
  // I guess in order to avoid adding extra layers of code, 
  // the tree extension might have to be coded in the custom tree class

  
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
  
  Int_t NChan_bbps = 52;
  Double_t gatewidth_bbps = 100.;
  Double_t gain_bbps = 2.e6;
  Double_t ped_bbps = 600.;//
  Double_t pedsigma_bbps = 3.;//
  Double_t trigoffset_bbps = 18.2;//
  Double_t ADCconv_bbps = 50.;
  Int_t ADCbits_bbps = 12;
  Double_t sigmapulse_bbps = 3.0;
    
  Int_t NChan_bbsh = 189;
  Double_t gatewidth_bbsh = 100.;
  Double_t gain_bbsh = 7.5e5;
  Double_t ped_bbsh = 500.;
  Double_t pedsigma_bbsh = 4.5;
  Double_t trigoffset_bbsh = 18.5;
  Double_t ADCconv_bbsh = 50.;
  Int_t ADCbits_bbsh = 12;
  Double_t sigmapulse_bbsh = 3.0;
    
  Int_t NChan_ecal = 2400;
  Double_t gatewidth_ecal = 100.;
  Double_t gain_ecal = 7.5e5;
  Double_t ped_ecal = 500.;
  Double_t pedsigma_ecal = 4.5;
  Double_t trigoffset_ecal = 18.5;
  Double_t ADCconv_ecal = 50.;
  Int_t ADCbits_ecal = 12;
  Double_t sigmapulse_ecal = 3.0;
    
  Int_t NChan_grinch = 510;
  Double_t gatewidth_grinch = 100.;
  Double_t gain_grinch = 7.e6;
  Double_t ped_grinch = 0.;
  Double_t pedsigma_grinch = 0.;
  Double_t trigoffset_grinch = 15.3;
  Double_t threshold_grinch = 3.e-3;
  Double_t ADCconv_grinch = 100;
  Int_t ADCbits_grinch = 12;
  Double_t TDCconv_grinch = 1.;
  Int_t TDCbits_grinch = 16;
  Double_t sigmapulse_grinch = 3.75;
 
  Int_t NChan_bbhodo = 180;
  Double_t gatewidth_bbhodo = 100.;
  Double_t gain_bbhodo = 1.e5;
  Double_t ped_bbhodo = 0.;
  Double_t pedsigma_bbhodo = 0.;
  Double_t trigoffset_bbhodo = 18.6;
  Double_t threshold_bbhodo = 3.e3;
  Double_t ADCconv_bbhodo = 100.;
  Int_t ADCbits_bbhodo = 12;
  Double_t TDCconv_bbhodo = 0.1;
  Int_t TDCbits_bbhodo = 19;
  Double_t sigmapulse_bbhodo = 1.6;
  
  Int_t NChan_cdet = 2352;
  Double_t gatewidth_cdet = 100.;
  Double_t gain_cdet = 1.e5;
  Double_t ped_cdet = 0.;
  Double_t pedsigma_cdet = 0.;
  Double_t trigoffset_cdet = 18.6;
  Double_t threshold_cdet = 3.e3;
  Double_t ADCconv_cdet = 100.;
  Int_t ADCbits_cdet = 12;
  Double_t TDCconv_cdet = 0.1;
  Int_t TDCbits_cdet = 19;
  Double_t sigmapulse_cdet = 1.6; 
  
  Int_t NChan_hcal = 288;
  Double_t gatewidth_hcal = 80;
  Double_t gain_hcal = 1.e6;
  Double_t ped_hcal = 0.;
  Double_t pedsigma_hcal = 0.;
  Double_t trigoffset_hcal = 81.;
  Double_t threshold_hcal = 3.e-3;
  Double_t ADCconv_hcal = 1.;
  Double_t TDCconv_hcal = 0.12;
  Int_t TDCbits_hcal = 16;
  Int_t FADC_ADCbits = 12;
  Double_t FADC_sampsize = 4.0;
 
  Int_t NPlanes_bbgem = 32;// number of planes/modules/readout
  Double_t gatewidth_bbgem = 400.;
  Double_t ZsupThr_bbgem = 240.;
  Int_t Nlayers_bbgem = 5;
  std::vector<Double_t> bbgem_layer_z;
  Int_t* layer_bbgem;
  Int_t* nstrips_bbgem;
  Double_t* offset_bbgem;
  Double_t* RO_angle_bbgem;
  Double_t* triggeroffset_bbgem;
  Double_t* commonmode_array_bbgem;
  UShort_t nAPV_bbgem = 0;
  
  Int_t NPlanes_ft = 36;// number of planes/modules/readout
  Double_t gatewidth_ft = 400.;
  Double_t ZsupThr_ft = 240.;
  Int_t Nlayers_ft = 6;
  std::vector<Double_t> ft_layer_z;
  Int_t* layer_ft;
  Int_t* nstrips_ft;
  Double_t* offset_ft;
  Double_t* RO_angle_ft;
  Double_t* triggeroffset_ft;
  Double_t* commonmode_array_ft;
  UShort_t nAPV_ft = 0;

  Int_t NPlanes_fpp1 = 40;// number of planes/modules/readout
  Double_t gatewidth_fpp1 = 400.;
  Double_t ZsupThr_fpp1 = 240.;
  Int_t Nlayers_fpp1 = 5;
  std::vector<Double_t> fpp1_layer_z;
  Int_t* layer_fpp1;
  Int_t* nstrips_fpp1;
  Double_t* offset_fpp1;
  Double_t* RO_angle_fpp1;
  Double_t* triggeroffset_fpp1;
  Double_t* commonmode_array_fpp1;
  UShort_t nAPV_fpp1 = 0;
  
  Int_t NPlanes_fpp2 = 40;// number of planes/modules/readout
  Double_t gatewidth_fpp2 = 400.;
  Double_t ZsupThr_fpp2 = 240.;
  Int_t Nlayers_fpp2 = 5;
  std::vector<Double_t> fpp2_layer_z;
  Int_t* layer_fpp2;
  Int_t* nstrips_fpp2;
  Double_t* offset_fpp2;
  Double_t* RO_angle_fpp2;
  Double_t* triggeroffset_fpp2;
  Double_t* commonmode_array_fpp2;
  UShort_t nAPV_fpp2 = 0;
    
  // ** How to add a new subsystem **
  // Add param for new detectors there...
  Int_t NChan_polscint_bs = 180;
  Double_t gatewidth_polscint_bs = 100.;
  Double_t gain_polscint_bs = 1.e5;
  Double_t ped_polscint_bs = 0.;
  Double_t pedsigma_polscint_bs = 0.;
  Double_t trigoffset_polscint_bs = 18.6;
  Double_t threshold_polscint_bs = 3.e3;
  Double_t ADCconv_polscint_bs = 100.;
  Int_t ADCbits_polscint_bs = 12;
  Double_t TDCconv_polscint_bs = 0.1;
  Int_t TDCbits_polscint_bs = 19;
  Double_t sigmapulse_polscint_bs = 1.6;
  
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
      Int_t ntokens = 0;
      std::unique_ptr<TObjArray> tokens( currentline.Tokenize(", \t") );
      if( !tokens->IsEmpty() ) {
	ntokens = tokens->GetLast()+1;
      }
      //TObjArray *tokens = currentline.Tokenize(" ");//vg: def lost => versions prior to 6.06; should be fixed! ??? 
      //int ntokens = tokens->GetEntries();
      
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
	if(skey=="NChan_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_bbps = stemp.Atoi();
	}
	
	if(skey=="gatewidth_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_bbps = stemp.Atof();
	}
	
	if(skey=="gain_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_bbps = stemp.Atof();
	}
	
	if(skey=="ped_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_bbps = stemp.Atof();
	}
	
	if(skey=="pedsigma_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_bbps = stemp.Atof();
	}
	
	if(skey=="trigoffset_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_bbps = stemp.Atof();
	}
	
	if(skey=="ADCconv_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_bbps = stemp.Atof();
	}	

	if(skey=="ADCbits_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_bbps = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_bbps = stemp.Atof();
	}
	
	//BBSH
	if(skey=="NChan_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_bbsh = stemp.Atoi();
	}
	
	if(skey=="gatewidth_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_bbsh = stemp.Atof();
	}
	
	if(skey=="gain_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_bbsh = stemp.Atof();
	}
	
	if(skey=="ped_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_bbsh = stemp.Atof();
	}
	
	if(skey=="pedsigma_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_bbsh = stemp.Atof();
	}
	
	if(skey=="trigoffset_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_bbsh = stemp.Atof();
	}
	
	if(skey=="ADCconv_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_bbsh = stemp.Atof();
	}	

	if(skey=="ADCbits_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_bbsh = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_bbsh = stemp.Atof();
	}
	
	//ECAL
	if(skey=="NChan_ecal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_ecal = stemp.Atoi();
	}
	
	if(skey=="gatewidth_ecal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_ecal = stemp.Atof();
	}
	
	if(skey=="gain_ecal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_ecal = stemp.Atof();
	}
	
	if(skey=="ped_ecal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_ecal = stemp.Atof();
	}
	
	if(skey=="pedsigma_ecal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_ecal = stemp.Atof();
	}
	
	if(skey=="trigoffset_ecal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_ecal = stemp.Atof();
	}
	
	if(skey=="ADCconv_ecal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_ecal = stemp.Atof();
	}	

	if(skey=="ADCbits_ecal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_ecal = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_ecal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_ecal = stemp.Atof();
	}
	
	//GRINCH
	if(skey=="NChan_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_grinch = stemp.Atoi();
	}
	
	if(skey=="gatewidth_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_grinch = stemp.Atof();
	}
	
	if(skey=="gain_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_grinch = stemp.Atof();
	}
	
	if(skey=="ped_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_grinch = stemp.Atof();
	}
	
	if(skey=="pedsigma_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_grinch = stemp.Atof();
	}
	
	if(skey=="trigoffset_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_grinch = stemp.Atof();
	}
	
	if(skey=="threshold_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_grinch = stemp.Atof();
	}
	
	if(skey=="ADCconv_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_grinch = stemp.Atof();
	}	

	if(skey=="ADCbits_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_grinch = stemp.Atoi();
	}	
	
	if(skey=="TDCconv_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_grinch = stemp.Atof();
	}	
	
	if(skey=="TDCbits_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_grinch = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_grinch = stemp.Atof();
	}
	
	//BBHODO
	if(skey=="NChan_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_bbhodo = stemp.Atoi();
	}
	
	if(skey=="gatewidth_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_bbhodo = stemp.Atof();
	}
	
	if(skey=="gain_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_bbhodo = stemp.Atof();
	}
	
	if(skey=="ped_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_bbhodo = stemp.Atof();
	}
	
	if(skey=="pedsigma_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_bbhodo = stemp.Atof();
	}
	
	if(skey=="trigoffset_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_bbhodo = stemp.Atof();
	}
	
	if(skey=="threshold_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_bbhodo = stemp.Atof();
	}
	
	if(skey=="ADCconv_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_bbhodo = stemp.Atof();
	}	

	if(skey=="ADCbits_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_bbhodo = stemp.Atoi();
	}	
	
	if(skey=="TDCconv_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_bbhodo = stemp.Atof();
	}	
	
	if(skey=="TDCbits_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_bbhodo = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_bbhodo = stemp.Atof();
	}
	
	//CDET
	if(skey=="NChan_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_cdet = stemp.Atoi();
	}
	
	if(skey=="gatewidth_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_cdet = stemp.Atof();
	}
	
	if(skey=="gain_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_cdet = stemp.Atof();
	}
	
	if(skey=="ped_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_cdet = stemp.Atof();
	}
	
	if(skey=="pedsigma_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_cdet = stemp.Atof();
	}
	
	if(skey=="trigoffset_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_cdet = stemp.Atof();
	}
	
	if(skey=="threshold_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_cdet = stemp.Atof();
	}
	
	if(skey=="ADCconv_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_cdet = stemp.Atof();
	}	

	if(skey=="ADCbits_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_cdet = stemp.Atoi();
	}	
	
	if(skey=="TDCconv_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_cdet = stemp.Atof();
	}	
	
	if(skey=="TDCbits_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_cdet = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_cdet = stemp.Atof();
	}
	
	//HCal
	if(skey=="NChan_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_hcal = stemp.Atoi();
	}
	
	if(skey=="gatewidth_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_hcal = stemp.Atof();
	}
	
	if(skey=="gain_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_hcal = stemp.Atof();
	}
	
	if(skey=="ped_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_hcal = stemp.Atof();
	}
	
	if(skey=="pedsigma_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_hcal = stemp.Atof();
	}
	
	if(skey=="trigoffset_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_hcal = stemp.Atof();
	}
	
	if(skey=="threshold_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_hcal = stemp.Atof();
	}
	
	if(skey=="ADCconv_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_hcal = stemp.Atof();
	}	
	
	if(skey=="TDCconv_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_hcal = stemp.Atof();
	}	
	
	if(skey=="TDCbits_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_hcal = stemp.Atoi();
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
	if(skey=="NChan_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_polscint_bs = stemp.Atoi();
	}
	
	if(skey=="gatewidth_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_polscint_bs = stemp.Atof();
	}
	
	if(skey=="gain_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_polscint_bs = stemp.Atof();
	}
	
	if(skey=="ped_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_polscint_bs = stemp.Atof();
	}
	
	if(skey=="pedsigma_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_polscint_bs = stemp.Atof();
	}
	
	if(skey=="trigoffset_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_polscint_bs = stemp.Atof();
	}
	
	if(skey=="threshold_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_polscint_bs = stemp.Atof();
	}
	
	if(skey=="ADCconv_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_polscint_bs = stemp.Atof();
	}	

	if(skey=="ADCbits_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_polscint_bs = stemp.Atoi();
	}	
	
	if(skey=="TDCconv_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_polscint_bs = stemp.Atof();
	}	
	
	if(skey=="TDCbits_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_polscint_bs = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_polscint_bs = stemp.Atof();
	}
	
	//GEMs
	if(skey=="NPlanes_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NPlanes_bbgem = stemp.Atoi();
	  
	  layer_bbgem = new Int_t[NPlanes_bbgem];
	  nstrips_bbgem = new Int_t[NPlanes_bbgem];
	  offset_bbgem = new Double_t[NPlanes_bbgem];
	  RO_angle_bbgem = new Double_t[NPlanes_bbgem];
	  triggeroffset_bbgem = new Double_t[NPlanes_bbgem/2];
	}
	
	if(skey=="gatewidth_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_bbgem = stemp.Atof();
	}
		
	if(skey=="ZsupThr_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ZsupThr_bbgem = stemp.Atof();
	}

	if(skey=="nlayers_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Nlayers_bbgem = stemp.Atof();
	}
	
	if(skey=="bbgem_layer_z"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  if(ntokens==Nlayers_bbgem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      bbgem_layer_z.push_back(stemp.Atof());
	    }
	  }else{
	    cout << "number of entries for bbgem_layer_z = " << ntokens-1 
		 << " don't match nlayers = " << Nlayers_bbgem << endl;
	    cout << "fix your db " << endl;
	  }
	}
	
	if(skey=="layer_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_bbgem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      layer_bbgem[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for layer_bbgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_bbgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="nstrips_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_bbgem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      nstrips_bbgem[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for nstrips_bbgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_bbgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="offset_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_bbgem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      offset_bbgem[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for offset_bbgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_bbgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="RO_angle_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_bbgem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      RO_angle_bbgem[k-1] = stemp.Atof()*TMath::DegToRad();
	    }
	  }else{
	    cout << "number of entries for RO_angle_bbgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_bbgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="triggeroffset_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_bbgem/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      triggeroffset_bbgem[k-1] = stemp.Atof();
	      // if this is affected at k-1, program crashes... it should not and that's confusing...
	    }
	  }else{
	    cout << "number of entries for triggeroffset_bbgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_bbgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="commonmode_array_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  commonmode_array_bbgem = new Double_t[ntokens-1];
	  for(int k = 1; k<ntokens; k++){
	    nAPV_bbgem++;
	    TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	    commonmode_array_bbgem[k-1] = stemp.Atof();
	  }
	}

	//FT
	if(skey=="NPlanes_ft"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NPlanes_ft = stemp.Atoi();
	  
	  layer_ft = new Int_t[NPlanes_ft];
	  nstrips_ft = new Int_t[NPlanes_ft];
	  offset_ft = new Double_t[NPlanes_ft];
	  RO_angle_ft = new Double_t[NPlanes_ft];
	  triggeroffset_ft = new Double_t[NPlanes_ft/2];
	}
	
	if(skey=="gatewidth_ft"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_ft = stemp.Atof();
	}
	
	if(skey=="ZsupThr_ft"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ZsupThr_ft = stemp.Atof();
	}

	if(skey=="nlayers_ft"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Nlayers_ft = stemp.Atof();
	}
	
	if(skey=="ft_layer_z"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  if(ntokens==Nlayers_ft+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      ft_layer_z.push_back(stemp.Atof());
	    }
	  }else{
	    cout << "number of entries for ft_layer_z = " << ntokens-1 
		 << " don't match nlayers = " << Nlayers_ft << endl;
	    cout << "fix your db " << endl;
	  }
	}
		
	if(skey=="layer_ft"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_ft+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      layer_ft[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for layer_ft = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_ft << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="nstrips_ft"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_ft+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      nstrips_ft[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for nstrips_ft = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_ft << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="offset_ft"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_ft+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      offset_ft[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for offset_ft = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_ft << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="RO_angle_ft"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_ft+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      RO_angle_ft[k-1] = stemp.Atof()*TMath::DegToRad();
	    }
	  }else{
	    cout << "number of entries for RO_angle_ft = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_ft << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="triggeroffset_ft"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_ft/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      triggeroffset_ft[k-1] = stemp.Atof();
	      // if this is affected at k-1, program crashes... it should not and that's confusing...
	    }
	  }else{
	    cout << "number of entries for triggeroffset_ft = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_ft << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="commonmode_array_ft"){
	  cout << "reading " << skey.Data() << endl;
	  commonmode_array_ft = new Double_t[ntokens-1];
	  for(int k = 1; k<ntokens; k++){
	    nAPV_ft++;
	    TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	    commonmode_array_ft[k-1] = stemp.Atof();
	  }
	}
	
	//FPP1
	if(skey=="NPlanes_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NPlanes_fpp1 = stemp.Atoi();
	  
	  layer_fpp1 = new Int_t[NPlanes_fpp1];
	  nstrips_fpp1 = new Int_t[NPlanes_fpp1];
	  offset_fpp1 = new Double_t[NPlanes_fpp1];
	  RO_angle_fpp1 = new Double_t[NPlanes_fpp1];
	  triggeroffset_fpp1 = new Double_t[NPlanes_fpp1/2];
	}
	
	if(skey=="gatewidth_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_fpp1 = stemp.Atof();
	}
	
	if(skey=="ZsupThr_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ZsupThr_fpp1 = stemp.Atof();
	}
	
	if(skey=="nlayers_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Nlayers_fpp1 = stemp.Atof();
	}
	
	if(skey=="fpp1_layer_z"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  if(ntokens==Nlayers_fpp1+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      fpp1_layer_z.push_back(stemp.Atof());
	    }
	  }else{
	    cout << "number of entries for fpp1_layer_z = " << ntokens-1 
		 << " don't match nlayers = " << Nlayers_fpp1 << endl;
	    cout << "fix your db " << endl;
	  }
	}
	
	if(skey=="layer_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp1+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      layer_fpp1[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for layer_fpp1 = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_fpp1 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="nstrips_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp1+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      nstrips_fpp1[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for nstrips_fpp1 = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_fpp1 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="offset_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp1+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      offset_fpp1[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for offset_fpp1 = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_fpp1 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="RO_angle_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp1+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      RO_angle_fpp1[k-1] = stemp.Atof()*TMath::DegToRad();
	    }
	  }else{
	    cout << "number of entries for RO_angle_fpp1 = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_fpp1 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="triggeroffset_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp1/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      triggeroffset_fpp1[k-1] = stemp.Atof();
	      // if this is affected at k-1, program crashes... it should not and that's confusing...
	    }
	  }else{
	    cout << "number of entries for triggeroffset_fpp1 = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_fpp1 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="commonmode_array_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  commonmode_array_fpp1 = new Double_t[ntokens-1];
	  for(int k = 1; k<ntokens; k++){
	    nAPV_fpp1++;
	    TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	    commonmode_array_fpp1[k-1] = stemp.Atof();
	  }
	}
	
	//FPP2
	if(skey=="NPlanes_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NPlanes_fpp2 = stemp.Atoi();
	  
	  layer_fpp2 = new Int_t[NPlanes_fpp2];
	  nstrips_fpp2 = new Int_t[NPlanes_fpp2];
	  offset_fpp2 = new Double_t[NPlanes_fpp2];
	  RO_angle_fpp2 = new Double_t[NPlanes_fpp2];
	  triggeroffset_fpp2 = new Double_t[NPlanes_fpp2/2];
	}
	
	if(skey=="gatewidth_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_fpp2 = stemp.Atof();
	}
	
	if(skey=="ZsupThr_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ZsupThr_fpp2 = stemp.Atof();
	}
	
	if(skey=="nlayers_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Nlayers_fpp2 = stemp.Atof();
	}
	
	if(skey=="fpp2_layer_z"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  if(ntokens==Nlayers_fpp2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      fpp2_layer_z.push_back(stemp.Atof());
	    }
	  }else{
	    cout << "number of entries for fpp2_layer_z = " << ntokens-1 
		 << " don't match nlayers = " << Nlayers_fpp2 << endl;
	    cout << "fix your db " << endl;
	  }
	}
	
	if(skey=="layer_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      layer_fpp2[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for layer_fpp2 = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_fpp2 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="nstrips_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      nstrips_fpp2[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for nstrips_fpp2 = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_fpp2 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="offset_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      offset_fpp2[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for offset_fpp2 = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_fpp2 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="RO_angle_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      RO_angle_fpp2[k-1] = stemp.Atof()*TMath::DegToRad();
	    }
	  }else{
	    cout << "number of entries for RO_angle_fpp2 = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_fpp2 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="triggeroffset_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp2/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      triggeroffset_fpp2[k-1] = stemp.Atof();
	      // if this is affected at k-1, program crashes... it should not and that's confusing...
	    }
	  }else{
	    cout << "number of entries for triggeroffset_fpp2 = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_fpp2 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="commonmode_array_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  commonmode_array_fpp2 = new Double_t[ntokens-1];
	  for(int k = 1; k<ntokens; k++){
	    nAPV_fpp2++;
	    TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	    commonmode_array_fpp2[k-1] = stemp.Atof();
	  }
	}
	
      }//end if( ntokens >= 2 )
      tokens->~TObjArray();// ineffective... :(
    }//end if( !currentline.BeginsWith("#"))
  }//end while
  
  //-----------------------------
  //  Declare detectors
  //-----------------------------
  cout << " declaring detectors " << endl;
  for(int k = 0; k<detectors_list.size(); k++){
    cout << "detector: " << detectors_list[k].Data() << "... " << endl;
    if(detectors_list[k] == "bbgem"){
      SBSDigGEMDet* bbgem = new SBSDigGEMDet(BBGEM_UNIQUE_DETID, NPlanes_bbgem, layer_bbgem, nstrips_bbgem, offset_bbgem, RO_angle_bbgem, 6, ZsupThr_bbgem);
      SBSDigGEMSimDig* gemdig = new SBSDigGEMSimDig(NPlanes_bbgem/2, triggeroffset_bbgem, ZsupThr_bbgem, nAPV_bbgem, commonmode_array_bbgem);
      for(int m = 0; m<Nlayers_bbgem; m++){
	bbgem->fZLayer.push_back(bbgem_layer_z[m]);
      }
      bbgem->fGateWidth = gatewidth_bbgem;
      
      GEMdetectors.push_back(bbgem);
      gemdetmap.push_back(BBGEM_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "ft"){
      SBSDigGEMDet* ft = new SBSDigGEMDet(FT_UNIQUE_DETID, NPlanes_ft, layer_ft, nstrips_ft, offset_ft, RO_angle_ft, 6, ZsupThr_ft);
      SBSDigGEMSimDig* gemdig = new SBSDigGEMSimDig(NPlanes_ft/2, triggeroffset_ft, ZsupThr_ft, nAPV_ft, commonmode_array_ft);
      for(int m = 0; m<Nlayers_ft; m++){
	ft->fZLayer.push_back(ft_layer_z[m]);
      }
      ft->fGateWidth = gatewidth_ft;
      
      GEMdetectors.push_back(ft);
      gemdetmap.push_back(FT_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "fpp1"){
      SBSDigGEMDet* fpp1 = new SBSDigGEMDet(FPP1_UNIQUE_DETID, NPlanes_fpp1, layer_fpp1, nstrips_fpp1, offset_fpp1, RO_angle_fpp1, 6, ZsupThr_fpp1);
      SBSDigGEMSimDig* gemdig = new SBSDigGEMSimDig(NPlanes_fpp1/2, triggeroffset_fpp1, ZsupThr_fpp1, nAPV_fpp1, commonmode_array_fpp1);
      for(int m = 0; m<Nlayers_fpp1; m++){
	fpp1->fZLayer.push_back(fpp1_layer_z[m]);
      }
      fpp1->fGateWidth = gatewidth_fpp1;
      
      GEMdetectors.push_back(fpp1);
      gemdetmap.push_back(FPP1_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "fpp2"){
      SBSDigGEMDet* fpp2 = new SBSDigGEMDet(FPP2_UNIQUE_DETID, NPlanes_fpp2, layer_fpp2, nstrips_fpp2, offset_fpp2, RO_angle_fpp2, 6, ZsupThr_fpp2);
      SBSDigGEMSimDig* gemdig = new SBSDigGEMSimDig(NPlanes_fpp2/2, triggeroffset_fpp2, ZsupThr_fpp2, nAPV_fpp2, commonmode_array_fpp2);
      for(int m = 0; m<Nlayers_fpp2; m++){
	fpp2->fZLayer.push_back(fpp2_layer_z[m]);
      }
      fpp2->fGateWidth = gatewidth_fpp2;
      
      GEMdetectors.push_back(fpp2);
      gemdetmap.push_back(FPP2_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "bbps"){
      SBSDigPMTDet* bbps = new SBSDigPMTDet(BBPS_UNIQUE_DETID, NChan_bbps, gain_bbps*qe, sigmapulse_bbps, gatewidth_bbps);

      bbps->fGain = gain_bbps;
      bbps->fPedestal = ped_bbps;
      bbps->fPedSigma = pedsigma_bbps;
      bbps->fTrigOffset = trigoffset_bbps;
      bbps->fGateWidth = gatewidth_bbps;
      bbps->fADCconv = ADCconv_bbps;
      bbps->fADCbits = ADCbits_bbps;
      
      PMTdetectors.push_back(bbps);
      detmap.push_back(BBPS_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "bbsh"){
      SBSDigPMTDet* bbsh = new SBSDigPMTDet(BBSH_UNIQUE_DETID, NChan_bbsh, gain_bbsh*qe, sigmapulse_bbsh, gatewidth_bbsh);
      
      bbsh->fGain = gain_bbsh;
      bbsh->fPedestal = ped_bbsh;
      bbsh->fPedSigma = pedsigma_bbsh;
      bbsh->fTrigOffset = trigoffset_bbsh;
      bbsh->fGateWidth = gatewidth_bbsh;
      bbsh->fADCconv = ADCconv_bbsh;
      bbsh->fADCbits = ADCbits_bbsh;    
      
      PMTdetectors.push_back(bbsh);
      detmap.push_back(BBSH_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "ecal"){
      SBSDigPMTDet* ecal = new SBSDigPMTDet(ECAL_UNIQUE_DETID, NChan_ecal, gain_ecal*qe, sigmapulse_ecal, gatewidth_ecal);
      
      ecal->fGain = gain_ecal;
      ecal->fPedestal = ped_ecal;
      ecal->fPedSigma = pedsigma_ecal;
      ecal->fTrigOffset = trigoffset_ecal;
      ecal->fGateWidth = gatewidth_ecal;
      ecal->fADCconv = ADCconv_ecal;
      ecal->fADCbits = ADCbits_ecal;    
      
      PMTdetectors.push_back(ecal);
      detmap.push_back(ECAL_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "grinch"){
      SBSDigPMTDet* grinch = new SBSDigPMTDet(GRINCH_UNIQUE_DETID, NChan_grinch, gain_grinch*qe, sigmapulse_grinch, gatewidth_grinch);
  
      grinch->fGain = gain_grinch;
      grinch->fPedestal = ped_grinch;
      grinch->fPedSigma = pedsigma_grinch;
      grinch->fTrigOffset = trigoffset_grinch;
      grinch->fThreshold = threshold_grinch*spe_unit/ROimpedance;
      grinch->fGateWidth = gatewidth_grinch;
      grinch->fADCconv = ADCconv_grinch;
      grinch->fADCbits = ADCbits_grinch;
      grinch->fTDCconv = TDCconv_grinch;
      grinch->fTDCbits = TDCbits_grinch;
      
      PMTdetectors.push_back(grinch);
      detmap.push_back(GRINCH_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "bbhodo"){
      SBSDigPMTDet* bbhodo = new SBSDigPMTDet(HODO_UNIQUE_DETID, NChan_bbhodo, gain_bbhodo*qe, sigmapulse_bbhodo, gatewidth_bbhodo);
      
      bbhodo->fGain = gain_bbhodo;
      bbhodo->fPedestal = ped_bbhodo;
      bbhodo->fPedSigma = pedsigma_bbhodo;
      bbhodo->fTrigOffset = trigoffset_bbhodo;
      bbhodo->fThreshold = threshold_bbhodo*spe_unit/ROimpedance;
      bbhodo->fGateWidth = gatewidth_bbhodo;
      bbhodo->fADCconv = ADCconv_bbhodo;
      bbhodo->fADCbits = ADCbits_bbhodo;
      bbhodo->fTDCconv = TDCconv_bbhodo;
      bbhodo->fTDCbits = TDCbits_bbhodo; 
      
      PMTdetectors.push_back(bbhodo);
      detmap.push_back(HODO_UNIQUE_DETID);
      cout << " set up! " << endl;
    }

    if(detectors_list[k] == "h_cdet" || detectors_list[k] == "e_cdet"){
      SBSDigPMTDet* cdet = new SBSDigPMTDet(CDET_UNIQUE_DETID, NChan_cdet, gain_cdet*qe, sigmapulse_cdet, gatewidth_cdet);
      
      cdet->fGain = gain_cdet;
      cdet->fPedestal = ped_cdet;
      cdet->fPedSigma = pedsigma_cdet;
      cdet->fTrigOffset = trigoffset_cdet;
      cdet->fThreshold = threshold_cdet*spe_unit/ROimpedance;
      cdet->fGateWidth = gatewidth_cdet;
      cdet->fADCconv = ADCconv_cdet;
      cdet->fADCbits = ADCbits_cdet;
      cdet->fTDCconv = TDCconv_cdet;
      cdet->fTDCbits = TDCbits_cdet; 
      
      PMTdetectors.push_back(cdet);
      detmap.push_back(CDET_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "hcal"){
      SBSDigPMTDet* hcal = new SBSDigPMTDet(HCAL_UNIQUE_DETID, NChan_hcal);
      
      hcal->fGain = gain_hcal;
      hcal->fPedestal = ped_hcal;
      hcal->fPedSigma = pedsigma_hcal;
      hcal->fTrigOffset = trigoffset_hcal;
      hcal->fThreshold = threshold_hcal*spe_unit/ROimpedance;
      hcal->fGateWidth = gatewidth_hcal;
      hcal->fADCconv = ADCconv_hcal;
      hcal->fADCbits = FADC_ADCbits;
      hcal->fTDCconv = TDCconv_hcal;
      hcal->fTDCbits = TDCbits_hcal; 
      hcal->SetSamples(FADC_sampsize);
      
      //ordered by increasing uinque id
      PMTdetectors.push_back(hcal);
      detmap.push_back(HCAL_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    // ** How to add a new subsystem **
    // Add the new detector here!
    if(detectors_list[k] == "prpolscint_bs"){
      SBSDigPMTDet* polscint_bs = new SBSDigPMTDet(PRPOLBS_SCINT_UNIQUE_DETID, NChan_polscint_bs, gain_polscint_bs*qe, sigmapulse_polscint_bs, gatewidth_polscint_bs);
      
      polscint_bs->fGain = gain_polscint_bs;
      polscint_bs->fPedestal = ped_polscint_bs;
      polscint_bs->fPedSigma = pedsigma_polscint_bs;
      polscint_bs->fTrigOffset = trigoffset_polscint_bs;
      polscint_bs->fThreshold = threshold_polscint_bs*spe_unit/ROimpedance;
      polscint_bs->fGateWidth = gatewidth_polscint_bs;
      polscint_bs->fADCconv = ADCconv_polscint_bs;
      polscint_bs->fADCbits = ADCbits_polscint_bs;
      polscint_bs->fTDCconv = TDCconv_polscint_bs;
      polscint_bs->fTDCbits = TDCbits_polscint_bs; 
      
      PMTdetectors.push_back(polscint_bs);
      detmap.push_back(PRPOLBS_SCINT_UNIQUE_DETID);
      cout << " set up! " << endl;
    } 
  }
  
  //restart reading to get background info
  while( currentline.ReadLine(in_db) && !currentline.BeginsWith("end_bkgdinfo")){
    if( !currentline.BeginsWith("#") ){
      Int_t ntokens = 0;
      std::unique_ptr<TObjArray> tokens( currentline.Tokenize(", \t") );
      if( !tokens->IsEmpty() ) {
	ntokens = tokens->GetLast()+1;
      }
      if(ntokens==3){
	inputbkgdfile = ( (TObjString*) (*tokens)[0] )->GetString();
	TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	BkgdTimeWindow = stemp.Atof();
	stemp = ( (TObjString*) (*tokens)[2] )->GetString();
	LumiFrac = stemp.Atof();
      }
    }
  }
  
  if(LumiFrac){
    cout << "Includes background from file: " << inputbkgdfile.c_str() 
	 << " (integrated on " << BkgdTimeWindow << " ns time window);" 
	 << endl << " assuming " << LumiFrac*100 << "% luminosity."<< endl;
  }
  
  TFile* f_bkgd;
  SBSDigBkgdGen* BkgdGenerator;
  if(LumiFrac>0){
    f_bkgd = TFile::Open(inputbkgdfile.c_str());
    if(f_bkgd->IsZombie()){
      LumiFrac = 0;
    }else{
      BkgdGenerator = new SBSDigBkgdGen(f_bkgd, BkgdTimeWindow);
    }
  }
  //f_bkgd->Close();
  
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
  
  //G4SBSRunData* run_data;
  
  double Theta_SBS, D_HCal;
  
  //gmn_tree *T_s_;//, *T_b;
  //g4sbs_tree *T_s;
  
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
    //run_data = (G4SBSRunData*)f_s.Get("run_data");
    G4SBSRunData* run_data = (G4SBSRunData*)f_s.Get("run_data");
    Theta_SBS = run_data->fSBStheta;
    D_HCal = run_data->fHCALdist;
    //TFile fs_c(Form("digitized/simdigtest_%d.root", i_fs), "UPDATE");
    //f_s.Cp(Form("digitized/simdigtest_%d.root", i_fs));
    //if(fs_c.IsOpen())cout << "copy of file is open" << endl;
    //cout << fs_c->ReOpen("UPDATE") << endl;
    //C_s = (TChain*)fs_c.Get("T");
    C_s = (TChain*)f_s.Get("T");
    //T_s = new gmn_tree(C_s);
    //T_s = new g4sbs_tree(C_s, detectors_list);//vg: def lost
    g4sbs_tree *T_s = new g4sbs_tree(C_s, detectors_list);
    //g4sbs_tree T_s(C_s, detectors_list);
    
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
      //for(int i = 0; i<NPlanes_bbgem; i++){
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
      //T_s.ClearDigBranches();
      //T_s.GetEntry(ev_s);
      
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
      //T_s.FillDigBranches();
      //T_s->fChain->Fill();
    }// end loop on signal events 
    
    for(int k = 0; k<GEMdetectors.size(); k++){
      GEMsimDig[k]->write_histos();
    }
    /**/
    if(LumiFrac>0)BkgdGenerator->WriteXCHistos();
    T_s->fChain->Write("", TObject::kOverwrite);
    //T_s.fChain->Write("", TObject::kOverwrite);
    //fs_c.Write();
    //fs_c.Close();
    f_s.Write();
    f_s.Close();
    //T_s->~g4sbs_tree();
    //T_s.~g4sbs_tree();
    i_fs++;
  }// end loop on signal files
  
  
  exit(0);
}
