//includes: standard
#include <iostream>
#include <fstream>
#include <string>

//includes: root
#include <TROOT.h>
#include "TString.h"
#include "TChain.h"
#include "TChainElement.h"

//includes: specific
//#include "G4SBSRunData.hh"
#include "gmn_tree.h"

#ifdef __APPLE__
#include "unistd.h"
#endif

using namespace std;
//____________________________________________________
int main(int argc, char** argv){

  // Step 0: read out arguments
  string inputsigfile, inputbkgdfile = "";//sources of files
  ULong64_t Nentries = -1;//number of events to process
  UShort_t Nbkgd = 0;//number of background files to add to each event
  
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
    cout << " Background input files from: " << inputbkgdfile << endl;
    Nbkgd = atoi(argv[4]);
    cout << " Number of background files to superimpose to signal = " << Nbkgd << endl;
  }
  
  // ------------------- // dev notes // ------------------- //
  // The loop on the input signal and background chains 
  // is going to happen here in the main I guess.
  //
  // First, we want to extend the input tree (for signal only!!!)
  // I guess in order to avoid adding extra layers of code, 
  // the tree extension might have to be coded in the custom tree class
  
  // Step 1: read input files build the input chains
  TString currentline;
  
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
  
  // build background chain
  ifstream beam_inputfile(inputbkgdfile);
  TChain *C_b = new TChain("T");
  if(Nbkgd!=0){
    while( currentline.ReadLine(beam_inputfile) && !currentline.BeginsWith("endlist") ){
      if( !currentline.BeginsWith("#") ){
	C_b->Add(currentline.Data());
      }
    }
  }
  TObjArray *fileElements_b=C_b->GetListOfFiles();
  TIter next_b(fileElements_b);
  TChainElement *chEl_b=0;
  
  
  gmn_tree *T_s, *T_b;
  Long64_t Nev_fs, Nev_fb;
  Long64_t ev_s, ev_b;
  
  Long64_t NEventsTotal = 0;
  UShort_t nbkgd = 0;
  
  while (( chEl_s=(TChainElement*)next_s() )) {
    if(NEventsTotal>=Nentries)break;
    
    TFile f_s(chEl_s->GetTitle());
    C_s = (TChain*)f_s.Get("T");
    T_s = new gmn_tree(C_s);
    Nev_fs = C_s->GetEntries();
    
    cout << chEl_s->GetTitle() << ": " << Nev_fs << " entries" << endl;
    
    // Expend tree here! (again, signal only!!!)
    
    for(ev_s = 0; ev_s<Nev_fs; ev_s++, NEventsTotal++){
      if(NEventsTotal>=Nentries)break;
      if(NEventsTotal%1000==0)cout << NEventsTotal << "/" << Nentries << endl;
      
      T_s->GetEntry(ev_s);
      
      // unfold the thing then...
      
      // loop here for background
      nbkgd = 0;
      while (( chEl_b=(TChainElement*)next_b() )) {
	if(nbkgd>=Nbkgd)break;
	
	TFile f_b(chEl_b->GetTitle());
	C_b = (TChain*)f_b.Get("T");
	T_b = new gmn_tree(C_b);
	Nev_fb = C_b->GetEntries();
	
	cout << chEl_b->GetTitle() << ": " << Nev_fb << " entries" << endl;
	
	for(ev_b = 0; ev_b<Nev_fb; ev_b++){
	  
	  
	  T_b->GetEntry(ev_b);
	}// end loop on background events
	nbkgd++;
      }// end loop on background files
      
    }// end loop on signal events 
    
  }// end loop on signal files
  
  
}

