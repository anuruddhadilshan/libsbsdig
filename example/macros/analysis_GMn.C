#include "TMatrixD.h"
#include "TVectorD.h"
#include "TFile.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "TString.h"
#include "TArray.h"
#include "TDecompSVD.h"
#include "TCut.h"
#include "TEventList.h"
#include "digtree.C"
//#include "gep_tree_with_spin.C"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TMath.h"
#include "TF2.h"
//#include "HistLoader.h"
#include "TChainElement.h"
#include "TRandom3.h"
#include "TGraph.h"

//This macro is for the moment oversimplified, and is just to check roughly the addition of background for the digitization in the GRINCH and hosdoscope.
//It will be expanded in a very near future.


void analysis_GMn(const char *inputfilename, 
		  const char *root_outfile)
{
  // load the files that we want to analyze
  cout << "reading input files" << endl;
  
  ifstream inputfile(inputfilename);
  TFile *fout = new TFile( root_outfile, "RECREATE" );
  TChain *C = new TChain("T");
  
  set<TString> files;
  
  TString currentline;
  while( currentline.ReadLine(inputfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline.Data());

      files.insert( currentline );
    }
  }
  
  inputfile.close();
  
  //Declaring histos
  TH1D* h1_grinchMCtime = new TH1D("h1_grinchMCtime", "Grinch;MC_time (ns)", 200, -100, 100);
  TH1D* h1_hodoMCtime = new TH1D("h1_hodoMCtime", "Hodoscope;MC_time (ns)", 200, -100, 100);

  int file_num = 0;
    
  // --------------------------
  // Loop on files
  // --------------------------
  TObjArray *fileElements=C->GetListOfFiles();
  TIter next(fileElements);
  TChainElement *chEl=0;
  while (( chEl=(TChainElement*)next() )) {
    //TFile *f = new TFile(chEl->GetTitle());
    TFile f(chEl->GetTitle());
    
    TChain *C1 = (TChain*)f.Get("digtree");
    digtree *T1 = new digtree(C1);
    
    Long64_t NEvts1 = C1->GetEntries();
    
    cout << chEl->GetTitle() << ": " << NEvts1 << " entries" << endl;
    
    // --------------------------
    // Loop on events
    // --------------------------
    for(Long64_t nevent = 0; nevent<NEvts1; nevent++){
      if( nevent%1000 == 0 ){
	cout << nevent << endl;
      }
      T1->GetEntry(nevent);
      //Plot MC time for GRINCH as an example
      if(T1->NSimData_bb_grinch){
	     
	for(int i = 0; i<T1->NSimData_bb_grinch; i++){
	  if(T1->SimData_bb_grinch_Type->at(i)==2){
	    for(int j = 0; j<T1->SimData_bb_grinch_Ndata->at(j); j++){
	      h1_grinchMCtime->Fill( (T1->SimData_bb_grinch_Data->at(i))[j] );
	    }
	  }
	}
      }

      if(T1->NSimData_bb_hodo){
	
	for(int i = 0; i<T1->NSimData_bb_hodo; i++){
	  if(T1->SimData_bb_hodo_Type->at(i)==2){
	    for(int j = 0; j<T1->SimData_bb_hodo_Ndata->at(j); j++){
	      h1_hodoMCtime->Fill( (T1->SimData_bb_hodo_Data->at(i))[j] );
	    }
	  }
	}
      }
    }//end loop on events
    
    
    
    file_num++;
    
    //if(file_num>=Nfiles)break;
    
  }//end loop on files
  
  TCanvas* C1 = new TCanvas("C1", "C1", 1000, 500);
  C1->Divide(2,1);
  C1->cd(1);
  h1_grinchMCtime->Draw("");
  C1->cd(2);
  h1_hodoMCtime->Draw("");
  
  fout->Write();
  
  C->Delete();
  //histloader1.CloseAll();
}
