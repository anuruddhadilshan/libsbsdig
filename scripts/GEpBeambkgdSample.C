#include "TMatrixD.h"
#include "TVectorD.h"
#include "TFile.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
#include "TDecompSVD.h"
#include "TCut.h"
#include "TEventList.h"
#include "G4SBSRunData.hh"
#include "gep_deftree.C"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TMath.h"
#include "TF2.h"
//#include "HistLoader.h" do we need this?
#include "TChainElement.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TSystem.h"


#include <set>
// This script is needed to generate histograms of distribution 
// of energy deposit, position, etc... for beam induced background hits.

void BeamBackground(const char *inputfilename = "list_gep1_beambkgd.txt", 
		    const char *root_outfile = "gep1_beambkgd.root")
		    
{
  cout << " This file uses root generated class gep_deftree. " << endl
       << " if you can't find it or your file seems to be inadapted to your file " << endl
       << " (missing branches warnings, unidentified crashes)" << endl
       << " you may generate file gmn_deftree.C and .h yourself " << endl
       << " by opening one of your g4sbs files and type: " 
       << " T->MakeClass(\"gmn_deftree\") " << endl;
  
  // reads and stores the list of files to process
  cout << "reading input files" << endl;
  
  ifstream inputfile(inputfilename);
  //TFile *fout = new TFile( outputfilename, "RECREATE" );
  TChain *C = new TChain("T");

  set<TString> files;
  
  TString currentline;
  while( currentline.ReadLine(inputfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline.Data());

      files.insert( currentline );
    }
  }
  
  //const double Thr_Hodo = 8.0e-3;// threshold to ensure the slat with max energy deposit is recorded
  
  //string part_str[12] = {"#pi^{0}", "K^{-}", "#pi^{-}", "#mu^{-}", "e^{-}", "#gamma",
  //"e^{+}", "#mu^{+}", "#pi^{+}", "K^{+}", "n", "p"};
  
  //log binning for momentum: 
  double pbins_log10[102];
  for(int k = 0; k<102; k++){
    pbins_log10[k] = pow(10.0, (double)k*0.05-4.0);
  }
  
  double edepbins_log10[121];
  for(int k = 0; k<121; k++){
    edepbins_log10[k] = pow(10.0, (double)k*0.05-6.0);
  }

  double pebins_log10[81];
  for(int k = 0; k<81; k++){
    pebins_log10[k] = pow(10.0, (double)k*0.05);
  }
  
  // declare the output files and the histograms
  // 
  TFile *fout = new TFile( root_outfile, "RECREATE" );
  
  TH1D *h1_Ntries = new TH1D("h1_Ntries","number of total generated events",1, 0, 1);
  
  //GEMs
  TH1D* h1_FT_nhits_[8];
  
  TH2D* h1_FT_yVsx_[8];
  TH2D* h1_FT_dyVsdx_[8];
  TH1D* h1_FT_Edep_[8];
  TH1D* h1_FT_Edep_log_[8];
  
  TH1D* h1_FPP1_nhits_[8];
  
  TH2D* h1_FPP1_yVsx_[8];
  TH2D* h1_FPP1_dyVsdx_[8];
  TH1D* h1_FPP1_Edep_[8];
  TH1D* h1_FPP1_Edep_log_[8];
  
  for(int m = 0; m<8; m++){
    h1_FT_nhits_[m] = new TH1D(Form("h1_FT_nhits_%d",m), 
				  Form("BB GEM plane %d energy deposit", m+1),
				  500, 0.0, 500);
    h1_FT_yVsx_[m] = new TH2D(Form("h1_FT_yVsx_%d",m), 
  				 Form("BB GEM plane %d y vs x coordinate", m+1),
  				 202, -1.01, 1.01, 62, -0.31, 0.31);
    h1_FT_dyVsdx_[m] = new TH2D(Form("h1_FT_dyVsdx_%d",m), 
				   Form("BB GEM plane %d dy vs dx coordinate", m+1),
				   100, -0.1, 0.1, 100, -0.1, 0.1);
    h1_FT_Edep_[m] = new TH1D(Form("h1_FT_Edep_%d",m), 
  				 Form("BB GEM plane %d energy deposit;MeV", m+1),
  				 100, 0.0, 1.e-1);
    h1_FT_Edep_log_[m] = new TH1D(Form("h1_FT_Edep_log_%d",m), 
				     Form("BB GEM plane %d energy deposit;MeV", m+1),
				     110, edepbins_log10);

    h1_FPP1_nhits_[m] = new TH1D(Form("h1_FPP1_nhits_%d",m), 
				  Form("BB GEM plane %d energy deposit", m+1),
				  500, 0.0, 500);
    h1_FPP1_yVsx_[m] = new TH2D(Form("h1_FPP1_yVsx_%d",m), 
  				 Form("BB GEM plane %d y vs x coordinate", m+1),
  				 202, -1.01, 1.01, 62, -0.31, 0.31);
    h1_FPP1_dyVsdx_[m] = new TH2D(Form("h1_FPP1_dyVsdx_%d",m), 
				   Form("BB GEM plane %d dy vs dx coordinate", m+1),
				   100, -0.1, 0.1, 100, -0.1, 0.1);
    h1_FPP1_Edep_[m] = new TH1D(Form("h1_FPP1_Edep_%d",m), 
  				 Form("BB GEM plane %d energy deposit;MeV", m+1),
  				 100, 0.0, 1.e-1);
    h1_FPP1_Edep_log_[m] = new TH1D(Form("h1_FPP1_Edep_log_%d",m), 
				     Form("BB GEM plane %d energy deposit;MeV", m+1),
				     110, edepbins_log10);
  }

  //HCal
  TH2D *h1_HCal_nhitsVsChan = new TH2D("h1_HCal_nhitsVsChan", "", 288, 0, 288, 100, 0.0, 100);
  
  TH2D *h1_HCal_EdepHitVsChan = new TH2D("h1_HCal_EdepHitVsChan", ";GeV", 288, 0, 288, 100, 0.0+1.0e-3, 1.0+1.0e-3);
  TH2D *h1_HCal_EdepHitVsChan_log = new TH2D("h1_HCal_EdepHitVsChan_log", ";GeV", 288, 0, 288, 110, edepbins_log10);
  TH2D *h1_HCal_zHitVsChan = new TH2D("h1_HCal_zHitVsChan", ";m", 288, 0, 288, 1000, 0.0, 1.0);
  
  TH2D *h1_HCal_EdepTotVsChan = new TH2D("h1_HCal_EdepTotVsChan;GeV", ";GeV", 288, 0, 288, 1000, 0.0, 2.0);
  TH2D *h1_HCal_EdepTotVsChan_log = new TH2D("h1_HCal_EdepTotVsChan_log", ";GeV", 288, 0, 288, 110, edepbins_log10);
  
  //CDet
  TH2D *h1_CDet_nhitsVsSlat = new TH2D("h1_CDet_nhitsVsSlat", "", 2352, 0, 2352, 100, 0, 100);
  
  TH2D *h1_CDet_xhitVsSlat = new TH2D("h1_CDet_xhitVsSlat", "", 2352, 0, 2352, 60, -0.3, 0.3);
  
  TH2D *h1_CDet_EdepHitVsSlat = new TH2D("h1_CDet_EdepHitVsSlat", ";GeV", 2352, 0, 2352, 200, 0.0, 0.2);
  TH2D *h1_CDet_EdepHitVsSlat_log = new TH2D("h1_CDet_EdepHitVsSlat_log", ";GeV", 2352, 0, 2352, 110, edepbins_log10);
  
  TH2D *h1_CDet_EdepTotVsSlat = new TH2D("h1_CDet_EdepTotVsSlat", ";GeV", 2352, 0, 2352, 250, 0.0, 0.5);
  TH2D *h1_CDet_EdepTotVsSlat_log = new TH2D("h1_CDet_EdepTotVsSlat_log", ";GeV", 2352, 0, 2352, 110, edepbins_log10);

  //ECal
  TH2D *h1_ECal_nhitsVsChan = new TH2D("h1_ECal_nhitsVsChan", "", 1656, 0, 1656, 100, 0.0, 100);
  
  TH2D *h1_ECal_EdepHitVsChan = new TH2D("h1_ECal_EdepHitVsChan", ";GeV", 1656, 0, 1656, 100, 0.0+1.0e-3, 1.0+1.0e-3);
  TH2D *h1_ECal_EdepHitVsChan_log = new TH2D("h1_ECal_EdepHitVsChan_log", ";GeV", 1656, 0, 1656, 110, edepbins_log10);
  TH2D *h1_ECal_zHitVsChan = new TH2D("h1_ECal_zHitVsChan", ";m", 1656, 0, 1656, 1000, 0.0, 1.0);
  
  TH2D *h1_ECal_EdepTotVsChan = new TH2D("h1_ECal_EdepTotVsChan;GeV", ";GeV", 1656, 0, 1656, 1000, 0.0, 2.0);
  TH2D *h1_ECal_EdepTotVsChan_log = new TH2D("h1_ECal_EdepTotVsChan_log", ";GeV", 1656, 0, 1656, 110, edepbins_log10);

  
  int FileCounter = 0;
  
  // declare arrays for the number of hits
  int nhits_FT[8];
  int nhits_FPP1[8];
  int nhits_HCal[288];
  int nhits_CDet[2352];
  int nhits_ECal[1656];
  
  memset(nhits_FT, 0, 8*sizeof(int));
  memset(nhits_FPP1, 0, 8*sizeof(int));
  memset(nhits_HCal, 0, 288*sizeof(int));
  memset(nhits_CDet, 0, 2352*sizeof(int));
  memset(nhits_ECal, 0, 1656*sizeof(int));

  Long64_t N1_tot = 0;
  int nplane;
  int chan;

  double theta_sbs, d_hcal;
  double x_ref, z_ref;
  double z_hit;
  
  // loop on files
  TObjArray *fileElements=C->GetListOfFiles();
  TIter next(fileElements);
  TChainElement *chEl=0;
  while (( chEl=(TChainElement*)next() )) {
    //TFile *f = new TFile(chEl->GetTitle());
    TFile f(chEl->GetTitle());
     
    // for each file, get the G4SBSRunData object...
    G4SBSRunData *RD = (G4SBSRunData*)f.Get("run_data");
    N1_tot+= (double)RD->fNtries;
    theta_sbs = RD->fSBStheta;
    d_hcal = RD->fHCALdist;
    x_ref = -d_hcal*sin(theta_sbs);
    z_ref = d_hcal*cos(theta_sbs);
    
    TChain *C1 = (TChain*)f.Get("T");
    
    // and the tree...
    // the tree has to be generated by: TTree::MakeClass(...) method.
    // eg. the tree class below was generated with:
    // root gep_beambkgd_job_0.root
    // root [0] T->MakeClass("gep_deftree");
    gep_deftree *T1 = new gep_deftree(C1);
    
    FileCounter++;
  
    Long64_t MaxEvt = C1->GetEntries();
    cout << chEl->GetTitle() << ": " << MaxEvt << " entries" << endl;
    // loop on events:
    for(Long64_t nevent = 0; nevent<MaxEvt; nevent++){
      //while( T1->GetEntry(nevent++) && nevent < 10000 ){
      if( nevent%1000 == 0){// && nevent!=0){
	cout << nevent << endl;
      }
      
      T1->GetEntry(nevent);
      
      // fill HCal hit energy deposit and position histograms, and increment number of hits
      if(T1->Harm_HCalScint_hit_nhits){
	for(int i = 0; i<T1->Harm_HCalScint_hit_nhits; i++){
	  if(T1->Harm_HCalScint_hit_sumedep->at(i)>=1.e-3){
	    nhits_HCal[T1->Harm_HCalScint_hit_cell->at(i)]++;
	  
	    h1_HCal_EdepHitVsChan->Fill(T1->Harm_HCalScint_hit_cell->at(i), T1->Harm_HCalScint_hit_sumedep->at(i));
	    h1_HCal_EdepHitVsChan_log->Fill(T1->Harm_HCalScint_hit_cell->at(i), T1->Harm_HCalScint_hit_sumedep->at(i));
	    
	    z_hit = -(T1->Harm_HCalScint_hit_xhitg->at(i)-x_ref)*sin(theta_sbs)+(T1->Harm_HCalScint_hit_zhitg->at(i)-z_ref)*cos(theta_sbs);
	    
	    h1_HCal_zHitVsChan->Fill(T1->Harm_HCalScint_hit_cell->at(i), z_hit);
	  }
	  // Edep_HCal+= T1->Harm_HCalScint_hit_sumedep->at(i);
	  // h1_HCal_blocks_Edep_rates->Fill(T1->Harm_HCalScint_hit_cell->at(i), T1->Harm_HCalScint_hit_sumedep->at(i));  
	  // //in MeV
	  // Edep_tot_50ns+= T1->Harm_HCalScint_hit_sumedep->at(i)*1.e3;
	  // Edep_tot_blocks_50ns[T1->Harm_HCalScint_hit_cell->at(i)]+= T1->Harm_HCalScint_hit_sumedep->at(i)*1.e3;
	  
	  // h1_HCal_Edep_tot_XC->Fill(T1->Harm_HCalScint_hit_cell->at(i), T1->Harm_HCalScint_hit_sumedep->at(i));
	}
      }
      
      // fill ECal energy deposit histograms, and increment number of hits
      if(T1->Earm_ECalTF1_hit_nhits){
	for(int i = 0; i<T1->Earm_ECalTF1_hit_nhits; i++){
	  nhits_ECal[T1->Earm_ECalTF1_hit_cell->at(i)]++;
	  h1_ECal_EdepHitVsChan->Fill(T1->Earm_ECalTF1_hit_cell->at(i), T1->Earm_ECalTF1_hit_sumedep->at(i));
	  h1_ECal_EdepHitVsChan_log->Fill(T1->Earm_ECalTF1_hit_cell->at(i), T1->Earm_ECalTF1_hit_sumedep->at(i));
	}
      }
      
      // fill CDET energy deposit and position histograms, and increment number of hits
      if(T1->Earm_CDET_Scint_hit_nhits){
	for(int i = 0; i<T1->Earm_CDET_Scint_hit_nhits; i++){
	  nhits_CDet[T1->Earm_CDET_Scint_hit_cell->at(i)]++;
	  h1_CDet_xhitVsSlat->Fill(T1->Earm_CDET_Scint_hit_cell->at(i), T1->Earm_CDET_Scint_hit_xhit->at(i));
	  h1_CDet_EdepHitVsSlat->Fill(T1->Earm_CDET_Scint_hit_cell->at(i), T1->Earm_CDET_Scint_hit_sumedep->at(i));
	  h1_CDet_EdepHitVsSlat_log->Fill(T1->Earm_CDET_Scint_hit_cell->at(i), T1->Earm_CDET_Scint_hit_sumedep->at(i));
	}
      }
      
      // fill GEMs energy deposit, position, and drift exit-entrance position difference histograms, and increment number of hits
      //cout << "FT: " << T1->Earm_FT_hit_nhits << endl;
      if(T1->Harm_FT_hit_nhits){
	for(int i = 0; i<T1->Harm_FT_hit_nhits; i++){
	  nplane = T1->Harm_FT_hit_plane->at(i)-1;
	  
	  nhits_FT[nplane]++;
	  // 
	  
	  // //cout << T1->Harm_FT_hit_tx->at(i) << " " << T1->Harm_FT_hit_ty->at(i) << endl;
	  
	  h1_FT_yVsx_[nplane]->Fill(T1->Harm_FT_hit_xin->at(i), 
	   			       T1->Harm_FT_hit_yin->at(i));
	  h1_FT_dyVsdx_[nplane]->Fill(T1->Harm_FT_hit_xout->at(i)-T1->Harm_FT_hit_xin->at(i), 
					 T1->Harm_FT_hit_yout->at(i)-T1->Harm_FT_hit_yin->at(i));
	  
	  h1_FT_Edep_[nplane]->Fill(T1->Harm_FT_hit_edep->at(i)*1.e3);
	  h1_FT_Edep_log_[nplane]->Fill(T1->Harm_FT_hit_edep->at(i)*1.e3);
	  
	}
      }
      
      if(T1->Harm_FPP1_hit_nhits){
	for(int i = 0; i<T1->Harm_FPP1_hit_nhits; i++){
	  nplane = T1->Harm_FPP1_hit_plane->at(i)-1;
	  
	  nhits_FPP1[nplane]++;
	  // 
	  
	  // //cout << T1->Harm_FPP1_hit_tx->at(i) << " " << T1->Harm_FPP1_hit_ty->at(i) << endl;
	  
	  h1_FPP1_yVsx_[nplane]->Fill(T1->Harm_FPP1_hit_xin->at(i), 
	   			       T1->Harm_FPP1_hit_yin->at(i));
	  h1_FPP1_dyVsdx_[nplane]->Fill(T1->Harm_FPP1_hit_xout->at(i)-T1->Harm_FPP1_hit_xin->at(i), 
					 T1->Harm_FPP1_hit_yout->at(i)-T1->Harm_FPP1_hit_yin->at(i));
	  
	  h1_FPP1_Edep_[nplane]->Fill(T1->Harm_FPP1_hit_edep->at(i)*1.e3);
	  h1_FPP1_Edep_log_[nplane]->Fill(T1->Harm_FPP1_hit_edep->at(i)*1.e3);
	  
	}
      }
    }//end for(Long_64t...)
    
    cout << " Background sample per file: " << RD->fNtries << " events per file " << endl;
    for(int j = 0; j<288; j++){
      h1_HCal_nhitsVsChan->Fill(j, nhits_HCal[j]);
    }
    
    for(int j = 0; j<2352; j++){
      h1_CDet_nhitsVsSlat->Fill(j, nhits_CDet[j]);
    }
    
    for(int j = 0; j<1656; j++){
      h1_ECal_nhitsVsChan->Fill(j, nhits_ECal[j]);
    }
    
    memset(nhits_HCal, 0, 288*sizeof(int));
    memset(nhits_CDet, 0, 2352*sizeof(int));
    memset(nhits_ECal, 0, 1656*sizeof(int));

    for(int j = 0; j<8; j++){
      h1_FT_nhits_[j]->Fill(nhits_FT[j]);
      h1_FPP1_nhits_[j]->Fill(nhits_FPP1[j]);
    }
    memset(nhits_FT, 0, 8*sizeof(int));
    memset(nhits_FPP1, 0, 8*sizeof(int));
    
    C1->Delete();
    T1->~gep_deftree();
  }//end while (( chEl= )) 
  h1_Ntries->Fill(0.5, N1_tot);
  
  
  //h1_HCal_Edep_rates->Scale(I_exp/1.602e-19/N1_tot);

  C->Delete();
  fout->Write();
  //return;
}
