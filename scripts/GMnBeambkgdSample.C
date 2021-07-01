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
#include "gmn_deftree.C"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TMath.h"
#include "TF2.h"
#include "HistLoader.h"
#include "TChainElement.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TSystem.h"

// This script is needed to generate histograms of distribution 
// of energy deposit, position, etc... for beam induced background hits.

void BeamBackground(const char *inputfilename, 
		    const char *root_outfile)
		    
{
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
  TH1D* h1_BBGEM_nhits_[5];
  
  TH2D* h1_BBGEM_yVsx_[5];
  TH2D* h1_BBGEM_dyVsdx_[5];
  TH1D* h1_BBGEM_Edep_[5];
  TH1D* h1_BBGEM_Edep_log_[5];
  
  for(int m = 0; m<5; m++){
    h1_BBGEM_nhits_[m] = new TH1D(Form("h1_BBGEM_nhits_%d",m), 
				  Form("BB GEM plane %d energy deposit", m+1),
				  500, 0.0, 500);
    h1_BBGEM_yVsx_[m] = new TH2D(Form("h1_BBGEM_yVsx_%d",m), 
  				 Form("BB GEM plane %d y vs x coordinate", m+1),
  				 202, -1.01, 1.01, 62, -0.31, 0.31);
    h1_BBGEM_dyVsdx_[m] = new TH2D(Form("h1_BBGEM_dyVsdx_%d",m), 
				   Form("BB GEM plane %d dy vs dx coordinate", m+1),
				   100, -0.1, 0.1, 100, -0.1, 0.1);
    h1_BBGEM_Edep_[m] = new TH1D(Form("h1_BBGEM_Edep_%d",m), 
  				 Form("BB GEM plane %d energy deposit;MeV", m+1),
  				 100, 0.0, 1.e-1);
    h1_BBGEM_Edep_log_[m] = new TH1D(Form("h1_BBGEM_Edep_log_%d",m), 
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
  
  //Hodoscope
  TH2D *h1_BBHodo_nhitsVsSlat = new TH2D("h1_BBHodo_nhitsVsSlat", "", 90, 0, 90, 100, 0, 100);
  
  TH2D *h1_BBHodo_xhitVsSlat = new TH2D("h1_BBHodo_xhitVsSlat", "", 90, 0, 90, 60, -0.3, 0.3);
  
  TH2D *h1_BBHodo_EdepHitVsSlat = new TH2D("h1_BBHodo_EdepHitVsSlat", ";GeV", 90, 0, 90, 200, 0.0, 0.2);
  TH2D *h1_BBHodo_EdepHitVsSlat_log = new TH2D("h1_BBHodo_EdepHitVsSlat_log", ";GeV", 90, 0, 90, 110, edepbins_log10);
  
  TH2D *h1_BBHodo_EdepTotVsSlat = new TH2D("h1_BBHodo_EdepTotVsSlat", ";GeV", 90, 0, 90, 250, 0.0, 0.5);
  TH2D *h1_BBHodo_EdepTotVsSlat_log = new TH2D("h1_BBHodo_EdepTotVsSlat_log", ";GeV", 90, 0, 90, 110, edepbins_log10);

  //PS
  TH2D *h1_BBPS_nhitsVsChan = new TH2D("h1_BBPS_nhitsVsChan", "", 52, 0, 52, 150, 0, 150);
  
  TH2D *h1_BBPS_EdepHitVsChan = new TH2D("h1_BBPS_EdepHitVsChan", ";GeV", 52, 0, 52, 250, 0.0+1.0e-3, 2.5+1.0e-3);
  TH2D *h1_BBPS_EdepHitVsChan_log = new TH2D("h1_BBPS_EdepHitVsChan_log", ";GeV", 52, 0, 52, 110, edepbins_log10);
  
  //SH
  TH2D *h1_BBSH_nhitsVsChan = new TH2D("h1_BBSH_nhitsVsChan", "", 189, 0, 189, 100, 0, 100);
  
  TH2D *h1_BBSH_EdepHitVsChan = new TH2D("h1_BBSH_EdepHitVsChan", ";GeV", 189, 0, 189, 250+1.0e-3, 0.0, 2.5+1.0e-3);
  TH2D *h1_BBSH_EdepHitVsChan_log = new TH2D("h1_BBSH_EdepHitVsChan_log", ";GeV", 189, 0, 189, 110, edepbins_log10);
  
  //GRINCH
  TH2D *h1_GRINCH_nhitsVsChan = new TH2D("h1_GRINCH_nhitsVsChan", "", 510, 0, 510, 20, 0, 20);
   
  TH2D *h1_GRINCH_NpeVsChan = new TH2D("h1_GRINCH_NpeVsChan", "", 510, 0, 510, 100, 0, 100);
  
  int FileCounter = 0;
  
  // declare arrays for the number of hits
  int nhits_GEM[5];
  int nhits_HCal[288];
  int nhits_Hodo[90];
  int nhits_PS[52];
  int nhits_SH[189];
  int nhits_GRINCH[510];
    
  double HCal_Edep_blocks_400ns[288];
  double BBHodo_Edep_slats_400ns[90];
  
  memset(nhits_GEM, 0, 5*sizeof(int));
  memset(nhits_HCal, 0, 288*sizeof(int));
  memset(nhits_Hodo, 0, 90*sizeof(int));
  memset(nhits_PS, 0, 52*sizeof(int));
  memset(nhits_SH, 0, 189*sizeof(int));
  memset(nhits_GRINCH, 0, 510*sizeof(int));

  memset(HCal_Edep_blocks_400ns, 0, 288*sizeof(double));
  memset(BBHodo_Edep_slats_400ns, 0, 90*sizeof(double));

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
    //N1_tot+= (double)RD->fNtries;
    theta_sbs = RD->fSBStheta;
    d_hcal = RD->fHCALdist;
    x_ref = -d_hcal*sin(theta_sbs);
    z_ref = d_hcal*cos(theta_sbs);
    
    TChain *C1 = (TChain*)f.Get("T");
    
    // and the tree...
    gmn_deftree *T1 = new gmn_deftree(C1);
    
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
	  HCal_Edep_blocks_400ns[T1->Harm_HCalScint_hit_cell->at(i)]+= T1->Harm_HCalScint_hit_sumedep->at(i);
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
      
      // fill BB hodoscope energy deposit and position histograms, and increment number of hits
      if(T1->Earm_BBHodoScint_hit_nhits){
	for(int i = 0; i<T1->Earm_BBHodoScint_hit_nhits; i++){
	  nhits_Hodo[T1->Earm_BBHodoScint_hit_cell->at(i)]++;
	  BBHodo_Edep_slats_400ns[T1->Earm_BBHodoScint_hit_cell->at(i)]+= T1->Earm_BBHodoScint_hit_sumedep->at(i);
	  
	  h1_BBHodo_xhitVsSlat->Fill(T1->Earm_BBHodoScint_hit_cell->at(i), T1->Earm_BBHodoScint_hit_xhit->at(i));
	  h1_BBHodo_EdepHitVsSlat->Fill(T1->Earm_BBHodoScint_hit_cell->at(i), T1->Earm_BBHodoScint_hit_sumedep->at(i));
	  h1_BBHodo_EdepHitVsSlat_log->Fill(T1->Earm_BBHodoScint_hit_cell->at(i), T1->Earm_BBHodoScint_hit_sumedep->at(i));
	  // h1_BBHodo_slat_rates_nothr->Fill(T1->Earm_BBHodoScint_hit_cell->at(i));
	  // if(T1->Earm_BBHodoScint_hit_sumedep->at(i)>=Thr_Hodo){
	  //   h1_BBHodo_slat_rates->Fill(T1->Earm_BBHodoScint_hit_cell->at(i));
	  // }
	  // h1_BBHodo_slat_Edep_rates->Fill(T1->Earm_BBHodoScint_hit_cell->at(i), T1->Earm_BBHodoScint_hit_sumedep->at(i));
	  
	  // Nph_hodo_1 = Npe_E_GeV_bbhs*T1->Earm_BBHodoScint_hit_sumedep->at(i)*LCE_0*exp(-(0.3+(T1->Earm_BBHodoScint_hit_xhit->at(i)-T1->Earm_BBHodoScint_hit_xcellg->at(i))/cos(theta_BB))/LCE_lambda);
	  // Nph_hodo_2 = Npe_E_GeV_bbhs*T1->Earm_BBHodoScint_hit_sumedep->at(i)*LCE_0*exp(-(0.3-(T1->Earm_BBHodoScint_hit_xhit->at(i)-T1->Earm_BBHodoScint_hit_xcellg->at(i))/cos(theta_BB))/LCE_lambda);
	  
	  // Npe_hodo_1 = R.Poisson(Nph_hodo_1*0.24);
	  // //0.28 is the QE, 0.1 is a facto to take into account of the variable photon path due to the solid angle (very conservative)
	  // Npe_hodo_2 = R.Poisson(Nph_hodo_2*0.24);
	  
	  // if(Npe_hodo_1)h1_BBHodo_slat_pe_rates->Fill(2*T1->Earm_BBHodoScint_hit_cell->at(i), Npe_hodo_1);
	  // if(Npe_hodo_2)h1_BBHodo_slat_pe_rates->Fill(2*T1->Earm_BBHodoScint_hit_cell->at(i)+1, Npe_hodo_2);
	}
      }
      
      // fill BigBite GEMs energy deposit, position, and drift exit-entrance position difference histograms, and increment number of hits
      //cout << "BBGEM: " << T1->Earm_BBGEM_hit_nhits << endl;
      if(T1->Earm_BBGEM_hit_nhits){
	for(int i = 0; i<T1->Earm_BBGEM_hit_nhits; i++){
	  nplane = T1->Earm_BBGEM_hit_plane->at(i)-1;
	  
	  nhits_GEM[nplane]++;
	  // 
	  
	  // //cout << T1->Earm_BBGEM_hit_tx->at(i) << " " << T1->Earm_BBGEM_hit_ty->at(i) << endl;
	  
	  h1_BBGEM_yVsx_[nplane]->Fill(T1->Earm_BBGEM_hit_xin->at(i), 
	   			       T1->Earm_BBGEM_hit_yin->at(i));
	  h1_BBGEM_dyVsdx_[nplane]->Fill(T1->Earm_BBGEM_hit_xout->at(i)-T1->Earm_BBGEM_hit_xin->at(i), 
					 T1->Earm_BBGEM_hit_yout->at(i)-T1->Earm_BBGEM_hit_yin->at(i));
	  
	  h1_BBGEM_Edep_[nplane]->Fill(T1->Earm_BBGEM_hit_edep->at(i)*1.e3);
	  h1_BBGEM_Edep_log_[nplane]->Fill(T1->Earm_BBGEM_hit_edep->at(i)*1.e3);
	  
	}
      }
      
      // fill Preshower energy deposit histograms, and increment number of hits
      if(T1->Earm_BBPSTF1_hit_nhits){
	for(int i = 0; i<T1->Earm_BBPSTF1_hit_nhits; i++){
	  if(T1->Earm_BBPSTF1_hit_sumedep->at(i)>=1.e-3){
	    nhits_PS[T1->Earm_BBPSTF1_hit_cell->at(i)]++;
	    
	    h1_BBPS_EdepHitVsChan->Fill(T1->Earm_BBPSTF1_hit_cell->at(i), T1->Earm_BBPSTF1_hit_sumedep->at(i));
	    h1_BBPS_EdepHitVsChan_log->Fill(T1->Earm_BBPSTF1_hit_cell->at(i), T1->Earm_BBPSTF1_hit_sumedep->at(i));
	  }
	  // h1_BBPS_Edep->Fill(T1->Earm_BBPSTF1_hit_row->at(i), 
	  // 		     T1->Earm_BBPSTF1_hit_col->at(i), 
	  // 		     T1->Earm_BBPSTF1_hit_sumedep->at(i));
	  // h1_BBPS_EdepSpectrum->Fill(T1->Earm_BBPSTF1_hit_cell->at(i), T1->Earm_BBPSTF1_hit_sumedep->at(i));
	}
      }
      //cout << "BBSHTF1: " << T1->Earm_BBSHTF1_hit_nhits << endl;
      // fill Shower energy deposit histograms
      if(T1->Earm_BBSHTF1_hit_nhits){
	for(int i = 0; i<T1->Earm_BBSHTF1_hit_nhits; i++){
	  if(T1->Earm_BBSHTF1_hit_sumedep->at(i)>=1.e-3){
	    nhits_SH[T1->Earm_BBSHTF1_hit_cell->at(i)]++;
	    
	    h1_BBSH_EdepHitVsChan->Fill(T1->Earm_BBSHTF1_hit_cell->at(i), T1->Earm_BBSHTF1_hit_sumedep->at(i));
	    h1_BBSH_EdepHitVsChan_log->Fill(T1->Earm_BBSHTF1_hit_cell->at(i), T1->Earm_BBSHTF1_hit_sumedep->at(i));
	  }
	  // h1_BBSH_Edep->Fill(T1->Earm_BBSHTF1_hit_row->at(i), 
	  // 		     T1->Earm_BBSHTF1_hit_col->at(i), 
	  // 		     T1->Earm_BBSHTF1_hit_sumedep->at(i));
	  // h1_BBSH_EdepSpectrum->Fill(T1->Earm_BBSHTF1_hit_cell->at(i), T1->Earm_BBSHTF1_hit_sumedep->at(i));
	}
      }
      
      // fill GRINCH histograms, and increment number of hits
      if(T1->Earm_GRINCH_hit_nhits){
	for(int i = 0; i<T1->Earm_GRINCH_hit_nhits; i++){
	  chan = int(T1->Earm_GRINCH_hit_PMT->at(i)/5);
	  nhits_GRINCH[chan]++;

	  h1_GRINCH_NpeVsChan->Fill(chan, T1->Earm_GRINCH_hit_NumPhotoelectrons->at(i));
	  // h1_BBSH_Edep->Fill(T1->Earm_BBSHTF1_hit_row->at(i), 
	  // 		     T1->Earm_BBSHTF1_hit_col->at(i), 
	  // 		     T1->Earm_BBSHTF1_hit_sumedep->at(i));
	  // h1_BBSH_EdepSpectrum->Fill(T1->Earm_BBSHTF1_hit_cell->at(i), T1->Earm_BBSHTF1_hit_sumedep->at(i));
	}
      }
     
      
    }//end for(Long_64t...)
    
    //if(FileCounter==5){
    // fill the hits histograms and reinitialize the number of hits 
    // WARNING: the following assumes that each beam-on-target file has 5 million generated beam electrons
    // TODO: replace this with something smarter (e.g. set a "time window" and a beam current).
    if(FileCounter%3==0){
      cout << "File Counter = 3 => 15 M evts beam-on-target" << endl;
      for(int j = 0; j<288; j++){
    	h1_HCal_nhitsVsChan->Fill(j, nhits_HCal[j]);
    	h1_HCal_EdepTotVsChan->Fill(j, HCal_Edep_blocks_400ns[j]);
    	h1_HCal_EdepTotVsChan_log->Fill(j, HCal_Edep_blocks_400ns[j]);
      }
      
      for(int j = 0; j<90; j++){
    	h1_BBHodo_nhitsVsSlat->Fill(j, nhits_Hodo[j]);
    	h1_BBHodo_EdepTotVsSlat->Fill(j, BBHodo_Edep_slats_400ns[j]);
    	h1_BBHodo_EdepTotVsSlat_log->Fill(j, BBHodo_Edep_slats_400ns[j]);
      }

      for(int j = 0; j<52; j++){
    	h1_BBPS_nhitsVsChan->Fill(j, nhits_PS[j]);
      }
      
      for(int j = 0; j<189; j++){
    	h1_BBSH_nhitsVsChan->Fill(j, nhits_SH[j]);
      }
      
      for(int j = 0; j<510; j++){
    	h1_GRINCH_nhitsVsChan->Fill(j, nhits_GRINCH[j]);
      }
      memset(nhits_HCal, 0, 288*sizeof(int));
      memset(nhits_Hodo, 0, 90*sizeof(int));
      memset(nhits_PS, 0, 52*sizeof(int));
      memset(nhits_SH, 0, 189*sizeof(int));
      memset(nhits_GRINCH, 0, 510*sizeof(int));
      
      memset(HCal_Edep_blocks_400ns, 0, 288*sizeof(double));
      memset(BBHodo_Edep_slats_400ns, 0, 90*sizeof(double));
      //}
      //if(FileCounter%15==0){
      for(int j = 0; j<5; j++){
	h1_BBGEM_nhits_[j]->Fill(nhits_GEM[j]);
      }
      memset(nhits_GEM, 0, 5*sizeof(int));
    }
    
    C1->Delete();
    T1->~gmn_deftree();
    //T1->~gmn_nope_tree();
    
  }//end while (( chEl= )) 
  //h1_Ntries->Fill(0.5, N1_tot);
  
  
  //h1_HCal_Edep_rates->Scale(I_exp/1.602e-19/N1_tot);

  C->Delete();
  fout->Write();
  //return;
}
