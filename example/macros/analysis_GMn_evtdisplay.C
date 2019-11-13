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

void chan2pos_GC(int chan, double& xpos, double& ypos){
  int scol = chan%17;
  int srow = (chan-scol)/17;
  int row = srow*2+scol/9;
  int col = scol%9;
  //cout << chan << " " << srow << " " << row << " " << scol << " " << col << endl;
  xpos = 12.0-3.0*col;
  if(row%2==1)xpos = 10.5-3.0*col;
  ypos = 88.5-3.0*row;
  
}

void chan2pos_PS(int chan, double& xpos, double& ypos){
  int row = chan%27;
  int col = (chan-row)/27;
  xpos = 18.5-37.0*col;
  ypos = 110.5-8.5*row;
}

void chan2pos_TH(int chan, double& xpos, double& ypos){
  int col = chan%2;
  int row = (chan-col)/2;
  xpos = -60+120.0*col;
  ypos = 111.25-2.5*row;
}

void chan2pos_SH(int chan, double& xpos, double& ypos){
  int row = chan%27;
  int col = (chan-row)/27;
  xpos = 25.5-8.5*col;
  ypos = 110.5-8.5*row;
}

void EvtDisplay_GMn(const char *inputfilename, 
		    const int Evt2display, const char* filesave = "")
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(55);
  // load the files that we want to analyze
  cout << "reading input files" << endl;
  
  ifstream inputfile(inputfilename);
  //TFile *fout = new TFile( root_outfile, "RECREATE" );
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
  //TH1D* h1_grinchMCtime = new TH1D("h1_grinchMCtime", "Grinch;MC_time (ns)", 200, -100, 100);
  //TH1D* h1_hodoMCtime = new TH1D("h1_hodoMCtime", "Hodoscope;MC_time (ns)", 200, -100, 100);
  //TH1D* h1_cdetMCtime = new TH1D("h1_cdetMCtime", "CDet;MC_time (ns)", 200, -100, 100);
  TH2D* h1_GRINCH = new TH2D("h1_GRINCH", "GRINCH PMT matrix;x (cm);y (cm)", 19, -14.25, 14.25, 122, -91.5, 91.5);
  TH2D* h1_PS = new TH2D("h1_PS", "Preshower;x (cm);y (cm)", 2, -37.0, +37.0, 27, -114.75, +114.75);
  TH2D* h1_TH = new TH2D("h1_TH", "Hodoscope;x (cm);y (cm)", 3, -90.0, +90.0, 90, -112.5, +112.5);
  TH2D* h1_SH = new TH2D("h1_SH", "Shower;x (cm);y (cm)", 7, -29.75, +29.75, 27, -114.75, +114.75);

  TH2D* h1_GRINCH_time = new TH2D("h1_GRINCH_time", "GRINCH;chan;time (ns)", 510, 0, 510, 200, -100, 100);
  TH2D* h1_TH_time = new TH2D("h1_TH_time", "Hodoscope;chan;time (ns)", 180, 0, 180, 200, -100, 100);
  
  double xpos, ypos;
  double tmean;
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
      if(T1->EvtID!=Evt2display)continue;
      
      if(T1->bb_grinch_Nsimhits){
	for(int i = 0; i<T1->bb_grinch_Nsimhits; i++){
	  chan2pos_GC(T1->bb_grinch_simhit_chan->at(i), xpos, ypos);
	  h1_GRINCH->Fill(xpos, ypos, T1->bb_grinch_simhit_npe->at(i));
	  h1_GRINCH_time->Fill(T1->bb_grinch_simhit_chan->at(i), T1->bb_grinch_simhit_t_lead->at(i));
	  h1_GRINCH_time->Fill(T1->bb_grinch_simhit_chan->at(i), T1->bb_grinch_simhit_t_trail->at(i));
	}
      }
      
      if(T1->bb_hodo_Nsimhits){
	for(int i = 0; i<T1->bb_hodo_Nsimhits; i++){
	  chan2pos_TH(T1->bb_hodo_simhit_chan->at(i), xpos, ypos);
	  h1_TH->Fill(xpos, ypos, T1->bb_hodo_simhit_npe->at(i));
	  h1_TH_time->Fill(T1->bb_hodo_simhit_chan->at(i), T1->bb_hodo_simhit_t_lead->at(i));
	  h1_TH_time->Fill(T1->bb_hodo_simhit_chan->at(i), T1->bb_hodo_simhit_t_trail->at(i));
	}
      }
      
      if(T1->bb_ps_Nsimhits){
	for(int i = 0; i<T1->bb_ps_Nsimhits; i++){
	  chan2pos_PS(T1->bb_ps_simhit_chan->at(i), xpos, ypos);
	  h1_PS->Fill(xpos, ypos, T1->bb_ps_simhit_npe->at(i));
	}
      }
      
      if(T1->bb_sh_Nsimhits){
	for(int i = 0; i<T1->bb_sh_Nsimhits; i++){
	  chan2pos_SH(T1->bb_sh_simhit_chan->at(i), xpos, ypos);
	  h1_SH->Fill(xpos, ypos, T1->bb_sh_simhit_npe->at(i));
	}
      }
      
      
      
    }//end loop on events

  }//end loop on files
  
  TCanvas* C1 = new TCanvas("C1", "C1", 1000, 800);
  C1->Divide(4,1);
  C1->cd(1);
  h1_GRINCH->Draw("colz, text");
  C1->cd(2);
  h1_PS->Draw("colz, text");
  C1->cd(3);
  h1_TH->Draw("colz, text");
  C1->cd(4);
  h1_SH->Draw("colz, text");
  
  // h1_GRINCH->SetMinimum(-100);
  // h1_GRINCH->SetMaximum(100);
  h1_GRINCH->SetMarkerSize(2.5);
  h1_GRINCH->GetXaxis()->SetTitleSize(0.07);
  h1_GRINCH->GetXaxis()->SetTitleOffset(0.5);
  h1_GRINCH->GetXaxis()->CenterTitle();
  h1_GRINCH->GetXaxis()->SetLabelSize(0.05);
  h1_GRINCH->GetYaxis()->SetTitleSize(0.07);
  h1_GRINCH->GetYaxis()->SetTitleOffset(0.5);
  h1_GRINCH->GetYaxis()->CenterTitle();
  h1_GRINCH->GetYaxis()->SetLabelSize(0.05);
  
  // h1_PS->SetMinimum(0);
  // h1_PS->SetMaximum(2000);
  h1_PS->SetMarkerSize(2.5);
  h1_PS->GetXaxis()->SetTitleSize(0.07);
  h1_PS->GetXaxis()->SetTitleOffset(0.5);
  h1_PS->GetXaxis()->CenterTitle();
  h1_PS->GetXaxis()->SetLabelSize(0.05);
  h1_PS->GetYaxis()->SetTitleSize(0.07);
  h1_PS->GetYaxis()->SetTitleOffset(0.5);
  h1_PS->GetYaxis()->CenterTitle();
  h1_PS->GetYaxis()->SetLabelSize(0.05);
  
  // h1_TH->SetMinimum(-100);
  // h1_TH->SetMaximum(100);
  h1_TH->SetMarkerSize(2.5);
  h1_TH->GetXaxis()->SetTitleSize(0.07);
  h1_TH->GetXaxis()->SetTitleOffset(0.5);
  h1_TH->GetXaxis()->CenterTitle();
  h1_TH->GetXaxis()->SetLabelSize(0.05);
  h1_TH->GetYaxis()->SetTitleSize(0.07);
  h1_TH->GetYaxis()->SetTitleOffset(0.5);
  h1_TH->GetYaxis()->CenterTitle();
  h1_TH->GetYaxis()->SetLabelSize(0.05);
  
  // h1_SH->SetMinimum(0);
  // h1_SH->SetMaximum(2000);
  h1_SH->SetMarkerSize(2.5);
  h1_SH->GetXaxis()->SetTitleSize(0.07);
  h1_SH->GetXaxis()->SetTitleOffset(0.5);
  h1_SH->GetXaxis()->CenterTitle();
  h1_SH->GetXaxis()->SetLabelSize(0.05);
  h1_SH->GetYaxis()->SetTitleSize(0.07);
  h1_SH->GetYaxis()->SetTitleOffset(0.5);
  h1_SH->GetYaxis()->CenterTitle();
  h1_SH->GetYaxis()->SetLabelSize(0.05);
  
  TCanvas* C2 = new TCanvas("C2", "C2", 1000, 800);
  C2->Divide(2,1);
  C2->cd(1);
  h1_GRINCH_time->Draw("colz");
  C2->cd(2);
  h1_TH_time->Draw("colz");
  
  h1_GRINCH_time->GetXaxis()->SetTitleSize(0.07);
  h1_GRINCH_time->GetXaxis()->SetTitleOffset(0.5);
  h1_GRINCH_time->GetXaxis()->CenterTitle();
  h1_GRINCH_time->GetXaxis()->SetLabelSize(0.05);
  h1_GRINCH_time->GetYaxis()->SetTitleSize(0.07);
  h1_GRINCH_time->GetYaxis()->SetTitleOffset(0.5);
  h1_GRINCH_time->GetYaxis()->CenterTitle();
  h1_GRINCH_time->GetYaxis()->SetLabelSize(0.05);

  h1_TH_time->GetXaxis()->SetTitleSize(0.07);
  h1_TH_time->GetXaxis()->SetTitleOffset(0.5);
  h1_TH_time->GetXaxis()->CenterTitle();
  h1_TH_time->GetXaxis()->SetLabelSize(0.05);
  h1_TH_time->GetYaxis()->SetTitleSize(0.07);
  h1_TH_time->GetYaxis()->SetTitleOffset(0.5);
  h1_TH_time->GetYaxis()->CenterTitle();
  h1_TH_time->GetYaxis()->SetLabelSize(0.05);
  
  if(strcmp(filesave, "")!=0)C1->SaveAs(Form("%s_1.pdf", filesave));
  if(strcmp(filesave, "")!=0)C2->SaveAs(Form("%s_2.pdf", filesave));
  
  //fout->Write();
  
  //C->Delete();
  //histloader1.CloseAll();
}
