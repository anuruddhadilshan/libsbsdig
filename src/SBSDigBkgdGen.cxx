#include "g4sbs_types.h"
#include "SBSDigBkgdGen.h"

#include "TH2D.h"
#include "TF1.h"
#include "TMath.h"

using namespace std;

SBSDigBkgdGen::SBSDigBkgdGen()
{
  //Initialization of arrays and histograms
  NhitsBBGEMs = new Double_t[5];
  h_xhitBBGEMs = new TH1D*[5];
  h_yhitBBGEMs = new TH1D*[5];
  h_dxhitBBGEMs = new TH1D*[5];
  h_dyhitBBGEMs = new TH1D*[5];
  NhitsHCal = new Double_t[288];
  NhitsBBPS = new Double_t[52];
  NhitsBBSH = new Double_t[189];
  NhitsBBHodo = new Double_t[90];
  NhitsPRPOLBS_SCINT = new Double_t[24];
  NhitsPRPOLFS_SCINT = new Double_t[90];
  NhitsACTIVEANA = new Double_t[90];
  P1hitGRINCH = new Double_t[510];
  P2hitsGRINCH = new Double_t[510];
  NhitsECal = new Double_t[1770];
  NhitsCDet = new Double_t[2352];
  NhitsFT = new Double_t[6];
  h_xhitFT = new TH1D*[5];
  h_yhitFT = new TH1D*[5];
  h_dxhitFT = new TH1D*[5];
  h_dyhitFT = new TH1D*[5];
  NhitsFPP1 = new Double_t[5];
  h_xhitFPP1 = new TH1D*[5];
  h_yhitFPP1 = new TH1D*[5];
  h_dxhitFPP1 = new TH1D*[5];
  h_dyhitFPP1 = new TH1D*[5];
  NhitsFPP2 = new Double_t[5];
  h_xhitFPP2 = new TH1D*[5];
  h_yhitFPP2 = new TH1D*[5];
  h_dxhitFPP2 = new TH1D*[5];
  h_dyhitFPP2 = new TH1D*[5];

  NhitsSBSGEMs = new Double_t[8];
  h_xhitSBSGEMs = new TH1D*[8];
  h_yhitSBSGEMs = new TH1D*[8];
  h_dxhitSBSGEMs = new TH1D*[8];
  h_dyhitSBSGEMs = new TH1D*[8];
  
  NhitsCEPOL_GEMFRONTs = new Double_t[4];
  h_xhitCEPOL_GEMFRONTs = new TH1D*[4];
  h_yhitCEPOL_GEMFRONTs = new TH1D*[4];
  h_dxhitCEPOL_GEMFRONTs = new TH1D*[4];
  h_dyhitCEPOL_GEMFRONTs = new TH1D*[4];

  NhitsCEPOL_GEMREARs = new Double_t[4];
  h_xhitCEPOL_GEMREARs = new TH1D*[4];
  h_yhitCEPOL_GEMREARs = new TH1D*[4];
  h_dxhitCEPOL_GEMREARs = new TH1D*[4];
  h_dyhitCEPOL_GEMREARs = new TH1D*[4];
  
  NhitsPRPOLBS_GEMs = new Double_t[4];
  h_xhitPRPOLBS_GEMs = new TH1D*[4];
  h_yhitPRPOLBS_GEMs = new TH1D*[4];
  h_dxhitPRPOLBS_GEMs = new TH1D*[4];
  h_dyhitPRPOLBS_GEMs = new TH1D*[4];
 
  NhitsPRPOLFS_GEMs = new Double_t[4];
  h_xhitPRPOLFS_GEMs = new TH1D*[4];
  h_yhitPRPOLFS_GEMs = new TH1D*[4];
  h_dxhitPRPOLFS_GEMs = new TH1D*[4];
  h_dyhitPRPOLFS_GEMs = new TH1D*[4];
}

SBSDigBkgdGen::SBSDigBkgdGen(TFile* f_bkgd, std::vector<TString> det_list, double timewindow, bool pmtbkgddig)
{
  //Initialization of arrays and histograms
  fTimeWindow = timewindow;
  NhitsBBGEMs = new Double_t[5];
  h_xhitBBGEMs = new TH1D*[5];
  h_yhitBBGEMs = new TH1D*[5];
  h_dxhitBBGEMs = new TH1D*[5];
  h_dyhitBBGEMs = new TH1D*[5];
  NhitsHCal = new Double_t[288];
  NhitsBBPS = new Double_t[52];
  NhitsBBSH = new Double_t[189];
  NhitsBBHodo = new Double_t[90];
  NhitsPRPOLBS_SCINT = new Double_t[24];
  NhitsPRPOLFS_SCINT = new Double_t[90];
  NhitsACTIVEANA = new Double_t[90];
  P1hitGRINCH = new Double_t[510];
  P2hitsGRINCH = new Double_t[510];
  NhitsECal = new Double_t[1770];
  NhitsCDet = new Double_t[2352];
  NhitsFT = new Double_t[6];
  h_xhitFT = new TH1D*[6];
  h_yhitFT = new TH1D*[6];
  h_dxhitFT = new TH1D*[6];
  h_dyhitFT = new TH1D*[6];
  NhitsFPP1 = new Double_t[5];
  h_xhitFPP1 = new TH1D*[5];
  h_yhitFPP1 = new TH1D*[5];
  h_dxhitFPP1 = new TH1D*[5];
  h_dyhitFPP1 = new TH1D*[5];
  NhitsFPP2 = new Double_t[5];
  h_xhitFPP2 = new TH1D*[5];
  h_yhitFPP2 = new TH1D*[5];
  h_dxhitFPP2 = new TH1D*[5];
  h_dyhitFPP2 = new TH1D*[5];

  NhitsSBSGEMs = new Double_t[8];
  h_xhitSBSGEMs = new TH1D*[8];
  h_yhitSBSGEMs = new TH1D*[8];
  h_dxhitSBSGEMs = new TH1D*[8];
  h_dyhitSBSGEMs = new TH1D*[8];
  
  NhitsCEPOL_GEMFRONTs = new Double_t[4];
  h_xhitCEPOL_GEMFRONTs = new TH1D*[4];
  h_yhitCEPOL_GEMFRONTs = new TH1D*[4];
  h_dxhitCEPOL_GEMFRONTs = new TH1D*[4];
  h_dyhitCEPOL_GEMFRONTs = new TH1D*[4];

  NhitsCEPOL_GEMREARs = new Double_t[4];
  h_xhitCEPOL_GEMREARs = new TH1D*[4];
  h_yhitCEPOL_GEMREARs = new TH1D*[4];
  h_dxhitCEPOL_GEMREARs = new TH1D*[4];
  h_dyhitCEPOL_GEMREARs = new TH1D*[4];
 
  NhitsPRPOLBS_GEMs = new Double_t[2];
  h_xhitPRPOLBS_GEMs = new TH1D*[2];
  h_yhitPRPOLBS_GEMs = new TH1D*[2];
  h_dxhitPRPOLBS_GEMs = new TH1D*[2];
  h_dyhitPRPOLBS_GEMs = new TH1D*[2];
 
  NhitsPRPOLFS_GEMs = new Double_t[2];
  h_xhitPRPOLFS_GEMs = new TH1D*[2];
  h_yhitPRPOLFS_GEMs = new TH1D*[2];
  h_dxhitPRPOLFS_GEMs = new TH1D*[2];
  h_dyhitPRPOLFS_GEMs = new TH1D*[2];
 
  fPMTBkgdDig = pmtbkgddig;
  
  Initialize(f_bkgd, det_list);
}

SBSDigBkgdGen::~SBSDigBkgdGen()
{
}

void SBSDigBkgdGen::Initialize(TFile* f_bkgd, std::vector<TString> det_list)
{
  double mu, sigma;
  
  // BB GEMs
  TH1D* h1_BBGEM_nhits_[5];
  TF1* f1_bbgemnhits_[5];
  TH2D* h1_BBGEM_yVsx_[5];
  TH2D* h1_BBGEM_dyVsdx_[5];
  TH1D* h1_BBGEM_Edep_[5];
  
  TH1D* h1_SBSGEM_nhits_[8];
  TF1* f1_sbsgemnhits_[8];
  TH2D* h1_SBSGEM_yVsx_[8];
  TH2D* h1_SBSGEM_dyVsdx_[8];
  TH1D* h1_SBSGEM_Edep_[8];
  
  //CEPOL_GEMFRONT
  TH1D* h1_CEPOL_GEMFRONT_nhits_[4];
  TF1* f1_CEPOL_GEMFRONTnhits_[4];
  TH2D* h1_CEPOL_GEMFRONT_yVsx_[4];
  TH2D* h1_CEPOL_GEMFRONT_dyVsdx_[4];
  TH1D* h1_CEPOL_GEMFRONT_Edep_[4];
  
 //CEPOL_GEMREAR
  TH1D* h1_CEPOL_GEMREAR_nhits_[4];
  TF1* f1_CEPOL_GEMREARnhits_[4];
  TH2D* h1_CEPOL_GEMREAR_yVsx_[4];
  TH2D* h1_CEPOL_GEMREAR_dyVsdx_[4];
  TH1D* h1_CEPOL_GEMREAR_Edep_[4];
  
 //PRPOLBS_GEM
  TH1D* h1_PRPOLBS_GEM_nhits_[2];
  TF1* f1_PRPOLBS_GEMnhits_[2];
  TH2D* h1_PRPOLBS_GEM_yVsx_[2];
  TH2D* h1_PRPOLBS_GEM_dyVsdx_[2];
  TH1D* h1_PRPOLBS_GEM_Edep_[2];
  
  //PRPOLFS_GEM
  TH1D* h1_PRPOLFS_GEM_nhits_[2];
  TF1* f1_PRPOLFS_GEMnhits_[2];
  TH2D* h1_PRPOLFS_GEM_yVsx_[2];
  TH2D* h1_PRPOLFS_GEM_dyVsdx_[2];
  TH1D* h1_PRPOLFS_GEM_Edep_[2];
  
  /*
  //cross-check histos
  h_NhitsBBGEMs_XC = new TH1D*[5];
  h_EdephitBBGEMs_XC = new TH1D("h_EdephitBBGEMs_XC", "", 1000, 0.0, 0.1);
  h_xhitBBGEMs_XC = new TH1D*[5];
  h_yhitBBGEMs_XC = new TH1D*[5];
  h_modBBGEMs_XC = new TH1D("h_modBBGEMs_XC", "", 36, 0, 36);
  //TH1D* h_dxhitBBGEMs_XC[5];
  //TH1D* h_dyhitBBGEMs_XC[5];
  
  h_NhitsHCal_XC = new TH2D("h_NhitsHCal_XC", "", 288, 0, 288, 100, 0, 100);
  h_EdephitHCal_XC = new TH1D("h_EdephitHCal_XC", "", 100, 0.0+1.0e-3, 1.0+1.0e-3);
  h_zhitHCal_XC = new TH1D("h_zhitHCal_XC", "", 100, 0., 1.);
  
  h_NhitsBBPS_XC = new TH2D("h_NhitsBBPS_XC", "", 52, 0, 52, 150, 0, 150);
  h_EdephitBBPS_XC = new TH1D("h_EdephitBBPS_XC", "", 150, 0.0+1.0e-3, 1.5+1.0e-3);
  
  h_NhitsBBSH_XC = new TH2D("h_NhitsBBSH_XC", "", 189, 0, 189, 100, 0, 100);
  h_EdephitBBSH_XC = new TH1D("h_EdephitBBSH_XC", "", 150, 0.0+1.0e-3, 1.5+1.0e-3);
  
  h_NhitsPRPOLFS_SCINT_XC = new TH2D("h_NhitsBBHodo_XC", "", 90, 0, 90, 100, 0, 100);
  h_EdephitBBHodo_XC = new TH1D("h_EdephitBBHodo_XC", "", 250, 0., 0.5);
  h_xhitBBHodo_XC = new TH1D("h_xhitBBHodo_XC", "", 60, -0.3, 0.3);
  
  h_NhitsGRINCH_XC = new TH2D("h_NhitsGRINCH_XC", "", 510, 0, 510, 20, 0, 20);
  h_NpeGRINCH_XC = new TH1D("h_NpeGRINCH_XC", "", 100, 0, 100);
  */
  
  // Initialization of BBGEM histograms:
  // most histograms are 1D projections the 2D histograms stored in the input file.
  // for the energy deposit, the input file histograms (one per plane) are consolidated into one
  for(size_t k = 0; k<det_list.size(); k++){
    if(det_list[k]=="bbgem"){
      cout << "BB GEMs" << endl;
      for(int m = 0; m<5; m++){
	// fit of the hits mulitplicity distribution.
	h1_BBGEM_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_BBGEM_nhits_%d",m));
	f1_bbgemnhits_[m] = new TF1(Form("f1_bbgemnhits_%d", m), "gaus", 0., 400.);
	h1_BBGEM_nhits_[m]->Fit(f1_bbgemnhits_[m], "QRN");
	mu = f1_bbgemnhits_[m]->GetParameter(1);
	sigma = f1_bbgemnhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_bbgemnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_BBGEM_nhits_[m]->Fit(f1_bbgemnhits_[m], "QRN");
	}
	NhitsBBGEMs[m] = max(1.0, f1_bbgemnhits_[m]->GetParameter(1));
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_BBGEM_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_BBGEM_yVsx_%d",m));
	h_xhitBBGEMs[m] = h1_BBGEM_yVsx_[m]->ProjectionX(Form("h1_xhitBBGEM_%d",m));
	h_yhitBBGEMs[m] = h1_BBGEM_yVsx_[m]->ProjectionY(Form("h1_yhitBBGEM_%d",m));
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_BBGEM_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_BBGEM_dyVsdx_%d",m));
	h_dxhitBBGEMs[m] = h1_BBGEM_dyVsdx_[m]->ProjectionX(Form("h1_dxhitBBGEM_%d",m));
	h_dyhitBBGEMs[m] = h1_BBGEM_dyVsdx_[m]->ProjectionY(Form("h1_dyhitBBGEM_%d",m));
	
	h1_BBGEM_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_BBGEM_Edep_%d",m));
	//h1_BBGEM_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_BBGEM_Edep_log_%d",m));
	
	if(m==0){
	  h_EdephitBBGEMs = h1_BBGEM_Edep_[m];
	}else{
	  h_EdephitBBGEMs->Add(h1_BBGEM_Edep_[m]);
	}
	
	/*
	//cross-check histos
	h_NhitsBBGEMs_XC[m] = new TH1D(Form("h_NhitsBBGEMs_XC_%d", m), "", 250, 0, 1000);
	h_xhitBBGEMs_XC[m] = new TH1D(Form("h_xhitBBGEMs_XC_%d", m), "", 205, -1.025, 1.025);
	h_yhitBBGEMs_XC[m] = new TH1D(Form("h_yhitBBGEMs_XC_%d", m), "", 62, -0.31, 0.31);
	//h_dxhitBBGEMs_XC[m] = new TH1D(Form("h_dxhitBBGEMs_XC_%d", m), "", 100, 0.05, 0.05);
	//h_dyhitBBGEMs_XC[m] = new TH1D(Form("h_dyhitBBGEMs_XC_%d", m), "", 100, -0.05, 0.05);
	*/
	
	
	//cout << m << " " << NhitsBBGEMs[m] << endl;
      }
    }
  // Initialization of SBSGEM histograms:
  // most histograms are 1D projections the 2D histograms stored in the input file.
  // for the energy deposit, the input file histograms (one per plane) are consolidated into one
    if(det_list[k]=="sbsgem"){
      cout << "SBS GEMs" << endl;
      for(int m = 0; m<8; m++){
	// fit of the hits mulitplicity distribution.
	h1_SBSGEM_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_SBSGEM_nhits_%d",m));
	f1_sbsgemnhits_[m] = new TF1(Form("f1_sbsgemnhits_%d", m), "gaus", 0., 400.);
	h1_SBSGEM_nhits_[m]->Fit(f1_sbsgemnhits_[m], "QRN");
	mu = f1_sbsgemnhits_[m]->GetParameter(1);
	sigma = f1_sbsgemnhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_sbsgemnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_SBSGEM_nhits_[m]->Fit(f1_sbsgemnhits_[m], "QRN");
	}
	NhitsSBSGEMs[m] = max(1.0, f1_sbsgemnhits_[m]->GetParameter(1));
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_SBSGEM_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_SBSGEM_yVsx_%d",m));
	h_xhitSBSGEMs[m] = h1_SBSGEM_yVsx_[m]->ProjectionX(Form("h1_xhitSBSGEM_%d",m));
	h_yhitSBSGEMs[m] = h1_SBSGEM_yVsx_[m]->ProjectionY(Form("h1_yhitSBSGEM_%d",m));
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_SBSGEM_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_SBSGEM_dyVsdx_%d",m));
	h_dxhitSBSGEMs[m] = h1_SBSGEM_dyVsdx_[m]->ProjectionX(Form("h1_dxhitSBSGEM_%d",m));
	h_dyhitSBSGEMs[m] = h1_SBSGEM_dyVsdx_[m]->ProjectionY(Form("h1_dyhitSBSGEM_%d",m));
	
	h1_SBSGEM_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_SBSGEM_Edep_%d",m));
	//h1_SBSGEM_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_SBSGEM_Edep_log_%d",m));
	
	if(m==0){
	  h_EdephitSBSGEMs = h1_SBSGEM_Edep_[m];
	}else{
	  h_EdephitSBSGEMs->Add(h1_SBSGEM_Edep_[m]);
	}
	
	/*
	//cross-check histos
	h_NhitsSBSGEMs_XC[m] = new TH1D(Form("h_NhitsSBSGEMs_XC_%d", m), "", 250, 0, 1000);
	h_xhitSBSGEMs_XC[m] = new TH1D(Form("h_xhitSBSGEMs_XC_%d", m), "", 205, -1.025, 1.025);
	h_yhitSBSGEMs_XC[m] = new TH1D(Form("h_yhitSBSGEMs_XC_%d", m), "", 62, -0.31, 0.31);
	//h_dxhitSBSGEMs_XC[m] = new TH1D(Form("h_dxhitSBSGEMs_XC_%d", m), "", 100, 0.05, 0.05);
	//h_dyhitSBSGEMs_XC[m] = new TH1D(Form("h_dyhitSBSGEMs_XC_%d", m), "", 100, -0.05, 0.05);
	*/
	
	
	//cout << m << " " << NhitsSBSGEMs[m] << endl;
      }
    }
    
    if(det_list[k]=="cepol_front"){
      //CEPOL_GEMFRONT
      
      cout << "CEPOL_GEMFRONT GEMs" << endl;
      
      for(int m = 0; m<4; m++){
	// fit of the hits mulitplicity distribution.
	h1_CEPOL_GEMFRONT_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_CEPOL_GEMFRONT_nhits_%d",m));
	f1_CEPOL_GEMFRONTnhits_[m] = new TF1(Form("f1_CEPOL_GEMFRONTnhits_%d", m), "gaus", 0., 400.);
	cout << m << " " << f1_CEPOL_GEMFRONTnhits_[m] << " " << h1_CEPOL_GEMFRONT_nhits_[m] << endl;
	h1_CEPOL_GEMFRONT_nhits_[m]->Fit(f1_CEPOL_GEMFRONTnhits_[m], "QRN");
	mu = f1_CEPOL_GEMFRONTnhits_[m]->GetParameter(1);
	sigma = f1_CEPOL_GEMFRONTnhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_CEPOL_GEMFRONTnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_CEPOL_GEMFRONT_nhits_[m]->Fit(f1_CEPOL_GEMFRONTnhits_[m], "QRN");
	}
	NhitsCEPOL_GEMFRONTs[m] = max(1.0, f1_CEPOL_GEMFRONTnhits_[m]->GetParameter(1));
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_CEPOL_GEMFRONT_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_CEPOL_GEMFRONT_yVsx_%d",m));
	h_xhitCEPOL_GEMFRONTs[m] = h1_CEPOL_GEMFRONT_yVsx_[m]->ProjectionX(Form("h1_xhitCEPOL_GEMFRONT_%d",m));
	h_yhitCEPOL_GEMFRONTs[m] = h1_CEPOL_GEMFRONT_yVsx_[m]->ProjectionY(Form("h1_yhitCEPOL_GEMFRONT_%d",m));
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_CEPOL_GEMFRONT_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_CEPOL_GEMFRONT_dyVsdx_%d",m));
	h_dxhitCEPOL_GEMFRONTs[m] = h1_CEPOL_GEMFRONT_dyVsdx_[m]->ProjectionX(Form("h1_dxhitCEPOL_GEMFRONT_%d",m));
	h_dyhitCEPOL_GEMFRONTs[m] = h1_CEPOL_GEMFRONT_dyVsdx_[m]->ProjectionY(Form("h1_dyhitCEPOL_GEMFRONT_%d",m));
	
	h1_CEPOL_GEMFRONT_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_CEPOL_GEMFRONT_Edep_%d",m));
	//h1_CEPOL_GEMFRONT_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_CEPOL_GEMFRONT_Edep_log_%d",m));
	
	if(m==0){
	  h_EdephitCEPOL_GEMFRONTs = h1_CEPOL_GEMFRONT_Edep_[m];
	}else{
	  h_EdephitCEPOL_GEMFRONTs->Add(h1_CEPOL_GEMFRONT_Edep_[m]);
	}
	
	//cout << m << " " << NhitsBBGEMs[m] << endl;
      }
    }
    
    if(det_list[k]=="cepol_rear"){
      cout << "CEPOL_GEMREAR GEMs" << endl;
      
      //CEPOL_GEMREAR
      for(int m = 0; m<4; m++){
	// fit of the hits mulitplicity distribution.
	h1_CEPOL_GEMREAR_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_CEPOL_GEMREAR_nhits_%d",m));
	f1_CEPOL_GEMREARnhits_[m] = new TF1(Form("f1_CEPOL_GEMREARnhits_%d", m), "gaus", 0., 400.);
	h1_CEPOL_GEMREAR_nhits_[m]->Fit(f1_CEPOL_GEMREARnhits_[m], "QRN");
	mu = f1_CEPOL_GEMREARnhits_[m]->GetParameter(1);
	sigma = f1_CEPOL_GEMREARnhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_CEPOL_GEMREARnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_CEPOL_GEMREAR_nhits_[m]->Fit(f1_CEPOL_GEMREARnhits_[m], "QRN");
	}
	NhitsCEPOL_GEMREARs[m] = max(1.0, f1_CEPOL_GEMREARnhits_[m]->GetParameter(1));
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_CEPOL_GEMREAR_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_CEPOL_GEMREAR_yVsx_%d",m));
	h_xhitCEPOL_GEMREARs[m] = h1_CEPOL_GEMREAR_yVsx_[m]->ProjectionX(Form("h1_xhitCEPOL_GEMREAR_%d",m));
	h_yhitCEPOL_GEMREARs[m] = h1_CEPOL_GEMREAR_yVsx_[m]->ProjectionY(Form("h1_yhitCEPOL_GEMREAR_%d",m));
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_CEPOL_GEMREAR_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_CEPOL_GEMREAR_dyVsdx_%d",m));
	h_dxhitCEPOL_GEMREARs[m] = h1_CEPOL_GEMREAR_dyVsdx_[m]->ProjectionX(Form("h1_dxhitCEPOL_GEMREAR_%d",m));
	h_dyhitCEPOL_GEMREARs[m] = h1_CEPOL_GEMREAR_dyVsdx_[m]->ProjectionY(Form("h1_dyhitCEPOL_GEMREAR_%d",m));
	
	h1_CEPOL_GEMREAR_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_CEPOL_GEMREAR_Edep_%d",m));
	//h1_CEPOL_GEMREAR_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_CEPOL_GEMREAR_Edep_log_%d",m));
	
	if(m==0){
	  h_EdephitCEPOL_GEMREARs = h1_CEPOL_GEMREAR_Edep_[m];
	}else{
	  h_EdephitCEPOL_GEMREARs->Add(h1_CEPOL_GEMREAR_Edep_[m]);
	}
	
	//cout << m << " " << NhitsBBGEMs[m] << endl;
      }
    }
    
    if(det_list[k]=="prpolbs_gem"){
      cout << "PRPOLBS_GEM GEMs" << endl;
      
      //PRPOLBS_GEM
      for(int m = 0; m<2; m++){
	// fit of the hits mulitplicity distribution.
	h1_PRPOLBS_GEM_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_PRPOLBS_GEM_nhits_%d",m));
	f1_PRPOLBS_GEMnhits_[m] = new TF1(Form("f1_PRPOLBS_GEMnhits_%d", m), "gaus", 0., 400.);
	h1_PRPOLBS_GEM_nhits_[m]->Fit(f1_PRPOLBS_GEMnhits_[m], "QRN");
	mu = f1_PRPOLBS_GEMnhits_[m]->GetParameter(1);
	sigma = f1_PRPOLBS_GEMnhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_PRPOLBS_GEMnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_PRPOLBS_GEM_nhits_[m]->Fit(f1_PRPOLBS_GEMnhits_[m], "QRN");
	}
	NhitsPRPOLBS_GEMs[m] = max(1.0, f1_PRPOLBS_GEMnhits_[m]->GetParameter(1));
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_PRPOLBS_GEM_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_PRPOLBS_GEM_yVsx_%d",m));
	h_xhitPRPOLBS_GEMs[m] = h1_PRPOLBS_GEM_yVsx_[m]->ProjectionX(Form("h1_xhitPRPOLBS_GEM_%d",m));
	h_yhitPRPOLBS_GEMs[m] = h1_PRPOLBS_GEM_yVsx_[m]->ProjectionY(Form("h1_yhitPRPOLBS_GEM_%d",m));
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_PRPOLBS_GEM_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_PRPOLBS_GEM_dyVsdx_%d",m));
	h_dxhitPRPOLBS_GEMs[m] = h1_PRPOLBS_GEM_dyVsdx_[m]->ProjectionX(Form("h1_dxhitPRPOLBS_GEM_%d",m));
	h_dyhitPRPOLBS_GEMs[m] = h1_PRPOLBS_GEM_dyVsdx_[m]->ProjectionY(Form("h1_dyhitPRPOLBS_GEM_%d",m));
	
	h1_PRPOLBS_GEM_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_PRPOLBS_GEM_Edep_%d",m));
	//h1_PRPOLBS_GEM_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_PRPOLBS_GEM_Edep_log_%d",m));
	
	if(m==0){
	  h_EdephitPRPOLBS_GEMs = h1_PRPOLBS_GEM_Edep_[m];
	}else{
	  h_EdephitPRPOLBS_GEMs->Add(h1_PRPOLBS_GEM_Edep_[m]);
	}
	
	//cout << m << " " << NhitsBBGEMs[m] << endl;
      }
    }
    
    if(det_list[k]=="prpolfs_gem"){
      cout << "PRPOLFS_GEM GEMs" << endl;
      
      //PRPOLFS_GEM
      for(int m = 0; m<2; m++){
	// fit of the hits mulitplicity distribution.
	h1_PRPOLFS_GEM_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_PRPOLFS_GEM_nhits_%d",m));
	f1_PRPOLFS_GEMnhits_[m] = new TF1(Form("f1_PRPOLFS_GEMnhits_%d", m), "gaus", 0., 400.);
	h1_PRPOLFS_GEM_nhits_[m]->Fit(f1_PRPOLFS_GEMnhits_[m], "QRN");
	mu = f1_PRPOLFS_GEMnhits_[m]->GetParameter(1);
	sigma = f1_PRPOLFS_GEMnhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_PRPOLFS_GEMnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_PRPOLFS_GEM_nhits_[m]->Fit(f1_PRPOLFS_GEMnhits_[m], "QRN");
	}
	NhitsPRPOLFS_GEMs[m] = max(1.0, f1_PRPOLFS_GEMnhits_[m]->GetParameter(1));
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_PRPOLFS_GEM_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_PRPOLFS_GEM_yVsx_%d",m));
	h_xhitPRPOLFS_GEMs[m] = h1_PRPOLFS_GEM_yVsx_[m]->ProjectionX(Form("h1_xhitPRPOLFS_GEM_%d",m));
	h_yhitPRPOLFS_GEMs[m] = h1_PRPOLFS_GEM_yVsx_[m]->ProjectionY(Form("h1_yhitPRPOLFS_GEM_%d",m));
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_PRPOLFS_GEM_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_PRPOLFS_GEM_dyVsdx_%d",m));
	h_dxhitPRPOLFS_GEMs[m] = h1_PRPOLFS_GEM_dyVsdx_[m]->ProjectionX(Form("h1_dxhitPRPOLFS_GEM_%d",m));
	h_dyhitPRPOLFS_GEMs[m] = h1_PRPOLFS_GEM_dyVsdx_[m]->ProjectionY(Form("h1_dyhitPRPOLFS_GEM_%d",m));
	
	h1_PRPOLFS_GEM_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_PRPOLFS_GEM_Edep_%d",m));
	//h1_PRPOLFS_GEM_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_PRPOLFS_GEM_Edep_log_%d",m));
	
	if(m==0){
	  h_EdephitPRPOLFS_GEMs = h1_PRPOLFS_GEM_Edep_[m];
	}else{
	  h_EdephitPRPOLFS_GEMs->Add(h1_PRPOLFS_GEM_Edep_[m]);
	}
	
	//cout << m << " " << NhitsBBGEMs[m] << endl;
      }
    }
    
    if(det_list[k]=="hcal"){
      // Initialization of HCal histograms:
      // most histograms are 1D projections the 2D histograms stored in the input file.
      cout << "HCal" << endl;
      TH2D *h1_HCal_nhitsVsChan = (TH2D*)f_bkgd->Get("h1_HCal_nhitsVsChan");
      TH1D* h1_HCal_nhits_[288];
      TF1* f1_hcalnhits_[288];
      TH2D *h1_HCal_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_HCal_EdepHitVsChan");
      //TH2D *h1_HCal_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_HCal_EdepHitVsChan_log");
      TH2D *h1_HCal_zHitVsChan = (TH2D*)f_bkgd->Get("h1_HCal_zHitVsChan");
      
      for(int m = 0; m<288; m++){
	// fit of the hits mulitplicity distribution.
	h1_HCal_nhits_[m] = h1_HCal_nhitsVsChan->ProjectionY(Form("h1_HCal_nhits_%d", m), m+1, m+1);
	f1_hcalnhits_[m] = new TF1(Form("f1_hcalnhits_%d", m), "gaus", 0, 50);
	h1_HCal_nhits_[m]->Fit(f1_hcalnhits_[m], "QR0");
	mu = f1_hcalnhits_[m]->GetParameter(1);
	sigma = f1_hcalnhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_hcalnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_HCal_nhits_[m]->Fit(f1_hcalnhits_[m], "QR0");
	}
	NhitsHCal[m] = max(1.0, f1_hcalnhits_[m]->GetParameter(1));
	//cout << m << " " << NhitsHCal[m] << endl;
      }
      // copy of 2D position histograms, then projection to obtain 1D histograms
      h_EdephitHCal = h1_HCal_EdepHitVsChan->ProjectionY("h_EdephitHCal");
      h_zhitHCal = h1_HCal_zHitVsChan->ProjectionY("h_zhitHCal");
    }
    
    if(det_list[k]=="bbps"){
      //PS
      cout << "PS" << endl;
      TH2D *h1_BBPS_nhitsVsChan = (TH2D*)f_bkgd->Get("h1_BBPS_nhitsVsChan");
      TH1D* h1_BBPS_nhits_[52];
      TF1* f1_bbpsnhits_[52];
      TH2D *h1_BBPS_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_BBPS_EdepHitVsChan");
      //TH2D *h1_BBPS_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_BBPS_EdepHitVsChan_log");
      
      for(int m = 0; m<52; m++){
	h1_BBPS_nhits_[m] = h1_BBPS_nhitsVsChan->ProjectionY(Form("h1_BBPS_nhits_%d", m), m+1, m+1);
	f1_bbpsnhits_[m] = new TF1(Form("f1_bbpsnhits_%d", m), "gaus", 0, 150);
	h1_BBPS_nhits_[m]->Fit(f1_bbpsnhits_[m], "QR0");
	mu = f1_bbpsnhits_[m]->GetParameter(1);
	sigma = f1_bbpsnhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_bbpsnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_BBPS_nhits_[m]->Fit(f1_bbpsnhits_[m], "QR0");
	}
	
	NhitsBBPS[m] = max(1.0, f1_bbpsnhits_[m]->GetParameter(1));
	//cout << m << " " << NhitsBBPS[m] << endl;
      }
      h_EdephitBBPS = h1_BBPS_EdepHitVsChan->ProjectionY("h_EdephitBBPS");
    }
    
    if(det_list[k]=="bbsh"){
      //SH
      cout << "SH" << endl;
      TH2D *h1_BBSH_nhitsVsChan = (TH2D*)f_bkgd->Get("h1_BBSH_nhitsVsChan");
      TH1D* h1_BBSH_nhits_[189];
      TF1* f1_bbshnhits_[189];
      TH2D *h1_BBSH_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_BBSH_EdepHitVsChan");
      //TH2D *h1_BBSH_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_BBSH_EdepHitVsChan_log");
      
      for(int m = 0; m<189; m++){
	h1_BBSH_nhits_[m] = h1_BBSH_nhitsVsChan->ProjectionY(Form("h1_BBSH_nhits_%d", m), m+1, m+1);
	f1_bbshnhits_[m] = new TF1(Form("f1_bbshnhits_%d", m), "gaus", 0, 150);
	h1_BBSH_nhits_[m]->Fit(f1_bbshnhits_[m], "QR0");
	mu = f1_bbshnhits_[m]->GetParameter(1);
	sigma = f1_bbshnhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_bbshnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_BBSH_nhits_[m]->Fit(f1_bbshnhits_[m], "QR0");
	}
	NhitsBBSH[m] = max(1.0, f1_bbshnhits_[m]->GetParameter(1));
	//cout << m << " " << NhitsBBSH[m] << endl;
      }
      h_EdephitBBSH = h1_BBSH_EdepHitVsChan->ProjectionY("h_EdephitBBSH");
    }
    
    if(det_list[k]=="bbhodo"){
      //BB Hodo
      cout << "Hodo" << endl;
      TH2D *h1_BBHodo_nhitsVsSlat = (TH2D*)f_bkgd->Get("h1_BBHodo_nhitsVsSlat");
      TH1D* h1_BBHodo_nhits_[90];
      TF1* f1_bbhodonhits_[90];
      
      TH2D *h1_BBHodo_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_BBHodo_EdepHitVsSlat");
      //TH2D *h1_BBHodo_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_BBHodo_EdepHitVsSlat_log");
      TH2D *h1_BBHodo_xhitVsSlat = (TH2D*)f_bkgd->Get("h1_BBHodo_xhitVsSlat");
      
      for(int m = 0; m<90; m++){
	h1_BBHodo_nhits_[m] = h1_BBHodo_nhitsVsSlat->ProjectionY(Form("h1_BBHodo_nhits_%d", m), m+1, m+1);
	f1_bbhodonhits_[m] = new TF1(Form("f1_bbhodonhits_%d", m), "gaus", 0, 100);
	h1_BBHodo_nhits_[m]->Fit(f1_bbhodonhits_[m], "QR0");
	mu = f1_bbhodonhits_[m]->GetParameter(1);
	sigma = f1_bbhodonhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_bbhodonhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_BBHodo_nhits_[m]->Fit(f1_bbhodonhits_[m], "QR0");
	}
	NhitsBBHodo[m] = max(1.0, f1_bbhodonhits_[m]->GetParameter(1));
	//cout << m << " " << NhitsBBHodo[m] << endl;
      }
      
      h_EdephitBBHodo = h1_BBHodo_EdepHitVsSlat->ProjectionY("h_EdephitBBHodo");
      h_xhitBBHodo = h1_BBHodo_xhitVsSlat->ProjectionY("h_xhitBBHodo");
    }
    
    if(det_list[k]=="prpolscint_bs"){
      //GEN-rp PRPOLBS_SCINT
      cout << "PRPOLBS_SCINT" << endl;
      TH2D *h1_PRPOLBS_SCINT_nhitsVsSlat = (TH2D*)f_bkgd->Get("h1_PRPOLBS_SCINT_nhitsVsSlat");
      TH1D* h1_PRPOLBS_SCINT_nhits_[24];
      TF1* f1_PRPOLBS_SCINTnhits_[24];
      
      TH2D *h1_PRPOLBS_SCINT_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_PRPOLBS_SCINT_EdepHitVsSlat");
      //TH2D *h1_PRPOLBS_SCINT_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_PRPOLBS_SCINT_EdepHitVsSlat_log");
      TH2D *h1_PRPOLBS_SCINT_xhitVsSlat = (TH2D*)f_bkgd->Get("h1_PRPOLBS_SCINT_xhitVsSlat");
      
      for(int m = 0; m<24; m++){
	h1_PRPOLBS_SCINT_nhits_[m] = h1_PRPOLBS_SCINT_nhitsVsSlat->ProjectionY(Form("h1_PRPOLBS_SCINT_nhits_%d", m), m+1, m+1);
	f1_PRPOLBS_SCINTnhits_[m] = new TF1(Form("f1_PRPOLBS_SCINTnhits_%d", m), "gaus", 0, 100);
	h1_PRPOLBS_SCINT_nhits_[m]->Fit(f1_PRPOLBS_SCINTnhits_[m], "QR0");
	mu = f1_PRPOLBS_SCINTnhits_[m]->GetParameter(1);
	sigma = f1_PRPOLBS_SCINTnhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_PRPOLBS_SCINTnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_PRPOLBS_SCINT_nhits_[m]->Fit(f1_PRPOLBS_SCINTnhits_[m], "QR0");
	}
	NhitsPRPOLBS_SCINT[m] = max(1.0, f1_PRPOLBS_SCINTnhits_[m]->GetParameter(1));
	//cout << m << " " << NhitsPRPOLBS_SCINT[m] << endl;
      }
      
      h_EdephitPRPOLBS_SCINT = h1_PRPOLBS_SCINT_EdepHitVsSlat->ProjectionY("h_EdephitPRPOLBS_SCINT");
      h_xhitPRPOLBS_SCINT = h1_PRPOLBS_SCINT_xhitVsSlat->ProjectionY("h_xhitPRPOLBS_SCINT");
    }
    
    if(det_list[k]=="prpolscint_fs"){
      //GEN-rp PRPOLFS_SCINT
      cout << "PRPOLFS_SCINT" << endl;
      TH2D *h1_PRPOLFS_SCINT_nhitsVsSlat = (TH2D*)f_bkgd->Get("h1_PRPOLFS_SCINT_nhitsVsSlat");
      TH1D* h1_PRPOLFS_SCINT_nhits_[90];
      TF1* f1_PRPOLFS_SCINTnhits_[90];
      
      TH2D *h1_PRPOLFS_SCINT_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_PRPOLFS_SCINT_EdepHitVsSlat");
      //TH2D *h1_PRPOLFS_SCINT_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_PRPOLFS_SCINT_EdepHitVsSlat_log");
      TH2D *h1_PRPOLFS_SCINT_xhitVsSlat = (TH2D*)f_bkgd->Get("h1_PRPOLFS_SCINT_xhitVsSlat");
      
      for(int m = 0; m<90; m++){
	h1_PRPOLFS_SCINT_nhits_[m] = h1_PRPOLFS_SCINT_nhitsVsSlat->ProjectionY(Form("h1_PRPOLFS_SCINT_nhits_%d", m), m+1, m+1);
	f1_PRPOLFS_SCINTnhits_[m] = new TF1(Form("f1_PRPOLFS_SCINTnhits_%d", m), "gaus", 0, 100);
	h1_PRPOLFS_SCINT_nhits_[m]->Fit(f1_PRPOLFS_SCINTnhits_[m], "QR0");
	mu = f1_PRPOLFS_SCINTnhits_[m]->GetParameter(1);
	sigma = f1_PRPOLFS_SCINTnhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_PRPOLFS_SCINTnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_PRPOLFS_SCINT_nhits_[m]->Fit(f1_PRPOLFS_SCINTnhits_[m], "QR0");
	}
	NhitsPRPOLFS_SCINT[m] = max(1.0, f1_PRPOLFS_SCINTnhits_[m]->GetParameter(1));
	//cout << m << " " << NhitsPRPOLFS_SCINT[m] << endl;
      }
      
      h_EdephitPRPOLFS_SCINT = h1_PRPOLFS_SCINT_EdepHitVsSlat->ProjectionY("h_EdephitPRPOLFS_SCINT");
      h_xhitPRPOLFS_SCINT = h1_PRPOLFS_SCINT_xhitVsSlat->ProjectionY("h_xhitPRPOLFS_SCINT");
    }

    if(det_list[k]=="activeana"){
      //GEN-rp ACTIVEANA
      cout << "ACTIVEANA" << endl;
      TH2D *h1_ACTIVEANA_nhitsVsSlat = (TH2D*)f_bkgd->Get("h1_ACTIVEANA_nhitsVsSlat");
      TH1D* h1_ACTIVEANA_nhits_[90];
      TF1* f1_ACTIVEANAnhits_[90];
      
      TH2D *h1_ACTIVEANA_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_ACTIVEANA_EdepHitVsSlat");
      //TH2D *h1_ACTIVEANA_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_ACTIVEANA_EdepHitVsSlat_log");
      TH2D *h1_ACTIVEANA_xhitVsSlat = (TH2D*)f_bkgd->Get("h1_ACTIVEANA_xhitVsSlat");
      
      for(int m = 0; m<90; m++){
	h1_ACTIVEANA_nhits_[m] = h1_ACTIVEANA_nhitsVsSlat->ProjectionY(Form("h1_ACTIVEANA_nhits_%d", m), m+1, m+1);
	f1_ACTIVEANAnhits_[m] = new TF1(Form("f1_ACTIVEANAnhits_%d", m), "gaus", 0, 100);
	h1_ACTIVEANA_nhits_[m]->Fit(f1_ACTIVEANAnhits_[m], "QR0");
	mu = f1_ACTIVEANAnhits_[m]->GetParameter(1);
	sigma = f1_ACTIVEANAnhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_ACTIVEANAnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_ACTIVEANA_nhits_[m]->Fit(f1_ACTIVEANAnhits_[m], "QR0");
	}
	NhitsACTIVEANA[m] = max(1.0, f1_ACTIVEANAnhits_[m]->GetParameter(1));
	//cout << m << " " << NhitsACTIVEANA[m] << endl;
      }
      
      h_EdephitACTIVEANA = h1_ACTIVEANA_EdepHitVsSlat->ProjectionY("h_EdephitACTIVEANA");
      h_xhitACTIVEANA = h1_ACTIVEANA_xhitVsSlat->ProjectionY("h_xhitACTIVEANA");
    }
    
    if(det_list[k]=="grinch"){
      //GRINCH
      cout << "GRINCH" << endl;
      TH2D *h1_GRINCH_nhitsVsChan = (TH2D*)f_bkgd->Get("h1_GRINCH_nhitsVsChan");
      // TH1D* h1_GRINCH_Chan_1hit = h1_GRINCH_nhitsVsChan->ProjectionX("h1_GRINCH_Chan_1hit", 2, 2);
      // TH1D* h1_GRINCH_Chan_2hits = h1_GRINCH_nhitsVsChan->ProjectionX("h1_GRINCH_Chan_2hit", 3, -1);
      TH2D *h1_GRINCH_NpeVsChan = (TH2D*)f_bkgd->Get("h1_GRINCH_NpeVsChan");
      
      for(int m = 1; m<=510; m++){
	P1hitGRINCH[m] = h1_GRINCH_nhitsVsChan->Integral(m, m, 2, 2)/h1_GRINCH_nhitsVsChan->Integral(m, m, 0, -1);
	P2hitsGRINCH[m] = h1_GRINCH_nhitsVsChan->Integral(m, m, 3, -1)/h1_GRINCH_nhitsVsChan->Integral(m, m, 0, -1);
      }
      h_NpeGRINCH = h1_GRINCH_NpeVsChan->ProjectionY("h1_GRINCH_Npe");
    }
    
    
    if(det_list[k]=="ft"){
      TH1D* h1_FT_nhits_[6];
      TF1* f1_ftnhits_[6];
      TH2D* h1_FT_yVsx_[6];
      TH2D* h1_FT_dyVsdx_[6];
      TH1D* h1_FT_Edep_[6];
      
      cout << "FT" << endl;
      
      for(int m = 0; m<6; m++){
	//Nhits
	h1_FT_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_FT_nhits_%d",m));
	f1_ftnhits_[m] = new TF1(Form("f1_ftnhits_%d", m), "gaus", 0., 400.);
	h1_FT_nhits_[m]->Fit(f1_ftnhits_[m], "QRN");
	mu = f1_ftnhits_[m]->GetParameter(1);
	sigma = f1_ftnhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_ftnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_FT_nhits_[m]->Fit(f1_ftnhits_[m], "QRN");
	}
	NhitsFT[m] = max(1.0, f1_ftnhits_[m]->GetParameter(1));
	
	h1_FT_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_FT_yVsx_%d",m));
	h_xhitFT[m] = h1_FT_yVsx_[m]->ProjectionX(Form("h1_xhitFT_%d",m));
	h_yhitFT[m] = h1_FT_yVsx_[m]->ProjectionY(Form("h1_yhitFT_%d",m));
	
	h1_FT_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_FT_dyVsdx_%d",m));
	h_dxhitFT[m] = h1_FT_dyVsdx_[m]->ProjectionX(Form("h1_dxhitFT_%d",m));
	h_dyhitFT[m] = h1_FT_dyVsdx_[m]->ProjectionY(Form("h1_dyhitFT_%d",m));
	
	h1_FT_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_FT_Edep_%d",m));
	//h1_FT_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_FT_Edep_log_%d",m));
	
	if(m==0){
	  h_EdephitFT = h1_FT_Edep_[m];
	}else{
	  h_EdephitFT->Add(h1_FT_Edep_[m]);
	}
	//cout << m << " " << NhitsFT[m] << endl;
      }
    }
  
    if(det_list[k]=="fpp1"){
      TH1D* h1_FPP1_nhits_[5];
      TF1* f1_fpp1nhits_[5];
      TH2D* h1_FPP1_yVsx_[5];
      TH2D* h1_FPP1_dyVsdx_[5];
      TH1D* h1_FPP1_Edep_[5];
      
      cout << "FPP1" << endl;
      
      for(int m = 0; m<6; m++){
	//Nhits
	h1_FPP1_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_FPP1_nhits_%d",m));
	f1_fpp1nhits_[m] = new TF1(Form("f1_fpp1nhits_%d", m), "gaus", 0., 400.);
	h1_FPP1_nhits_[m]->Fit(f1_fpp1nhits_[m], "QRN");
	mu = f1_fpp1nhits_[m]->GetParameter(1);
	sigma = f1_fpp1nhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_fpp1nhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_FPP1_nhits_[m]->Fit(f1_fpp1nhits_[m], "QRN");
	}
	NhitsFPP1[m] = max(1.0, f1_fpp1nhits_[m]->GetParameter(1));
	
	h1_FPP1_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_FPP1_yVsx_%d",m));
	h_xhitFPP1[m] = h1_FPP1_yVsx_[m]->ProjectionX(Form("h1_xhitFPP1_%d",m));
	h_yhitFPP1[m] = h1_FPP1_yVsx_[m]->ProjectionY(Form("h1_yhitFPP1_%d",m));
	
	h1_FPP1_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_FPP1_dyVsdx_%d",m));
	h_dxhitFPP1[m] = h1_FPP1_dyVsdx_[m]->ProjectionX(Form("h1_dxhitFPP1_%d",m));
	h_dyhitFPP1[m] = h1_FPP1_dyVsdx_[m]->ProjectionY(Form("h1_dyhitFPP1_%d",m));
	
	h1_FPP1_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_FPP1_Edep_%d",m));
	//h1_FPP1_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_FPP1_Edep_log_%d",m));
	
	if(m==0){
	  h_EdephitFPP1 = h1_FPP1_Edep_[m];
	}else{
	  h_EdephitFPP1->Add(h1_FPP1_Edep_[m]);
	}
	//cout << m << " " << NhitsFPP1[m] << endl;
      }
    }
   
    if(det_list[k]=="fpp2"){
      TH1D* h1_FPP2_nhits_[5];
      TF1* f1_fpp2nhits_[5];
      TH2D* h1_FPP2_yVsx_[5];
      TH2D* h1_FPP2_dyVsdx_[5];
      TH1D* h1_FPP2_Edep_[5];
      
      cout << "FPP2" << endl;
      
      for(int m = 0; m<6; m++){
	//Nhits
	h1_FPP2_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_FPP2_nhits_%d",m));
	f1_fpp2nhits_[m] = new TF1(Form("f1_fpp2nhits_%d", m), "gaus", 0., 400.);
	h1_FPP2_nhits_[m]->Fit(f1_fpp2nhits_[m], "QRN");
	mu = f1_fpp2nhits_[m]->GetParameter(1);
	sigma = f1_fpp2nhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_fpp2nhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_FPP2_nhits_[m]->Fit(f1_fpp2nhits_[m], "QRN");
	}
	NhitsFPP2[m] = max(1.0, f1_fpp2nhits_[m]->GetParameter(1));
	
	h1_FPP2_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_FPP2_yVsx_%d",m));
	h_xhitFPP2[m] = h1_FPP2_yVsx_[m]->ProjectionX(Form("h1_xhitFPP2_%d",m));
	h_yhitFPP2[m] = h1_FPP2_yVsx_[m]->ProjectionY(Form("h1_yhitFPP2_%d",m));
	
	h1_FPP2_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_FPP2_dyVsdx_%d",m));
	h_dxhitFPP2[m] = h1_FPP2_dyVsdx_[m]->ProjectionX(Form("h1_dxhitFPP2_%d",m));
	h_dyhitFPP2[m] = h1_FPP2_dyVsdx_[m]->ProjectionY(Form("h1_dyhitFPP2_%d",m));
	
	h1_FPP2_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_FPP2_Edep_%d",m));
	//h1_FPP2_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_FPP2_Edep_log_%d",m));
	
	if(m==0){
	  h_EdephitFPP2 = h1_FPP2_Edep_[m];
	}else{
	  h_EdephitFPP2->Add(h1_FPP2_Edep_[m]);
	}
	//cout << m << " " << NhitsFPP2[m] << endl;
      }
    }

    if(det_list[k]=="ecal"){
      //ECal
      cout << "ECal" << endl;
      TH2D *h1_ECal_nhitsVsChan = (TH2D*)f_bkgd->Get("h1_ECal_nhitsVsChan");
      TH1D* h1_ECal_nhits_[1770];
      TF1* f1_ecalnhits_[1770];
      TH2D *h1_ECal_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_ECal_EdepHitVsChan");
      //TH2D *h1_ECal_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_ECal_EdepHitVsChan_log");
      
      for(int m = 0; m<1770; m++){
	h1_ECal_nhits_[m] = h1_ECal_nhitsVsChan->ProjectionY(Form("h1_ECal_nhits_%d", m), m+1, m+1);
	f1_ecalnhits_[m] = new TF1(Form("f1_ecalnhits_%d", m), "gaus", 0, 150);
	h1_ECal_nhits_[m]->Fit(f1_ecalnhits_[m], "QR0");
	mu = f1_ecalnhits_[m]->GetParameter(1);
	sigma = f1_ecalnhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_ecalnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_ECal_nhits_[m]->Fit(f1_ecalnhits_[m], "QR0");
	}
	
	NhitsECal[m] = max(1.0, f1_ecalnhits_[m]->GetParameter(1));
	//cout << m << " " << NhitsECal[m] << endl;
      }
      h_EdephitECal = h1_ECal_EdepHitVsChan->ProjectionY("h_EdephitECal");
    }

    if(det_list[k]=="cdet"){
      //CDet
      cout << "CDet" << endl;
      TH2D *h1_CDet_nhitsVsSlat = (TH2D*)f_bkgd->Get("h1_CDet_nhitsVsSlat");
      TH1D* h1_CDet_nhits_[2352];
      TF1* f1_cdetnhits_[2352];
      
      TH2D *h1_CDet_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_CDet_EdepHitVsSlat");
      //TH2D *h1_CDet_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_CDet_EdepHitVsSlat_log");
      TH2D *h1_CDet_xhitVsSlat = (TH2D*)f_bkgd->Get("h1_CDet_xhitVsSlat");
      
      for(int m = 0; m<90; m++){
	h1_CDet_nhits_[m] = h1_CDet_nhitsVsSlat->ProjectionY(Form("h1_CDet_nhits_%d", m), m+1, m+1);
	f1_cdetnhits_[m] = new TF1(Form("f1_cdetnhits_%d", m), "gaus", 0, 100);
	h1_CDet_nhits_[m]->Fit(f1_cdetnhits_[m], "QR0");
	mu = f1_cdetnhits_[m]->GetParameter(1);
	sigma = f1_cdetnhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_cdetnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_CDet_nhits_[m]->Fit(f1_cdetnhits_[m], "QR0");
	}
	NhitsCDet[m] = max(1.0, f1_cdetnhits_[m]->GetParameter(1));
	//cout << m << " " << NhitsCDet[m] << endl;
      }
      
      h_EdephitCDet = h1_CDet_EdepHitVsSlat->ProjectionY("h_EdephitCDet");
      h_xhitCDet = h1_CDet_xhitVsSlat->ProjectionY("h_xhitCDet");
    }
  }//end loop on detector list
}

  
void SBSDigBkgdGen::GenerateBkgd(TRandom3* R, 
				 std::vector<SBSDigPMTDet*> pmtdets,
				 std::vector<int> detmap, 
				 std::vector<SBSDigGEMDet*> gemdets, 
				 std::vector<int> gemmap, 
				 double lumifrac)
{
  // this function performs basically the same action as SBSDigAuxi::UnfoldData, except that instead of taking hits from an input tree, it generates hits from the histograms
  // see examples with BigBite GEMs and BigBite Hodoscope
  int nhits;
  int mod;
  double edep;
  int Npe;
  double t, p;
  double z_hit, Npe_Edep_ratio, sigma_tgen;
  double beta, sin2thetaC;
  double x_hit, y_hit;
  int idet = 0;
  
  if(fPMTBkgdDig){//
  
  //ordering by increasing unique det ID
  while(detmap[idet]!=HCAL_UNIQUE_DETID && idet<(int)detmap.size())idet++;
  if(idet>=detmap.size())idet = -1;
  if(idet>=0){
    //cout << "hcal" << endl;
    for(int m = 0; m<288; m++){
      nhits = R->Poisson(NhitsHCal[m]*lumifrac*pmtdets[idet]->fGateWidth/fTimeWindow);
      //cout << m << " " << NhitsHCal[m]*lumifrac << " " << nhits << endl;*
      //h_NhitsHCal_XC->Fill(m, nhits);
      for(int i = 0; i<nhits; i++){
	edep = h_EdephitHCal->GetRandom();//R);
	z_hit = h_zhitHCal->GetRandom();//R); //R->Uniform(-0.91, 0.91);//for the time being
	//h_EdephitHCal_XC->Fill(edep);
	//h_zhitHCal_XC->Fill(z_hit);
	
	//cout << " " << i << " " << edep << " " << z_hit << endl;
	Npe_Edep_ratio = 5.242+11.39*z_hit+10.41*pow(z_hit, 2);
	Npe = R->Poisson(Npe_Edep_ratio*edep*1.0e3);
	sigma_tgen = 0.4244+11380/pow(Npe+153.4, 2);
	
	t = R->Uniform(-pmtdets[idet]->fGateWidth/2., pmtdets[idet]->fGateWidth/2.);
	
	//if(edep>1.e-3)
	pmtdets[idet]->PMTmap[m].Fill_FADCmode1(Npe, pmtdets[idet]->fThreshold, t, sigma_tgen, 1);// edep > 1 MeV
      }
    }
  }
  
  while(detmap[idet]!=BBPS_UNIQUE_DETID && idet<(int)detmap.size())idet++;
  if(idet>=detmap.size())idet = -1;
  if(idet>=0){
    //cout << "ps" << endl;
    for(int m = 0; m<52; m++){
      nhits = R->Poisson(NhitsBBPS[m]*lumifrac*pmtdets[idet]->fGateWidth/fTimeWindow);
      //cout << m << " " << NhitsBBPS[m]*lumifrac << " " << nhits << endl;
      //h_NhitsBBPS_XC->Fill(m, nhits);
      for(int i = 0; i<nhits; i++){
	edep = h_EdephitBBPS->GetRandom();//R);
	
	//h_EdephitBBPS_XC->Fill(edep);
	
	if(edep<1.e-4)continue;
	//check probability to generate p.e. yield
	bool genpeyield = true;
	if(edep<1.e-2)
	  genpeyield = R->Uniform(0, 1)<=1-exp(0.29-950.*edep);
	//if we're go, let's generate
	if(genpeyield){
	  beta = sqrt( pow(m_e+edep, 2)-m_e*m_e )/(m_e + edep);
	  sin2thetaC = TMath::Max(1.-1./pow(n_leadglass*beta, 2), 0.);
	  //1500. Used to be 454.: just wrong
	  Npe = R->Poisson(300.0*edep*sin2thetaC/(1.-1./(n_leadglass*n_leadglass)) );
	  t = R->Uniform(-pmtdets[idet]->fGateWidth/2., pmtdets[idet]->fGateWidth/2.);
	  
	  //cout << " " << i << " " << edep << " " << Npe << endl;
	  //if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	  //if(edep>1.e-3)
	  pmtdets[idet]->PMTmap[m].Fill(pmtdets[idet]->fRefPulse, Npe, 0, t, 1);// edep > 1 MeV
	}
      }
    }
  }
  
  while(detmap[idet]!=BBSH_UNIQUE_DETID && idet<(int)detmap.size())idet++;
  if(idet>=detmap.size())idet = -1;
  if(idet>=0){
    //cout << "sh" << endl;
    for(int m = 0; m<189; m++){
      nhits = R->Poisson(NhitsBBSH[m]*lumifrac*pmtdets[idet]->fGateWidth/fTimeWindow);
      //h_NhitsBBSH_XC->Fill(m, nhits);
      for(int i = 0; i<nhits; i++){
	edep = h_EdephitBBSH->GetRandom();//R);
	
	//h_EdephitBBSH_XC->Fill(edep);
	
	if(edep<1.e-4)continue;
	//check probability to generate p.e. yield
	bool genpeyield = true;
	if(edep<1.e-2)
	  genpeyield = R->Uniform(0, 1)<=1-exp(0.29-950.*edep);
	//if we're go, let's generate
	if(genpeyield){
	  beta = sqrt( pow(m_e+edep, 2)-m_e*m_e )/(m_e + edep);
	  sin2thetaC = TMath::Max(1.-1./pow(n_leadglass*beta, 2), 0.);
	  //1500. Used to be 454.: just wrong
	  Npe = R->Poisson(360.0*edep*sin2thetaC/(1.-1./(n_leadglass*n_leadglass)) );
	  t = R->Uniform(-pmtdets[idet]->fGateWidth/2., pmtdets[idet]->fGateWidth/2.);
	  //if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	  //if(edep>1.e-3)
	  pmtdets[idet]->PMTmap[m].Fill(pmtdets[idet]->fRefPulse, Npe, 0, t, 1);// edep > 1 MeV
	}
      }
    }
  }
  
  while(detmap[idet]!=ECAL_UNIQUE_DETID && idet<(int)detmap.size())idet++;
  if(idet>=detmap.size())idet = -1;
  if(idet>=0){
    //cout << "sh" << endl;
    for(int m = 0; m<189; m++){
      nhits = R->Poisson(NhitsECal[m]*lumifrac*pmtdets[idet]->fGateWidth/fTimeWindow);
      //h_NhitsECal_XC->Fill(m, nhits);
      for(int i = 0; i<nhits; i++){
	edep = h_EdephitECal->GetRandom();//R);
	
	//h_EdephitECal_XC->Fill(edep);
	
	if(edep<1.e-4)continue;
	//check probability to generate p.e. yield
	bool genpeyield = true;
	if(edep<1.e-2)
	  genpeyield = R->Uniform(0, 1)<=1-exp(0.29-950.*edep);
	//if we're go, let's generate
	if(genpeyield){
	  beta = sqrt( pow(m_e+edep, 2)-m_e*m_e )/(m_e + edep);
	  sin2thetaC = TMath::Max(1.-1./pow(n_leadglass*beta, 2), 0.);
	  //1500. Used to be 454.: just wrong
	  Npe = R->Poisson(360.0*edep*sin2thetaC/(1.-1./(n_leadglass*n_leadglass)) );
	  t = R->Uniform(-pmtdets[idet]->fGateWidth/2., pmtdets[idet]->fGateWidth/2.);
	  //if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	  //if(edep>1.e-3)
	  pmtdets[idet]->PMTmap[m].Fill(pmtdets[idet]->fRefPulse, Npe, 0, t, 1);// edep > 1 MeV
	}
      }
    }
  }
  
  while(detmap[idet]!=GRINCH_UNIQUE_DETID && idet<(int)detmap.size())idet++;
  if(idet>=detmap.size())idet = -1;
  if(idet>=0){
    //cout << "grinch" << endl;
    for(int m = 0; m<510; m++){
      p = R->Uniform(0, 1);
      if(p<P2hitsGRINCH[m]*lumifrac*pmtdets[idet]->fGateWidth/fTimeWindow){
	nhits = 2;
      }else if(p<P1hitGRINCH[m]*lumifrac*pmtdets[idet]->fGateWidth/fTimeWindow)nhits = 1;
      
      //h_NhitsGRINCH_XC->Fill(m, nhits);
      
      for(int i = 0; i<nhits; i++){
	t = R->Uniform(-pmtdets[idet]->fGateWidth/2., pmtdets[idet]->fGateWidth/2.);
	Npe = h_NpeGRINCH->GetRandom();
	
	pmtdets[idet]->PMTmap[m].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, 1);
	
	//h_NpeGRINCH_XC->Fill(Npe);
      }
    }
  }
  
  // the block of code below is similar to the code that unfolds the data from the BigBite hodoscope in SBSDigAuxi::UnfoldData(...)
  while(detmap[idet]!=HODO_UNIQUE_DETID && idet<(int)detmap.size())idet++;
  if(idet>=detmap.size())idet = -1;
  //cout << " " << idet;
  // Process hodoscope data
  if(idet>=0){
    //cout << "hodo" << endl;
    for(int m = 0; m<90; m++){
      // determine the number of hits to generate, then loop on this number of hits
      nhits = R->Poisson(NhitsBBHodo[m]*lumifrac*pmtdets[idet]->fGateWidth/fTimeWindow);
      
      //h_NhitsBBHodo_XC->Fill(m, nhits);
            
      for(int i = 0; i<nhits; i++){
	// energy deposit, hit position
	// generated from sampling the histograms with function GetRandom();
	edep =  h_EdephitBBHodo->GetRandom();//*1.e6;
	//if(edep<0.002)continue;
	x_hit =  h_xhitBBHodo->GetRandom();
	
	//h_EdephitBBHodo_XC->Fill(edep);
	//h_xhitBBHodo_XC->Fill(x_hit);
	
	//p = R->Uniform(-50.,50.);
	// from that point, it's almost the same code as in 
	// SBSDigAuxi::UnfoldData(...), with one exception:
	// * t is generated uniformly within the DAQ time window.
	for(int j = 0; j<2; j++){//j = 0: close beam PMT, j = 1: far beam PMT
	  // Evaluation of number of photoelectrons and time from energy deposit documented at:
	  // https://hallaweb.jlab.org/dvcslog/SBS/170711_172759/BB_hodoscope_restudy_update_20170711.pdf
	  Npe = R->Poisson(1.0e7*edep*0.113187*exp(-(0.3+pow(-1, j)*x_hit)/1.03533)* 0.24);
	  t = R->Uniform(-pmtdets[idet]->fGateWidth/2., pmtdets[idet]->fGateWidth/2.);//+p+(0.55+pow(-1, j)*x_hit)/0.15;
	  //T->Earm_BBHodoScint_hit_sumedep->at(i);
	  //if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	  pmtdets[idet]->PMTmap[m*2+j].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, 1);
	}
      }
    }
  }
 //GEN-rp PRPOLBS_SCINT 
  while(detmap[idet]!=HODO_UNIQUE_DETID && idet<(int)detmap.size())idet++;
  if(idet>=detmap.size())idet = -1;
  //cout << " " << idet;
  // Process hodoscope data
  if(idet>=0){
    //cout << "hodo" << endl;
    for(int m = 0; m<90; m++){
      nhits = R->Poisson(NhitsPRPOLBS_SCINT[m]*lumifrac*pmtdets[idet]->fGateWidth/fTimeWindow);
      
      //h_NhitsPRPOLBS_SCINT_XC->Fill(m, nhits);
            
      for(int i = 0; i<nhits; i++){
	edep =  h_EdephitPRPOLBS_SCINT->GetRandom();//*1.e6;
	//if(edep<0.002)continue;
	x_hit =  h_xhitPRPOLBS_SCINT->GetRandom();
	
	//h_EdephitPRPOLBS_SCINT_XC->Fill(edep);
	//h_xhitPRPOLBS_SCINT_XC->Fill(x_hit);
	
	//p = R->Uniform(-50.,50.);
	for(int j = 0; j<2; j++){//j = 0: close beam PMT, j = 1: far beam PMT
	  // Evaluation of number of photoelectrons and time from energy deposit documented at:
	  // https://hallaweb.jlab.org/dvcslog/SBS/170711_172759/BB_hodoscope_restudy_update_20170711.pdf
	  Npe = R->Poisson(1.0e7*edep*0.113187*exp(-(0.3+pow(-1, j)*x_hit)/1.03533)* 0.24);
	  t = R->Uniform(-pmtdets[idet]->fGateWidth/2., pmtdets[idet]->fGateWidth/2.);//+p+(0.55+pow(-1, j)*x_hit)/0.15;
	  //T->Earm_PRPOLBS_SCINTScint_hit_sumedep->at(i);
	  //if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	  pmtdets[idet]->PMTmap[m*2+j].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, 1);
	}
      }
    }
  }
  
 //GEN-rp PRPOLFS_SCINT 
  while(detmap[idet]!=HODO_UNIQUE_DETID && idet<(int)detmap.size())idet++;
  if(idet>=detmap.size())idet = -1;
  //cout << " " << idet;
  // Process hodoscope data
  if(idet>=0){
    //cout << "hodo" << endl;
    for(int m = 0; m<90; m++){
      nhits = R->Poisson(NhitsPRPOLFS_SCINT[m]*lumifrac*pmtdets[idet]->fGateWidth/fTimeWindow);
      
      //h_NhitsPRPOLFS_SCINT_XC->Fill(m, nhits);
            
      for(int i = 0; i<nhits; i++){
	edep =  h_EdephitPRPOLFS_SCINT->GetRandom();//*1.e6;
	//if(edep<0.002)continue;
	x_hit =  h_xhitPRPOLFS_SCINT->GetRandom();
	
	//h_EdephitPRPOLFS_SCINT_XC->Fill(edep);
	//h_xhitPRPOLFS_SCINT_XC->Fill(x_hit);
	
	//p = R->Uniform(-50.,50.);
	for(int j = 0; j<2; j++){//j = 0: close beam PMT, j = 1: far beam PMT
	  // Evaluation of number of photoelectrons and time from energy deposit documented at:
	  // https://hallaweb.jlab.org/dvcslog/SBS/170711_172759/BB_hodoscope_restudy_update_20170711.pdf
	  Npe = R->Poisson(1.0e7*edep*0.113187*exp(-(0.3+pow(-1, j)*x_hit)/1.03533)* 0.24);
	  t = R->Uniform(-pmtdets[idet]->fGateWidth/2., pmtdets[idet]->fGateWidth/2.);//+p+(0.55+pow(-1, j)*x_hit)/0.15;
	  //T->Earm_PRPOLFS_SCINTScint_hit_sumedep->at(i);
	  //if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	  pmtdets[idet]->PMTmap[m*2+j].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, 1);
	}
      }
    }
  }
  
 //GEN-rp ACTIVEANA 
  while(detmap[idet]!=HODO_UNIQUE_DETID && idet<(int)detmap.size())idet++;
  if(idet>=detmap.size())idet = -1;
  //cout << " " << idet;
  // Process hodoscope data
  if(idet>=0){
    //cout << "hodo" << endl;
    for(int m = 0; m<90; m++){
      nhits = R->Poisson(NhitsACTIVEANA[m]*lumifrac*pmtdets[idet]->fGateWidth/fTimeWindow);
      
      //h_NhitsACTIVEANA_XC->Fill(m, nhits);
            
      for(int i = 0; i<nhits; i++){
	edep =  h_EdephitACTIVEANA->GetRandom();//*1.e6;
	//if(edep<0.002)continue;
	x_hit =  h_xhitACTIVEANA->GetRandom();
	
	//h_EdephitACTIVEANA_XC->Fill(edep);
	//h_xhitACTIVEANA_XC->Fill(x_hit);
	
	//p = R->Uniform(-50.,50.);
	for(int j = 0; j<2; j++){//j = 0: close beam PMT, j = 1: far beam PMT
	  // Evaluation of number of photoelectrons and time from energy deposit documented at:
	  // https://hallaweb.jlab.org/dvcslog/SBS/170711_172759/BB_hodoscope_restudy_update_20170711.pdf
	  Npe = R->Poisson(1.0e7*edep*0.113187*exp(-(0.3+pow(-1, j)*x_hit)/1.03533)* 0.24);
	  t = R->Uniform(-pmtdets[idet]->fGateWidth/2., pmtdets[idet]->fGateWidth/2.);//+p+(0.55+pow(-1, j)*x_hit)/0.15;
	  //T->Earm_ACTIVEANAScint_hit_sumedep->at(i);
	  //if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	  pmtdets[idet]->PMTmap[m*2+j].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, 1);
	}
      }
    }
  }
  
  while(detmap[idet]!=CDET_UNIQUE_DETID && idet<(int)detmap.size())idet++;
  if(idet>=detmap.size())idet = -1;
  // Process hodoscope data
  if(idet>=0){
    //cout << "hodo" << endl;
    for(int m = 0; m<2352; m++){
      nhits = R->Poisson(NhitsCDet[m]*lumifrac*pmtdets[idet]->fGateWidth/fTimeWindow);
      for(int i = 0; i<nhits; i++){
	edep =  h_EdephitCDet->GetRandom();//*1.e6;
	//if(edep<0.002)continue;
	x_hit =  h_xhitCDet->GetRandom();
	
	Npe = R->Poisson( edep*5.634e3 );
	t = R->Uniform(-pmtdets[idet]->fGateWidth/2., pmtdets[idet]->fGateWidth/2.);
	
	pmtdets[idet]->PMTmap[m].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, 1);
	
      }
    }
  }
  
  
  }//end if(fDetailedDig) 

  // the block of code below is similar to the code that unfolds the data from the BigBite GEMs in SBSDigAuxi::UnfoldData(...)
  idet = 0;
  while(gemmap[idet]!=BBGEM_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
  if(idet>=gemmap.size())idet = -1;
  
  if(idet>=0){
    //    cout << "bbgems" << endl;
    // loop on the GEM layers
    for(int m = 0; m<5; m++){
      // determine the number of hits to generate, then loop on this number of hits
      nhits = R->Poisson(NhitsBBGEMs[m]*lumifrac*gemdets[idet]->fGateWidth/fTimeWindow);
      //h_NhitsBBGEMs_XC[m]->Fill(nhits);
      for(int i = 0; i<nhits; i++){
	// energy deposit, hit position (at entrance of drift) 
	// generated from sampling the histograms with function GetRandom();
	edep =  h_EdephitBBGEMs->GetRandom();
	x_hit =  h_xhitBBGEMs[m]->GetRandom();
	y_hit =  h_yhitBBGEMs[m]->GetRandom();

	//h_EdephitBBGEMs_XC->Fill(edep);
	//h_xhitBBGEMs_XC[m]->Fill(x_hit);
	//h_yhitBBGEMs_XC[m]->Fill(y_hit);
	
	//x_hit =  h_dxhitBBGEMs[m]->GetRandom();
	//y_hit =  h_dyhitBBGEMs[m]->GetRandom();
	SBSDigGEMDet::gemhit hit; 
	hit.source = 1;
	
	mod = 0;
	// from that point, it's almost the same code as in 
	// SBSDigAuxi::UnfoldData(...), with a few exceptions:
	// * z_in, z_out are assumed to be -/+ 0.0015m respectively; 
	//   (that's not too big of an approximation)
	// * x(y)_out assumed to be x(y)_in+dx(y) 
	//   where dx(y) is sampled from the distribution of x(y)_out-x(y)_in.
	// * t is generated uniformly within the DAQ time window.
	while(mod<gemdets[idet]->fNPlanes/2){
	  if( (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5)<=x_hit && x_hit<=(gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) && m+1==gemdets[idet]->GEMPlanes[mod*2].Layer() )break;	  mod++;
	}//that does the job, but maybe can be optimized???
	if(mod==gemdets[idet]->fNPlanes/2)continue;
	//h_modBBGEMs_XC->Fill(mod);
	
	if(fabs(y_hit)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.)continue;
	
	//cout << x_hit << " " << y_hit << " " << mod << endl;
	
	hit.module = mod; 
	hit.edep = edep*1.0e6;//already in MeV for some reasons...
	
	hit.xin = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	hit.yin = y_hit;
	hit.zin = -0.0015;
	hit.xout = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset()+h_dxhitBBGEMs[m]->GetRandom();// 
	hit.yout = y_hit+h_dyhitBBGEMs[m]->GetRandom();
	//if(fabs(hit.yout)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.){
	//hit.yout *= gemdets[idet]->GEMPlanes[mod*2+1].dX()/2./hit.yout;
	//}
	hit.zout = 0.0015;
	hit.t = R->Uniform(-gemdets[idet]->fGateWidth/2.-50., gemdets[idet]->fGateWidth/2.-50.);
	
	gemdets[idet]->fGEMhits.push_back(hit);
      }
    }
  }
  
  while(gemmap[idet]!=SBSGEM_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
  if(idet>=gemmap.size())idet = -1;
  
  if(idet>=0){
    //    cout << "sbsgems" << endl;
    // loop on the GEM layers
    for(int m = 0; m<8; m++){
      // determine the number of hits to generate, then loop on this number of hits
      nhits = R->Poisson(NhitsSBSGEMs[m]*lumifrac*gemdets[idet]->fGateWidth/fTimeWindow);
      //h_NhitsSBSGEMs_XC[m]->Fill(nhits);
      for(int i = 0; i<nhits; i++){
	// energy deposit, hit position (at entrance of drift) 
	// generated from sampling the histograms with function GetRandom();
	edep =  h_EdephitSBSGEMs->GetRandom();
	x_hit =  h_xhitSBSGEMs[m]->GetRandom();
	y_hit =  h_yhitSBSGEMs[m]->GetRandom();

	//h_EdephitSBSGEMs_XC->Fill(edep);
	//h_xhitSBSGEMs_XC[m]->Fill(x_hit);
	//h_yhitSBSGEMs_XC[m]->Fill(y_hit);
	
	//x_hit =  h_dxhitSBSGEMs[m]->GetRandom();
	//y_hit =  h_dyhitSBSGEMs[m]->GetRandom();
	SBSDigGEMDet::gemhit hit; 
	hit.source = 1;
	
	mod = 0;
	// from that point, it's almost the same code as in 
	// SBSDigAuxi::UnfoldData(...), with a few exceptions:
	// * z_in, z_out are assumed to be -/+ 0.0015m respectively; 
	//   (that's not too big of an approximation)
	// * x(y)_out assumed to be x(y)_in+dx(y) 
	//   where dx(y) is sampled from the distribution of x(y)_out-x(y)_in.
	// * t is generated uniformly within the DAQ time window.
	while(mod<gemdets[idet]->fNPlanes/2){
	  if( (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5)<=x_hit && x_hit<=(gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) && m+1==gemdets[idet]->GEMPlanes[mod*2].Layer() )break;	  mod++;
	}//that does the job, but maybe can be optimized???
	if(mod==gemdets[idet]->fNPlanes/2)continue;
	//h_modSBSGEMs_XC->Fill(mod);
	
	if(fabs(y_hit)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.)continue;
	
	//cout << x_hit << " " << y_hit << " " << mod << endl;
	
	hit.module = mod; 
	hit.edep = edep*1.0e6;//already in MeV for some reasons...
	
	hit.xin = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	hit.yin = y_hit;
	hit.zin = -0.0015;
	hit.xout = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset()+h_dxhitSBSGEMs[m]->GetRandom();// 
	hit.yout = y_hit+h_dyhitSBSGEMs[m]->GetRandom();
	//if(fabs(hit.yout)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.){
	//hit.yout *= gemdets[idet]->GEMPlanes[mod*2+1].dX()/2./hit.yout;
	//}
	hit.zout = 0.0015;
	hit.t = R->Uniform(-gemdets[idet]->fGateWidth/2.-50., gemdets[idet]->fGateWidth/2.-50.);
	
	gemdets[idet]->fGEMhits.push_back(hit);
      }
    }
  }
    
  ///GEN-rp GEMs CEPOL_GEMFRONT
 idet = 0;
  while(gemmap[idet]!=CEPOL_GEMFRONT_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
  if(idet>=gemmap.size())idet = -1;
  
  if(idet>=0){
    //    cout << "CEPOL_GEMFRONTs" << endl;
    for(int m = 0; m<4; m++){
      nhits = R->Poisson(NhitsCEPOL_GEMFRONTs[m]*lumifrac*gemdets[idet]->fGateWidth/fTimeWindow);
      //h_NhitsCEPOL_GEMFRONTs_XC[m]->Fill(nhits);
      for(int i = 0; i<nhits; i++){
	edep =  h_EdephitCEPOL_GEMFRONTs->GetRandom();
	x_hit =  h_xhitCEPOL_GEMFRONTs[m]->GetRandom();
	y_hit =  h_yhitCEPOL_GEMFRONTs[m]->GetRandom();

	//h_EdephitCEPOL_GEMFRONTs_XC->Fill(edep);
	//h_xhitCEPOL_GEMFRONTs_XC[m]->Fill(x_hit);
	//h_yhitCEPOL_GEMFRONTs_XC[m]->Fill(y_hit);
	
	//x_hit =  h_dxhitCEPOL_GEMFRONTs[m]->GetRandom();
	//y_hit =  h_dyhitCEPOL_GEMFRONTs[m]->GetRandom();
	SBSDigGEMDet::gemhit hit; 
	hit.source = 1;
	
	mod = 0;
	
	while(mod<gemdets[idet]->fNPlanes/2){
	  if( (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5)<=x_hit && x_hit<=(gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) && m+1==gemdets[idet]->GEMPlanes[mod*2].Layer() )break;	  mod++;
	}//that does the job, but maybe can be optimized???
	if(mod==gemdets[idet]->fNPlanes/2)continue;
	//h_modCEPOL_GEMFRONTs_XC->Fill(mod);
	
	if(fabs(y_hit)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.)continue;
	
	//cout << x_hit << " " << y_hit << " " << mod << endl;
	
	hit.module = mod; 
	hit.edep = edep*1.0e6;//already in MeV for some reasons...
	
	hit.xin = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	hit.yin = y_hit;
	hit.zin = -0.0015;
	hit.xout = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset()+h_dxhitCEPOL_GEMFRONTs[m]->GetRandom();
	hit.yout = y_hit+h_dyhitCEPOL_GEMFRONTs[m]->GetRandom();
	//if(fabs(hit.yout)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.){
	//hit.yout *= gemdets[idet]->GEMPlanes[mod*2+1].dX()/2./hit.yout;
	//}
	hit.zout = 0.0015;
	hit.t = R->Uniform(-gemdets[idet]->fGateWidth/2.-50., gemdets[idet]->fGateWidth/2.-50.);
	
	gemdets[idet]->fGEMhits.push_back(hit);
      }
    }
  }
  
  ///GEN-rp GEMs CEPOL_GEMREAR
 idet = 0;
  while(gemmap[idet]!=CEPOL_GEMREAR_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
  if(idet>=gemmap.size())idet = -1;
  
  if(idet>=0){
    //    cout << "CEPOL_GEMREARs" << endl;
    for(int m = 0; m<4; m++){
      nhits = R->Poisson(NhitsCEPOL_GEMREARs[m]*lumifrac*gemdets[idet]->fGateWidth/fTimeWindow);
      //h_NhitsCEPOL_GEMREARs_XC[m]->Fill(nhits);
      for(int i = 0; i<nhits; i++){
	edep =  h_EdephitCEPOL_GEMREARs->GetRandom();
	x_hit =  h_xhitCEPOL_GEMREARs[m]->GetRandom();
	y_hit =  h_yhitCEPOL_GEMREARs[m]->GetRandom();

	//h_EdephitCEPOL_GEMREARs_XC->Fill(edep);
	//h_xhitCEPOL_GEMREARs_XC[m]->Fill(x_hit);
	//h_yhitCEPOL_GEMREARs_XC[m]->Fill(y_hit);
	
	//x_hit =  h_dxhitCEPOL_GEMREARs[m]->GetRandom();
	//y_hit =  h_dyhitCEPOL_GEMREARs[m]->GetRandom();
	SBSDigGEMDet::gemhit hit; 
	hit.source = 1;
	
	mod = 0;
	
	while(mod<gemdets[idet]->fNPlanes/2){
	  if( (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5)<=x_hit && x_hit<=(gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) && m+1==gemdets[idet]->GEMPlanes[mod*2].Layer() )break;	  mod++;
	}//that does the job, but maybe can be optimized???
	if(mod==gemdets[idet]->fNPlanes/2)continue;
	//h_modCEPOL_GEMREARs_XC->Fill(mod);
	
	if(fabs(y_hit)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.)continue;
	
	//cout << x_hit << " " << y_hit << " " << mod << endl;
	
	hit.module = mod; 
	hit.edep = edep*1.0e6;//already in MeV for some reasons...
	
	hit.xin = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	hit.yin = y_hit;
	hit.zin = -0.0015;
	hit.xout = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset()+h_dxhitCEPOL_GEMREARs[m]->GetRandom();
	hit.yout = y_hit+h_dyhitCEPOL_GEMREARs[m]->GetRandom();
	//if(fabs(hit.yout)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.){
	//hit.yout *= gemdets[idet]->GEMPlanes[mod*2+1].dX()/2./hit.yout;
	//}
	hit.zout = 0.0015;
	hit.t = R->Uniform(-gemdets[idet]->fGateWidth/2.-50., gemdets[idet]->fGateWidth/2.-50.);
	
	gemdets[idet]->fGEMhits.push_back(hit);
      }
    }
  }
  

  ///GEN-rp GEMs PRPOLBS_GEM
 idet = 0;
  while(gemmap[idet]!=PRPOLBS_GEM_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
  if(idet>=gemmap.size())idet = -1;
  
  if(idet>=0){
    //    cout << "PRPOLBS_GEMs" << endl;
    for(int m = 0; m<2; m++){
      nhits = R->Poisson(NhitsPRPOLBS_GEMs[m]*lumifrac*gemdets[idet]->fGateWidth/fTimeWindow);
      //h_NhitsPRPOLBS_GEMs_XC[m]->Fill(nhits);
      for(int i = 0; i<nhits; i++){
	edep =  h_EdephitPRPOLBS_GEMs->GetRandom();
	x_hit =  h_xhitPRPOLBS_GEMs[m]->GetRandom();
	y_hit =  h_yhitPRPOLBS_GEMs[m]->GetRandom();

	//h_EdephitPRPOLBS_GEMs_XC->Fill(edep);
	//h_xhitPRPOLBS_GEMs_XC[m]->Fill(x_hit);
	//h_yhitPRPOLBS_GEMs_XC[m]->Fill(y_hit);
	
	//x_hit =  h_dxhitPRPOLBS_GEMs[m]->GetRandom();
	//y_hit =  h_dyhitPRPOLBS_GEMs[m]->GetRandom();
	SBSDigGEMDet::gemhit hit; 
	hit.source = 1;
	
	mod = 0;
	
	while(mod<gemdets[idet]->fNPlanes/2){
	  if( (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5)<=x_hit && x_hit<=(gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) && m+1==gemdets[idet]->GEMPlanes[mod*2].Layer() )break;	  mod++;
	}//that does the job, but maybe can be optimized???
	if(mod==gemdets[idet]->fNPlanes/2)continue;
	//h_modPRPOLBS_GEMs_XC->Fill(mod);
	
	if(fabs(y_hit)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.)continue;
	
	//cout << x_hit << " " << y_hit << " " << mod << endl;
	
	hit.module = mod; 
	hit.edep = edep*1.0e6;//already in MeV for some reasons...
	
	hit.xin = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	hit.yin = y_hit;
	hit.zin = -0.0015;
	hit.xout = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset()+h_dxhitPRPOLBS_GEMs[m]->GetRandom();
	hit.yout = y_hit+h_dyhitPRPOLBS_GEMs[m]->GetRandom();
	//if(fabs(hit.yout)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.){
	//hit.yout *= gemdets[idet]->GEMPlanes[mod*2+1].dX()/2./hit.yout;
	//}
	hit.zout = 0.0015;
	hit.t = R->Uniform(-gemdets[idet]->fGateWidth/2.-50., gemdets[idet]->fGateWidth/2.-50.);
	
	gemdets[idet]->fGEMhits.push_back(hit);
      }
    }
  }
  
  ///GEN-rp GEMs PRPOLFS_GEM
 idet = 0;
  while(gemmap[idet]!=PRPOLFS_GEM_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
  if(idet>=gemmap.size())idet = -1;
  
  if(idet>=0){
    //    cout << "PRPOLFS_GEMs" << endl;
    for(int m = 0; m<2; m++){
      nhits = R->Poisson(NhitsPRPOLFS_GEMs[m]*lumifrac*gemdets[idet]->fGateWidth/fTimeWindow);
      //h_NhitsPRPOLFS_GEMs_XC[m]->Fill(nhits);
      for(int i = 0; i<nhits; i++){
	edep =  h_EdephitPRPOLFS_GEMs->GetRandom();
	x_hit =  h_xhitPRPOLFS_GEMs[m]->GetRandom();
	y_hit =  h_yhitPRPOLFS_GEMs[m]->GetRandom();

	//h_EdephitPRPOLFS_GEMs_XC->Fill(edep);
	//h_xhitPRPOLFS_GEMs_XC[m]->Fill(x_hit);
	//h_yhitPRPOLFS_GEMs_XC[m]->Fill(y_hit);
	
	//x_hit =  h_dxhitPRPOLFS_GEMs[m]->GetRandom();
	//y_hit =  h_dyhitPRPOLFS_GEMs[m]->GetRandom();
	SBSDigGEMDet::gemhit hit; 
	hit.source = 1;
	
	mod = 0;
	
	while(mod<gemdets[idet]->fNPlanes/2){
	  if( (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5)<=x_hit && x_hit<=(gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) && m+1==gemdets[idet]->GEMPlanes[mod*2].Layer() )break;	  mod++;
	}//that does the job, but maybe can be optimized???
	if(mod==gemdets[idet]->fNPlanes/2)continue;
	//h_modPRPOLFS_GEMs_XC->Fill(mod);
	
	if(fabs(y_hit)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.)continue;
	
	//cout << x_hit << " " << y_hit << " " << mod << endl;
	
	hit.module = mod; 
	hit.edep = edep*1.0e6;//already in MeV for some reasons...
	
	hit.xin = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	hit.yin = y_hit;
	hit.zin = -0.0015;
	hit.xout = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset()+h_dxhitPRPOLFS_GEMs[m]->GetRandom();
	hit.yout = y_hit+h_dyhitPRPOLFS_GEMs[m]->GetRandom();
	//if(fabs(hit.yout)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.){
	//hit.yout *= gemdets[idet]->GEMPlanes[mod*2+1].dX()/2./hit.yout;
	//}
	hit.zout = 0.0015;
	hit.t = R->Uniform(-gemdets[idet]->fGateWidth/2.-50., gemdets[idet]->fGateWidth/2.-50.);
	
	gemdets[idet]->fGEMhits.push_back(hit);
      }
    }
  }
  

  idet = 0;
  while(gemmap[idet]!=FT_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
  if(idet>=gemmap.size())idet = -1;
  
  if(idet>=0){
    //    cout << "ft" << endl;
    for(int m = 0; m<6; m++){
      nhits = R->Poisson(NhitsFT[m]*lumifrac*gemdets[idet]->fGateWidth/fTimeWindow);
      //h_NhitsFT_XC[m]->Fill(nhits);
      for(int i = 0; i<nhits; i++){
	edep =  h_EdephitFT->GetRandom();
	x_hit =  h_xhitFT[m]->GetRandom();
	y_hit =  h_yhitFT[m]->GetRandom();

	//h_EdephitFT_XC->Fill(edep);
	//h_xhitFT_XC[m]->Fill(x_hit);
	//h_yhitFT_XC[m]->Fill(y_hit);
	
	//x_hit =  h_dxhitFT[m]->GetRandom();
	//y_hit =  h_dyhitFT[m]->GetRandom();
	SBSDigGEMDet::gemhit hit; 
	hit.source = 1;
	
	mod = 0;
	
	while(mod<gemdets[idet]->fNPlanes/2){
	  if( (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5)<=x_hit && x_hit<=(gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) && m+1==gemdets[idet]->GEMPlanes[mod*2].Layer() )break;	  mod++;
	}//that does the job, but maybe can be optimized???
	if(mod==gemdets[idet]->fNPlanes/2)continue;
	//h_modFT_XC->Fill(mod);
	
	if(fabs(y_hit)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.)continue;
	
	//cout << x_hit << " " << y_hit << " " << mod << endl;
	
	hit.module = mod; 
	hit.edep = edep*1.0e6;//already in MeV for some reasons...
	
	hit.xin = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	hit.yin = y_hit;
	hit.zin = -0.0015;
	hit.xout = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset()+h_dxhitFT[m]->GetRandom();
	hit.yout = y_hit+h_dyhitFT[m]->GetRandom();
	//if(fabs(hit.yout)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.){
	//hit.yout *= gemdets[idet]->GEMPlanes[mod*2+1].dX()/2./hit.yout;
	//}
	hit.zout = 0.0015;
	hit.t = R->Uniform(-gemdets[idet]->fGateWidth/2.-50., gemdets[idet]->fGateWidth/2.-50.);
	
	gemdets[idet]->fGEMhits.push_back(hit);
      }
    }
  }
  
  idet = 0;
  while(gemmap[idet]!=FPP1_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
  if(idet>=gemmap.size())idet = -1;
  
  if(idet>=0){
    //    cout << "fpp1" << endl;
    for(int m = 0; m<5; m++){
      nhits = R->Poisson(NhitsFPP1[m]*lumifrac*gemdets[idet]->fGateWidth/fTimeWindow);
      //h_NhitsFPP1_XC[m]->Fill(nhits);
      for(int i = 0; i<nhits; i++){
	edep =  h_EdephitFPP1->GetRandom();
	x_hit =  h_xhitFPP1[m]->GetRandom();
	y_hit =  h_yhitFPP1[m]->GetRandom();

	//h_EdephitFPP1_XC->Fill(edep);
	//h_xhitFPP1_XC[m]->Fill(x_hit);
	//h_yhitFPP1_XC[m]->Fill(y_hit);
	
	//x_hit =  h_dxhitFPP1[m]->GetRandom();
	//y_hit =  h_dyhitFPP1[m]->GetRandom();
	SBSDigGEMDet::gemhit hit; 
	hit.source = 1;
	
	mod = 0;
	
	while(mod<gemdets[idet]->fNPlanes/2){
	  if( (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5)<=x_hit && x_hit<=(gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) && m+1==gemdets[idet]->GEMPlanes[mod*2].Layer() )break;	  mod++;
	}//that does the job, but maybe can be optimized???
	if(mod==gemdets[idet]->fNPlanes/2)continue;
	//h_modFPP1_XC->Fill(mod);
	
	if(fabs(y_hit)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.)continue;
	
	//cout << x_hit << " " << y_hit << " " << mod << endl;
	
	hit.module = mod; 
	hit.edep = edep*1.0e6;//already in MeV for some reasons...
	
	hit.xin = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	hit.yin = y_hit;
	hit.zin = -0.0015;
	hit.xout = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset()+h_dxhitFPP1[m]->GetRandom();
	hit.yout = y_hit+h_dyhitFPP1[m]->GetRandom();
	//if(fabs(hit.yout)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.){
	//hit.yout *= gemdets[idet]->GEMPlanes[mod*2+1].dX()/2./hit.yout;
	//}
	hit.zout = 0.0015;
	hit.t = R->Uniform(-gemdets[idet]->fGateWidth/2.-50., gemdets[idet]->fGateWidth/2.-50.);
	
	gemdets[idet]->fGEMhits.push_back(hit);
      }
    }
  }
  
  idet = 0;
  while(gemmap[idet]!=FPP2_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
  if(idet>=gemmap.size())idet = -1;
  
  if(idet>=0){
    //    cout << "fpp2" << endl;
    for(int m = 0; m<6; m++){
      nhits = R->Poisson(NhitsFPP2[m]*lumifrac*gemdets[idet]->fGateWidth/fTimeWindow);
      //h_NhitsFPP2_XC[m]->Fill(nhits);
      for(int i = 0; i<nhits; i++){
	edep =  h_EdephitFPP2->GetRandom();
	x_hit =  h_xhitFPP2[m]->GetRandom();
	y_hit =  h_yhitFPP2[m]->GetRandom();

	//h_EdephitFPP2_XC->Fill(edep);
	//h_xhitFPP2_XC[m]->Fill(x_hit);
	//h_yhitFPP2_XC[m]->Fill(y_hit);
	
	//x_hit =  h_dxhitFPP2[m]->GetRandom();
	//y_hit =  h_dyhitFPP2[m]->GetRandom();
	SBSDigGEMDet::gemhit hit; 
	hit.source = 1;
	
	mod = 0;
	
	while(mod<gemdets[idet]->fNPlanes/2){
	  if( (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5)<=x_hit && x_hit<=(gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) && m+1==gemdets[idet]->GEMPlanes[mod*2].Layer() )break;	  mod++;
	}//that does the job, but maybe can be optimized???
	if(mod==gemdets[idet]->fNPlanes/2)continue;
	//h_modFPP2_XC->Fill(mod);
	
	if(fabs(y_hit)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.)continue;
	
	//cout << x_hit << " " << y_hit << " " << mod << endl;
	
	hit.module = mod; 
	hit.edep = edep*1.0e6;//already in MeV for some reasons...
	
	hit.xin = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	hit.yin = y_hit;
	hit.zin = -0.0015;
	hit.xout = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset()+h_dxhitFPP2[m]->GetRandom();
	hit.yout = y_hit+h_dyhitFPP2[m]->GetRandom();
	//if(fabs(hit.yout)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.){
	//hit.yout *= gemdets[idet]->GEMPlanes[mod*2+1].dX()/2./hit.yout;
	//}
	hit.zout = 0.0015;
	hit.t = R->Uniform(-gemdets[idet]->fGateWidth/2.-50., gemdets[idet]->fGateWidth/2.-50.);
	
	gemdets[idet]->fGEMhits.push_back(hit);
      }
    }
  }

  
}

void SBSDigBkgdGen::WriteXCHistos()
{
  /*
  h_EdephitBBGEMs_XC->Write();
  for(int m = 0; m<5; m++){
    h_NhitsBBGEMs_XC[m]->Write();
    h_xhitBBGEMs_XC[m]->Write();
    h_yhitBBGEMs_XC[m]->Write();
  }
  h_modBBGEMs_XC->Write();
  
  h_NhitsHCal_XC->Write();
  h_EdephitHCal_XC->Write();
  h_zhitHCal_XC->Write();
  
  h_NhitsBBPS_XC->Write();
  h_EdephitBBPS_XC->Write();
  
  h_NhitsBBSH_XC->Write();
  h_EdephitBBSH_XC->Write();
  
  h_NhitsBBHodo_XC->Write();
  h_EdephitBBHodo_XC->Write();
  h_xhitBBHodo_XC->Write();
  
  h_NhitsGRINCH_XC->Write();
  h_NpeGRINCH_XC->Write();
  */
}
