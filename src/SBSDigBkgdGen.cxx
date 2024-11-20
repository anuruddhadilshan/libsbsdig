#include "g4sbs_types.h"
#include "SBSDigBkgdGen.h"

#include "TH2D.h"
#include "TF1.h"
#include "TMath.h"

using namespace std;

const UInt_t NHCalElements = 288;
const UInt_t NBBGEMPlanes = 5;
const UInt_t NBBPSElements = 52;
const UInt_t NBBSHElements = 189;
const UInt_t NBBHodoElements = 90;
const UInt_t NGRINCHElements = 510;
const UInt_t NECalElements = 1656;
const UInt_t NCDETElements = 2352;
const UInt_t NFTPlanes = 8;
const UInt_t NFPP1Planes = 8;
const UInt_t NSBSGEMPlanes = 8;
//const UInt_t NFPP2Planes = 5;//Legacy...
const UInt_t NPrPolBS_ScintElements = 24;
const UInt_t NPrPolFS_ScintElements = 24;
const UInt_t NActiveAnaElements = 32;
const UInt_t NCEPol_GEMFrontPlanes = 4;
const UInt_t NCEPol_GEMRearPlanes = 4;
const UInt_t NPrPolBS_GEMPlanes = 2;
const UInt_t NPrPolFS_GEMPlanes = 2;

SBSDigBkgdGen::SBSDigBkgdGen()
{
  //Initialization of arrays and histograms
  NhitsBBGEM = new Double_t[NBBGEMPlanes];
  h_xhitBBGEM = new TH1D*[NBBGEMPlanes];
  h_yhitBBGEM = new TH1D*[NBBGEMPlanes];
  h_dxhitBBGEM = new TH1D*[NBBGEMPlanes];
  h_dyhitBBGEM = new TH1D*[NBBGEMPlanes];
  NhitsHCal = new Double_t[NHCalElements];
  NhitsBBPS = new Double_t[NBBPSElements];
  NhitsBBSH = new Double_t[NBBSHElements];
  NhitsBBHodo = new Double_t[NBBHodoElements];
  NhitsPrPolBS_Scint = new Double_t[NPrPolBS_ScintElements];
  NhitsPrPolFS_Scint = new Double_t[NPrPolFS_ScintElements];
  NhitsActiveAna = new Double_t[NActiveAnaElements];
  P1hitGRINCH = new Double_t[NGRINCHElements];
  P2hitsGRINCH = new Double_t[NGRINCHElements];
  NhitsECal = new Double_t[NECalElements];
  NhitsCDet = new Double_t[NCDETElements];
  NhitsFT = new Double_t[NFTPlanes];
  h_xhitFT = new TH1D*[NFTPlanes];
  h_yhitFT = new TH1D*[NFTPlanes];
  h_dxhitFT = new TH1D*[NFTPlanes];
  h_dyhitFT = new TH1D*[NFTPlanes];
  NhitsFPP1 = new Double_t[NFPP1Planes];
  h_xhitFPP1 = new TH1D*[NFPP1Planes];
  h_yhitFPP1 = new TH1D*[NFPP1Planes];
  h_dxhitFPP1 = new TH1D*[NFPP1Planes];
  h_dyhitFPP1 = new TH1D*[NFPP1Planes];
  // NhitsFPP2 = new Double_t[5];
  // h_xhitFPP2 = new TH1D*[5];
  // h_yhitFPP2 = new TH1D*[5];
  // h_dxhitFPP2 = new TH1D*[5];
  // h_dyhitFPP2 = new TH1D*[5];

  NhitsSBSGEM = new Double_t[NSBSGEMPlanes];
  h_xhitSBSGEM = new TH1D*[NSBSGEMPlanes];
  h_yhitSBSGEM = new TH1D*[NSBSGEMPlanes];
  h_dxhitSBSGEM = new TH1D*[NSBSGEMPlanes];
  h_dyhitSBSGEM = new TH1D*[NSBSGEMPlanes];
  
  NhitsCEPol_GEMFront = new Double_t[NCEPol_GEMFrontPlanes];
  h_xhitCEPol_GEMFront = new TH1D*[NCEPol_GEMFrontPlanes];
  h_yhitCEPol_GEMFront = new TH1D*[NCEPol_GEMFrontPlanes];
  h_dxhitCEPol_GEMFront = new TH1D*[NCEPol_GEMFrontPlanes];
  h_dyhitCEPol_GEMFront = new TH1D*[NCEPol_GEMFrontPlanes];

  NhitsCEPol_GEMRear = new Double_t[NCEPol_GEMRearPlanes];
  h_xhitCEPol_GEMRear = new TH1D*[NCEPol_GEMRearPlanes];
  h_yhitCEPol_GEMRear = new TH1D*[NCEPol_GEMRearPlanes];
  h_dxhitCEPol_GEMRear = new TH1D*[NCEPol_GEMRearPlanes];
  h_dyhitCEPol_GEMRear = new TH1D*[NCEPol_GEMRearPlanes];
  
  NhitsPrPolBS_GEM = new Double_t[NPrPolBS_GEMPlanes];
  h_xhitPrPolBS_GEM = new TH1D*[NPrPolBS_GEMPlanes];
  h_yhitPrPolBS_GEM = new TH1D*[NPrPolBS_GEMPlanes];
  h_dxhitPrPolBS_GEM = new TH1D*[NPrPolBS_GEMPlanes];
  h_dyhitPrPolBS_GEM = new TH1D*[NPrPolBS_GEMPlanes];
 
  NhitsPrPolFS_GEM = new Double_t[NPrPolFS_GEMPlanes];
  h_xhitPrPolFS_GEM = new TH1D*[NPrPolFS_GEMPlanes];
  h_yhitPrPolFS_GEM = new TH1D*[NPrPolFS_GEMPlanes];
  h_dxhitPrPolFS_GEM = new TH1D*[NPrPolFS_GEMPlanes];
  h_dyhitPrPolFS_GEM = new TH1D*[NPrPolFS_GEMPlanes];
}

SBSDigBkgdGen::SBSDigBkgdGen(TFile* f_bkgd, std::vector<TString> det_list, double timewindow, bool pmtbkgddig)
{
  //Initialization of arrays and histograms
  fTimeWindow = timewindow;
  NhitsBBGEM = new Double_t[NBBGEMPlanes];
  h_xhitBBGEM = new TH1D*[NBBGEMPlanes];
  h_yhitBBGEM = new TH1D*[NBBGEMPlanes];
  h_dxhitBBGEM = new TH1D*[NBBGEMPlanes];
  h_dyhitBBGEM = new TH1D*[NBBGEMPlanes];
  NhitsHCal = new Double_t[NHCalElements];
  NhitsBBPS = new Double_t[NBBPSElements];
  NhitsBBSH = new Double_t[NBBSHElements];
  NhitsBBHodo = new Double_t[NBBHodoElements];
  NhitsPrPolBS_Scint = new Double_t[NPrPolBS_ScintElements];
  NhitsPrPolFS_Scint = new Double_t[NPrPolFS_ScintElements];
  NhitsActiveAna = new Double_t[NActiveAnaElements];
  P1hitGRINCH = new Double_t[NGRINCHElements];
  P2hitsGRINCH = new Double_t[NGRINCHElements];
  NhitsECal = new Double_t[NECalElements];
  NhitsCDet = new Double_t[NCDETElements];
  NhitsFT = new Double_t[NFTPlanes];
  h_xhitFT = new TH1D*[NFTPlanes];
  h_yhitFT = new TH1D*[NFTPlanes];
  h_dxhitFT = new TH1D*[NFTPlanes];
  h_dyhitFT = new TH1D*[NFTPlanes];
  NhitsFPP1 = new Double_t[NFPP1Planes];
  h_xhitFPP1 = new TH1D*[NFPP1Planes];
  h_yhitFPP1 = new TH1D*[NFPP1Planes];
  h_dxhitFPP1 = new TH1D*[NFPP1Planes];
  h_dyhitFPP1 = new TH1D*[NFPP1Planes];
  // NhitsFPP2 = new Double_t[5];
  // h_xhitFPP2 = new TH1D*[5];
  // h_yhitFPP2 = new TH1D*[5];
  // h_dxhitFPP2 = new TH1D*[5];
  // h_dyhitFPP2 = new TH1D*[5];

  NhitsSBSGEM = new Double_t[NSBSGEMPlanes];
  h_xhitSBSGEM = new TH1D*[NSBSGEMPlanes];
  h_yhitSBSGEM = new TH1D*[NSBSGEMPlanes];
  h_dxhitSBSGEM = new TH1D*[NSBSGEMPlanes];
  h_dyhitSBSGEM = new TH1D*[NSBSGEMPlanes];
  
  NhitsCEPol_GEMFront = new Double_t[NCEPol_GEMFrontPlanes];
  h_xhitCEPol_GEMFront = new TH1D*[NCEPol_GEMFrontPlanes];
  h_yhitCEPol_GEMFront = new TH1D*[NCEPol_GEMFrontPlanes];
  h_dxhitCEPol_GEMFront = new TH1D*[NCEPol_GEMFrontPlanes];
  h_dyhitCEPol_GEMFront = new TH1D*[NCEPol_GEMFrontPlanes];

  NhitsCEPol_GEMRear = new Double_t[NCEPol_GEMRearPlanes];
  h_xhitCEPol_GEMRear = new TH1D*[NCEPol_GEMRearPlanes];
  h_yhitCEPol_GEMRear = new TH1D*[NCEPol_GEMRearPlanes];
  h_dxhitCEPol_GEMRear = new TH1D*[NCEPol_GEMRearPlanes];
  h_dyhitCEPol_GEMRear = new TH1D*[NCEPol_GEMRearPlanes];
 
  NhitsPrPolBS_GEM = new Double_t[NPrPolBS_GEMPlanes];
  h_xhitPrPolBS_GEM = new TH1D*[NPrPolBS_GEMPlanes];
  h_yhitPrPolBS_GEM = new TH1D*[NPrPolBS_GEMPlanes];
  h_dxhitPrPolBS_GEM = new TH1D*[NPrPolBS_GEMPlanes];
  h_dyhitPrPolBS_GEM = new TH1D*[NPrPolBS_GEMPlanes];
 
  NhitsPrPolFS_GEM = new Double_t[NPrPolFS_GEMPlanes];
  h_xhitPrPolFS_GEM = new TH1D*[NPrPolFS_GEMPlanes];
  h_yhitPrPolFS_GEM = new TH1D*[NPrPolFS_GEMPlanes];
  h_dxhitPrPolFS_GEM = new TH1D*[NPrPolFS_GEMPlanes];
  h_dyhitPrPolFS_GEM = new TH1D*[NPrPolFS_GEMPlanes];
 
  

  fPMTBkgdDig = pmtbkgddig;
  
  Initialize(f_bkgd, det_list);
}

SBSDigBkgdGen::~SBSDigBkgdGen()
{
}

void SBSDigBkgdGen::Initialize(TFile* f_bkgd, std::vector<TString> det_list)
{
  double mu, sigma;
  
  // BB GEM
  TH1D* h1_BBGEM_nhits_[NBBGEMPlanes];
  //TF1* f1_bbgemnhits_[NBBGEMPlanes];
  TH2D* h1_BBGEM_yVsx_[NBBGEMPlanes];
  TH2D* h1_BBGEM_dyVsdx_[NBBGEMPlanes];
  TH1D* h1_BBGEM_Edep_[NBBGEMPlanes];
  
  TH1D* h1_SBSGEM_nhits_[NFPP1Planes];
  //TF1* f1_sbsgemnhits_[NFPP1Planes];
  TH2D* h1_SBSGEM_yVsx_[NFPP1Planes];
  TH2D* h1_SBSGEM_dyVsdx_[NFPP1Planes];
  TH1D* h1_SBSGEM_Edep_[NFPP1Planes];
  
  //CEPol_GEMFRONT
  TH1D* h1_CEPol_GEMFront_nhits_[NCEPol_GEMFrontPlanes];
  //TF1* f1_CEPol_GEMFrontnhits_[NCEPol_GEMFrontPlanes];
  TH2D* h1_CEPol_GEMFront_yVsx_[NCEPol_GEMFrontPlanes];
  TH2D* h1_CEPol_GEMFront_dyVsdx_[NCEPol_GEMFrontPlanes];
  TH1D* h1_CEPol_GEMFront_Edep_[NCEPol_GEMFrontPlanes];
  
  //CEPol_GEMRear
  TH1D* h1_CEPol_GEMRear_nhits_[NCEPol_GEMRearPlanes];
  //TF1* f1_CEPol_GEMRearnhits_[NCEPol_GEMRearPlanes];
  TH2D* h1_CEPol_GEMRear_yVsx_[NCEPol_GEMRearPlanes];
  TH2D* h1_CEPol_GEMRear_dyVsdx_[NCEPol_GEMRearPlanes];
  TH1D* h1_CEPol_GEMRear_Edep_[NCEPol_GEMRearPlanes];
  
  //PrPolBS_GEM
  TH1D* h1_PrPolBS_GEM_nhits_[NPrPolBS_GEMPlanes];
  //TF1* f1_PrPolBS_GEMnhits_[NPrPolBS_GEMPlanes];
  TH2D* h1_PrPolBS_GEM_yVsx_[NPrPolBS_GEMPlanes];
  TH2D* h1_PrPolBS_GEM_dyVsdx_[NPrPolBS_GEMPlanes];
  TH1D* h1_PrPolBS_GEM_Edep_[NPrPolBS_GEMPlanes];
  
  //PrPolFS_GEM
  TH1D* h1_PrPolFS_GEM_nhits_[NPrPolFS_GEMPlanes];
  //TF1* f1_PrPolFS_GEMnhits_[NPrPolFS_GEMPlanes];
  TH2D* h1_PrPolFS_GEM_yVsx_[NPrPolFS_GEMPlanes];
  TH2D* h1_PrPolFS_GEM_dyVsdx_[NPrPolFS_GEMPlanes];
  TH1D* h1_PrPolFS_GEM_Edep_[NPrPolFS_GEMPlanes];

  //FT GEM
  TH1D* h1_FT_nhits_[NFTPlanes];
  //TF1* f1_ftnhits_[NFTPlanes];
  TH2D* h1_FT_yVsx_[NFTPlanes];
  TH2D* h1_FT_dyVsdx_[NFTPlanes];
  TH1D* h1_FT_Edep_[NFTPlanes];
  
  //FPP1 GEM
  TH1D* h1_FPP1_nhits_[NSBSGEMPlanes];
  //TF1* f1_fpp1nhits_[NSBSGEMPlanes];
  TH2D* h1_FPP1_yVsx_[NSBSGEMPlanes];
  TH2D* h1_FPP1_dyVsdx_[NSBSGEMPlanes];
  TH1D* h1_FPP1_Edep_[NSBSGEMPlanes];
  
  /*
  //cross-check histos
  h_NhitsBBGEM_XC = new TH1D*[NBBGEMPlanes];
  h_EdephitBBGEM_XC = new TH1D("h_EdephitBBGEM_XC", "", 1000, 0.0, 0.1);
  h_xhitBBGEM_XC = new TH1D*[NBBGEMPlanes];
  h_yhitBBGEM_XC = new TH1D*[NBBGEMPlanes];
  h_modBBGEM_XC = new TH1D("h_modBBGEM_XC", "", 36, 0, 36);
  //TH1D* h_dxhitBBGEM_XC[NBBGEMPlanes];
  //TH1D* h_dyhitBBGEM_XC[NBBGEMPlanes];
  
  h_NhitsHCal_XC = new TH2D("h_NhitsHCal_XC", "", NHCalElements, 0, NHCalElements, 100, 0, 100);
  h_EdephitHCal_XC = new TH1D("h_EdephitHCal_XC", "", 100, 0.0+1.0e-3, 1.0+1.0e-3);
  h_zhitHCal_XC = new TH1D("h_zhitHCal_XC", "", 100, 0., 1.);
  
  h_NhitsBBPS_XC = new TH2D("h_NhitsBBPS_XC", "", NBBPSElements, 0, NBBPSElements, 150, 0, 150);
  h_EdephitBBPS_XC = new TH1D("h_EdephitBBPS_XC", "", 150, 0.0+1.0e-3, 1.5+1.0e-3);
  
  h_NhitsBBSH_XC = new TH2D("h_NhitsBBSH_XC", "", NBBSHElements, 0, NBBSHElements, 100, 0, 100);
  h_EdephitBBSH_XC = new TH1D("h_EdephitBBSH_XC", "", 150, 0.0+1.0e-3, 1.5+1.0e-3);
  
  h_NhitsPrPolFS_Scint_XC = new TH2D("h_NhitsBBHodo_XC", "", 90, 0, 90, 100, 0, 100);
  h_EdephitBBHodo_XC = new TH1D("h_EdephitBBHodo_XC", "", 250, 0., 0.5);
  h_xhitBBHodo_XC = new TH1D("h_xhitBBHodo_XC", "", 60, -0.3, 0.3);
  
  h_NhitsGRINCH_XC = new TH2D("h_NhitsGRINCH_XC", "", NGRINCHElements, 0, NGRINCHElements, 20, 0, 20);
  h_NpeGRINCH_XC = new TH1D("h_NpeGRINCH_XC", "", 100, 0, 100);
  */
  
  // Initialization of BBGEM histograms:
  // most histograms are 1D projections the 2D histograms stored in the input file.
  // for the energy deposit, the input file histograms (one per plane) are consolidated into one
  for(size_t k = 0; k<det_list.size(); k++){
    if(det_list[k]=="bbgem"){
      cout << "BB GEM" << endl;
      for(int m = 0; m<NBBGEMPlanes; m++){
	// // fit of the hits mulitplicity distribution.
	// h1_BBGEM_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_BBGEM_nhits_%d",m));
	// f1_bbgemnhits_[m] = new TF1(Form("f1_bbgemnhits_%d", m), "gaus", 0., 400.);
	// h1_BBGEM_nhits_[m]->Fit(f1_bbgemnhits_[m], "QRN");
	// mu = f1_bbgemnhits_[m]->GetParameter(1);
	// sigma = f1_bbgemnhits_[m]->GetParameter(2);
	// if(mu>=0){
	//   f1_bbgemnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	//   h1_BBGEM_nhits_[m]->Fit(f1_bbgemnhits_[m], "QRN");
	// }
	//Get mean number of hits per layer in background time window from histogram entries:
	
	//	NhitsBBGEM[m] = max(1.0, f1_bbgemnhits_[m]->GetParameter(1));
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_BBGEM_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_BBGEM_yVsx_%d",m));
	h_xhitBBGEM[m] = h1_BBGEM_yVsx_[m]->ProjectionX(Form("h1_xhitBBGEM_%d",m));
	h_yhitBBGEM[m] = h1_BBGEM_yVsx_[m]->ProjectionY(Form("h1_yhitBBGEM_%d",m));

	NhitsBBGEM[m] = h1_BBGEM_yVsx_[m]->GetEntries();
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_BBGEM_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_BBGEM_dyVsdx_%d",m));
	h_dxhitBBGEM[m] = h1_BBGEM_dyVsdx_[m]->ProjectionX(Form("h1_dxhitBBGEM_%d",m));
	h_dyhitBBGEM[m] = h1_BBGEM_dyVsdx_[m]->ProjectionY(Form("h1_dyhitBBGEM_%d",m));
	
	//h1_BBGEM_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_BBGEM_Edep_%d",m));
	h1_BBGEM_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_BBGEM_Edep_log_%d",m)); //this is log( edep (keV))
	
	if(m==0){
	  h_EdephitBBGEM = h1_BBGEM_Edep_[m];
	}else{
	  h_EdephitBBGEM->Add(h1_BBGEM_Edep_[m]);
	}
	
	/*
	//cross-check histos
	h_NhitsBBGEM_XC[m] = new TH1D(Form("h_NhitsBBGEM_XC_%d", m), "", 250, 0, 1000);
	h_xhitBBGEM_XC[m] = new TH1D(Form("h_xhitBBGEM_XC_%d", m), "", 205, -1.025, 1.025);
	h_yhitBBGEM_XC[m] = new TH1D(Form("h_yhitBBGEM_XC_%d", m), "", 62, -0.31, 0.31);
	//h_dxhitBBGEM_XC[m] = new TH1D(Form("h_dxhitBBGEM_XC_%d", m), "", 100, 0.05, 0.05);
	//h_dyhitBBGEM_XC[m] = new TH1D(Form("h_dyhitBBGEM_XC_%d", m), "", 100, -0.05, 0.05);
	*/
	
	
	//cout << m << " " << NhitsBBGEM[m] << endl;
      }
    }
  // Initialization of SBSGEM histograms:
  // most histograms are 1D projections the 2D histograms stored in the input file.
  // for the energy deposit, the input file histograms (one per plane) are consolidated into one
    if(det_list[k]=="sbsgem"){
      cout << "SBS GEM" << endl;
      for(int m = 0; m<NSBSGEMPlanes; m++){
	// fit of the hits mulitplicity distribution.
	// h1_SBSGEM_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_SBSGEM_nhits_%d",m));
	// f1_sbsgemnhits_[m] = new TF1(Form("f1_sbsgemnhits_%d", m), "gaus", 0., 400.);
	// h1_SBSGEM_nhits_[m]->Fit(f1_sbsgemnhits_[m], "QRN");
	// mu = f1_sbsgemnhits_[m]->GetParameter(1);
	// sigma = f1_sbsgemnhits_[m]->GetParameter(2);
	// if(mu>=0){
	//   f1_sbsgemnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	//   h1_SBSGEM_nhits_[m]->Fit(f1_sbsgemnhits_[m], "QRN");
	// }
	// NhitsSBSGEM[m] = max(1.0, f1_sbsgemnhits_[m]->GetParameter(1));
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_SBSGEM_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_SBSGEM_yVsx_%d",m));
	h_xhitSBSGEM[m] = h1_SBSGEM_yVsx_[m]->ProjectionX(Form("h1_xhitSBSGEM_%d",m));
	h_yhitSBSGEM[m] = h1_SBSGEM_yVsx_[m]->ProjectionY(Form("h1_yhitSBSGEM_%d",m));

	NhitsSBSGEM[m] = h1_SBSGEM_yVsx_[m]->GetEntries();
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_SBSGEM_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_SBSGEM_dyVsdx_%d",m));
	h_dxhitSBSGEM[m] = h1_SBSGEM_dyVsdx_[m]->ProjectionX(Form("h1_dxhitSBSGEM_%d",m));
	h_dyhitSBSGEM[m] = h1_SBSGEM_dyVsdx_[m]->ProjectionY(Form("h1_dyhitSBSGEM_%d",m));
	
	h1_SBSGEM_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_SBSGEM_Edep_%d",m));
	//h1_SBSGEM_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_SBSGEM_Edep_log_%d",m));
	
	if(m==0){
	  h_EdephitSBSGEM = h1_SBSGEM_Edep_[m];
	}else{
	  h_EdephitSBSGEM->Add(h1_SBSGEM_Edep_[m]);
	}
	
	/*
	//cross-check histos
	h_NhitsSBSGEM_XC[m] = new TH1D(Form("h_NhitsSBSGEM_XC_%d", m), "", 250, 0, 1000);
	h_xhitSBSGEM_XC[m] = new TH1D(Form("h_xhitSBSGEM_XC_%d", m), "", 205, -1.025, 1.025);
	h_yhitSBSGEM_XC[m] = new TH1D(Form("h_yhitSBSGEM_XC_%d", m), "", 62, -0.31, 0.31);
	//h_dxhitSBSGEM_XC[m] = new TH1D(Form("h_dxhitSBSGEM_XC_%d", m), "", 100, 0.05, 0.05);
	//h_dyhitSBSGEM_XC[m] = new TH1D(Form("h_dyhitSBSGEM_XC_%d", m), "", 100, -0.05, 0.05);
	*/
	
	
	//cout << m << " " << NhitsSBSGEM[m] << endl;
      }
    }
    
    if(det_list[k]=="cepol_front"){
      //CEPol_GEMFRONT
      
      cout << "CEPol_GEMFront GEM" << endl;
      
      for(int m = 0; m<NCEPol_GEMFrontPlanes; m++){
	// fit of the hits mulitplicity distribution.
	// h1_CEPol_GEMFront_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_CEPol_GEMFront_nhits_%d",m));
	// f1_CEPol_GEMFrontnhits_[m] = new TF1(Form("f1_CEPol_GEMFrontnhits_%d", m), "gaus", 0., 400.);
	// cout << m << " " << f1_CEPol_GEMFrontnhits_[m] << " " << h1_CEPol_GEMFront_nhits_[m] << endl;
	// h1_CEPol_GEMFront_nhits_[m]->Fit(f1_CEPol_GEMFrontnhits_[m], "QRN");
	// mu = f1_CEPol_GEMFrontnhits_[m]->GetParameter(1);
	// sigma = f1_CEPol_GEMFrontnhits_[m]->GetParameter(2);
	// if(mu>=0){
	//   f1_CEPol_GEMFrontnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	//   h1_CEPol_GEMFront_nhits_[m]->Fit(f1_CEPol_GEMFrontnhits_[m], "QRN");
	// }
	// NhitsCEPol_GEMFront[m] = max(1.0, f1_CEPol_GEMFrontnhits_[m]->GetParameter(1));
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_CEPol_GEMFront_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_CEPol_GEMFront_yVsx_%d",m));
	h_xhitCEPol_GEMFront[m] = h1_CEPol_GEMFront_yVsx_[m]->ProjectionX(Form("h1_xhitCEPol_GEMFront_%d",m));
	h_yhitCEPol_GEMFront[m] = h1_CEPol_GEMFront_yVsx_[m]->ProjectionY(Form("h1_yhitCEPol_GEMFront_%d",m));
	
	NhitsCEPol_GEMRear[m] = h1_CEPol_GEMFront_yVsx_[m]->GetEntries();
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_CEPol_GEMFront_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_CEPol_GEMFront_dyVsdx_%d",m));
	h_dxhitCEPol_GEMFront[m] = h1_CEPol_GEMFront_dyVsdx_[m]->ProjectionX(Form("h1_dxhitCEPol_GEMFront_%d",m));
	h_dyhitCEPol_GEMFront[m] = h1_CEPol_GEMFront_dyVsdx_[m]->ProjectionY(Form("h1_dyhitCEPol_GEMFront_%d",m));
	
	h1_CEPol_GEMFront_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_CEPol_GEMFront_Edep_%d",m));
	//h1_CEPol_GEMFront_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_CEPol_GEMFront_Edep_log_%d",m));
	
	if(m==0){
	  h_EdephitCEPol_GEMFront = h1_CEPol_GEMFront_Edep_[m];
	}else{
	  h_EdephitCEPol_GEMFront->Add(h1_CEPol_GEMFront_Edep_[m]);
	}
	
	//cout << m << " " << NhitsBBGEM[m] << endl;
      }
    }
    
    if(det_list[k]=="cepol_rear"){
      cout << "CEPol_GEMRear GEM" << endl;
      
      //CEPol_GEMRear
      for(int m = 0; m<NCEPol_GEMRearPlanes; m++){
	// fit of the hits mulitplicity distribution.
	// h1_CEPol_GEMRear_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_CEPol_GEMRear_nhits_%d",m));
	// f1_CEPol_GEMRearnhits_[m] = new TF1(Form("f1_CEPol_GEMRearnhits_%d", m), "gaus", 0., 400.);
	// h1_CEPol_GEMRear_nhits_[m]->Fit(f1_CEPol_GEMRearnhits_[m], "QRN");
	// mu = f1_CEPol_GEMRearnhits_[m]->GetParameter(1);
	// sigma = f1_CEPol_GEMRearnhits_[m]->GetParameter(2);
	// if(mu>=0){
	//   f1_CEPol_GEMRearnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	//   h1_CEPol_GEMRear_nhits_[m]->Fit(f1_CEPol_GEMRearnhits_[m], "QRN");
	// }
	// NhitsCEPol_GEMRear[m] = max(1.0, f1_CEPol_GEMRearnhits_[m]->GetParameter(1));
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_CEPol_GEMRear_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_CEPol_GEMRear_yVsx_%d",m));
	h_xhitCEPol_GEMRear[m] = h1_CEPol_GEMRear_yVsx_[m]->ProjectionX(Form("h1_xhitCEPol_GEMRear_%d",m));
	h_yhitCEPol_GEMRear[m] = h1_CEPol_GEMRear_yVsx_[m]->ProjectionY(Form("h1_yhitCEPol_GEMRear_%d",m));
	
	NhitsCEPol_GEMRear[m] = h1_CEPol_GEMRear_yVsx_[m]->GetEntries();
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_CEPol_GEMRear_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_CEPol_GEMRear_dyVsdx_%d",m));
	h_dxhitCEPol_GEMRear[m] = h1_CEPol_GEMRear_dyVsdx_[m]->ProjectionX(Form("h1_dxhitCEPol_GEMRear_%d",m));
	h_dyhitCEPol_GEMRear[m] = h1_CEPol_GEMRear_dyVsdx_[m]->ProjectionY(Form("h1_dyhitCEPol_GEMRear_%d",m));
	
	h1_CEPol_GEMRear_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_CEPol_GEMRear_Edep_%d",m));
	//h1_CEPol_GEMRear_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_CEPol_GEMRear_Edep_log_%d",m));
	
	if(m==0){
	  h_EdephitCEPol_GEMRear = h1_CEPol_GEMRear_Edep_[m];
	}else{
	  h_EdephitCEPol_GEMRear->Add(h1_CEPol_GEMRear_Edep_[m]);
	}
	
	//cout << m << " " << NhitsBBGEM[m] << endl;
      }
    }
    
    if(det_list[k]=="prpolbs_gem"){
      cout << "PrPolBS_GEM GEM" << endl;
      
      //PrPolBS_GEM
      for(int m = 0; m<NPrPolBS_GEMPlanes; m++){
	// fit of the hits mulitplicity distribution.
	// h1_PrPolBS_GEM_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_PrPolBS_GEM_nhits_%d",m));
	// f1_PrPolBS_GEMnhits_[m] = new TF1(Form("f1_PrPolBS_GEMnhits_%d", m), "gaus", 0., 400.);
	// h1_PrPolBS_GEM_nhits_[m]->Fit(f1_PrPolBS_GEMnhits_[m], "QRN");
	// mu = f1_PrPolBS_GEMnhits_[m]->GetParameter(1);
	// sigma = f1_PrPolBS_GEMnhits_[m]->GetParameter(2);
	// if(mu>=0){
	//   f1_PrPolBS_GEMnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	//   h1_PrPolBS_GEM_nhits_[m]->Fit(f1_PrPolBS_GEMnhits_[m], "QRN");
	// }
	// NhitsPrPolBS_GEM[m] = max(1.0, f1_PrPolBS_GEMnhits_[m]->GetParameter(1));
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_PrPolBS_GEM_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_PrPolBS_GEM_yVsx_%d",m));
	h_xhitPrPolBS_GEM[m] = h1_PrPolBS_GEM_yVsx_[m]->ProjectionX(Form("h1_xhitPrPolBS_GEM_%d",m));
	h_yhitPrPolBS_GEM[m] = h1_PrPolBS_GEM_yVsx_[m]->ProjectionY(Form("h1_yhitPrPolBS_GEM_%d",m));
	
	NhitsPrPolBS_GEM[m] = h1_PrPolBS_GEM_yVsx_[m]->GetEntries();

	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_PrPolBS_GEM_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_PrPolBS_GEM_dyVsdx_%d",m));
	h_dxhitPrPolBS_GEM[m] = h1_PrPolBS_GEM_dyVsdx_[m]->ProjectionX(Form("h1_dxhitPrPolBS_GEM_%d",m));
	h_dyhitPrPolBS_GEM[m] = h1_PrPolBS_GEM_dyVsdx_[m]->ProjectionY(Form("h1_dyhitPrPolBS_GEM_%d",m));
	
	h1_PrPolBS_GEM_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_PrPolBS_GEM_Edep_%d",m));
	//h1_PrPolBS_GEM_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_PrPolBS_GEM_Edep_log_%d",m));
	
	if(m==0){
	  h_EdephitPrPolBS_GEM = h1_PrPolBS_GEM_Edep_[m];
	}else{
	  h_EdephitPrPolBS_GEM->Add(h1_PrPolBS_GEM_Edep_[m]);
	}
	
	//cout << m << " " << NhitsBBGEM[m] << endl;
      }
    }
    
    if(det_list[k]=="prpolfs_gem"){
      cout << "PrPolFS_GEM GEM" << endl;
      
      //PrPolFS_GEM
      for(int m = 0; m<NPrPolFS_GEMPlanes; m++){
	// // fit of the hits mulitplicity distribution.
	// h1_PrPolFS_GEM_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_PrPolFS_GEM_nhits_%d",m));
	// f1_PrPolFS_GEMnhits_[m] = new TF1(Form("f1_PrPolFS_GEMnhits_%d", m), "gaus", 0., 400.);
	// h1_PrPolFS_GEM_nhits_[m]->Fit(f1_PrPolFS_GEMnhits_[m], "QRN");
	// mu = f1_PrPolFS_GEMnhits_[m]->GetParameter(1);
	// sigma = f1_PrPolFS_GEMnhits_[m]->GetParameter(2);
	// if(mu>=0){
	//   f1_PrPolFS_GEMnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	//   h1_PrPolFS_GEM_nhits_[m]->Fit(f1_PrPolFS_GEMnhits_[m], "QRN");
	// }
	// NhitsPrPolFS_GEM[m] = max(1.0, f1_PrPolFS_GEMnhits_[m]->GetParameter(1));
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_PrPolFS_GEM_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_PrPolFS_GEM_yVsx_%d",m));
	h_xhitPrPolFS_GEM[m] = h1_PrPolFS_GEM_yVsx_[m]->ProjectionX(Form("h1_xhitPrPolFS_GEM_%d",m));
	h_yhitPrPolFS_GEM[m] = h1_PrPolFS_GEM_yVsx_[m]->ProjectionY(Form("h1_yhitPrPolFS_GEM_%d",m));
	
	NhitsPrPolFS_GEM[m] = h1_PrPolFS_GEM_yVsx_[m]->GetEntries();
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_PrPolFS_GEM_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_PrPolFS_GEM_dyVsdx_%d",m));
	h_dxhitPrPolFS_GEM[m] = h1_PrPolFS_GEM_dyVsdx_[m]->ProjectionX(Form("h1_dxhitPrPolFS_GEM_%d",m));
	h_dyhitPrPolFS_GEM[m] = h1_PrPolFS_GEM_dyVsdx_[m]->ProjectionY(Form("h1_dyhitPrPolFS_GEM_%d",m));
	
	h1_PrPolFS_GEM_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_PrPolFS_GEM_Edep_%d",m));
	//h1_PrPolFS_GEM_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_PrPolFS_GEM_Edep_log_%d",m));
	
	if(m==0){
	  h_EdephitPrPolFS_GEM = h1_PrPolFS_GEM_Edep_[m];
	}else{
	  h_EdephitPrPolFS_GEM->Add(h1_PrPolFS_GEM_Edep_[m]);
	}
	
	//cout << m << " " << NhitsBBGEM[m] << endl;
      }
    }
    
    if(det_list[k]=="ft"){
      cout << "FT GEM" << endl;
      for(int m = 0; m<NFTPlanes; m++){
	// fit of the hits mulitplicity distribution.
	// h1_FT_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_FT_nhits_%d",m));
	// NhitsFT[m] = h1_FT_nhits_[m]->GetMean();
	
	// if(NhitsFT[m]>h1_FT_nhits_[m]->GetRMS()*5){
	//   f1_ftnhits_[m] = new TF1(Form("f1_ftnhits_%d", m), "gaus", 0., 400.);
	  
	//   h1_FT_nhits_[m]->Fit(f1_ftnhits_[m], "QRN");
	//   mu = f1_ftnhits_[m]->GetParameter(1);
	//   sigma = f1_ftnhits_[m]->GetParameter(2);
	//   if(mu>=0){
	//     f1_ftnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	//     h1_FT_nhits_[m]->Fit(f1_ftnhits_[m], "QRN");
	//   }
	//   NhitsFT[m] = f1_ftnhits_[m]->GetParameter(1);
	// }
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_FT_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_FT_yVsx_%d",m));
	h_xhitFT[m] = h1_FT_yVsx_[m]->ProjectionX(Form("h1_xhitFT_%d",m));
	h_yhitFT[m] = h1_FT_yVsx_[m]->ProjectionY(Form("h1_yhitFT_%d",m));

	NhitsFT[m] = h1_FT_yVsx_[m]->GetEntries();
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
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
	
	/*
	//cross-check histos
	h_NhitsFT_XC[m] = new TH1D(Form("h_NhitsFT_XC_%d", m), "", 250, 0, 1000);
	h_xhitFT_XC[m] = new TH1D(Form("h_xhitFT_XC_%d", m), "", 205, -1.025, 1.025);
	h_yhitFT_XC[m] = new TH1D(Form("h_yhitFT_XC_%d", m), "", 62, -0.31, 0.31);
	//h_dxhitFT_XC[m] = new TH1D(Form("h_dxhitFT_XC_%d", m), "", 100, 0.05, 0.05);
	//h_dyhitFT_XC[m] = new TH1D(Form("h_dyhitFT_XC_%d", m), "", 100, -0.05, 0.05);
	*/
	
	
	//cout << m << " " << NhitsFT[m] << endl;
      }
    }
    
    if(det_list[k]=="fpp1"){
      cout << "FPP1 GEM" << endl;
      for(int m = 0; m<NFPP1Planes; m++){
	// fit of the hits mulitplicity distribution.
	// h1_FPP1_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_FPP1_nhits_%d",m));
	// NhitsFPP1[m] = h1_FPP1_nhits_[m]->GetMean();
	
	// if(NhitsFPP1[m]>h1_FPP1_nhits_[m]->GetRMS()*5){
	//   f1_ftnhits_[m] = new TF1(Form("f1_ftnhits_%d", m), "gaus", 0., 400.);
	  
	//   h1_FPP1_nhits_[m]->Fit(f1_ftnhits_[m], "QRN");
	//   mu = f1_ftnhits_[m]->GetParameter(1);
	//   sigma = f1_ftnhits_[m]->GetParameter(2);
	//   if(mu>=0){
	//     f1_ftnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	//     h1_FPP1_nhits_[m]->Fit(f1_ftnhits_[m], "QRN");
	//   }
	//   NhitsFPP1[m] = f1_ftnhits_[m]->GetParameter(1);
	// }
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
	h1_FPP1_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_FPP1_yVsx_%d",m));
	h_xhitFPP1[m] = h1_FPP1_yVsx_[m]->ProjectionX(Form("h1_xhitFPP1_%d",m));
	h_yhitFPP1[m] = h1_FPP1_yVsx_[m]->ProjectionY(Form("h1_yhitFPP1_%d",m));
	
	NhitsFPP1[m] = h1_FT_yVsx_[m]->GetEntries();
	
	// copy of 2D position histograms, then projection to obtain 1D histograms
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
	
	/*
	//cross-check histos
	h_NhitsFPP1_XC[m] = new TH1D(Form("h_NhitsFPP1_XC_%d", m), "", 250, 0, 1000);
	h_xhitFPP1_XC[m] = new TH1D(Form("h_xhitFPP1_XC_%d", m), "", 205, -1.025, 1.025);
	h_yhitFPP1_XC[m] = new TH1D(Form("h_yhitFPP1_XC_%d", m), "", 62, -0.31, 0.31);
	//h_dxhitFPP1_XC[m] = new TH1D(Form("h_dxhitFPP1_XC_%d", m), "", 100, 0.05, 0.05);
	//h_dyhitFPP1_XC[m] = new TH1D(Form("h_dyhitFPP1_XC_%d", m), "", 100, -0.05, 0.05);
	*/
	
	
	//cout << m << " " << NhitsFPP1[m] << endl;
      }
    }
    
    if(det_list[k]=="hcal" && fPMTBkgdDig ){
      // Initialization of HCal histograms:
      // most histograms are 1D projections the 2D histograms stored in the input file.
      cout << "HCal" << endl;
      TH2D *h1_HCal_nhitsVsChan = (TH2D*)f_bkgd->Get("h1_HCal_nhitsVsChan");
      TH1D* h1_HCal_nhits_[NHCalElements];
      TF1* f1_hcalnhits_[NHCalElements];
      TH2D *h1_HCal_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_HCal_EdepHitVsChan");
      //TH2D *h1_HCal_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_HCal_EdepHitVsChan_log");
      TH2D *h1_HCal_zHitVsChan = (TH2D*)f_bkgd->Get("h1_HCal_zHitVsChan");
      
      for(int m = 0; m<NHCalElements; m++){
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
    
    if(det_list[k]=="bbps" && fPMTBkgdDig){
      //PS
      cout << "PS" << endl;
      TH2D *h1_BBPS_nhitsVsChan = (TH2D*)f_bkgd->Get("h1_BBPS_nhitsVsChan");
      TH1D* h1_BBPS_nhits_[NBBPSElements];
      TF1* f1_bbpsnhits_[NBBPSElements];
      TH2D *h1_BBPS_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_BBPS_EdepHitVsChan");
      //TH2D *h1_BBPS_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_BBPS_EdepHitVsChan_log");
      
      for(int m = 0; m<NBBPSElements; m++){
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
    
    if(det_list[k]=="bbsh" && fPMTBkgdDig){
      //SH
      cout << "SH" << endl;
      TH2D *h1_BBSH_nhitsVsChan = (TH2D*)f_bkgd->Get("h1_BBSH_nhitsVsChan");
      TH1D* h1_BBSH_nhits_[NBBSHElements];
      TF1* f1_bbshnhits_[NBBSHElements];
      TH2D *h1_BBSH_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_BBSH_EdepHitVsChan");
      //TH2D *h1_BBSH_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_BBSH_EdepHitVsChan_log");
      
      for(int m = 0; m<NBBSHElements; m++){
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
    
    if(det_list[k]=="bbhodo" && fPMTBkgdDig){
      //BB Hodo
      cout << "Hodo" << endl;
      TH2D *h1_BBHodo_nhitsVsSlat = (TH2D*)f_bkgd->Get("h1_BBHodo_nhitsVsSlat");
      TH1D* h1_BBHodo_nhits_[NBBHodoElements];
      TF1* f1_bbhodonhits_[NBBHodoElements];
      
      TH2D *h1_BBHodo_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_BBHodo_EdepHitVsSlat");
      //TH2D *h1_BBHodo_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_BBHodo_EdepHitVsSlat_log");
      TH2D *h1_BBHodo_xhitVsSlat = (TH2D*)f_bkgd->Get("h1_BBHodo_xhitVsSlat");
      
      for(int m = 0; m<NBBHodoElements; m++){
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
      //GEN-rp PrPolBS_Scint
      cout << "PrPolBS_Scint" << endl;
      TH2D *h1_PrPolBS_Scint_nhitsVsSlat = (TH2D*)f_bkgd->Get("h1_PrPolBS_Scint_nhitsVsSlat");
      TH1D* h1_PrPolBS_Scint_nhits_[NPrPolBS_ScintElements];
      TF1* f1_PrPolBS_Scintnhits_[NPrPolBS_ScintElements];
      
      TH2D *h1_PrPolBS_Scint_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_PrPolBS_Scint_EdepHitVsSlat");
      //TH2D *h1_PrPolBS_Scint_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_PrPolBS_Scint_EdepHitVsSlat_log");
      TH2D *h1_PrPolBS_Scint_xhitVsSlat = (TH2D*)f_bkgd->Get("h1_PrPolBS_Scint_xhitVsSlat");
      
      for(int m = 0; m<NPrPolBS_ScintElements; m++){
	h1_PrPolBS_Scint_nhits_[m] = h1_PrPolBS_Scint_nhitsVsSlat->ProjectionY(Form("h1_PrPolBS_Scint_nhits_%d", m), m+1, m+1);
	f1_PrPolBS_Scintnhits_[m] = new TF1(Form("f1_PrPolBS_Scintnhits_%d", m), "gaus", 0, 100);
	h1_PrPolBS_Scint_nhits_[m]->Fit(f1_PrPolBS_Scintnhits_[m], "QR0");
	mu = f1_PrPolBS_Scintnhits_[m]->GetParameter(1);
	sigma = f1_PrPolBS_Scintnhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_PrPolBS_Scintnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_PrPolBS_Scint_nhits_[m]->Fit(f1_PrPolBS_Scintnhits_[m], "QR0");
	}
	NhitsPrPolBS_Scint[m] = max(1.0, f1_PrPolBS_Scintnhits_[m]->GetParameter(1));
	//cout << m << " " << NhitsPrPolBS_Scint[m] << endl;
      }
      
      h_EdephitPrPolBS_Scint = h1_PrPolBS_Scint_EdepHitVsSlat->ProjectionY("h_EdephitPrPolBS_Scint");
      h_xhitPrPolBS_Scint = h1_PrPolBS_Scint_xhitVsSlat->ProjectionY("h_xhitPrPolBS_Scint");
    }
    
    if(det_list[k]=="prpolscint_fs"){
      //GEN-rp PrPolFS_Scint
      cout << "PrPolFS_Scint" << endl;
      TH2D *h1_PrPolFS_Scint_nhitsVsSlat = (TH2D*)f_bkgd->Get("h1_PrPolFS_Scint_nhitsVsSlat");
      TH1D* h1_PrPolFS_Scint_nhits_[NPrPolFS_ScintElements];
      TF1* f1_PrPolFS_Scintnhits_[NPrPolFS_ScintElements];
      
      TH2D *h1_PrPolFS_Scint_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_PrPolFS_Scint_EdepHitVsSlat");
      //TH2D *h1_PrPolFS_Scint_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_PrPolFS_Scint_EdepHitVsSlat_log");
      TH2D *h1_PrPolFS_Scint_xhitVsSlat = (TH2D*)f_bkgd->Get("h1_PrPolFS_Scint_xhitVsSlat");
      
      for(int m = 0; m<NPrPolFS_ScintElements; m++){
	h1_PrPolFS_Scint_nhits_[m] = h1_PrPolFS_Scint_nhitsVsSlat->ProjectionY(Form("h1_PrPolFS_Scint_nhits_%d", m), m+1, m+1);
	f1_PrPolFS_Scintnhits_[m] = new TF1(Form("f1_PrPolFS_Scintnhits_%d", m), "gaus", 0, 100);
	h1_PrPolFS_Scint_nhits_[m]->Fit(f1_PrPolFS_Scintnhits_[m], "QR0");
	mu = f1_PrPolFS_Scintnhits_[m]->GetParameter(1);
	sigma = f1_PrPolFS_Scintnhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_PrPolFS_Scintnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_PrPolFS_Scint_nhits_[m]->Fit(f1_PrPolFS_Scintnhits_[m], "QR0");
	}
	NhitsPrPolFS_Scint[m] = max(1.0, f1_PrPolFS_Scintnhits_[m]->GetParameter(1));
	//cout << m << " " << NhitsPrPolFS_Scint[m] << endl;
      }
      
      h_EdephitPrPolFS_Scint = h1_PrPolFS_Scint_EdepHitVsSlat->ProjectionY("h_EdephitPrPolFS_Scint");
      h_xhitPrPolFS_Scint = h1_PrPolFS_Scint_xhitVsSlat->ProjectionY("h_xhitPrPolFS_Scint");
    }

    if(det_list[k]=="activeana"){
      //GEN-rp ActiveAna
      cout << "ActiveAna" << endl;
      TH2D *h1_ActiveAna_nhitsVsSlat = (TH2D*)f_bkgd->Get("h1_ActiveAna_nhitsVsSlat");
      TH1D* h1_ActiveAna_nhits_[NActiveAnaElements];
      TF1* f1_ActiveAnanhits_[NActiveAnaElements];
      
      TH2D *h1_ActiveAna_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_ActiveAna_EdepHitVsSlat");
      //TH2D *h1_ActiveAna_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_ActiveAna_EdepHitVsSlat_log");
      TH2D *h1_ActiveAna_xhitVsSlat = (TH2D*)f_bkgd->Get("h1_ActiveAna_xhitVsSlat");
      
      for(int m = 0; m<NActiveAnaElements; m++){
	h1_ActiveAna_nhits_[m] = h1_ActiveAna_nhitsVsSlat->ProjectionY(Form("h1_ActiveAna_nhits_%d", m), m+1, m+1);
	f1_ActiveAnanhits_[m] = new TF1(Form("f1_ActiveAnanhits_%d", m), "gaus", 0, 100);
	h1_ActiveAna_nhits_[m]->Fit(f1_ActiveAnanhits_[m], "QR0");
	mu = f1_ActiveAnanhits_[m]->GetParameter(1);
	sigma = f1_ActiveAnanhits_[m]->GetParameter(2);
	if(mu>=0){
	  f1_ActiveAnanhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
	  h1_ActiveAna_nhits_[m]->Fit(f1_ActiveAnanhits_[m], "QR0");
	}
	NhitsActiveAna[m] = max(1.0, f1_ActiveAnanhits_[m]->GetParameter(1));
	//cout << m << " " << NhitsActiveAna[m] << endl;
      }
      
      h_EdephitActiveAna = h1_ActiveAna_EdepHitVsSlat->ProjectionY("h_EdephitActiveAna");
      h_xhitActiveAna = h1_ActiveAna_xhitVsSlat->ProjectionY("h_xhitActiveAna");
    }
    
    if(det_list[k]=="grinch" && fPMTBkgdDig){
      //GRINCH
      cout << "GRINCH" << endl;
      TH2D *h1_GRINCH_nhitsVsChan = (TH2D*)f_bkgd->Get("h1_GRINCH_nhitsVsChan");
      // TH1D* h1_GRINCH_Chan_1hit = h1_GRINCH_nhitsVsChan->ProjectionX("h1_GRINCH_Chan_1hit", 2, 2);
      // TH1D* h1_GRINCH_Chan_2hits = h1_GRINCH_nhitsVsChan->ProjectionX("h1_GRINCH_Chan_2hit", 3, -1);
      TH2D *h1_GRINCH_NpeVsChan = (TH2D*)f_bkgd->Get("h1_GRINCH_NpeVsChan");
      
      for(int m = 1; m<=NGRINCHElements; m++){
	P1hitGRINCH[m] = h1_GRINCH_nhitsVsChan->Integral(m, m, 2, 2)/h1_GRINCH_nhitsVsChan->Integral(m, m, 0, -1);
	P2hitsGRINCH[m] = h1_GRINCH_nhitsVsChan->Integral(m, m, 3, -1)/h1_GRINCH_nhitsVsChan->Integral(m, m, 0, -1);
      }
      h_NpeGRINCH = h1_GRINCH_NpeVsChan->ProjectionY("h1_GRINCH_Npe");
    }
    
    
    // if(det_list[k]=="ft"){
    //   TH1D* h1_FT_nhits_[NFTPlanes];
    //   TF1* f1_ftnhits_[NFTPlanes];
    //   TH2D* h1_FT_yVsx_[NFTPlanes];
    //   TH2D* h1_FT_dyVsdx_[NFTPlanes];
    //   TH1D* h1_FT_Edep_[NFTPlanes];
      
    //   cout << "FT" << endl;
      
    //   for(int m = 0; m<6; m++){
    // 	//Nhits
    // 	h1_FT_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_FT_nhits_%d",m));
    // 	f1_ftnhits_[m] = new TF1(Form("f1_ftnhits_%d", m), "gaus", 0., 400.);
    // 	h1_FT_nhits_[m]->Fit(f1_ftnhits_[m], "QRN");
    // 	mu = f1_ftnhits_[m]->GetParameter(1);
    // 	sigma = f1_ftnhits_[m]->GetParameter(2);
    // 	if(mu>=0){
    // 	  f1_ftnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
    // 	  h1_FT_nhits_[m]->Fit(f1_ftnhits_[m], "QRN");
    // 	}
    // 	NhitsFT[m] = max(1.0, f1_ftnhits_[m]->GetParameter(1));
	
    // 	h1_FT_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_FT_yVsx_%d",m));
    // 	h_xhitFT[m] = h1_FT_yVsx_[m]->ProjectionX(Form("h1_xhitFT_%d",m));
    // 	h_yhitFT[m] = h1_FT_yVsx_[m]->ProjectionY(Form("h1_yhitFT_%d",m));
	
    // 	h1_FT_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_FT_dyVsdx_%d",m));
    // 	h_dxhitFT[m] = h1_FT_dyVsdx_[m]->ProjectionX(Form("h1_dxhitFT_%d",m));
    // 	h_dyhitFT[m] = h1_FT_dyVsdx_[m]->ProjectionY(Form("h1_dyhitFT_%d",m));
	
    // 	h1_FT_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_FT_Edep_%d",m));
    // 	//h1_FT_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_FT_Edep_log_%d",m));
	
    // 	if(m==0){
    // 	  h_EdephitFT = h1_FT_Edep_[m];
    // 	}else{
    // 	  h_EdephitFT->Add(h1_FT_Edep_[m]);
    // 	}
    // 	//cout << m << " " << NhitsFT[m] << endl;
    //   }
    // }
  
    // if(det_list[k]=="fpp1"){
    //   TH1D* h1_FPP1_nhits_[5];
    //   TF1* f1_fpp1nhits_[5];
    //   TH2D* h1_FPP1_yVsx_[5];
    //   TH2D* h1_FPP1_dyVsdx_[5];
    //   TH1D* h1_FPP1_Edep_[5];
      
    //   cout << "FPP1" << endl;
      
    //   for(int m = 0; m<6; m++){
    // 	//Nhits
    // 	h1_FPP1_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_FPP1_nhits_%d",m));
    // 	f1_fpp1nhits_[m] = new TF1(Form("f1_fpp1nhits_%d", m), "gaus", 0., 400.);
    // 	h1_FPP1_nhits_[m]->Fit(f1_fpp1nhits_[m], "QRN");
    // 	mu = f1_fpp1nhits_[m]->GetParameter(1);
    // 	sigma = f1_fpp1nhits_[m]->GetParameter(2);
    // 	if(mu>=0){
    // 	  f1_fpp1nhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
    // 	  h1_FPP1_nhits_[m]->Fit(f1_fpp1nhits_[m], "QRN");
    // 	}
    // 	NhitsFPP1[m] = max(1.0, f1_fpp1nhits_[m]->GetParameter(1));
	
    // 	h1_FPP1_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_FPP1_yVsx_%d",m));
    // 	h_xhitFPP1[m] = h1_FPP1_yVsx_[m]->ProjectionX(Form("h1_xhitFPP1_%d",m));
    // 	h_yhitFPP1[m] = h1_FPP1_yVsx_[m]->ProjectionY(Form("h1_yhitFPP1_%d",m));
	
    // 	h1_FPP1_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_FPP1_dyVsdx_%d",m));
    // 	h_dxhitFPP1[m] = h1_FPP1_dyVsdx_[m]->ProjectionX(Form("h1_dxhitFPP1_%d",m));
    // 	h_dyhitFPP1[m] = h1_FPP1_dyVsdx_[m]->ProjectionY(Form("h1_dyhitFPP1_%d",m));
	
    // 	h1_FPP1_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_FPP1_Edep_%d",m));
    // 	//h1_FPP1_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_FPP1_Edep_log_%d",m));
	
    // 	if(m==0){
    // 	  h_EdephitFPP1 = h1_FPP1_Edep_[m];
    // 	}else{
    // 	  h_EdephitFPP1->Add(h1_FPP1_Edep_[m]);
    // 	}
    // 	//cout << m << " " << NhitsFPP1[m] << endl;
    //   }
    // }
   
    // if(det_list[k]=="fpp2"){
    //   TH1D* h1_FPP2_nhits_[5];
    //   TF1* f1_fpp2nhits_[5];
    //   TH2D* h1_FPP2_yVsx_[5];
    //   TH2D* h1_FPP2_dyVsdx_[5];
    //   TH1D* h1_FPP2_Edep_[5];
      
    //   cout << "FPP2" << endl;
      
    //   for(int m = 0; m<6; m++){
    // 	//Nhits
    // 	h1_FPP2_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_FPP2_nhits_%d",m));
    // 	f1_fpp2nhits_[m] = new TF1(Form("f1_fpp2nhits_%d", m), "gaus", 0., 400.);
    // 	h1_FPP2_nhits_[m]->Fit(f1_fpp2nhits_[m], "QRN");
    // 	mu = f1_fpp2nhits_[m]->GetParameter(1);
    // 	sigma = f1_fpp2nhits_[m]->GetParameter(2);
    // 	if(mu>=0){
    // 	  f1_fpp2nhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
    // 	  h1_FPP2_nhits_[m]->Fit(f1_fpp2nhits_[m], "QRN");
    // 	}
    // 	NhitsFPP2[m] = max(1.0, f1_fpp2nhits_[m]->GetParameter(1));
	
    // 	h1_FPP2_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_FPP2_yVsx_%d",m));
    // 	h_xhitFPP2[m] = h1_FPP2_yVsx_[m]->ProjectionX(Form("h1_xhitFPP2_%d",m));
    // 	h_yhitFPP2[m] = h1_FPP2_yVsx_[m]->ProjectionY(Form("h1_yhitFPP2_%d",m));
	
    // 	h1_FPP2_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_FPP2_dyVsdx_%d",m));
    // 	h_dxhitFPP2[m] = h1_FPP2_dyVsdx_[m]->ProjectionX(Form("h1_dxhitFPP2_%d",m));
    // 	h_dyhitFPP2[m] = h1_FPP2_dyVsdx_[m]->ProjectionY(Form("h1_dyhitFPP2_%d",m));
	
    // 	h1_FPP2_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_FPP2_Edep_%d",m));
    // 	//h1_FPP2_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_FPP2_Edep_log_%d",m));
	
    // 	if(m==0){
    // 	  h_EdephitFPP2 = h1_FPP2_Edep_[m];
    // 	}else{
    // 	  h_EdephitFPP2->Add(h1_FPP2_Edep_[m]);
    // 	}
    // 	//cout << m << " " << NhitsFPP2[m] << endl;
    //   }
    // }

    if(det_list[k]=="ecal" && fPMTBkgdDig ){
      //ECal
      cout << "ECal" << endl;
      TH2D *h1_ECal_nhitsVsChan = (TH2D*)f_bkgd->Get("h1_ECal_nhitsVsChan");
      TH1D* h1_ECal_nhits_[NECalElements];
      TF1* f1_ecalnhits_[NECalElements];
      TH2D *h1_ECal_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_ECal_EdepHitVsChan");
      //TH2D *h1_ECal_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_ECal_EdepHitVsChan_log");
      
      for(int m = 0; m<NECalElements; m++){
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

    if(det_list[k]=="cdet" && fPMTBkgdDig ){
      //CDet
      cout << "CDet" << endl;
      TH2D *h1_CDet_nhitsVsSlat = (TH2D*)f_bkgd->Get("h1_CDet_nhitsVsSlat");
      TH1D* h1_CDet_nhits_[NCDETElements];
      TF1* f1_cdetnhits_[NCDETElements];
      
      TH2D *h1_CDet_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_CDet_EdepHitVsSlat");
      //TH2D *h1_CDet_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_CDet_EdepHitVsSlat_log");
      TH2D *h1_CDet_xhitVsSlat = (TH2D*)f_bkgd->Get("h1_CDet_xhitVsSlat");
      
      for(int m = 0; m<NCDETElements; m++){
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
  // see examples with BigBite GEM and BigBite Hodoscope
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
      for(int m = 0; m<NHCalElements; m++){
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
      for(int m = 0; m<NBBPSElements; m++){
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
	    sin2thetaC = TMath::Max(1.-1./(pmtdets[idet]->fRefIndex*pmtdets[idet]->fRefIndex * beta*beta), 0.);
	    //1500. Used to be 454.: just wrong
	    Npe = R->Poisson(300.0*edep*sin2thetaC/(1.-1./(pmtdets[idet]->fRefIndex*pmtdets[idet]->fRefIndex)) );
	    t = R->Uniform(-pmtdets[idet]->fGateWidth/2., pmtdets[idet]->fGateWidth/2.);
	  
	    //cout << " " << i << " " << edep << " " << Npe << endl;
	    //if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	    //if(edep>1.e-3)
	    //pmtdets[idet]->PMTmap[m].Fill(pmtdets[idet]->fRefPulse, Npe, 0, t, 1);// edep > 1 MeV
	    pmtdets[idet]->PMTmap[m].Fill_FADCmode1(Npe, pmtdets[idet]->fThreshold, t, pmtdets[idet]->fSigmaPulse, 1);
	  }
	}
      }
    }
  
    while(detmap[idet]!=BBSH_UNIQUE_DETID && idet<(int)detmap.size())idet++;
    if(idet>=detmap.size())idet = -1;
    if(idet>=0){
      //cout << "sh" << endl;
      for(int m = 0; m<NBBSHElements; m++){
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
	    sin2thetaC = TMath::Max(1.-1./(pmtdets[idet]->fRefIndex*pmtdets[idet]->fRefIndex * beta*beta), 0.);
	    //1500. Used to be 454.: just wrong
	    Npe = R->Poisson(360.0*edep*sin2thetaC/(1.-1./(pmtdets[idet]->fRefIndex*pmtdets[idet]->fRefIndex)) );
	    t = R->Uniform(-pmtdets[idet]->fGateWidth/2., pmtdets[idet]->fGateWidth/2.);
	    //if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	    //if(edep>1.e-3)
	    //pmtdets[idet]->PMTmap[m].Fill(pmtdets[idet]->fRefPulse, Npe, 0, t, 1);// edep > 1 MeV
	    pmtdets[idet]->PMTmap[m].Fill_FADCmode1(Npe, pmtdets[idet]->fThreshold, t, pmtdets[idet]->fSigmaPulse, 1);
	  }
	}
      }
    }
  
    while(detmap[idet]!=ECAL_UNIQUE_DETID && idet<(int)detmap.size())idet++;
    if(idet>=detmap.size())idet = -1;
    if(idet>=0){
      //cout << "sh" << endl;
      for(int m = 0; m<NECalElements; m++){
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
	    sin2thetaC = TMath::Max(1.-1./(pmtdets[idet]->fRefIndex*pmtdets[idet]->fRefIndex * beta*beta), 0.);
	    //1500. Used to be 454.: just wrong
	    Npe = R->Poisson(360.0*edep*sin2thetaC/(1.-1./(pmtdets[idet]->fRefIndex*pmtdets[idet]->fRefIndex)) );
	    t = R->Uniform(-pmtdets[idet]->fGateWidth/2., pmtdets[idet]->fGateWidth/2.);
	    //if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	    //if(edep>1.e-3)
	    //pmtdets[idet]->PMTmap[m].Fill(pmtdets[idet]->fRefPulse, Npe, 0, t, 1);// edep > 1 MeV
	    pmtdets[idet]->PMTmap[m].Fill_FADCmode1(Npe, pmtdets[idet]->fThreshold, t, pmtdets[idet]->fSigmaPulse, 1);
	  }
	}
      }
    }
  
    while(detmap[idet]!=GRINCH_UNIQUE_DETID && idet<(int)detmap.size())idet++;
    if(idet>=detmap.size())idet = -1;
    if(idet>=0){
      //cout << "grinch" << endl;
      for(int m = 0; m<NGRINCHElements; m++){
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
      for(int m = 0; m<NBBHodoElements; m++){
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
    //GEN-rp PrPolBS_Scint 
    while(detmap[idet]!=PRPOLBS_SCINT_UNIQUE_DETID && idet<(int)detmap.size())idet++;
    if(idet>=detmap.size())idet = -1;
    //cout << " " << idet;
    // Process hodoscope data
    if(idet>=0){
      //cout << "hodo" << endl;
      for(int m = 0; m<NPrPolBS_ScintElements; m++){
	nhits = R->Poisson(NhitsPrPolBS_Scint[m]*lumifrac*pmtdets[idet]->fGateWidth/fTimeWindow);
      
	//h_NhitsPrPolBS_Scint_XC->Fill(m, nhits);
            
	for(int i = 0; i<nhits; i++){
	  edep =  h_EdephitPrPolBS_Scint->GetRandom();//*1.e6;
	  //if(edep<0.002)continue;
	  x_hit =  h_xhitPrPolBS_Scint->GetRandom();
	
	  //h_EdephitPrPolBS_Scint_XC->Fill(edep);
	  //h_xhitPrPolBS_Scint_XC->Fill(x_hit);
	
	  //p = R->Uniform(-50.,50.);
	  for(int j = 0; j<2; j++){//j = 0: close beam PMT, j = 1: far beam PMT
	    // Evaluation of number of photoelectrons and time from energy deposit documented at:
	    // https://hallaweb.jlab.org/dvcslog/SBS/170711_172759/BB_hodoscope_restudy_update_20170711.pdf
	    Npe = R->Poisson(1.0e7*edep*0.113187*exp(-(0.3+pow(-1, j)*x_hit)/1.03533)* 0.24);
	    t = R->Uniform(-pmtdets[idet]->fGateWidth/2., pmtdets[idet]->fGateWidth/2.);//+p+(0.55+pow(-1, j)*x_hit)/0.15;
	    //T->Earm_PrPolBS_ScintScint_hit_sumedep->at(i);
	    //if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	    pmtdets[idet]->PMTmap[m*2+j].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, 1);
	  }
	}
      }
    }
  
    //GEN-rp PrPolFS_Scint 
    while(detmap[idet]!=PRPOLFS_SCINT_UNIQUE_DETID && idet<(int)detmap.size())idet++;
    if(idet>=detmap.size())idet = -1;
    //cout << " " << idet;
    // Process hodoscope data
    if(idet>=0){
      //cout << "hodo" << endl;
      for(int m = 0; m<NPrPolFS_ScintElements; m++){
	nhits = R->Poisson(NhitsPrPolFS_Scint[m]*lumifrac*pmtdets[idet]->fGateWidth/fTimeWindow);
      
	//h_NhitsPrPolFS_Scint_XC->Fill(m, nhits);
            
	for(int i = 0; i<nhits; i++){
	  edep =  h_EdephitPrPolFS_Scint->GetRandom();//*1.e6;
	  //if(edep<0.002)continue;
	  x_hit =  h_xhitPrPolFS_Scint->GetRandom();
	
	  //h_EdephitPrPolFS_Scint_XC->Fill(edep);
	  //h_xhitPrPolFS_Scint_XC->Fill(x_hit);
	
	  //p = R->Uniform(-50.,50.);
	  for(int j = 0; j<2; j++){//j = 0: close beam PMT, j = 1: far beam PMT
	    // Evaluation of number of photoelectrons and time from energy deposit documented at:
	    // https://hallaweb.jlab.org/dvcslog/SBS/170711_172759/BB_hodoscope_restudy_update_20170711.pdf
	    Npe = R->Poisson(1.0e7*edep*0.113187*exp(-(0.3+pow(-1, j)*x_hit)/1.03533)* 0.24);
	    t = R->Uniform(-pmtdets[idet]->fGateWidth/2., pmtdets[idet]->fGateWidth/2.);//+p+(0.55+pow(-1, j)*x_hit)/0.15;
	    //T->Earm_PrPolFS_ScintScint_hit_sumedep->at(i);
	    //if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	    pmtdets[idet]->PMTmap[m*2+j].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, 1);
	  }
	}
      }
    }
  
    //GEN-rp ActiveAna 
    while(detmap[idet]!=ACTIVEANA_UNIQUE_DETID && idet<(int)detmap.size())idet++;
    if(idet>=detmap.size())idet = -1;
    //cout << " " << idet;
    // Process hodoscope data
    if(idet>=0){
      //cout << "hodo" << endl;
      for(int m = 0; m<NActiveAnaElements; m++){
	nhits = R->Poisson(NhitsActiveAna[m]*lumifrac*pmtdets[idet]->fGateWidth/fTimeWindow);
      
	//h_NhitsActiveAna_XC->Fill(m, nhits);
            
	for(int i = 0; i<nhits; i++){
	  edep =  h_EdephitActiveAna->GetRandom();//*1.e6;
	  //if(edep<0.002)continue;
	  x_hit =  h_xhitActiveAna->GetRandom();
	
	  //h_EdephitActiveAna_XC->Fill(edep);
	  //h_xhitActiveAna_XC->Fill(x_hit);
	
	  //p = R->Uniform(-50.,50.);
	  for(int j = 0; j<2; j++){//j = 0: close beam PMT, j = 1: far beam PMT
	    // Evaluation of number of photoelectrons and time from energy deposit documented at:
	    // https://hallaweb.jlab.org/dvcslog/SBS/170711_172759/BB_hodoscope_restudy_update_20170711.pdf
	    Npe = R->Poisson(1.0e7*edep*0.113187*exp(-(0.3+pow(-1, j)*x_hit)/1.03533)* 0.24);
	    t = R->Uniform(-pmtdets[idet]->fGateWidth/2., pmtdets[idet]->fGateWidth/2.);//+p+(0.55+pow(-1, j)*x_hit)/0.15;
	    //T->Earm_ActiveAnaScint_hit_sumedep->at(i);
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
      for(int m = 0; m<NCDETElements; m++){
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

  // the block of code below is similar to the code that unfolds the data from the BigBite GEM in SBSDigAuxi::UnfoldData(...)
  idet = 0;
  while(gemmap[idet]!=BBGEM_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
  if(idet>=gemmap.size())idet = -1;
  
  if(idet>=0){
    //    cout << "bbgems" << endl;
    // loop on the GEM layers
    for(int m = 0; m<NBBGEMPlanes; m++){
      // determine the number of hits to generate, then loop on this number of hits
      nhits = R->Poisson(NhitsBBGEM[m]*lumifrac*gemdets[idet]->fGateWidth/fTimeWindow);

      // cout << "layer, NhitsBBGEM[layer], gate width, time window, num bkgd hits = " << m << ", " 
      //  	   << NhitsBBGEM[m] << ", " << gemdets[idet]->fGateWidth << ", " << fTimeWindow << ", " << nhits << endl;
      //h_NhitsBBGEM_XC[m]->Fill(nhits);
      for(int i = 0; i<nhits; i++){
	// energy deposit, hit position (at entrance of drift) 
	// generated from sampling the histograms with function GetRandom();

	//now the edep histogram is actually log( Edep (keV) );

	double logedep_keV = h_EdephitBBGEM->GetRandom();
	edep = exp( logedep_keV )/1000.; //energy deposit in MeV!

	//cout << "hit, edep (MeV) = " << i << ", " << edep << endl;
	
	x_hit =  h_xhitBBGEM[m]->GetRandom();
	y_hit =  h_yhitBBGEM[m]->GetRandom();

	//h_EdephitBBGEM_XC->Fill(edep);
	//h_xhitBBGEM_XC[m]->Fill(x_hit);
	//h_yhitBBGEM_XC[m]->Fill(y_hit);
	
	//x_hit =  h_dxhitBBGEM[m]->GetRandom();
	//y_hit =  h_dyhitBBGEM[m]->GetRandom();
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
	//h_modBBGEM_XC->Fill(mod);
	
	if(fabs(y_hit)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.)continue;
	
	//cout << "(x_hit,y_hit,mod) = (" << x_hit << ", " << y_hit << ", " << mod << ")" << endl;
	
	hit.module = mod; 
	hit.edep = edep*1.0e6;//already in MeV for some reasons...
	
	hit.xin = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	hit.yin = y_hit;
	hit.zin = -0.0015;
	hit.xout = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset()+h_dxhitBBGEM[m]->GetRandom();// 
	hit.yout = y_hit+h_dyhitBBGEM[m]->GetRandom();
	//if(fabs(hit.yout)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.){
	//hit.yout *= gemdets[idet]->GEMPlanes[mod*2+1].dX()/2./hit.yout;
	//}
	hit.zout = 0.0015;
	hit.t = R->Uniform(-gemdets[idet]->fGateWidth/2.-50., gemdets[idet]->fGateWidth/2.-50.);

	//cout << "Adding new background hit " << i << endl;
	
	gemdets[idet]->fGEMhits.push_back(hit);
      }
    }
  }
  
  while(gemmap[idet]!=SBSGEM_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
  if(idet>=gemmap.size())idet = -1;
  
  if(idet>=0){
    //    cout << "sbsgems" << endl;
    // loop on the GEM layers
    for(int m = 0; m<NSBSGEMPlanes; m++){
      // determine the number of hits to generate, then loop on this number of hits
      nhits = R->Poisson(NhitsSBSGEM[m]*lumifrac*gemdets[idet]->fGateWidth/fTimeWindow);
      //h_NhitsSBSGEM_XC[m]->Fill(nhits);
      for(int i = 0; i<nhits; i++){
	// energy deposit, hit position (at entrance of drift) 
	// generated from sampling the histograms with function GetRandom();
	edep =  h_EdephitSBSGEM->GetRandom();
	x_hit =  h_xhitSBSGEM[m]->GetRandom();
	y_hit =  h_yhitSBSGEM[m]->GetRandom();

	//h_EdephitSBSGEM_XC->Fill(edep);
	//h_xhitSBSGEM_XC[m]->Fill(x_hit);
	//h_yhitSBSGEM_XC[m]->Fill(y_hit);
	
	//x_hit =  h_dxhitSBSGEM[m]->GetRandom();
	//y_hit =  h_dyhitSBSGEM[m]->GetRandom();
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
	//h_modSBSGEM_XC->Fill(mod);
	
	if(fabs(y_hit)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.)continue;
	
	//cout << x_hit << " " << y_hit << " " << mod << endl;
	
	hit.module = mod; 
	hit.edep = edep*1.0e6;//already in MeV for some reasons...
	
	hit.xin = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	hit.yin = y_hit;
	hit.zin = -0.0015;
	hit.xout = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset()+h_dxhitSBSGEM[m]->GetRandom();// 
	hit.yout = y_hit+h_dyhitSBSGEM[m]->GetRandom();
	//if(fabs(hit.yout)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.){
	//hit.yout *= gemdets[idet]->GEMPlanes[mod*2+1].dX()/2./hit.yout;
	//}
	hit.zout = 0.0015;
	hit.t = R->Uniform(-gemdets[idet]->fGateWidth/2.-50., gemdets[idet]->fGateWidth/2.-50.);
	
	gemdets[idet]->fGEMhits.push_back(hit);
      }
    }
  }
    
  ///GEN-rp GEM CEPol_GEMFRONT
 idet = 0;
  while(gemmap[idet]!=CEPOL_GEMFRONT_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
  if(idet>=gemmap.size())idet = -1;
  
  if(idet>=0){
    //    cout << "CEPol_GEMFront" << endl;
    for(int m = 0; m<NCEPol_GEMFrontPlanes; m++){
      nhits = R->Poisson(NhitsCEPol_GEMFront[m]*lumifrac*gemdets[idet]->fGateWidth/fTimeWindow);
      //h_NhitsCEPol_GEMFront_XC[m]->Fill(nhits);
      for(int i = 0; i<nhits; i++){
	edep =  h_EdephitCEPol_GEMFront->GetRandom();
	x_hit =  h_xhitCEPol_GEMFront[m]->GetRandom();
	y_hit =  h_yhitCEPol_GEMFront[m]->GetRandom();

	//h_EdephitCEPol_GEMFront_XC->Fill(edep);
	//h_xhitCEPol_GEMFront_XC[m]->Fill(x_hit);
	//h_yhitCEPol_GEMFront_XC[m]->Fill(y_hit);
	
	//x_hit =  h_dxhitCEPol_GEMFront[m]->GetRandom();
	//y_hit =  h_dyhitCEPol_GEMFront[m]->GetRandom();
	SBSDigGEMDet::gemhit hit; 
	hit.source = 1;
	
	mod = 0;
	
	while(mod<gemdets[idet]->fNPlanes/2){
	  if( (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5)<=x_hit && x_hit<=(gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) && m+1==gemdets[idet]->GEMPlanes[mod*2].Layer() )break;	  mod++;
	}//that does the job, but maybe can be optimized???
	if(mod==gemdets[idet]->fNPlanes/2)continue;
	//h_modCEPol_GEMFront_XC->Fill(mod);
	
	if(fabs(y_hit)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.)continue;
	
	//cout << x_hit << " " << y_hit << " " << mod << endl;
	
	hit.module = mod; 
	hit.edep = edep*1.0e6;//already in MeV for some reasons...
	
	hit.xin = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	hit.yin = y_hit;
	hit.zin = -0.0015;
	hit.xout = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset()+h_dxhitCEPol_GEMFront[m]->GetRandom();
	hit.yout = y_hit+h_dyhitCEPol_GEMFront[m]->GetRandom();
	//if(fabs(hit.yout)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.){
	//hit.yout *= gemdets[idet]->GEMPlanes[mod*2+1].dX()/2./hit.yout;
	//}
	hit.zout = 0.0015;
	hit.t = R->Uniform(-gemdets[idet]->fGateWidth/2.-50., gemdets[idet]->fGateWidth/2.-50.);
	
	gemdets[idet]->fGEMhits.push_back(hit);
      }
    }
  }
  
  ///GEN-rp GEM CEPol_GEMRear
 idet = 0;
  while(gemmap[idet]!=CEPOL_GEMREAR_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
  if(idet>=gemmap.size())idet = -1;
  
  if(idet>=0){
    //    cout << "CEPol_GEMRear" << endl;
    for(int m = 0; m<NCEPol_GEMRearPlanes; m++){
      nhits = R->Poisson(NhitsCEPol_GEMRear[m]*lumifrac*gemdets[idet]->fGateWidth/fTimeWindow);
      //h_NhitsCEPol_GEMRear_XC[m]->Fill(nhits);
      for(int i = 0; i<nhits; i++){
	edep =  h_EdephitCEPol_GEMRear->GetRandom();
	x_hit =  h_xhitCEPol_GEMRear[m]->GetRandom();
	y_hit =  h_yhitCEPol_GEMRear[m]->GetRandom();

	//h_EdephitCEPol_GEMRear_XC->Fill(edep);
	//h_xhitCEPol_GEMRear_XC[m]->Fill(x_hit);
	//h_yhitCEPol_GEMRear_XC[m]->Fill(y_hit);
	
	//x_hit =  h_dxhitCEPol_GEMRear[m]->GetRandom();
	//y_hit =  h_dyhitCEPol_GEMRear[m]->GetRandom();
	SBSDigGEMDet::gemhit hit; 
	hit.source = 1;
	
	mod = 0;
	
	while(mod<gemdets[idet]->fNPlanes/2){
	  if( (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5)<=x_hit && x_hit<=(gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) && m+1==gemdets[idet]->GEMPlanes[mod*2].Layer() )break;	  mod++;
	}//that does the job, but maybe can be optimized???
	if(mod==gemdets[idet]->fNPlanes/2)continue;
	//h_modCEPol_GEMRear_XC->Fill(mod);
	
	if(fabs(y_hit)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.)continue;
	
	//cout << x_hit << " " << y_hit << " " << mod << endl;
	
	hit.module = mod; 
	hit.edep = edep*1.0e6;//already in MeV for some reasons...
	
	hit.xin = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	hit.yin = y_hit;
	hit.zin = -0.0015;
	hit.xout = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset()+h_dxhitCEPol_GEMRear[m]->GetRandom();
	hit.yout = y_hit+h_dyhitCEPol_GEMRear[m]->GetRandom();
	//if(fabs(hit.yout)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.){
	//hit.yout *= gemdets[idet]->GEMPlanes[mod*2+1].dX()/2./hit.yout;
	//}
	hit.zout = 0.0015;
	hit.t = R->Uniform(-gemdets[idet]->fGateWidth/2.-50., gemdets[idet]->fGateWidth/2.-50.);
	
	gemdets[idet]->fGEMhits.push_back(hit);
      }
    }
  }
  

  ///GEN-rp GEM PrPolBS_GEM
 idet = 0;
  while(gemmap[idet]!=PRPOLBS_GEM_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
  if(idet>=gemmap.size())idet = -1;
  
  if(idet>=0){
    //    cout << "PrPolBS_GEM" << endl;
    for(int m = 0; m<NPrPolBS_GEMPlanes; m++){
      nhits = R->Poisson(NhitsPrPolBS_GEM[m]*lumifrac*gemdets[idet]->fGateWidth/fTimeWindow);
      //h_NhitsPrPolBS_GEM_XC[m]->Fill(nhits);
      for(int i = 0; i<nhits; i++){
	edep =  h_EdephitPrPolBS_GEM->GetRandom();
	x_hit =  h_xhitPrPolBS_GEM[m]->GetRandom();
	y_hit =  h_yhitPrPolBS_GEM[m]->GetRandom();

	//h_EdephitPrPolBS_GEM_XC->Fill(edep);
	//h_xhitPrPolBS_GEM_XC[m]->Fill(x_hit);
	//h_yhitPrPolBS_GEM_XC[m]->Fill(y_hit);
	
	//x_hit =  h_dxhitPrPolBS_GEM[m]->GetRandom();
	//y_hit =  h_dyhitPrPolBS_GEM[m]->GetRandom();
	SBSDigGEMDet::gemhit hit; 
	hit.source = 1;
	
	mod = 0;
	
	while(mod<gemdets[idet]->fNPlanes/2){
	  if( (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5)<=x_hit && x_hit<=(gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) && m+1==gemdets[idet]->GEMPlanes[mod*2].Layer() )break;	  mod++;
	}//that does the job, but maybe can be optimized???
	if(mod==gemdets[idet]->fNPlanes/2)continue;
	//h_modPrPolBS_GEM_XC->Fill(mod);
	
	if(fabs(y_hit)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.)continue;
	
	//cout << x_hit << " " << y_hit << " " << mod << endl;
	
	hit.module = mod; 
	hit.edep = edep*1.0e6;//already in MeV for some reasons...
	
	hit.xin = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	hit.yin = y_hit;
	hit.zin = -0.0015;
	hit.xout = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset()+h_dxhitPrPolBS_GEM[m]->GetRandom();
	hit.yout = y_hit+h_dyhitPrPolBS_GEM[m]->GetRandom();
	//if(fabs(hit.yout)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.){
	//hit.yout *= gemdets[idet]->GEMPlanes[mod*2+1].dX()/2./hit.yout;
	//}
	hit.zout = 0.0015;
	hit.t = R->Uniform(-gemdets[idet]->fGateWidth/2.-50., gemdets[idet]->fGateWidth/2.-50.);
	
	gemdets[idet]->fGEMhits.push_back(hit);
      }
    }
  }
  
  ///GEN-rp GEM PrPolFS_GEM
 idet = 0;
  while(gemmap[idet]!=PRPOLFS_GEM_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
  if(idet>=gemmap.size())idet = -1;
  
  if(idet>=0){
    //    cout << "PrPolFS_GEM" << endl;
    for(int m = 0; m<NPrPolFS_GEMPlanes; m++){
      nhits = R->Poisson(NhitsPrPolFS_GEM[m]*lumifrac*gemdets[idet]->fGateWidth/fTimeWindow);
      //h_NhitsPrPolFS_GEM_XC[m]->Fill(nhits);
      for(int i = 0; i<nhits; i++){
	edep =  h_EdephitPrPolFS_GEM->GetRandom();
	x_hit =  h_xhitPrPolFS_GEM[m]->GetRandom();
	y_hit =  h_yhitPrPolFS_GEM[m]->GetRandom();

	//h_EdephitPrPolFS_GEM_XC->Fill(edep);
	//h_xhitPrPolFS_GEM_XC[m]->Fill(x_hit);
	//h_yhitPrPolFS_GEM_XC[m]->Fill(y_hit);
	
	//x_hit =  h_dxhitPrPolFS_GEM[m]->GetRandom();
	//y_hit =  h_dyhitPrPolFS_GEM[m]->GetRandom();
	SBSDigGEMDet::gemhit hit; 
	hit.source = 1;
	
	mod = 0;
	
	while(mod<gemdets[idet]->fNPlanes/2){
	  if( (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5)<=x_hit && x_hit<=(gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) && m+1==gemdets[idet]->GEMPlanes[mod*2].Layer() )break;	  mod++;
	}//that does the job, but maybe can be optimized???
	if(mod==gemdets[idet]->fNPlanes/2)continue;
	//h_modPrPolFS_GEM_XC->Fill(mod);
	
	if(fabs(y_hit)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.)continue;
	
	//cout << x_hit << " " << y_hit << " " << mod << endl;
	
	hit.module = mod; 
	hit.edep = edep*1.0e6;//already in MeV for some reasons...
	
	hit.xin = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	hit.yin = y_hit;
	hit.zin = -0.0015;
	hit.xout = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset()+h_dxhitPrPolFS_GEM[m]->GetRandom();
	hit.yout = y_hit+h_dyhitPrPolFS_GEM[m]->GetRandom();
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
    for(int m = 0; m<NFTPlanes; m++){
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
    for(int m = 0; m<NFPP1Planes; m++){
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
  
  
  // idet = 0;
  // while(gemmap[idet]!=FPP2_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
  // if(idet>=gemmap.size())idet = -1;
  
  // if(idet>=0){
  //   //    cout << "fpp2" << endl;
  //   for(int m = 0; m<6; m++){
  //     nhits = R->Poisson(NhitsFPP2[m]*lumifrac*gemdets[idet]->fGateWidth/fTimeWindow);
  //     //h_NhitsFPP2_XC[m]->Fill(nhits);
  //     for(int i = 0; i<nhits; i++){
  // 	edep =  h_EdephitFPP2->GetRandom();
  // 	x_hit =  h_xhitFPP2[m]->GetRandom();
  // 	y_hit =  h_yhitFPP2[m]->GetRandom();

  // 	//h_EdephitFPP2_XC->Fill(edep);
  // 	//h_xhitFPP2_XC[m]->Fill(x_hit);
  // 	//h_yhitFPP2_XC[m]->Fill(y_hit);
	
  // 	//x_hit =  h_dxhitFPP2[m]->GetRandom();
  // 	//y_hit =  h_dyhitFPP2[m]->GetRandom();
  // 	SBSDigGEMDet::gemhit hit; 
  // 	hit.source = 1;
	
  // 	mod = 0;
	
  // 	while(mod<gemdets[idet]->fNPlanes/2){
  // 	  if( (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5)<=x_hit && x_hit<=(gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) && m+1==gemdets[idet]->GEMPlanes[mod*2].Layer() )break;	  mod++;
  // 	}//that does the job, but maybe can be optimized???
  // 	if(mod==gemdets[idet]->fNPlanes/2)continue;
  // 	//h_modFPP2_XC->Fill(mod);
	
  // 	if(fabs(y_hit)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.)continue;
	
  // 	//cout << x_hit << " " << y_hit << " " << mod << endl;
	
  // 	hit.module = mod; 
  // 	hit.edep = edep*1.0e6;//already in MeV for some reasons...
	
  // 	hit.xin = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
  // 	hit.yin = y_hit;
  // 	hit.zin = -0.0015;
  // 	hit.xout = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset()+h_dxhitFPP2[m]->GetRandom();
  // 	hit.yout = y_hit+h_dyhitFPP2[m]->GetRandom();
  // 	//if(fabs(hit.yout)>gemdets[idet]->GEMPlanes[mod*2+1].dX()/2.){
  // 	//hit.yout *= gemdets[idet]->GEMPlanes[mod*2+1].dX()/2./hit.yout;
  // 	//}
  // 	hit.zout = 0.0015;
  // 	hit.t = R->Uniform(-gemdets[idet]->fGateWidth/2.-50., gemdets[idet]->fGateWidth/2.-50.);
	
  // 	gemdets[idet]->fGEMhits.push_back(hit);
  //     }
  //   }
  // }

  
}

void SBSDigBkgdGen::WriteXCHistos()
{
  /*
  h_EdephitBBGEM_XC->Write();
  for(int m = 0; m<5; m++){
    h_NhitsBBGEM_XC[m]->Write();
    h_xhitBBGEM_XC[m]->Write();
    h_yhitBBGEM_XC[m]->Write();
  }
  h_modBBGEM_XC->Write();
  
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
