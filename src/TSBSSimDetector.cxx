#include "TSBSSimDetector.h"
#include "TSBSDBManager.h"

TSBSSimDetector::TSBSSimDetector() : fHasData(false)
{
  fDBmanager = TSBSDBManager::GetInstance();
  
}

TSBSSimDetector::~TSBSSimDetector()
{
  fDBmanager->Delete();
}

TSPEModel::TSPEModel(DigInfo diginfo, const char* detname):
  fDigInfo(diginfo)//, qe(1.602e-19), unit(1e-9)
{
  fScale = fDigInfo.fROImpedance*qe/spe_unit;
  
  if(fDigInfo.fGain.size()==1){
    fScale*= fDigInfo.fGain[0];
  }
  
  //start_t = -12.5;
  double mint = fDigInfo.fTriggerOffset-fDigInfo.fGateWidth/2.0;
  double maxt = fDigInfo.fTriggerOffset+fDigInfo.fGateWidth/2.0;
  // test values
  double tau = fDigInfo.fSPEtau;
  double sig = fDigInfo.fSPEsig;
  double t0 = fDigInfo.fSPEtransittime;
  
  TF1 fFunc1(Form("fFunc1%s",detname),
	     TString::Format("TMath::Max(0.,(x/%g)*TMath::Exp(-x/(%g)))",
			     tau*tau,tau),
	     mint,maxt);
  TF1 fFunc2(Form("fFunc2%s",detname),
	     TString::Format("%g*TMath::Exp(-((x-%g)**2)/(%g))",
			     1./TMath::Sqrt(2*TMath::Pi()*sig),t0,sig*sig),
	     mint,maxt);
  TF1Convolution fConvolution(&fFunc1,&fFunc2);
  
  model = new TF1(Form("fSignal%s",detname),fConvolution,mint,maxt, fConvolution.GetNpar());
}

double TSPEModel::Eval(double t, int chan)
{
  if(fDigInfo.fGain.size()>1){
    if(fDigInfo.fGain.size()<=chan){
      cout << "warning: requested channel number " << chan << "larger than number of channel size " << fDigInfo.fGain.size() << " check database ! " << endl;
      exit(-1);
    }
    fScale*= fDigInfo.fGain[chan];
  }
  
  return fScale*model->Eval(t);
  //return model->Eval(t);
  //return 1.0;
}
