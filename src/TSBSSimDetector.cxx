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

SPEModel::SPEModel(DigInfo diginfo, const char* detname):
  fDigInfo(diginfo), qe(1.602e-19), unit(1e-9)
{
  gain_pmt = fDigInfo.fGain;  
  resistance = fDigInfo.fROImpedance;
  scale = gain_pmt*resistance*qe/unit;
  //start_t = -12.5;
  mint = fDigInfo.fTriggerOffset-fDigInfo.fGateWidth/2.0;
  maxt = fDigInfo.fTriggerOffset+fDigInfo.fGateWidth/2.0;
  // test values
  tao = 2.08;
  sig = 2.20;
  t0 = 5.0;
  fFunc1 = new TF1(Form("fFunc1%s",detname),
		   TString::Format("TMath::Max(0.,(x/%g)*TMath::Exp(-x/(%g)))",
				   tao*tao,tao),
		   mint,maxt);
  fFunc2 = new TF1(Form("fFunc2%s",detname),
		   TString::Format("%g*TMath::Exp(-((x-%g)**2)/(%g))",
				   1./TMath::Sqrt(2*TMath::Pi()*sig),t0,sig*sig),
		   mint,maxt);
  fConvolution = new TF1Convolution(fFunc1,fFunc2);
  
  model = new TF1(Form("fSignal%s",detname),*fConvolution,mint,maxt, fConvolution->GetNpar());
}

double SPEModel::Eval(double t)
{
  return scale*model->Eval(t);
  //return model->Eval(t);
  //return 1.0;
}
