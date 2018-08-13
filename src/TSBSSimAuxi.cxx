#include "TSBSSimAuxi.h"
#include "TSBSDBManager.h"

//
// Class TNPEModel
//
TSPEModel::TSPEModel(const char* detname, 
		     double tau, double sigma, 
		     double t0, double tmin, double tmax)
{
  TF1 fFunc1(Form("fFunc1%s",detname),
	     TString::Format("TMath::Max(0.,(x/%g)*TMath::Exp(-x/(%g)))", 
			     tau*tau, tau),
	     tmin ,tmax);
  TF1 fFunc2(Form("fFunc2%s",detname),
	     TString::Format("%g*TMath::Exp(-((x-%g)**2)/(%g))", 
			     1./TMath::Sqrt(2*TMath::Pi()*sigma), t0, sigma*sigma),
	     tmin, tmax);
  TF1Convolution fConvolution(&fFunc1, &fFunc2);
  
  fPulseModel = new TF1(Form("fSignal%s",detname), fConvolution, tmin, tmax, fConvolution.GetNpar());
  
}

bool TSPEModel::PulseOverThr(double charge, double thr)
{
  if(fPulseModel->GetMaximum()<thr/charge){
    return false;
  }else{
    return true;
  }
};
 
void TSPEModel::FindLeadTrailTime(double charge, double thr, double &t_lead, double &t_trail)
{
  if(!PulseOverThr(charge, thr)){
    t_lead = 1.0e38;
    t_trail = 1.0e38;
  }else{
    double xmax = fPulseModel->GetMaximumX();
    t_lead = fPulseModel->GetX(thr/charge, fPulseModel->GetXmin(), xmax);
    t_trail = fPulseModel->GetX(thr/charge, xmax, fPulseModel->GetXmax());
  }
}

//
// Class TPMTSignal
//
TPMTSignal::TPMTSignal()
  : fSumEdep(0), fNpe(0), fNpeChargeConv(1.0), fADC(0), fEventTime(0)
{ 
  fLeadTimes.clear();
  fTrailTimes.clear();
  fTDCs.clear();
}

TPMTSignal::TPMTSignal(double npechargeconv)
  : fSumEdep(0), fNpe(0), fNpeChargeConv(npechargeconv), fADC(0), fEventTime(0)
{ 
  fLeadTimes.clear();
  fTrailTimes.clear();
  fTDCs.clear();
}

void TPMTSignal::Fill(TSPEModel *model, int npe, double thr, double evttime, bool signal)
{
  if(signal)fEventTime = evttime;
  fNpe+= npe;
  //if(model->PulseOverThr(fCharge, thr))fNpe//fADC = model->GetCharge()*model->GetADCconversion();
  
  //determine lead and trail times
  double t_lead, t_trail;
  // find the lead and trail time for *this* pulse, not the total pulse
  model->FindLeadTrailTime(npe*fNpeChargeConv, thr, t_lead, t_trail);
  fLeadTimes.push_back(evttime+t_lead);
  fTrailTimes.push_back(evttime+t_trail);
  
}

void TPMTSignal::Digitize(DigInfo diginfo)
{
  if(fNpe<=0)
    return;
  
  fADC = fNpe*fNpeChargeConv*diginfo.fADCconversion;
  //for(
  
  
  /*
  int braw = 0;
  double max = 0;
  sum = 0;
  for(int bs = 0; bs < nbins; bs++) {
    max = 0;
    for(int br = 0; br < dnraw; br++) {
      if(samples_raw[br+braw] > max)
        max = samples_raw[br+braw];
    }
    if(max>2)
      max = 2;
    //samples[bs] =int((max/2.);// *4095);
    samples[bs] =int((max/2.)*4095);
    //samples[bs] = samples[bs] > 4095 ? 4095 : samples[bs];
    braw += dnraw;
    sum += samples[bs];
  }

  // Also digitize the sumedep
  sumedep *= 1e9; // To store in eV
  */
}

void TPMTSignal::Clear()
{
  fSumEdep = 0;
  fNpe = 0;
  fNpeChargeConv = 0;
  fADC = 0;
  
  fEventTime = 0;
  fLeadTimes.clear();
  fTrailTimes.clear();
  fTDCs.clear();
}
 
TPMTSignal::~TPMTSignal()
{
  Clear();
}


//
// Class TNPEModel
//
TNPEModel::TNPEModel(DigInfo diginfo, const char* detname, int npe)
  : fDigInfo(diginfo), fNpe(npe)
{
  fScale = fDigInfo.fROImpedance*qe/spe_unit;
  
  if(fDigInfo.fGain.size()==1){
    fScale*= fDigInfo.fGain[0];
  }
  
  double mint = -fDigInfo.fGateWidth/2.0;
  double maxt = +fDigInfo.fGateWidth/2.0;
  // test values
  double tau = fDigInfo.fSPEtau;
  double sig = fDigInfo.fSPEsig;
  double t0 = fDigInfo.fSPEtransittime-fDigInfo.fTriggerOffset;
  fStartTime = t0;
  
  TF1 fFunc1(Form("fFunc1%s",detname),
	     TString::Format("TMath::Max(0.,(x/%g)*TMath::Exp(-x/(%g)))",
			     tau*tau,tau),
	     mint,maxt);
  TF1 fFunc2(Form("fFunc2%s",detname),
	     TString::Format("%g*TMath::Exp(-((x-%g)**2)/(%g))",
			     1./TMath::Sqrt(2*TMath::Pi()*sig),t0,sig*sig),
	     mint,maxt);
  TF1Convolution fConvolution(&fFunc1,&fFunc2);
  
  fModel = new TF1(Form("fSignal%s",detname),fConvolution,mint,maxt, fConvolution.GetNpar());
}

// return the total charge: comes in handy to evaluate ADC value quickly.
double TNPEModel::GetCharge(int chan)
{
  double totalcharge = fNpe*fScale;
  if(fDigInfo.fGain.size()>1){//if not, fScale already includes the gain
    if(fDigInfo.fGain.size()<=chan){
      cout << "warning: requested channel number " << chan << "larger than number of channel size " << fDigInfo.fGain.size() << " check your code ! " << endl;
      exit(-1);
    }
    totalcharge*= fDigInfo.fGain[chan];
  }
  
  return totalcharge;
}

//
double TNPEModel::Eval(int chan, double t)
{
  double totalscale = fNpe*fScale;
  if(fDigInfo.fGain.size()>1){//if not, fScale already includes the gain
    if(fDigInfo.fGain.size()<=chan){
      cout << "warning: requested channel number " << chan << "larger than number of channel size " << fDigInfo.fGain.size() << " check your code ! " << endl;
      exit(-1);
    }
    totalscale*= fDigInfo.fGain[chan];
  }
  
  return totalscale*fModel->Eval(t);
}

// 
void TNPEModel::FindLeadTrailTime(int chan, double &t_lead, double &t_trail)
{
  double thresh = fDigInfo.fThreshold[0];
  if(fDigInfo.fThreshold.size()>1){
    if(fDigInfo.fThreshold.size()<=chan){
      cout << "warning: requested channel number " << chan << "larger than number of channel size " << fDigInfo.fThreshold.size() << " check your code ! " << endl;
      exit(-1);
    }
    thresh = fDigInfo.fThreshold[chan];
  }
  double totalscale = fNpe*fScale;
  if(fDigInfo.fGain.size()>1){
    totalscale*= fDigInfo.fGain[chan];
  }
  //Since we cannot seem to scale the pulse, we'll scale the threshold
  thresh = thresh/totalscale;
  
  if(PulseOverThr()){
    t_lead = 1.0e38;
    t_trail = 1.0e38;
  }else{
    double xmax = fModel->GetMaximumX();
    t_lead = fModel->GetX(thresh, -fDigInfo.fGateWidth, xmax);
    t_trail = fModel->GetX(thresh, xmax, +fDigInfo.fGateWidth);
    
  }
}

// 
bool TNPEModel::PulseOverThr(int chan)
{
  double thresh = fDigInfo.fThreshold[0];
  if(fDigInfo.fThreshold.size()>1){
    if(fDigInfo.fThreshold.size()<=chan){
      cout << "warning: requested channel number " << chan << "larger than number of channel size " << fDigInfo.fThreshold.size() << " check your code ! " << endl;
      exit(-1);
    }
    thresh = fDigInfo.fThreshold[chan];
  }
  double totalscale = fNpe*fScale;
  if(fDigInfo.fGain.size()>1){
    totalscale*= fDigInfo.fGain[chan];
  }
  //Since we cannot seem to scale the pulse, we'll scale the threshold
  thresh = thresh/totalscale;
  
  if(fModel->GetMaximum(-fDigInfo.fGateWidth, +fDigInfo.fGateWidth)<thresh){
    return false;
  }else{
    return true;
  }
};
