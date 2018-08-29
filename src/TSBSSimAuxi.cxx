#include "TSBSSimAuxi.h"
#include "TSBSDBManager.h"
#define DEBUG 0

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
  TF1 fFunc2(Form("fFunc2%s",detname), "gaus", 
  // 	     TString::Format("%g*TMath::Exp(-((x-%g)**2)/(%g))", 
  // 			     1./TMath::Sqrt(2*TMath::Pi()*sigma), t0, sigma*sigma),
   	     tmin, tmax);
  fFunc2.SetParameters(1.0e0/10.0/(sqrt(2*TMath::Pi())*sigma), 0, sigma);//do it this way for thaat purpose
  
  //TF1Convolution is deemed too premature to be used - and indeed pretty slow
  //TF1Convolution fConvolution(&fFunc1, &fFunc2);
  //fPulseModel = new TF1(Form("fSignal%s",detname), fConvolution, tmin, tmax, fConvolution.GetNpar());
  
  const int NbinsTotal = int(tmax-tmin)*10;// 10 bins/ns should do... since we will extrapolate after...
  fPulseHisto = new TH1D(Form("fPulseHisto%s",detname), "", NbinsTotal, tmin, tmax);
  double t_i, t_j;
  double ps_i, g_j;
  for(int i = 1; i<=NbinsTotal; i++){
    t_i = fPulseHisto->GetBinCenter(i+1);
    ps_i = fFunc1.Eval(t_i);
    if(sigma>0){
      fFunc2.SetParameter(1, t_i);
      for(int j = 1; j<=NbinsTotal; j++){
	t_j = fPulseHisto->GetBinCenter(j+1);
	g_j = fFunc2.Eval(t_j);
	fPulseHisto->Fill(t_j, ps_i*g_j);
      }
    }else{
      fPulseHisto->Fill(t_i, ps_i);
    }
  }
  
  cout << endl<< detname << " pulse histo built" << endl;
#if DEBUG>2
  for(int i = 0; i<NbinsTotal; i++){
    if(fPulseHisto->GetBinContent(i)>0)cout << fPulseHisto->GetBinContent(i) << " ";
  }
  cout << endl;
#endif
  
}

bool TSPEModel::PulseOverThr(double charge, double thr)
{
#if DEBUG>0
  cout << "unnormalized pulse max (ns -1) " << fPulseModel->GetMaximum() << ", threshold (C/ns) " << thr << ", charge (C) " << charge << endl;
#endif
  //if(fPulseModel->GetMaximum()<thr/charge){
  if(fPulseHisto->GetMaximum()<thr/charge){
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
    /*
    double xmax = fPulseModel->GetMaximumX();
    t_lead = fPulseModel->GetX(thr/charge, fPulseModel->GetXmin(), xmax);
    t_trail = fPulseModel->GetX(thr/charge, xmax, fPulseModel->GetXmax());
    */
    double xmax = fPulseHisto->GetBinCenter(fPulseHisto->GetMaximumBin());
    t_lead = GetHistoX(thr/charge, fPulseHisto->GetBinLowEdge(1), xmax);
    t_trail = GetHistoX(thr/charge, xmax, fPulseHisto->GetBinLowEdge(fPulseHisto->GetNbinsX()+1));
  }
}

double TSPEModel::GetHistoX(double y, double x1, double x2)
{
  double splineslope;
  for(int k = fPulseHisto->FindBin(x1); k<=fPulseHisto->FindBin(x2); k++){
    if(  ( (fPulseHisto->GetBinContent(k+1)-y)*(fPulseHisto->GetBinContent(k)-y) )<0 ){
      // threshold crossed if diff(TH1::GetBinContent-y) changes sign
      splineslope = (fPulseHisto->GetBinContent(k+1)-fPulseHisto->GetBinContent(k))/(fPulseHisto->GetBinCenter(k+1)-fPulseHisto->GetBinCenter(k));
      
      return fPulseHisto->GetBinCenter(k)+(y-fPulseHisto->GetBinContent(k))/splineslope;
    }
  }
  return 1.0e38;
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
  
  //cout << "TPMTSignal::Fill : fNpeChargeConv = " << fNpeChargeConv << endl;
  
  //determine lead and trail times
  double t_lead, t_trail;
  // find the lead and trail time for *this* pulse, not the total pulse
  model->FindLeadTrailTime(npe*fNpeChargeConv, thr, t_lead, t_trail);
  if(t_lead<1e30 && t_trail<1e30){
    fLeadTimes.push_back(evttime+t_lead);
    fTrailTimes.push_back(evttime+t_trail);
  }
}

void TPMTSignal::Digitize(TDigInfo diginfo, int chan)
{
  if(fNpe<=0)
    return;
  
#if DEBUG>0
  cout << "Charge (C) " << Charge() << " (fC) " << Charge()*1.0e15 << ", ADC conversion (fC/ch) " << diginfo.ADCConversion();
#endif
  
  fADC = TMath::Nint(Charge()*1.0e15/diginfo.ADCConversion()+diginfo.GenPedestal(chan));
  //if ADC value bigger than number of ADC bits, ADC saturates
  if( fADC>TMath::Nint( TMath::Power(2, diginfo.ADCBits()) ) ){
    fADC = TMath::Nint( TMath::Power(2, diginfo.ADCBits()) );
  }
  //cout << "TPMTSignal::Digitize():  " << fLeadTimes.size() << " - " << fTrailTimes.size() << endl;
  
#if DEBUG>0
  cout << " => ADC = " << fADC << endl;
#endif
  
  UInt_t tdc_value;
  
  // For the sake of going forward, we assume that the signal is the first entry of each vector
  if(fLeadTimes.size() && fTrailTimes.size()){
#if DEBUG>0
    cout << " fLeadTimes.at(0) " << fLeadTimes.at(0) << " fTrailTimes.at(0) " << fTrailTimes.at(0) << endl;
#endif
    
    // trim "all" bits that are above the number of TDC bits - a couple to speed it up
    // (since TDC have a revolving clock, as far as I understand)
    tdc_value = TMath::Nint(fLeadTimes.at(0)/diginfo.TDCConversion());
    for(int i = diginfo.TDCBits()+2; i>=diginfo.TDCBits(); i--){
      tdc_value ^= ( -0 ^ tdc_value) & ( 1 << (i) );
    }
    tdc_value ^= ( -0 ^ tdc_value) & ( 1 << (31) );
    fTDCs.insert(fTDCs.begin()+0, TMath::Nint(fLeadTimes.at(0)*diginfo.TDCConversion()));
    // also mark the traling time with setting bin 31 to 1
    tdc_value = TMath::Nint(fTrailTimes.at(0)/diginfo.TDCConversion());
    for(int i = diginfo.TDCBits()+2; i>=diginfo.TDCBits(); i--){
      tdc_value ^= ( -0 ^ tdc_value) & ( 1 << (i) );
    }
    tdc_value ^= ( -1 ^ tdc_value) & ( 1 << (31) );
    fTDCs.insert(fTDCs.begin()+1, tdc_value);
    
#if DEBUG>0
    cout << " fTDCs.at(0) " << fTDCs.at(0) << " fTDCs.at(1) " << fTDCs.at(1) << endl;
#endif
  }
  
  /*
  // TDCs: select only lead and trail times not between a lead and a trail time.
  // too complicated for the moment. 
  int minsize = min(fLeadTimes.size(), fTrailTimes.size());
  UInt_t LeadTDC;
  UInt_t TrailTDC;
  for(int i = 1; i<minsize; i++){
    LeadTimeBoxed = false;
    TrailTimeBoxed = false;
    LeadTDC = fLeadTimes.at(i)*diginfo.fTDCconversion;
    TrailTDC = fTrailTimes.at(i)*diginfo.fTDCconversion;
    
    for(j = 0; j<fTDCs.size(); j+=2){
      //if(fTDCs.at(j)<=LeadTDC && LeadTDC<=fTDCs.at(j+1))LeadTimeBoxed = true;
      //if(fTDCs.at(j)<=TrailTDC && TrailTDC<=fTDCs.at(j+1))TrailTimeBoxed = true;
      
      if(LeadTDC<fTDCs.at(j)){//current leading time before "recorded" TDC leading time
	if(TrailTDC>=fTDCs.at(j)){//current trailing time after: replace leading time with current
	  fTDCs.erase(fTDCs.begin()+j);
	  fTDCs.insert(fTDCs.begin()+j, LeadTDC);
	}else{// current trailing time before: insert a new "pair"
	  fTDCs.insert(fTDCs.begin()+j, LeadTDC);
	  fTDCs.insert(fTDCs.begin()+j+1, TrailTDC);
	}
	break;
      }
      if(fTDCs.at(j)<=LeadTDC && LeadTDC<=fTDCs.at(j+1)){
	if(fTDCs.at(j)<=TrailTDC && TrailTDC<=fTDCs.at(j+1)){
	  break;
	}else{
	  fTDCs.erase(fTDCs.begin()+j+1);
	  fTDCs.insert(fTDCs.begin()+j+1, TrailTDC);
	}
      }
      //   if(LeadTDC<fTDCs.at(j)){
    // 	fTDCs.insert(fTDCs.begin()+j-1);
    // 	if(TrailTDC>=fTDCs.at(j)){
	  
    // 	}
    //   }
    // 	&& TrailTDC>=fTDCs.at(j)){
    // 	fTDCs.erase(fTDCs.begin()+j);
    // 	fTDCs.insert(fTDCs.begin()+j);
    //   }
      
      // 
    }
    //fTDCs.push_back(TMath::Nint(fLeadTimes.at(i)*diginfo.fTDCconversion));
  }
  // for(int i = 0; i<fLeadTimes.size(); i++){
  //   fTDCs.push_back(TMath::Nint(fLeadTimes.at(i)*diginfo.fTDCconversion));
  // }
  */
  
  fSumEdep*=1.0e9;// store in eV.
}

void TPMTSignal::Clear()
{
  //cout << " TPMTSignal::Clear() " << endl;
  
  fSumEdep = 0;
  fNpe = 0;
  fADC = 0;
  
  fEventTime = 0;
  fLeadTimes.clear();
  fTrailTimes.clear();
  fTDCs.clear();
}

//
// Class TDigSlot
//
TDigSlot::TDigSlot() : fCrate(-1), fSlot(-1), fNchan(-1), fChanLo(-1),
  fChanHi(-1)
{
};

TDigSlot::TDigSlot(Int_t crate, Int_t slot, Int_t lo,
    Int_t hi) : fCrate(crate), fSlot(slot), fChanLo(lo), fChanHi(hi)
{
  fNchan = fChanHi - fChanLo;
}

TDigSlot::~TDigSlot()
{
}

Int_t TDigSlot::GetChanNumber(Int_t ch)
{
  Int_t lch = fChanLo + ch;
  if (lch > fChanHi)
    return -1;

  return lch;
}

//
// Class TDigDetMap
//
/*
Int_t TDigDetMap::Fill(std::vector<Int_t> vals)
{
  // TODO: Check if some of the data does not repeat.
  // As a first try, just pre-fill all data the user provides
  for(size_t k = 0; k < vals.size(); k+=4) {
    TDigSlot slot(vals[k],vals[k+1],vals[k+2],vals[k+3]);
    fSlots.push_back(slot);
  }

  return 0;
}
*/


//
// Class TDetInfo
//
TDetInfo::TDetInfo()
{
  fModSlots.clear();
  fNmodules.clear();
  fGeoInfo.clear();
}

TDetInfo::TDetInfo(const std::string detname)
{
  fDetName = detname;
  fNmodules.clear();
  fGeoInfo.clear();
}

TDetInfo::~TDetInfo()
{
  fNmodules.clear();
  fGeoInfo.clear();
}

Int_t TDetInfo::AddSlot(Int_t crate, Int_t slot, Int_t lo, Int_t hi)
{
  TDigSlot modslot(crate,slot,lo,hi);
  fModSlots.push_back(modslot);
  return modslot.GetNchan();
}


TDigChannelInfo TDetInfo::FindLogicalChannelSlot(Int_t lch)
{
  TDigChannelInfo info;
  info.ch = -1;
  info.slot = -1;
  info.crate = -1;
  // If we have a detector map, then use that
  if(!fDetMap.empty()) {
    std::map<int,std::pair<int,int> >::iterator it = fDetMap.find(lch);
    if(it != fDetMap.end() ) {
      TDigSlot &sl = fModSlots[it->second.first];
      info.ch = sl.GetChanNumber(it->second.second);
      info.slot = sl.GetSlot();
      info.crate = sl.GetCrate();
      return info;
    }
  } else if(! fModSlots.empty()) {
    // Otherwise, loop through all the modules and find the one we want
    for(std::vector<TDigSlot>::iterator it = fModSlots.begin();
        it != fModSlots.end(); it++) {
      if ( (*it).GetNchan() < lch ) {
        info.ch = (*it).GetChanNumber(lch);
        info.slot = (*it).GetSlot();
        info.crate = (*it).GetCrate();
        return info;
      }
      lch -= (*it).GetNchan();
    }
  } else {
    // No map of any kind, so then come up with the channel
    // based on fFirstSlot and fFirstCrate
    info.ch = lch%fChanPerSlot;
    info.slot = ((lch-info.ch)/fChanPerSlot)%fSlotPerCrate+fFirstSlot;
    info.crate = (lch-info.slot*fChanPerSlot-info.ch)/fSlotPerCrate+fFirstCrate;
  }

  // Well, if we got to here, then clearly we didn't find it, so instead
  // make it up based on the firstSlot and lastSlot
  return info;

}


void TDetInfo::LoadChannelMap(std::vector<int> chanmap)
{
  // Assume DBManager has checked it for proper size and proceed blindly
  // accepting the format.
  int lch = 0;
  int nmods = fModSlots.size();
  int nch = 0;
  for(int i = 0; i < nmods; i++) {
    nch = fModSlots[i].GetNchan();
    for(int k = 0; k < nch; k++) {
      fDetMap[lch++] = std::pair<int,int>(i,k);
    }
  }
}

//
// Class TDigInfo
//
TDigInfo::TDigInfo()
{
  fRN = new TRandom3(0);
  fGain.clear();
  fPedestal.clear();
  fPedNoise.clear();
  fThreshold.clear();
}

TDigInfo::~TDigInfo()
{
  fGain.clear();
  fPedestal.clear();
  fPedNoise.clear();
  fThreshold.clear();
}

double TDigInfo::Gain(uint chan)
{
  if(fGain.size()>1){
    if(fGain.size()<=chan){
      printf("warning: requested channel number %ud larger than channel size %ld. Check code where this method (DigInfo.Gain(int)) is employed and /or database!\n", 
	     chan, fGain.size());
      exit(-1);
    }
    return fGain.at(chan);
  }else{
    return fGain.at(0);
  }
};

double TDigInfo::Pedestal(uint chan)
{
  if(fPedestal.size()>1){
    if(fPedestal.size()<=chan){
      printf("warning: requested channel number %ud larger than channel size %ld. Check code where this method (DigInfo.Pedestal(int)) is employed and /or database!\n", 
	     chan, fPedestal.size());
      exit(-1);
    }
    return fPedestal.at(chan);
  }else{
    return fPedestal.at(0);
  }
};

double TDigInfo::PedestalNoise(uint chan)
{
  if(fPedNoise.size()>1){
    if(fPedNoise.size()<=chan){
      printf("warning: requested channel number %ud larger than channel size %ld. Check code where this method (DigInfo.PedestalNoise(int)) is employed and /or database!\n", 
	     chan, fPedNoise.size());
      exit(-1);
    }
    return fPedNoise.at(chan);
  }else{
    return fPedNoise.at(0);
  }
};

double TDigInfo::Threshold(uint chan)//returns threshold in C/ns !!!
{
  if(fThreshold.size()>1){
    if(fThreshold.size()<=chan){
      printf("warning: requested channel number %ud larger than channel size %ld. Check code where this method (DigInfo.ThresholdNoise(int)) is employed and /or database!\n", 
	     chan, fThreshold.size());
      exit(-1);
    }
    //Thr(V)/Omega = Thr(A); Thr(A)*unit = Thr(C/ns)
    return fThreshold.at(chan)*spe_unit/fROimpedance;
  }else{
    return fThreshold.at(0)*spe_unit/fROimpedance;
  }
};






/*
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
*/
