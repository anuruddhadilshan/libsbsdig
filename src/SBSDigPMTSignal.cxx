#include "SBSDigPMTSignal.h"
#include "TMath.h"

using namespace std;

//
// Class SPEModel
//
SPEModel::SPEModel()
{
  fPulseHisto = new TH1D("fPulseHisto", "", 1000, -50, 50);
}

SPEModel::SPEModel(UShort_t uniqueid, double sigma, 
		   double t0, double tmin, double tmax)
{
  TF1 fFunc("fFunc", "landaun", tmin, tmax);
  const int NbinsTotal = int(tmax-tmin)*10;// 10 bins/ns should do... since we will extrapolate after...
  fPulseHisto = new TH1D(Form("fPulseHisto_%d", uniqueid), "", NbinsTotal, tmin, tmax);
  double t_i;//, t_j;
  double ps_i;//, g_j;
  for(int i = 1; i<=NbinsTotal; i++){
    t_i = fPulseHisto->GetBinCenter(i+1);
    ps_i = fFunc.Eval(t_i);
    fPulseHisto->Fill(t_i, ps_i);
  }
}

SPEModel::~SPEModel()
{
  fPulseHisto->Delete();
}

bool SPEModel::PulseOverThr(double charge, double thr)
{
  if(fPulseHisto->GetMaximum()<thr/charge){
    return false;
  }else{
    return true;
  }
};
 
bool SPEModel::FindLeadTrailTime(double charge, double thr, double &t_lead, double &t_trail)
{
  if(!PulseOverThr(charge, thr)){
    t_lead = 1.0e38;
    t_trail = 1.0e38;
    return false;
  }else{
    double xmax = fPulseHisto->GetBinCenter(fPulseHisto->GetMaximumBin());
    if(fPulseHisto->GetBinContent(1)<thr/charge){
      t_lead = GetHistoX(thr/charge, fPulseHisto->GetBinLowEdge(1), xmax);
    }else{
      t_lead = 1.0e38;
      return false;
    }
    if(fPulseHisto->GetBinContent(fPulseHisto->GetNbinsX())<thr/charge){
      t_trail = GetHistoX(thr/charge, xmax, fPulseHisto->GetBinLowEdge(fPulseHisto->GetNbinsX()+1));
    }else{
      t_trail = 1.0e38;
      return false;
    }
    //why does this function seem to oblitarate fMCHit containers when t_trail is calculated to be 1.e38... GetHistoX?
    if(t_trail>1.e30) cout << "thr/charge" << thr/charge << " >? " <<  fPulseHisto->GetBinContent(fPulseHisto->GetNbinsX()) << " or " << fPulseHisto->GetBinContent(fPulseHisto->GetNbinsX()+1) << endl;
    return true;
  }
}

double SPEModel::GetHistoX(double y, double x1, double x2)
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
// Class PMTSignal
//
PMTSignal::PMTSignal()
  : fSumEdep(0), fNpe(0), fNpeChargeConv(1.0), fADC(0), fEventTime(0), fNADCSamps(0), fNSamps(0)
{ 
  fLeadTimes.clear();
  fTrailTimes.clear();
  fTDCs.clear();
  
  f1 = 0;
  //R = TRndmManager::GetInstance();
}

PMTSignal::PMTSignal(double npechargeconv)
  : fSumEdep(0), fNpe(0), fNpeChargeConv(npechargeconv), fADC(0), fEventTime(0), fNADCSamps(0), fNSamps(0)
{ 
  fLeadTimes.clear();
  fTrailTimes.clear();
  fTDCs.clear();
  
  f1 = 0;
}

void PMTSignal::Fill(SPEModel *model, int npe, double thr, double evttime, int signal)
{
  if(signal==0)fEventTime = evttime;
  fNpe+= npe;
  //if(model->PulseOverThr(fCharge, thr))fNpe//fADC = model->GetCharge()*model->GetADCconversion();
  
  //determine lead and trail times
  double t_lead, t_trail;
  // find the lead and trail time for *this* pulse, not the total pulse
  bool goodtime = model->FindLeadTrailTime(npe*fNpeChargeConv, thr, t_lead, t_trail);
  
  t_lead+=evttime;
  t_trail+=evttime;
  
  if(goodtime){
    //Filter here the lead and trail times
    if(fLeadTimes.size()>0){
      // Check if the lead and trail times are inside an existing lead time- trail time pair
      // we assume here that fLeadTimes and fTrailTimes are same size 
      // *if we do things correctly, that should be the case*
      // we shall keep lead-trail times pair in timing order
      // we neglect pulse overlaps ftm.
      if(fLeadTimes.size()!=fTrailTimes.size()){
	cout << " B - Warning: size of lead times container: " << fLeadTimes.size() 
	     << " != size of trail times container: " << fTrailTimes.size() << endl;
      }
      
      for(size_t i = 0; i<fLeadTimes.size(); i++){
	// possibility of the current pair straddling with others.... :/
	// treat those separately to simplify...
	// tL < tT_i-1 < tL_i < tT
	if(i>0){
	  if(t_lead < fTrailTimes.at(i-1) && fLeadTimes.at(i) < t_trail){
	    //do necessary substitutions first, then erase
	    // tL < tL_i-1 => tL *replaces* tL_i-1
	    if(t_lead < fLeadTimes.at(i-1)){
	      fLeadTimes.erase(fLeadTimes.begin()+i-1);
	      fLeadTimes.insert(fLeadTimes.begin()+i-1, t_lead);
	    }
	    // tT_i < tT => tT *replaces* tT_i
	    if(fTrailTimes.at(i) < t_trail){
	      fTrailTimes.erase(fTrailTimes.begin()+i);
	      fTrailTimes.insert(fTrailTimes.begin()+i, t_trail);
	    }
	    fLeadTimes.erase(fLeadTimes.begin()+i);
	    fTrailTimes.erase(fTrailTimes.begin()+i-1);
	    break;
	  }
	}//end if(i>0)
	// tL < tT_i < tL_i+1 < tT
	if(i<fLeadTimes.size()-1){
	  if(t_lead < fTrailTimes.at(i) && fLeadTimes.at(i+1) < t_trail){
	    //do necessary substitutions first, then erase
	    // tL < tL_i => tL *replaces* tL_i
	    if(t_lead < fLeadTimes.at(i)){
	      fLeadTimes.erase(fLeadTimes.begin()+i);
	      fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	    }
	    // tT_i+1 < tT => tT *replaces* tT_i+1
	    if(fTrailTimes.at(i+1) < t_trail){
	      fTrailTimes.erase(fTrailTimes.begin()+i+1);
	      fTrailTimes.insert(fTrailTimes.begin()+i+1, t_trail);
	    }
	    fLeadTimes.erase(fLeadTimes.begin()+i+1);
	    fTrailTimes.erase(fTrailTimes.begin()+i);
	    break;
	  }
	}//end if(i<fLeadTimes.size()-1)
	
	// if not, 6 cases to consider:
	// tL < tT < tL_i < tT_i => both inserted *before* existing pair 
	if(t_trail < fLeadTimes.at(i)){
	  fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	  fTrailTimes.insert(fTrailTimes.begin()+i, t_trail); 
	  break;
	}
	// tL < tL_i < tT < tT_i => tL *replaces* tL_i
	if(t_lead < fLeadTimes.at(i) && fLeadTimes.at(i) < t_trail && t_trail < fTrailTimes.at(i)){
	  fLeadTimes.erase(fLeadTimes.begin()+i);
	  fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	  break;
	}
	// tL_i < tL < tT < tT_i => tL *replaces* tL_i AND tT *replaces* tT_i
	if(t_lead < fLeadTimes.at(i) && fTrailTimes.at(i) < t_trail){
	  fLeadTimes.erase(fLeadTimes.begin()+i);
	  fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	  fTrailTimes.erase(fTrailTimes.begin()+i);
	  fTrailTimes.insert(fTrailTimes.begin()+i, t_trail);
	  break;
	}
	if(fLeadTimes.at(i) < t_lead && t_trail < fTrailTimes.at(i)){
	  break;
	}
	if(fLeadTimes.at(i) < t_lead && t_lead < fTrailTimes.at(i) && fTrailTimes.at(i) < t_trail){
	  fTrailTimes.erase(fTrailTimes.begin()+i);
	  fTrailTimes.insert(fTrailTimes.begin()+i, t_trail);
	}
	if(fTrailTimes.at(i) < t_lead){
	  fLeadTimes.insert(fLeadTimes.begin()+i+1, t_lead);
	  fTrailTimes.insert(fTrailTimes.begin()+i+1, t_trail);
	  break;
	}
	if(fLeadTimes.size()!=fTrailTimes.size()){
	  cout << " A - Warning: size of lead times container: " << fLeadTimes.size() 
	       << " != size of trail times container: " << fTrailTimes.size() << endl;
	}
	if(i>=fLeadTimes.size())cout << "Warning: i = " << i << " >= size of containers = " << fLeadTimes.size() << endl;
      }
    }else{
      //of course, if initial size was 0, just psuh it back
      //hopefully it will be the case most of the time
      fLeadTimes.push_back(t_lead);
      fTrailTimes.push_back(t_trail);
    }
    if(fLeadTimes.size()!=fTrailTimes.size()){
      cout << " A - Warning: size of lead times container: " << fLeadTimes.size() 
	   << " != size of trail times container: " << fTrailTimes.size() << endl;
    }
  }//end if(t_lead && t_trail<30)
}

void PMTSignal::Fill(int npe, double thr, double evttime, double sigmatime, int signal)
{
  if(signal==0)fEventTime = evttime;
  fNpe+= npe;
  //if(model->PulseOverThr(fCharge, thr))fNpe//fADC = model->GetCharge()*model->GetADCconversion();
  
  if(evttime>=fTmin+fNSamps*fSampSize)return;
  
  f1->SetParameters(npe*fNpeChargeConv, evttime, sigmatime);
  //determine lead and trail times
  double t_lead, t_trail;
  bool goodtime = false;//model->FindLeadTrailTime(npe*fNpeChargeConv, thr, t_lead, t_trail);
  if(fNSamps){
    fSamples[0]+= f1->Eval(fTmin+(0.5)*fSampSize)*fSampSize;
    for(int i = 1; i<fNSamps; i++){
      fSamples[i]+= f1->Eval(fTmin+(i+0.5)*fSampSize)*fSampSize;
      //if(i>0){
      if(fSamples[i-1]<=thr && thr<fSamples[i]){
	t_lead = fTmin+(i-0.5)*fSampSize+fSampSize*(thr-fSamples[i-1])/(fSamples[i]-fSamples[i-1]);
	goodtime = true;
      }
      if(fSamples[i-1]>=thr && thr>fSamples[i]){
	t_trail = fTmin+(i-0.5)*fSampSize+fSampSize*(thr-fSamples[i-1])/(fSamples[i]-fSamples[i-1]);
	goodtime = true;
      }
    }
  }
    
  // find the lead and trail time for *this* pulse, not the total pulse
  //t_lead+=evttime;
  //t_trail+=evttime;
  
  if(goodtime){
    //Filter here the lead and trail times
    if(fLeadTimes.size()>0){
      // Check if the lead and trail times are inside an existing lead time- trail time pair
      // we assume here that fLeadTimes and fTrailTimes are same size 
      // *if we do things correctly, that should be the case*
      // we shall keep lead-trail times pair in timing order
      // we neglect pulse overlaps ftm.
      if(fLeadTimes.size()!=fTrailTimes.size()){
	cout << " B - Warning: size of lead times container: " << fLeadTimes.size() 
	     << " != size of trail times container: " << fTrailTimes.size() << endl;
      }
      
      for(size_t i = 0; i<fLeadTimes.size(); i++){
	// possibility of the current pair straddling with others.... :/
	// treat those separately to simplify...
	// tL < tT_i-1 < tL_i < tT
	if(i>0){
	  if(t_lead < fTrailTimes.at(i-1) && fLeadTimes.at(i) < t_trail){
	    //do necessary substitutions first, then erase
	    // tL < tL_i-1 => tL *replaces* tL_i-1
	    if(t_lead < fLeadTimes.at(i-1)){
	      fLeadTimes.erase(fLeadTimes.begin()+i-1);
	      fLeadTimes.insert(fLeadTimes.begin()+i-1, t_lead);
	    }
	    // tT_i < tT => tT *replaces* tT_i
	    if(fTrailTimes.at(i) < t_trail){
	      fTrailTimes.erase(fTrailTimes.begin()+i);
	      fTrailTimes.insert(fTrailTimes.begin()+i, t_trail);
	    }
	    fLeadTimes.erase(fLeadTimes.begin()+i);
	    fTrailTimes.erase(fTrailTimes.begin()+i-1);
	    break;
	  }
	}//end if(i>0)
	// tL < tT_i < tL_i+1 < tT
	if(i<fLeadTimes.size()-1){
	  if(t_lead < fTrailTimes.at(i) && fLeadTimes.at(i+1) < t_trail){
	    //do necessary substitutions first, then erase
	    // tL < tL_i => tL *replaces* tL_i
	    if(t_lead < fLeadTimes.at(i)){
	      fLeadTimes.erase(fLeadTimes.begin()+i);
	      fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	    }
	    // tT_i+1 < tT => tT *replaces* tT_i+1
	    if(fTrailTimes.at(i+1) < t_trail){
	      fTrailTimes.erase(fTrailTimes.begin()+i+1);
	      fTrailTimes.insert(fTrailTimes.begin()+i+1, t_trail);
	    }
	    fLeadTimes.erase(fLeadTimes.begin()+i+1);
	    fTrailTimes.erase(fTrailTimes.begin()+i);
	    break;
	  }
	}//end if(i<fLeadTimes.size()-1)
	
	// if not, 6 cases to consider:
	// tL < tT < tL_i < tT_i => both inserted *before* existing pair 
	if(t_trail < fLeadTimes.at(i)){
	  fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	  fTrailTimes.insert(fTrailTimes.begin()+i, t_trail); 
	  break;
	}
	// tL < tL_i < tT < tT_i => tL *replaces* tL_i
	if(t_lead < fLeadTimes.at(i) && fLeadTimes.at(i) < t_trail && t_trail < fTrailTimes.at(i)){
	  fLeadTimes.erase(fLeadTimes.begin()+i);
	  fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	  break;
	}
	// tL_i < tL < tT < tT_i => tL *replaces* tL_i AND tT *replaces* tT_i
	if(t_lead < fLeadTimes.at(i) && fTrailTimes.at(i) < t_trail){
	  fLeadTimes.erase(fLeadTimes.begin()+i);
	  fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	  fTrailTimes.erase(fTrailTimes.begin()+i);
	  fTrailTimes.insert(fTrailTimes.begin()+i, t_trail);
	  break;
	}
	if(fLeadTimes.at(i) < t_lead && t_trail < fTrailTimes.at(i)){
	  break;
	}
	if(fLeadTimes.at(i) < t_lead && t_lead < fTrailTimes.at(i) && fTrailTimes.at(i) < t_trail){
	  fTrailTimes.erase(fTrailTimes.begin()+i);
	  fTrailTimes.insert(fTrailTimes.begin()+i, t_trail);
	}
	if(fTrailTimes.at(i) < t_lead){
	  fLeadTimes.insert(fLeadTimes.begin()+i+1, t_lead);
	  fTrailTimes.insert(fTrailTimes.begin()+i+1, t_trail);
	  break;
	}
	if(fLeadTimes.size()!=fTrailTimes.size()){
	  cout << " A - Warning: size of lead times container: " << fLeadTimes.size() 
	       << " != size of trail times container: " << fTrailTimes.size() << endl;
	}
	if(i>=fLeadTimes.size())cout << "Warning: i = " << i << " >= size of containers = " << fLeadTimes.size() << endl;
      }
    }else{
      //of course, if initial size was 0, just psuh it back
      //hopefully it will be the case most of the time
      fLeadTimes.push_back(t_lead);
      fTrailTimes.push_back(t_trail);
    }
    if(fLeadTimes.size()!=fTrailTimes.size()){
      cout << " A - Warning: size of lead times container: " << fLeadTimes.size() 
	   << " != size of trail times container: " << fTrailTimes.size() << endl;
    }
  }//end if(t_lead && t_trail<30)
}



void PMTSignal::Digitize()//TDigInfo diginfo, int chan)
{
  /*
  if(fNpe<=0){
    fADC = R->Gaus(diginfo.Pedestal(chan), diginfo.PedestalNoise(chan));
    return;
  }
    
#if DEBUG>0
  cout << "Charge (C) " << Charge() << " (fC) " << Charge()*1.0e15 << ", ADC conversion (fC/ch) " << diginfo.ADCConversion();
#endif
  
  fADC = TMath::Nint(Charge()*1.0e15/diginfo.ADCConversion()+R->Gaus(diginfo.Pedestal(chan), diginfo.PedestalNoise(chan)));
  //if ADC value bigger than number of ADC bits, ADC saturates
  if( fADC>UInt_t(TMath::Nint( TMath::Power(2, diginfo.ADCBits()) )) ){
    fADC = TMath::Nint( TMath::Power(2, diginfo.ADCBits()) );
  }
  //cout << "PMTSignal::Digitize():  " << fLeadTimes.size() << " - " << fTrailTimes.size() << endl;
  
#if DEBUG>0
  cout << " => ADC = " << fADC << endl;
#endif
  
  UInt_t tdc_value;
  
  // For the sake of going forward, we assume that the signal is the first entry of each vector
  fTDCData.time.clear();
  if(fLeadTimes.size() && fTrailTimes.size()){
    for(size_t i = 0; i<fLeadTimes.size(); i++){
#if DEBUG>0
      cout << " fLeadTimes.at(" << i << ") " << fLeadTimes.at(i) 
	   << " fTrailTimes.at(" << i << ") " << fTrailTimes.at(i) << endl;
#endif
      // trim "all" bits that are above the number of TDC bits - a couple to speed it up
      // (since TDC have a revolving clock, as far as I understand)
      // let's use an arbitrary reference time offset of 1us before the trigger
      tdc_value = TMath::Nint((fLeadTimes.at(i)+1.e3)/diginfo.TDCConversion());
      for(int j = 30; j>=diginfo.TDCBits(); j--){
	tdc_value ^= ( -0 ^ tdc_value) & ( 1 << (j) );
      }
      tdc_value ^= ( -0 ^ tdc_value) & ( 1 << (31) );
      fTDCData.time.push_back(tdc_value);
      //fTDCs.insert(fTDCs.begin()+0, TMath::Nint(fLeadTimes.at(0)*diginfo.TDCConversion()));//bug!!!!
      fTDCs.push_back(tdc_value);//they're already sorted in order, presumably
      // also mark the traling time with setting bin 31 to 1 // need to reconvert then
      tdc_value = TMath::Nint((fTrailTimes.at(i)+1.e3)/diginfo.TDCConversion());
      for(int j = 30; j>=diginfo.TDCBits(); j--){
	tdc_value ^= ( -0 ^ tdc_value) & ( 1 << (j) );
      }
      tdc_value ^= ( -1 ^ tdc_value) & ( 1 << (31) );
      fTDCs.push_back(tdc_value);
      fTDCData.time.push_back(tdc_value);
      
#if DEBUG>0
      cout << " fTDCs.at(0) " << fTDCs.at(0) << " fTDCs.at(1) " << fTDCs.at(1) << endl;
#endif
    }
  }
  */
  fSumEdep*=1.0e9;// store in eV.
}

void PMTSignal::SetSamples(double tmin, double tmax, double sampsize)
{
  fTmin = tmin;
  fADCSampSize = sampsize;
  fSampSize = sampsize*2.5;
  fNADCSamps = round((tmax-tmin)/fADCSampSize);
  fNSamps = round((tmax-tmin)/fSampSize);
  fADCSamples = new double[fNADCSamps];
  fSamples = new double[fNSamps];
  
  memset(fSamples, 0, fNSamps*sizeof(double));
  memset(fADCSamples, 0, fNADCSamps*sizeof(double));
  f1 = new TF1("f1", "landaun", tmin, tmax);
}

void PMTSignal::Clear(bool dosamples)
{
  //cout << " PMTSignal::Clear() " << endl;
  
  fSumEdep = 0;
  fNpe = 0;
  fADC = 0;
  
  fEventTime = 0;
  fLeadTimes.clear();
  fTrailTimes.clear();
  fTDCs.clear();
  
  if(dosamples){
    memset(fSamples, 0, fNSamps*sizeof(double));
    memset(fADCSamples, 0, fNSamps*sizeof(double));
  }
}



//ClassImp(SPEModel)

