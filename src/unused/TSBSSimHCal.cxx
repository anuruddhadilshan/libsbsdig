#include "TSBSSimHCal.h"
#include <iostream>
#include <TSBSSimData.h>
//#include "sbs_types.h"
//#include <TF1.h>
//#include <TF1Convolution.h>
//#include <TTree.h>
//#include <TFile.h>
#include <TSBSSimEvent.h>
#include "TSBSDBManager.h"
#include "SBSSimDataEncoder.h"

#define HCAL_TDC_THRESH 0.1

TSBSSimHCal::TSBSSimHCal(const char* name, short id)
{
  fName = name;
  SetUniqueDetID(id);
  Init();
  fHasFADC = false;
  if(fEncoderADC && fEncoderADC->IsFADC()) {
    fHasFADC = true;
  }
}

TSBSSimHCal::~TSBSSimHCal()
{
}

void TSBSSimHCal::Init()
{
  TSBSSimDetector::Init();
  if(fDebug>=1)
    cout << "HCal detector with UniqueDetID = " << UniqueDetID() << ": TSBSSimHCal::Init() " << endl;
  
  fDetInfo = fDBmanager->GetDetInfo(fName.Data());
  
  double tau = fDetInfo.DigInfo().SPE_Tau();
  double sigma = fDetInfo.DigInfo().SPE_Sigma();
  double tmin = -fDetInfo.DigInfo().GateWidth()/2.0;
  double tmax = +fDetInfo.DigInfo().GateWidth()/2.0;
  double t0 = 0.0;//+fDetInfo.DigInfo().SPE_TransitTime()-fDetInfo.DigInfo().TriggerOffset();
  
  fSPE = new TSPEModel(fName.Data(), tau, sigma, t0, tmin, tmax);

  /*
  fSignals.resize(fDetInfo.NChan());
  for(size_t i_ch = 0; i_ch<fDetInfo.NChan(); i_ch++){
    fSignals[i_ch].Clear();
    fSignals[i_ch].SetNpeChargeConv(fDetInfo.DigInfo().NpeChargeConv(i_ch));
  }
  */
  fSignals.resize(fDetInfo.NChan());
  for(uint i = 0; i<fDetInfo.NChan(); i++){
    fSignals[i].mint = tmin;
    fSignals[i].maxt = tmax;
  }
}


void TSBSSimHCal::LoadEventData(const std::vector<g4sbshitdata*> &evbuffer)
{
  Clear();
  LoadAccumulateData(evbuffer);
}


void TSBSSimHCal::LoadAccumulateData(const std::vector<g4sbshitdata*> &evbuffer)
{
  //Double_t mint = 1e9;
  //bool signal = false;
  int chan = 0;
  int type = 0;
  double data = 0;
  for(std::vector<g4sbshitdata*>::const_iterator it = evbuffer.begin(); it!= evbuffer.end(); ++it ) {
    g4sbshitdata* ev = (*it);
    // Only get detector data for HCAL
    if(ev->GetDetUniqueID() == UniqueDetID()) {
      fSignals[chan].mc_source = ev->GetData(0)==0;
      chan  = ev->GetData(1);
      type = ev->GetData(2);
      data = ev->GetData(3);
      if(type == 0) {
        //std::cout << "Filling data for chan: " << ev->GetData(0) << ", t=" << 
        // ev->GetData(1) - 60. << std::endl;
        //if(ev->GetData(1)<mint)
        //  mint = ev->GetData(1);
        // Data is time, so use info from the configuration
        data += fDetInfo.DigInfo().SPE_TransitTime()-fDetInfo.DigInfo().TriggerOffset() + fTimeZero;//+fDetInfo.DigInfo().TriggerJitter()
        //pulsenorm = fDetInfo.DigInfo().Gain(chan)*fDetInfo.DigInfo().ROImpedance()*qe/spe_unit;
        //fSignals[chan].Fill(fSPE, data-75.);
	//cout << "spe time "<< data << endl;
        fSignals[chan].Fill(data);
        //fSignals[chan].Fill(fSPE, pulsenorm,data);
      } else if (type == 1) { // sumedep data
        fSignals[chan].sumedep = data;
      }
    }
  }
  //std::cout << "Mint = " << mint << std::endl;
}

void TSBSSimHCal::Signal::Digitize(TSPEModel *model, double pulsenorm,
    double toffset, double max_val)
{
  if(npe <= 0)
    return;

  // First, make the "scope" picture
  double t = mint;
  for(int bin = 0; bin < nbins_times; bin++) {
    if(times_histo[bin]>0) {
      FillNPE(model,pulsenorm*times_histo[bin],t,toffset);
    }
    t+=dx_raw_time;
  }

  // Now digitize
  int braw = 0;
  double max = 0;
  sum = 0;
  for(int bs = 0; bs < nbins; bs++) {
    max = 0;
    for(int br = 0; br < dnraw; br++) {
      if(samples_raw[br+braw] > max)
        max = samples_raw[br+braw];
      if( !met_tdc_thresh) {
        tdc_time += dx_raw;
        if(samples_raw[br+braw] >= HCAL_TDC_THRESH ) {
          met_tdc_thresh = true;
        }
      }
    }
    fadc.samples[bs] =int((max/2.)*max_val);
    if(fadc.samples[bs]>max_val)
      fadc.samples[bs]=max_val;
    braw += dnraw;
    sum += fadc.samples[bs];
  }

  // Also digitize the sumedep
  sumedep *= 1e9; // To store in eV

  // If we have TDC threshold met, let's set the time in a format
  // suitable for the F1 TDC
  // The F1 TDC in high resolution mode has a resolution of 60 ps LSB
  // and a range of 3.9 us (16 bits).
  // TODO: I should also stop hard coding this! Use Eric's methods
  // in TPMTSignal instead!
  if(met_tdc_thresh) {
    tdc_time -= mint; // Make sure tdc_time is always positive
    tdc.time.push_back(int((tdc_time/3.9e3)*65535));
    //tdc_time = int((tdc_time/3.9e3)*65535);
  }

}

void TSBSSimHCal::Digitize(TSBSSimEvent &event)
{
  if(fDebug>=2)cout << "TSBSSimHCal::Digitize() : Unique Det ID " << UniqueDetID() << " signal size = " << fSignals.size() << endl;
  
  bool any_events = false;
  double pulsenorm = 0;
  double max_val = pow(2,fDetInfo.DigInfo().ADCBits());
  //TSBSSimEvent::DetectorData data;
  std::vector<uint32_t> data_mod;
  std::vector<uint32_t> data;
  std::vector<int> data_dec;
  UInt_t tdcval;
  int mult = 0;
  Int_t ADCsum = 0;
  for(size_t m = 0; m < fSignals.size(); m++) {
    data.clear();
    data_dec.clear();
    data_mod.clear();
    if(fDebug>=3 && fSignals[m].npe)cout << fSignals[m].npe << endl;
    if(fSignals[m].npe > 0) {
      pulsenorm = fDetInfo.DigInfo().Gain(m)*fDetInfo.DigInfo().ROImpedance()
        *qe/spe_unit;
      if(fDebug>=4)cout << pulsenorm << endl;
      fSignals[m].Digitize(fSPE,pulsenorm,0.0,max_val);
      any_events = true;
      //data.fDetID = UniqueDetID();
      //data.fChannel = m;
      if(fDebug>=4)cout << "signal " << m << " digitization done" << endl;
      
      event.fSimHitMCOutData[fDetInfo.DetFullName()].fNSimHits++;
      event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimSource.push_back(fSignals[m].mc_source);
      event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimTRID.push_back(fSignals[m].trid);
      event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimPID.push_back(fSignals[m].pid);
      event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimChannel.push_back(Short_t(m));
      event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimEdep.push_back(fSignals[m].sumedep);
      event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimNpe.push_back(fSignals[m].npe);
      event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimTime.push_back(fSignals[m].tdc_time);
      
      if(fDebug>=4)cout << " encode TDC ? " << fEncoderTDC << endl; 
      if(fEncoderTDC){
	if(fSignals[m].tdc.time.size()){
	  if(fDebug>=4)cout << "tdc size " << fSignals[m].tdc.time.size() << endl;
	  bool tlfill = false;
	  bool ttfill = false;
	  for(uint i = 0; i<fSignals[m].tdc.time.size(); i++){
	    if(fSignals[m].tdc.getEdge(i)==0){
	      event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimLeadTime.push_back(fSignals[m].tdc.getTime(i));
	      tlfill = true;
	    }else{
	      event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimTrailTime.push_back(fSignals[m].tdc.getTime(i));
	      ttfill = true;
	    }
	  }
	  if(!tlfill)event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimLeadTime.push_back(-1000000);
	  if(!ttfill)event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimTrailTime.push_back(-1000000);
	}else{
	  event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimLeadTime.push_back(-1000000);
	  event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimTrailTime.push_back(-1000000);
	}
      }//end if(encoder)
      
      if(fDebug>=4)cout << " encode ADC " << endl; 
      mult = 0;
      fEncoderADC->EncodeFADC(fSignals[m].fadc,fEncBuffer,
          fNEncBufferWords);
      CopyEncodedData(fEncoderADC,mult++,data_mod);//.fData);

      if(fDebug>=4)cout << GetName() << " ADC size " << data_mod.size() << endl;

      ADCsum = 0;
      data_dec.clear();//data_dec.push_back(0);
      data.clear();//data.push_back(0);
      //First, the header!
      event.fSimDigSampOutData[fDetInfo.DetFullName()].fNHits++;
      event.fSimDigSampOutData[fDetInfo.DetFullName()].fChannel.push_back(-1000);
      event.fSimDigSampOutData[fDetInfo.DetFullName()].fDataWord.push_back(data_mod[0]);
      event.fSimDigSampOutData[fDetInfo.DetFullName()].fADC.push_back(-1000000);
      event.fSimDigSampOutData[fDetInfo.DetFullName()].fNsamps.push_back(0);
      event.fSimDigSampOutData[fDetInfo.DetFullName()].fADC_samps.push_back(data_dec);
      event.fSimDigSampOutData[fDetInfo.DetFullName()].fDataWord_samps.push_back(data);
      if(fEncoderTDC){
	event.fSimDigSampOutData[fDetInfo.DetFullName()].fTDC_L.push_back(-1000000);
	event.fSimDigSampOutData[fDetInfo.DetFullName()].fTDC_T.push_back(-1000000);
      }
      if(fDebug>=4)cout <<"ADC header: "<< data_mod[0] << endl;
      
      data.clear();
      data_dec.clear();
      for(uint i = 1; i<data_mod.size(); i++){
	if(fDebug>=4)cout << i << "/" << data_mod.at(i) << endl;
	data.push_back(data_mod[i]);
      }
      for(uint i = 0; i<fSignals[m].fadc.samples.size(); i++){
	if(fDebug>=4)cout << i << "/" << fSignals[m].fadc.samples.at(i) << endl;
	ADCsum+= fSignals[m].fadc.samples.at(i)-fDetInfo.DigInfo().Pedestal(m);
	data_dec.push_back(fSignals[m].fadc.samples.at(i)-fDetInfo.DigInfo().Pedestal(m));
      }
      

      event.fSimDigSampOutData[fDetInfo.DetFullName()].fNHits++;
      event.fSimDigSampOutData[fDetInfo.DetFullName()].fChannel.push_back(Short_t(m));
      //this is a bit of abuse: I use the dataword in this case to store the number of words
      event.fSimDigSampOutData[fDetInfo.DetFullName()].fDataWord.push_back(data.size());
      event.fSimDigSampOutData[fDetInfo.DetFullName()].fADC.push_back(ADCsum);
      event.fSimDigSampOutData[fDetInfo.DetFullName()].fNsamps.push_back(fSignals[m].fadc.samples.size());
      event.fSimDigSampOutData[fDetInfo.DetFullName()].fADC_samps.push_back(data_dec);
      event.fSimDigSampOutData[fDetInfo.DetFullName()].fDataWord_samps.push_back(data);
      if(fEncoderTDC){
	event.fSimDigSampOutData[fDetInfo.DetFullName()].fTDC_L.push_back(-1000000);
	event.fSimDigSampOutData[fDetInfo.DetFullName()].fTDC_T.push_back(-1000000);
      }
      data.clear();//data.push_back(0);
      data_dec.clear();//data_dec.push_back(0);
      data_mod.clear();
      
      // Now add the TDC if the threshold was met
      if(fSignals[m].met_tdc_thresh && 
	 fEncoderTDC->EncodeTDC(fSignals[m].tdc,fEncBuffer,fNEncBufferWords)){
        CopyEncodedData(fEncoderTDC,mult++,data_mod);
	if(fDebug>=4)cout << GetName() << " TDC size " << data.size() << endl;
 
	for(uint i = 0; i<data_mod.size(); i++){
	  if(fDebug>=4)cout << i << "/" << data_mod.at(i) << endl;
	  event.fSimDigSampOutData[fDetInfo.DetFullName()].fNHits++;
	  event.fSimDigSampOutData[fDetInfo.DetFullName()].fDataWord.push_back(data_mod[i]);
	  event.fSimDigSampOutData[fDetInfo.DetFullName()].fADC.push_back(-1000000);
	  event.fSimDigSampOutData[fDetInfo.DetFullName()].fNsamps.push_back(0);
	  event.fSimDigSampOutData[fDetInfo.DetFullName()].fADC_samps.push_back(data_dec);
	  event.fSimDigSampOutData[fDetInfo.DetFullName()].fDataWord_samps.push_back(data);
	  if(i==0){//header
	    event.fSimDigSampOutData[fDetInfo.DetFullName()].fChannel.push_back(-1000);
	    event.fSimDigSampOutData[fDetInfo.DetFullName()].fTDC_L.push_back(-1000000);
	    event.fSimDigSampOutData[fDetInfo.DetFullName()].fTDC_T.push_back(-1000000);
	  }else{
	    event.fSimDigSampOutData[fDetInfo.DetFullName()].fChannel.push_back(Short_t(m)+fDetInfo.NChan());//Shall we store here the "logical" chan ????
	    if( fSignals[m].tdc.getTime(i-1) & ( 1 << (31) ) ){
	      tdcval = fSignals[m].tdc.getTime(i-1);
	      tdcval ^= ( -0 ^ tdcval) & ( 1 << (31) );
	      event.fSimDigSampOutData[fDetInfo.DetFullName()].fTDC_L.push_back(-1000000);
	      event.fSimDigSampOutData[fDetInfo.DetFullName()].fTDC_T.push_back(tdcval-1.e3/fDetInfo.DigInfo().TDCConversion());
	    }else{
	      event.fSimDigSampOutData[fDetInfo.DetFullName()].fTDC_L.push_back(fSignals[m].tdc.getTime(i-1)-1.e3/fDetInfo.DigInfo().TDCConversion());
	      event.fSimDigSampOutData[fDetInfo.DetFullName()].fTDC_T.push_back(-1000000);
	    } 
	  }
	}
	data.clear();
      }
      

      /*
      //event.DetID.push_back(Short_t(UniqueDetID()));
      event.DetChannel[fDetInfo.DetFullName()].push_back(Short_t(m));
      event.DetNData[fDetInfo.DetFullName()].push_back(Short_t(data.size()));
      event.DetData[fDetInfo.DetFullName()].push_back(data);
      event.NDetData[fDetInfo.DetFullName()]++;

      //event.fDetectorData.push_back(data);
      //data.fData.clear();
      data.clear();
      */
    }
  }
  if(!event.fSimDigSampOutData[fDetInfo.DetFullName()].CheckSize(bool(fEncoderTDC), fDebug>=1)){
    cout << "Warning: output vectors for" << fDetInfo.DetFullName() << " don't have the same size! (any events ?" << any_events << ")" << endl;
  }
  if(!event.fSimHitMCOutData[fDetInfo.DetFullName()].CheckSize(true, true, bool(fEncoderTDC), fDebug>=1)){
  //if(!event.fSimHitMCOutData[fDetInfo.DetFullName()].CheckSize(true, true, false, fDebug>=1)){
    cout << "Warning:  output MC vectors for " << fDetInfo.DetFullName() << " don't have the same size! " << endl;
  }

  SetHasDataFlag(any_events);
}

TSBSSimHCal::Signal::Signal() : sumedep(0.0), mint(-40.0), maxt(40.0), npe(0), dnraw(10), dx_samples(4.0)
  //mint(0.0), maxt(50.), nbins(50),
{
  // hard coded, 'cause' why not? :D
  dx_raw_time = 0.120;// TDC bin (in ns)
  nbins_times = (maxt-mint)/dx_raw_time;
  times_histo.resize(nbins_times);
  nbins = (maxt-mint)/dx_samples;
  dx_raw = dx_samples/double(dnraw);
  nbins_raw= (maxt-mint)/dx_raw;
  fadc.samples.resize(nbins);
  samples_raw.resize(nbins_raw);
  //cout << maxt - mint << " " << nbins_times << " " << nbins << " " << dx_raw << " " << nbins_raw << endl;
  Clear();
}

void TSBSSimHCal::Signal::Fill(double t)
{
  int bin = (t-mint)/dx_raw_time;
  if(bin < 0 || bin > nbins_times)
    return;
  times_histo[bin]++;
  npe++;
}

void TSBSSimHCal::Signal::FillNPE(TSPEModel *model, double pulsenorm, double t, double toffset)
{
  int start_bin = 0;
  if( mint > t )
    toffset -= (mint-t);
  else
    start_bin = (t-mint)/dx_raw;

  if(start_bin > nbins_raw)
    return; // Way outside our window anyways

  // Now digitize this guy into the raw_bin (scope)
  //double tt = model->start_t-toffset;
  double tt = -12.5-toffset;//model->start_t hardcoded anyway
  //std::cout << "t=" << t << ", tt=" << tt << std::endl;
  for(int bin = start_bin; bin < nbins_raw; bin++) {
    samples_raw[bin] += pulsenorm*model->Eval(tt);
    tt += dx_raw;
  }
}

/* 
//This is the old function TSBSSimHCal::Signal::Fill, which uses struct SPEModel instead of class TSPEModel:
void TSBSSimHCal::Signal::Fill(SPEModel *model,double t, double toffset)
{
  int start_bin = 0;
  if( mint > t )
    toffset -= (mint-t);
  else
    start_bin = (t-mint)/dx_raw;

  if(start_bin > nbins_raw)
    return; // Way outside our window anyways

  // Now digitize this guy into the raw_bin (scope)
  double tt = model->start_t-toffset;
  //std::cout << "t=" << t << ", tt=" << tt << std::endl;
  for(int bin = start_bin; bin < nbins_raw; bin++) {
    samples_raw[bin] += model->Eval(tt);
    tt += dx_raw;
  }
  npe++;
}

TSBSSimHCal::SPEModel::SPEModel() :
  gain_pmt(1e6), resistance(50.0)//,qe(1.602e-19), unit(1e-9)
{
  scale = gain_pmt*resistance*qe/spe_unit;
  start_t = -12.5;
  mint = -25;
  maxt = 75.0;
  tao = 2.08*5; //ns
  sig = 2.20*5; //ns
  // test values
  tao = 2.08;
  sig = 2.20;
  t0 = 5.0;
  fFunc1 = new TF1("fFunc1HCal",TString::Format("TMath::Max(0.,"
        "(x/%g)*TMath::Exp(-x/(%g)))",tao*tao,tao),mint,maxt);
  fFunc2 = new TF1("fFunc2HCal",TString::Format(
        "%g*TMath::Exp(-((x-%g)**2)/(%g))",
        1./TMath::Sqrt(2*TMath::Pi()*sig),t0,sig*sig),mint,maxt);
  fConvolution = new TF1Convolution(fFunc1,fFunc2);

  model = new TF1("fHCalSignal",*fConvolution,mint,maxt,
        fConvolution->GetNpar());
}

double TSBSSimHCal::SPEModel::Eval(double t)
{
  return scale*model->Eval(t);
  //return model->Eval(t);
  //return 1.0;
}
*/

void TSBSSimHCal::Clear(Option_t*)
{
  for(size_t i = 0; i < fSignals.size(); i++ ) {
    fSignals[i].Clear();
  }
}

void TSBSSimHCal::Signal::Clear()
{
  for(size_t i = 0; i < fadc.samples.size(); i++) {
    fadc.samples[i] = 0;
  }
  for(size_t i = 0; i < samples_raw.size(); i++) {
    samples_raw[i] = 0;
  }
  for(size_t i = 0; i < times_histo.size(); i++) {
    times_histo[i] = 0;
  }

  npe = 0;
  met_tdc_thresh = false;
  tdc_time = mint-dx_raw;
  tdc.time.clear();
}
ClassImp(TSBSSimHCal) // Implements TSBSSimHCal
