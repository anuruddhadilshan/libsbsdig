#include "TSBSSimHCal.h"
#include <iostream>
#include <TSBSSimData.h>
//#include <TF1.h>
//#include <TF1Convolution.h>
//#include <TTree.h>
//#include <TFile.h>
#include <TSBSSimEvent.h>
#include "TSBSDBManager.h"

#define HCAL_TDC_THRESH 0.1

TSBSSimHCal::TSBSSimHCal(const char* name, short id)
{
  fName = name;
  SetUniqueDetID(id);
  Init();
}

TSBSSimHCal::~TSBSSimHCal()
{
}

void TSBSSimHCal::Init()
{
  if(fDebug>=1)
    cout << "HCal detector with UniqueDetID = " << UniqueDetID() << ": TSBSSimHCal::Init() " << endl;
  
  fDetInfo = fDBmanager->GetDetInfo(fName.Data());
  
  double tau = fDetInfo.DigInfo().SPE_Tau();
  double sigma = fDetInfo.DigInfo().SPE_Sigma();
  double tmin = -fDetInfo.DigInfo().GateWidth()/2.0;
  double tmax = +fDetInfo.DigInfo().GateWidth()/2.0;
  double t0 = 0.0;//+fDetInfo.DigInfo().SPE_TransitTime()-fDetInfo.DigInfo().TriggerOffset();
  
  fSPE = new TSPEModel(fName.Data(), tau, sigma, t0, tmin, tmax);
  
  fSignals.resize(fDetInfo.NChan());
  /*
    fSPE = new SPEModel();
    //fSPE = new SPEModel( new TF1("fHCalSignal",*fConvolution,mint,maxt,
    //    fConvolution->GetNpar()));
    fSignals.resize(288); // TODO: Don't hard code this!!!
  */
  //fFileOut = new TFile("rootfiles/testout.root","RECREATE");
  //fTreeOut = new TTree("TTest","");
  /*for(int m = 0; m < int(fSignals.size()); m++) {
    fTreeOut->Branch(TString::Format("m%d.npe",m),&(fSignals[m].npe));
    fTreeOut->Branch(TString::Format("m%d.sum",m),&(fSignals[m].sum));
    fTreeOut->Branch(TString::Format("m%d.samples",m),&(fSignals[m].samples));
    }
  */
}


void TSBSSimHCal::LoadEventData(const std::vector<g4sbshitdata*> &evbuffer)
{
  Clear();
  LoadAccumulateData(evbuffer);
  /*
  // Just make HCAL be 288 modules for now to make it easier....
  //Double_t mint = 1e9;
  int mod = 0;
  int type = 0;
  double data = 0;
  double pulsenorm = 0;
  // for( const g4sbshitdata *ev: evbuffer) {
  for(std::vector<g4sbshitdata*>::const_iterator it = evbuffer.begin(); it!= evbuffer.end(); ++it ) {
    g4sbshitdata* ev = (*it);
    // Only get detector data for HCAL
    // TODO: Don't hard code DetID here!!!
    if(ev->GetDetType() == kHCal) {
      mod  = ev->GetData(0);
      type = ev->GetData(1);
      data = ev->GetData(2);
      if(type == 0) {
        //std::cout << "Filling data for mod: " << ev->GetData(0) << ", t=" << 
        // ev->GetData(1) - 60. << std::endl;
        //if(ev->GetData(1)<mint)
        //  mint = ev->GetData(1);
	pulsenorm = fDetInfo.DigInfo().Gain(mod)*fDetInfo.DigInfo().ROImpedance()*qe/spe_unit;
        //fSignals[mod].Fill(fSPE, data-75.);
	fSignals[mod].Fill(fSPE, pulsenorm,data-75.);
      } else if (type == 1) { // sumedep data
        fSignals[mod].sumedep = data;
      }
    }
  }
  //std::cout << "Mint = " << mint << std::endl;
  */
}


void TSBSSimHCal::LoadAccumulateData(const std::vector<g4sbshitdata*> &evbuffer)
{
  // Just make HCAL be 288 modules for now to make it easier....
  //Double_t mint = 1e9;
  int mod = 0;
  int type = 0;
  double data = 0;
  double pulsenorm = 0;
  // for( const g4sbshitdata *ev: evbuffer) {
  for(std::vector<g4sbshitdata*>::const_iterator it = evbuffer.begin(); it!= evbuffer.end(); ++it ) {
    g4sbshitdata* ev = (*it);
    // Only get detector data for HCAL
    // TODO: Don't hard code DetID here!!!
    if(ev->GetDetType() == kHCal) {
      mod  = ev->GetData(0);
      type = ev->GetData(1);
      data = ev->GetData(2);
      if(type == 0) {
        //std::cout << "Filling data for mod: " << ev->GetData(0) << ", t=" << 
        // ev->GetData(1) - 60. << std::endl;
        //if(ev->GetData(1)<mint)
        //  mint = ev->GetData(1);
	pulsenorm = fDetInfo.DigInfo().Gain(mod)*fDetInfo.DigInfo().ROImpedance()*qe/spe_unit;
        //fSignals[mod].Fill(fSPE, data-75.);
	fSignals[mod].Fill(fSPE, pulsenorm,data-75.);
      } else if (type == 1) { // sumedep data
        fSignals[mod].sumedep = data;
      }
    }
  }
  //std::cout << "Mint = " << mint << std::endl;
}

void TSBSSimHCal::Signal::Digitize()
{
  if(npe <= 0)
    return;

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
    if(max>2)
      max = 2;
    //samples[bs] =int((max/2.);//*4095);
    samples[bs] =int((max/2.)*4095);
    //samples[bs] = samples[bs] > 4095 ? 4095 : samples[bs];
    braw += dnraw;
    sum += samples[bs];
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
    tdc_time = int((tdc_time/3.9e3)*65535);
  }
}

void TSBSSimHCal::Digitize(TSBSSimEvent &event)
{
  bool any_events = false;
  for(size_t m = 0; m < fSignals.size(); m++) {
    fSignals[m].Digitize();
    if(fSignals[m].npe > 0) {
      any_events = true;
      TSBSSimEvent::DetectorData data;
      data.fDetID = HCAL_UNIQUE_DETID; // 2 for fADC data
      data.fChannel = m;
      //data.fData.push_back(data.fChannel);
      //data.fData.push_back(m);
      data.fData.push_back(0); // For samples data
      data.fData.push_back(fSignals[m].samples.size());
      //std::cout << "Module : " << m << " npe=" << fSignals[m].npe;
      for(size_t j = 0; j < fSignals[m].samples.size(); j++) {
        //std::cout << " " << fSignals[m].samples[j];
        data.fData.push_back(fSignals[m].samples[j]);
      }
      event.fDetectorData.push_back(data);
      // Now add the sum (or edep)
      data.fData.clear();
      //data.fData.push_back(m);
      // Since it's still uncertain if we can populate the sumedet
      // parts, I'll leave this out for now...
      //data.fData.push_back(1);
      //data.fData.push_back(1);
      //data.fData.push_back(fSignals[m].sumedep);
      //event.fDetectorData.push_back(data);

      // Now add the TDC if the threshold was met
      if( fSignals[m].met_tdc_thresh) {
        data.fData.clear();
        data.fData.push_back(1);
        data.fData.push_back(1);
        data.fData.push_back(fSignals[m].tdc_time);
        event.fDetectorData.push_back(data);
      }
    }
  }
  SetHasDataFlag(any_events);
}

TSBSSimHCal::Signal::Signal() : sumedep(0.0), mint(0.0), maxt(50.0), npe(0), dnraw(10), dx_samples(1.0)
  //mint(0.0), maxt(50.), nbins(50),
{
  nbins = (maxt-mint)/dx_samples;
  dx_raw = dx_samples/double(dnraw);
  nbins_raw= (maxt-mint)/dx_raw;
  samples.resize(nbins);
  samples_raw.resize(nbins_raw);
}

void TSBSSimHCal::Signal::Fill(TSPEModel *model, double pulsenorm, double t, double toffset)
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
  npe++;
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

void TSBSSimHCal::Clear()
{
  for(size_t i = 0; i < fSignals.size(); i++ ) {
    fSignals[i].Clear();
  }
}

void TSBSSimHCal::Signal::Clear()
{
  for(size_t i = 0; i < samples.size(); i++) {
    samples[i] = 0;
  }
  for(size_t i = 0; i < samples_raw.size(); i++) {
    samples_raw[i] = 0;
  }
  npe = 0;
  met_tdc_thresh = false;
  tdc_time = mint-dx_raw;
}
ClassImp(TSBSSimHCal) // Implements TSBSSimHCal
