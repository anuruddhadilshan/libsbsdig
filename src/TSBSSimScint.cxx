#include "TSBSSimScint.h"
#include <iostream>
#include <TSBSSimData.h>
#include <TF1.h>
#include <TF1Convolution.h>
#include <TTree.h>
#include <TFile.h>
#include <TSBSSimEvent.h>
#include "TSBSDBManager.h"

TSBSSimScint::TSBSSimScint(const char* name, short id)
{
  fName = name;
  SetUniqueID(id);
  Init();
}

TSBSSimScint::~TSBSSimScint()
{
}

void TSBSSimScint::Init()
{
  //fNPE = new NPEModel( new TF1("fHCalSignal",*fConvolution,mint,maxt,
  //    fConvolution->GetNpar()));
  //fSignals.resize(180); // TODO: Don't hard code this!!!
  fDetInfo = fDBmanager->GetDetInfo(fName.Data());
  //fNPE = new TNPEModel(fDetInfo.fDigInfo, fName.Data());
  //fNPE->DigInfo = fDetInfo.fDigInfo;
  
  double tau = fDetInfo.DigInfo().SPE_Tau();
  double sigma = fDetInfo.DigInfo().SPE_Sigma();
  double tmin = -fDetInfo.DigInfo().GateWidth()/2.0;
  double tmax = +fDetInfo.DigInfo().GateWidth()/2.0;
  double t0 = +fDetInfo.DigInfo().SPE_TransitTime()-fDetInfo.DigInfo().TriggerOffset();
  
  fSPE = new TSPEModel(fName.Data(), tau, sigma, t0, tmin, tmax);
  
  fSignals.resize(fDetInfo.NChan());
  for(int i_ch = 0; i_ch<fDetInfo.NChan(); i_ch++){
    fSignals[i_ch].SetNpeChargeConv(fDetInfo.DigInfo().NpeChargeConv(i_ch));
  }
  //fFileOut = new TFile("rootfiles/testout.root","RECREATE");
  //fTreeOut = new TTree("TTest","");
  /*for(int m = 0; m < int(fSignals.size()); m++) {
    fTreeOut->Branch(TString::Format("m%d.npe",m),&(fSignals[m].npe));
    fTreeOut->Branch(TString::Format("m%d.sum",m),&(fSignals[m].sum));
    fTreeOut->Branch(TString::Format("m%d.samples",m),&(fSignals[m].samples));
  }
  */
}


void TSBSSimScint::LoadEventData(const std::vector<g4sbshitdata*> &evbuffer)
{
  Clear();
  
  bool signal = false;
  int chan = 0;
  int type = 0;
  double time = 0;
  double data = 0;
  
  for( const g4sbshitdata *ev: evbuffer) {
    // Only get detector data for Scintillator
    // new detector ID convetion proposal: UniqueDetID = 10*DetType+DetID
    if(ev->GetDetType() == kScint) {
      signal = (ev->GetData(0)==0);
      chan = ev->GetData(1);
      type = ev->GetData(2);
      time = ev->GetData(3)+
	fDetInfo.DigInfo().SPE_TransitTime()-fDetInfo.DigInfo().TriggerOffset()+fDetInfo.DigInfo().TriggerJitter();//add 
      data = ev->GetData(4);
      
      if(fDebug>=3)cout << "Evt Time: "<< ev->GetData(3) << " " << time << endl;
      
      if(type == 0) {
        //std::cout << "Filling data for chan: " << ev->GetData(0) << ", t=" << 
        // ev->GetData(1) - 60. << std::endl;
        //if(ev->GetData(1)<mint)
        //  mint = ev->GetData(1);
	//fSPE->SetNpe(data);
        //fSignals[chan].Fill(chan, fNPE, data);//
	fSignals[chan].Fill(fSPE, data, fDetInfo.DigInfo().Threshold(chan), time, 1);//
      } else if (type == 1) { // sumedep data
        fSignals[chan].AddSumEdep(data);
      }
    }
  }
  //std::cout << "Mint = " << mint << std::endl;
}

void TSBSSimScint::Digitize(TSBSSimEvent &event)
{
  
  bool any_events = false;
  for(size_t m = 0; m < fSignals.size(); m++) {
    fSignals[m].Digitize(fDetInfo.DigInfo(), m);
    if(fSignals[m].Npe() > 0) {
      any_events = true;
      TSBSSimEvent::DetectorData data;
      cout << UniqueID() << endl;
      data.fDetID = UniqueID();//need 
      data.fChannel = m;
      //for scintillators, we only need to 
      data.fData.push_back(0);//Digitized data
      data.fData.push_back(fSignals[m].TDCSize());
      for(int i = 0; i<fSignals[m].TDCSize(); i++){
	data.fData.push_back(fSignals[m].TDC(i));
	if(fDebug>=3)cout << i << " " << fSignals[m].TDC(i) << endl;
      }
      //data.fData.push_back(m);
      //data.fData.push_back(0); // For samples data
      //data.fData.push_back(fSignals[m].samples.size());
      //std::cout << "Module : " << m << " npe=" << fSignals[m].npe;
      // for(size_t j = 0; j < fSignals[m].samples.size(); j++) {
      //   //std::cout << " " << fSignals[m].samples[j];
      //   data.fData.push_back(fSignals[m].samples[j]);
      // }
      event.fDetectorData.push_back(data);
      TSBSSimEvent::SimDetectorData simdata;
      simdata.fDetID = UniqueID();//need 
      simdata.fChannel = m;
      simdata.fData.push_back(1+fSignals[m].LeadTimesSize()+fSignals[m].TrailTimesSize());
      data.fData.push_back(fSignals[m].Sumedep());
      for(int i = 0; i<fSignals[m].LeadTimesSize(); i++){
	simdata.fData.push_back(fSignals[m].LeadTime(i));
      }
      for(int i = 0; i<fSignals[m].TrailTimesSize(); i++){
	simdata.fData.push_back(fSignals[m].TrailTime(i));
      }
      event.fSimDetectorData.push_back(simdata);
    }
  }
  SetHasDataFlag(any_events);
  
}

/*
TSBSSimScint::Signal::Signal() : mint(0.0), maxt(50.0), dx_samples(1.0), npe(0),
  dnraw(10), sumedep(0.0)
  //mint(0.0), maxt(50.), nbins(50),
{
  nbins = (maxt-mint)/dx_samples;
  dx_raw = dx_samples/double(dnraw);
  nbins_raw= (maxt-mint)/dx_raw;
  samples.resize(nbins);
  samples_raw.resize(nbins_raw);
}
*/

/*
void TSBSSimScint::Signal::Fill(NPEModel *model,double t, double toffset)
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
*/


void TSBSSimScint::Clear()
{
  for(size_t i = 0; i < fSignals.size(); i++ ) {
    fSignals[i].Clear();
  }
}
/*
void TSBSSimScint::Signal::Clear()
{
  // for(size_t i = 0; i < samples.size(); i++) {
  //   samples[i] = 0;
  // }
  // for(size_t i = 0; i < samples_raw.size(); i++) {
  //   samples_raw[i] = 0;
  // }
  npe = 0;
}
*/
