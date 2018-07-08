#include "TSBSSimScint.h"
#include <iostream>
#include <TSBSSimData.h>
#include <TF1.h>
#include <TF1Convolution.h>
#include <TTree.h>
#include <TFile.h>
#include <TSBSSimEvent.h>
#include "TSBSDBManager.h"

TSBSSimScint::TSBSSimScint(const char* name)
{
  fName = name;
  Init();
}

TSBSSimScint::~TSBSSimScint()
{
}

void TSBSSimScint::Init()
{
  //fSPE = new SPEModel( new TF1("fHCalSignal",*fConvolution,mint,maxt,
  //    fConvolution->GetNpar()));
  //fSignals.resize(180); // TODO: Don't hard code this!!!
  
  fDetInfo = fDBmanager->GetDetInfo(fName.Data());
  fSPE = new SPEModel(fDetInfo.fDigInfo, fName.Data());
  //fSPE->DigInfo = fDetInfo.fDigInfo;
    
  fSignals.resize(fDetInfo.fNChan);
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
  /*
  Clear();
  // Just make HCAL be 288 modules for now to make it easier....
  //Double_t mint = 1e9;
  int mod = 0;
  int type = 0;
  double data = 0;
  for( const g4sbshitdata *ev: evbuffer) {
    // Only get detector data for Scintillator
    // new detector ID convetion proposal: UniqueDetID = 10*DetType+DetID
    if(ev->GetDetType() == kScint) {
      mod  = ev->GetData(0);
      type = ev->GetData(1);
      data = ev->GetData(2);
      if(type == 0) {
        //std::cout << "Filling data for mod: " << ev->GetData(0) << ", t=" << 
        // ev->GetData(1) - 60. << std::endl;
        //if(ev->GetData(1)<mint)
        //  mint = ev->GetData(1);
        fSignals[mod].Fill(fSPE,data-75.);
      } else if (type == 1) { // sumedep data
        fSignals[mod].sumedep = data;
      }
    }
  }
  //std::cout << "Mint = " << mint << std::endl;
  */
}

void TSBSSimScint::Signal::Digitize()
{
  /*
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

void TSBSSimScint::Digitize(TSBSSimEvent &event)
{
  /*
  bool any_events = false;
  for(size_t m = 0; m < fSignals.size(); m++) {
    fSignals[m].Digitize();
    if(fSignals[m].npe > 0) {
      any_events = true;
      TSBSSimEvent::DetectorData data;
      data.fDetID = 2; // 2 for fADC data
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
      data.fData.clear();
      //data.fData.push_back(m);
      data.fData.push_back(1);
      data.fData.push_back(1);
      data.fData.push_back(fSignals[m].sumedep);
      event.fDetectorData.push_back(data);
    }
  }
  SetHasDataFlag(any_events);
  */
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
void TSBSSimScint::Signal::Fill(SPEModel *model,double t, double toffset)
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
