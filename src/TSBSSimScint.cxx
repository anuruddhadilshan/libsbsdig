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
  SetUniqueDetID(id);
  Init();
}

TSBSSimScint::~TSBSSimScint()
{
}

void TSBSSimScint::Init()
{
  if(fDebug>=1)
    cout << "Scintillator detector with UniqueDetID = " << UniqueDetID() << ": TSBSSimScint::Init() " << endl;
  
  fDetInfo = fDBmanager->GetDetInfo(fName.Data());
  
  double tau = fDetInfo.DigInfo().SPE_Tau();
  double sigma = fDetInfo.DigInfo().SPE_Sigma();
  double tmin = -fDetInfo.DigInfo().GateWidth()/2.0;
  double tmax = +fDetInfo.DigInfo().GateWidth()/2.0;
  double t0 = 0.0;//+fDetInfo.DigInfo().SPE_TransitTime()-fDetInfo.DigInfo().TriggerOffset();
  
  fSPE = new TSPEModel(fName.Data(), tau, sigma, t0, tmin, tmax);
  
  fSignals.resize(fDetInfo.NChan());
  for(int i_ch = 0; i_ch<fDetInfo.NChan(); i_ch++){
    fSignals[i_ch].SetNpeChargeConv(fDetInfo.DigInfo().NpeChargeConv(i_ch));
  }
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
    if(ev->GetDetUniqueID() == UniqueDetID()) {
      signal = (ev->GetData(0)==0);
      chan = ev->GetData(1);
      type = ev->GetData(2);
      time = ev->GetData(3)+fDetInfo.DigInfo().SPE_TransitTime()-fDetInfo.DigInfo().TriggerOffset()+fDetInfo.DigInfo().TriggerJitter();//add 
      data = ev->GetData(4);
      
      if(fDebug>=3)
	cout << " chan = " << chan << " Evt Time: "<< ev->GetData(3) << " " << time << endl;
      
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
	if(fDebug>=3)
	  cout << "chan " << chan << " data " << data 
	       << " fSignals[chan].SumEdep() " << fSignals[chan].SumEdep() << endl;
      }
    }
  }
}

void TSBSSimScint::Digitize(TSBSSimEvent &event)
{
  bool any_events = false;
  for(size_t m = 0; m < fSignals.size(); m++) {
    fSignals[m].Digitize(fDetInfo.DigInfo(), m);
    if(fSignals[m].Npe() > 0) {
      any_events = true;
      TSBSSimEvent::DetectorData data;
      data.fDetID = UniqueDetID();
      data.fChannel = m;
      
      data.fData.push_back(0);//Digitized data
      data.fData.push_back(fSignals[m].TDCSize());
      if(fDebug>=3)cout << "TSBSSimScint::Digitize() : Unique Det ID " << UniqueDetID()  
			<< " = > fSignals[m].TDCSize() " << fSignals[m].TDCSize() << endl;
      for(int i = 0; i<fSignals[m].TDCSize(); i++){
	data.fData.push_back(fSignals[m].TDC(i));
	if(fDebug>=3)cout << " TDC " << i << " = " << fSignals[m].TDC(i) << endl;
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
      simdata.fDetID = UniqueDetID();
      simdata.fChannel = m;
      simdata.fData.push_back(1);
      simdata.fData.push_back(fSignals[m].LeadTimesSize()+fSignals[m].TrailTimesSize());
      if(fDebug>=3){
	cout << "SumEdep = " << fSignals[m].SumEdep() 
	     << ", Charge " << fSignals[m].Charge() 
	     << ", Npe = " << fSignals[m].Npe() << endl;
      }
      //data.fData.push_back(fSignals[m].SumEdep());
      //data.fData.push_back(fSignals[m].Npe());
      for(int i = 0; i<fSignals[m].LeadTimesSize(); i++){
	simdata.fData.push_back(fSignals[m].LeadTime(i));
	if(fDebug>=3)cout << " leadtime " << i << " = " << fSignals[m].LeadTime(i) << endl;
      }
      for(int i = 0; i<fSignals[m].TrailTimesSize(); i++){
	simdata.fData.push_back(fSignals[m].TrailTime(i));
	if(fDebug>=3)cout << " trail time " << i << " = " << fSignals[m].TrailTime(i) << endl;;
      }
      event.fSimDetectorData.push_back(simdata);
    }
  }
  SetHasDataFlag(any_events);
  
}

// Clear signals in array
void TSBSSimScint::Clear()
{
  for(size_t i = 0; i < fSignals.size(); i++ ) {
    fSignals[i].Clear();
  }
}

