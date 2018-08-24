#include "TSBSSimCher.h"
#include <iostream>
#include <TSBSSimData.h>
#include <TF1.h>
#include <TF1Convolution.h>
#include <TTree.h>
#include <TFile.h>
#include <TSBSSimEvent.h>
#include "TSBSDBManager.h"

TSBSSimCher::TSBSSimCher(const char* name, short id)
{
  fName = name;
  SetUniqueDetID(id);
  Init();
}

TSBSSimCher::~TSBSSimCher()
{
}

void TSBSSimCher::Init()
{
  if(fDebug>=1)
    cout << "Cherenkov detector with UniqueDetID = " << UniqueDetID() << ": TSBSSimCher::Init() " << endl;
  
  // Get the Detector info
  fDetInfo = fDBmanager->GetDetInfo(fName.Data());
  
  // Get all necessary info to parameterize the PMT pulse shape.
  double tau = fDetInfo.DigInfo().SPE_Tau();
  double sigma = fDetInfo.DigInfo().SPE_Sigma();
  double tmin = -fDetInfo.DigInfo().GateWidth()/2.0;
  double tmax = +fDetInfo.DigInfo().GateWidth()/2.0;
  double t0 = 0.0;//+fDetInfo.DigInfo().SPE_TransitTime()-fDetInfo.DigInfo().TriggerOffset();
  
  fSPE = new TSPEModel(fName.Data(), tau, sigma, t0, tmin, tmax);
  
  //Configure the PMT signals array
  fSignals.resize(fDetInfo.NChan());
  for(int i_ch = 0; i_ch<fDetInfo.NChan(); i_ch++){
    fSignals[i_ch].SetNpeChargeConv(fDetInfo.DigInfo().NpeChargeConv(i_ch));
  }
}


void TSBSSimCher::LoadEventData(const std::vector<g4sbshitdata*> &evbuffer)
{
  Clear();
  
  bool signal = false;
  int chan = 0;
  int type = 0;
  double time = 0;
  double data = 0;
  
  for( const g4sbshitdata *ev: evbuffer) {
    // Only get detector data for Cherenkov
    // new detector ID convetion proposal: UniqueDetID = 10*DetType+DetID
    if(ev->GetDetUniqueID() == UniqueDetID()) {
      signal = (ev->GetData(0)==0);
      chan = ev->GetData(1);
      type = ev->GetData(2);
      time = ev->GetData(3)+fDetInfo.DigInfo().SPE_TransitTime()-fDetInfo.DigInfo().TriggerOffset()+fDetInfo.DigInfo().TriggerJitter();//add 
      data = ev->GetData(4);
      
      if(fDebug>=3)
	cout << "Detector " << UniqueDetID() << " chan = " << chan << " Evt Time: "<< ev->GetData(3) << " " << time << endl;
      
      if(type == 0) {
        //std::cout << "Filling data for chan: " << ev->GetData(0) << ", t=" << 
        // ev->GetData(1) - 60. << std::endl;
        //if(ev->GetData(1)<mint)
        //  mint = ev->GetData(1);
	//fSPE->SetNpe(data);
        //fSignals[chan].Fill(chan, fNPE, data);//
	fSignals[chan].Fill(fSPE, data, fDetInfo.DigInfo().Threshold(chan), time, 1);//
      }
    }
  }
}

void TSBSSimCher::Digitize(TSBSSimEvent &event)
{
  bool any_events = false;

  TSBSSimEvent::DetectorData data;
  TSBSSimEvent::SimDetectorData simdata;
  
  UInt_t VETROCword;
  
  bool header[8] = {0, 0, 0, 0, 0, 0, 1, 1};
  bool channel[8];
  bool tdc[16];
  bool trail;
  short edgebitpos = 26;
  
  for(size_t m = 0; m < fSignals.size(); m++) {
    data.fData.clear();
    simdata.fData.clear();
    fSignals[m].Digitize(fDetInfo.DigInfo(), m);
    if(fSignals[m].Npe() > 0) {
      any_events = true;
      data.fDetID = UniqueDetID();
      data.fChannel = m;
      
      //define convention for type:
      // 0: ADC
      // 1: TDC
      // push back a different word for ADC and TDC ?
      // // Fill ADC 
      // data.fData.push_back(0);//ADC data flag
      // data.fData.push_back(1);//ADC data size
      // data.fData.push_back(fSignals[m].ADC());//ADC data
      // simdata.fData.clear();
      // Fill TDC 
      data.fData.push_back(1);//Digitized data
      data.fData.push_back(fSignals[m].TDCSize());
      if(fDebug>=3)cout << "TSBSSimCher::Digitize() : Unique Det ID " << UniqueDetID()  
			<< " = > fSignals[m].TDCSize() " << fSignals[m].TDCSize() << endl;
      for(int i = 0; i<fSignals[m].TDCSize(); i++){
	if(fDebug>=3)cout << " TDC " << i << " = " << fSignals[m].TDC(i) << endl;
	// Build here the vetroc word:
	// we can afford to do very ad-hoc code because 
	// all Cherenkov detectors are presumably going to use this.
	// code bits one by one... a bit tedious (and slow...)
	for(int j = 0; j<fDetInfo.DigInfo().TDCBits(); j++){
	  if(j<8){
	    //cout << j+24 << " " << header[j] << endl;
	    //cout << j+16 << " " << channel[j] << endl;
	    VETROCword ^= (-header[j] ^ VETROCword) & (1 << (j+24));
	    channel[j] = (m >> j) & 1;
	    VETROCword ^= (-channel[j] ^ VETROCword) & (1 << (j+16));
	  }
	  tdc[j] = (fSignals[m].TDC(i) >> j) & 1;
	  //cout << j << " " << tdc[j] << endl;
	  VETROCword ^= (-tdc[j] ^ VETROCword) & (1UL << j);
	}
	trail  = (fSignals[m].TDC(i) >> 31) & 1;
	VETROCword ^= (-trail ^ VETROCword) & (1UL << edgebitpos);
	//data.fData.push_back(fSignals[m].TDC(i));
	if(fDebug>=3){
	  cout << "channel " << m << " TDC " << i << " = " << fSignals[m].TDC(i) << endl;
	  if(fDebug>=5){
	    cout << "signal tdc: " << endl;
	    for(int j = 31; j>=0; j--){
	      bool bit = (fSignals[m].TDC(i) >> j) & 1;
	      cout << bit;
	    }
	    cout << endl << "vetroc word: " << endl;
	    for(int j = 31; j>=0; j--){
	      bool bit = (VETROCword >> j) & 1;
	      cout << bit;
	    }
	    cout << endl;
	  }
	}
	//Then feed here the vetroc word to the data vector
	data.fData.push_back(VETROCword);
      }
      event.fDetectorData.push_back(data);
      
      //Now take care of simulated data
      simdata.fDetID = UniqueDetID();
      simdata.fChannel = m;
      
      //define convention for type:
      // 0: SumEdep
      // 1: Npe
      // 2: Time
      // Fill SumEdep
      simdata.fData.push_back(0);
      simdata.fData.push_back(1);
      simdata.fData.push_back(fSignals[m].SumEdep());
      event.fSimDetectorData.push_back(simdata);
      simdata.fData.clear();
      // Fill Npe
      simdata.fData.push_back(1);
      simdata.fData.push_back(1);
      simdata.fData.push_back(fSignals[m].Npe());
      event.fSimDetectorData.push_back(simdata);
      simdata.fData.clear();
      // Fill Times
      simdata.fData.push_back(2);
      simdata.fData.push_back(fSignals[m].LeadTimesSize()+fSignals[m].TrailTimesSize());
      if(fDebug>=3){
	cout << "SumEdep = " << fSignals[m].SumEdep() 
	     << ", Charge " << fSignals[m].Charge() 
	     << ", Npe = " << fSignals[m].Npe() << endl;
      }
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
void TSBSSimCher::Clear()
{
  for(size_t i = 0; i < fSignals.size(); i++ ) {
    fSignals[i].Clear();
  }
}
