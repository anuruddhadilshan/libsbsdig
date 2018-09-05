#include "TSBSSimECal.h"
#include <iostream>
#include <TSBSSimData.h>
#include <TSBSSimEvent.h>
#include "TSBSDBManager.h"

TSBSSimECal::TSBSSimECal(const char* name, short id)
{
  fName = name;
  SetUniqueDetID(id);
  Init();
}

TSBSSimECal::~TSBSSimECal()
{
}

void TSBSSimECal::Init()
{
  if(fDebug>=1)
    cout << "ECal detector with UniqueDetID = " << UniqueDetID() << ": TSBSSimECal::Init() " << endl;
  
  fDetInfo = fDBmanager->GetDetInfo(fName.Data());
  
  double tau = fDetInfo.DigInfo().SPE_Tau();
  double sigma = fDetInfo.DigInfo().SPE_Sigma();
  double tmin = -fDetInfo.DigInfo().GateWidth()/2.0;
  double tmax = +fDetInfo.DigInfo().GateWidth()/2.0;
  double t0 = 0.0;//+fDetInfo.DigInfo().SPE_TransitTime()-fDetInfo.DigInfo().TriggerOffset();
  
  fSPE = new TSPEModel(fName.Data(), tau, sigma, t0, tmin, tmax);
  
  fSignals.resize(fDetInfo.NChan());
  for(size_t i_ch = 0; i_ch<fDetInfo.NChan(); i_ch++){
    fSignals[i_ch].SetNpeChargeConv(fDetInfo.DigInfo().NpeChargeConv(i_ch));
  }
}


void TSBSSimECal::LoadEventData(const std::vector<g4sbshitdata*> &evbuffer)
{
  Clear();
  LoadAccumulateData(evbuffer);
  /*
  bool signal = false;
  int chan = 0;
  int type = 0;
  double time = 0;
  double data = 0;
  
  // for( const g4sbshitdata *ev: evbuffer) {
  for(std::vector<g4sbshitdata*>::const_iterator it = evbuffer.begin(); it!= evbuffer.end(); ++it ) {
    g4sbshitdata* ev = (*it);
    // Only get detector data for ECal
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
	fSignals[chan].Fill(fSPE, data, fDetInfo.DigInfo().Threshold(chan), time, signal);//
	if(fDebug>=3)
	  cout << "chan " << chan << " data " << data 
	       << " fSignals[chan].Npe() " << fSignals[chan].Npe() << endl;
      } else if (type == 1) { // sumedep data
        fSignals[chan].AddSumEdep(data);
	if(fDebug>=3)
	  cout << "chan " << chan << " data " << data 
	       << " fSignals[chan].SumEdep() " << fSignals[chan].SumEdep() << endl;
      }
    }
  }
  */
}

void TSBSSimECal::LoadAccumulateData(const std::vector<g4sbshitdata*> &evbuffer)
{
  bool signal = false;
  int chan = 0;
  int type = 0;
  double time = 0;
  double data = 0;
  
  // for( const g4sbshitdata *ev: evbuffer) {
  for(std::vector<g4sbshitdata*>::const_iterator it = evbuffer.begin(); it!= evbuffer.end(); ++it ) {
    g4sbshitdata* ev = (*it);
    // Only get detector data for ECal
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
	fSignals[chan].Fill(fSPE, data, fDetInfo.DigInfo().Threshold(chan), time, signal);//
	if(fDebug>=3)
	  cout << "chan " << chan << " data " << data 
	       << " fSignals[chan].Npe() " << fSignals[chan].Npe() << endl;
      } else if (type == 1) { // sumedep data
        fSignals[chan].AddSumEdep(data);
	if(fDebug>=3)
	  cout << "chan " << chan << " data " << data 
	       << " fSignals[chan].SumEdep() " << fSignals[chan].SumEdep() << endl;
      }
    }
  }
}

void TSBSSimECal::Digitize(TSBSSimEvent &event)
{
  bool any_events = false;
  
  TSBSSimEvent::DetectorData data;
  TSBSSimEvent::SimDetectorData simdata;
  
  UInt_t TDCword;
  
  bool header[8] = {0, 0, 0, 0, 0, 0, 0, 0};//bits ...-31 
  // we ignore (for the moment) all other informations for the ADC.
  bool channel[7];//Common features between V1190A and 1877: 7 digit channel
  bool tdc[19];
  bool trail;
  //hardcoded, but add all ADC/TDC bit coding in DB shall be cumbersome
  short edgebitpos = 16;
  short nheaderbits = 32-8-fDetInfo.DigInfo().TDCBits();
  if(UniqueDetID()==30){
    edgebitpos = 26;
  }
  short chanfirstbit = fDetInfo.DigInfo().TDCBits()+Short_t(edgebitpos==fDetInfo.DigInfo().TDCBits());
  
  for(size_t m = 0; m < fSignals.size(); m++) {
    data.fData.clear();
    simdata.fData.clear();
    fSignals[m].Digitize(fDetInfo.DigInfo(), m);
    if(fSignals[m].Npe() > 0) {
      any_events = true;
      data.fDetID = UniqueDetID();
      data.fChannel = m;
      
      if(fDebug>=3)cout << "TSBSSimECal::Digitize() : Unique Det ID " << UniqueDetID() 
			<< " = > fSignals[m].ADC() " << fSignals[m].ADC() << endl;
      //define convention for type:
      // 0: ADC
      // 1: TDC
      // push back a different word for ADC and TDC ?
      // Fill ADC 
      data.fData.push_back(0);//ADC data flag
      data.fData.push_back(1);//ADC data size
      data.fData.push_back(fSignals[m].ADC());//ADC data
      event.fDetectorData.push_back(data);
      data.fData.clear();
      
      // Fill TDC 
      if(fDetInfo.DigInfo().TDCBits()>0 && fDetInfo.DigInfo().TDCConversion()>0){
	data.fData.push_back(1);//TDCs data
	data.fData.push_back(fSignals[m].TDCSize());//TDC data size
	if(fDebug>=3)cout << "TSBSSimECal::Digitize() : Unique Det ID " << UniqueDetID()  
			  << " = > fSignals[m].TDCSize() " << fSignals[m].TDCSize() << endl;
	for(size_t i = 0; i<fSignals[m].TDCSize(); i++){
	  if(fDebug>=3)cout << " TDC " << i << " = " << fSignals[m].TDC(i) << endl;
	  
	  // Build here the TDC word:
	  //code bits one by one... a bit tedious (and slow...)
	  for(int j = 0; j<fDetInfo.DigInfo().TDCBits(); j++){
	    if(j<nheaderbits){
	      //cout << j+24 << " " << header[j] << endl;
	    //cout << j+16 << " " << channel[j] << endl;
	      TDCword ^= (-header[j] ^ TDCword) & (1 << (j+32-nheaderbits));
	    }
	    if(j<7){
	      channel[j] = (m >> j) & 1;
	      TDCword ^= (-channel[j] ^ TDCword) & (1 << (j+chanfirstbit));
	  }
	    tdc[j] = (fSignals[m].TDC(i) >> j) & 1;
	    //cout << j << " " << tdc[j] << endl;
	    TDCword ^= (-tdc[j] ^ TDCword) & (1UL << j);
	  }
	  trail  = (fSignals[m].TDC(i) >> 31) & 1;
	  TDCword ^= (-trail ^ TDCword) & (1UL << edgebitpos);
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
		bool bit = (TDCword >> j) & 1;
		cout << bit;
	      }
	      cout << endl;
	    }
	  }
	  //Then feed here the TDC word to the data vector
	  data.fData.push_back(TDCword);
	}
	event.fDetectorData.push_back(data);
	data.fData.clear();
      }
      
      //Now take care of simulated data
      simdata.fDetID = UniqueDetID();
      simdata.fChannel = m;
      
      //define convention for type:
      // 0: SumEdep
      // 1: Npe
      // 2: Time
      // Fill SumEdep
      simdata.fDataType = 0;
      simdata.fNdata = 1;
      // simdata.fData.push_back(0);
      // simdata.fData.push_back(1);
      simdata.fData.push_back(fSignals[m].SumEdep());
      event.fSimDetectorData.push_back(simdata);
      simdata.fData.clear();
      // Fill Npe
      simdata.fDataType = 1;
      simdata.fNdata = 1;
      // simdata.fData.push_back(1);
      // simdata.fData.push_back(1);
      simdata.fData.push_back(fSignals[m].Npe());
      event.fSimDetectorData.push_back(simdata);
      simdata.fData.clear();
      /*
      // Fill Times
      simdata.fData.push_back(2);
      simdata.fData.push_back(fSignals[m].LeadTimesSize()+fSignals[m].TrailTimesSize());
      if(fDebug>=3){
	cout << "SumEdep = " << fSignals[m].SumEdep() 
	     << ", Charge " << fSignals[m].Charge() 
	     << ", Npe = " << fSignals[m].Npe() << endl;
      }
      //data.fData.push_back(fSignals[m].SumEdep());
      //data.fData.push_back(fSignals[m].Npe());
      
      for(size_t i = 0; i<fSignals[m].LeadTimesSize(); i++){
	simdata.fData.push_back(fSignals[m].LeadTime(i));
	if(fDebug>=3)cout << " leadtime " << i << " = " << fSignals[m].LeadTime(i) << endl;
      }
      for(size_t i = 0; i<fSignals[m].TrailTimesSize(); i++){
	simdata.fData.push_back(fSignals[m].TrailTime(i));
	if(fDebug>=3)cout << " trail time " << i << " = " << fSignals[m].TrailTime(i) << endl;;
      }
      event.fSimDetectorData.push_back(simdata);
      simdata.fData.clear();
      */
    }
  }
  SetHasDataFlag(any_events);
  
}

// Clear signals in array
void TSBSSimECal::Clear()
{
  for(size_t i = 0; i < fSignals.size(); i++ ) {
    fSignals[i].Clear();
  }
}

ClassImp(TSBSSimECal) // Implements TSBSSimECal
