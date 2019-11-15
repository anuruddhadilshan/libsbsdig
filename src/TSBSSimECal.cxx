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
  TSBSSimDetector::Init();
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
      time = ev->GetData(3)+fDetInfo.DigInfo().SPE_TransitTime()-fDetInfo.DigInfo().TriggerOffset() + fTimeZero;//+fDetInfo.DigInfo().TriggerJitter()
      data = ev->GetData(4);
      
      //if(fabs(time)>fDetInfo.DigInfo().GateWidth()/2.0)continue;
      
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

  if(fDebug>=3)cout << "TSBSSimECal::Digitize() : Unique Det ID " << UniqueDetID() << " signal size = " << fSignals.size() << endl;
  
  //TSBSSimEvent::DetectorData data;
  //TSBSSimEvent::SimDetectorData simdata;
  
  /*
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
 */ 
  int mult = 0;
  SimEncoder::adc_data adc_data;
  std::vector<uint32_t> data;
  //std::vector<double> simdata;
  for(size_t m = 0; m < fSignals.size(); m++) {
    //data.fData.clear();
    //simdata.fData.clear();
    data.clear();
    //simdata.clear();
    fSignals[m].Digitize(fDetInfo.DigInfo(), m);
    if(fSignals[m].Npe() > 0) {
      any_events = true;
      //data.fDetID = UniqueDetID();
      //data.fChannel = m;
      
      if(fDebug>=4)cout << " = > fSignals[" << m << "].ADC() " << fSignals[m].ADC() << endl;
      
      /*
      //Fill the simulation data
      //event.SimDetID.push_back(Short_t(UniqueDetID()));
      event.SimDetChannel[fDetInfo.DetFullName()].push_back(Short_t(m));
      event.SimDetDataType[fDetInfo.DetFullName()].push_back(0);
      event.SimDetNData[fDetInfo.DetFullName()].push_back(1);
      simdata.push_back(fSignals[m].SumEdep());
      event.SimDetData[fDetInfo.DetFullName()].push_back(simdata);
      event.NSimDetData[fDetInfo.DetFullName()]++;
      simdata.clear();
      
      //event.SimDetID.push_back(Short_t(UniqueDetID()));
      event.SimDetChannel[fDetInfo.DetFullName()].push_back(Short_t(m));
      event.SimDetDataType[fDetInfo.DetFullName()].push_back(1);
      event.SimDetNData[fDetInfo.DetFullName()].push_back(1);
      simdata.push_back(fSignals[m].Npe());
      event.SimDetData[fDetInfo.DetFullName()].push_back(simdata);
      event.NSimDetData[fDetInfo.DetFullName()]++;
      simdata.clear();
      
      event.SimDetChannel[fDetInfo.DetFullName()].push_back(Short_t(m));
      event.SimDetDataType[fDetInfo.DetFullName()].push_back(2);
      event.SimDetNData[fDetInfo.DetFullName()].push_back(1);
      simdata.push_back(fSignals[m].EventTime());
      event.SimDetData[fDetInfo.DetFullName()].push_back(simdata);
      event.NSimDetData[fDetInfo.DetFullName()]++;
      simdata.clear();
      */
      if(fDebug>=4){
	cout << fDetInfo.DetFullName() << ": check MC vec size" << endl;
	fSignals[m].check_vec_size();
      }
      for(uint i_mc = 0; i_mc<fSignals[m].MCHitSize(); i_mc++){
	event.NSimDetHits[fDetInfo.DetFullName()]++;
	event.SimDetChannel[fDetInfo.DetFullName()].push_back(Short_t(m));
	event.SimDetEdep[fDetInfo.DetFullName()].push_back(fSignals[m].MCHitEdep(i_mc));
	event.SimDetNpe[fDetInfo.DetFullName()].push_back(fSignals[m].MCHitNpe(i_mc));
	event.SimDetTime[fDetInfo.DetFullName()].push_back(fSignals[m].MCHitTime(i_mc));
	if(fEncoderTDC){
	  event.SimDetLeadTime[fDetInfo.DetFullName()].push_back(fSignals[m].MCHitLeadTime(i_mc));
	  event.SimDetTrailTime[fDetInfo.DetFullName()].push_back(fSignals[m].MCHitTrailTime(i_mc));
	}
      }
      
      //define convention for type:
      // 0: ADC
      // 1: TDC
      mult = 0;
      
      // Fill ADC 
      if(fEncoderADC) {
        adc_data.integral=fSignals[m].ADC();
        fEncoderADC->EncodeADC(adc_data,fEncBuffer,fNEncBufferWords);
        CopyEncodedData(fEncoderADC,mult++,data);//.fData);

	for(uint i = 0; i<data.size(); i++){
	  event.NDetHits[fDetInfo.DetFullName()]++;
	  event.DetChannel[fDetInfo.DetFullName()].push_back(Short_t(m));
	  event.DetDataWord[fDetInfo.DetFullName()].push_back(data.at(i));
	  // if(i==0){//header
	  //   event.DetADC[fDetInfo.DetFullName()].push_back(-1);
	  //   event.DetTDC[fDetInfo.DetFullName()].push_back(-1);
	  // }else{
	  //   event.DetADC[fDetInfo.DetFullName()].push_back(fSignals[m].ADC());
	  //   if(fEncoderTDC)event.DetTDC[fDetInfo.DetFullName()].push_back(-1);
	  // }
	}
	data.clear();
      }
      
      // Fill TDC
      if(fEncoderTDC) {
        fEncoderTDC->EncodeTDC(fSignals[m].TDCData(),fEncBuffer,
            fNEncBufferWords);
        CopyEncodedData(fEncoderTDC,mult++,data);//.fData);

	for(uint i = 0; i<data.size(); i++){
	  event.NDetHits[fDetInfo.DetFullName()]++;
	  event.DetChannel[fDetInfo.DetFullName()].push_back(Short_t(m));
	  event.DetDataWord[fDetInfo.DetFullName()].push_back(data.at(i));
	  // if(i==0){//header
	  //   event.DetADC[fDetInfo.DetFullName()].push_back(-1);
	  //   event.DetTDC[fDetInfo.DetFullName()].push_back(-1);
	  // }else{
	  //   if(fEncoderADC)event.DetADC[fDetInfo.DetFullName()].push_back(-1);
	  //   event.DetTDC[fDetInfo.DetFullName()].push_back(fSignals[m].TDC(i-1));
	  // }
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
  SetHasDataFlag(any_events);
  
}

// Clear signals in array
void TSBSSimECal::Clear(Option_t*)
{
  for(size_t i = 0; i < fSignals.size(); i++ ) {
    fSignals[i].Clear();
  }
}

ClassImp(TSBSSimECal) // Implements TSBSSimECal
