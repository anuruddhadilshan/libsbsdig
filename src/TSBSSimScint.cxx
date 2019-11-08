#include "TSBSSimScint.h"
#include <iostream>
#include <TSBSSimData.h>
#include <TSBSSimEvent.h>
#include "TSBSDBManager.h"
#include <TSBSSimDataEncoder.h>

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
  TSBSSimDetector::Init();
  if(fDebug>=1)
    cout << "Scintillator detector with UniqueDetID = " << UniqueDetID() << ": TSBSSimScint::Init() " << endl;
  
  // Get the Detector info
  fDetInfo = fDBmanager->GetDetInfo(fName.Data());
  
  double tau = fDetInfo.DigInfo().SPE_Tau();
  double sigma = fDetInfo.DigInfo().SPE_Sigma();
  double tmin = -fDetInfo.DigInfo().GateWidth()/2.0;
  double tmax = +fDetInfo.DigInfo().GateWidth()/2.0;
  double t0 = 0.0;//+fDetInfo.DigInfo().SPE_TransitTime()-fDetInfo.DigInfo().TriggerOffset();
  
  // Get all necessary info to parameterize the PMT pulse shape.
  fSPE = new TSPEModel(fName.Data(), tau, sigma, t0, tmin, tmax);
  
  //Configure the PMT signals array
  fSignals.resize(fDetInfo.NChan());
  for(size_t i_ch = 0; i_ch<fDetInfo.NChan(); i_ch++){
    fSignals[i_ch].SetNpeChargeConv(fDetInfo.DigInfo().NpeChargeConv(i_ch));
  }
}


void TSBSSimScint::LoadEventData(const std::vector<g4sbshitdata*> &evbuffer)
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
    // Only get detector data for Scintillator
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

void TSBSSimScint::LoadAccumulateData(const std::vector<g4sbshitdata*> &evbuffer)
{
  bool signal = false;
  int chan = 0;
  int type = 0;
  double time = 0;
  double data = 0;
  
  // for( const g4sbshitdata *ev: evbuffer) {
  for(std::vector<g4sbshitdata*>::const_iterator it = evbuffer.begin(); it!= evbuffer.end(); ++it ) {
    g4sbshitdata* ev = (*it);
    // Only get detector data for Scintillator
    // new detector ID convetion proposal: UniqueDetID = 10*DetType+DetID
    if(ev->GetDetUniqueID() == UniqueDetID()) {
      signal = (ev->GetData(0)==0);
      chan = ev->GetData(1);
      type = ev->GetData(2);
      time = ev->GetData(3)+fDetInfo.DigInfo().SPE_TransitTime()-fDetInfo.DigInfo().TriggerOffset()+fDetInfo.DigInfo().TriggerJitter() + fTimeZero;//add 
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

void TSBSSimScint::Digitize(TSBSSimEvent &event)
{
  bool any_events = false;
  
  if(fDebug>=3)cout << "TSBSSimScint::Digitize() : Unique Det ID " << UniqueDetID()  << endl;
  
  //TSBSSimEvent::DetectorData data;
  //TSBSSimEvent::SimDetectorData simdata;

  //UInt_t TDCword;
  
  //bool header[8] = {0, 0, 0, 0, 0, 0, 0, 0};//bits ...-31 
  // we ignore (for the moment) all other informations for the ADC.
  //bool channel[7];//Common features between V1190A and 1877: 7 digit channel
  //bool tdc[19];
  //bool trail;
  //hardcoded, but add all ADC/TDC bit coding in DB shall be cumbersome
  //short edgebitpos = 16;
  //short nheaderbits = 32-8-fDetInfo.DigInfo().TDCBits();
  //if(UniqueDetID()==30){
  //  edgebitpos = 26;
  //}
  //short chanfirstbit = fDetInfo.DigInfo().TDCBits()+Short_t(edgebitpos==fDetInfo.DigInfo().TDCBits());

  SimEncoder::adc_data adc_data;
  std::vector<uint32_t> data;
  std::vector<double> simdata;
  short mult = 0; // logical channel multiplier (in case of ADC + TDC together)

  for(size_t m = 0; m < fSignals.size(); m++) {
    //data.fData.clear();
    //simdata.fData.clear();
    fSignals[m].Digitize(fDetInfo.DigInfo(), m);
    if(fSignals[m].Npe() > 0) {
      any_events = true;
      //data.fDetID = UniqueDetID();
      //data.fChannel = m;

      if(fDebug>=4)cout << " = > fSignals[" << m << "].TDCSize() " << fSignals[m].TDCSize() << endl;
      
      /*
      //event.SimDetID.push_back(Short_t(UniqueDetID()));
      event.SimDetChannel[fDetInfo.DetFullName()].push_back(Short_t(m));
      event.SimDetDataType[fDetInfo.DetFullName()].push_back(1);
      event.SimDetNData[fDetInfo.DetFullName()].push_back(1);
      simdata.push_back(fSignals[m].Npe());
      event.SimDetData[fDetInfo.DetFullName()].push_back(simdata);
      event.NSimDetData[fDetInfo.DetFullName()]++;
      simdata.clear();
      
      //event.SimDetID.push_back(Short_t(UniqueDetID()));
      event.SimDetChannel[fDetInfo.DetFullName()].push_back(Short_t(m));
      event.SimDetDataType[fDetInfo.DetFullName()].push_back(2);
      event.SimDetNData[fDetInfo.DetFullName()].push_back(Short_t(fSignals[m].LeadTimesSize()+fSignals[m].TrailTimesSize()));
      for(size_t i = 0; i<fSignals[m].LeadTimesSize(); i++){
	simdata.push_back(fSignals[m].LeadTime(i));
	
	if(fDebug>=3)cout << " leadtime " << i << " = " << fSignals[m].LeadTime(i) << endl;
      }
      for(size_t i = 0; i<fSignals[m].TrailTimesSize(); i++){
	simdata.push_back(fSignals[m].TrailTime(i));
	if(fDebug>=3)cout << " trail time " << i << " = " << fSignals[m].TrailTime(i) << endl;;
      }
      event.SimDetData[fDetInfo.DetFullName()].push_back(simdata);
      event.NSimDetData[fDetInfo.DetFullName()]++;
      simdata.clear();
      */
      
      for(int i_mc = 0; i_mc<fSignals[m].MCHitSize(); i_mc++){
	event.NSimDetHits[fDetInfo.DetFullName()]++;
	event.SimDetChannel[fDetInfo.DetFullName()].push_back(Short_t(m));
	event.SimDetEdep[fDetInfo.DetFullName()].push_back(fSignals[m].MCHitEdep(i_mc));
	event.SimDetNpe[fDetInfo.DetFullName()].push_back(fSignals[m].MCHitNpe(i_mc));
	event.SimDetTime[fDetInfo.DetFullName()].push_back(fSignals[m].MCHitTime(i_mc));
	event.SimDetLeadTime[fDetInfo.DetFullName()].push_back(fSignals[m].MCHitLeadTime(i_mc));
	event.SimDetTrailTime[fDetInfo.DetFullName()].push_back(fSignals[m].MCHitTrailTime(i_mc));
      }
      
      
      //define convention for type:
      // 0: ADC
      // 1: TDC
      // push back a different word for ADC and TDC ?
      mult = 0;
      //if(fDetInfo.DigInfo().ADCBits()>0 && fDetInfo.DigInfo().ADCConversion()>0){
      if(fEncoderADC) {
        adc_data.integral=fSignals[m].ADC();
        fEncoderADC->EncodeADC(adc_data,fEncBuffer,fNEncBufferWords);
        CopyEncodedData(fEncoderADC,mult++,data);//.fData);

	//simdata.fData.clear();
      }
      
      // Fill ADC
      for(int i = 0; i<data.size(); i++){
	event.NDetHits[fDetInfo.DetFullName()]++;
	event.DetChannel[fDetInfo.DetFullName()].push_back(Short_t(m));
	event.DetDataWord[fDetInfo.DetFullName()].push_back(data.at(i));
	event.DetADC[fDetInfo.DetFullName()].push_back(fSignals[m].ADC());
	event.DetTDC[fDetInfo.DetFullName()].push_back(-1);
      }
      data.clear();
      
      
      if(fEncoderTDC) {
        fEncoderTDC->EncodeTDC(fSignals[m].TDCData(),fEncBuffer,
            fNEncBufferWords);
        CopyEncodedData(fEncoderTDC,mult++,data);//.fData);
      }
      // Fill TDC
      for(int i = 0; i<data.size(); i++){
	event.NDetHits[fDetInfo.DetFullName()]++;
	event.DetChannel[fDetInfo.DetFullName()].push_back(Short_t(m));
	event.DetDataWord[fDetInfo.DetFullName()].push_back(data.at(i));
	event.DetADC[fDetInfo.DetFullName()].push_back(-1);
	event.DetTDC[fDetInfo.DetFullName()].push_back(fSignals[m].TDC(i));
      }
      data.clear();
      
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
void TSBSSimScint::Clear(Option_t*)
{
  for(size_t i = 0; i < fSignals.size(); i++ ) {
    fSignals[i].Clear();
  }
}

ClassImp(TSBSSimScint) // Implements TSBSSimScint
