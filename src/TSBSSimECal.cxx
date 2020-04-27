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
    fSignals[i_ch].Clear();
    fSignals[i_ch].SetNpeChargeConv(fDetInfo.DigInfo().NpeChargeConv(i_ch));
  }
}


void TSBSSimECal::LoadEventData(const std::vector<g4sbshitdata*> &evbuffer)
{
  // New event - we want to clear all arrays before storing some new data
  Clear();
  LoadAccumulateData(evbuffer);
}

void TSBSSimECal::LoadAccumulateData(const std::vector<g4sbshitdata*> &evbuffer)
{
  bool signal = false;
  int chan = 0;
  //int type = 0;
  double time = 0;
  double npe = 0;
  double edep = 0;
  
  // for( const g4sbshitdata *ev: evbuffer) {
  for(std::vector<g4sbshitdata*>::const_iterator it = evbuffer.begin(); it!= evbuffer.end(); ++it ) {
    g4sbshitdata* ev = (*it);
    // Only get detector data for ECal
    // new detector ID convetion proposal: UniqueDetID = 10*DetType+DetID
    if(ev->GetDetUniqueID() == UniqueDetID()) {
      signal = ev->GetData(0);
      chan = ev->GetData(1);
      //type = ev->GetData(2);
      time = ev->GetData(2)+fDetInfo.DigInfo().SPE_TransitTime()-fDetInfo.DigInfo().TriggerOffset() + fTimeZero;//+fDetInfo.DigInfo().TriggerJitter()
      npe = ev->GetData(3);
      edep = ev->GetData(4);
      
      if(fabs(time)>fDetInfo.DigInfo().GateWidth()/2.0)continue;
      
      if(fDebug>=3)
	cout << "Detector " << UniqueDetID() << " chan = " << chan << " Evt Time: "<< ev->GetData(3) << " " << time << endl;
      
      if(!fSignals[chan].check_vec_size())cout << "Warning: *before filling* MC vector container sizes for det " << fDetInfo.DetFullName() << " chan " << chan << " dont check out!" << endl;
      fSignals[chan].Fill(fSPE, npe, fDetInfo.DigInfo().Threshold(chan), time, signal);
      fSignals[chan].AddSumEdep(edep);
      if(!fSignals[chan].check_vec_size())cout << "Warning: Size of MC info containers for " << fDetInfo.DetFullName() << " chan " << chan << " don't check out!!!" << endl;
    }
  }//end loop on evbuffer
}

void TSBSSimECal::Digitize(TSBSSimEvent &event)
{
  bool any_events = false;

  if(fDebug>=2)cout << "TSBSSimECal::Digitize() : Unique Det ID " << UniqueDetID() << " signal size = " << fSignals.size() << endl;
  
  int mult = 0;
  SimEncoder::adc_data adc_data;
  std::vector<uint32_t> data;
  UInt_t tdcval;
  //std::vector<double> simdata;
  for(size_t m = 0; m < fSignals.size(); m++) {
    //data.fData.clear();
    //simdata.fData.clear();
    data.clear();
    //simdata.clear();
    if(fDebug>=3)cout << "digitize channel " << m << endl;
    fSignals[m].Digitize(fDetInfo.DigInfo(), m);
    if(fSignals[m].Npe() >= 0) {// we want the blocks with no signal, for the pedestals
      any_events = true;
      //data.fDetID = UniqueDetID();
      //data.fChannel = m;
      
      if(fDebug>=4)cout << " = > fSignals[" << m << "].ADC() " << fSignals[m].ADC() << endl;
      
      if(fDebug>=4){
	cout << fDetInfo.DetFullName() << ": check MC vec size" << endl;
	fSignals[m].check_vec_size();
      }
      for(uint i_mc = 0; i_mc<fSignals[m].MCHitSize(); i_mc++){
	event.fSimHitMCOutData[fDetInfo.DetFullName()].fNSimHits++;
	event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimSource.push_back(fSignals[m].MCHitSource(i_mc));
	//event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimTRID.push_back(fSignals[m].MCHitSource(i_mc));//dummy values for the moment
	//event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimPID.push_back(fSignals[m].MCHitSource(i_mc));//dummy values for the moment
	event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimChannel.push_back(Short_t(m));
	event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimEdep.push_back(fSignals[m].MCHitEdep(i_mc));
	event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimNpe.push_back(fSignals[m].MCHitNpe(i_mc));
	event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimTime.push_back(fSignals[m].MCHitTime(i_mc));
	
	if(fEncoderTDC){
	  if(fDebug>=4)cout << "MC times " << fSignals[m].MCHitLeadTime(i_mc) << " " << fSignals[m].MCHitTrailTime(i_mc) << endl;
	  event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimLeadTime.push_back(fSignals[m].MCHitLeadTime(i_mc));
	  event.fSimHitMCOutData[fDetInfo.DetFullName()].fSimTrailTime.push_back(fSignals[m].MCHitTrailTime(i_mc));
	}
      }
      
      //define convention for type:
      // 0: ADC
      // 1: TDC
      mult = 0;
      
      // Fill ADC
      if(fEncoderADC) {
	if(fDebug>=4)cout << " dig ADC " << fSignals[m].ADC() << endl;
        adc_data.integral=fSignals[m].ADC();
        fEncoderADC->EncodeADC(adc_data,fEncBuffer,fNEncBufferWords);
        CopyEncodedData(fEncoderADC,mult++,data);
	
	for(uint i = 0; i<data.size(); i++){
	  event.fSimDigOutData[fDetInfo.DetFullName()].fNHits++;
	  event.fSimDigOutData[fDetInfo.DetFullName()].fDataWord.push_back(data.at(i));
	  
	  if(i==0){//header
	    event.fSimDigOutData[fDetInfo.DetFullName()].fChannel.push_back(-1000);
	    event.fSimDigOutData[fDetInfo.DetFullName()].fADC.push_back(-1000000);
	    if(fEncoderTDC){
	      event.fSimDigOutData[fDetInfo.DetFullName()].fTDC_L.push_back(-1000000);
	      event.fSimDigOutData[fDetInfo.DetFullName()].fTDC_T.push_back(-1000000);
	    }
	  }else{
	    event.fSimDigOutData[fDetInfo.DetFullName()].fChannel.push_back(Short_t(m));
	    event.fSimDigOutData[fDetInfo.DetFullName()].fADC.push_back(fSignals[m].ADC()-fDetInfo.DigInfo().Pedestal(m));
	    if(fEncoderTDC){
	      event.fSimDigOutData[fDetInfo.DetFullName()].fTDC_L.push_back(-1000000);
	      event.fSimDigOutData[fDetInfo.DetFullName()].fTDC_T.push_back(-1000000);
	    }
	  }
	}
	data.clear();
      }
      
      // Fill TDC 
      if(fEncoderTDC && fSignals[m].TDCSize()) {
 	if(fDebug>=4)cout << fSignals[m].TDC(0) << " " << fSignals[m].TDC(1) << endl;
	fEncoderTDC->EncodeTDC(fSignals[m].TDCData(),fEncBuffer,
            fNEncBufferWords);
        CopyEncodedData(fEncoderTDC,mult++,data);
	
	for(uint i = 0; i<data.size(); i++){
	  event.fSimDigOutData[fDetInfo.DetFullName()].fNHits++;
	  event.fSimDigOutData[fDetInfo.DetFullName()].fDataWord.push_back(data.at(i));
	  if(fEncoderADC)event.fSimDigOutData[fDetInfo.DetFullName()].fADC.push_back(-1000000);

	  if(i==0){//header
	    event.fSimDigOutData[fDetInfo.DetFullName()].fChannel.push_back(-1000);
	    event.fSimDigOutData[fDetInfo.DetFullName()].fTDC_L.push_back(-1000000);
	    event.fSimDigOutData[fDetInfo.DetFullName()].fTDC_T.push_back(-1000000);
	  }else{
	    event.fSimDigOutData[fDetInfo.DetFullName()].fNHits++;
	    event.fSimDigOutData[fDetInfo.DetFullName()].fChannel.push_back(Short_t(m));
	    event.fSimDigOutData[fDetInfo.DetFullName()].fDataWord.push_back(data.at(i));
	    if(fEncoderADC)event.fSimDigOutData[fDetInfo.DetFullName()].fADC.push_back(-1000000);
	    if( fSignals[m].TDC(i-1) & ( 1 << (31) ) ){
	      tdcval = fSignals[m].TDC(i-1);
	      if(fDebug>=4)cout << " T: " << tdcval << " => ";
	      tdcval ^= ( -0 ^ tdcval) & ( 1 << (31) );
	      if(fDebug>=4)cout << tdcval << endl;
	      event.fSimDigOutData[fDetInfo.DetFullName()].fTDC_L.push_back(-1000000);
	      event.fSimDigOutData[fDetInfo.DetFullName()].fTDC_T.push_back(tdcval-1.e3/fDetInfo.DigInfo().TDCConversion());
	    }else{
	      event.fSimDigOutData[fDetInfo.DetFullName()].fTDC_L.push_back(fSignals[m].TDC(i-1)-1.e3/fDetInfo.DigInfo().TDCConversion());
	      if(fDebug>=4)cout << " L: " << fSignals[m].TDC(i-1) << endl;
	      event.fSimDigOutData[fDetInfo.DetFullName()].fTDC_T.push_back(-1000000);
	    }
	  }
	}
	data.clear();
      }
    }//end if fSignals.Npe
  }//end loop on signals
  //if(fDebug>=2)
  if(!event.fSimDigOutData[fDetInfo.DetFullName()].CheckSize(bool(fEncoderADC), bool(fEncoderTDC), fDebug>=2)){
    cout << "Warning: output vectors for" << fDetInfo.DetFullName() << " don't have the same size! (any events ?" << any_events << ")" << endl;
  }
  if(!event.fSimHitMCOutData[fDetInfo.DetFullName()].CheckSize(true, true, bool(fEncoderTDC), fDebug>=0)){
    cout << "Warning:  output MC vectors for " << fDetInfo.DetFullName() << " don't have the same size! " << endl;
  }
  SetHasDataFlag(any_events);
}

// Clear signals in array
void TSBSSimECal::Clear(Option_t*)
{
  for(size_t i = 0; i < fDetInfo.NChan(); i++ ) {
    fSignals[i].Clear();
  }
}

ClassImp(TSBSSimECal) // Implements TSBSSimECal
