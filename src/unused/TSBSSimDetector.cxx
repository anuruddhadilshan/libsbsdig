#include "TSBSSimDetector.h"
#include "TSBSDBManager.h"
#include "TSBSSimData.h"
#include "TSBSSimEvent.h"


int TSBSSimDetector::fEvNum = 0;

TSBSSimDetector::TSBSSimDetector() : fEncoderADC(0),
  fEncoderTDC(0), fHasData(false)
{
  fDBmanager = TSBSDBManager::GetInstance();
}

void TSBSSimDetector::Init()
{
  fDetInfo = fDBmanager->GetDetInfo(fName.Data());
  // Find out what encoders this detector has available
  fEncoderADC = fDetInfo.DigInfo().GetEncoderADC();
  fEncoderTDC = fDetInfo.DigInfo().GetEncoderTDC();
}

TSBSSimDetector::~TSBSSimDetector()
{
  fDBmanager->Delete();
}

void TSBSSimDetector::LoadMCTrackData(const std::vector<g4sbsgendata*> &evbuffer, 
				      TSBSSimEvent &event)
{
  for(std::vector<g4sbsgendata*>::const_iterator it = evbuffer.begin(); it!= evbuffer.end(); ++it ) {
    g4sbsgendata* ev = (*it);
    if(ev->GetDetUniqueID() == UniqueDetID()) {
      if(fDetInfo.DetFullName().find("gem")!=std::string::npos){//GEMs
	event.fMCTrackOutData[fDetInfo.DetFullName()].fNMCTracks++;
	event.fMCTrackOutData[fDetInfo.DetFullName()].fMCTrackSource.push_back(ev->GetSource());
	event.fMCTrackOutData[fDetInfo.DetFullName()].fMCTrackTRID.push_back(ev->GetTRID());
	event.fMCTrackOutData[fDetInfo.DetFullName()].fMCTrackPID.push_back(ev->GetPID());
	
	event.fMCTrackOutData[fDetInfo.DetFullName()].fMCTrackX.push_back(ev->GetX());
	event.fMCTrackOutData[fDetInfo.DetFullName()].fMCTrackY.push_back(ev->GetY());
	event.fMCTrackOutData[fDetInfo.DetFullName()].fMCTrackT.push_back(ev->GetT());
	event.fMCTrackOutData[fDetInfo.DetFullName()].fMCTrackP.push_back(ev->GetP());
	event.fMCTrackOutData[fDetInfo.DetFullName()].fMCTrackdX.push_back(ev->GetdX());
	event.fMCTrackOutData[fDetInfo.DetFullName()].fMCTrackdY.push_back(ev->GetdY());
	event.fMCTrackOutData[fDetInfo.DetFullName()].fMCTrackXv.push_back(ev->GetVertexAtTarget().X());
	event.fMCTrackOutData[fDetInfo.DetFullName()].fMCTrackYv.push_back(ev->GetVertexAtTarget().Y());
	event.fMCTrackOutData[fDetInfo.DetFullName()].fMCTrackZv.push_back(ev->GetVertexAtTarget().Z());
	event.fMCTrackOutData[fDetInfo.DetFullName()].fMCTrackPXv.push_back(ev->GetMomentumAtTarget().X());
	event.fMCTrackOutData[fDetInfo.DetFullName()].fMCTrackPYv.push_back(ev->GetMomentumAtTarget().Y());
	event.fMCTrackOutData[fDetInfo.DetFullName()].fMCTrackPZv.push_back(ev->GetMomentumAtTarget().Z());
	event.fMCTrackOutData[fDetInfo.DetFullName()].fMCTrackWeight.push_back(ev->GetWeight());
      }else{
	event.fTrackMCHitOutData[fDetInfo.DetFullName()].fNTrackMCHits++;
	event.fTrackMCHitOutData[fDetInfo.DetFullName()].fTrackMCSource.push_back(ev->GetSource());
	event.fTrackMCHitOutData[fDetInfo.DetFullName()].fTrackMCTRID.push_back(ev->GetTRID());
	event.fTrackMCHitOutData[fDetInfo.DetFullName()].fTrackMCPID.push_back(ev->GetPID());
	event.fTrackMCHitOutData[fDetInfo.DetFullName()].fTrackMCXhit.push_back(ev->GetX());
	event.fTrackMCHitOutData[fDetInfo.DetFullName()].fTrackMCYhit.push_back(ev->GetY());
	event.fTrackMCHitOutData[fDetInfo.DetFullName()].fTrackMCThit.push_back(ev->GetT()+fTimeZero);
	event.fTrackMCHitOutData[fDetInfo.DetFullName()].fTrackMCE.push_back(ev->GetE());
	event.fTrackMCHitOutData[fDetInfo.DetFullName()].fTrackMCWeight.push_back(ev->GetWeight());
      }

    }
  }//end loop on g4sbsgendata
}

void TSBSSimDetector::CopyEncodedData(SBSSimDataEncoder *enc,
    unsigned short mult, std::vector<unsigned int> &dat)
{
  dat.push_back(SBSSimDataEncoder::EncodeHeader(enc->GetId(),
        mult,fNEncBufferWords));
  for(unsigned short n = 0; n < fNEncBufferWords; n++) {
    dat.push_back(fEncBuffer[n]);
  }
  fNEncBufferWords = 0;
}

ClassImp(TSBSSimDetector)
/*
//
// Class TSPEModel
//
TSPEModel::TSPEModel(DigInfo diginfo, const char* detname):
  fDigInfo(diginfo)//, qe(1.602e-19), unit(1e-9)
{
  fScale = fDigInfo.fROImpedance*qe/spe_unit;
  
  if(fDigInfo.fGain.size()==1){
    fScale*= fDigInfo.fGain[0];
  }
  
  //fStartTime = ;
  //start_t = -12.5;
  double mint = -fDigInfo.fGateWidth/2.0;
  double maxt = +fDigInfo.fGateWidth/2.0;
  // test values
  double tau = fDigInfo.fSPEtau;
  double sig = fDigInfo.fSPEsig;
  double t0 = fDigInfo.fSPEtransittime-fDigInfo.fTriggerOffset;
  
  TF1 fFunc1(Form("fFunc1%s",detname),
	     TString::Format("TMath::Max(0.,(x/%g)*TMath::Exp(-x/(%g)))",
			     tau*tau,tau),
	     mint,maxt);
  TF1 fFunc2(Form("fFunc2%s",detname),
	     TString::Format("%g*TMath::Exp(-((x-%g)**2)/(%g))",
			     1./TMath::Sqrt(2*TMath::Pi()*sig),t0,sig*sig),
	     mint,maxt);
  TF1Convolution fConvolution(&fFunc1,&fFunc2);
  
  model = new TF1(Form("fSignal%s",detname),fConvolution,mint,maxt, fConvolution.GetNpar());
}

double TSPEModel::Eval(double t, int chan)
{
  if(fDigInfo.fGain.size()>1){
    if(fDigInfo.fGain.size()<=chan){
      cout << "warning: requested channel number " << chan << "larger than number of channel size " << fDigInfo.fGain.size() << " check database ! " << endl;
      exit(-1);
    }
    fScale*= fDigInfo.fGain[chan];
  }
  
  return fScale*model->Eval(t);
  //return model->Eval(t);
  //return 1.0;
}

//
// Class TPMTSignal
//

TPMTSignal::TPMTSignal()
  : fNpe(0), fADC(0)
{  
  leadtimes.clear();
  trailtimes.clear();
}

void TPMTSignal::Fill(TSPEModel *model,double t, double toffset)
{
  // int start_bin = 0;
  // if( mint > t )
  //   toffset -= (mint-t);
  // else
  //   start_bin = (t-mint)/dx_raw;

  // if(start_bin > nbins_raw)
  //   return; // Way outside our window anyways

  // // Now digitize this guy into the raw_bin (scope)
  //double tt = model->GetStartTime-toffset;
  // //std::cout << "t=" << t << ", tt=" << tt << std::endl;
  // for(int bin = start_bin; bin < nbins_raw; bin++) {
  //   samples_raw[bin] += model->Eval(tt);
  //   tt += dx_raw;
  // }
  // npe++;
}

void TPMTSignal::Clear()
{
  fNpe = 0;
  fADC = 0;
  leadtimes.clear();
  trailtimes.clear();
}
*/
ClassImp(TSPEModel) // Implements TSPEModel
ClassImp(TPMTSignal) // Implements TPMTSignal
