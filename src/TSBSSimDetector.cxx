#include "TSBSSimDetector.h"
#include "TSBSDBManager.h"

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

void TSBSSimDetector::CopyEncodedData(TSBSSimDataEncoder *enc,
    unsigned short mult, std::vector<unsigned int> &dat)
{
  dat.push_back(TSBSSimDataEncoder::EncodeHeader(enc->GetId(),
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
