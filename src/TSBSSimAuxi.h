#ifndef _TSBSSIMAUXI_H
#define _TSBSSIMAUXI_H

#include <vector>
#include <map>
#include "g4sbs_types.h"
#include "THaAnalysisObject.h"
#include "g4sbs_types.h"
#include "TF1.h"
//#include "TF1Convolution.h"
#include "TH1D.h"
#include "TRandom3.h"

//
// Classes for DB information
//
//__________________________________
class TSignalInfo : public TObject {
 public:
  TSignalInfo() {};
 TSignalInfo(int apid, int atid):fPID(apid), fTID(atid) {};
  ~TSignalInfo() {};
  
  double PID(){return fPID;};
  double TID(){return fTID;};

  void SetPID(int apid){fPID = apid;};
  void SetTID(int atid){fTID = atid;};
  
 private:
  int fPID;
  int fTID;
  
  ClassDef(TSignalInfo, 1);
};

//__________________________________
class TSpectroInfo : public TObject{
 public:
  TSpectroInfo() {};
  TSpectroInfo(int ndets):fNDets(ndets) {};
  ~TSpectroInfo() {
    fMCsignalInfo.clear();
    fDetNames.clear();
  };
  
  double MCAngle(){return fMCangle;};
  UInt_t NDets(){return fNDets;};
  UInt_t SignalInfoSize(){return fMCsignalInfo.size();};
  
  TSignalInfo MCSignalInfo(int i){return fMCsignalInfo.at(i);};
  std::string DetName(int i){return fDetNames.at(i);};
  
  void SetMCAngle(double ang){fMCangle = ang;};
  void SetNDets(UInt_t ndets){fNDets = ndets;};
  void AddMCSignalInfo(TSignalInfo siginfo){fMCsignalInfo.push_back(siginfo);};
  void AddDetName(std::string detname){fDetNames.push_back(detname);};
  
 private:
  double fMCangle;
  UInt_t fNDets;
  std::vector<TSignalInfo> fMCsignalInfo;
  std::vector<std::string> fDetNames;
  
  ClassDef(TSpectroInfo, 1);
};

//______________________________
class TGeoInfo : public TObject{
 public:
  TGeoInfo() {};
  ~TGeoInfo() {};
  
  UInt_t NRows(){return fNrows;};
  UInt_t NCols(){return fNcols;};
  double XSize(){return fXsize;};
  double YSize(){return fYsize;};
  double ZPos(){return fZpos;};
  double XOffset(){return fXoffset;};
  double YOffset(){return fYoffset;};
  
  void SetNRows(UInt_t nrows){fNrows = nrows;};
  void SetNCols(UInt_t ncols){fNcols = ncols;};
  void SetXSize(double xsize){fXsize = xsize;};
  void SetYSize(double ysize){fYsize = ysize;};
  void SetZPos(double zpos){fZpos = zpos;};
  void SetXOffset(double xoffset){fXoffset = xoffset;};
  void SetYOffset(double yoffset){fYoffset = yoffset;};
  
 private:
  UInt_t  fNrows;      // number of rows
  UInt_t  fNcols;      // number of columns
  double  fXsize;      // detector X size (in transport coordinates)
  double  fYsize;      // detector Y size (in transport coordinates)
  double  fZpos;       // detector position on spectrometer axis
  double  fXoffset;    // detector X offset - handy for detectors with many modules (in transport coordinates)
  double  fYoffset;    // detector Y offset - handy for detectors with many modules (in transport coordinates)
  
  ClassDef(TGeoInfo, 1);
};

//______________________________
class TDigInfo : public TObject{
 public:
  TDigInfo();
  ~TDigInfo();

  double ROImpedance(){return fROimpedance;};
  double ADCConversion(){return fADCconversion;};
  Int_t  ADCBits(){return fADCbits;};
  double TDCConversion(){return fTDCconversion;};
  Int_t  TDCBits(){return fTDCbits;};
  UInt_t GainSize(){return fGain.size();};
  double Gain(UInt_t chan);
  UInt_t PedestalSize(){return fPedestal.size();};
  double Pedestal(UInt_t chan);
  UInt_t PedestalNoiseSize(){return fPedNoise.size();};
  double PedestalNoise(UInt_t chan);
  double GenPedestal(UInt_t chan){return fRN->Gaus(Pedestal(chan), PedestalNoise(chan));};
  UInt_t ThresholdSize(){return fThreshold.size();};
  double Threshold(UInt_t chan);
  double TriggerJitter(){return fTriggerJitter;};
  double TriggerOffset(){return fTriggerOffset;};
  double GateWidth(){return fGateWidth;};
  double SPE_Tau(){return fSPE_tau;};
  double SPE_Sigma(){return fSPE_sigma;};
  double SPE_TransitTime(){return fSPE_transittime;};
    
  double NpeChargeConv(int chan){return Gain(chan)*qe;};//charge in Coulomb
  
  void SetROImpedance(double roimp){fROimpedance = roimp;};
  void SetADCConversion(double adcconv){fADCconversion = adcconv;};
  void SetADCBits(int adcbits){fADCbits = adcbits;};
  void SetTDCConversion(double tdcconv){fTDCconversion = tdcconv;};
  void SetTDCBits(int tdcbits){fTDCbits = tdcbits;};
  void AddGain(double gain){fGain.push_back(gain);};
  void AddPedestal(double ped){fPedestal.push_back(ped);};
  void AddPedestalNoise(double pednoise){fPedNoise.push_back(pednoise);};
  void AddThreshold(double thr){fThreshold.push_back(thr);};
  void SetTriggerJitter(double triggerjitter){fTriggerJitter = triggerjitter;};
  void SetTriggerOffset(double triggeroffset){fTriggerOffset = triggeroffset;};
  void SetGateWidth(double gatewidth){fGateWidth = gatewidth;};
  void SetSPE_Tau(double spe_tau){fSPE_tau = spe_tau;};
  void SetSPE_Sigma(double spe_sig){fSPE_sigma = spe_sig;};
  void SetSPE_TransitTime(double spe_transit){fSPE_transittime = spe_transit;};

 private:
  double  fROimpedance;     // readout impedance
  double  fADCconversion;    // charge/ADC channel conversion
  Int_t   fADCbits;           // number of bits in ADC
  double  fTDCconversion;      // time/TDC channel conversion
  Int_t   fTDCbits;             // number of bits in ADC
  std::vector<double>  fGain;    // Gain 
  std::vector<double>  fPedestal; // Pedestal value (adc value)
  std::vector<double>  fPedNoise;  // Pedestal noise (adc value)
  std::vector<double>  fThreshold; //Channel threshold
  double  fTriggerJitter;         // trigger jitter (ns)
  double  fTriggerOffset;        // trigger offset (ns)
  double  fGateWidth;           // gate width (ns)
  double  fSPE_tau;            // tau param for SPE
  double  fSPE_sigma;         // sigma param for SPE
  double  fSPE_transittime;  // pmt transit time param for SPE
  
  TRandom3* fRN;
  
  ClassDef(TDigInfo, 1);
};

//______________________________
//
// A "digital" module in a "digital crate" with corresponding
// slot#, and number of channels, etc...
class TDigSlot : public TObject {
public:
  TDigSlot();
  TDigSlot(Int_t crate, Int_t slot, Int_t lo, Int_t hi);
  virtual ~TDigSlot();

  Int_t GetCrate()  { return fCrate; }
  Int_t GetSlot()   { return fSlot; }
  Int_t GetNchan()  { return fNchan; }
  Int_t GetChanLo() { return fChanLo; }
  Int_t GetChanHi() { return fChanHi; }

  // Get corresponding channel number
  Int_t GetChanNumber(Int_t ch);

private:
  Int_t fCrate;
  Int_t fSlot;
  Int_t fNchan;
  Int_t fChanLo;
  Int_t fChanHi;

  ClassDef(TDigSlot,0);
};

//______________________________
//
// Contains information about the logical channel, such as "digital" slot
// crate, slot#, channel number in the current TDigSlot, etc..
struct TDigChannelInfo {
  int ch;
  int slot;
  int crate;
};


//______________________________
class TDetInfo : public TObject{
 public:
  TDetInfo();
  TDetInfo(const std::string detname);
  ~TDetInfo();
  
  std::string DetName(){return fDetName;};
  det_type DetType(){return fDetType;};
  UInt_t   NChan(){return fNchan;};
  UInt_t   NLogChan(){return fNlogChan;};
  UInt_t   ChanPerSlot(){return fChanPerSlot;};
  UInt_t   SlotPerCrate(){return fSlotPerCrate;};
  UInt_t   FirstSlot(){return fFirstSlot;};
  UInt_t   FirstCrate(){return fFirstCrate;};
  UInt_t   NPlanes(){return fNplanes;};
  UInt_t   NModulesSize(){return fNmodules.size();};
  UInt_t   NModules(int i){return fNmodules.at(i);};
  
  int GeoInfoSize(){return fGeoInfo.size();};
  TGeoInfo GeoInfo(int i){return fGeoInfo.at(i);};
  TDigInfo DigInfo(){return fDigInfo;};
  
  void SetDetName(std::string detname){fDetName = detname;};
  void SetDetType(det_type type){fDetType = type;};
  void SetNChan(UInt_t nchan){fNchan = nchan;};
  void SetNLogChan(UInt_t nchan){fNlogChan = nchan;};
  void SetChanPerSlot(UInt_t chanperslot){fChanPerSlot = chanperslot;};
  void SetSlotPerCrate(UInt_t slotpercrate){fSlotPerCrate = slotpercrate;};
  void SetFirstSlot(UInt_t firstslot){fFirstSlot = firstslot;};
  void SetFirstCrate(UInt_t firstcrate){fFirstCrate = firstcrate;};
  void SetNPlanes(UInt_t nplanes){fNplanes = nplanes;};
  void AddNModules(UInt_t nmodules){fNmodules.push_back(nmodules);};
  
  void AddGeoInfo(TGeoInfo geoinfo){fGeoInfo.push_back(geoinfo);};
  void SetDigInfo(TDigInfo diginfo){fDigInfo = diginfo;};

  Int_t AddSlot(Int_t crate, Int_t slot, Int_t lo, Int_t hi);
  TDigChannelInfo FindLogicalChannelSlot(Int_t lch);
  void LoadChannelMap(std::vector<int> chanmap);

 private:
  std::string fDetName;      // Detector name
  det_type     fDetType;      // DetectorType
  UInt_t        fNchan;        // Total number of channels over all detector
  UInt_t         fNlogChan;        // Total number of *logical* channels over all detector
  UInt_t          fChanPerSlot;  // Number of channels per slot
  UInt_t           fSlotPerCrate; // Number of slots per crate
  UInt_t            fFirstSlot;  // first slot in detmap
  UInt_t             fFirstCrate; // first crate in detmpa
  UInt_t              fNplanes;      // Number of planes // useful e.g. GEM, CDet
  std::vector<UInt_t>  fNmodules;     // Number of modules per plane // useful e.g. GEM, CDet
  
  std::vector<TGeoInfo> fGeoInfo;
  TDigInfo fDigInfo;
  std::vector<TDigSlot> fModSlots;
  std::map<int,std::pair<int,int> > fDetMap;
  
  ClassDef(TDetInfo, 1);
};

//
// classes for signal digitization
//
//________________________________
class TSPEModel : public TObject {
 public:
  TSPEModel(const char* detname, double tau, double sigma, double t0 = 0, double tmin = -100, double tmax = +100);
  double Eval(double t){return fPulseModel->Eval(t);};
  bool   PulseOverThr(double charge, double thr);
  void   FindLeadTrailTime(double charge, double thr, double &t_lead, double &t_trail);
  
 private:
  TF1 *fPulseModel;// We'll have to ditch that stuff... *sigh*
  TH1D *fPulseHisto;//At least we'll have to setup and use one per detector
  double GetHistoX(double y, double x1, double x2);
  //void  BuildHisto(double tau, double sigma);
  
  ClassDef(TSPEModel,1);
};

//_________________________________
class TPMTSignal : public TObject {
 public:
  TPMTSignal();
  TPMTSignal(double npechargeconv);
  void Fill(TSPEModel *model, int npe, double thr, double evttime, bool signal);
  void Digitize(TDigInfo diginfo, int chan);
  void Clear();
  ~TPMTSignal(){Clear();};
  
  void AddSumEdep(double edep){fSumEdep+= edep;};
  void SetNpeChargeConv(double npechargeconv){fNpeChargeConv = npechargeconv;};
    
  double SumEdep(){return fSumEdep;};
  UInt_t Npe(){return fNpe;};
  double Charge(){return fNpe*fNpeChargeConv;};
  Int_t  ADC(){return fADC;};

  double EventTime(){return fEventTime;};
  UInt_t LeadTimesSize(){return fLeadTimes.size();};
  double LeadTime(int i){return fLeadTimes.at(i);};
  UInt_t TrailTimesSize(){return fTrailTimes.size();};
  double TrailTime(int i){return fTrailTimes.at(i);};
  UInt_t TDCSize(){return fTDCs.size();};
  UInt_t TDC(int i){return fTDCs.at(i);};
  
  
 private:
  double fSumEdep;//Not forced to use it for everything
  UInt_t fNpe;
  double fNpeChargeConv;
  Int_t  fADC;// One unique ADC value ?

  double fEventTime;
  //TDCs: multiple values possible.
  std::vector<double> fLeadTimes;
  std::vector<double> fTrailTimes;
  std::vector<UInt_t> fTDCs;
  
  ClassDef(TPMTSignal,1);
};


/*
//
// That class was too complicated... about to scrap...
//
class TNPEModel : public TObject {
 public:
  TNPEModel(DigInfo diginfo, const char* detname, int npe = 1);
  //TNPEModel(const char* detname, double);
  double GetNpe(){return fNpe;};
  double GetCharge(int chan);
  void   SetNpe(double npe){fNpe = npe;};
  double GetStartTime(){return fStartTime;};
  void   SetStartTime(double t){fStartTime = t;};
  // double GetADCconversion(){return fDigInfo.fADCconversion;};
  // double GetTDCconversion(){return fDigInfo.fTDCconversion;};

  double Eval(int chan, double t);
  void   FindLeadTrailTime(int chan, double &t_lead, double &t_trail);
  bool   PulseOverThr(int chan = 0);
    
  ~TNPEModel(){};
  
 private:
  DigInfo fDigInfo;// too much info/memory - better off not having it - 
  // specially since it will be used in detector classes which have the info.
  TF1 *fModel;
  double fScale;// Npe charge scale.
  double fNpe;
  double fStartTime;
  
  / *
  double gain_pmt;
  double resistance; //ohm
  double qe; //
  double unit;
  double scale;
  TF1 *fFunc1;
  TF1 *fFunc2;
  TF1Convolution *fConvolution;
  double mint;
  double start_t;
  double maxt;
  double tau;
  double sig;
  double t0;
  * /
  ClassDef(TNPEModel,1);
};
*/


#endif // _TSBSSIMAUXI_H

