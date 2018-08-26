#ifndef _TSBSSIMAUXI_H
#define _TSBSSIMAUXI_H

#include <vector>
#include "g4sbs_types.h"
#include "THaAnalysisObject.h"
#include "g4sbs_types.h"
#include "TF1.h"
#include "TF1Convolution.h"

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
  int    NDets(){return fNDets;};
  
  TSignalInfo MCSignalInfo(int i){return fMCsignalInfo.at(i);};
  std::string DetName(int i){return fDetNames.at(i);};
  
  void SetMCAngle(double ang){fMCangle = ang;};
  void SetNDets(double ndets){fNDets = ndets;};
  void AddMCSignalInfo(TSignalInfo siginfo){fMCsignalInfo.push_back(siginfo);};
  void AddDetName(std::string detname){fDetNames.push_back(detname);};
  
 private:
  double fMCangle;
  int fNDets;
  std::vector<TSignalInfo> fMCsignalInfo;
  std::vector<std::string> fDetNames;
  
  ClassDef(TSpectroInfo, 1);
};

//______________________________
class TGeoInfo : public TObject{
 public:
  TGeoInfo() {};
  ~TGeoInfo() {};
  
  int    NRows(){return fNrows;};
  int    NCols(){return fNcols;};
  double XSize(){return fXsize;};
  double YSize(){return fYsize;};
  double ZPos(){return fZpos;};
  double XOffset(){return fXoffset;};
  double YOffset(){return fYoffset;};
  
  void SetNRows(int nrows){fNrows = nrows;};
  void SetNCols(int ncols){fNcols = ncols;};
  void SetXSize(double xsize){fXsize = xsize;};
  void SetYSize(double ysize){fYsize = ysize;};
  void SetZPos(double zpos){fZpos = zpos;};
  void SetXOffset(double xoffset){fXoffset = xoffset;};
  void SetYOffset(double yoffset){fYoffset = yoffset;};
  
 private:
  int     fNrows;      // number of rows
  int     fNcols;      // number of columns
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
  int    ADCBits(){return fADCbits;};
  double TDCConversion(){return fTDCconversion;};
  int    TDCBits(){return fTDCbits;};
  int    GainSize(){return fGain.size();};
  double Gain(uint chan);
  int    PedestalSize(){return fPedestal.size();};
  double Pedestal(uint chan);
  int    PedestalNoiseSize(){return fPedNoise.size();};
  double PedestalNoise(uint chan);
  int    ThresholdSize(){return fThreshold.size();};
  double Threshold(uint chan);
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
  int     fADCbits;           // number of bits in ADC
  double  fTDCconversion;      // time/TDC channel conversion
  int     fTDCbits;             // number of bits in ADC
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
  
  ClassDef(TDigInfo, 1);
};

//______________________________
class TDetInfo : public TObject{
 public:
  TDetInfo();
  TDetInfo(const std::string detname);
  ~TDetInfo();
  
  std::string DetName(){return fDetName;};
  det_type   DetType(){return fDetType;};
  int       NChan(){return fNchan;};
  int      ChanPerSlot(){return fChanPerSlot;};
  int     SlotPerCrate(){return fSlotPerCrate;};
  int    NPlanes(){return fNplanes;};
  int   NModulesSize(){return fNmodules.size();};
  int  NModules(int i){return fNmodules.at(i);};
  
  int GeoInfoSize(){return fGeoInfo.size();};
  TGeoInfo GeoInfo(int i){return fGeoInfo.at(i);};
  TDigInfo DigInfo(){return fDigInfo;};
  
  void SetDetName(std::string detname){fDetName = detname;};
  void SetDetType(det_type type){fDetType = type;};
  void SetNChan(int nchan){fNchan = nchan;};
  void SetChanPerSlot(int chanperslot){fChanPerSlot = chanperslot;};
  void SetSlotPerCrate(int slotpercrate){fSlotPerCrate = slotpercrate;};
  void SetNPlanes(int nplanes){fNplanes = nplanes;};
  void AddNModules(int nmodules){fNmodules.push_back(nmodules);};
  
  void AddGeoInfo(TGeoInfo geoinfo){fGeoInfo.push_back(geoinfo);};
  void SetDigInfo(TDigInfo diginfo){fDigInfo = diginfo;};
  
 private:
  std::string fDetName;      // Detector name
  det_type     fDetType;      // DetectorType
  int           fNchan;        // Total number of channels over all detector
  int            fChanPerSlot;  // Number of channels per slot
  int             fSlotPerCrate; // Number of slots per crate
  int              fNplanes;      // Number of planes // useful e.g. GEM, CDet
  std::vector<int>  fNmodules;     // Number of modules per plane // useful e.g. GEM, CDet
  
  std::vector<TGeoInfo> fGeoInfo;
  TDigInfo fDigInfo;
  
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
  TF1 *fPulseModel;
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
  UInt_t ADC(){return fADC;};

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
  UInt_t fADC;// One unique ADC value ?

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

