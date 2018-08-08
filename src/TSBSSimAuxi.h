#ifndef _TSBSSIMAUXI_H
#define _TSBSSIMAUXI_H

#include <vector>
#include "g4sbs_types.h"
#include "THaAnalysisObject.h"
#include "g4sbs_types.h"
#include "TF1.h"
#include "TF1Convolution.h"

class TNPEModel : public TObject {
 public:
  TNPEModel(DigInfo diginfo, const char* detname, int npe = 1);
  double GetNpe(){return fNpe;};
  void   SetNpe(double npe){fNpe = npe;};
  double GetStartTime(){return fStartTime;};
  void   SetStartTime(double t){fStartTime = t;};
  double GetADCconversion(){return fDigInfo.fADCconversion;};
  double GetTDCconversion(){return fDigInfo.fTDCconversion;};

  double Eval(int chan, double t);
  void   FindLeadTrailTime(int chan, double &t_lead, double &t_trail);
  bool   PulseOverThr(int chan = 0);
    
  ~TNPEModel(){};
  
 private:
  DigInfo fDigInfo;
  TF1 *fModel;
  double fScale;
  double fNpe;
  double fStartTime;
  
  /*
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
  */
  ClassDef(TNPEModel,1);
};

class TPMTSignal : public TObject {
 public:
  /*
    std::vector<double> samples;
    std::vector<double> samples_raw;
    double sumedep;
    double mint;
    double maxt;
    int nbins;
    int nbins_raw;
    int npe;
    int sum;
    int dnraw;
    double dx_samples;
    double dx_raw;
  */
  double fSumedep;//Not forced to use it for everything
  //double fADC;
  int fNpe;
  int fADC;// One unique ADC value ?
  //TDCs: multiple values possible.
  std::vector<double> fLeadTimes;
  std::vector<double> fTrailTimes;
  
  TPMTSignal();
  void Fill(int chan, TNPEModel *model, double t, double toffset = 0.0);
  void Digitize();
  void Clear();
  ~TPMTSignal();
  
  ClassDef(TPMTSignal,1);
};

//
// Classes for DB information: TO BE COMPLETED LATER! 
// structs already exist for DB and code is functional. 
//
class TSignalInfo : public TObject {
 public:
  TSignalInfo() {};
 TSignalInfo(int apid, int atid):fPID(apid), fTID(atid) {};
  ~TSignalInfo() {};
  
  double GetPID(){return fPID;};
  double GetTID(){return fTID;};

  void SetPID(int apid){fPID = apid;};
  void SetTID(int atid){fTID = atid;};
  
 private:
  int fPID;
  int fTID;
  
  ClassDef(TSignalInfo, 1);
};

class TSpectroInfo : public TObject{
 public:
  TSpectroInfo() {};
  TSpectroInfo(int ndets):fNDets(ndets) {};
  ~TSpectroInfo() {};

  //TSignalInfo
  
  
 private:
  
  double fMCangle;
  int fNDets;
  std::vector<TSignalInfo> fMCsignalInfo;
  std::vector<std::string> fDetNames;
  
  ClassDef(TSpectroInfo, 1);
};

/*
class  : public TObject{
 public:

  ClassDef(, 1);
};

class  : public TObject{
 public:

  ClassDef(, 1);
};

class  : public TObject{
 public:

  ClassDef(, 1);
};
*/
/*
struct TSignalInfo{
  int pid;
  int tid;
  
  // Int_t fillBitsGEM;
  // Int_t fillBitsEC;
  // Int_t signalSector; //used if map sector
  // Double_t ECEDep;
  // Double_t momentum;
  // Double_t R;
  
  SignalInfo() {}
  SignalInfo(int apid, int atid):pid(apid), tid(atid) {}
  ~SignalInfo() {}
};

struct SpectroInfo{
  double fMCangle;
  std::vector<SignalInfo> MCsignalInfo;
  int fNDets;
  std::vector<std::string> fDetNames;
  
  SpectroInfo() {}
  SpectroInfo(int ndets):fNDets(ndets) {}
  ~SpectroInfo() {}
};

struct GeoInfo{
  // I don't think this can fit all detectors. 
  // Might have to do a class or something so it can be inherited ?
  int     fNrows;      // number of rows
  int     fNcols;      // number of columns
  double  fXsize;      // detector X size (in transport coordinates)
  double  fYsize;      // detector Y size (in transport coordinates)
  double  fZpos;       // detector position on spectrometer axis
};

struct DigInfo{
  // same remark as for GeoInfo
  double  fROImpedance;          // readout impedance
  std::vector<double>  fGain;     // Gain 
  std::vector<double>  fPedestal; // Pedestal value (adc value)
  std::vector<double>  fPedNoise; // Pedestal noise (adc value)
  double  fTriggerJitter;         // trigger jitter (ns)
  double  fTriggerOffset;        // trigger offset (ns)
  double  fGateWidth;           // gate width (ns)
  double  fSPEtau;            // tau param for SPE
  double  fSPEsig;           // sigma param for SPE
  double  fSPEtransittime; // pmt transit time param for SPE
};

struct DetInfo{
  std::string fDetName;      // Detector name
  det_type     fDetType;      // DetectorType
  int           fNChan;        // Total number of channels over all detector
  int            fChanPerSlot;  // Number of channels per slot
  int             fSlotPerCrate; // Number of slots per crate
  int              fNPlanes;      // Number of planes // useful e.g. GEM, CDet
  std::vector<int>  fNModules;     // Number of modules per plane // useful e.g. GEM, CDet
  
  std::vector<GeoInfo> fGeoInfo;
  DigInfo fDigInfo;
  DetInfo()
  {
    fNModules.clear();
    fGeoInfo.clear();
  }
  DetInfo(const std::string detname)
  {
    fDetName = detname;
    fNModules.clear();
    fGeoInfo.clear();
  }
  ~DetInfo(){}

};
*/

#endif // _TSBSSIMAUXI_H

