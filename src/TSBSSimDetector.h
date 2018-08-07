#ifndef _TSBSSIMDETECTOR_H
#define _TSBSSIMDETECTOR_H

#include <vector>
#include "g4sbs_types.h"
#include "THaAnalysisObject.h"
#include "g4sbs_types.h"
#include "TF1.h"
#include "TF1Convolution.h"

class g4sbshitdata;
class TSBSSimEvent;
class TSBSDBManager;

class TSBSSimDetector : public THaAnalysisObject {
public:
  TSBSSimDetector();
  virtual ~TSBSSimDetector();
  virtual void LoadEventData(const std::vector<g4sbshitdata*> &evbuffer) = 0;
  virtual void Digitize(TSBSSimEvent &event) = 0;
  virtual bool HasData() { return fHasData; }
protected:
  void SetHasDataFlag(bool has_data) { fHasData = has_data; }

  TSBSDBManager* fDBmanager;
  DetInfo fDetInfo;
private:
  bool fHasData;

  ClassDef(TSBSSimDetector,1)
};

/*
class TSPEModel : public TObject {
 public:
  TSPEModel(DigInfo diginfo, const char* detname);
  double Eval(double t, int chan = 0);
  double GetStartTime(){return fStartTime;};
  void SetStartTime(double t){fStartTime = t;};
  
 private:
  DigInfo fDigInfo;
  TF1 *model;
  double fScale;
  double fStartTime;
  
  // double gain_pmt;
  // double resistance; //ohm
  // double qe; //
  // double unit;
  // double scale;
  // TF1 *fFunc1;
  // TF1 *fFunc2;
  // TF1Convolution *fConvolution;
  // double mint;
  // double start_t;
  // double maxt;
  // double tau;
  // double sig;
  // double t0;
  
  ClassDef(TSPEModel,1)
};

class TPMTSignal : public TObject {
 public:
  // std::vector<double> samples;
  // std::vector<double> samples_raw;
  // double sumedep;
  // double mint;
  // double maxt;
  // int nbins;
  // int nbins_raw;
  // int npe;
  // int sum;
  // int dnraw;
  // double dx_samples;
  // double dx_raw;
  
  //double sumedep;
  int fNpe;
  int fADC;// One unique ADC value ?
  //TDCs: multiple values possible.
  std::vector<double> leadtimes;
  std::vector<double> trailtimes;
  
  

  TPMTSignal();
  void Fill(TSPEModel *model, double t, double toffset = 0.0);
  void Digitize();
  void Clear();
  
  ClassDef(TPMTSignal,1)
};
*/

#endif // _TSBSSIMDETECTOR_H
