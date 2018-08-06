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

class TSPEModel {
 public:
  TSPEModel(DigInfo diginfo, const char* detname);
  double Eval(double t, int chan = 0);
  double GetStartTime(){return start_t;};
  void SetStartTime(double t){start_t = t;};
  
 private:
  DigInfo fDigInfo;
  TF1 *model;
  double fScale;
  double start_t;
  
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
  ClassDef(TSPEModel,1)
};

#endif // _TSBSSIMDETECTOR_H
