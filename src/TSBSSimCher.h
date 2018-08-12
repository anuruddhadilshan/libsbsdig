#ifndef _TSBSSIMCHER_H
#define _TSBSSIMCHER_H

#include "TSBSSimDetector.h"
#include "TSBSSimAuxi.h"

class TF1;
class TF1Convolution;
class TTree;
class TFile;

class TSBSSimCher : public TSBSSimDetector {
public:
  TSBSSimCher();
  virtual ~TSBSSimCher();
  // This loads the simulation event data
  virtual void LoadEventData(const std::vector<g4sbshitdata*> &evbuffer);
  virtual void Digitize(TSBSSimEvent &event);

  virtual void Clear();

  // Initialize
  void Init();
  struct Signal {
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
    Signal();
    void Fill(TSPEModel *model,double t, double toffset = 0.0);
    void Digitize();
    void Clear();
  };
private:
  TSPEModel *fSPE;
  std::vector<Signal> fSignals;
  
  ClassDef(TSBSSimCher,1)
};

#endif //_TSBSSIMHCAL_H
