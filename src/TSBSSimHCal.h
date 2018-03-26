#ifndef _TSBSSIMHCAL_H
#define _TSBSSIMHCAL_H

#include "TSBSSimDetector.h"

class TF1;
class TF1Convolution;
class TTree;
class TFile;

class TSBSSimHCal : public TSBSSimDetector {
public:
  TSBSSimHCal();
  virtual ~TSBSSimHCal();
  // This loads the simulation event data
  virtual void LoadEventData(const std::vector<g4sbshitdata*> &evbuffer);
  virtual void Digitize(TSBSSimEvent &event);

  virtual void Clear();

  // Initialize
  void Init();
  struct SPEModel {
    double gain_pmt;
    double resistance; //ohm
    double qe; //
    double unit;
    double scale;
    TF1 *model;
    SPEModel();
    double Eval(double t);
    TF1 *fFunc1;
    TF1 *fFunc2;
    TF1Convolution *fConvolution;
    double mint;
    double start_t;
    double maxt;
    double tao;
    double sig;
    double t0;
  };
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
    void Fill(SPEModel *model,double t, double toffset = 0.0);
    void Digitize();
    void Clear();
  };
private:
  SPEModel *fSPE;
  std::vector<Signal> fSignals;
};

#endif //_TSBSSIMHCAL_H
