#ifndef _TSBSSIMHCAL_H
#define _TSBSSIMHCAL_H

#include "TSBSSimDetector.h"
#include "TSBSSimAuxi.h"

// class TF1;
// class TF1Convolution;
// class TTree;
// class TFile;

class TSBSSimHCal : public TSBSSimDetector {
public:
  TSBSSimHCal(const char* name, short id);
  virtual ~TSBSSimHCal();
  // This loads the simulation event data
  virtual void LoadEventData(const std::vector<g4sbshitdata*> &evbuffer);
  //The following function is to accumulate data without reinitializing the event.
  virtual void LoadAccumulateData(const std::vector<g4sbshitdata*> &evbuffer);
  virtual void Digitize(TSBSSimEvent &event);

  virtual void Clear(Option_t *op = "");

  // Silence compiler warnings about Init from parent class
  using THaAnalysisObject::Init;
  // Initialize
  void Init();

  struct Signal {
    SimEncoder::fadc_data fadc;// encoded FADC samples
    SimEncoder::tdc_data tdc;// encoded TDC values
    //std::vector<double> samples;
    std::vector<double> samples_raw;// "unencoded" FADC samples 
    std::vector<double> times_histo;// "unencoded" TDC values 
    int nbins_times;// number of TDC bins!
    double sumedep;// Energy deposit
    double mint;// DAQ window min 
    double maxt;// DAQ window max 
    int nbins;// number of ADC bins
    int nbins_raw;// number of ADC subdivisions
    int npe;// Number of photoelectrons
    double tdc_time;//tdc time to feed SimEncoder::tdc_data tdc
    bool met_tdc_thresh;// true if signal is large enough to trigger TDC
    int sum;// ADC sum 
    int dnraw;// number of samples read off the FADC
    double dx_samples;// size of FADC sample 
    double dx_raw;// subdivision of FADC samples
    double dx_raw_time;// TDC bin size
    
    short mc_source;// MC source: sig ? /bkgd ?
    int trid;// original particle track ID 
    int pid;// original particle PID
    Signal();
    void FillNPE(TSPEModel *model, double pulsenorm, double t, double toffset = 0.0);
    void Fill(double t);
    void Digitize(TSPEModel *model, double pulsenorm, double toffset,
        double max_val);
    void DigitizeOld();
    void Clear();
  };
private:
  //SPEModel *fSPE;
  TSPEModel *fSPE;
  std::vector<Signal> fSignals;
  bool fHasFADC;
  // TODO: Try to use the standard TPMTSignal class (but must have
  // the ability to provide samples)
  //std::vector<TPMTSignal> fSignals;

  ClassDef(TSBSSimHCal,1)
};

#endif //_TSBSSIMHCAL_H
