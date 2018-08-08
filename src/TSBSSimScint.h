#ifndef _TSBSSIMSCINT_H
#define _TSBSSIMSCINT_H

#include "TSBSSimDetector.h"
#include "TSBSSimAuxi.h"

class TF1;
class TF1Convolution;
class TTree;
class TFile;

class TSBSSimScint : public TSBSSimDetector {
public:
  TSBSSimScint(const char* name);
  virtual ~TSBSSimScint();
  // This loads the simulation event data
  virtual void LoadEventData(const std::vector<g4sbshitdata*> &evbuffer);
  virtual void Digitize(TSBSSimEvent &event);
  
  virtual void Clear();
  
  // Initialize
  void Init();
  
 private:
  
  TNPEModel *fNPE;
  std::vector<TPMTSignal> fSignals;
  
  ClassDef(TSBSSimScint,1)
};

#endif //_TSBSSIMHCAL_H
