#ifndef _TSBSSIMCHER_H
#define _TSBSSIMCHER_H

#include "TSBSSimDetector.h"
#include "TSBSSimAuxi.h"

class TSBSSimCher : public TSBSSimDetector {
public:
  TSBSSimCher(const char* name, short id);
  virtual ~TSBSSimCher();
  // This loads the simulation event data
  virtual void LoadEventData(const std::vector<g4sbshitdata*> &evbuffer);
  virtual void Digitize(TSBSSimEvent &event);
  
  virtual void Clear();
  
  // Initialize
  void Init();
  
 private:
  
  TSPEModel *fSPE;
  std::vector<TPMTSignal> fSignals;
  
  ClassDef(TSBSSimCher,1)
};

#endif //_TSBSSIMCHER_H
