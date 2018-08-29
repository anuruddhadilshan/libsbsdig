#ifndef _TSBSSIMSCINT_H
#define _TSBSSIMSCINT_H

#include "TSBSSimDetector.h"
#include "TSBSSimAuxi.h"

class TSBSSimScint : public TSBSSimDetector {
public:
  TSBSSimScint(const char* name, short id);
  virtual ~TSBSSimScint();
  // This loads the simulation event data
  virtual void LoadEventData(const std::vector<g4sbshitdata*> &evbuffer);
  virtual void Digitize(TSBSSimEvent &event);
  
  virtual void Clear();
  
  // Initialize
  void Init();
  
 private:
  
  TSPEModel *fSPE;
  std::vector<TPMTSignal> fSignals;
  
  ClassDef(TSBSSimScint,1)
};

#endif //_TSBSSIMSCINT_H
