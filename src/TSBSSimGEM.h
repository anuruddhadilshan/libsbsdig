#ifndef _TSBSSIMGEM_H
#define _TSBSSIMGEM_H

#include "TSBSSimDetector.h"
#include "TSBSSimAuxi.h"

class TGEMSBSSimDigitization;
class TGEMSBSGEMChamber;
class TGEMSBSSpec;
class TGEMSBSDBManager;

class TSBSSimGEM : public TSBSSimDetector {
public:
  TSBSSimGEM(const char* name, short id);
  virtual ~TSBSSimGEM();
  // This loads the simulation event data
  virtual void LoadEventData(const std::vector<g4sbshitdata*> &evbuffer);
  //The following function is to accumulate data without reinitializing the event.
  virtual void LoadAccumulateData(const std::vector<g4sbshitdata*> &evbuffer);
  virtual void Digitize(TSBSSimEvent &event);

  virtual void Clear(Option_t *opt = "");
  virtual void EventStart();

  // Silence compiler warnings about Init from parent class
  using THaAnalysisObject::Init;
  // Initialize
  void Init();

  ClassDef(TSBSSimGEM,1)
private:
    TGEMSBSSimDigitization *fGEMDigi;
    TGEMSBSSpec *fSpec; ///< A bit much to create our own spectrometer, but oh well...
    TGEMSBSDBManager *fManager;
};

#endif //_TSBSSIMGEM_H
