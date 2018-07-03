#ifndef _TSBSSIMDETECTOR_H
#define _TSBSSIMDETECTOR_H

#include <vector>
#include "g4sbs_types.h"
#include "THaAnalysisObject.h"

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
};

#endif // _TSBSSIMDETECTOR_H
