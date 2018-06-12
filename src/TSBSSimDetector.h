#ifndef _TSBSSIMDETECTOR_H
#define _TSBSSIMDETECTOR_H

#include <vector>
#include "g4sbs_types.h"

class g4sbshitdata;
class TSBSSimEvent;
class TSBSDBManager;

class TSBSSimDetector {
public:
  TSBSSimDetector();
  virtual ~TSBSSimDetector();
  virtual void LoadEventData(const std::vector<g4sbshitdata*> &evbuffer) = 0;
  virtual void Digitize(TSBSSimEvent &event) = 0;
  virtual bool HasData() { return fHasData; }
protected:
  void SetHasDataFlag(bool has_data) { fHasData = has_data; }
private:
  bool fHasData;
  det_type   fDetType;//detector type (see g4sbs_types.h) no need right now, but it might come in handy
  TSBSDBManager* fDBmanager;
};

#endif // _TSBSSIMDETECTOR_H
