#ifndef SBSDIGPMTDET_H
#define SBSDIGPMTDET_H

#include <iostream>
#include <vector>
#include <map>
#include "SBSDigPMTSignal.h"

//________________________________
class SBSDigPMTDet {
 public:
  SBSDigPMTDet();
  SBSDigPMTDet(int nchan, double npechargeconv);
  virtual ~SBSDigPMTDet();
  void Clear();
  
 private:
  int fNChan;
  std::map<int, PMTSignal> PMTmap;
  
};

#endif
