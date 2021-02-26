#ifndef SBSDIGSCAL_H
#define SBSDIGSCAL_H

#include <iostream>
#include <vector>
#include <map>
#include "SBSDigPMTSignal.h"

//________________________________
class SBSDigSCal {
 public:
  SBSDigSCal();
  SBSDigSCal(int nchan, double npechargeconv);
  virtual ~SBSDigSCal();
  void Clear();
  
 private:
  int fNChan;
  std::map<int, PMTSignal> PMTmap;
  
};

#endif
