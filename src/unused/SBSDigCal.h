#ifndef SBSDIGCAL_H
#define SBSDIGCAL_H

#include <iostream>
#include <vector>
#include <map>
#include "SBSDigPMTSignal.h"

//________________________________
class SBSDigCal {
 public:
  SBSDigCal();
  SBSDigCal(int nchan, double npechargeconv);
  virtual ~SBSDigCal();
  void Clear();
  
 private:
  int fNChan;
  std::map<int, PMTSignal> PMTmap;
  
};

#endif
