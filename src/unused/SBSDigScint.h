#ifndef SBSDIGSCINT_H
#define SBSDIGSCINT_H

#include <iostream>
#include <vector>
#include <map>
#include "SBSDigPMTSignal.h"

//________________________________
class SBSDigScint {
 public:
  SBSDigScint();
  SBSDigScint(int nchan, double npechargeconv);
  virtual ~SBSDigScint();
  void Clear();
  
  //private:
  int fNChan;
  std::map<int, PMTSignal> PMTmap;
  
};

#endif
