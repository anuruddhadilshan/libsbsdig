#ifndef SBSDIGCHER_H
#define SBSDIGCHER_H

#include <iostream>
#include <vector>
#include <map>
#include "SBSDigPMTSignal.h"

//________________________________
class SBSDigCher {
 public:
  SBSDigCher();
  SBSDigCher(int nchan, double npechargeconv);
  virtual ~SBSDigCher();
  void Clear();
  
  //private:
  int fNChan;
  std::map<int, PMTSignal> PMTmap;
  
};

#endif
