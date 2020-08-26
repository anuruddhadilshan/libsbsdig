#ifndef SBSDIGGEMDET_H
#define SBSDIGGEMDET_H

#include <iostream>
#include <vector>
#include <map>
#include "SBSDigGEMPlane.h"

//________________________________
class SBSDigGEMDet {
 public:
  SBSDigGEMDet();
  SBSDigGEMDet(int nplanes, double* nstrips, int nsamp, double zsup_thr);
  virtual ~SBSDigGEMDet();
  void Clear();
  
 private:
  int fNPlanes;
  std::map<int, SBSDigGEMPlane> GEMPlanes;
  
};

#endif
