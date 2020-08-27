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
  SBSDigGEMDet(UInt_t nplanes, double* nstrips, int nsamp, double zsup_thr);
  virtual ~SBSDigGEMDet();
  void Clear();
  
 private:
  UInt_t fNPlanes;
  std::map<int, SBSDigGEMPlane> GEMPlanes;
};

#endif
