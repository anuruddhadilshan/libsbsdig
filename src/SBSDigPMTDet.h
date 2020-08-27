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
  SBSDigPMTDet(UInt_t nchan);
  SBSDigPMTDet(UInt_t nchan, double NpeChargeConv, double sigmapulse, double gatewidth);
  virtual ~SBSDigPMTDet();
  void Clear(bool dosamples = false);
  void SetSamples();
  
  //private:
  UInt_t fNChan;
  std::map<int, PMTSignal> PMTmap;
  SPEModel* fRefPulse;
  
  double fGain;
  double fPedestal;
  double fPedSigma;
  double fTrigOffset;
  double fThreshold;
  double fGateWidth;
  double fADCconv;
  double fADCbits;
  double fTDCconv;
  double fTDCbits;
};

#endif
