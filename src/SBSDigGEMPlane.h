#ifndef SBSDIGGEMPLANE_H
#define SBSDIGGEMPLANE_H

#include <iostream>
#include <vector>
#include <map>
#include <TROOT.h>
//#include "TH1D.h"
//#include "TRandom3.h"

class SBSDigGEMPlane {
 public:
  SBSDigGEMPlane();
  SBSDigGEMPlane(int nstrips, int nsamples = 6, double thr = 100);
  virtual ~SBSDigGEMPlane();
  void Clear();
  
  //void SetStripThreshold(double thr){StripThr = ;};
  
  UShort_t GetADC(int strip, int samp){return fStripADC[strip*fNSamples+samp];};
  UInt_t GetADCSum(int strip){return fStripADCsum[strip];};
  
 private:
  // ADC sampled value of strip array of each axis
  Int_t fNStrips;
  Int_t fNSamples;
  Double_t fStripThr;//threshold for ADC sum
  UInt_t* fStripADCsum;
  UShort_t* fStripADC;
  
  //ClassDef(SBSDigGEMPlane, 1)
};
#endif
