#include "SBSDigGEMPlane.h"
#include "TMath.h"

using namespace std;

//
// Class SBSDigGEMPlane
//
SBSDigGEMPlane::SBSDigGEMPlane() :
  fNStrips(3840), fNSamples(6), fStripThr(100), fXoffset(0), fROangle(0)
{
}

SBSDigGEMPlane::SBSDigGEMPlane(int nstrips, int nsamples, double thr, double offset, double roangle) :
  fNStrips(nstrips), fNSamples(nsamples), fXoffset(offset), fStripThr(thr), fROangle(roangle)
{
  fdX = fNStrips*4.e-4;
  
  fStripADCsum = new UInt_t[fNStrips];
  fStripADC = new UShort_t[fNStrips*fNSamples];
  Clear();
}

SBSDigGEMPlane::~SBSDigGEMPlane()
{
  Clear();
}


void SBSDigGEMPlane::Clear()
{
  memset(fStripADCsum, 0, fNStrips*sizeof(UInt_t));
  memset(fStripADC, 0, fNStrips*fNSamples*sizeof(UShort_t));
}

//ClassImp(SBSDigGEMPlane);

