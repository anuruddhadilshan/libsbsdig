#include "SBSDigGEMPlane.h"
#include "TMath.h"

using namespace std;

//
// Class SBSDigGEMPlane
//
SBSDigGEMPlane::SBSDigGEMPlane() :
  fNStrips(3840), fNSamples(6), fStripThr(100)
{
}

SBSDigGEMPlane::SBSDigGEMPlane(int nstrips, int nsamples, double thr) :
  fNStrips(nstrips), fNSamples(nsamples), fStripThr(thr)
{
  fStripADCsum = new UInt_t[fNStrips];
  fStripADC = new UShort_t[fNStrips*fNSamples];
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

