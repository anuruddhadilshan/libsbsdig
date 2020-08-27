#include "SBSDigGEMDet.h"
#include "TMath.h"
#define DEBUG 0

using namespace std;

SBSDigGEMDet::SBSDigGEMDet()
{
}

SBSDigGEMDet::SBSDigGEMDet(UInt_t nplanes, double* nstrips, int nsamp, double zsup_thr):
  fNPlanes(nplanes)
{
  for(uint i = 0; i<fNPlanes; i++)GEMPlanes[i] = SBSDigGEMPlane(nstrips[i], nsamp, zsup_thr);
}

SBSDigGEMDet::~SBSDigGEMDet()
{
  
}

void SBSDigGEMDet::Clear()
{
  for(uint i = 0; i<fNPlanes; i++)GEMPlanes[i].Clear();
}
