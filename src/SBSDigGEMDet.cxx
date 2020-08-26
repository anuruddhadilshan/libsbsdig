#include "SBSDigGEMDet.h"
#include "TMath.h"
#define DEBUG 0

using namespace std;

SBSDigGEMDet::SBSDigGEMDet()
{
}

SBSDigGEMDet::SBSDigGEMDet(int nplanes, double* nstrips, int nsamp, double zsup_thr)
{
  fNPlanes = nplanes;
  for(int i = 0; i<fNPlanes; i++)GEMPlanes[i] = SBSDigGEMPlane(nstrips[i], nsamp, zsup_thr);
}

SBSDigGEMDet::~SBSDigGEMDet()
{
  
}

void SBSDigGEMDet::Clear()
{
  for(int i = 0; i<fNPlanes; i++)GEMPlanes[i].Clear();
}
