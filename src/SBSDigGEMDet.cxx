#include "SBSDigGEMDet.h"
#include "TMath.h"
#define DEBUG 0

using namespace std;

SBSDigGEMDet::SBSDigGEMDet()
{
}

SBSDigGEMDet::SBSDigGEMDet(UShort_t uniqueid, UInt_t nplanes, double* nstrips, int nsamp, double zsup_thr):
  fUniqueID(uniqueid), fNPlanes(nplanes)
{
  //for(uint i = 0; i<fNPlanes; i++)GEMPlanes[i] = SBSDigGEMPlane(nstrips[i], nsamp, zsup_thr);
  for(uint i = 0; i<fNPlanes; i++)GEMPlanes.push_back(SBSDigGEMPlane(nstrips[i], nsamp, zsup_thr));
}

SBSDigGEMDet::~SBSDigGEMDet()
{
  
}

void SBSDigGEMDet::Clear()
{
  for(uint i = 0; i<fNPlanes; i++)GEMPlanes[i].Clear();
}
