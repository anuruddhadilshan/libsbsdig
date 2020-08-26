#include "SBSDigPMTDet.h"
#include "TMath.h"
#define DEBUG 0

using namespace std;

SBSDigPMTDet::SBSDigPMTDet()
{
}

SBSDigPMTDet::SBSDigPMTDet(int nchan, double npechargeconv)
{
  fNChan = nchan;
  for(int i = 0; i<fNChan; i++)PMTmap[i] = PMTSignal(npechargeconv);
}

SBSDigPMTDet::~SBSDigPMTDet()
{
  
}

void SBSDigPMTDet::Clear()
{
  for(int i = 0; i<fNChan; i++)PMTmap[i].Clear();
}
