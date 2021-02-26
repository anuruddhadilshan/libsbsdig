#include "SBSDigCher.h"
#include "TMath.h"
#define DEBUG 0

using namespace std;

SBSDigCher::SBSDigCher()
{
}

SBSDigCher::SBSDigCher(int nchan, double npechargeconv)
{
  fNChan = nchan;
  for(int i = 0; i<fNChan; i++)PMTmap[i] = PMTSignal(npechargeconv);
}

SBSDigCher::~SBSDigCher()
{
  
}

void SBSDigCher::Clear()
{
  for(int i = 0; i<fNChan; i++)PMTmap[i].Clear();
}

