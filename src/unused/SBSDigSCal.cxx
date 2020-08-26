#include "SBSDigSCal.h"
#include "TMath.h"
#define DEBUG 0

using namespace std;

SBSDigSCal::SBSDigSCal()
{
}

SBSDigCal::SBSDigSCal(int nchan, double npechargeconv)
{
  fNChan = nchan;
  for(int i = 0; i<fNChan; i++)PMTmap[i] = PMTSignal(npechargeconv);
}

SBSDigCal::~SBSDigSCal()
{
  
}

void SBSDigSCal::Clear()
{
  for(int i = 0; i<fNChan; i++)PMTmap[i].Clear();
}
