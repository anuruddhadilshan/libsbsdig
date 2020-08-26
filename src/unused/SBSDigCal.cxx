#include "SBSDigCal.h"
#include "TMath.h"
#define DEBUG 0

using namespace std;

SBSDigCal::SBSDigCal()
{
}

SBSDigCal::SBSDigCal(int nchan, double npechargeconv)
{
  fNChan = nchan;
  for(int i = 0; i<fNChan; i++)PMTmap[i] = PMTSignal(npechargeconv);
}

SBSDigCal::~SBSDigCal()
{
  
}

void SBSDigCal::Clear()
{
  for(int i = 0; i<fNChan; i++)PMTmap[i].Clear();
}
