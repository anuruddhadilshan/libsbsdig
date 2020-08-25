#include "SBSDigScint.h"
#include "TMath.h"
#define DEBUG 0

using namespace std;

SBSDigScint::SBSDigScint()
{
}

SBSDigScint::SBSDigScint(int nchan, double npechargeconv)
{
  fNChan = nchan;
  for(int i = 0; i<fNChan; i++)PMTmap[i] = PMTSignal(npechargeconv);
}

SBSDigScint::~SBSDigScint()
{
  
}
