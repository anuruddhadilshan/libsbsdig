#include "SBSDigPMTDet.h"
#include "TMath.h"
#define DEBUG 0

using namespace std;

SBSDigPMTDet::SBSDigPMTDet()
{
}

SBSDigPMTDet::SBSDigPMTDet(UInt_t nchan):
  fNChan(nchan)
{
  for(int i = 0; i<fNChan; i++)PMTmap[i] = PMTSignal();
}

SBSDigPMTDet::SBSDigPMTDet(UInt_t nchan, double NpeChargeConv, double sigmapulse, double gatewidth):
  fNChan(nchan)
{
  for(int i = 0; i<fNChan; i++)PMTmap[i] = PMTSignal(NpeChargeConv);
  fRefPulse = new SPEModel(sigmapulse, 0, -gatewidth/2., gatewidth/2.);
}

SBSDigPMTDet::~SBSDigPMTDet()
{
  
}

void SBSDigPMTDet::SetSamples()
{
  for(int i = 0; i<fNChan; i++)PMTmap[i].SetSamples(-fGateWidth/2, fGateWidth/2, 4.0);
}

void SBSDigPMTDet::Clear(bool dosamples)
{
  for(int i = 0; i<fNChan; i++)PMTmap[i].Clear(dosamples);
}
