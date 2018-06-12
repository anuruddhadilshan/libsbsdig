#include <iostream>

#include "TSBSCherData.h"

using namespace std;

/*
TSBSCherData::TSBSCherData (UInt_t h)
{
  fHitData = new TClonesArray("CherHitData", h);
}

TSBSCherData::~TSBSCherData()
{
  fHitData->Delete();
}

void
TSBSCherData::ClearEvent ()
{
  fHitData->Clear();
}

void
TSBSCherData::InitEvent (UInt_t h)
{
  if (h <= 0)
    return;

  fHitData = new TClonesArray("CherHitData", h);
}
*/

TSBSCherData::TSBSCherData (UInt_t h)
{
  fHitData.reserve(h);
}

TSBSCherData::~TSBSCherData()
{
}

void
TSBSCherData::ClearEvent() 
{
  fHitData.clear();
};

void
TSBSCherData::InitEvent (UInt_t h)
{
  if (h <= 0)
    return;

  fHitData.resize(h);
}

void 
TSBSCherData::PrintHit (UInt_t k) const
{
  cout << "  Event " << GetEvent() << ":" << endl;
  cout << "    Hit time: " << GetHitTime(k) << " ns, RMS " 
       <<  GetHitTimeRMS(k) << " ns " << endl;
  cout << "    Hit PE yield: " << GetHitPEyield(k) << endl;
  cout << "    Hit PMT: " << GetHitPMTID(k) << endl;
  cout << "    PMT X_transport: " << GetHitXPMT(k) << endl;
  cout << "    PMT Y_transport: " << GetHitYPMT(k) << endl;
  cout << "    Particle type: " << GetParticleType(k) << endl;
  cout << "    Photon detection position: " << GetPositionDet(k).X()
       << " " << GetPositionDet(k).Y() 
       << " " << GetPositionDet(k).Z() 
       << " (mag. " << GetPositionDet(k).Mag() << ")"
       << " eV" << endl;
  if(GetOriginVolFlag(k)==2){
    cout << "    MCtrack PID: " << GetMCtrackPID(k) << endl;
    cout << "    MC track momentum: " << GetMCtrackMomentum(k).X()
	 << " " << GetMCtrackMomentum(k).Y() 
	 << " " << GetMCtrackMomentum(k).Z() 
	 << " (mag. " << GetMCtrackMomentum(k).Mag() << ")"
	 << " MeV" << endl;
    cout << "    MC track vertex: " << GetMCtrackVertex(k).X()
	 << " " << GetMCtrackVertex(k).Y() 
	 << " " << GetMCtrackVertex(k).Z() 
	 << " m" << endl;
  }
   cout << "    Photon volume origin flag (1: aerogel, 2: gas, 3: vessel, 4: PMT window):  " << GetOriginVolFlag(k) << endl;
}
