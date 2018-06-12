#include <iostream>

#include "TSBSDetData.h"

using namespace std;

TSBSDetData::TSBSDetData ()
{
  fRunID = 0; 
  fEvtID = 0; 
  fSource = 0;
  //fHitData->Expand(h);
}

TSBSDetData::~TSBSDetData()
{
  //fHitData->Delete();
}

void
TSBSDetData::ClearEvent() 
{
  //fHitData->Clear();
};

void
TSBSDetData::InitEvent (UInt_t h)
{
  /*
  if (h <= 0)
    return;

  fHitData->Expand(h);
  */
}

void 
TSBSDetData::Print() const
{
  cout << "Run " << GetRun() << " Event " << GetEvent() 
    //<< " " << GetNHit() << " hits" 
       << endl;
}

void 
TSBSDetData::PrintHit (UInt_t k) const
{
}
