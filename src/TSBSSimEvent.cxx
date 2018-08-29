//*-- Author :    Ole Hansen (ole@jlab.org)    9-Dec-2011

/////////////////////////////////////////////////////////////////////
//
//   TSBSSimEvent
//
//   Common class definitions for SoLID simulation decoder
//
/////////////////////////////////////////////////////////////////////

#include "TSBSSimEvent.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TMath.h"

#include <iostream>

using namespace std;

//-----------------------------------------------------------------------------
TSBSSimEvent::TSBSSimEvent()
  : fRunID(0), fEvtID(0), fNSignal(0)//, fMCTracks(0)
{
}

//-----------------------------------------------------------------------------
void TSBSSimEvent::Clear( const Option_t* opt )
{
  // Clear the event in preparation for reading next tree entry

  TString sopt(opt);

  fNSignal = 0;
  fDetectorData.clear();
  fSimDetectorData.clear();
}

//-----------------------------------------------------------------------------
void TSBSSimEvent::Print( const Option_t* opt ) const
{
  // Print current event info

  cout << ">>>>> =====================================" << endl;
  cout << "Run number:               " << fRunID << endl;
  cout << "Event number:               " << fEvtID << endl;
  cout << "Event weight:               " << fWeight << endl;
  
  // TODO: rewrite this Print function
  // TString sopt(opt);
  // bool do_all   = sopt.Contains("all",   TString::kIgnoreCase);
  // bool do_hit   = sopt.Contains("hit",   TString::kIgnoreCase) || do_all;
  // bool do_clust = sopt.Contains("clust", TString::kIgnoreCase) || do_all;
  // bool do_track = sopt.Contains("track", TString::kIgnoreCase) || do_all;
}

//-----------------------------------------------------------------------------
ClassImp(TSBSSimEvent)
