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
TSBSSimTrack::TSBSSimTrack( Int_t number, Int_t pid,
			    const TVector3& vertex, const TVector3& momentum )
  : Podd::MCTrack( number, pid, vertex, momentum )
{
}

//-----------------------------------------------------------------------------
TSBSSimTrack::TSBSSimTrack() : Podd::MCTrack()
{
}

//Special function for debugging
Double_t TSBSSimTrack::MCFitX_print() const
{
  printf("MCFit[0-3]: %f %f %f %f \n", fMCFitPar[0], fMCFitPar[1], fMCFitPar[2], fMCFitPar[3]);
  printf("MCFit[4-8]: %f %f %f %f %f \n", fMCFitPar[4], fMCFitPar[5], fMCFitPar[6], fMCFitPar[7], fMCFitPar[8]);
  // printf("RcFit[0-3]: %f %f %f %f \n", fRcFitPar[0], fRcFitPar[1], fRcFitPar[2], fRcFitPar[3]);
  // printf("RcFit[4-8]: %f %f %f %f %f \n", fRcFitPar[4], fRcFitPar[5], fRcFitPar[6], fRcFitPar[7], fRcFitPar[8]);
  return fMCFitPar[0];
}

//-----------------------------------------------------------------------------
TSBSSimEvent::TSBSSimEvent()
  : fRunID(0), fEvtID(0), fNSignal(0), fMCTracks(0)
{
}

//-----------------------------------------------------------------------------
TSBSSimEvent::TSBSSimEvent( UInt_t ntracks )
  : fRunID(0), fEvtID(0), fNSignal(0)
{
  if( ntracks == 0 ) ntracks = 1;
  fMCTracks = new TClonesArray( "TSBSSimTrack", ntracks );
}

//-----------------------------------------------------------------------------
TSBSSimEvent::~TSBSSimEvent()
{
  delete fMCTracks;
}

//-----------------------------------------------------------------------------
TSBSSimTrack* TSBSSimEvent::AddTrack( Int_t number, Int_t pid,
				      const TVector3& vertex,
				      const TVector3& momentum )
{
  // Add a physics track with the given parameters

  return
    new( (*fMCTracks)[GetNtracks()] ) TSBSSimTrack( number, pid,
						    vertex, momentum );
}

//-----------------------------------------------------------------------------
void TSBSSimEvent::Clear( const Option_t* opt )
{
  // Clear the event in preparation for reading next tree entry

  TString sopt(opt);

  fNSignal = 0;
  fPMTHits.clear();
  fDetectorData.clear();
  //fMCClusterHitID.clear();
}

//-----------------------------------------------------------------------------
Int_t TSBSSimEvent::GetNtracks() const
{
  // Get number of physics tracks

  return fMCTracks->GetLast()+1;
}

//-----------------------------------------------------------------------------
void TSBSSimEvent::Print( const Option_t* opt ) const
{
  // Print current event info

  cout << ">>>>> =====================================" << endl;
  cout << "Run number:               " << fRunID << endl;
  cout << "Event number:               " << fEvtID << endl;
  cout << "Event weight:               " << fWeight << endl;
  cout << "Number of hits:             " << fPMTHits.size()   << endl;
  
  TString sopt(opt);
  bool do_all   = sopt.Contains("all",   TString::kIgnoreCase);
  bool do_hit   = sopt.Contains("hit",   TString::kIgnoreCase) || do_all;
  bool do_clust = sopt.Contains("clust", TString::kIgnoreCase) || do_all;
  bool do_track = sopt.Contains("track", TString::kIgnoreCase) || do_all;
  
  if( do_track && fMCTracks ) {
    for( Int_t i=0; i<GetNtracks(); ++i ) {
      fMCTracks->UncheckedAt(i)->Print(opt);
    }
  }
  /*
  if( do_clust ) {
    cout << "Cluster list size = " << fMCClusterHitID.size() << endl;
    for(int i_ = 0; i_<fMCClusterHitID.size(); i_++){
      cout << " i_ = " << i_ << "; MC track ID " << fMCClusterHitID[i_].first 
	   << " cluster size " << fMCClusterHitID[i_].second.size() << endl;
      for(int j_ = 0; j_<fMCClusterHitID[i_].second.size(); j_++){
	cout << " j_ = " << " PMT hit ID " << fMCClusterHitID[i_].second[j_] << endl;
      }
    }
  }
  */
  if( do_hit ) {
    for( vector<PMTHit>::const_iterator ih = fPMTHits.begin();
	 ih != fPMTHits.end(); ++ih ) {
      const PMTHit& h = *ih;
      cout << "hit = " << h.fID 
	   << ", type = "    << h.fType 
	   << ", src = "    << h.fSource 
	   << endl;
      cout << "PMT row = "   << h.fPMTrow 
	   << ", PMT col = "   << h.fPMTcol 
	   << ", X_PMT = "   << h.fXPMT 
	   << ", Y_PMT = "   << h.fYPMT 
	   << ", number of photoelectrons = " << h.fNpe 
	   << ", time = " << h.fTime 
	   << " ns; TDC rising time = " << h.fTDCtime[0] 
	   << " ns, TDC falling time = " << h.fTDCtime[1] 
	   << endl;
      cout << " Channel = " << h.fChannel << endl
	   << "VETROC number " << h.fVETROCID
	   << ", TDC 'words' for rising time: " << h.fTDC[0] 
	   << "; for falling time: " << h.fTDC[1] 
	   << endl;
    }
  }
}

//-----------------------------------------------------------------------------
ClassImp(TSBSSimEvent)
