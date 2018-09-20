//*-- Author :    Ole Hansen (ole@jlab.org)    9-Dec-2011

/////////////////////////////////////////////////////////////////////
//
//   TGEMSBSSimEvent
//
//   Common class definitions for SoLID simulation decoder
//
/////////////////////////////////////////////////////////////////////

#include "TGEMSBSSimEvent.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TMath.h"

#include <iostream>

using namespace std;

//-----------------------------------------------------------------------------
TGEMSBSSimTrack::TGEMSBSSimTrack( Int_t number, Int_t pid,
			    const TVector3& vertex, const TVector3& momentum, const TVector3& vertexAtTarget, const TVector3& momentumAtTarget  )
  : Podd::MCTrack( number, pid, vertex, momentum )
{
  vertex_target = vertexAtTarget;
  momentum_target= momentumAtTarget;
}

//-----------------------------------------------------------------------------
TGEMSBSSimTrack::TGEMSBSSimTrack() : Podd::MCTrack()
{
}

//Special function for debugging
Double_t TGEMSBSSimTrack::MCFitX_print() const
{
  printf("MCFit[0-3]: %f %f %f %f \n", fMCFitPar[0], fMCFitPar[1], fMCFitPar[2], fMCFitPar[3]);
  printf("MCFit[4-8]: %f %f %f %f %f \n", fMCFitPar[4], fMCFitPar[5], fMCFitPar[6], fMCFitPar[7], fMCFitPar[8]);
  // printf("RcFit[0-3]: %f %f %f %f \n", fRcFitPar[0], fRcFitPar[1], fRcFitPar[2], fRcFitPar[3]);
  // printf("RcFit[4-8]: %f %f %f %f %f \n", fRcFitPar[4], fRcFitPar[5], fRcFitPar[6], fRcFitPar[7], fRcFitPar[8]);
  return fMCFitPar[0];
}
// Those below are not useful for SBS, which needs X, Y, Xdir, Ydir (unless otherwise demonstrated)
// refer to comment in TGEMSBSSimEvent.h l. 30-32
/*
//-----------------------------------------------------------------------------
Double_t TGEMSBSSimTrack::MCFitR() const
{
  return TMath::Sqrt(fMCFitPar[0]*fMCFitPar[0] + fMCFitPar[2]*fMCFitPar[2] );
}
//-----------------------------------------------------------------------------
Double_t TGEMSBSSimTrack::MCFitPhi() const
{
  return TVector3( fMCFitPar[0], fMCFitPar[2], 0 ).Phi();
}
//-----------------------------------------------------------------------------
Double_t TGEMSBSSimTrack::MCFitThetaDir() const
{
  return TVector3( fMCFitPar[1], fMCFitPar[3], 1.0 ).Theta();
}
//-----------------------------------------------------------------------------
Double_t TGEMSBSSimTrack::MCFitPhiDir() const
{
  return TVector3( fMCFitPar[1], fMCFitPar[3], 1.0 ).Phi();
}
//-----------------------------------------------------------------------------
Double_t TGEMSBSSimTrack::RcFitR() const
{
  return TMath::Sqrt(fRcFitPar[0]*fRcFitPar[0] + fRcFitPar[2]*fRcFitPar[2] );
}
//-----------------------------------------------------------------------------
Double_t TGEMSBSSimTrack::RcFitPhi() const
{
  return TVector3( fRcFitPar[0], fRcFitPar[2], 0 ).Phi();
}
//-----------------------------------------------------------------------------
Double_t TGEMSBSSimTrack::RcFitThetaDir() const
{
  return TVector3( fRcFitPar[1], fRcFitPar[3], 1.0 ).Theta();
}
//-----------------------------------------------------------------------------
Double_t TGEMSBSSimTrack::RcFitPhiDir() const
{
  return TVector3( fRcFitPar[1], fRcFitPar[3], 1.0 ).Phi();
}
*/

//-----------------------------------------------------------------------------
TGEMSBSECalCluster::TGEMSBSECalCluster()
  : fEnergy(0), fXPos(0), fYPos(0), fTime(0)
{
}

//-----------------------------------------------------------------------------
TGEMSBSECalCluster::TGEMSBSECalCluster(double E, double X, double Y, double t)
  : fEnergy(E), fXPos(X), fYPos(Y), fTime(t)
{
}

//-----------------------------------------------------------------------------
TGEMSBSSimEvent::TGEMSBSSimEvent()
  : fRunID(0), fEvtID(0), fWeight(1.), fMCTracks(0), fNSignal(0),
    fSectorsMapped(false), fSignalSector(0)
{
}

//-----------------------------------------------------------------------------
TGEMSBSSimEvent::TGEMSBSSimEvent( UInt_t ntracks )
  : fRunID(0), fEvtID(0), fWeight(1.), fNSignal(0), fSectorsMapped(false),
    fSignalSector(0)
{
  if( ntracks == 0 ) ntracks = 1;
  fMCTracks = new TClonesArray( "TGEMSBSSimTrack", ntracks );
}

//-----------------------------------------------------------------------------
TGEMSBSSimEvent::~TGEMSBSSimEvent()
{
  delete fMCTracks;
}

//-----------------------------------------------------------------------------
TGEMSBSSimTrack* TGEMSBSSimEvent::AddTrack( Int_t number, Int_t pid,
				      const TVector3& vertex,
				      const TVector3& momentum, 
				      const TVector3& vertexAtTarget, 
				      const TVector3& momentumAtTarget )
{
  // Add a physics track with the given parameters

  return
    new( (*fMCTracks)[GetNtracks()] ) TGEMSBSSimTrack( number, pid,
						    vertex, momentum, vertexAtTarget, momentumAtTarget );
}

//-----------------------------------------------------------------------------
void TGEMSBSSimEvent::Clear( Option_t* opt )
{
  // Clear the event in preparation for reading next tree entry

  TString sopt(opt);

  fNSignal = 0;
  fGEMClust.clear();
  fGEMStrips.clear();
 
  if( sopt.Contains("all",TString::kIgnoreCase) ) {
    fECalClusters.clear();   
    if( fMCTracks ) {
      fMCTracks->Clear(opt);
    }
  }
}

//-----------------------------------------------------------------------------
Int_t TGEMSBSSimEvent::GetNtracks() const
{
  // Get number of physics tracks

  return fMCTracks->GetLast()+1;
}

//-----------------------------------------------------------------------------
void TGEMSBSSimEvent::Print( Option_t* opt ) const
{
  // Print current event info

  cout << ">>>>> =====================================" << endl;
  cout << "Event number:               " << fEvtID << endl;
  cout << "Number of hits:             " << fGEMClust.size()   << endl;
  cout << "Number of fired GEM strips: " << fGEMStrips.size()  << endl;

  TString sopt(opt);
  bool do_all    = sopt.Contains("all",   TString::kIgnoreCase);
  bool do_hit    = sopt.Contains("hit",   TString::kIgnoreCase) || do_all;
  bool do_clust  = sopt.Contains("clust", TString::kIgnoreCase) || do_all;
  bool do_track  = sopt.Contains("track", TString::kIgnoreCase) || do_all;
  // bool do_signal = sopt.Contains("signal", TString::kIgnoreCase);
  // bool do_backgr = sopt.Contains("backgr", TString::kIgnoreCase);
  // bool do_pileup = sopt.Contains("pileup", TString::kIgnoreCase);

  if( do_track && fMCTracks ) {
    for( Int_t i=0; i<GetNtracks(); ++i ) {
      fMCTracks->UncheckedAt(i)->Print(opt);
    }
  }
  if( do_clust ) {
    for( vector<GEMCluster>::const_iterator ic = fGEMClust.begin();
	 ic != fGEMClust.end(); ++ic ) {
      const GEMCluster& c = *ic;
      cout << "hit = " << c.fID
	   << ", src = "    << c.fSource
	   << ", sect = "   << c.fRealSector;
      if( fSectorsMapped ) cout << "->" << c.fSector;
      cout << ", plane = "  << c.fPlane
	   << ", GID = "    << c.fType
	   << ", TRID = "   << c.fTRID
	   << ", PID = "    << c.fPID
	   << ", chrg = "   << c.fCharge
	   << ", time = "   << c.fTime
	   << ", size = ("  << c.fSize[0]  << ", " << c.fSize[1]  << ")"
	   << ", start = (" << c.fStart[0] << ", " << c.fStart[1] << ")"
	   << ", coord = (" << c.fXProj[0] << ", " << c.fXProj[1] << ")"
	   << endl;
      cout << " mom lab  = "; c.fP.Print();
      cout << " mom spec = "; c.fPspec.Print();
      cout << " mcpos  = "; c.fMCpos.Print();
      cout << " hitpos = "; c.fHitpos.Print();
    }
    if( do_hit && !fGEMClust.empty() )
      cout << "-------------" << endl;
  }

  if( do_hit ) {
    UInt_t i = 0;
    for( vector<DigiGEMStrip>::const_iterator is = fGEMStrips.begin();
	 is != fGEMStrips.end(); ++is ) {
      const DigiGEMStrip& s = *is;
      cout << "strip = " << i++
	   << ", sect = "   << s.fSector
	   << ", plane = "  << s.fPlane
	   << ", proj = "   << s.fProj
	   << ", strip = "  << s.fChan
	   << ", type = "   << s.fSigType
	   << ", chrg = "   << s.fCharge
	   << ", adc = ";
      for( int k=0; k<s.fNsamp; k++ ) {
	cout << s.fADC[k];
	if( k+1 != s.fNsamp ) cout << ", ";
      }
      cout << ", hits = ";
      for( Int_t isc = 0; isc < s.fClusters.GetSize(); isc++ ) {
	cout << s.fClusters[isc];
	if( isc+1 != s.fClusters.GetSize() ) cout << ", ";
      }
      cout << endl;
    }
  }
}

//-----------------------------------------------------------------------------
ClassImp(TGEMSBSSimEvent)
ClassImp(TGEMSBSSimTrack)
ClassImp(TGEMSBSECalCluster)
