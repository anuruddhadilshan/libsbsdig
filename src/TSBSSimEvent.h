/////////////////////////////////////////////////////////////////////
//
//   TSBSSimEvent
//
//   Common class definitions (event, MC track, etc.) for SoLID
//   simulation decoder
//
//   Ole Hansen (ole@jlab.org)  December 2011
//
/////////////////////////////////////////////////////////////////////

#ifndef __TSBSSim_h
#define __TSBSSim_h

#include "SimDecoder.h"
#include "TVector3.h"
//#include "TArrayS.h"
#include <vector>
#include <stdint.h>

class TClonesArray;

//-----------------------------------------------------------------------------
class TSBSSimTrack : public Podd::MCTrack {
public:
  TSBSSimTrack( Int_t number, Int_t pid,
		const TVector3& vertex, const TVector3& momentum );
  TSBSSimTrack();
  
  //Special function for debugging
  Double_t MCFitX_print()     const;

  // Accessors for SBS-specific parameters
  // EFuchey: 2017/01/24:
  // Those below are not useful for SBS, which needs X, Y, Xdir, Ydir (unless otherwise demonstrated)
  // Now I'm wondering if TSBSSimTrack is useful at all, but in case, we will stick to that
  /* Double_t MCFitR()     const; */
  /* Double_t MCFitPhi()     const; */
  /* Double_t MCFitThetaDir()  const; */
  /* Double_t MCFitPhiDir()  const; */
  /* Double_t RcFitR()     const; */
  /* Double_t RcFitPhi()     const; */
  /* Double_t RcFitThetaDir()  const; */
  /* Double_t RcFitPhiDir()  const; */
  
  ClassDef(TSBSSimTrack,3)  // A MC physics track in SBS
};

//-----------------------------------------------------------------------------
// Kludgy hardcoded parameters necessary because I can't get ROOT to allocate
// arrays dynamically via TTree::GetEntry

#define treeName "digtree"
#define eventBranchName "event"

class TSBSSimEvent : public TObject {
public:
  TSBSSimEvent();                 // Default constructor, for ROOT I/O
  virtual ~TSBSSimEvent();
  TSBSSimEvent( UInt_t ntracks ); // Construct and initialize fMCTracks
  
  virtual void Clear( const Option_t* opt="" );
  virtual void Print( const Option_t* opt="" ) const;
  TSBSSimTrack* AddTrack( Int_t number, Int_t pid,
			  const TVector3& vertex, const TVector3& momentum );
  
  Int_t GetNPMTHits()  const { return fPMTHits.size(); }
  Int_t GetNtracks() const;
  
  // Event identification
  Int_t     fRunID;               // Run number
  Int_t     fEvtID;               // Event number

  Double_t  fWeight;              // Event weight
  Int_t     fNSignal;             // Number of clusters from trigger track (signal)
  
  // MC tracks
  TClonesArray*   fMCTracks;      //-> Physics tracks
  //std::vector< std::pair< Int_t, std::vector<Short_t> > > fMCClusterHitID;
  // collections of hit IDs coming from the same MC track ID
  // save only clusters with size >= 2

  struct PMTHit {
    // MC hit data
    Short_t  fID;          // Hit number
    Int_t    fSource;      // MC data set source (0 = signal, >0 background)
    Short_t  fType;        // GEANT particle type (1 = primary)
    Int_t    fMCtrackID;   // GEANT track ID
    Int_t    fMCtrackPID;  // GEANT particle ID (defined only if particle is primary)
    Short_t  fOrigVolFlag; // 
    Float_t  fXPMT;        // X coordinate of the PMT in transport coordinates
    Float_t  fYPMT;        // Y coordinate of the PMT in transport coordinates
    Float_t  fNpe;         // Number of photoelectrons
    Double_t fTime;        // Arrival time at electronics
    Double_t fTDCtime[2];  // TDC rising and falling times values wrt trigger time
    // Digitization results for this hit
    Short_t  fDetID;       // Detector ID
    Short_t  fChannel;     // Channel number
    Short_t  fPMTrow;      // Row number: cross reference to Channel number
    Short_t  fPMTcol;      // Column number: cross reference to Channel number
    Short_t  fVETROCID;    // VETROC ID
    uint32_t fTDC[2];      // TDC VETROC words
  };
  std::vector<PMTHit> fPMTHits;  // All MC-generated clusters in the GEMs

  struct DetectorData {
    Short_t fDetID;  // Source detector number
    Short_t fChannel; // Channel number for this hit
    std::vector<uint32_t> fData;
  };
  std::vector<DetectorData> fDetectorData;

  ClassDef(TSBSSimEvent, 6) // Simulated data for one event
};

#endif
