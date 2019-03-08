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
// Kludgy hardcoded parameters necessary because I can't get ROOT to allocate
// arrays dynamically via TTree::GetEntry

#define treeName "digtree"
#define eventBranchName "event"

class TSBSSimEvent : public TObject {
public:
  TSBSSimEvent();                 // Default constructor, for ROOT I/O
  virtual ~TSBSSimEvent(){};
  
  virtual void Clear( const Option_t* opt="" );
  virtual void Print( const Option_t* opt="" ) const;
  
  // Event identification
  Int_t     fRunID;               // Run number
  Int_t     fEvtID;               // Event number

  Double_t  fWeight;              // Event weight
  Int_t     fNSignal;             // Number of clusters from trigger track (signal)
  
  
  struct SimDetectorData {
    Short_t fDetID;  // Source detector number
    Short_t fChannel; // Channel number for this hit
    Short_t fDataType; // Data type for this hit: 0: Npe, 1: Time, 2: SumEdep
    Short_t fNdata;  // number of data stored
    std::vector<double> fData;
  };
  struct DetectorData {
    UShort_t fDetID;  // Source detector number
    UShort_t fChannel; // Channel number for this hit
    //not sure we want to use that yet - there might be a reason why it is what it is
    // Short_t fDataType; // Data type for this hit: 0: Npe, 1: Time, 2: SumEdep
    // Short_t fNdata;  // number of data stored
    std::vector<uint32_t> fData;
  };
  std::vector<DetectorData> fDetectorData;
  std::vector<SimDetectorData> fSimDetectorData;
  /**/
  Int_t NSimDetData;
  std::vector<Short_t> SimDetID;
  std::vector<Short_t> SimDetChannel;
  std::vector<Short_t> SimDetDataType;
  std::vector<Short_t> SimDetNData;
  std::vector< std::vector<double> > SimDetData;
  
  Int_t NDetData;
  std::vector<Short_t> DetID;
  std::vector<Short_t> DetChannel;
  //std::vector<Short_t> DetDataType;
  std::vector<Short_t> DetNData;
  std::vector< std::vector<uint32_t> > DetData;
  
  ClassDef(TSBSSimEvent, 6) // Simulated data for one event
};

#endif
