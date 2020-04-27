
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
#include "TSBSSimData.h"
//#include "TArrayS.h"
#include <vector>
#include <stdint.h>

class TClonesArray;
class TSBSDBManager;
class simdig_outdata;
class simgemdig_outdata;

//-----------------------------------------------------------------------------
// Kludgy hardcoded parameters necessary because I can't get ROOT to allocate
// arrays dynamically via TTree::GetEntry

#define treeName "digtree"
#define eventBranchName "event"

class TSBSSimEvent : public TObject {
public:
  TSBSSimEvent();                 // Default constructor, for ROOT I/O
  TSBSSimEvent(int run, int evt, int signal = 0);                 // Default constructor, for ROOT I/O
  virtual ~TSBSSimEvent(){};
  
  virtual void Clear( const Option_t* opt="" );
  virtual void Print( const Option_t* opt="" ) const;
  
  // Event identification
  Int_t     fRunID;               // Run number
  Int_t     fEvtID;               // Event number

  Double_t  fWeight;              // Event weight
  Int_t     fNSignal;             // Number of clusters from trigger track (signal)
  
  /*
  struct SimDetectorData {
    //Short_t fDetID;  // Source detector number
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

  std::map< std::string, Int_t >         NSimDetHits;
  std::map< std::string, std::vector<Short_t> >  SimDetSource;
  std::map< std::string, std::vector<Short_t> >  SimDetChannel;
  std::map< std::string, std::vector<Double_t> > SimDetEdep;
  std::map< std::string, std::vector<Int_t> >    SimDetNpe;
  std::map< std::string, std::vector<Double_t> > SimDetTime;
  std::map< std::string, std::vector<Double_t> > SimDetLeadTime;
  std::map< std::string, std::vector<Double_t> > SimDetTrailTime;
  */
  
  std::map< std::string, trackmchit_outdata > fTrackMCHitOutData;
  std::map< std::string, simhitmc_outdata > fSimHitMCOutData;
  std::map< std::string, simgemhitmc_outdata > fSimGEMHitMCOutData;
  std::map< std::string, simdig_outdata > fSimDigOutData;
  std::map< std::string, simdigsamp_outdata > fSimDigSampOutData;
  
  TSBSDBManager *fManager;
  
  ClassDef(TSBSSimEvent, 6) // Simulated data for one event
};

#endif
