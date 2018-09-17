#ifndef __TSBSDBMANAGER_H
#define __TSBSDBMANAGER_H

#include <iostream>
#include <fstream>
#include "THaAnalysisObject.h"
#include <vector>
#include <string>
#include <sstream>
#include "Rtypes.h"
#include "VarDef.h"
#include "TMath.h"
#include "g4sbs_types.h"
#include "TSBSSimAuxi.h"

using namespace std;

class TSBSDBManager : public THaAnalysisObject {
public:
  ~TSBSDBManager();
  static TSBSDBManager* GetInstance() {
    if (fManager == NULL) fManager = new TSBSDBManager();
    return fManager;
  }
  
  Int_t LoadGenInfo(const string& fileName);
  Int_t LoadDetInfo(const string& specname, const string& detname);
  
  exp_type GetExpType(){return(fSBSExpType);};
  
  // This is left as a patch, to make the program compile, but we'll need to replace them eventually,
  // cause we'll need these information for each detector.
  // const int    &   GetChanPerSlot()  { return 0; }//fChanPerSlot;  }
  // const int    &   GetSlotPerCrate() { return 0; }//fSlotPerCrate; }
  
  const TDetInfo & GetDetInfo(const char* detname);
  bool IsDetInfoAvailable(const char* detname); // Check if a detector is defined
  
 protected:
  TSBSDBManager();
  
  //New database parameters:
  exp_type fSBSExpType;
  UInt_t fNSpecs;
  UInt_t fNextGeneratedCrate;
  std::vector<string> fSpecNames;
  std::vector<TSpectroInfo> fSpectroInfo;
  std::vector<TDetInfo> fDetInfo;

  // int    LoadDB(ifstream& inp, DBRequest* request, const string& prefix);
  // string FindKey( ifstream& inp, const string& key );
  
  static TSBSDBManager* fManager;
  
  /*
  const int    &   GetNDetectors()   { return fNDetectors;   }
  const int    &   Getg4sbsExpType() { return fg4sbsExpType; }
  const int    &   Getg4sbsDetType() { return fg4sbsDetType; }
  
  const int    &   GetSigPID(unsigned int i);
  const int    &   GetSigTID(unsigned int i);
  
  const double &   GetZCkovIn(int i);
  const double &   GetNradiator(int i);
  const double &   GetLradiator(int i);
  const int    &   GetNPMTs(int i);
  const int    &   GetNPMTrows(int i);
  const int    &   GetNPMTcolsMax(int i);
  const double &   GetPMTmatrixHext(int i);
  const double &   GetPMTmatrixVext(int i);
  const double &   GetPMTdistX(int i);
  const double &   GetPMTdistY(int i);
  const double &   GetX_TCPMTs(int i);
  const double &   GetY_TCPMTs(int i);
    
 protected:
  TSBSDBManager();
  int    LoadDB(ifstream& inp, DBRequest* request, const string& prefix);
  string FindKey( ifstream& inp, const string& key );
  bool   CheckIndex(int i);
  
  static TSBSDBManager* fManager;
  
  //variable for data base information
  int fNDetectors;  // number of Cherenkov detectors in arm (usually 1...)
  int fChanPerSlot;  // number of PMTs per VETROC
  int fSlotPerCrate;  // number of VETROC per crate
  
  // Parameters for TSBSGeant4File
  int fg4sbsExpType;// experiment flag. Choices are
  // 1 - GMn;
  // 2 - GEn;
  // 3 - GEp;
  // 4 - SIDIS;
  // 5 - A1n;
  int fg4sbsDetType;// flag to determine which type of GEM should be read. Choices are:
  // 1 - GRINCH
  // 2 - RICH
  
  // Parameters for signal particles
  int fNSigParticle; // number of signal particles
  vector<int>    fSigPID;
  vector<int>    fSigTID;
  
  //map< int, vector<GeoInfo> > fGeoInfo;
  vector<GeoInfo> fGeoInfo;
  */
  
  int    fErrID;
  double fErrVal;
  
  ClassDef(TSBSDBManager,1)
};

#endif
