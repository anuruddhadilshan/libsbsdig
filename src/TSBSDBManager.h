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
  
  const std::vector<TDetInfo>& GetAllDetInfo() { return fDetInfo; }
  const std::vector<TSpectroInfo>& GetAllSpectroInfo() { return fSpectroInfo; }
  const TSpectroInfo & GetSpectroInfo(const char* specname);
  const TDetInfo & GetDetInfo(const char* detname);
  const TDetInfo& GetDetInfoById(Int_t id);
  bool IsDetInfoAvailable(const char* detname); // Check if a detector is defined
  bool IsDetInfoAvailableById(Int_t id); // Check if a detector is defined by ID
  
  double GetBkgdSpreadTimeWindowHW(){return fBkgdSpreadTimeWindowHW;};
  double GetTriggerJitter(){return fTriggerJitter;};
  //void SetBkgdSpreadTimeWindowHW(int bstwhw){fBkgdSpreadTimeWindowHW = bstwhw;};
  
  //stuff for GEMs...
  
 protected:
  TSBSDBManager();
  
  //New database parameters:
  exp_type fSBSExpType;
  UInt_t fNSpecs;
  std::vector<string> fSpecNames;
  std::vector<TSpectroInfo> fSpectroInfo;
  std::vector<TDetInfo> fDetInfo;

  static TSBSDBManager* fManager;
    
  int    fErrID;
  double fErrVal;
  
  Double_t fBkgdSpreadTimeWindowHW;
  Double_t fTriggerJitter;
  
  TRndmManager *fRN;
  
  //stuff for GEMs...
  
  
  ClassDef(TSBSDBManager,1)
};

#endif
