#ifndef __TSBSDBMANAGER_H
#define __TSBSDBMANAGER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "Rtypes.h"
#include "VarDef.h"
#include "TMath.h"
#include "types.h"

using namespace std;

class TSBSDBManager {
public:
  ~TSBSDBManager();
  static TSBSDBManager* GetInstance() {
    if (fManager == NULL) fManager = new TSBSDBManager();
    return fManager;
  }
  
  void LoadGeneralInfo(const string& fileName);
  void LoadGeoInfo(const string& prefix);
  
  const int    &   GetNDetectors()   { return fNDetectors;   }
  const int    &   Getg4sbsExpType() { return fg4sbsExpType; }
  const int    &   Getg4sbsDetType() { return fg4sbsDetType; }
  const int    &   GetChanPerSlot()  { return fChanPerSlot;  }
  const int    &   GetSlotPerCrate() { return fSlotPerCrate; }
  
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
  
  int    fErrID;
  double fErrVal;
};

#endif
