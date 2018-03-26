#ifndef __TYPES_
#define __TYPES_

/*
enum GEMDir_t { kGEMX, kGEMY, kGEMR, kGEMPhi };

enum StripSignalType{
  kPrimaryStrip = 0,
  kSecondaryStrip,
  kInducedStrip
};
*/

struct SignalInfo{
  Int_t pid;
  Int_t tid;
  /*
    Int_t fillBitsGEM;
    Int_t fillBitsEC;
    Int_t signalSector; //used if map sector
    Double_t ECEDep;
    Double_t momentum;
    Double_t R;
  */
  SignalInfo() {}
  SignalInfo(Int_t apid, Int_t atid):pid(apid), tid(atid) {}
  ~SignalInfo() {}
};

struct GeoInfo{
  int    fZCkovIn;       // Z of the entrance window in the spectrometer central ray;
  int    fNradiator;     // radiator index of refraction;
  int    fLradiator;     // radiator length on central ray;
  //int    fNquartz;       // quartz window index of refraction;
  int    fNPMTs;         // number of PMTs
  int    fNPMTrows;      // number of PMT rows
  int    fNPMTcolsMax;   // max number of PMT columns 
  double fPMTmatrixHext; // horizontal extension, in m, of the PMT matrix (from lower PMT center to higher PMT center)
  double fPMTmatrixVext; // vertical extension, in m, of the PMT matrix (from left PMT center to right PMT center)
  double fPMTdistX;      // projected X distance between the center of 2 PMT tubes in consecutive rows, in m
  double fPMTdistY;      // Y distance between the center of 2 PMT tubes in consecutive columns, in m
  double fX_TCPMT;       // X position of the top close PMT center in the PMT matrix (transport coord)
  double fY_TCPMT;       // Y position of the top close PMT center in the PMT matrix (transport coord)
};

#endif//__TYPES_
