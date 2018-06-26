#ifndef __G4SBS_TYPES_H
#define __G4SBS_TYPES_H

#include <vector>
#include <string>

////////////////////////////////////////////////////////
//  Data for extracting things from GEMC
//
//  we'll hardcode them here, but it would be nice to
//  maybe get them into a database
//  I guess this could also be done through a mysql
//  interface, but I think that makes it more complicated
//  and breakable


#define NBANKS 1

// Tag numbers associated in the GEMC banks
#define __GENERATED_TAG  10
#define __CER_TAG  110

#define __GENERATED_SIZE 7

static int __g4sbs_types_datasize[NBANKS] = {21};

// FIXME:  Need to do this better,  map?
static int data_size(int tag){
    if( tag == __CER_TAG ) {
	return __g4sbs_types_datasize[0];
    }
    return 0;
}

enum exp_type{
  kGMn, kGEp, kGEn, 
  kSIDIS, kA1n, kTDIS, kDVCS
};

enum det_type{
  kHCal, kECal,
  kCher, kScint,
  kGEM
};

struct SignalInfo{
  int pid;
  int tid;
  /*
    Int_t fillBitsGEM;
    Int_t fillBitsEC;
    Int_t signalSector; //used if map sector
    Double_t ECEDep;
    Double_t momentum;
    Double_t R;
  */
  SignalInfo() {}
  SignalInfo(int apid, int atid):pid(apid), tid(atid) {}
  ~SignalInfo() {}
};

struct SpectroInfo{
  std::vector<SignalInfo> MCsignalInfo;
  int fNDets;
  std::vector<std::string> fDetNames;
  
  SpectroInfo() {}
  SpectroInfo(int ndets):fNDets(ndets) {}
  ~SpectroInfo() {}
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

#endif//__GEMC_TYPES_H
