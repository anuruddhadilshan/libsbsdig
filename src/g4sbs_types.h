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
  // I don't think this can fit all detectors. 
  // Might have to do a class or something so it can be inherited ?
  int     fNrows;      // number of rows
  int     fNcols;      // number of columns
  double  fXsize;      // detector X size (in transport coordinates)
  double  fYsize;      // detector Y size (in transport coordinates)
  double  fZpos;       // detector position on spectrometer axis
};

struct DigInfo{
  // same remark as for GeoInfo
  double  fROImpedance;   // readout impedance
  double  fGain;          // Gain 
  double  fPedestal;      // Pedestal value (adc value)
  double  fPedNoise;      // Pedestal noise (adc value)
  double  fTriggerJitter; // trigger jitter (ns)
  double  fTriggerOffset; // trigger offset (ns)
  double  fGateWidth;     // gate width (ns)
};

struct DetInfo{
  std::string fDetName;      // Detector name
  det_type     fDetType;      // DetectorType
  int           fNChan;        // Total number of channels over all detector
  int            fChanPerSlot;  // Number of channels per slot
  int             fSlotPerCrate; // Number of slots per crate
  int              fNPlanes;      // Number of planes // useful e.g. GEM, CDet
  std::vector<int>  fNModules;     // Number of modules per plane // useful e.g. GEM, CDet
  
  std::vector<GeoInfo> fGeoInfo;
  DigInfo fDigInfo;
  DetInfo()
  {
    fNModules.clear();
    fGeoInfo.clear();
  }
  DetInfo(const std::string detname)
  {
    fDetName = detname;
    fNModules.clear();
    fGeoInfo.clear();
  }
  ~DetInfo(){}

};

#endif//__GEMC_TYPES_H
