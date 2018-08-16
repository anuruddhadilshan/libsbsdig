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

#define qe 1.602e-19
#define spe_unit 1.0e-9

#define NBANKS 1

// Tag numbers associated in the GEMC banks
#define __GENERATED_TAG  10
#define __CER_TAG  110

#define __GENERATED_SIZE 7

// List of detector unique IDs: 
// by (proposed) convention: DetUniqueID = DetType*10+DetID
// DetType of type det_type defined in g4sbs_types: kHCal(0), kECal(1), kCher(2), kScint(3), kGEM(4);
#define HCAL_UNIQUE_DETID 0
#define BBPS_UNIQUE_DETID 10
#define BBSH_UNIQUE_DETID 11
#define ECAL_UNIQUE_DETID 12
#define GRINCH_UNIQUE_DETID 20
#define RICH_UNIQUE_DETID 21
#define HODO_UNIQUE_DETID 30
#define CDET_UNIQUE_DETID 31
#define BBGEM_UNIQUE_DETID 40
#define SBSGEM_UNIQUE_DETID 41
#define FT_UNIQUE_DETID 42
#define FPP1_UNIQUE_DETID 43
#define FPP2_UNIQUE_DETID 44
//#define CHER_HIT_ID 0

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


/*
struct SignalInfo{
  int pid;
  int tid;
  // Int_t fillBitsGEM;
  // Int_t fillBitsEC;
  // Int_t signalSector; //used if map sector
  // Double_t ECEDep;
  // Double_t momentum;
  // Double_t R;
  SignalInfo() {}
  SignalInfo(int apid, int atid):pid(apid), tid(atid) {}
  ~SignalInfo() {}
};

struct SpectroInfo{
  double fMCangle;
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
  double  fROImpedance;        // readout impedance
  double  fADCconversion;        // charge/ADC channel conversion
  int     fADCbits;               // number of bits in ADC
  double  fTDCconversion;          // time/TDC channel conversion
  int     fTDCbits;               // number of bits in ADC
  std::vector<double>  fGain;      // Gain 
  std::vector<double>  fPedestal;  // Pedestal value (adc value)
  std::vector<double>  fPedNoise;  // Pedestal noise (adc value)
  std::vector<double>  fThreshold; //Channel threshold
  double  fTriggerJitter;         // trigger jitter (ns)
  double  fTriggerOffset;        // trigger offset (ns)
  double  fGateWidth;           // gate width (ns)
  double  fSPEtau;            // tau param for SPE
  double  fSPEsig;           // sigma param for SPE
  double  fSPEtransittime; // pmt transit time param for SPE
  // just for all those, clearly, we'll have to switch to a class... TODO :S
  double NpeChargeConv(int chan){return Gain(chan)*fROImpedance*qe/spe_unit;}
  double Gain(uint chan){
    if(fGain.size()>1){
      if(fGain.size()<=chan){
	printf("warning: requested channel number %ud larger than channel size %ld. Check code where this method (DigInfo.Gain(int)) is employed and /or database!\n", 
	       chan, fGain.size());
	exit(-1);
      }
      return fGain.at(chan);
    }else{
      return fGain.at(0);
    }
  };
  //
  double Pedestal(uint chan){
    if(fPedestal.size()>1){
      if(fPedestal.size()<=chan){
	printf("warning: requested channel number %ud larger than channel size %ld. Check code where this method (DigInfo.Pedestal(int)) is employed and /or database!\n", 
	       chan, fPedestal.size());
	exit(-1);
      }
      return fPedestal.at(chan);
    }else{
      return fPedestal.at(0);
    }
  };
  double PedestalNoise(uint chan){
    if(fPedNoise.size()>1){
      if(fPedNoise.size()<=chan){
	printf("warning: requested channel number %ud larger than channel size %ld. Check code where this method (DigInfo.PedestalNoise(int)) is employed and /or database!\n", 
	       chan, fPedNoise.size());
	exit(-1);
      }
      return fPedNoise.at(chan);
    }else{
      return fPedNoise.at(0);
    }
  };
  double Threshold(uint chan){
    if(fThreshold.size()>1){
      if(fThreshold.size()<=chan){
	printf("warning: requested channel number %ud larger than channel size %ld. Check code where this method (DigInfo.ThresholdNoise(int)) is employed and /or database!\n", 
	       chan, fThreshold.size());
	exit(-1);
      }
      return fThreshold.at(chan);
    }else{
      return fThreshold.at(0);
    }
  };
  
  
  
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
*/

#endif//__GEMC_TYPES_H
