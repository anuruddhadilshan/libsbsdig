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
#define spe_unit 1.0e-9 //to convert ns to s...

#define NBANKS 1

// Tag numbers associated in the GEMC banks

#define __GENERATED_SIZE 7

// List of detector unique IDs: 
// by (proposed) convention: DetUniqueID = DetType*10+DetID
// DetType of type det_type defined in g4sbs_types: kHCal(0), kECal(1), kCher(2), kScint(3), kGEM(4);
// TODO: also put those in the detector DB, and have DB manager ensure no detector share and indentical unique ID
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

enum exp_type{
  kGMn, kGEp, kGEn, 
  kSIDIS, kA1n, kTDIS, kDVCS
};

enum det_type{
  kHCal, kECal,
  kCher, kScint,
  kGEM
};

#endif//__GEMC_TYPES_H
