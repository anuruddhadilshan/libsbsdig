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

//Replacing these C-style #define statements with C++ const global variable declarations
const double qe = 1.602e-19;
const double spe_unit = 1.0e-9; //to convert ns to s...
const double m_e = 511.e-6; //mass electron in GeV
const double n_lg  = 1.68;    //WTF is this?
const double ROimpedance = 50.0; //Ohm

// List of detector unique IDs: 
// by (proposed) convention: DetUniqueID = DetType*10+DetID
// DetType of type det_type defined in g4sbs_types: kHCal(0), kECal(1), kCher(2), kScint(3), kGEM(4);
// TODO: also put those in the detector DB, and have DB manager ensure no detector share and indentical unique ID
const int  HCAL_UNIQUE_DETID = 0;
const int  BBPS_UNIQUE_DETID = 10;
const int  BBSH_UNIQUE_DETID = 11;
const int  ECAL_UNIQUE_DETID = 12;
const int  GRINCH_UNIQUE_DETID = 20;
const int  RICH_UNIQUE_DETID = 21;
const int  HODO_UNIQUE_DETID = 30;
const int  CDET_UNIQUE_DETID = 31;
const int  ACTIVEANA_UNIQUE_DETID = 32;
const int  PRPOLBS_SCINT_UNIQUE_DETID = 33;
const int  PRPOLFS_SCINT_UNIQUE_DETID = 34;
const int  BBGEM_UNIQUE_DETID = 40;
const int  SBSGEM_UNIQUE_DETID = 41;
const int  FT_UNIQUE_DETID = 42;
const int  FPP1_UNIQUE_DETID = 43;
const int  FPP2_UNIQUE_DETID = 44;
const int  CEPOL_GEMFRONT_UNIQUE_DETID = 45;
const int  CEPOL_GEMREAR_UNIQUE_DETID = 46;
const int  PRPOLBS_GEM_UNIQUE_DETID = 47;
const int  PRPOLFS_GEM_UNIQUE_DETID = 48;

/*
enum exp_type{
  kNeutronExp, kGEp, kGEnRP, 
  kSIDIS, kA1n, kTDIS, kDVCS
};

enum det_type{
  kHCal, kECal,
  kCher, kScint,
  kGEM
};
*/

const std::string kProj_str[2] = {"x", "y"};
//db ??? YES!
//const double bbgem_z[5] = {0.85, 1.0, 1.15, 1.30, 2.33223544};
//const double ft_z[6] = {1.835510, 1.925510, 2.015510, 2.105510, 2.195510, 2.285510};
//const double fpp1_z[5] = {3.05906, 3.15906, 3.25906, 3.35906, 3.45906};
//const double fpp2_z[5] = {4.19731, 4.29731, 4.39731, 4.49731, 4.59731};


#endif//__GEMC_TYPES_H
