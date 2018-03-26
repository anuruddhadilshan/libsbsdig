#ifndef __TSBSDET_H
#define __TSBSDET_H

//#include "THaDetector.h"
#include "THaCherenkov.h"

class THaEvData;

// TSBSDet also inherits form THaDetector, which grants it all the functions from its class
// (see http://hallaweb.jlab.org/podd/doc/html_v16/ClassIndex.html for more info).

class TSBSDet : public THaDetector {
 public:
  //Constructors and destructor
  TSBSDet(const char *name, const char *desc);//It is recommended to use this constructor
  TSBSDet() : fNChannels(0), fZPos(0) {} // for ROOT RTTI

  virtual ~TSBSDet();
  
  //Read the geometry for the TSBSBox in the data base
  Int_t ReadDatabase (const TDatime& date);
  Int_t ReadGeometry (FILE* file, const TDatime& date, Bool_t required);// This class *must* be overloaded 
  const char* GetDBFileName() const;

  Int_t Decode( const THaEvData & );// This class *must* be overloaded 
  
  // Access to the info of TSBSGEMPlane which is regarded as a subdetector of TSBSDet.
  // (see comments in the code of class TSBSGEMPlane)
  UInt_t GetNChannels() const {return fNChannels;};
  Double_t GetZpos() const {return fZPos;};
  
  void Print ();// This class *must* be overloaded 
  
 protected:
  UInt_t   fNChannels; // number of channels
  Double_t fZPos;      // position of the detector from the spectrometer bending point
  
  ClassDef(TSBSDet,0)

};

#endif//__TSBSCHER_H
