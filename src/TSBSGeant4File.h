#ifndef __TSBSGEANT4FILE_H
#define __TSBSGEANT4FILE_H

// Put prototypes here first so that it doens't freak out
#ifndef __CINT__
//#include "evioUtil.hxx"
#endif//__CINT__

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "g4sbs_tree.h"
#include "TSBSCherData.h"
#include "TSBSDBManager.h"


class g4sbshitdata;
class g4sbsgendata;

///////////////////////////////////////////////////////////////////////////////


class TSBSGeant4File {

 public:
  //constructor may be inputed a data file to input some of the paramaters used by this class
  //NB: if the second file path does not select a valid file, default parameters will be used.
  TSBSGeant4File();// Default constructor
  TSBSGeant4File( const char *name);// Constructor with input file name: recommanded
  virtual ~TSBSGeant4File();// Default destructor

  //void ReadGasData(const char* filename); // NB: See comment lines 128-129 
  
  // Standard getters and setters
  void  SetFilename( const char *name );
  void  SetSource( Int_t i ) { fSource = i; }
  void  Clear();
  Int_t Open();
  Int_t Close();
  
  Long64_t GetEntries(){
    return fTree->GetEntries();
  };
  
  const char* GetFileName() const { return fFilename; }
  Int_t GetSource() const { return fSource; }
  
  // This is actually where the data is read: 
  Int_t ReadNextEvent(int d_flag = 0);

  // Returns a pointer to the data vector
  std::vector<g4sbshitdata*>& GetDataVector() { return fg4sbsHitData; }
  
  //return the size of the hit arrays
  UInt_t GetNData() const { return fg4sbsHitData.size(); }
  UInt_t GetNGen() const { return fg4sbsGenData.size(); }
  
  UInt_t GetEvNum() const { return fEvNum; }
  
  //get one hit from the hit data arrays
  g4sbshitdata *GetHitData(Int_t i) const { return fg4sbsHitData[i]; }
  g4sbsgendata *GetGenData(Int_t i) const { return fg4sbsGenData[i]; }
  
  //get GEM data
  TSBSCherData *GetCherData();
  void GetCherData(TSBSCherData* gd);
  
 private:
  // Members
  char  fFilename[255];
  TFile *fFile;
  g4sbs_tree *fTree;// needed to easily unfold root file data
  Int_t fSource;   // User-defined source ID (e.g. MC run number)  
  
  //hit data arrays
  vector<g4sbshitdata *> fg4sbsHitData;
  vector<g4sbsgendata *> fg4sbsGenData;
  
  unsigned int fEvNum;// global event incrementer

  TSBSDBManager *fManager;
};



#endif//__TSOLSIMG4SBSFILE_H
