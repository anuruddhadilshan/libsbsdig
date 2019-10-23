#ifndef _TSBSSIMDIGITIZER_H
#define _TSBSSIMDIGITIZER_H

#include "THaAnalysisObject.h"
#include <vector>
#include <set>
#include <map>
//#include "TRandom3.h"

class TSBSGeant4File;
class TSBSSimDetector;
class TSBSSimEvent;
class TFile;
class TChain;
class TTree;
class THaAnalysisObject;
class TSBSDBManager;
class TRndmManager;

class TSBSSimDigitizer : public THaAnalysisObject {
public:
  TSBSSimDigitizer(const char* outputfilename = "digitized/simdig_test.root");
  virtual ~TSBSSimDigitizer();

  // Parsses the passed file event by event and digitizes
  int AddFileToEvent(TSBSGeant4File *file);
  // File superposition: Procession of a stack of files ?
  // make the file stack a member of TSBSSimDigitizer, and just "Process" it with function below 
  // -> not functional yet :/
  int Process(ULong_t max_events = -1);// Process the signal chain
  
  // Add a new detector to the list
  void AddDetector(TSBSSimDetector* detector);

  // Add a new file to the file stack
  //void AddInputFile(TSBSGeant4File* file, Int_t weight = 1);
  void AddInputFile(const char* filename, Int_t source, Int_t weight = 1);
  
private:
  std::vector<TSBSSimDetector*> fDetectors;
  TSBSSimEvent *fEvent;
  TFile *fOutFile;
  TTree *fOutTree;
    
  // maps of files and weights with source to perform additive digitization.
  std::vector<Int_t> fSources;
  std::map<Int_t, TChain*> fSourceChainMap;
  std::map<Int_t, Int_t> fSourceWeightMap;
  
  TSBSDBManager *fManager;
  
  TRndmManager *fRN;
  
  ClassDef(TSBSSimDigitizer,1)
};

#endif // _TSBSSIMDIGITIZER_H

