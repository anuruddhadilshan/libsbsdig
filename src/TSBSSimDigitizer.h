#ifndef _TSBSSIMDIGITIZER_H
#define _TSBSSIMDIGITIZER_H

#include "THaAnalysisObject.h"
#include <vector>

class TSBSGeant4File;
class TSBSSimDetector;
class TSBSSimEvent;
class TFile;
class TTree;
class THaAnalysisObject;

class TSBSSimDigitizer : public THaAnalysisObject {
public:
  TSBSSimDigitizer(const char* outputfilename = "digitized/simdig_test.root");
  virtual ~TSBSSimDigitizer();

  // Parsses the passed file event by event and digitizes
  int Process(TSBSGeant4File *file, int max_events = 0);
  // File superposition: Procession of a stack of files ?
  // make the file stack a member of TSBSSimDigitizer, and just "Process" it with function below 
  // -> not functional yet :/
  int Process(int max_events = 0);// Process the mmeber file stack
  
  // Add a new detector to the list
  void AddDetector(TSBSSimDetector* detector);

  // Add a new file to the file stack
  void AddInputFile(TSBSGeant4File* file, UInt_t weight = 1);
  
private:
  std::vector<TSBSSimDetector*> fDetectors;
  TSBSSimEvent *fEvent;
  TFile *fOutFile;
  TTree *fOutTree;
  
  // Lists of files with weight to perform additive digitization.
  // Files and weights are added at the same time, and cannot be accessed from the outside... 
  std::vector< TSBSGeant4File* > fG4FileStack_;
  // vector of vector if strings: 
  // the global vector contains the type of file, the inner vector contains the list of names of files.
  std::vector< std::vector<TString> > fG4FileStack;
  std::vector< UInt_t > fG4FileWeights;
  
  
  ClassDef(TSBSSimDigitizer,1)
};

#endif // _TSBSSIMDIGITIZER_H

