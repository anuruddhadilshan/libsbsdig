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
  TSBSSimDigitizer();
  virtual ~TSBSSimDigitizer();

  // Parsses the passed file event by event and digitizes
  int Process(TSBSGeant4File *file, int max_events = 0);

  // Add a new detector to the list
  int Add(TSBSSimDetector* detector);
private:
  std::vector<TSBSSimDetector*> fDetectors;
  TSBSSimEvent *fEvent;
  TFile *fOutFile;
  TTree *fOutTree;
  
  ClassDef(TSBSSimDigitizer,1)
};

#endif // _TSBSSIMDIGITIZER_H

