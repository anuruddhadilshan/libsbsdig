#ifndef _TSBSSIMDIGITIZER_H
#define _TSBSSIMDIGITIZER_H

#include <vector>

class TSBSGeant4File;
class TSBSSimDetector;
class TSBSSimEvent;
class TFile;
class TTree;

class TSBSSimDigitizer {
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
};

#endif // _TSBSSIMDIGITIZER_H

