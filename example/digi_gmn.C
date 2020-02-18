// Example "replay" script
//#define DEBUG 1
// Only include these headers if we are compiling this macro
// (there is a strange dictionary behavior otherwise...)
#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "TSystem.h"
#include "TDatime.h"
#include "TSBSSimECal.h"
#include "TSBSGeant4File.h"
#include "TSBSSimHCal.h"
#include "TSBSSimScint.h"
#include "TSBSSimCher.h"
#include "TSBSSimGEM.h"
#include "TSBSDBManager.h"
#include "TSBSSimDigitizer.h"
#include "THaAnalysisObject.h"
#include "TSBSSimGEM.h"
#endif

//a bit more complex macro, where you have to provide the paths to the g4sbs files via  input_sigfile and input_bkgdfile (if you want background)

void digi_gmn(ULong64_t nentries, const char* input_sigfile, int nbkgd = 0, const char* input_bkgdfile = "", int debuglevel = 1)
{
  printf("\n** This gets called with 'analyzer' and not 'root' **\n");
  printf("** If you're getting missing symbol errors, this is likely the cause **\n\n");

  TDatime run_time = 991231;

  gSystem->AddDynamicPath("${SBS_ANALYSIS}");
  gSystem->Load("../libsbsdig.so");

  ////////////////////////////////////////////////////////////////
  
  TSBSDBManager* manager = TSBSDBManager::GetInstance();
  manager->SetDebug(debuglevel);
  
  if(debuglevel>=1)cout << "About to read database " << endl;

  manager->LoadGenInfo("db_geninfo_gmn.dat");
  
  if(debuglevel>=1)cout << "Setup digitizer " << endl;
  
  // Create the SBS Digitizer (will control the digitization process)
  TSBSSimDigitizer *digitizer = new TSBSSimDigitizer("digitized/simdig_test.root");
  digitizer->SetDebug(debuglevel);
  
  if(debuglevel>=1)cout << "Setup input file " << endl;
  
  ifstream sig_inputfile(input_sigfile);
  TString currentline;
  while( currentline.ReadLine(sig_inputfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      digitizer->AddInputFile(currentline.Data(), 0, 1);
    }
  }
  
  if(nbkgd){
    ifstream beam_inputfile(input_bkgdfile);
    while( currentline.ReadLine(beam_inputfile) && !currentline.BeginsWith("endlist") ){
      if( !currentline.BeginsWith("#") ){
	digitizer->AddInputFile(currentline.Data(), 1, -nbkgd);
      }
    }
  }
  
  // It is recommended  to declare the detector with its unique ID (second parameter)
  // See list of unique det IDs defined in src/g4sbs_types.h
  
  if(debuglevel>=1)cout << "Declare detectors and add them to digitizer " << endl;
  
  TSBSSimHCal *hcal = new TSBSSimHCal("hcal", 0);
  hcal->SetDebug(debuglevel);
  digitizer->AddDetector(hcal);
  
  TSBSSimScint *hodo = new TSBSSimScint("hodo", 30);
  hodo->SetDebug(debuglevel);
  digitizer->AddDetector(hodo);
  
  TSBSSimCher *grinch = new TSBSSimCher("grinch", 20);
  grinch->SetDebug(debuglevel);
  digitizer->AddDetector(grinch);
  
  TSBSSimECal *ps = new TSBSSimECal("ps", 10);
  ps->SetDebug(debuglevel);
  digitizer->AddDetector(ps);
  
  TSBSSimECal *sh = new TSBSSimECal("sh", 11);
  sh->SetDebug(debuglevel);
  digitizer->AddDetector(sh);
  
  TSBSSimGEM *bbgem = new TSBSSimGEM("gem", 40);
  bbgem->SetDebug(debuglevel);
  digitizer->AddDetector(bbgem);
  
  if(debuglevel>=1)cout << "About to process digitization for " << nentries << "events " << endl;
  
  digitizer->Process(nentries);
}
