// Example "replay" script
//#define DEBUG 1
#include "TSystem.h"
#include "TDatime.h"
#include "TSBSGeant4File.h"
#include "TSBSSimHCal.h"
#include "TSBSSimScint.h"
#include "TSBSSimCher.h"
#include "TSBSDBManager.h"
#include "TSBSSimDigitizer.h"
#include "THaAnalysisObject.h"

/*
#include "TSBSCher.h"
#include "TSystem.h"
#include "TSBSSpec.h"
#include "TSBSSimCherDigitization.h"
#include "TSBSDBManager.h"
#include "TSBSGeant4File.h"
#include "TSBSSimEvent.h"
*/

R__LOAD_LIBRARY(../libsbsdig)

void digi_all_test(int nentries = 100, int debuglevel = 1)
{
  printf("\n** This gets called with 'analyzer' and not 'root' **\n");
  printf("** If you're getting missing symbol errors, this is likely the cause **\n\n");

  TDatime run_time = 991231;

  gSystem->Load("../libsbsdig.so");

  ////////////////////////////////////////////////////////////////
  
  TSBSDBManager* manager = TSBSDBManager::GetInstance();
  manager->SetDebug(debuglevel);
  //manager->LoadGeneralInfo(Form("%s/db_generalinfo_grinch.dat",gSystem->Getenv("DB_DIR")));
  //manager->LoadGeoInfo("g4sbs_grinch");
  manager->LoadGenInfo("db_geninfo_gmn.dat");
  
  // Create the SBS Digitizer (will control the digitization process)
  TSBSSimDigitizer *digitizer = new TSBSSimDigitizer("digitized/simdig_test.root");
  digitizer->SetDebug(debuglevel+1);
  
  // First load the input root file
  TSBSGeant4File *f = new TSBSGeant4File("/work/halla/sbs/efuchey/gmn13.5_elastic_sig_20180709_22/elastic_0.root");
  // TSBSGeant4File *f_b = new TSBSGeant4File("/volatile/halla/sbs/efuchey/gmn13.5_beam_bkgd_20180718_14/beam_bkgd_0.root");
  
  // digitizer->AddInputFile(f, 1);
  // digitizer->AddInputFile(f_b, 100);
  
  // It is recommended  to declare the detector with its unique ID (second parameter)
  // See list of unique det IDs defined in src/g4sbs_types.h
  TSBSSimScint *hodo = new TSBSSimScint("hodo", 30);
  hodo->SetDebug(debuglevel);
  digitizer->AddDetector(hodo);
  
  TSBSSimScint *cdet = new TSBSSimScint("cdet", 31);
  cdet->SetDebug(debuglevel);
  digitizer->AddDetector(cdet);
  
  TSBSSimCher *grinch = new TSBSSimCher("grinch", 20);
  grinch->SetDebug(debuglevel);
  digitizer->AddDetector(grinch);
  
  TSBSSimECal *ps = new TSBSSimECal("ps", 10);
  ps->SetDebug(debuglevel);
  digitizer->AddDetector(ps);
  
  TSBSSimECal *sh = new TSBSSimECal("sh", 11);
  sh->SetDebug(debuglevel);
  digitizer->AddDetector(sh);
  
  digitizer->Process(f, nentries);
  //digitizer->Process(nentries);
  
  /*
    TSBSSimHCal *hcal = new TSBSSimHCal();
    digitizer->Add(hcal);
    // Now start the digitization on this file
    delete hcal;
  */
  
  //cout << "delete detectors" << endl;
  //delete hodo;
  //delete cdet;
  //delete grinch;
}
