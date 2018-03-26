// Example "replay" script
//#define DEBUG 1
#include "TSystem.h"
#include "TDatime.h"
#include "TSBSGeant4File.h"
#include "TSBSSimHCal.h"
#include "TSBSDBManager.h"
#include "TSBSSimDigitizer.h"

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

void digi_hcal_test()
{
  printf("\n** This gets called with 'analyzer' and not 'root' **\n");
  printf("** If you're getting missing symbol errors, this is likely the cause **\n\n");

  TDatime run_time = 991231;

  gSystem->Load("../libsbsdig.so");

  ////////////////////////////////////////////////////////////////

  TSBSDBManager* manager = TSBSDBManager::GetInstance();
  manager->LoadGeneralInfo(Form("%s/db_generalinfo_grinch.dat",gSystem->Getenv("DB_DIR")));
  manager->LoadGeoInfo("g4sbs_grinch");


  // Create the SBS Digitizer (will control the digitization process)
  TSBSSimDigitizer *digitizer = new TSBSSimDigitizer();

  // First load the input root file
  TSBSGeant4File *f = new TSBSGeant4File("data/sbsin_931.root");

  TSBSSimHCal *hcal = new TSBSSimHCal();
  digitizer->Add(hcal);

  // Now start the digitization on this file
  digitizer->Process(f);


  delete hcal;
}
