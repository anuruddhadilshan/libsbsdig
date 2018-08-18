// Example "replay" script
//#define DEBUG 1
#include "TSystem.h"
#include "TDatime.h"
#include "TSBSGeant4File.h"
#include "TSBSSimHCal.h"
#include "TSBSSimScint.h"
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

void digi_scint_test(int nentries = 100)
{
  printf("\n** This gets called with 'analyzer' and not 'root' **\n");
  printf("** If you're getting missing symbol errors, this is likely the cause **\n\n");

  TDatime run_time = 991231;

  gSystem->Load("../libsbsdig.so");

  ////////////////////////////////////////////////////////////////
  
  TSBSDBManager* manager = TSBSDBManager::GetInstance();
  manager->SetDebug(1);
  //manager->LoadGeneralInfo(Form("%s/db_generalinfo_grinch.dat",gSystem->Getenv("DB_DIR")));
  //manager->LoadGeoInfo("g4sbs_grinch");
  manager->LoadGenInfo("db_geninfo_gmn.dat");
  
  // Create the SBS Digitizer (will control the digitization process)
  TSBSSimDigitizer *digitizer = new TSBSSimDigitizer();
  
  // First load the input root file
  TSBSGeant4File *f = new TSBSGeant4File("/work/halla/sbs/efuchey/gmn13.5_elastic_sig_20180709_22/elastic_0.root");
  
  TSBSSimScint *hodo = new TSBSSimScint("hodo", 30);
  hodo->SetDebug(3);
  digitizer->Add(hodo);

  digitizer->Process(f, nentries);
  
  /*

  TSBSSimHCal *hcal = new TSBSSimHCal();
  digitizer->Add(hcal);

  // Now start the digitization on this file


  delete hcal;
  */
  
  cout << "delete hodoscope" << endl;
  
  delete hodo;
}
