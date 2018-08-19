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
  TSBSSimDigitizer *digitizer = new TSBSSimDigitizer();
  digitizer->SetDebug(debuglevel);
  
  // First load the input root file
  TSBSGeant4File *f = new TSBSGeant4File("/work/halla/sbs/efuchey/gmn13.5_elastic_sig_20180709_22/elastic_0.root");
  
  // It is recommended  to declare the detector with its unique ID (second parameter)
  // See list of unique det IDs defined in src/g4sbs_types.h
  TSBSSimScint *hodo = new TSBSSimScint("hodo", 30);
  hodo->SetDebug(debuglevel);
  digitizer->Add(hodo);
  
  /*
  TSBSSimScint *cdet = new TSBSSimScint("cdet", 31);
  hodo->SetDebug(debuglevel);
  digitizer->Add(cdet);
  
  TSBSSimCher *grinch = new TSBSSimCher("grinch", 20);
  hodo->SetDebug(debuglevel);
  digitizer->Add(grinch);
  */

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
