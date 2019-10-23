// Example "replay" script
//#define DEBUG 1
#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "TSystem.h"
#include "TDatime.h"
#include "TSBSGeant4File.h"
#include "TSBSSimHCal.h"
#include "TSBSSimECal.h"
#include "TSBSSimScint.h"
#include "TSBSSimCher.h"
#include "TSBSDBManager.h"
#include "TSBSSimDigitizer.h"
#include "THaAnalysisObject.h"
#endif

void digi_all_test(int nentries = 100, int nbkgd = 0, int debuglevel = 1)
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

  //manager->LoadGeneralInfo(Form("%s/db_generalinfo_grinch.dat",gSystem->Getenv("DB_DIR")));
  //manager->LoadGeoInfo("g4sbs_grinch");
  manager->LoadGenInfo("db_geninfo_gmn.dat");
  
  if(debuglevel>=1)cout << "Setup digitizer " << endl;
  
  // Create the SBS Digitizer (will control the digitization process)
  TSBSSimDigitizer *digitizer = new TSBSSimDigitizer("digitized/simdig_test.root");
  digitizer->SetDebug(debuglevel);
  
  if(debuglevel>=1)cout << "Setup input file " << endl;
  
  // First load the input root file
  // TSBSGeant4File *f = new TSBSGeant4File("/volatile/halla/sbs/efuchey/gmn13.5_elastic_sig_20190725_15/elastic_0.root");
  // f->SetSource(0);
  // if(!f->Open()){
  //   exit(-1);
  // }
  // if(debuglevel>=2)cout << "Add to digitizer file " << f->GetName() << endl;
  // digitizer->AddInputFile(f, 1);
  
  digitizer->AddInputFile("/volatile/halla/sbs/efuchey/gmn13.5_elastic_sig_20190725_15/elastic_0.root", 0, 1);
  
  int nmiss = 0;
  for(int i = 0; i<nentries*nbkgd+nmiss; i++){
    // TSBSGeant4File *f_b = new TSBSGeant4File(Form("/volatile/halla/sbs/efuchey/gmn13.5_beam_bkgd_blsh_20190724_01/beam_bkgd_%d.root", i));
    // if(debuglevel>=2)cout << "Add to digitizer file " << f_b->GetName() << endl;
    // if(!f_b->Open()){
    //   continue;
    //   nmiss++;
    // }
    //f_b->SetSource(1);
    digitizer->AddInputFile(Form("/volatile/halla/sbs/efuchey/gmn13.5_beam_bkgd_blsh_20190724_01/beam_bkgd_%d.root", i), 1, -nbkgd);
  }
  
  // It is recommended  to declare the detector with its unique ID (second parameter)
  // See list of unique det IDs defined in src/g4sbs_types.h
  
  if(debuglevel>=1)cout << "Declare detectors and add them to digitizer " << endl;
  
  TSBSSimHCal *hcal = new TSBSSimHCal("hcal", 0);
  digitizer->AddDetector(hcal);
  
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
  
  if(debuglevel>=1)cout << "About to process digitization for " << nentries << "events " << endl;
  
  //digitizer->Process(f, nentries);
  digitizer->Process(nentries);
  
  //cout << "delete detectors" << endl;
  //delete hodo;
  //delete cdet;
  //delete grinch;
}
