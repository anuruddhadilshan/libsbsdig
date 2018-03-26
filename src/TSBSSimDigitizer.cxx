#include "TSBSSimDigitizer.h"
#include "TSBSGeant4File.h"
#include "TSBSSimDetector.h"
#include "TSBSSimEvent.h"
#include <TTree.h>
#include <TFile.h>

TSBSSimDigitizer::TSBSSimDigitizer()
{
  fEvent = new TSBSSimEvent(5);
  fOutFile = new TFile("rootfiles/simout_test.root","RECREATE");
  //fOutTree = new TTree("TSBSDigi","");
  fOutTree = new TTree("digtree","");
  //fOutTree->Branch("fEvents",&fEvent);
  fOutTree->Branch("event",&fEvent);
}


TSBSSimDigitizer::~TSBSSimDigitizer()
{
}

int TSBSSimDigitizer::Process(TSBSGeant4File *f, int max_events)
{
  if(!f)
    return 0;

  // Can we open the file?
  f->SetSource(0);
  int res = f->Open();
  if( res != 1) {
    std::cerr << "Failed to open g4sbs rootfile. Failed with error code: "
      << res << std::endl;
    return 0;
  }

  if ( max_events <= 0 || max_events > f->GetEntries() )
    max_events = f->GetEntries();

  // Now loop through the file and digitize entries
  int d_flag_readevent = 0;
  int nevent = 0;
  int ngood = 0;
  bool has_data;
  while( f->ReadNextEvent(d_flag_readevent) && nevent<max_events ) {
    fEvent->Clear();
    has_data = false;

    // Loop through all detectors and have them parse data vector
    for(size_t det = 0; det < fDetectors.size(); det++) {
      fDetectors[det]->LoadEventData(f->GetDataVector());
    }
    // Now digitize all detectors
    for(size_t det = 0; det < fDetectors.size(); det++) {
      fDetectors[det]->Digitize(*fEvent);
    }
    fEvent->fEvtID = nevent;

    // Fill in tree if any of the detectors have data to fill
    for(size_t det = 0; det < fDetectors.size() && !has_data; det++) {
      has_data = fDetectors[det]->HasData();
    }
    if(has_data) {
      // Write to the tree
      fOutTree->Fill();
      std::cout << "Have data for event: " << nevent << std::endl;
    }

    nevent++;
  }

  fOutFile->Write();
  return 0;
}

int TSBSSimDigitizer::Add(TSBSSimDetector* detector)
{
  fDetectors.push_back(detector);
}
