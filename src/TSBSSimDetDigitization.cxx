#include "TSBSSimDetDigitization.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TClonesArray.h"

#include "TSBSGeant4File.h"  // needed for g4sbsgendata class def
#include "TSBSDetData.h"
#include "TSBSSpec.h"
#include "TSBSSimEvent.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <utility>

#define q_e 1.602176487e-19

using namespace std;

static TSBSDBManager* manager = TSBSDBManager::GetInstance();

TSBSSimDetDigitization::TSBSSimDetDigitization( const TSBSSpec& spect,
						const char* name)
  : THaAnalysisObject(name, "Generic simulation digitizer"),
    // fDoMapSector(false), fSignalSector(0), fDP(0), fdh(0), fNChambers(0), fNPlanes(0),
    fOFile(0), fOTree(0), fEvent(0)
{
  Init();
  Initialize (spect);
  
  fEvent = new TSBSSimEvent(5);
}

TSBSSimDetDigitization::~TSBSSimDetDigitization()
{
  DeleteObjects();
}

void TSBSSimDetDigitization::DeleteObjects()
{
  delete fOFile;      fOFile = 0;
  delete fOTree;      fOTree = 0;
  delete fEvent;      fEvent = 0;
  
  /*
  fTDCArrays.clear();
  fTDCtimeArrays.clear();
  fPulseHeight.clear();
  fPulseTimes.clear();
  */
}

void
TSBSSimDetDigitization::Initialize(const TSBSSpec& spect)
{
  // Initialize digitization structures based on parameters from given
  // spectrometer

  // Avoid memory leaks in case of reinitialization
  DeleteObjects();
  
  /*
  fNDetectors = manager->GetNDetectors();
  fNPMTs.resize(fNDetectors);
  for(int i = 0; i<fNDetectors; i++){
    fNPMTs[i] = manager->GetNPMTs(i);
    TArrayI TDCarray1(fNPMTs[i]);
    TArrayI TDCarray2(fNPMTs[i]);
    fTDCArrays.push_back(std::make_pair(TDCarray1, TDCarray2));
    TArrayD TDCtimearray1(fNPMTs[i]);
    TArrayD TDCtimearray2(fNPMTs[i]);
    fTDCtimeArrays.push_back(std::make_pair(TDCtimearray1, TDCtimearray2));
  }
  
  if(strcmp(fPMTtimesDataFileNames.Data(),"")!=0)ReadPMTtimesDataFiles();
  */
}

Int_t
TSBSSimDetDigitization::ReadDatabase (const TDatime& date)
{
  /*
  FILE* file = OpenFile (date);
  if (!file) return kFileError;

  const DBRequest request[] =
    {
      { "pmtgain",                   &fPMTGain,                   kDouble },
      //{ "pmtpulseshapetau",          &fPMTPulseShapeTau,          kDouble },
      { "pmttransittime",            &fPMTTransitTime,            kDouble },
      // { "pmtrisetime",               &fPMTRiseTime,               kDouble },
      // { "pmtjitter",                 &fPMTJitter,                 kDouble },
      // { "pmtfwhm",                   &fPMTFWHM,                   kDouble },
      { "pmttimesdatafilenames",     &fPMTtimesDataFileNames,     kTString, 0, 1},
      { "tdctimeconversion",         &fTDCTimeConversion,         kDouble },
      { "tdcbits",                   &fTDCbits,                   kInt    },
      { "tdcthreshold",              &fTDCthreshold,              kDouble },
      { "tdcresolution",             &fTDCresolution,             kDouble },
      { "triggeroffset",             &fTriggerOffset,             kDouble },
      { "triggerjitter",             &fTriggerJitter,             kDouble },
      { "gatewidth",                 &fGateWidth,                 kDouble },
      { "pulsenoiseconst",           &fPulseNoiseConst,           kDouble },
      { "pulsenoisesigma",           &fPulseNoiseSigma,           kDouble },
      { "do_crosstalk",              &fDoCrossTalk,               kInt    },
      { "crosstalk_mean",            &fCrossTalkMean,             kDouble },
      { "crosstalk_sigma",           &fCrossTalkSigma,            kDouble },
      { "crosstalk_chan_apart",      &fCrossTalkChanApart,        kInt    },
      { "readoutimpedance",          &fReadOutImpedance,          kDouble },
      { 0 }
    };
  
  Int_t err = LoadDB (file, date, request, fPrefix);
  fclose(file);
  if (err)
    return kInitError; 
  
  // fPMTTransitTime*= 1.e-9;
  // fPMTRiseTime*= 1.e-9;
  // fPMTJitter*= 1.e-9;
  // fPMTFWHM*= 1.e-9;
  // fTDCresolution*= 1.e-9;
  */
  
  /*
  if( fEleSamplingPoints < 0 || fEleSamplingPoints > 10 )
    fEleSamplingPoints = 10;
  if( fADCbits < 1 || fADCbits > MAX_ADCBITS ) {
    Error("ReadDatabase", "Invalid parameter adcbits = %d", fADCbits );
    return kInitError;
  }
  fAvalancheFiducialBand = TMath::Abs(fAvalancheFiducialBand);
  */
  
  return kOK;
}


Int_t
TSBSSimDetDigitization::Digitize (const TSBSDetData& gdata, const TSBSSpec& spect)
{
  // Digitize event after clearing all previous digitization results.

  fEvent->Clear();
  
  return AdditiveDigitize( gdata, spect );
}

Int_t
TSBSSimDetDigitization::AdditiveDigitize (const TSBSDetData& chdata, const TSBSSpec& spect)
{

  return 0;
}

void
TSBSSimDetDigitization::NoDigitize (const TSBSDetData& chdata, const TSBSSpec& spect) // do not digitize event, just fill the tree
{
  //  if (!fEvCleared)  //?
  fEvent->Clear();
}

//___________________________________________________________________________________
void
TSBSSimDetDigitization::Print() const
{
}

// Tree methods
void
TSBSSimDetDigitization::InitTree (const TSBSSpec& spect, const TString& ofile)
{
  fOFile = new TFile( ofile, "RECREATE");

  if (fOFile == 0 || fOFile->IsZombie() )
    {
      cerr << "Error: cannot open output file " << ofile << endl;
      delete fOFile; fOFile = 0;
      return;
    }

  fOTree = new TTree( treeName, "Tree of digitized values");

  // create the tree variables

  fOTree->Branch( eventBranchName, "TSBSSimEvent", &fEvent );

}

void
TSBSSimDetDigitization::SetTreeEvent (const TSBSDetData& tscd,
				       const TSBSGeant4File& f, Int_t evnum )
{
  // Set overall event info.
  fEvent->Clear("all");
  fEvent->fRunID = tscd.GetRun();
  // FIXME: still makes sense if background added?
  if(tscd.GetSource()==0){
    fEvent->fEvtID = (evnum < 0) ? tscd.GetEvent() : evnum;
  }
  for( UInt_t i=0; i<f.GetNGen(); ++i ) {
    const g4sbsgendata* gd = f.GetGenData(i);
    fEvent->AddTrack( gd->GetTRID(), gd->GetPID(),
     		      gd->GetV(), // Vertex coordinates in [m]
     		      gd->GetP()  // Momentum in [GeV]
     		      );
  }
  // FIXED: one GenData per event: signal, primary particle
  if( f.GetNGen() > 0 )
    fEvent->fWeight = f.GetGenData(0)->GetWeight();
  
}

Short_t
TSBSSimDetDigitization::SetTreeHit (const UInt_t ih,
				     const TSBSSpec& spect,
				     //TSolGEMVStrip* const *dh,
				     const TSBSDetData& tscd,
				     Double_t t0 )
{
  return 0;
}


void
TSBSSimDetDigitization::FillTree ()
{
  //cout << "Fill tree " << fOFile << " " << fOTree << endl;
  //fOFile = fOTree->GetCurrentFile();//CHECK ?
  if (fOFile && fOTree
      // added this line to not write events where there are no entries

      // Remove for background study
      //      && fEvent->fGEMStrips.size() > 0 && fEvent->fGEMHit.size() > 0

      )
    {
      //CleanClusterList();
      fOFile->cd();
      //fEvent->Print("all");
      fOTree->Fill();
    }
}

void
TSBSSimDetDigitization::WriteTree () const
{
  //cout << "write tree " << fOFile << " " << fOTree << endl;
  
  if (fOFile && fOTree) {
    fOFile->cd();
    fOTree->Write();
  }
}

void
TSBSSimDetDigitization::CloseTree () const
{
  if (fOFile) fOFile->Close();
}

