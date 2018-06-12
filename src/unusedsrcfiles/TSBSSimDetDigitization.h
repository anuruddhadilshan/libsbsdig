
// Class handling digitization of GEMs

#ifndef __TSBSSimDetDigitization__
#define __TSBSSimDetDigitization__

#include "TRandom3.h"
#include "TVector3.h"
#include "TArrayI.h"
#include "TArrayD.h"

#include "THaAnalysisObject.h"

#include <vector>
#include <map>

class TFile;
class TTree;

class TSBSDetData;
class TSBSSpec;
class TSBSSimEvent;
class TSBSGeant4File;


class TSBSSimDetDigitization: public THaAnalysisObject
{
 public:
  //Constructor and destructor
  TSBSSimDetDigitization( const TSBSSpec& spect,
			  const char* name = "ratedig");
  virtual ~TSBSSimDetDigitization();
  
  //full initialization of all parameters with database
  void Initialize(const TSBSSpec& spect);
  
  Double_t GetGateWidth(){ return fGateWidth; }
  Double_t GetTriggerOffset(){ return fTriggerOffset; }
  Double_t GetTriggerJitter(){ return fTriggerJitter; }
  Double_t GetPulseNoiseSigma(){ return fPulseNoiseSigma; }
  
  //bool GetTDCtimes(int detnum, double C, double t0, double& t1, double& t2);
  
  // Tree methods
  // To write a tree with digitization results:
  //   Call InitTree before main loop;
  //   Call SetTreeEvent in main loop (before or after Digitize)
  //   Call FillTree in main loop (after Digitize and SetTreeEvent)
  // Call WriteTree and CloseTree after main loop

  void InitTree (const TSBSSpec& spect, const TString& ofile);
  
  void FillTree ();
  void WriteTree () const;
  void CloseTree () const;

  TSBSSimEvent* GetEvent() const { return fEvent; }

 protected:
  //Everything which goes here MUST be overloaded
  void Print() const;// print info

  Int_t ReadDatabase (const TDatime& date);
  
  //This is in those three functions that the job is done, more specifically in AddititveDigitize
  Int_t Digitize (const TSBSDetData& gdata, const TSBSSpec& spect); // digitize event
  Int_t AdditiveDigitize (const TSBSDetData& gdata, const TSBSSpec& spect); // add more digitized data, e.g. background
  void NoDigitize (const TSBSDetData& gdata, const TSBSSpec& spect); // do not digitize event, just fill tree

  void SetTreeEvent (const TSBSDetData& tsgd,
		     const TSBSGeant4File& f,
		     Int_t evnum = -1);
  Short_t SetTreeHit (const UInt_t ih,
		      const TSBSSpec& spect,
		      //TSolGEMVStrip* const *dh,
		      const TSBSDetData& tsgd,
		      Double_t t0 ); // called from Digitization
  
  TFile* fOFile;          // Output ROOT file
  TTree* fOTree;          // Output tree
  TSBSSimEvent* fEvent;   // Output event structure, written to tree

  void MakePrefix() { THaAnalysisObject::MakePrefix(0); }
  void DeleteObjects();

  Double_t fGateWidth;
  Double_t fTriggerOffset;
  Double_t fTriggerJitter;
  Double_t fPulseNoiseSigma;
  
  TRandom3 fTrnd;     // time randomizer
  
 private:


  ClassDef (TSBSSimDetDigitization, 0)
};

#endif

