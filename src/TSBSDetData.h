// manage the input data (output of the GEANT4 simulation)
// for the gem part

#ifndef __TSBSDetData__
#define __TSBSDetData__

#include <TRandom.h>
#include <TClonesArray.h>
#include <TVector3.h>

class TSBSDetData 
{
 public:

  TSBSDetData ();
  virtual ~TSBSDetData();

  void ClearEvent();
  
  //void SetNHit (UInt_t h)   { fHitData->Expand(h); }
  void SetEvent (UInt_t id) { fEvtID = id; }
  void SetRun (UInt_t r)    { fRunID = r; }
  void SetSource (Int_t s)  { fSource = s; }
  
  //UInt_t GetNHit()   const { return fHitData->GetSize(); }
  UInt_t GetEvent()  const { return fEvtID; }
  UInt_t GetRun()    const { return fRunID; }
  Int_t  GetSource() const { return fSource; }
  
  void Print() const;

 protected:
  //TClonesArray* fHitData;
  void InitEvent (UInt_t h);
  void PrintHit (UInt_t k) const;// This class must be overloaded
  
 private:

  UInt_t fRunID, fEvtID;
  Int_t  fSource; // MC source file ID (0 = signal, >0 background)

};
#endif
