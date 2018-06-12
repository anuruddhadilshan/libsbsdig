// manage the input data (output of the GEANT4 simulation)
// for the Cherenkov (Cher) detectors

#ifndef __TSBSCherData__
#define __TSBSCherData__

#include <TRandom.h>
#include <TVector3.h>
#include "TSBSDetData.h"
#include <vector>

/*
class CherHitData : public TObject {
 public:
  
  CherHitData(){
    fDetID = 0; fPMTID = 0; fXPMT = 0.; fYPMT = 0.; fNpe = 0.; fTime = 0.; fTimeRMS = 0.; fType = 0; fPos_det = TVector3(0., 0., 0.); fPID_MCtrack = 0; fMom_MCtrack = TVector3(0., 0., 0.); fVtx_MCtrack = TVector3(0., 0., 0.); fOriginVolFlag = 0;
  }
  
  
  ~CherHitData(){}
  
  void SetDetID (UInt_t detid)             { fDetID = detid; }
  void SetPMTID (UInt_t pmtid)             { fPMTID = pmtid; }
  void SetXPMT (Double_t xpmt)             { fXPMT = xpmt; }
  void SetYPMT (Double_t ypmt)             { fYPMT = ypmt; }
  void SetPEyield (Double_t npe)           { fNpe = npe; }
  void SetTime(Double_t t)                 { fTime = t; }
  void SetTimeRMS(Double_t trms)           { fTimeRMS = trms; }
  void SetParticleType (UInt_t type)          { fType = type; }
  void SetPositionDet (const TVector3& X)     { fPos_det = X; }
  void SetMCtrackPID (Int_t pid)              { fPID_MCtrack = pid; }
  void SetMCtrackMomentum (const TVector3& p) { fMom_MCtrack = p; }
  void SetMCtrackVertex (const TVector3& X)   { fVtx_MCtrack = X; }
  void SetOriginVolFlag (UInt_t volflag)      { fOriginVolFlag = volflag; }
  
  UInt_t          GetDetID ()          const { return fDetID; }
  UInt_t          GetPMTID ()          const { return fPMTID; }
  Double_t        GetXPMT ()           const { return fXPMT; }
  Double_t        GetYPMT ()           const { return fYPMT; }
  Double_t        GetPEyield ()        const { return fNpe; }
  Double_t        GetTime ()           const { return fTime; }
  Double_t        GetTimeRMS ()        const { return fTimeRMS; }
  UInt_t          GetParticleType ()      const { return fType; }
  const TVector3& GetPositionDet ()       const { return fPos_det; }
  UInt_t          GetMCtrackPID ()        const { return fPID_MCtrack; }
  const TVector3& GetMCtrackMomentum ()   const { return fMom_MCtrack; }
  const TVector3& GetMCtrackVertex ()     const { return fVtx_MCtrack; }
  UInt_t          GetOriginVolFlag ()     const { return fOriginVolFlag; }
    
 private:
  UInt_t    fDetID;          // Hit detector ID
  UInt_t    fPMTID;          // Hit PMT
  Double_t  fXPMT;           // X coordinate of the PMT in transport coordinates
  Double_t  fYPMT;           // Y coordinate of the PMT in transport coordinates
  Double_t  fNpe;            // Number of photoelectrons
  Double_t  fTime;           // hit time
  Double_t  fTimeRMS;        // hit time RMS: should be small (<0.1 ns)
  UInt_t    fType;           // particle type: 1 = primary, >1 secondaries
  TVector3  fPos_det;        // position of the optical photon at detection
  Int_t     fPID_MCtrack;    // Particle ID of the MC particle which produced the hit (*)
  TVector3  fMom_MCtrack;    // momentum of the MC particle which produced the hit (*)
  TVector3  fVtx_MCtrack;    // vertex of the MC particle which produced the hit (*)
  UInt_t    fOriginVolFlag;  // origin volume of the photon
  
  ClassDef(CherHitData, 0);
};
*/

class TSBSCherData : public TSBSDetData {
 public:
  
  TSBSCherData (UInt_t h = 0);
  virtual ~TSBSCherData();

  void ClearEvent();
  void InitEvent (UInt_t h);

  // Positions are in mm
  // Time is in ns
  // Momenta are in MeV

  /*  
  void SetNHit (UInt_t h)   { fHitData->Expand(h); }
  void SetHitDetID (UInt_t k, UInt_t detid)             { ((CherHitData*)fHitData->At(k))->SetDetID(detid); }
  void SetHitPMTID (UInt_t k, UInt_t pmtid)             { ((CherHitData*)fHitData->At(k))->SetPMTID (pmtid); }
  void SetHitXPMT (UInt_t k, Double_t xpmt)             { ((CherHitData*)fHitData->At(k))->SetXPMT (xpmt); }
  void SetHitYPMT (UInt_t k, Double_t ypmt)             { ((CherHitData*)fHitData->At(k))->SetYPMT (ypmt); }
  void SetHitPEyield (UInt_t k, Double_t npe)           { ((CherHitData*)fHitData->At(k))->SetPEyield(npe); }
  void SetHitTime(UInt_t k, Double_t t)                 { ((CherHitData*)fHitData->At(k))->SetTime(t); }
  void SetHitTimeRMS(UInt_t k, Double_t trms)           { ((CherHitData*)fHitData->At(k))->SetTimeRMS(trms); }
  void SetParticleType (UInt_t k, UInt_t type)          { ((CherHitData*)fHitData->At(k))->SetParticleType(type); }
  void SetPositionDet (UInt_t k, const TVector3& X)     { ((CherHitData*)fHitData->At(k))->SetPositionDet(X); }
  void SetMCtrackPID (UInt_t k, Int_t pid)              { ((CherHitData*)fHitData->At(k))->SetMCtrackPID(pid); }
  void SetMCtrackMomentum (UInt_t k, const TVector3& p) { ((CherHitData*)fHitData->At(k))->SetMCtrackMomentum(p); }
  void SetMCtrackVertex (UInt_t k, const TVector3& X)   { ((CherHitData*)fHitData->At(k))->SetMCtrackVertex(X); }
  void SetOriginVolFlag (UInt_t k, UInt_t volflag)      { ((CherHitData*)fHitData->At(k))->SetOriginVolFlag(volflag); }
  
  UInt_t GetNHit()   const { return fHitData->GetSize(); }
  UInt_t          GetHitDetID (UInt_t k)        const { return ((CherHitData*)fHitData->At(k))->GetDetID(); }
  UInt_t          GetHitPMTID (UInt_t k)        const { return ((CherHitData*)fHitData->At(k))->GetPMTID(); }
  Double_t        GetHitXPMT (UInt_t k)         const { return ((CherHitData*)fHitData->At(k))->GetXPMT(); }
  Double_t        GetHitYPMT (UInt_t k)         const { return ((CherHitData*)fHitData->At(k))->GetYPMT(); }
  Double_t        GetHitPEyield (UInt_t k)      const { return ((CherHitData*)fHitData->At(k))->GetPEyield(); }
  Double_t        GetHitTime (UInt_t k)         const { return ((CherHitData*)fHitData->At(k))->GetTime(); }
  Double_t        GetHitTimeRMS (UInt_t k)      const { return ((CherHitData*)fHitData->At(k))->GetTimeRMS(); }
  UInt_t          GetParticleType (UInt_t k)    const { return ((CherHitData*)fHitData->At(k))->GetParticleType(); }
  const TVector3& GetPositionDet (UInt_t k)     const { return ((CherHitData*)fHitData->At(k))->GetPositionDet(); }
  UInt_t          GetMCtrackPID (Int_t k)       const { return ((CherHitData*)fHitData->At(k))->GetMCtrackPID(); }
  const TVector3& GetMCtrackMomentum (UInt_t k) const { return ((CherHitData*)fHitData->At(k))->GetMCtrackMomentum(); }
  const TVector3& GetMCtrackVertex (UInt_t k)   const { return ((CherHitData*)fHitData->At(k))->GetMCtrackVertex(); }
  UInt_t          GetOriginVolFlag (UInt_t k)   const { return ((CherHitData*)fHitData->At(k))->GetOriginVolFlag(); }
  */
    
  void SetNHit (UInt_t h)   { fHitData.resize(h); }
  void SetHitDetID (UInt_t k, UInt_t detid)             { fHitData[k].fDetID  = detid; }
  void SetHitPMTID (UInt_t k, UInt_t pmtid)             { fHitData[k].fPMTID  = pmtid; }
  void SetHitXPMT (UInt_t k, Double_t xpmt)             { fHitData[k].fXPMT  = xpmt; }
  void SetHitYPMT (UInt_t k, Double_t ypmt)             { fHitData[k].fYPMT  = ypmt; }
  void SetHitPEyield (UInt_t k, Double_t npe)           { fHitData[k].fNpe = npe; }
  void SetHitTime(UInt_t k, Double_t t)                 { fHitData[k].fTime = t; }
  void SetHitTimeRMS(UInt_t k, Double_t trms)           { fHitData[k].fTimeRMS = trms; }
  void SetParticleType (UInt_t k, UInt_t type)          { fHitData[k].fType = type; }
  void SetPositionDet (UInt_t k, const TVector3& X)     { fHitData[k].fPos_det = X; }
  void SetMCtrackPID (UInt_t k, Int_t pid)              { fHitData[k].fPID_MCtrack = pid; }
  void SetMCtrackMomentum (UInt_t k, const TVector3& p) { fHitData[k].fMom_MCtrack = p; }
  void SetMCtrackVertex (UInt_t k, const TVector3& X)   { fHitData[k].fVtx_MCtrack = X; }
  void SetOriginVolFlag (UInt_t k, UInt_t volflag)      { fHitData[k].fOriginVolFlag = volflag; }
  
  UInt_t GetNHit()   const { return fHitData.size(); }
  UInt_t          GetHitDetID (UInt_t k)        const { return fHitData[k].fDetID; }
  UInt_t          GetHitPMTID (UInt_t k)        const { return fHitData[k].fPMTID; }
  Double_t        GetHitXPMT (UInt_t k)         const { return fHitData[k].fXPMT; }
  Double_t        GetHitYPMT (UInt_t k)         const { return fHitData[k].fYPMT; }
  Double_t        GetHitPEyield (UInt_t k)      const { return fHitData[k].fNpe; }
  Double_t        GetHitTime (UInt_t k)         const { return fHitData[k].fTime; }
  Double_t        GetHitTimeRMS (UInt_t k)      const { return fHitData[k].fTimeRMS; }
  UInt_t          GetParticleType (UInt_t k)    const { return fHitData[k].fType; }
  const TVector3& GetPositionDet (UInt_t k)     const { return fHitData[k].fPos_det; }
  UInt_t          GetMCtrackPID (Int_t k)       const { return fHitData[k].fPID_MCtrack; }
  const TVector3& GetMCtrackMomentum (UInt_t k) const { return fHitData[k].fMom_MCtrack; }
  const TVector3& GetMCtrackVertex (UInt_t k)   const { return fHitData[k].fVtx_MCtrack; }
  UInt_t GetOriginVolFlag (UInt_t k) const { return fHitData[k].fOriginVolFlag; }
  
  // Hit data
  // moved in "public" to allow it to compile with Root6/CentOS7
  struct CherHitData_t {
    UInt_t    fDetID;          // Hit detector ID
    UInt_t    fPMTID;          // Hit PMT
    Double_t  fXPMT;           // X coordinate of the PMT in transport coordinates
    Double_t  fYPMT;           // Y coordinate of the PMT in transport coordinates
    Double_t  fNpe;            // Number of photoelectrons
    Double_t  fTime;           // hit time
    Double_t  fTimeRMS;        // hit time RMS: should be small (<0.1 ns)
    UInt_t    fType;           // particle type: 1 = primary, >1 secondaries
    TVector3  fPos_det;        // position of the optical photon at detection
    Int_t     fPID_MCtrack;    // Particle ID of the MC particle which produced the hit (*)
    TVector3  fMom_MCtrack;    // momentum of the MC particle which produced the hit (*)
    TVector3  fVtx_MCtrack;    // vertex of the MC particle which produced the hit (*)
    UInt_t    fOriginVolFlag;  // origin volume of the photon
  };//  (*) valid if identified

  void PrintHit (UInt_t k) const;
  
 private:
  
  std::vector<CherHitData_t> fHitData;
  //TClonesArray* fHitData;
};
#endif
