#include <iostream>

#include "TSBSCher.h"
#include "THaEvData.h"
#include "THaApparatus.h"
#include "TMath.h"
#include "ha_compiledata.h"

using namespace std;

//Recommanded constructor
TSBSCher::TSBSCher( const char *name, const char *desc )
  : TSBSDet (name, desc)
{
  fNPMTrows = 0;
  fNPMTcolsMax = 0;
  
  return;
}

TSBSCher::~TSBSCher()
{
}


const char* TSBSCher::GetDBFileName() const {
    THaApparatus *app = GetApparatus();
    if( app )
      return Form ("%s.", app->GetName());
    else
      return fPrefix;
}

Int_t
TSBSCher::ReadGeometry (FILE* file, const TDatime& date,
			      Bool_t required)
{

  Int_t err = THaDetector::ReadGeometry (file, date, required);
  if (err)
    return err;
  
  const DBRequest request[] =
    {
      //{"zckov_in",     &fZCkovIn,      kDouble, 0, 1},
      {"zckov_in",     &fZPos,         kDouble, 0, 1},
      {"n_radiator",   &fNradiator,    kDouble, 0, 1},
      {"l_radiator",   &fLradiator,    kDouble, 0, 1},
      //{"npmts",        &fNPMTs,    kInt,    0, 1},
      {"npmts",        &fNChannels,    kInt,    0, 1},
      {"npmtrows",     &fNPMTrows,     kInt,    0, 1},
      {"npmtcolsmax",  &fNPMTcolsMax,  kInt,    0, 1},
      {"pmtdistx",     &fPMTdistX,     kDouble, 0, 1},
      {"pmtdisty",     &fPMTdistY,     kDouble, 0, 1},
      {"x_tcpmt",      &fX_TCPMT,      kDouble, 0, 1},
      {"y_tcpmt",      &fY_TCPMT,      kDouble, 0, 1},
      {0}
    };
  err = LoadDB (file, date, request, fPrefix);

  if (err)
    return err;
  
  fPMTmatrixHext = (fNPMTcolsMax-1)*fPMTdistY;
  fPMTmatrixVext = (fNPMTrows-1)*fPMTdistX;
  
  return kOK;
}

//
Int_t
TSBSCher::Decode (const THaEvData& ed )
{
  return 0;
}

void
TSBSCher::Print ()
{
  //Print TSBSCher info
  cout << "I'm a Cherenkov detector named " << GetName() << ", located " << fZPos << " m away from the spectrometer bending point." << endl;
  cout << " Total number of PMTs = " << fNChannels << ", set in " << fNPMTrows << " rows of " 
       << fNPMTcolsMax << " or " << fNPMTcolsMax-1 << " PMTs each." << endl;
  cout << " Distance between 2 PMTs rows : " << fPMTdistX << " m" << endl
       << " Distance between 2 PMTs columns : " << fPMTdistY << " m" << endl
       << " horizontal/vertical extension of PMT matrix (m) : " 
       << fPMTmatrixHext << "/" << fPMTmatrixVext << endl;
  cout << "Position of top close PMT in PMT matrix (in transport coord): " << endl
       << "X = " << fX_TCPMT << endl
       << "Y = " << fY_TCPMT << endl;
  
}
