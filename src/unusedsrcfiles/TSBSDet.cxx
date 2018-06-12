#include <iostream>

#include "TSBSDet.h"
#include "THaEvData.h"
#include "THaApparatus.h"
#include "TMath.h"
#include "ha_compiledata.h"

using namespace std;

//Recommanded constructor
TSBSDet::TSBSDet( const char *name, const char *desc )
  : THaDetector (name, desc)
{
  fNChannels = 0;
  fZPos = 0;
  
  return;
}

TSBSDet::~TSBSDet()
{
}


const char* TSBSDet::GetDBFileName() const {
    THaApparatus *app = GetApparatus();
    if( app )
      return Form ("%s.", app->GetName());
    else
      return fPrefix;
}
  

Int_t
TSBSDet::ReadDatabase (const TDatime& date)
{
  FILE* file = OpenFile (date);
  if (!file) return kFileError;

  Int_t err = ReadGeometry (file, date, false);

  fclose(file);
  if (err)
    return err;

  return kOK;
}

Int_t
TSBSDet::ReadGeometry (FILE* file, const TDatime& date,
			      Bool_t required)
{
  return kOK;
}


//
Int_t
TSBSDet::Decode (const THaEvData& ed )
{
  return 0;
}


void
TSBSDet::Print ()
{
}
