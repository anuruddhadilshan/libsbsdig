#include "TSBSDBManager.h"
#include "TSBSSimDecoder.h"
#include <cassert>
#include <cmath>
#include "TMath.h"
#include "TVector2.h"
#include "TRandom3.h"
#include "TString.h"

TSBSDBManager * TSBSDBManager::fManager = NULL;

TSBSDBManager::TSBSDBManager() 
: fErrID(-999), fErrVal(-999.)
{
}
//______________________________________________________________
TSBSDBManager::~TSBSDBManager()
{
}

//______________________________________________________________
Int_t TSBSDBManager::LoadGeneralInfo(const string& fileName)
{  
  // Load the experiment/setup general info
  ifstream input(fileName.c_str());
  if (!input.is_open()){
    cout<<"cannot find general information file "<<fileName
	<<". Exiting the program"<<endl;
        exit(0);
  }
  
  const string prefix = "geninfo.";
  
  string exp_str;
  string specs_str;
  
  //first, load the experiment general info: expt type, number and names of spectrometers
  DBRequest request[] = {
    {"sbsexptype",         &exp_str,    kTString, 0, 1},
    {"nspectrometers",     &fNSpecs,    kInt,     0, 1},
    {"spectrometernames",  &specs_str,  kTString, 0, 1},
    { 0 }
  };

  int err = LoadDB( input, request,  prefix);
  
  if( err ) exit(2); 
  
  //assing the right exp_type value to the exp_type flag according to the expt name
  if(exp_str.compare("gmn")==0)fSBSExpType = kGMn;
  if(exp_str.compare("gep")==0)fSBSExpType = kGEp;
  if(exp_str.compare("gen")==0)fSBSExpType = kGEn;
  if(exp_str.compare("sidis")==0)fSBSExpType = kSIDIS;
  if(exp_str.compare("a1n")==0)fSBSExpType = kA1n;
  if(exp_str.compare("tdis")==0)fSBSExpType = kTDIS;
  if(exp_str.compare("ndvcs")==0)fSBSExpType = kDVCS;
  
  //split the full string to extract the individual spectrometer names
  fSpecNames = vsplit(specs_str);
  
  //Then, loop on the spectrometers to gather the detector number and names, and the MC signal of interest
  for(int i_spec = 0; i_spec<fNSpecs; i_spec++){
    SpectroInfo specinfo;
    int ndets;
    string dets_str;
    int nsig;

    string prefix2 = prefix+fSpecNames.at(i_spec)+".";
    
    std::vector<int>* pid = 0;
    std::vector<int>* tid = 0;
    
    try{
      pid = new vector<int>;
      tid = new vector<int>;
      
      DBRequest request[] = {
	{"nsignal",        &nsig,      kInt,      0, 1},
	{"signal.pid",     pid,        kIntV,     0, 1},
	{"signal.tid",     tid,        kIntV,     0, 1},
	{"ndetectors",     &ndets,     kInt,      0, 1},
	{"detectornames",  &dets_str,  kTString,  0, 1},
	{ 0 }
      };
      
      Int_t err = LoadDB (input, request, prefix);
      //input.close();
      if (err){
	input.close();
	return kInitError;
      }
      
      specinfo.fNDets = ndets;
      specinfo.fDetNames = vsplit(dets_str);
      
      for(int i_sig = 0; i_sig<nsig; i_sig++){
	SignalInfo siginfo(pid->at(i_sig), tid->at(i_sig));
	specinfo.MCsignalInfo.push_back(siginfo);
      }
      fSpectroInfos.push_back(specinfo);
	
      delete pid;
      delete tid;
    }  catch(...) {
      delete pid;
      delete tid;
      input.close();
      throw;
    }//end try / catch
    
    
  }
  input.close();
  return(kOK);
}
/*
//______________________________________________________________
void TSBSDBManager::LoadGeneralInfo(const string& fileName)
{  
  ifstream input(fileName.c_str());
  if (!input.is_open()){
    cout<<"cannot find general information file "<<fileName
	<<". Exiting the program"<<endl;
        exit(0);
  }
  const string prefix = "generalinfo.";
  DBRequest request[] = {
    {"g4sbs_exptype",  &fg4sbsExpType , kInt,    0, 1},
    {"g4sbs_dettype",  &fg4sbsDetType , kInt,    0, 1},
    {"ndetectors",     &fNDetectors   , kInt,    0, 1},
    {"chan_per_slot",  &fChanPerSlot  , kInt,    0, 1},
    {"slot_per_crate", &fSlotPerCrate , kInt,    0, 1},
    {"nsignal",        &fNSigParticle , kInt,    0, 1},
    { 0 }
  };
  int pid, tid;
  DBRequest signalRequest[] = {
    {"pid",                 &pid,                   kInt, 0, 1},
    {"tid",                 &tid,                   kInt, 0, 1},
    { 0 }
  };
  int err = LoadDB( input, request,  prefix);
  
  if( err ) exit(2); 
  
  for (int i=0; i<fNSigParticle; i++){
    ostringstream signal_prefix(prefix, ios_base::ate);
    signal_prefix<<"signal"<<i+1<<".";
    
    err = LoadDB(input, signalRequest, signal_prefix.str());
    
    fSigPID.push_back(pid);
    fSigTID.push_back(tid);
    
    if( err ) exit(2); 
  }
  
  for (int i=0; i<fNDetectors; i++){
    vector<GeoInfo> thisInfo;
    thisInfo.clear();
    fGeoInfo = thisInfo;
  }
}


//______________________________________________________________
void TSBSDBManager::LoadGeoInfo(const string& prefix)
{
  // Include DB_DIR (standard Hall A analyzer DB path in the search)
  std::string path = "";
  if(std::getenv("DB_DIR")) {
    path = std::string(std::getenv("DB_DIR")) + "/";
  }
  const string& fileName = path+"db_"+prefix+".dat";
    
  ifstream input(fileName.c_str());
  if (!input.is_open()){
    cout<<"cannot find geometry file "<<fileName
	<<". Exiting the program"<<endl;
    exit(0);
  }
  
  GeoInfo thisGeo;
  
  DBRequest request[] = {
    {"zckov_in",     &thisGeo.fZCkovIn,      kDouble, 0, 1},
    {"n_radiator",   &thisGeo.fNradiator,    kDouble, 0, 1},
    {"l_radiator",   &thisGeo.fLradiator,    kDouble, 0, 1},
    {"npmts",        &thisGeo.fNPMTs,        kInt,    0, 1},
    {"npmtrows",     &thisGeo.fNPMTrows,     kInt,    0, 1},
    {"npmtcolsmax",  &thisGeo.fNPMTcolsMax,  kInt,    0, 1},
    {"pmtdistx",     &thisGeo.fPMTdistX,     kDouble, 0, 1},
    {"pmtdisty",     &thisGeo.fPMTdistY,     kDouble, 0, 1},
    {"x_tcpmt",      &thisGeo.fX_TCPMT,      kDouble, 0, 1},
    {"y_tcpmt",      &thisGeo.fY_TCPMT,      kDouble, 0, 1},
    { 0 }
  };
  
  for (int i=0; i<fNDetectors; i++){
    //map<int, vector<GeoInfo> >::iterator it = fGeoInfo.find(i);
    //if (it == fGeoInfo.end()) { cout<<"unexpected detector id "<<i<<endl; }
    ostringstream detector_prefix(prefix, ios_base::ate);
    detector_prefix<<".cher"<<i<<".";
    
    int err = LoadDB(input, request,detector_prefix.str());
    if( err ) exit(2);
    
    thisGeo.fPMTmatrixHext = (thisGeo.fNPMTcolsMax-1)*thisGeo.fPMTdistY;
    thisGeo.fPMTmatrixVext = (thisGeo.fNPMTrows-1)*thisGeo.fPMTdistX;
    
    fGeoInfo.push_back(thisGeo);
  }
  
  
}

//______________________________________________________________
string TSBSDBManager::FindKey( ifstream& inp, const string& key )
{
  static const string empty("");
  string line;
  string::size_type keylen = key.size();
  inp.seekg(0); // could probably be more efficient, but it's fast enough
  while( getline(inp,line) ) {
    if( line.size() <= keylen )
      continue;
    if( line.compare(0,keylen,key) == 0 ) {
      if( keylen < line.size() ) {
	string::size_type pos = line.find_first_not_of(" \t=", keylen);
	if( pos != string::npos )
	  return line.substr(pos);
      }
      break;
    }
  }
  return empty;
}
//_________________________________________________________________
int TSBSDBManager::LoadDB( ifstream& inp, DBRequest* request, const string& prefix )
{
  DBRequest* item = request;
  while( item->name ) {
    ostringstream sn(prefix, ios_base::ate);
    sn << item->name;
    const string& key = sn.str();
    string val = FindKey(inp,key);
    if( !val.empty() ) {
      istringstream sv(val);
      switch(item->type){
        case kDouble:
          sv >> *((double*)item->var);
          break;
        case kInt:
          sv >> *((Int_t*)item->var);
          break;
        default:
          return 1;
        break;
      }
      if( !sv ) {
	cerr << "Error converting key/value = " << key << "/" << val << endl;
	return 1;
      }
    } else {
      cerr << "key \"" << key << "\" not found" << endl;
      return 2;
    }
    ++item;
  }
  return 0;
}
//_____________________________________________________________________
const int & TSBSDBManager::GetSigPID(unsigned int i)
{
    if ( i >= fSigPID.size() ){ 
        cout<<"only "<<fSigPID.size()<<" signal particle registered"<<endl;
        return fErrID;
    }
    return fSigPID[i];
}
//______________________________________________________________________
const int & TSBSDBManager::GetSigTID(unsigned int i)
{
    if ( i >= fSigPID.size() ){ 
        cout<<"only "<<fSigPID.size()<<" signal particle registered"<<endl;
        return fErrID;
    }
    return fSigTID[i];
}

//_________________________________________________________________________
bool TSBSDBManager::CheckIndex(int i)
{
    if (i >= fNDetectors || i < 0){
        cout<<"invalid Detector ID requested: "<<i<<endl;
        return false;
    }

    return true;
}

//_________________________________________________________________________
const double & TSBSDBManager::GetZCkovIn(int i)
{
  if (!CheckIndex(i)) return fErrVal;
  return fGeoInfo.at(i).fZCkovIn;
}

//_________________________________________________________________________
const double & TSBSDBManager::GetNradiator(int i)
{
  if (!CheckIndex(i)) return fErrVal;
  return fGeoInfo.at(i).fNradiator;
}

//_________________________________________________________________________
const double & TSBSDBManager::GetLradiator(int i)
{
  if (!CheckIndex(i)) return fErrVal;
  return fGeoInfo.at(i).fLradiator;
}

//_________________________________________________________________________
const int & TSBSDBManager::GetNPMTs(int i)
{
  if (!CheckIndex(i)) return fErrID;
  return fGeoInfo.at(i).fNPMTs;
}

//_________________________________________________________________________
const int & TSBSDBManager::GetNPMTrows(int i)
{
  if (!CheckIndex(i)) return fErrID;
  return fGeoInfo.at(i).fNPMTrows;
}

//_________________________________________________________________________
const int & TSBSDBManager::GetNPMTcolsMax(int i)
{
  if (!CheckIndex(i)) return fErrID;
  return fGeoInfo.at(i).fNPMTcolsMax;
}

//_________________________________________________________________________
const double & TSBSDBManager::GetPMTmatrixHext(int i)
{
  if (!CheckIndex(i)) return fErrVal;
  return fGeoInfo.at(i).fPMTmatrixHext;
}

//_________________________________________________________________________
const double & TSBSDBManager::GetPMTmatrixVext(int i)
{
  if (!CheckIndex(i)) return fErrVal;
  return fGeoInfo.at(i).fPMTmatrixVext;
}

//_________________________________________________________________________
const double & TSBSDBManager::GetPMTdistX(int i)
{
  if (!CheckIndex(i)) return fErrVal;
  return fGeoInfo.at(i).fPMTdistX;
}

//_________________________________________________________________________
const double & TSBSDBManager::GetPMTdistY(int i)
{
  if (!CheckIndex(i)) return fErrVal;
  return fGeoInfo.at(i).fPMTdistY;
}

//_________________________________________________________________________
const double & TSBSDBManager::GetX_TCPMTs(int i)
{
  if (!CheckIndex(i)) return fErrVal;
  return fGeoInfo.at(i).fX_TCPMT;
}

//_________________________________________________________________________
const double & TSBSDBManager::GetY_TCPMTs(int i)
{
  if (!CheckIndex(i)) return fErrVal;
  return fGeoInfo.at(i).fY_TCPMT;
}
*/
