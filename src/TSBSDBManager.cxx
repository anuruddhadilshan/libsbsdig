#include "TSBSDBManager.h"
#include <cassert>
#include <cmath>
#include "TMath.h"
#include "TVector2.h"
#include "TRandom3.h"
#include "TString.h"
#include "TGEMSBSDBManager.h"

TSBSDBManager * TSBSDBManager::fManager = NULL;

TSBSDBManager::TSBSDBManager() 
  : fErrID(-999), fErrVal(-999.), fBkgdSpreadTimeWindowHW(0.)
{
  fRN = TRndmManager::GetInstance();
}
//______________________________________________________________
TSBSDBManager::~TSBSDBManager()
{
}

//______________________________________________________________
Int_t TSBSDBManager::LoadGenInfo(const string& fileName)
{  
  // Load the experiment/setup general info
  /*
  ifstream input(fileName.c_str());
  if (!input.is_open()){
    cout<<"cannot find general information file "<<fileName
	<<". Exiting the program"<<endl;
        exit(0);
  }
  */
  std::string path = "../db/";
  if(std::getenv("SBS_DIGI_DB")) {
    path = std::string(std::getenv("SBS_DIGI_DB")) + "/";
  }
  const string& PathfileName = path+fileName;
  
  cout << "File name " << fileName << endl;
  
  FILE* file = OpenFile( PathfileName.c_str(), GetInitDate() );
  
  const string prefix = "geninfo.";
  
  int rnseed = 0;
  
  string exp_str;
  string specs_str;
  
  //first, load the experiment general info: expt type, number and names of spectrometers
  DBRequest request[] = {
    {"sbsexptype", &exp_str,   kString, 0, 0},
    {"nspecs",     &fNSpecs,   kInt,    0, 0},
    {"specnames",  &specs_str, kString, 0, 0},
    {"randomseed", &rnseed,    kInt,    0, 1},
    { 0 }
  };
  
  fRN->SetSeed(rnseed);
  
  //int err = LoadDB( input, request,  prefix);
  int err = LoadDB(file, GetInitDate(), request,  prefix.c_str());
  
  if(fDebug>=2){
    cout << "err " << err << endl;
    cout << "nspecs "<< fNSpecs << endl;
    cout << specs_str.c_str() << endl;
  }
  
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
  for(size_t i_spec = 0; i_spec<fNSpecs; i_spec++){
    TSpectroInfo specinfo;
    //specinfo.SetName(fSpecNames.at(i_spec));
    double mcangle;
    int ndets;
    string dets_str;
    int nsig;
    
    string prefix2 = prefix+fSpecNames.at(i_spec)+".";
    
    if(fDebug>=3){
      cout << "spec info prefix " << prefix2.c_str() << endl;
    }
    
    std::vector<int>* pid = 0;
    std::vector<int>* tid = 0;
    
    try{
      pid = new vector<int>;
      tid = new vector<int>;
      
      DBRequest request[] = {
	{"mcangle",        &mcangle,   kDouble,   0, 0},
	{"nsignal",        &nsig,      kInt,      0, 0},
	{"signal.pid",     pid,        kIntV,     0, 0},
	{"signal.tid",     tid,        kIntV,     0, 0},
	{"ndets",          &ndets,     kInt,      0, 0},
	{"detnames",       &dets_str,  kString,   0, 0},
	{ 0 }
      };
      
      //Int_t err = LoadDB (input, request, prefix2);
      Int_t err = LoadDB (file, GetInitDate(), request, prefix2.c_str());
      
      if (err){
	//input.close();
	fclose(file);
	return kInitError;
      }
      
      specinfo.SetMCAngle(mcangle);
      specinfo.SetNDets(ndets);
      std::vector<std::string> detnames = vsplit(dets_str);
      for(size_t i_str = 0; i_str<detnames.size(); i_str++){
	specinfo.AddDetName(detnames.at(i_str));
      }
      
      if(fDebug>=3){
	cout << " spec " << i_spec << ": ndetectors = " << specinfo.NDets() << endl;
      }
      
      for(Int_t i_sig = 0; i_sig<nsig; i_sig++){
	TSignalInfo siginfo(pid->at(i_sig), tid->at(i_sig));
	specinfo.AddMCSignalInfo(siginfo);
      }
      fSpectroInfo.push_back(specinfo);
	
      delete pid;
      delete tid;
    }  catch(...) {
      delete pid;
      delete tid;
      //input.close();
      fclose(file);
      throw;
    }//end try / catch
    
    if(fDebug>=3){
      cout << "after 'try': spec " << i_spec << ": ndetectors = " << specinfo.NDets() << endl;
    }
    
    // then loop on detectors
    for(size_t i_det = 0; i_det<specinfo.NDets(); i_det++){
      err = LoadDetInfo(fSpecNames.at(i_spec), specinfo.DetName(i_det));
      if(err) {
        exit(3);
      }
    }
    
  }// end spectrometer loop
  //input.close();
  cout << "Background spread time window half width = " << fBkgdSpreadTimeWindowHW << " ns." <<endl;
  
  fclose(file);
  return(kOK);
}

//______________________________________________________________
Int_t TSBSDBManager::LoadDetInfo(const string& specname, const string& detname)
{
  const char *here = "TSBSDBManager::LoadDetInfo()";
  TDetInfo detinfo(detname,specname);
  // Use DB_DIR (standard Hall A analyzer DB path) for
  // common values shared by both SBS-Offline and digitization library
  // and then use SBS_DIGI_DB to search for digitization specific
  // values.
  // Include DB_DIR (standard Hall A analyzer DB path in the search)
  // If any is not specified, the default path is ../db/
  std::string path = "../db/";
  std::string pathCommon = "../db/";
  if(std::getenv("SBS_DIGI_DB")) {
    path = std::string(std::getenv("SBS_DIGI_DB")) + "/";
  }
  if(std::getenv("DB_DIR")) {
    pathCommon = std::string(std::getenv("DB_DIR")) + "/";
  }
  const string& fileName = path+"db_"+specname+"."+detname+".dat";
  const string& fileNameCommon = pathCommon+"db_"+specname+"."+detname+".dat";

  const string prefix = specname+"."+detname+".";
  // First, open the common db file and parse info there, later, the
  // digitization specific db can be used to override any values
  FILE* fileCommon  = OpenFile( fileNameCommon.c_str(), GetInitDate() );
  std::vector<int> detmap,chanmap;
  int chanmap_start = 0;

  // variables from common db (in SBS-offline)
  DBRequest requestCommon[] = {
    {"detmap", &detmap, kIntV, 0, true}, ///< Optional
    {"chanmap", &chanmap, kIntV, 0, true}, ///< Optional
    {"chanmap_start", &chanmap_start, kInt, 0, true}, ///< Optional
    { 0 }
  };
  Int_t err = LoadDB (fileCommon, GetInitDate(), requestCommon, prefix.c_str());
  // Could close the common file already
  fclose(fileCommon);
  
  FILE* file = OpenFile( fileName.c_str(), GetInitDate() );

  string dettype_str;
  int det_id = -1;
  int nchan, nlog_chan = 0;
  int first_crate = -1, first_slot = 0;
  bool ignore_slotcrate = true;
  int slot_per_crate = 22, chan_per_slot = 16; //< Just an arbitrary values

  //First load the parameters which will be common to *all* detectors (including digitization parameters)
  // nlog_chan is number of "logical" channels, in case there
  // are multiple channels that correspond to one physical detector
  // component.  For example, HCAL will have one ADC channel and one TDC
  // channel for each block.
  DBRequest requestDetType[] = {
    {"dettype",        &dettype_str,    kString,  0, 0},
    ///< REQUIRED! See g4sbs_types.h for list of unique IDs
    {"unique_id",      &det_id,    kInt,  0, 0},
    {"detmap", &detmap, kIntV, 0, true}, ///< Optional (override detmap)
    {"chanmap", &chanmap, kIntV, 0, true}, ///< Optional (override)
    {"chanmap_start", &chanmap_start, kInt, 0, true}, ///< Optional (override)
    { 0 }
  };

  err = LoadDB (file, GetInitDate(), requestDetType, prefix.c_str());
  if(err)
    return kInitError;

  if(det_id < 0) {
    std::cerr << "ERROR: Detector ID cannot be negative." << std::endl;
    return kInitError;
  }
  if (IsDetInfoAvailableById(det_id)) { // Check that the det_id is unique
    TDetInfo conflict_det = GetDetInfoById(det_id);
    std::cerr << "Detector: "<< specname << "." << detname << " with ID: "
      << det_id << " is not unique! Conflicts with: "
      << conflict_det.DetFullName() << ", with ID: "
      << conflict_det.DetUniqueId() << std::endl;
    return kInitError;
  }
  detinfo.SetDetUniqueId(det_id);

  if(dettype_str.compare("HCal")==0)detinfo.SetDetType(kHCal);
  if(dettype_str.compare("ECal")==0)detinfo.SetDetType(kECal);
  if(dettype_str.compare("Cher")==0) detinfo.SetDetType(kCher);
  if(dettype_str.compare("Scint")==0) detinfo.SetDetType(kScint);
  // GEMs are their own thing for now...so process them differently
  if(dettype_str.compare("GEM")==0) {
    detinfo.SetDetType(kGEM);
    // Create a new DB manager for this detector
    TGEMSBSDBManager *gemdb = new TGEMSBSDBManager(specname.c_str(),
        detname.c_str());
    gemdb->SetDBFileName(fileName);
    gemdb->LoadGeneralInfo(fileName);
    gemdb->LoadGeoInfo(fileNameCommon);
    gemdb->InitializeGEMs();
    detinfo.SetGEMDB(gemdb);
    nchan = gemdb->GetNChan();
    // Now, process the channel map for these GEMs
    int crate,slot,mpd_id,gem_id,adc_id,i2c,gem_pos,gem_invert;
    std::vector<int> gem_map = gemdb->GetChanMap();
    for(size_t k = 0; k < gem_map.size(); k+=8) {
      crate      = gem_map[k  ];
      slot       = gem_map[k+1];
      mpd_id     = gem_map[k+2];
      gem_id     = gem_map[k+3];
      adc_id     = gem_map[k+4];
      i2c        = gem_map[k+5];
      gem_pos    = gem_map[k+6];
      gem_invert = gem_map[k+7];
      detinfo.AddGEMSlot(crate,slot,mpd_id,gem_id,adc_id,i2c,gem_pos,
          gem_invert);
    }
    fDetInfo.push_back(detinfo);
    return err;
  }

  if(detmap.empty()) { // If no detmap, then user MUST specify slot/crate info
    ignore_slotcrate = false;
  }

  DBRequest request[] = {
    {"nchan",          &nchan,          kInt,     0, 0}, ///<REQUIRED
    {"nlog_chan",      &nlog_chan,      kInt,     0, true},
    {"chan_per_slot",  &chan_per_slot,  kInt,     0, true},
    {"slot_per_crate", &slot_per_crate, kInt,     0, true},
    {"first_crate", &first_crate, kInt,     0, ignore_slotcrate},
    {"first_slot", &first_slot, kInt,     0, ignore_slotcrate},
    { 0 }
  };
  err = LoadDB (file, GetInitDate(), request, prefix.c_str());
  if(nlog_chan == 0) {
    nlog_chan = nchan;
  }

  detinfo.SetNChan(nchan);
  detinfo.SetNLogChan(nlog_chan);

  // Now see if we have a detector and channel map specified,
  // if so, we'll use that instead of the chan_per_slot or slot_per_crate
  // values.
  if(detmap.empty()) { // No detmap, so build simple one with values defined
    detinfo.SetChanPerSlot(chan_per_slot);
    detinfo.SetFirstSlot(first_slot);
    detinfo.SetSlotPerCrate(slot_per_crate);
    detinfo.SetFirstCrate(first_crate);
  } else {
    int crate,slot,ch_lo,ch_hi,chan_count, ch_count;
    chan_count = 0;
    for(size_t k = 0; k < detmap.size(); k+=4) {
      ch_count = 0;
      crate  = detmap[k];
      slot   = detmap[k+1];
      ch_lo  = detmap[k+2];
      ch_hi  = detmap[k+3];
      ch_count = 1 + (ch_hi-ch_lo);
      // Check to make sure numbers make sense
      if(ch_count <= 0) {
        Error(Here(here), "Cannot specify detmap where first channel (%d) is "
            "larger than last channel (%d)",ch_lo,ch_hi);
        err = kInitError;
      }
      detinfo.AddSlot(crate,slot,ch_lo,ch_hi);
      chan_count += ch_count;
    }
    if(chan_count != nlog_chan) {
      Error(Here(here), "Number of logical channels read in detmap (%d) does "
          "not match expected (%d)",chan_count,nlog_chan);
      err = kInitError;
    }
  }

  if(err)
    return err;

  // If the user specified a channel map, build that now
  if(!chanmap.empty()) {
    if (int(chanmap.size()) != nlog_chan) {
      Error(Here(here), "Number of logical channels read in chanmap (%d) does "
          "not match expected (%d)",int(chanmap.size()),nlog_chan);
      err = kInitError;
    } else {
      detinfo.LoadChannelMap(chanmap,chanmap_start);
    }
  }

  const string digprefix = "dig."+prefix;
  
  bool ignore_pmt = false;
  if(detinfo.DetType()==kGEM) ignore_pmt = true;

  bool ignore_adc = false;
  if(detinfo.DetType()==kScint || detinfo.DetType()==kCher) ignore_adc = true;
  ignore_adc = (ignore_pmt || ignore_adc);

  bool ignore_tdc = false;
  if(detinfo.DetType()==kECal) ignore_tdc = true;
  ignore_tdc = (ignore_pmt || ignore_tdc);
  
  TDigInfo diginfo;
  
  double roimp;
  double adcconv = -1;
  int    adcbits = -1;
  double tdcconv = -1;
  int    tdcbits = -1;
  std::vector<double>* gain = 0;
  std::vector<double>* pedestal = 0;
  std::vector<double>* pednoise = 0;
  std::vector<double>* threshold = 0;
  double triggerjitter;
  double triggeroffset;
  double gatewidth;
  double spe_tau;
  double spe_sig;
  double spe_transit;
  string tdc_encoder_str;
  string adc_encoder_str;
  
  try{
    gain = new vector<double>;
    pedestal = new vector<double>;
    pednoise = new vector<double>;
    threshold = new vector<double>;
    
    DBRequest request_dig[] = {
      {"roimpedance",   &roimp,          kDouble,  0, ignore_pmt},
      {"adcconversion", &adcconv,        kDouble,  0, ignore_adc},
      {"adcbits",       &adcbits,        kInt,     0, ignore_adc},
      {"adc_encoder",   &adc_encoder_str,kString,  0, ignore_adc},
      {"tdcconversion", &tdcconv,        kDouble,  0, ignore_tdc},
      {"tdcbits",       &tdcbits,        kInt,     0, ignore_tdc},
      {"tdc_encoder",   &tdc_encoder_str,kString,  0, ignore_tdc},
      {"gain",          gain,            kDoubleV, 0, 0},
      {"pedestal",      pedestal,        kDoubleV, 0, 0},
      {"pednoise",      pednoise,        kDoubleV, 0, 0},
      {"threshold",     threshold,       kDoubleV, 0, 0},
      {"triggerjitter", &triggerjitter,  kDouble,  0, 0},
      {"triggeroffset", &triggeroffset,  kDouble,  0, 0},
      {"gatewidth",     &gatewidth,      kDouble,  0, 0},
      {"spe_tau",       &spe_tau,        kDouble,  0, ignore_pmt},
      {"spe_sigma",     &spe_sig,        kDouble,  0, ignore_pmt},
      {"spe_transit",   &spe_transit,    kDouble,  0, ignore_pmt},
      { 0 }
    };
    
    if(fDebug>=3){
      cout << "spec " << specname.c_str() << " detector " << detname.c_str() << " dig prefix " << digprefix.c_str() << endl;
    }
    
    //err = LoadDB (input, request_dig, digprefix);
    err = LoadDB (file, GetInitDate(), request_dig, digprefix.c_str());
    if (err){
      //input.close();
      fclose(file);
      return kInitError;
    }
    
    if(gatewidth>fBkgdSpreadTimeWindowHW)fBkgdSpreadTimeWindowHW = gatewidth;
    
    //detinfo.DigInfo()
    diginfo.SetROImpedance(roimp);
    diginfo.SetADCConversion(adcconv);
    diginfo.SetADCBits(adcbits);
    diginfo.SetTDCConversion(tdcconv);
    diginfo.SetTDCBits(tdcbits);
    diginfo.SetTriggerJitter(triggerjitter);
    diginfo.SetTriggerOffset(triggeroffset);
    diginfo.SetGateWidth(gatewidth);
    diginfo.SetSPE_Tau(spe_tau);
    diginfo.SetSPE_Sigma(spe_sig);
    diginfo.SetSPE_TransitTime(spe_transit);
    TSBSSimDataEncoder *enc = 0;
    if(!adc_encoder_str.empty()) {
      enc = TSBSSimDataEncoder::GetEncoderByName(
            adc_encoder_str.c_str());
      if(enc && enc->IsADC()) {
        diginfo.SetEncoderADC(enc);
      } else {
        std::cerr << "Error: ADC encoder " << adc_encoder_str << " not found!"
          << std::endl;
        return kInitError;
      }
    }
    if(!tdc_encoder_str.empty()) {
      enc = TSBSSimDataEncoder::GetEncoderByName(
            tdc_encoder_str.c_str());
      if(enc && enc->IsTDC()) {
        diginfo.SetEncoderTDC(enc);
      } else {
        std::cerr << "Error: TDC encoder " << tdc_encoder_str << " not found!"
          << std::endl;
        return kInitError;
      }
    }
    
    if(fDebug>=3){
      cout << "roimp " << roimp << ", DigInfo ROinmpedance " << diginfo.ROImpedance() << endl;
      cout << "DigInfo Gain size " << diginfo.GainSize() << endl;
    }
   
    if(gain->size()!=detinfo.NChan()){
      cout << "warning: number of gains in input (" << gain->size() 
	   << ") does not match number of channels (" << detinfo.NChan()
	   << ")" << endl;
      cout << "First gain entry used for all channels. " << endl
	   << "If you want one gain value per channel, fix your DB" << endl<< endl;
      diginfo.AddGain(gain->at(0));
      if(fDebug>=3)cout << "Gain size (right after 'AddGain') = "<< diginfo.GainSize() << endl;
    }else{
      for(uint i__ = 0; i__<gain->size(); i__++){
	diginfo.AddGain(gain->at(i__));
      }
    }
    if(pedestal->size()!=detinfo.NChan()){
      cout << "warning: number of pedestals in input (" << pedestal->size() 
	   << ") does not match number of channels (" << detinfo.NChan()
	   << ")" << endl;
      cout << "First pedestal entry used for all channels. " << endl
	   << "If you want one pedestal value per channel, fix your DB" << endl<< endl;
      diginfo.AddPedestal(pedestal->at(0));
    }else{
      for(uint i__ = 0; i__<pedestal->size(); i__++){
	diginfo.AddPedestal(pedestal->at(i__));
      }
    }
    if(pednoise->size()!=detinfo.NChan()){
      cout << "warning: number of pedestal noises in input (" << pednoise->size() 
	   << ") does not match number of channels (" << detinfo.NChan()
	   << ")" << endl;
      cout << "First pedestal noise entry used for all channels. " << endl
	   << "If you want one pedestal noise value per channel, fix your DB" << endl<< endl;
      diginfo.AddPedestalNoise(pednoise->at(0));
    }else{
      for(uint i__ = 0; i__<pednoise->size(); i__++){
	diginfo.AddPedestalNoise(pednoise->at(i__));
      }
    }
    
    if(threshold->size()!=detinfo.NChan()){
      cout << "warning: number of thresholds in input (" << gain->size() 
	   << ") does not match number of channels (" << detinfo.NChan()
	   << ")" << endl;
      cout << "First threshold entry used for all channels. " << endl
	   << "If you want one threshold value per channel, fix your DB" << endl<< endl;
      diginfo.AddThreshold(threshold->at(0));
    }else{
      for(uint i__ = 0; i__<threshold->size(); i__++){
	diginfo.AddThreshold(threshold->at(i__));
      }
    }
    
    if(fDebug>=3){
      cout << "Dig vectors sizes: gain " << gain->size() 
	   << " DigInfo.GainSize() " << diginfo.GainSize() << endl
	   << ", pedestal " << pedestal->size() 
	   << " DigInfo.PedestalSize() " << diginfo.PedestalSize() << endl
	   << ", pedestal noise " << pednoise->size() 
	   << " DigInfo.PedestalNoiseSize() " << diginfo.PedestalNoiseSize() << endl
	   << ", threshold " << threshold->size() 
	   << " DigInfo.ThresholdSize() " << diginfo.ThresholdSize() << endl
	   << endl;
    }
    detinfo.SetDigInfo(diginfo);

    delete gain;
    delete pedestal;
    delete pednoise;
    delete threshold;
  }  catch(...) {
    delete gain;
    delete pedestal;
    delete pednoise;
    delete threshold;
    //input.close();
    fclose(file);
    throw;
  }//end try / catch

  if(fDebug>=3){
    cout << "After try: Dig vectors sizes: DigInfo.GainSize() " << detinfo.DigInfo().GainSize() << endl
	 << ", DigInfo.PedestalSize() " << detinfo.DigInfo().PedestalSize() << endl
	 << ", DigInfo.PedestalNoiseSize() " << detinfo.DigInfo().PedestalNoiseSize() << endl
	 << ", DigInfo.ThresholdSize() " << detinfo.DigInfo().ThresholdSize() << endl
	 << endl;
  }
    
  const string geoprefix = "geo."+prefix;
  
  if(detinfo.DetType()==kGEM || detinfo.DetType()==kScint){
    UInt_t nplanes;
    std::vector<Int_t>* nmodules = 0;
    
    try{
      nmodules = new vector<int>;
      DBRequest request[] = {
	{"nplanes",        &nplanes,        kInt,     0, 0},
	{"nmodules",       nmodules,        kIntV,    0, 0},
	{ 0 }
      };
       
      if(fDebug>=3){
	cout << " prefix.c_str() " << prefix.c_str() << endl;
      }
      
      //Int_t err = LoadDB (input, request, prefix);
      Int_t err = LoadDB (file, GetInitDate(), request, prefix.c_str());
      
      if (err){
	//input.close();
	fclose(file);
	return kInitError;
      }
      
      detinfo.SetNPlanes(nplanes);
      
      
      for(size_t i_pl = 0; i_pl<nplanes; i_pl++){
	detinfo.AddNModules(nmodules->at(i_pl));
	
	for(Int_t i_mod = 0; i_mod<nmodules->at(i_pl); i_mod++){
	  TGeoInfo thisGeo;
	  int nrows;
	  int ncols;
	  double xsize;
	  double ysize;
	  double zpos;
	  double xoffset;
	  double yoffset;
	  
	  DBRequest request_geo[] = {
	    {"nrows",     &nrows,      kInt,    0, 1},
	    {"ncols",     &ncols,      kInt,    0, 1},
	    {"xsize",     &xsize,      kDouble, 0, 1},
	    {"ysize",     &ysize,      kDouble, 0, 1},
	    {"zpos",      &zpos,       kDouble, 0, 1},
	    {"xoffset",   &xoffset,    kDouble, 0, 1},
	    {"yoffset",   &yoffset,    kDouble, 0, 1},
	    { 0 }
	  };
	  
	  if(fDebug>=3){
	    cout << " geoprefix.c_str() " << geoprefix.c_str() << endl;
	  }
	  
	  string geoprefix_ii = geoprefix;
	  char temp[100];
	  if(nplanes>1){
	    sprintf(temp, "%s%lu", geoprefix_ii.c_str(), (i_pl+1));
	    geoprefix_ii = std::string(temp)+".";
	  }
	  if(nmodules->at(i_pl)>1){
	    sprintf(temp, "%s%d", geoprefix_ii.c_str(), (i_mod+1));
	    geoprefix_ii = std::string(temp)+".";
	  }
	  // if(nplanes>1) geoprefix_ii = geoprefix_ii+std::to_string(i_pl+1)+".";
	  // if(nmodules->at(i_pl)>1) geoprefix_ii = geoprefix_ii+std::to_string(i_mod+1)+".";
	  	  
	  if(fDebug>=3){
	    cout << " geoprefix_ii.c_str() " << geoprefix_ii.c_str() << endl;
	  }
	  
	  //err = LoadDB (input, request_geo, geoprefix+"."+std::to_string(i_pl));
	  err = LoadDB (file, GetInitDate(), request_geo, geoprefix_ii.c_str());
	  if (err){
	    //input.close();
	    fclose(file);
	    return kInitError;
	  }
	  
	  thisGeo.SetNRows(nrows);
	  thisGeo.SetNCols(ncols);
	  thisGeo.SetXSize(xsize);
	  thisGeo.SetYSize(ysize);
	  thisGeo.SetZPos(zpos);
	  thisGeo.SetXOffset(xoffset);
	  thisGeo.SetYOffset(yoffset);
	  
	  detinfo.AddGeoInfo(thisGeo);
	  if(fDebug>=3)cout << "GeoInfo size " << detinfo.GeoInfoSize() << endl;
	  //detinfo.fGeoInfo.push_back(thisGeo);
	}
      }//end 
      delete nmodules;
    }  catch(...) {
      delete nmodules;
      //input.close();
      fclose(file);
      throw;
    }//end try / catch
    
    if(fDebug>=3){
      cout << "NModules size " << detinfo.NModulesSize() << endl;
    }
    
  }else{
    TGeoInfo thisGeo;
    int nrows;
    int ncols;
    double xsize;
    double ysize;
    double zpos;
    
    DBRequest request_geo[] = {
      {"nrows",     &nrows,      kInt,    0, 1},
      {"ncols",     &ncols,      kInt,    0, 1},
      {"xsize",     &xsize,      kDouble, 0, 1},
      {"ysize",     &ysize,      kDouble, 0, 1},
      {"zpos",      &zpos,       kDouble, 0, 1},
      { 0 }
    };
    
    err = LoadDB (file, GetInitDate(), request_geo, geoprefix.c_str());
    if (err){
      //input.close();
      fclose(file);
      return kInitError;
    }
    
    thisGeo.SetNRows(nrows);
    thisGeo.SetNCols(ncols);
    thisGeo.SetXSize(xsize);
    thisGeo.SetYSize(ysize);
    thisGeo.SetZPos(zpos);
    
    detinfo.AddGeoInfo(thisGeo);
    if(fDebug>=3)cout << "GeoInfo size " << detinfo.GeoInfoSize() << endl;
    //detinfo.fGeoInfo.push_back(thisGeo);
  }
  
  
  if(fDebug>=3){
    cout << "GeoInfo size " << detinfo.GeoInfoSize() << endl;
    cout << "GeoInfo ZPos " << detinfo.GeoInfo(0).ZPos() << endl;
    cout << "DigInfo ROinmpedance " << detinfo.DigInfo().ROImpedance() << endl;
    cout << "DigInfo Gain size " << detinfo.DigInfo().GainSize() << endl;
  }
   
  fDetInfo.push_back(detinfo);
  
  //input.close();
  fclose(file);
  return(kOK);
}

bool TSBSDBManager::IsDetInfoAvailableById(Int_t id)
{
  // Loop through all detectors to see if this one is available
  for(size_t i = 0; i<fDetInfo.size(); i++){
    if(fDetInfo.at(i).DetUniqueId()==id){
      return true;
    }
  }

  return false;
}

bool TSBSDBManager::IsDetInfoAvailable(const char* detname)
{
  // Loop through all detectors to see if this one is available
  for(size_t i = 0; i<fDetInfo.size(); i++){
    if(fDetInfo.at(i).DetName().compare(detname)==0 ||
        fDetInfo.at(i).DetFullName().compare(detname)==0){
      return true;
    }
  }

  return false;
}


const TDetInfo& TSBSDBManager::GetDetInfoById(Int_t id)
{
  //check if the detectors databases is loaded in the first place
  if(fDetInfo.size()==0){
    std::cerr << "Detector info has not been loaded by the manager yet; Exiting." << endl;
    std::cerr << "To avoid this, make sure you load the DB information with the DB manager before you declare any detector." << endl;
    exit(2);
  }
  // if so, loop on list of detectors.
  for(size_t i = 0; i<fDetInfo.size(); i++){
    if(fDetInfo.at(i).DetUniqueId() == id){
      return fDetInfo.at(i);
    }
  }

  // if no valid detectors have been found exit with error message
  cout << "No detector corresponding to id:" << id << " found in database. Check program or database" << endl;
  exit(2);
}

//function to retrieve the coorect detector information from the detector name
const TDetInfo & TSBSDBManager::GetDetInfo(const char* detname)
{
  //check if the detectors databases is loaded in the first place
  if(fDetInfo.size()==0){
    cout << "Detector info has not been loaded by the manager yet; Exiting." << endl;
    cout << "To avoid this, make sure you load the DB information with the DB manager before you declare any detector." << endl;
    exit(2);
  }
  
  // if so, loop on list of detectors.
  for(size_t i = 0; i<fDetInfo.size(); i++){
    if(fDetInfo.at(i).DetName().compare(detname)==0){
      return fDetInfo.at(i);
    }
  }
  
  // if no valid detectors have been found exit with error message
  cout << "No detector corresponding to " << detname << " found in database. Check program or database" << endl;
  exit(2);
}

ClassImp(TSBSDBManager)

/*
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
        case kInt:
          sv >> *((Int_t*)item->var);
          break;
        case kString:
          sv >> *((string*)item->var);
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
*/

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
  // Include SBS_DIGI_DB (standard Hall A analyzer DB path in the search)
  std::string path = "../db/";
  if(std::getenv("SBS_DIGI_DB")) {
    path = std::string(std::getenv("SBS_DIGI_DB")) + "/";
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
