#include "TSBSSimGEM.h"
#include <iostream>
#include <TSBSSimData.h>
#include <TSBSSimEvent.h>
#include "TSBSDBManager.h"
#include "SBSSimDataEncoder.h"
#include "TGEMSBSSimDigitization.h"
#include "TGEMSBSSpec.h"
#include "TGEMSBSGEMSimHitData.h"
#include "TGEMSBSGEMChamber.h"
#include "TGEMSBSDBManager.h"

ClassImp(TSBSSimGEM) // Implements TSBSSimGEM


TSBSSimGEM::TSBSSimGEM(const char* name, short id)
{
  fName = name;
  SetUniqueDetID(id);
  Init();
}

TSBSSimGEM::~TSBSSimGEM()
{
}

void TSBSSimGEM::Init()
{
  TSBSSimDetector::Init();

  if(fDebug>=1)
    cout << "GEM detector with UniqueDetID = " << UniqueDetID() << ": TSBSSimGEM::Init() " << endl;

  // Get the Detector info
  fDetInfo = fDBmanager->GetDetInfo(fName.Data());
  fManager = fDetInfo.GetGEMDB();
  fManager->SetDebug(fDebug);

  // Hard code the MPD encoder for GEMs
  fEncoderMPD = SBSSimDataEncoder::GetEncoderByName("mpd");

  // At this point, we should build all the GEM chambers that exist for
  // this 

  fGEMDigi = new TGEMSBSSimDigitization(fManager->GetSpec(),fName,fManager);
  fGEMDigi->SetDebug(fDebug);
}


void TSBSSimGEM::LoadEventData(const std::vector<g4sbshitdata*> &evbuffer)
{
  Clear();
  fGEMDigi->EventStart();//Put it here !!!
  LoadAccumulateData(evbuffer);
}

void TSBSSimGEM::LoadAccumulateData(const std::vector<g4sbshitdata*> &evbuffer)
{
  // Pack data into TGEMSBSGEMSimHitData
  //    printf("NEXT EVENT ---------------------------\n");

  TGEMSBSGEMSimHitData gd;
  gd.ClearEvent();
  gd.SetEvent(GetEventNum());

  // First, build a list of events that correspond to this detector only
  std::vector<g4sbshitdata*> evb;
  for(std::vector<g4sbshitdata*>::const_iterator it = evbuffer.begin();
      it!= evbuffer.end(); ++it ) {
    if((*it)->GetDetUniqueID() == UniqueDetID()) {
      evb.push_back(*it);
    }
  }

  if(evb.empty()) {
    return;
  }
  gd.InitEvent(evb.size());

  g4sbshitdata *h;//, *hs;
  //bool matchedstrip;
  bool first = true;
  unsigned int ngdata = 0;// j,
  double hit_data_temp[24];
  int module = 0;
  int sector = 0;
  int i = 0;
  for(std::vector<g4sbshitdata*>::const_iterator it = evb.begin();
      it!= evb.end(); ++it ) {
    h = (*it);
    
    if(h->GetDetUniqueID() == UniqueDetID()) {
      // Unpack the data in its entirety (for now)
      int      source  = h->GetData(0);
      int       plane   = h->GetData(1)-1;
      //double    strip   = h->GetData(2); 
      //double    x       = h->GetData(3); 
      //double    y       = h->GetData(4); 
      //double    z       = h->GetData(5); 
      //double    polx    = h->GetData(6); 
      //double    poly    = h->GetData(7); 
      //double    polz    = h->GetData(8); 
      //double    t       = h->GetData(9); 
      //double    trms    = h->GetData(10);
      double    tmin    = h->GetData(11);
      double    tmax    = h->GetData(12);
      double    tx      = h->GetData(13);
      double    ty      = h->GetData(14);
      double    txp     = h->GetData(15);
      double    typ     = h->GetData(16);
      //double    xg      = h->GetData(17);
      //double    yg      = h->GetData(18);
      //double    zg      = h->GetData(19);
      int       trid    = h->GetData(20);
      double    mid     = h->GetData(21);
      int       pid     = h->GetData(22);
      double    vx      = h->GetData(23);
      double    vy      = h->GetData(24);
      double    vz      = h->GetData(25);
      double    p       = h->GetData(26);
      double    edep    = h->GetData(27)*1.0e3;
      //double    beta    = h->GetData(28);
      double    xin     = h->GetData(29);
      double    yin     = h->GetData(30);
      double    zin     = h->GetData(31);
      double    xout    = h->GetData(32);
      double    yout    = h->GetData(33);
      double    zout    = h->GetData(34);
      int type = mid+1; //=1 if primary, >1 if secondary

      module = fManager->GetModuleIDFromPos(plane,tx);
      if(module==-1)continue;
      
      double pz = sqrt( pow(p, 2)/
			( pow(txp, 2) + 
			  pow(typ, 2) + 1.0) );

      TVector3 Mom = TVector3(txp*pz*1.0e3, // in MeV
			      typ*pz*1.0e3, // in MeV
			      pz*1.0e3);// in MeV

      TVector3 X_in = TVector3((xin-fManager->GetXOffset(plane, module))*1.0e3, // in mm
			       yin*1.0e3, // in mm
			       (zin+fManager->Getg4sbsZSpecOffset()-
				fManager->GetD0(plane, module))*1.0e3);// in mm

      TVector3 X_out = TVector3((xout-fManager->GetXOffset(plane, module))*1.0e3, // in mm
				yout*1.0e3, // in mm
				(zout+fManager->Getg4sbsZSpecOffset()-
				 fManager->GetD0(plane, module))*1.0e3);// in mm

      // we use X_in and X_out to extrapolate X_RO; 
      // not very clean, but since we don't use it in the digitization at all, it does not really matter...
      TVector3 X_RO = TVector3(X_in.X()+(xout-xin)*9.185/3.0 ,
			       // in mm
			       X_in.Y()+(yout-yin)*9.185/3.0 ,
			       // in mm
			       +7.685);// in mm
      // see comment lines 138-141 of TGEMSBSGeant4File.h

      if(fabs(X_out.X())>=fManager->GetDX(plane, module)*5.0e2){
	if(fDebug>1){
	  cout << "Warning: Evt " << fEvNum << ", hit " << i 
	       << ": X_out.X " << X_out.X() << " outside BBGEM plane " << plane
	       << " sector " << sector;
	} //D_FLAG
	double temp = fabs(X_out.X());
	X_out[0]*=fManager->GetDX(plane, module)*5.0e2/temp;
	if(fDebug>1){
	  cout  << "; set at limit: " << X_out.X() << " mm " << endl;
	  cout << "(X_in.X = " << X_in.X() << ",  " << tx*1.0e3 << " mm)" << endl;
	} //D_FLAG
	X_RO.SetX(X_out.X());
      }  
      if(fabs(X_out.Y())>=fManager->GetDY(plane, module)*5.0e2){
	if(fDebug>1){
	  cout << "Warning: Evt " << fEvNum << ", hit " << i 
	       << ": X_out.Y " << X_out.Y() << " outside FT plane " << plane << " sector " << sector;
	} //D_FLAG
	double temp = fabs(X_out.Y());
	X_out[1]*=fManager->GetDY(plane, module)*5.0e2/temp;
	if(fDebug>1){
	  cout  << "; set at limit: " << X_out.Y() << " mm " << endl;
	  cout << "(X_in.Y = " << X_in.Y() << ",  " << ty*1.0e3 << " mm)" << endl;
	} //D_FLAG
	X_RO.SetY(X_out.Y());
      }

      TVector3 Vtx = TVector3(vx*1.0e3, // in mm
			      vy*1.0e3, // in mm
			      vz*1.0e3);// in mm

      // Then copy them onto a temporary hit array
      //Filling hit_data temporary array...
      hit_data_temp[0] = (double)plane;
      hit_data_temp[1] = edep;
      hit_data_temp[8] = tmin;
      hit_data_temp[12] = tmax;
      hit_data_temp[13] = (double)type;
      hit_data_temp[17] = (double)trid;
      hit_data_temp[18] = (double)pid;
      hit_data_temp[19] = module;
      hit_data_temp[23] = sector;
      for(int k = 0; k<3; k++){
	hit_data_temp[k+2] = X_RO[k];
	hit_data_temp[k+5] = X_in[k];
	hit_data_temp[k+9] = X_out[k];
	hit_data_temp[k+14] = Vtx[k];
	hit_data_temp[k+20] = Mom[k];
      }

      ///  for(i=0; i<evbuffer.size(); i++){
      //    h = GetHitData(i);
      if(first) {
	gd.SetSource(source);
	first = false;
      }

      if( hit_data_temp[1]>0.0 ){ // edep
	
	// Chamber IDs are numbered as 
      	
	//if( h->GetDetID()%100 == __GEM_DRIFT_ID &&  hit_data_temp[1]>0.0 ){
	// Vector information
	TVector3 p(hit_data_temp[20], hit_data_temp[21], hit_data_temp[22]);
	gd.SetMomentum(ngdata, p);
	
	TVector3 li(hit_data_temp[5], hit_data_temp[6], hit_data_temp[7]);
	gd.SetHitEntrance(ngdata, li);
	
	TVector3 lo(hit_data_temp[9], hit_data_temp[10], hit_data_temp[11]);
	gd.SetHitExit(ngdata, lo);
	
	// Average over entrance and exit time
	gd.SetHitTime(ngdata, (hit_data_temp[8]+hit_data_temp[12])/2.0);
	
	TVector3 vert(hit_data_temp[14], hit_data_temp[15], hit_data_temp[16]);
	gd.SetVertex(ngdata, vert);
	
	TVector3 lr(hit_data_temp[2], hit_data_temp[3], hit_data_temp[4]);
	gd.SetHitReadout(ngdata, lr);
	// printf("%d %f %f\n", h->GetDetID()/100, li.X(), li.Y()  );
	
	gd.SetHitEnergy(ngdata, hit_data_temp[1]*1e6 ); // Gives eV
	gd.SetParticleType(ngdata, (UInt_t)hit_data_temp[13] );//  Track type (1 primary, >1 secondary) 
	gd.SetTrackID(ngdata, (UInt_t) hit_data_temp[17] );// track ID
	gd.SetParticleID(ngdata, hit_data_temp[18] );//  PID 
    
	// gd.SetHitChamber(ngdata, hit_data_temp[23]*fManager->GetNChamber()+hit_data_temp[0]);
	gd.SetHitChamber(ngdata, hit_data_temp[0]*3+hit_data_temp[19]);
	//cout<<(hit_data_temp[23]*fManager->GetNChamber()+hit_data_temp[0])<<" : "<<(hit_data_temp[0]*3+hit_data_temp[19])<<endl;

	gd.SetHitPlane(ngdata, hit_data_temp[0]);
	gd.SetHitModule(ngdata,hit_data_temp[19]);
      
	ngdata++;
      }
      i++;
    }
  }//end loop on g4sbshitdata
  gd.SetNHit(ngdata);

  //Add the hits in here ??? why not ?
  
  
  // Once this is done, now call fGEMDigi to actually process these hits
  fGEMDigi->SetTimeZero(fTimeZero);//override GEM time
  fGEMDigi->AdditiveDigitize (gd, fManager->GetSpec());
  gd.ClearEvent();
}

void TSBSSimGEM::Digitize(TSBSSimEvent &event)
{
  if(!fEncoderMPD) {
    std::cerr << "Wow! No MPD encoder!!" << std::endl;
    return;
  }
  int plane, module;
  std::string planename;
  fGEMDigi->SetTreeStrips();
  std::vector<uint32_t> data_mpd;
  std::vector<uint32_t> data;
  std::vector<Int_t> data_dec;
  Int_t ADCsum = 0;
  //TSBSSimEvent::DetectorData data;
  //data.fDetID = UniqueDetID();
  // Here data.fChannel corresponds to plane
  //data.fChannel = 0;
  SimEncoder::mpd_data mpd_data;
  mpd_data.channel = 0;
  Int_t adc = 0;
  UInt_t nstrip = 0;
  UInt_t strip, strip0_mpd, pl_strip;
  mpd_data.adc_id = 0; // For now, increment with each APV25
  // The other info will be encoded properly when the Decoder pass happens
  mpd_data.mpd_id = 0;
  mpd_data.gem_id = 0;
  mpd_data.i2c = 0;
  mpd_data.pos = 0;
  mpd_data.invert = 0;
  UInt_t idx = 0;
  SetHasDataFlag(false);
  if(fDebug>=3)cout << fGEMDigi->fGEMClust.size() << endl;
  for(uint i_mc = 0; i_mc<fGEMDigi->fGEMClust.size();i_mc++){
    event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].fNSimHits++;
    event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].fSimSource.push_back(fGEMDigi->fGEMClust[i_mc].fSource);
    event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].fSimTRID.push_back(fGEMDigi->fGEMClust[i_mc].fTRID);
    event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].fSimPID.push_back(fGEMDigi->fGEMClust[i_mc].fPID);
    event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].fSimChannel.push_back(-1);
    event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].fPlane.push_back(fGEMDigi->fGEMClust[i_mc].fPlane);
    event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].fModule.push_back(fGEMDigi->fGEMClust[i_mc].fModule);
    event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].fSimEdep.push_back(fGEMDigi->fGEMClust[i_mc].fCharge);
    event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].fSimTime.push_back(fGEMDigi->fGEMClust[i_mc].fTime);
    event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].fXpos.push_back(fGEMDigi->fGEMClust[i_mc].fHitpos.X());
    event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].fYpos.push_back(fGEMDigi->fGEMClust[i_mc].fHitpos.Y());
    event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].fPX.push_back(fGEMDigi->fGEMClust[i_mc].fPspec.X());
    event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].fPY.push_back(fGEMDigi->fGEMClust[i_mc].fPspec.Y());
    event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].fPZ.push_back(fGEMDigi->fGEMClust[i_mc].fPspec.Z());
    event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].fSizeX.push_back(fGEMDigi->fGEMClust[i_mc].fSize[0]);
    event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].fSizeY.push_back(fGEMDigi->fGEMClust[i_mc].fSize[1]);
    event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].fStartX.push_back(fGEMDigi->fGEMClust[i_mc].fStart[0]);
    event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].fStartY.push_back(fGEMDigi->fGEMClust[i_mc].fStart[1]);
  }

  // Here, Chamber is equivalent to a "Tracking-Plane" which is really
  // what the TGEMSBSSimDigitization uses
  
  for(UInt_t ich = 0; ich < fGEMDigi->GetNChambers(); ich++) {
    fManager->GetPMfromGlobalPlaneNum(ich, plane, module);
    if(fDebug>=3)cout << "ich " << ich << " plane " << plane << " module " << module << endl;
    for(UInt_t ip = 0; ip < fGEMDigi->GetNPlanes(ich); ip++) {
      mpd_data.nsamples = fGEMDigi->GetNSamples(ich,ip);
      // This is the total number of APV25's we'd need
      nstrip = fGEMDigi->GetNStrips(ich,ip);
      // two options on how to name the plane here: 
      // plane and module info specifically, and then actual strip number
      // planename = Form("%s.p%d.m%d.%s", 
      // 		       fDetInfo.DetFullName().c_str(), 
      // 		       plane+1, module+1, kProj_str[ip].c_str());
      // ... or just use plane, and then strip = strip + nstrips*nmodules
      planename = Form("%s.%d.%s", 
		       fDetInfo.DetFullName().c_str(), 
		       plane+1, kProj_str[ip].c_str());
      if(fDebug>=4)cout << planename.c_str() << " ich " << ich << " ip " << ip << " nstrip " << nstrip << endl;
      
      // Break up the data in number of strips that fit in an APV25
      // (128 channels). I found that TGEMSBSSimDigitization somehow gets one
      // extra strip that seems unreasonable, since I doubt we'd get one
      // APV25 chip just for one strip. Hence, I'm going to assume that's
      // not the intent and skip anything with only one strip left.
      // => Good call! EF
      strip = 0;
      while(nstrip > 1) {
        mpd_data.samples.clear();
        mpd_data.nstrips = nstrip >= SBS_APV25_NCH ? SBS_APV25_NCH : nstrip;
        nstrip -= mpd_data.nstrips;
        mpd_data.samples.resize(mpd_data.nsamples*mpd_data.nstrips);
        idx = 0;
	strip0_mpd = strip;
	//if I have to reloop anyway...
        for(UInt_t istrip = 0; istrip < mpd_data.nstrips; istrip++, strip++) {
	  for(UShort_t s = 0; s < mpd_data.nsamples; s++) {
	    adc = fGEMDigi->GetSimADC(ich,ip,strip,s);
	    // Negative values convert poorly to unsigned integers, so just
	    // set them to zero if the actual ADC is negative
	    // which should not be the case if common mode is switched on
	    mpd_data.samples[idx++] = (adc>0 ? adc : 0);
	  }
	}
	data_mpd.clear();
        fEncoderMPD->EncodeMPD(mpd_data,fEncBuffer,fNEncBufferWords);
        CopyEncodedData(fEncoderMPD,mpd_data.adc_id,data_mpd);
	
	if(fDebug>=5){
	  cout << data_mpd.size() << endl;
	  for(uint i = 0; i< data_mpd.size(); i++)cout << data_mpd[i] << " ";
	  cout << endl;
	}
	
	//Then I fill the tree :)
	// collection of headers first:
	// first header: data word; MPDheader: dataword_samps
	// perhaps it is worth doing one hit per header... //Nope
	data_dec.clear();//data_dec.push_back(0);
	data.clear();//data.push_back(0);
	int i_;
	for(i_ = 1; i_<3; i_++){
	  data.push_back(data_mpd[i_]);
	}
	event.fSimDigSampOutData[planename].fNHits++;
	event.fSimDigSampOutData[planename].fChannel.push_back(-1000);
	event.fSimDigSampOutData[planename].fDataWord.push_back(data_mpd[0]);
	event.fSimDigSampOutData[planename].fADC.push_back(-1000000);
	event.fSimDigSampOutData[planename].fNsamps.push_back(0);
	event.fSimDigSampOutData[planename].fADC_samps.push_back(data_dec);
	event.fSimDigSampOutData[planename].fDataWord_samps.push_back(data);
		
	strip = strip0_mpd;
	
        for(UInt_t istrip = 0; istrip < mpd_data.nstrips; istrip++, strip++) {
	  ADCsum = 0;
	  data.clear();
	  data_dec.clear();
	  pl_strip = fManager->GetGlobalStripPlane(strip, plane, module, ip);
	  
	  //zero suppression here:
	  if(fDebug>=5)cout << fGEMDigi->GetSimADCSum(ich,ip,strip) << " " << fGEMDigi->ZeroSupThreshold(idx) << endl;
	  if(fGEMDigi->GetSimADCSum(ich,ip,strip)>=fGEMDigi->ZeroSupThreshold(mpd_data.mpd_id)){
	    for(UShort_t s = 0; s < mpd_data.nsamples; s++) {
	      adc = fGEMDigi->GetSimADC(ich,ip,strip,s)
		-fGEMDigi->CommonMode(mpd_data.mpd_id);
	      ADCsum+=adc;
	      if(s%2==0)data.push_back(data_mpd[i_++]);
	      data_dec.push_back(adc);
	      if(adc>0)
		SetHasDataFlag(true);
	    }
	    
	    event.fSimDigSampOutData[planename].fNHits++;
	    //event.fSimDigSampOutData[planename].fChannel.push_back(strip);
	    event.fSimDigSampOutData[planename].fChannel.push_back(pl_strip);
	    //this is a bit of abuse: I use the dataword in this case to store the size of the vector samples
	    event.fSimDigSampOutData[planename].fDataWord.push_back(ceil(mpd_data.nsamples/2));//This way if nsamples is odd, we still save the singleton sample
	    assert(data.size()==ceil(mpd_data.nsamples/2));
	    event.fSimDigSampOutData[planename].fADC.push_back(ADCsum);
	    event.fSimDigSampOutData[planename].fNsamps.push_back(mpd_data.nsamples);
	    event.fSimDigSampOutData[planename].fADC_samps.push_back(data_dec);
	    event.fSimDigSampOutData[planename].fDataWord_samps.push_back(data);
	  }
	}
	//event.fDetectorData.push_back(data);
	mpd_data.adc_id++;
      }
      //data.fChannel++;
      //if(fDebug>=3)
      if(!event.fSimDigSampOutData[fDetInfo.DetFullName()].CheckSize(true, fDebug>=1)){
	cout << "Warning: output vectors for" << planename << " don't have the same size!" << endl;
      }
      if(!event.fSimGEMHitMCOutData[fDetInfo.DetFullName()].CheckSize(fDebug>=1)){
	cout << "Warning: MC output vectors for" << planename << " don't have the same size!" << endl;
      }
      
    }
  }
}


// Clear signals in array
void TSBSSimGEM::Clear(Option_t*)
{
  //for(size_t i = 0; i < fSignals.size(); i++ ) {
  //  fSignals[i].Clear();
  //}
}

void TSBSSimGEM::EventStart()
{
  // Tell fGEMDigi to clear out previous event (since now all calls to its
  // digitize are additive)
  fGEMDigi->EventStart();
}
