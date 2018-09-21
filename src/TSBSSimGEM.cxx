#include "TSBSSimGEM.h"
#include <iostream>
#include <TSBSSimData.h>
#include <TSBSSimEvent.h>
#include "TSBSDBManager.h"
#include "TSBSSimDataEncoder.h"
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

  // Hard code the MPD encoder for GEMs
  fEncoderADC = TSBSSimDataEncoder::GetEncoderByName("mpd");
  // Make spectrometer
  fSpec = new TGEMSBSSpec(fManager->GetPrefix().c_str(),
      Form("Temporary GEM spectrometer for %s",fName.Data()));
  // And make all the GEM chambers for this spectrometer
  TGEMSBSGEMChamber *dGEM;
  for(Int_t plane = 0; plane < fManager->GetNGEMPlane(); plane++) {
    for(Int_t mod = 0; mod < fManager->GetNModule(plane); mod++) {
      dGEM = new TGEMSBSGEMChamber(Form("plane%d.module%d",
            /*fManager->GetPrefix().c_str(),*/plane,mod),Form(
            "Test chamber for %s on Plane: %d, Module: %d",
              fName.Data(),plane,mod));
      dGEM->SetApparatus(fSpec);
      if(dGEM->Init()) { // true == error
          std::cerr << "ERROR!: TSBSSimGEM::Init(), UniqueDetID = "
          << UniqueDetID() << " error initializing GEM: " <<
          Form("%s.plane%d.module%d",fManager->GetPrefix().c_str(),plane,mod)
          << std::endl;
      } else {
        fSpec->AddGEM(dGEM);
      }
    }
  }

  // At this point, we should build all the GEM chambers that exist for
  // this 

  fGEMDigi = new TGEMSBSSimDigitization(*fSpec,fName,fManager);
}


void TSBSSimGEM::LoadEventData(const std::vector<g4sbshitdata*> &evbuffer)
{
  Clear();
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
  gd.SetNHit(ngdata);


  // Once this is done, now call fGEMDigi to actually process these hits
  fGEMDigi->AdditiveDigitize (gd, *fSpec);
  gd.ClearEvent();
}

void TSBSSimGEM::Digitize(TSBSSimEvent &event)
{
  fGEMDigi->SetTreeStrips();
  TSBSSimEvent::DetectorData data;
  data.fDetID = UniqueDetID();
  data.fChannel = 0;
  SimEncoder::fadc_data fadc_data;
  Int_t adc = 0;
  for(UInt_t ich = 0; ich < fGEMDigi->GetNChambers(); ich++) {
    for(UInt_t ip = 0; ip < fGEMDigi->GetNPlanes(ich); ip++) {
      for(UInt_t istrip = 0; istrip < fGEMDigi->GetNStrips(ich,ip); istrip++) {
        data.fData.clear();
        fadc_data.samples.clear();
        fadc_data.samples.resize(fGEMDigi->GetNSamples(ich,ip));
        for(UShort_t s = 0; s < fGEMDigi->GetNSamples(ich,ip); s++) {
          adc = fGEMDigi->GetSimADC(ich,ip,istrip,s);
          // Negative values convert poorly to unsigned integers, so just
          // set them to zero if the actual ADC is negative
          fadc_data.samples[s] = (adc>0 ? adc : 0);
        }
        fEncoderADC->EncodeFADC(fadc_data,fEncBuffer,fNEncBufferWords);
        CopyEncodedData(fEncoderADC,0,data.fData);
        event.fDetectorData.push_back(data);
        data.fChannel++;
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
