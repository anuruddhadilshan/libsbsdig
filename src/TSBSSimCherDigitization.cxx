#include "TSBSSimCherDigitization.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TClonesArray.h"

#include "TSBSGeant4File.h"  // needed for g4sbsgendata class def
#include "TSBSSimData.h"
#include "TSBSCherData.h"
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
/*
//for some reasons, if these parameters are declared as flags in the .h, it doesn't work...
Int_t    TSBSSimCherDigitization::fDoCrossTalk = 0;
Int_t    TSBSSimCherDigitization::fNCStripApart = 0;
Double_t TSBSSimCherDigitization::fCrossFactor = 0.;
Double_t TSBSSimCherDigitization::fCrossSigma = 0.;

// Chamber number -> sector/plane helper functions

inline
static void ChamberToSector( Short_t chamber, Short_t& sector, Short_t& plane )
{
  div_t d = div( chamber, manager->GetNChamber() );
  sector = d.quot;
  plane  = d.rem;
  //cout << "chamber " << chamber << ", sector " << sector << ", plane " << plane << endl;
}

inline
static UInt_t MapSector( UInt_t chamber )
{
  // Convert the true chamber index to one with sector = 0
  
  return manager->GetNSector() * UInt_t(chamber/manager->GetNSector());
}
*/

/*
// Auxiliary class

TSBSDigitizedPlane::TSBSDigitizedPlane (UShort_t nstrip,
					UShort_t nsample,
					Int_t threshold )
  : fNSamples(nsample), fNStrips(nstrip), fThreshold(threshold){
  // initialization of all arrays
  fType = new Short_t[nstrip];
  fCharge = new Float_t[nstrip];
  fTime = new Float_t[nstrip];
  fTotADC = new Int_t[nstrip];
  fOverThr = new Short_t[nstrip];

  fStripADC.Set(fNSamples*fNStrips);
  fStripClusters.resize(fNStrips);

  Clear();
};

TSBSDigitizedPlane::~TSBSDigitizedPlane()
{
  delete[] fType;
  delete[] fCharge;
  delete[] fTime;
  delete[] fTotADC;
  delete[] fOverThr;
};

void
TSBSDigitizedPlane::Clear()
{
  fStripADC.Reset();
  memset( fType,   0, fNStrips*sizeof(Short_t) );
  memset( fTotADC, 0, fNStrips*sizeof(Int_t) );
  memset( fCharge, 0, fNStrips*sizeof(Float_t) );

  for (Int_t i = 0; i < fNStrips; i++) {
    fTime[i] = 9999.;
    fStripClusters[i].clear();
  }
  fNOT = 0;
}

void
TSBSDigitizedPlane::Cumulate (const TSolGEMVStrip *vv, Short_t type,
			      Short_t clusterID )
{
  // cumulate hits (strips signals)
  if (vv) {
    for( Int_t j=0; j < vv->GetSize(); j++ ) {
      Int_t idx = vv->GetIdx(j);
      assert( idx >= 0 && idx < fNStrips );
      fType[idx] |= type;
      fTime[idx] = (fTime[idx] < vv->GetTime()) ? fTime[idx] : vv->GetTime();
      fCharge[idx] += vv->GetCharge(j);
      bool was_below = !( fTotADC[idx] > fThreshold );
      for( UInt_t k=0; k<fNSamples; k++ ) {
	Int_t nnn = vv->GetADC(j,k);
	//if(fTime[idx]>50.)cout << k << ":" << nnn << " ";
	assert( nnn >= 0 );
	if( nnn == 0 ) continue;
	Int_t iadc = idx*fNSamples+k;
	//cout << fStripADC[iadc] << " ";
	fStripADC[iadc] = fStripADC[iadc] + nnn;
	//if(fTime[idx]>50.)cout << k << ":" << fStripADC[iadc] << " ";
	//cout << fStripADC[iadc] << " ";
	fTotADC[idx] += nnn;
      }//if(fTime[idx]>50.)cout << endl;
      if( was_below && fTotADC[idx] > fThreshold ) {
	assert( fNOT < fNStrips );
	fOverThr[fNOT] = idx;
	++fNOT;
      }
      fStripClusters[idx].push_back(clusterID);
    }
    
    // //pedestal noise here ??? let's try
    // //printf("%d \n ", fNStrips);
    // double gain = 1.0;
    // double off = 0.0;
    // for(Int_t j = 0; j<fNStrips; j++){
    //   if(  vv->GetIdx(0) <= j && j <=  vv->GetIdx(vv->GetSize()-1) )continue;
    //   SETBIT(fType[j], kInducedStrip);
    //   fTime[j] = 0.0;
    //   fCharge[j] = 0.0;
    //   bool was_below = !( fTotADC[j] > fThreshold );
    //   for( UInt_t k=0; k<fNSamples; k++ ) {
    // 	double noise = (fRan.Gaus(30., 7.0)-off)/gain;
    // 	Int_t iadc = j*fNSamples+k;
    // 	Int_t vadc = TMath::Max(0, TMath::FloorNint(noise));
    // 	fStripADC[iadc] = vadc;
    // 	fTotADC[j] += vadc;
    //   }
    //   if( was_below && fTotADC[j] > fThreshold ) {
    // 	assert( fNOT < fNStrips );
    // 	fOverThr[fNOT] = j;
    // 	++fNOT;
    //   }
    // }//
    
    //do cross talk if requested, a big signal along the strips 
    //will induce a smaller signal as the bigger one going to the APV, 
    //the smaller signal will appear on strips that is 
    //about 32 channels away from the big signal
    if (!TSBSSimCherDigitization::fDoCrossTalk) return;
    Int_t isLeft = fRan.Uniform(1.) < 0.5 ? -1 : 1;
    Double_t factor = TSBSSimCherDigitization::fCrossFactor +
      fRan.Gaus(0., TSBSSimCherDigitization::fCrossSigma);
    if (factor <= 0.) return; //no induced signal
    
    for( Int_t j=0; j < vv->GetSize(); j++ ) {
      Int_t idx = vv->GetIdx(j);
      assert( idx >= 0 && idx < fNStrips );
      
      Int_t idxInduce = idx + isLeft*TSBSSimCherDigitization::fNCStripApart;
      if (idxInduce < 0 || idxInduce >= fNStrips ) continue; //outside the readout
      
      SETBIT(fType[idxInduce], kInducedStrip);
      //same time as the main signal strip
      fTime[idxInduce] = (fTime[idx] < vv->GetTime()) ? fTime[idx] : vv->GetTime();
      fCharge[idxInduce] += factor*vv->GetCharge(j);
      bool was_below = !( fTotADC[idxInduce] > fThreshold );
      for( UInt_t k=0; k<fNSamples; k++ ) {
	Int_t nnn = vv->GetADC(j,k);
	assert( nnn >= 0 );
	nnn *= factor;
	if( nnn == 0 ) continue;
	Int_t iadc = idxInduce*fNSamples+k;
	fStripADC[iadc] = fStripADC[iadc] + nnn;
	fTotADC[idxInduce] += nnn;
      }
      if( was_below && fTotADC[idxInduce] > fThreshold ) {
	assert( fNOT < fNStrips );
	fOverThr[fNOT] = idxInduce;
	++fNOT;
      }
    }
  }
};


UShort_t
TSBSDigitizedPlane::Threshold( Int_t thr )
{
  // Find number of strips over threshold 'thr'
  // and build index table for GetIdxOverThr.
  // This needs to be called only if one wants a change the threshold value.

  fNOT = 0;
  fThreshold = thr;

  for (UInt_t j = 0; j < fNStrips; j++)
    {
      if (fTotADC[j] > thr)
	{
	  fOverThr[fNOT] = j;
	  fNOT++;
	}
    }

  return fNOT;
};
*/

TSBSSimCherDigitization::TSBSSimCherDigitization( const TSBSSpec& spect,
						const char* name)
  : TSBSSimDetDigitization(spect, name)
{
  Init();
  Initialize (spect);
  
  fEvent = new TSBSSimEvent(5);
}

TSBSSimCherDigitization::~TSBSSimCherDigitization()
{
  DeleteObjects();
}

void
TSBSSimCherDigitization::Initialize(const TSBSSpec& spect)
{
  // Initialize digitization structures based on parameters from given
  // spectrometer

  // Avoid memory leaks in case of reinitialization
  DeleteObjects();
  
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
}

Int_t
TSBSSimCherDigitization::ReadDatabase (const TDatime& date)
{
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

void 
TSBSSimCherDigitization::ReadPMTtimesDataFiles()
{
  cout << "Reading PMT times data files " << endl;
  
  // retrieve here if there are several files... (TODO)
  TString* filenames[fNDetectors];

  Int_t index0 = 0;
  Int_t index1 = 0;
  
  for(int j = 0; j<fNDetectors; j++){
    index0 = index1;
    if(index0>0)index0++;
    index1 = fPMTtimesDataFileNames.Index(" ", index0);
    if(j==fNDetectors-1)index1 = fPMTtimesDataFileNames.Sizeof();
    filenames[j] = new TString(fPMTtimesDataFileNames(index0, index1-index0));
    
    //cout << i << "  --  " << index0 << " " << index1 << "  --  " << filenames[i]->Data() << endl;
    
    ifstream in(filenames[j]->Data());
    if(!in.is_open()){
      cout << filenames[j]->Data() << " does not exist: exit" << endl;
      return;
    }
    
    std::map< Double_t, std::pair<Double_t, double_t> > PMTtimes_detector;
    std::pair<Double_t, double_t> t_values;
    Double_t pulseheight, ne, tb, tr, tf;
    Double_t ne_prev, tr_prev, tf_prev;
    Double_t tr_calc, tf_calc;
    
    ne = tb = tr = tf = 1.0e38;
    
    in.ignore(100, ';');
    
    in >> pulseheight;
    
    fPulseHeight.push_back(pulseheight);
    
    in.ignore(100, ';');
    
    double i = 1.0;
    
    while(1){
      ne_prev = ne;
      tr_prev = tr;
      tf_prev = tf;
      
      in >> ne >> tb >> tr >> tf;
      if(!in.good())break;
      
      double ne_ = i*fPMTGain;
      
      if(ne_==ne){
	// cout  << ne_ << " ne " << ne << " tr " << tr << " tf " << tf << endl;
	
	t_values = std::make_pair(tr-tb, tf-tb);
	PMTtimes_detector.insert(std::pair<Double_t, std::pair<Double_t, Double_t> >(ne_,  t_values));
	i++;
      }else if(ne_prev<ne_ && ne_<ne){
	
	double b = (log(tr-tb)-log(tr_prev-tb))/(ne-ne_prev);
	tr_calc = (tr_prev-tb)*exp((ne_-ne_prev)*b);
	b = (log(30.0-(tf-tb))-log(30.0-(tf_prev-tb)))/(ne-ne_prev);
	tf_calc = 30.0-(30.0-(tf_prev-tb))*exp((ne_-ne_prev)*b);
	
	// cout << "ne " << ne_prev << " " << ne_ << " " << ne << endl;
	// cout << "tr " << tr_prev-tb << " " << tr_calc << " " << tr-tb << endl;
	// cout << "tf " << tf_prev-tb << " " << tf_calc << " " << tf-tb << endl;
	
	t_values = std::make_pair(tr_calc, tf_calc);
	PMTtimes_detector.insert(std::pair<Double_t, std::pair<Double_t, Double_t> >(ne_,  t_values));
	i++;
      }
      if(i>30)break;
    }
    fPulseTimes.push_back(PMTtimes_detector);
  }
}

Int_t
TSBSSimCherDigitization::Digitize (const TSBSCherData& gdata, const TSBSSpec& spect)
{
  // Digitize event after clearing all previous digitization results.

  fEvent->Clear();
  
  /*
  fSignalSector = 0;  // safe default, will normally be overridden in AdditiveDigitize

  for (UInt_t ic = 0; ic < fNChambers; ++ic) {
    for (UInt_t ip = 0; ip < fNPlanes[ic]; ++ip)
      fDP[ic][ip]->Clear();
  }
  fFilledStrips = true;
  */
  
  return AdditiveDigitize( gdata, spect );
}

Int_t
TSBSSimCherDigitization::AdditiveDigitize (const TSBSCherData& chdata, const TSBSSpec& spect)
{
  // Digitize event. Add results to any existing digitized data.
  
  UInt_t nh = chdata.GetNHit();
  
  // For signal data, determine the sector of the primary track
  bool is_background = chdata.GetSource() != 0;
  if( nh == 0 ) {
    //cout << "no hit, doing nothing " << endl;
    return 0;
  }

  UInt_t vsize = manager->GetNDetectors();
  vector<Float_t> event_time(vsize);
  vector<bool> time_set(vsize,false);
  UInt_t itime = 0;
  
  for (UInt_t ih = 0; ih < nh; ++ih) {  
    UInt_t idet = chdata.GetHitDetID (ih);
    UInt_t ipmt = chdata.GetHitPMTID (ih);
    if (ipmt >= fNPMTs[idet])
      continue;
    
    //Short_t itype = (chdata.GetParticleType(ih)==1) ? 1 : 2; // primary = 1, secondaries = 2
    if( !time_set[itime] ) {
    // Trigger time jitter, including an arbitrary offset to align signal timing
      Double_t trigger_jitter = fTrnd.Gaus(0, fTriggerJitter);
      //cout << "Offset, Jitter: " << fTriggerOffset << " " << fTriggerJitter << " => trig jitter = " << trigger_jitter << endl;
      //if( is_background ) {
	// For background data, uniformly randomize event time between
      //event_time[itime] = fTrnd.Uniform(fGateWidth) - fGateWidth/2.0 + trigger_jitter - fTriggerOffset;
	//cout << "GateWidth " << fGateWidth << ", sampling period " << fEleSamplingPeriod << endl;
      //} else {
	// Signal events occur at t = 0, smeared only by the trigger jitter
      event_time[itime] = trigger_jitter-fTriggerOffset;
      //}
    }
    time_set[itime] = true;
    
    // time_zero is defined as the beginning of the pulse being drawn out of the PMT anode
    Double_t time_zero = event_time[itime] + chdata.GetHitTime(ih) + fPMTTransitTime;
    if( is_background ) {
      time_zero+= fTrnd.Uniform(fGateWidth) - fGateWidth/2.0;
    }
    
    double totalpulsecharge = chdata.GetHitPEyield(ih)*fPMTGain;//*q_e; in electrons
    
    double t_TDC_1, t_TDC_2;
    bool TDCactive = GetTDCtimes(idet, totalpulsecharge, time_zero, t_TDC_1, t_TDC_2);
    
    //cout << t_TDC_1 << " " << t_TDC_2 << endl;
    
    if(TDCactive){
      // data flow test: it does not work that way...
      int TDCval_1 = TMath::Nint( TMath::Min( t_TDC_1/fTDCTimeConversion, pow(2, fTDCbits)-1 ) );
      int TDCval_2 = TMath::Nint( TMath::Min( t_TDC_2/fTDCTimeConversion, pow(2, fTDCbits)-1 ) );
      
      fTDCtimeArrays.at(idet).first[ipmt] = t_TDC_1;
      fTDCtimeArrays.at(idet).second[ipmt] = t_TDC_2;
      
      fTDCArrays.at(idet).first[ipmt] = TDCval_1;
      fTDCArrays.at(idet).second[ipmt] = TDCval_2;
      //fTDCArrays[idet][ipmt] = int(fabs(TDCval));
    }else{
      fTDCtimeArrays.at(idet).first[ipmt] = -1.e38;
      fTDCtimeArrays.at(idet).second[ipmt] = -1.e38;
      
      fTDCArrays.at(idet).first[ipmt] = 0;
      fTDCArrays.at(idet).second[ipmt] = 0;
    }
    
    //Short_t id = 
    SetTreeHit (ih, spect, chdata, time_zero);
  }
  return 0;
}

void
TSBSSimCherDigitization::NoDigitize (const TSBSCherData& chdata, const TSBSSpec& spect) // do not digitize event, just fill the tree
{
  //  if (!fEvCleared)  //?
  fEvent->Clear();
  UInt_t nh = chdata.GetNHit();
  
  for (UInt_t ih = 0; ih < nh; ++ih)
    {
      UInt_t idet = chdata.GetHitDetID(ih);
      UInt_t ipmt = chdata.GetHitPMTID(ih);
      
      if (ipmt >= fNPMTs[idet])
	continue;
      
      // Short_t id =
      //SetTreeHit (ih, spect, fdh, chdata, 0.0);
      SetTreeHit (ih, spect, chdata, 0.0);
    }
}

//___________________________________________________________________________________
bool 
TSBSSimCherDigitization::GetTDCtimes(int detnum, // detctor ID
				     double C,   // total pulse charge
				     double time_zero, //time of the beginning of the pulse
				     double& t1, // rising time 
				     double& t2  // falling time
				     )
{
  //dumm... test data flow
  if(C*fPulseHeight[detnum]<fTDCthreshold){
    cout << "Pulse amplitude = " << C*fPulseHeight[detnum]
	 << " V below threshold = " << fTDCthreshold << " V." << endl;
    // 	 << " PMT signal max voltage " << C*fReadOutImpedance/(fPMTFWHM*1.0e-9) << endl;
    return false;// return false if the pulse does not cross the set threshold
  }else{
    std::map< Double_t, std::pair<Double_t, Double_t> >::iterator it = fPulseTimes[detnum].find(C);
    std::pair<Double_t, Double_t> time_pair = it->second;
    t1 = time_zero+time_pair.first;
    t2 = time_zero+time_pair.second;
    
    return true;
  }
}


//___________________________________________________________________________________
void
TSBSSimCherDigitization::Print() const
{
  /*
  cout << "GEM digitization:" << endl;
  cout << "  Gas parameters:" << endl;
  cout << "    Gas ion width: " << fGasWion << endl;
  cout << "    Gas diffusion: " << fGasDiffusion << endl;
  cout << "    Gas drift velocity: " << fGasDriftVelocity << endl;
  cout << "    Avalanche fiducial band: " << fAvalancheFiducialBand << endl;
  cout << "    Avalanche charge statistics: " << fAvalancheChargeStatistics << endl;
  cout << "    Gain mean: " << fGainMean << endl;
  cout << "    Gain 0: " << fGain0 << endl;

  cout << "  Electronics parameters:" << endl;
  cout << "    Trigger offset: " << fTriggerOffset << endl;
  cout << "    Trigger jitter: " << fTriggerJitter << endl;
  cout << "    Sampling Period: " << fEleSamplingPeriod << endl;
  cout << "    Sampling Points: " << fEleSamplingPoints   << endl;
  cout << "    Pulse Noise width: " << fPulseNoiseSigma << endl;
  cout << "    ADC offset: " << fADCoffset << endl;
  cout << "    ADC gain: " << fADCgain << endl;
  cout << "    ADC bits: " << fADCbits << endl;
  cout << "    Gate width: " << fGateWidth << endl;

  cout << "  Pulse shaping parameters:" << endl;
  cout << "    Pulse shape tau0: " << fPulseShapeTau0 << endl;
  cout << "    Pulse shape tau1: " << fPulseShapeTau1 << endl;
  */
}

// Tree methods
void
TSBSSimCherDigitization::InitTree (const TSBSSpec& spect, const TString& ofile)
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
TSBSSimCherDigitization::SetTreeEvent (const TSBSCherData& tscd,
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
TSBSSimCherDigitization::SetTreeHit (const UInt_t ih,
				     const TSBSSpec& spect,
				     //TSolGEMVStrip* const *dh,
				     const TSBSCherData& tscd,
				     Double_t t0 )
{
  // Sets the variables in fEvent->fGEMClust describing a hit
  // This is later used to fill the tree.
  
  TSBSSimEvent::PMTHit hit;
  
  UInt_t idet = tscd.GetHitDetID(ih);
  UInt_t ipmt = tscd.GetHitPMTID(ih);
  
  hit.fID = fEvent->fPMTHits.size();//ih;
  hit.fSource = tscd.GetSource();
  hit.fType = (tscd.GetParticleType(ih)==1) ? 1 : 2;// primary = 1, secondaries = 2
  hit.fMCtrackID = tscd.GetParticleType(ih);
  hit.fMCtrackPID = tscd.GetMCtrackPID(ih);
  hit.fOrigVolFlag = tscd.GetOriginVolFlag(ih);
  hit.fXPMT = tscd.GetHitXPMT(ih);
  hit.fYPMT = tscd.GetHitYPMT(ih);
  hit.fNpe = tscd.GetHitPEyield(ih);
  hit.fTime = t0;
  hit.fTDCtime[0] = fTDCtimeArrays.at(idet).first[ipmt];
  hit.fTDCtime[1] = fTDCtimeArrays.at(idet).second[ipmt];
  // Digitization results for this hit
  hit.fDetID = idet;
  hit.fChannel = ipmt;
  hit.fPMTrow = round((hit.fXPMT+manager->GetPMTmatrixVext(hit.fDetID)/2.0)/manager->GetPMTdistX(hit.fDetID));
  double PMTcol_ = (hit.fYPMT+manager->GetPMTmatrixHext(hit.fDetID)/2.0)/manager->GetPMTdistY(hit.fDetID); 
  if( fabs(PMTcol_ - round(PMTcol_))<1.e-3 ){
    hit.fPMTcol = round(PMTcol_);
  }else{
    hit.fPMTcol = int(PMTcol_);
  }
  //hit.fTDC[0] = fTDCArrays.at(idet).first[ipmt];
  //hit.fTDC[1] = fTDCArrays.at(idet).second[ipmt];
  
  /*
  Bool_t newclus = true;
  //fill the 
  for(int i_ = 0; i_<fEvent->fMCClusterHitID.size(); i_++){
    if(tscd.GetParticleType(ih)==fEvent->fMCClusterHitID[i_].first){
      fEvent->fMCClusterHitID[i_].second.push_back(ih);
      newclus = false;
      break;
    }
  }
  if(newclus){
    Int_t trkID = tscd.GetParticleType(ih);
    std::vector<Short_t> ClusterPMTlist;
    ClusterPMTlist.push_back(ih);
    fEvent->fMCClusterHitID.push_back(make_pair(trkID, ClusterPMTlist));
  }
  */
  
  // cout << "PMT " << ipmt 
  //      << " tdc0 " << fTDCArrays.at(idet).first[ipmt] 
  //      << " tdc1 " << fTDCArrays.at(idet).second[ipmt] << endl;  
  
  div_t d = div( ipmt, manager->GetChanPerSlot() );
  hit.fVETROCID = d.quot;
  int i_chan = d.rem;
  
  uint32_t TDCvetrocWord0, TDCvetrocWord1;
  TDCvetrocWord0 = TDCvetrocWord1 = 0;
  
  bool header[8] = {0, 0, 0, 0, 0, 0, 1, 1};
  bool channel[8];
  bool tdc[16];
  
  //cout << "channel " << i_chan << " time_r " << fTDCArrays.at(idet).first[ipmt] << endl;
  //cout << "Vetroc words: init: " << TDCvetrocWord0 << " " << TDCvetrocWord1 << endl;
  for(int i = 0; i<16; i++){
    if(i<8){
      TDCvetrocWord0 ^= (-header[i] ^ TDCvetrocWord0) & (1 << (i+24));
      TDCvetrocWord1 ^= (-header[i] ^ TDCvetrocWord1) & (1 << (i+24));
      channel[i] = (i_chan >> i) & 1;
      TDCvetrocWord0 ^= (-channel[i] ^ TDCvetrocWord0) & (1 << (i+16));
      TDCvetrocWord1 ^= (-channel[i] ^ TDCvetrocWord1) & (1 << (i+16));
      
      // cout << "header bit " << i << " = " << header[i] 
      // 	   << ", VETROC word bit " << (i+24) << " = " << ((TDCvetrocWord0 >> (i+24)) & 1UL) << endl;
      // cout << "channel bit " << i << " = " << channel[i] 
      // 	   << ", VETROC word bit " << (i+16) << " = " << ((TDCvetrocWord0 >> (i+16)) & 1UL) << endl;
    }
    tdc[i] = (fTDCArrays.at(idet).first[ipmt] >> i) & 1;
    TDCvetrocWord0 ^= (-tdc[i] ^ TDCvetrocWord0) & (1UL << i);
    // cout << "tdc_r bit " << i << " = " << tdc[i] 
    // 	 << ", VETROC word bit " << i << " = " << ((TDCvetrocWord0 >> i) & 1UL) << endl;
    tdc[i] = (fTDCArrays.at(idet).second[ipmt] >> i) & 1;
    TDCvetrocWord1 ^= (-tdc[i] ^ TDCvetrocWord1) & (1UL << i);
  }
  
  TDCvetrocWord1 ^= (-1 ^ TDCvetrocWord1) & (1UL << 26);
  //cout << "Vetroc words: test: " << TDCvetrocWord0 << " " << TDCvetrocWord1 << endl;
  
  hit.fTDC[0] = TDCvetrocWord0;
  hit.fTDC[1] = TDCvetrocWord1;
  
  fEvent->fPMTHits.push_back( hit );
  
  if( hit.fType == 1 && hit.fSource == 0 )
    fEvent->fNSignal++;
  
  return hit.fID;
}

/*
void
TSBSSimCherDigitization::CleanClusterList ()
{
  //if(fEvent->fMCClusterHitID.size()>1)cout << "Cluster list size = " << fEvent->fMCClusterHitID.size() << endl;
  for(int i_ = fEvent->fMCClusterHitID.size()-1; i_>=0; i_--){
    // if(fEvent->fMCClusterHitID.size()>1){
    //   cout << " i_ = " << i_ << "; MC track ID " << fEvent->fMCClusterHitID[i_].first 
    // 	   << " cluster size " << fEvent->fMCClusterHitID[i_].second.size() << endl;
    //   for(int j_ = 0; j_<fEvent->fMCClusterHitID[i_].second.size(); j_++){
    // 	cout << " j_ = " << " PMT hit ID " << fEvent->fMCClusterHitID[i_].second[j_] << endl;
    //   }
    // }
    if(fEvent->fMCClusterHitID[i_].second.size()<2){
      fEvent->fMCClusterHitID.erase(fEvent->fMCClusterHitID.begin()+i_);
    }
  }
}
*/

void
TSBSSimCherDigitization::FillTree ()
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
TSBSSimCherDigitization::WriteTree () const
{
  //cout << "write tree " << fOFile << " " << fOTree << endl;
  
  if (fOFile && fOTree) {
    fOFile->cd();
    fOTree->Write();
  }
}

void
TSBSSimCherDigitization::CloseTree () const
{
  if (fOFile) fOFile->Close();
}

