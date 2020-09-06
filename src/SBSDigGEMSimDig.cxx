#include "SBSDigGEMSimDig.h"
#include "SBSDigGEMDet.h"
#include "gmn_tree.h"

//gas parameters
#define fGasWion 26             // eV
#define fGasDiffusion 1.e5       // mm2/s
#define fGasDriftVelocity 5.5e7   // mm/s
#define fAvalancheFiducialBand 10. // number of sigma defining the band around the avalanche in readout plane
#define fAvalancheChargeStatistics 1 // 0 Furry, 1 Gaussian
#define fGainMean 8.e3
#define fGain0 20.
#define fMaxNIon 1.e4               //maximum amount of ion pairs allowed in the digitization
  
#define fSNormNsigma 18.          //fSNormNsigma is an arbitrary multiplicative fact  
#define fAvaGain 20.
#define fLateralUncertainty 0.

//electronics parameters
#define fAPVTimeJitter 25.    // time jitter associated with the APV internal clock
  
#define fEleSamplingPoints 6
#define fEleSamplingPeriod 25. // ns
#define fADCoffset 0.         // ADC offset
#define fADCgain 1.          // ADC gain
#define fADCbits 12         // ADC resolutions in bits
#define fGateWidth 400.    // to be changed , ns - pulse shape width at ~1/10 max

//parameter for GEM pedestal noise
#define fPulseNoiseSigma 20.   // additional sigma term of the pedestal noise
#define fPulseNoisePeriod 200. // period of the pedestal noise, assuming sinusoidal function
#define fPulseNoiseAmpConst 0  // constant term of the pedestal noise amplitude
#define fPulseNoiseAmpSigma 0  // sigma term of the pedestal noise amplitude

//paramters for cross-talk
#define fNCStripApart 32  // # of strips the induced signal is away from the mean signal
#define fCrossFactor 0.1  //reduction factor for the induced signal
#define fCrossSigma 0.03  //uncertainty of the reduction factor

// Pulse shaping parameters
#define fPulseShapeTau 56.   // [ns] GEM model 0 = 50. in SiD model
//#define fPulseShapeTau1 0.   // [ns] GEM model only; if negative assume SiD model

#define fEntranceRef -1.5  // z position of the copper layer right before the first GEM gas layer,             // relative to the center of the GEM chamber
                         // which introduce additional uncertainty in the lateral direction
#define fRoutZ 9.185   // z-distance hit entrance to readout plane [mm]

//numerical integration parameters
#define fYIntegralStepsPerPitch 4
#define fXIntegralStepsPerPitch 4

#define fNROPlanes 2
#define fStripPitch 4.e-4

using namespace std;

//stupid but: let's do the defautl constructor:
SBSDigGEMSimDig::SBSDigGEMSimDig()// :
  //fGasWion(), fGasDiffusion(), fGasDriftVelocity(), fAvalancheFiducialBand(), fAvalancheChargeStatistics(), fGainMean(), fGain0(), fMaxNIon(), fSNormNsigma(), fAvaGain(), fAPVTimeJitter() 
{
}

SBSDigGEMSimDig::SBSDigGEMSimDig(int nchambers, double* trigoffset, double zsup_thr, int napv, double* commonmode_array) : fZeroSup(zsup_thr) 
{
  for(int i = 0; i<nchambers; i++)fTriggerOffset.push_back(trigoffset[i]);
  if(fZeroSup>0)fDoZeroSup = true;
  if(napv)fDoCommonMode = false;
  for(int i = 0; i<napv; i++)fCommonModeArray.push_back(commonmode_array[i]);
  fRIon.resize((int)fMaxNIon);
}


SBSDigGEMSimDig::~SBSDigGEMSimDig()
{
}


//.......................................................
// ionization Model
//
void
SBSDigGEMSimDig::IonModel(TRandom3* R,
			  const TVector3& xi,
			  const TVector3& xo,
			  const Double_t elost ) // eV
{
#define DBG_ION 0

  TVector3 vseg = xo-xi; // mm
  
  // ---- extract primary ions from Poisson
  fRNIon = R->Poisson(elost/fGasWion);

  if (fRNIon <=0)
    return;

#if DBG_ION > 0
  cout << "E lost = " << elost << ", " << fRNIon << " ions" << endl;
#endif
  if (fRNIon > fMaxNIon) {
#if DBG_ION > 0
    cout << __FUNCTION__ << ": WARNING: too many primary ions " << fRNIon << " limit to "
	 << fMaxNIon << endl;
#endif
    fRNIon = fMaxNIon;
  }

  fRSMax = 0.;
  fRTotalCharge = 0;
  fRTime0 = 999999.; // minimum time of drift

  for (UInt_t i=0;i<fRNIon;i++) { // first loop used to evaluate quantities
    IonPar_t ip;

    Double_t lion = R->Uniform(0.,1.); // position of the hit along the track segment (fraction)

    //In principle, the lateral uncertainty should have been put in the Ava model, but not here
    //But since we are not simulating the details of the avalanche, I think it is ok (Weizhi)
    ip.X = vseg.X()*lion+xi.X() + R->Gaus(0., fLateralUncertainty);
    ip.Y = vseg.Y()*lion+xi.Y() + R->Gaus(0., fLateralUncertainty);

    // Note the definition of fRoutZ is the distance from xi.Z() to xrout.Z():
    //        xi               xo   xrout
    // |<-LD->|<-----vseg----->|    |
    // |<-------fRoutZ---------|--->|
    // |      |<-lion*vseg->   |    |
    // |      |             <--LL-->|
    
    Double_t LD = TMath::Abs(xi.Z() - fEntranceRef);//usually should be 0,
                                            //unless particle is produced inside the gas layer

    Double_t LL = TMath::Abs(fRoutZ - LD - vseg.Z()*lion);
    Double_t ttime = LL/fGasDriftVelocity; // traveling time from the drift gap to the readout
    
    //cout << " rout Z  (mm?) " << fRoutZ << ", LD (mm?) " << LD << " vseg Z (mm?) " << vseg.Z()  << endl;
    //cout << " travelling length (mm?) " << LL << ", travelling time:  " <<  ttime << endl;
    
    fRTime0 = TMath::Min(ttime, fRTime0); // minimum traveling time [s]

    ip.SNorm = TMath::Sqrt(2.*fGasDiffusion*ttime); // spatial sigma on readout [mm]
    // cout<<"ip.SNorm: "<<ip.SNorm<<endl;
    //  cout<<"vseg: "<<vseg.X()<<" deltaX: "<<vseg.X()*lion<<endl;
    if( fAvalancheChargeStatistics == 1 ) {
      Double_t gnorm = fGainMean/TMath::Sqrt(fGain0); // overall gain TBC
      ip.Charge = R->Gaus(fGainMean, gnorm); // Gaussian distribution of the charge
    }
    else {
      ip.Charge = R->Exp(fGainMean); // Furry distribution
    }

    if( ip.Charge > 0 )
      fRTotalCharge += ip.Charge;
    else
      ip.Charge = 0;

    fRSMax = TMath::Max(ip.SNorm, fRSMax);

    // Derived quantities needed by the numerical integration in AvaModel
    ip.SNorm *= fSNormNsigma;
    ip.R2 = ip.SNorm * ip.SNorm;
    ip.ggnorm = ip.Charge * TMath::InvPi() / ip.R2; // normalized charge

#if DBG_ION > 1
    printf("z coords %f %f %f %f lion %f LL %lf\n",
	   xi.Z(), xo.Z(), vseg.Z(), lion, LL);
    printf("ttime = %e\n", ttime);
#endif
#if DBG_ION > 0
    cout << " x, y = " << ip.X << ", " << ip.Y << " snorm = "
	 << ip.SNorm/fSNormNsigma << " charge " << ip.Charge << endl;
    cout << "fRTime0 = " << fRTime0 << endl;
    cout << "fRion size " << fRIon.size() << " " << i << endl;
#endif
    
    fRIon[i] = ip;
  }
  return;
}

Short_t 
SBSDigGEMSimDig::ADCConvert(Double_t val, Double_t off, Double_t gain, Int_t bits)
{
  // Convert analog value 'val' to integer ADC reading with 'bits' resolution
  assert( bits >= 0 && bits <= 12 );

  if( val < 0. )
    val = 0.;
  Double_t vvv = (val - off)/gain;
  //std::cout<<val<<" : "<<vvv<<std::endl;
  //printf("offset = %1.3f, gain = %1.3f, input value = %1.3f , output value = %1.3f  \n", off, gain, val, vvv);
  Double_t saturation = static_cast<Double_t>( (1<<bits)-1 );
  if( vvv > saturation )
    vvv = saturation;

  Short_t dval =
    static_cast<Short_t>( TMath::Floor( (vvv>saturation) ? saturation : vvv ));

  //  cerr << val << " dval = " << dval << endl;
  if( dval < 0 ) dval = 0;
  return dval;

}

// Pulse Shape SiD model
// APV25 time function from  M. Friedl et al NIMA 572 (2007) pg 385-387 (APV25 for silicon!)
//

Double_t 
SBSDigGEMSimDig::PulseShape(Double_t t, 
			    Double_t C,  // normalization factor
			    Double_t Tp) // shaping time 
{

  Double_t v;
  Double_t x;
  x = t/Tp;
  v = C/Tp * x * TMath::Exp(-x);
  
  return ( v>0. ) ? v : 0.;

}



//.......................................................
// avalanche model
//

//TGEMSBSGEMHit **
void SBSDigGEMSimDig::AvaModel(const int ic,
			       SBSDigGEMDet* gemdet, 
			       TRandom3* R,
			       const TVector3& xi,
			       const TVector3& xo,
			       const Double_t t0)
{
#define DBG_AVA 0
#if DBG_AVA > 0
  cout << "Chamber " << ic << "----------------------------------" << endl;
  cout << "In  " << xi.X() << " " << xi.Y() << " " << xi.Z() << endl;
  cout << "Out " << xo.X() << " " << xo.Y() << " " << xo.Z() << endl;
#endif

  // xi, xo are in chamber frame, in mm

  Double_t nsigma = fAvalancheFiducialBand; // coverage factor

#if DBG_AVA > 0
  cout << "fRSMax, nsigma " << fRSMax << " " << nsigma << endl;
#endif

  Double_t x0,y0,x1,y1; // lower and upper corners of avalanche diffusion area
  
  if (xi.X()<xo.X()) {
    x0 = xi.X()-nsigma*fRSMax;
    x1 = xo.X()+nsigma*fRSMax;
  } else {
    x1 = xi.X()+nsigma*fRSMax;
    x0 = xo.X()-nsigma*fRSMax;
  }

  if (xi.Y()< xo.Y()) {
    y0 = xi.Y()-nsigma*fRSMax;
    y1 = xo.Y()+nsigma*fRSMax;
  } else {
    y1 = xi.Y()+nsigma*fRSMax;
    y0 = xo.Y()-nsigma*fRSMax;
  }
  
  // Check if any part of the avalanche region is in the active area of the sector.
  // Here, "active area" means the chamber's *bounding box*, which is
  // larger than the wedge's active area (section of a ring)
  
  //const TGEMSBSGEMChamber& chamber = spect.GetChamber(ic);
  Double_t glx = (-gemdet->GEMPlanes[ic*2].dX())/2.*1.e3;//+fStripPitch)/2.0*1000.;//(chamber.GetPlane(0).GetStripLowerEdge(0)+chamber.GetPlane(0).GetSPitch()/2.0) * 1000.0;
  Double_t gly = (-gemdet->GEMPlanes[ic*2+1].dX())/2.*1.e3;//+fStripPitch)/2.0*1000.;
  Double_t gux = (gemdet->GEMPlanes[ic*2].dX())/2.*1.e3;//+fStripPitch)/2.0*1000.;//(chamber.GetPlane(0).GetStripUpperEdge(chamber.GetPlane(0).GetNStrips()-1) -chamber.GetPlane(0).GetSPitch()/2.0) * 1000.0;
  Double_t guy = (gemdet->GEMPlanes[ic*2+1].dX())/2.*1.e3;//+fStripPitch)/2.0*1000.;//(chamber.GetPlane(1).GetStripUpperEdge(chamber.GetPlane(1).GetNStrips()-1) -chamber.GetPlane(1).GetSPitch()/2.0) * 1000.0;
  
  if (x1<glx || x0>gux ||
      y1<gly || y0>guy) { // out of the sector's bounding box
    cerr << __FILE__ << " " << __FUNCTION__ << ": out of sector, "
	 << "chamber " << ic << " sector " << ic/30 << " plane " << ic%30 << endl
	 << "Following relations should hold:" << endl
	 << "(x1 " << x1 << ">glx " << glx << ") (x0 " << x0 << "<gux " << gux << ")" << endl
	 << "(y1 " << y1 << ">gly " << gly << ") (y0 " << y0 << "<guy " << guy << ")" << endl;
    //return 0;
  }
  
  //bool bb_clipped = (x0<glx||y0<gly||x1>gux||y1>guy);
  if(x0<glx) x0=glx;
  if(y0<gly) y0=gly;
  if(x1>gux) x1=gux;
  if(y1>guy) y1=guy;

  // Loop over chamber planes
  double roangle, dx, xoffset;
  int nstrips;
  
  double xt_factor;
  int isLeft;
  //TGEMSBSGEMHit **virs;
  //virs = new TGEMSBSGEMHit *[fNROPlanes[ic]];
  for (UInt_t ipl = 0; ipl < fNROPlanes; ++ipl){
#if DBG_AVA > 0
     cout << "coordinate " << ipl << " =========================" << endl;
#endif

    xt_factor = fCrossFactor+R->Gaus(fCrossSigma);
    isLeft = R->Uniform(1.) < 0.5 ? -1 : 1;
    
    // Compute strips affected by the avalanche
    //const SBSDigGEMPlane& pl = gemdet->GEMPlanes[ic*2+ipl];
    roangle = gemdet->GEMPlanes[ic*2+ipl].ROangle();
    dx = gemdet->GEMPlanes[ic*2+ipl].dX();
    //xoffset = gemdet->GEMPlanes[ic*2+ipl].Xoffset();
    nstrips = gemdet->GEMPlanes[ic*2+ipl].GetNStrips();
    
    // Positions in strip frame
    Double_t xs0 = x0*cos(roangle) - y0*sin(roangle);
    Double_t ys0 = x0*sin(roangle) + y0*cos(roangle);
    Double_t xs1 = x1*cos(roangle) - y1*sin(roangle); 
    Double_t ys1 = y1*sin(roangle) + y1*cos(roangle);
#if DBG_AVA > 0
     cout << "glx gly gux guy " << glx << " " << gly << " " << gux << " " << guy << endl;
     cout << "xs0 ys0 xs1 ys1 " << xs0 << " " << ys0 << " " << xs1 << " " << ys1 << endl;
#endif

    Int_t iL = max(0, Int_t((xs0*1.e-3+dx/2.)/fStripPitch) );
    iL = min(iL, nstrips);
    //pl.GetStrip (xs0 * 1e-3, ys0 * 1e-3);
    Int_t iU = min(Int_t((xs1*1.e-3+dx/2.)/fStripPitch), nstrips);
    iU = max(0, iU);
    //pl.GetStrip (xs1 * 1e-3, ys1 * 1e-3);
     
    // Check for (part of) the avalanche area being outside of the strip region
    //if( (iL <= 0 && iU <= 0) || (iL>=pl.GetNStrips() && iU>=pl.GetNStrips()) ) {
      // All of the avalanche outside -> nothing to do
      // TODO: what if this happens for only one strip coordinate (ipl)?
// #if DBG_AVA > 0
//       cerr << __FILE__ << " " << __FUNCTION__ << ": out of active area, "
// 	   << "chamber " << ic << " sector " << ic%30 << " plane " << ic/30 << endl
// 	   << "iL_raw " << pl.GetStripUnchecked(xs0*1e-3) << " "
// 	   << "iU_raw " << pl.GetStripUnchecked(xs1*1e-3) << endl
// 	   << endl << endl;
// #endif
    if(iL==iU){//nothing to do
      return;
    }

    if(iU<iL)swap(iU, iL);


    //
    // Bounds of rectangular avalanche region, in strip frame
    //

    // Limits in x are low edge of first strip to high edge of last
    Double_t xl = (iL*fStripPitch-dx/2.)*1.e3;//pl.GetStripLowerEdge (iL) * 1000.0;
    Double_t xr = (iU*fStripPitch-dx/2.)*1.e3;//pl.GetStripUpperEdge (iU) * 1000.0;

#if DBG_AVA > 0
    cout << "iL gsle " << iL << " " << xl << endl;
    cout << "iU gsue " << iU << " " << xr << endl;
#endif
    
    // Limits in y are y limits of track plus some reasonable margin
    // We do this in units of strip pitch for convenience (even though
    // this is the direction orthogonal to the pitch direction)

    // Use y-integration step size of 1/10 of strip pitch (in mm)
    Double_t yq = fStripPitch * 1000.0 / fYIntegralStepsPerPitch;
    Double_t yb = ys0, yt = ys1;
    if (yb > yt)
      swap( yb, yt );
    yb = yq * TMath::Floor (yb / yq);
    yt = yq * TMath::Ceil  (yt / yq);

    // We should also allow x to have variable bin size based on the db
    // the new avalanche model (Cauchy-Lorentz) has a very sharp full width
    // half maximum, so if the bin size is too large, it can introduce
    // fairly large error on the charge deposition. Setting fXIntegralStepsPerPitch
    // to 1 will go back to the original version -- Weizhi Xiong

    Int_t nstrips = iU - iL + 1;
    Int_t nx = (iU - iL + 1) * fXIntegralStepsPerPitch;
    Int_t ny = TMath::Nint( (yt - yb)/yq );
#if DBG_AVA > 0
    cout << "xr xl yt yb nx ny "
	 << xr << " " << xl << " " << yt << " " << yb
	 << " " << nx << " " << ny << endl;
#endif
    assert( nx > 0 && ny > 0 );

    // define function, gaussian and sum of gaussian

    Double_t xbw = (xr - xl) / nx;
    Double_t ybw = (yt - yb) / ny;
#if DBG_AVA > 0
    cout << "xbw ybw " << xbw << " " << ybw << endl;
#endif
    
    Int_t sumASize = nx * ny;
#if DBG_AVA > 0
    cout<<nx<<" : "<<ny<< ", nx*ny " << sumASize <<" Nstrips: "<<nstrips<<endl;
#endif
    fSumA.resize(sumASize);
    memset (&fSumA[0], 0, fSumA.size() * sizeof (Double_t));
    //Double_t fSumA[sumASize];
    //memset (fSumA, 0, sumASize * sizeof (Double_t));
    //for(Int_t i=0; i<sumASize; i++){
    //fSumA[i] = 0;
    //}
    
    //    fSumA.resize(nx*ny);
    //memset (&fSumA[0], 0, fSumA.size() * sizeof (Double_t));
#if DBG_AVA > 0
    cout << fRNIon << " " << fRIon.size() << endl;
#endif
    
    for (UInt_t i = 0; i < fRNIon; i++){
      Double_t frxs = fRIon[i].X*cos(roangle) - fRIon[i].Y*sin(roangle);
      Double_t frys = fRIon[i].X*sin(roangle) + fRIon[i].Y*cos(roangle);
      //pl.PlaneToStrip (frxs, frys);
      //frxs *= 1e3; frys *= 1e3;
      //  cout<<"IonStrip: "<<pl.GetStrip(frxs*1e-3,frys*1e-3)<<endl;
      // bin containing center and # bins each side to process
      Int_t ix = (frxs-xl) / xbw;
      Int_t iy = (frys-yb) / ybw;
      Int_t dx = fRIon[i].SNorm / xbw  + 1;
      Int_t dy = fRIon[i].SNorm / ybw  + 1;
#if DBG_AVA > 1
       cout << "ix dx iy dy " << ix << " " << dx << " " << iy << " " << dy << endl;
#endif

      //
      // NL change:
      //
      // ggnorm is the avalance charge for the i^th ion, and R2 is the square of the radius of the diffusion 
      // circle, mutiplied by the kSNormNsigma factor: (ip.SNorm * ip.SNorm)*kSNormNsigma*kSNormNsigma. All 
      // strips falling within this circle are considered in charge summing. 
      //
      // The charge contribution to a given strip by the i^th ion is evaluated by a Lorentzian (or Gaussian)
      // distribution; the sigma for this distribution is eff_sigma, which is the actual avalance sigma. 
      //
      Double_t ggnorm = fRIon[i].ggnorm;
      Double_t r2 = fRIon[i].R2;
      Double_t eff_sigma_square = r2/(fSNormNsigma*fSNormNsigma);
      Double_t eff_sigma = TMath::Sqrt(eff_sigma_square);
      Double_t current_ion_amplitude = fAvaGain*ggnorm*(1./(TMath::Pi()*eff_sigma))*(eff_sigma*eff_sigma);
      
      
      // xc and yc are center of current bin
      
      // Loop over bins
      Int_t min_binNb_x = max(ix-dx,0);
      Int_t min_binNb_y = max(iy-dy,0);
      Int_t max_binNb_x = min(ix+dx+1,nx);
      Int_t max_binNb_y = min(iy+dy+1,ny);
      Int_t jx = min_binNb_x;
      Double_t xc = xl + (jx+0.5) * xbw;
      for (; jx < max_binNb_x; ++jx, xc += xbw){
	Double_t xd2 = frxs-xc; xd2 *= xd2;
	//if( xd2 > r2 ){
	//  if( (xc - frxs)>0 )
	//    break;
	//  else
	//    continue;
	//}
	Int_t jx_base = jx * ny;
	Int_t jy = min_binNb_y;
	Double_t yc = yb + (jy+0.5) * ybw;
	
	for (; jy < max_binNb_y; ++jy, yc += ybw){
	  Double_t yd2 = frys-yc; yd2 *= yd2;
	  //if( yd2 > r2 ){
	  //  if( (yc - frys)>0 )
	  //    break;
	  //  else
	  //    continue;
	  // }
	  if( xd2 + yd2 <= r2 ) {
	    //cout << current_ion_amplitude << " " << xd2+yd2 << " " << eff_sigma_square << endl;
	    fSumA[jx_base+jy] += current_ion_amplitude / ((xd2+yd2)+eff_sigma_square);
	  }
	}//cout<<endl;
      }//cout<<"##########################################################################"<<endl<<endl;getchar();
    }
    
#if DBG_AVA > 0
    cout << "t0 = " << t0 << " plane " << ipl 
	 << endl;
#endif

    //virs[ipl] = new TGEMSBSGEMHit(nx,fEleSamplingPoints);
    //virs[ipl]->SetTime(t0);
    //virs[ipl]->SetHitCharge(fRTotalCharge);
    
    Int_t ai=0;
    Double_t area = xbw * ybw;

//when we integrate in order to get the signal pulse, we want all charge
    //deposition on the area of a single strip -- Weizhi
    
    //cout << "number of strips: " << nstrips << ", number of samples " << fEleSamplingPoints << " area: " << area << endl;
    
    // if(nstrips>0){cout<<"nstrips: "<<nstrips<<" Nion: "<<fRNIon<<endl;}
    
    for (Int_t j = 0; j < nstrips; j++){
      //  cout<<"strip: "<<iL+j<<":    ";
      Int_t posflag = 0;
      Double_t us = 0.;
      for (UInt_t k=0; k<fXIntegralStepsPerPitch; k++){
	Double_t integralY_tmp = 0;
	int kx = (j * fXIntegralStepsPerPitch + k) * ny;
	for( Int_t jy = ny; jy != 0; --jy )
	  integralY_tmp += fSumA[kx++];
	
	us += integralY_tmp * area;
	//	us += IntegralY( fSumA, j * fXIntegralStepsPerPitch + k, nx, ny ) * area;
	//if(us>0)cout << "k " << k << ", us " << us << endl;
      }

#if DBG_AVA > 0
      cout << "strip " << j << " us " << us << endl;
#endif

      // cout <<setw(6)<< (Int_t)(us/100);
      //  cout<<iL+j<<" : "<<us<<endl;
      //generate the random pedestal phase and amplitude
      // Double_t phase = fTrnd.Uniform(0., fPulseNoisePeriod);
      // Double_t amp = fPulseNoiseAmpConst + fTrnd.Gaus(0., fPulseNoiseAmpSigma);
      
      for (Int_t b = 0; b < fEleSamplingPoints; b++){
	Double_t pulse = PulseShape (fEleSamplingPeriod * b - t0,
				     us,
				     fPulseShapeTau);
	//fPulseShapeTau0, fPulseShapeTau1 );
	
	Short_t dadc = ADCConvert( pulse,
				   0,// fADCoffset,
				   fADCgain,
				   fADCbits );
#if DBG_AVA > 0
	if(pulse>0)
	  cout << "strip number " << j << ", sampling number " << b << ", t0 = " << t0 << endl
	       << "pulse = " << pulse << ", (val - off)/gain = " 
	       << (pulse-fADCoffset)/fADCgain << ", dadc = " << dadc << endl;
#endif
	//fDADC[b] = dadc;
	gemdet->GEMPlanes[ic*2+ipl].AddADC(j, b, dadc);
	//posflag += dadc;
	//if(dadc>0)cout << t0 << " " << pulse << " " << dadc << endl;
	//cross talk here ?
	if(xt_factor>0){
	  if(j+isLeft*fNCStripApart>=0 && j+isLeft*fNCStripApart<nstrips){
	    gemdet->GEMPlanes[ic*2+ipl].AddADC(j+isLeft*fNCStripApart, b, dadc*xt_factor);
	  }
	  
	}
      }//cout <<"  "<<t0<< endl;

      

      //if (posflag > 0) { // store only strip with signal
      //for (Int_t b = 0; b < fEleSamplingPoints; b++)
	  //{virs[ipl]->AddSampleAt (fDADC[b], b, ai);}
	  //virs[ipl]->AddStripAt (iL+j, ai);
	  //virs[ipl]->AddChargeAt (us, ai);
      //ai++;
      //}

      //  cout<<endl;
    }//getchar();
    
    //cout << "number of strips with signal " << ai << endl;

    //virs[ipl]->SetSize(ai);
  }
  
  //return virs;
}


Int_t
SBSDigGEMSimDig::Digitize (SBSDigGEMDet* gemdet,
			   TRandom3* R)
{
  // Digitize event. Add results to any existing digitized data.

  //UInt_t nh = gdata.GetNHit();
  bool is_background = false;
  Float_t event_time=0,time_zero=0;
  Double_t trigger_jitter = R->Uniform(-fAPVTimeJitter/2, fAPVTimeJitter/2);
  
  // For signal data, determine the sector of the primary track
  
  for(int ih = 0; ih<gemdet->fGEMhits.size(); ih++){
    is_background = (gemdet->fGEMhits[ih].source==0);
    UInt_t igem = gemdet->fGEMhits[ih].module;
    //UInt_t igem = iplane/2;
    
    if(igem>=16)cout << igem << endl;
    //cout<<igem<<":"<<imodule<<":"<<iplane<<endl;
    //Short_t itype = (gdata.GetParticleType(ih)==1) ? 1 : 2; // primary = 1, secondaries = 2
    // if(gdata.GetParticleType(ih)!=1){cout<<"x"<<endl;getchar();}

    // These vectors are in the spec frame, we need them in the chamber frame
    TVector3 vv1(gemdet->fGEMhits[ih].xin*1.e3,
		 gemdet->fGEMhits[ih].yin*1.e3, 
		 gemdet->fGEMhits[ih].zin*1.e3);
    TVector3 vv2(gemdet->fGEMhits[ih].xout*1.e3,
		 gemdet->fGEMhits[ih].yout*1.e3, 
		 gemdet->fGEMhits[ih].zout*1.e3);
    
    if(abs(vv1.X()-vv2.X())>50 || abs(vv1.Y()-vv2.Y())>50){//in mm
      //cout<<abs(vv1.X()-vv2.X())<<endl;
      //getchar();
      continue;
    }
    IonModel (R, vv1, vv2, gemdet->fGEMhits[ih].edep );
    
    // Get Signal Start Time 'time_zero'
    if( is_background ) {
      // For background data, uniformly randomize event time between
      // -fGateWidth to +75 ns (assuming 3 useful 25 ns samples).
      // Not using HitTime from simulation file but randomize HitTime to cycle use background files
      //event_time = m(-fGateWidth, 6*fEleSamplingPeriod);
      event_time = fTimeZero;//fTrnd.Uniform(-fGateWidth/2.-fEleSamplingPeriod, fGateWidth-fEleSamplingPeriod);
      //event_time = fTrnd.Uniform((-fGateWidth+2*fEleSamplingPeriod), 8*fEleSamplingPeriod);
    } else {
      // Signal events occur at t = 0, 
      event_time = fTimeZero+gemdet->fGEMhits[ih].t;
    }
    //  cout<<event_time<<"  "<<ih<<endl;
    // Adding drift time and trigger_jitter
    time_zero = event_time - fTriggerOffset[igem] + fRTime0*1e9 - trigger_jitter;
    
    //cout << time_zero << " " << fTimeZero << " " << gemdet->fGEMhits[ih].t 
    //<< " " << trigger_jitter << " " << fRTime0*1e9 << endl;
    
#if DBG_AVA > 0
    if(time_zero>200.0)
      cout << "time_zero " << time_zero 
	   << "; evt time " << event_time 
	   << "; hit time " << gemdet->fGEMhits[ih].t
	   << "; drift time " << fRTime0*1e9
	   << endl;
#endif
    if (fRNIon > 0) {
      //cout << "AvaModel..." << endl;
      AvaModel (igem, gemdet, R, vv1, vv2, time_zero);
      //cout << "Done!" << endl;
    }
    
  }//end loop on hits
  //fFilledStrips = false;
  
  //Cumulate(gemdet, R);
    
  return 0;
}

//-------------------------------------------------------
// Helper functions for integration in AvaModel
inline static
Double_t IntegralY( Double_t* a, Int_t ix, Int_t nx, Int_t ny )
{
  register double sum = 0.;
  register int kx = ix*ny;
  for( Int_t jy = ny; jy != 0; --jy )
    sum += a[kx++];

  return sum;
}


void SBSDigGEMSimDig::CheckOut(SBSDigGEMDet* gemdet, 
			       TRandom3* R, 
			       gmn_tree* T)
{
  double commonmode = 0;
  int apv_ctr;
  for(size_t i = 0; i<gemdet->GEMPlanes.size(); i++){
    cout << i << " " << gemdet->GEMPlanes[i].GetNStrips() << endl;
    for(int j = 0; j<gemdet->GEMPlanes[i].GetNStrips(); j++){
      if(fDoCommonMode)
	if(j%128==0 && apv_ctr<fCommonModeArray.size())
	  commonmode = fCommonModeArray[apv_ctr++];

      if(gemdet->GEMPlanes[i].GetADCSum(j)>0){
	for(int k = 0; k<6; k++){
	  gemdet->GEMPlanes[j].AddADC(j, k, R->Gaus(commonmode, fPulseNoiseSigma));
	  //handle saturation
	  if(gemdet->GEMPlanes[j].GetADC(j, k)>pow(2, fADCbits) )gemdet->GEMPlanes[j].SetADC(j, k, pow(2, fADCbits) );
	}
	if(fDoZeroSup){
	  if(gemdet->GEMPlanes[i].GetADCSum(j)-commonmode*6>fZeroSup){
	    FillBBGEMTree(gemdet->GEMPlanes[i], T, j);
	  }
	}else{
	  FillBBGEMTree(gemdet->GEMPlanes[i], T, j); 
	}
	
      }
    }
  }  
}

void SBSDigGEMSimDig::FillBBGEMTree(SBSDigGEMPlane pl, gmn_tree* T, int j)
{
  short strip;
  if(pl.Module()<3){
    strip = j+pl.GetNStrips()*pl.Module();
    if(pl.ROangle()==0){
      T->Earm_BBGEM_1x_dighit_nstrips++;
      T->Earm_BBGEM_1x_dighit_strip->push_back(strip);
      T->Earm_BBGEM_1x_dighit_adc_0->push_back(pl.GetADC(strip, 0));
      T->Earm_BBGEM_1x_dighit_adc_1->push_back(pl.GetADC(strip, 1));
      T->Earm_BBGEM_1x_dighit_adc_2->push_back(pl.GetADC(strip, 2));
      T->Earm_BBGEM_1x_dighit_adc_3->push_back(pl.GetADC(strip, 3));
      T->Earm_BBGEM_1x_dighit_adc_4->push_back(pl.GetADC(strip, 4));
      T->Earm_BBGEM_1x_dighit_adc_5->push_back(pl.GetADC(strip, 5));
    }else{
      T->Earm_BBGEM_1y_dighit_nstrips++;
      T->Earm_BBGEM_1y_dighit_strip->push_back(strip);
      T->Earm_BBGEM_1y_dighit_adc_0->push_back(pl.GetADC(strip, 0));
      T->Earm_BBGEM_1y_dighit_adc_1->push_back(pl.GetADC(strip, 1));
      T->Earm_BBGEM_1y_dighit_adc_2->push_back(pl.GetADC(strip, 2));
      T->Earm_BBGEM_1y_dighit_adc_3->push_back(pl.GetADC(strip, 3));
      T->Earm_BBGEM_1y_dighit_adc_4->push_back(pl.GetADC(strip, 4));
      T->Earm_BBGEM_1y_dighit_adc_5->push_back(pl.GetADC(strip, 5));
    }
  }else if(pl.Module()<6){
    strip = j+pl.GetNStrips()*(pl.Module()-3);
    if(pl.ROangle()==0){
      T->Earm_BBGEM_2x_dighit_nstrips++;
      T->Earm_BBGEM_2x_dighit_strip->push_back(strip);
      T->Earm_BBGEM_2x_dighit_adc_0->push_back(pl.GetADC(strip, 0));
      T->Earm_BBGEM_2x_dighit_adc_1->push_back(pl.GetADC(strip, 1));
      T->Earm_BBGEM_2x_dighit_adc_2->push_back(pl.GetADC(strip, 2));
      T->Earm_BBGEM_2x_dighit_adc_3->push_back(pl.GetADC(strip, 3));
      T->Earm_BBGEM_2x_dighit_adc_4->push_back(pl.GetADC(strip, 4));
      T->Earm_BBGEM_2x_dighit_adc_5->push_back(pl.GetADC(strip, 5));
    }else{
      T->Earm_BBGEM_2y_dighit_nstrips++;
      T->Earm_BBGEM_2y_dighit_strip->push_back(strip);
      T->Earm_BBGEM_2y_dighit_adc_0->push_back(pl.GetADC(strip, 0));
      T->Earm_BBGEM_2y_dighit_adc_1->push_back(pl.GetADC(strip, 1));
      T->Earm_BBGEM_2y_dighit_adc_2->push_back(pl.GetADC(strip, 2));
      T->Earm_BBGEM_2y_dighit_adc_3->push_back(pl.GetADC(strip, 3));
      T->Earm_BBGEM_2y_dighit_adc_4->push_back(pl.GetADC(strip, 4));
      T->Earm_BBGEM_2y_dighit_adc_5->push_back(pl.GetADC(strip, 5));
    }
  }else if(pl.Module()<9){
    strip = j+pl.GetNStrips()*(pl.Module()-6);
    if(pl.ROangle()==0){
      T->Earm_BBGEM_3x_dighit_nstrips++;
      T->Earm_BBGEM_3x_dighit_strip->push_back(strip);
      T->Earm_BBGEM_3x_dighit_adc_0->push_back(pl.GetADC(strip, 0));
      T->Earm_BBGEM_3x_dighit_adc_1->push_back(pl.GetADC(strip, 1));
      T->Earm_BBGEM_3x_dighit_adc_2->push_back(pl.GetADC(strip, 2));
      T->Earm_BBGEM_3x_dighit_adc_3->push_back(pl.GetADC(strip, 3));
      T->Earm_BBGEM_3x_dighit_adc_4->push_back(pl.GetADC(strip, 4));
      T->Earm_BBGEM_3x_dighit_adc_5->push_back(pl.GetADC(strip, 5));
    }else{
      T->Earm_BBGEM_3y_dighit_nstrips++;
      T->Earm_BBGEM_3y_dighit_strip->push_back(strip);
      T->Earm_BBGEM_3y_dighit_adc_0->push_back(pl.GetADC(strip, 0));
      T->Earm_BBGEM_3y_dighit_adc_1->push_back(pl.GetADC(strip, 1));
      T->Earm_BBGEM_3y_dighit_adc_2->push_back(pl.GetADC(strip, 2));
      T->Earm_BBGEM_3y_dighit_adc_3->push_back(pl.GetADC(strip, 3));
      T->Earm_BBGEM_3y_dighit_adc_4->push_back(pl.GetADC(strip, 4));
      T->Earm_BBGEM_3y_dighit_adc_5->push_back(pl.GetADC(strip, 5));
    }
  }else if(pl.Module()<12){
    strip = j+pl.GetNStrips()*(pl.Module()-9);
    if(pl.ROangle()==0){
      T->Earm_BBGEM_4x_dighit_nstrips++;
      T->Earm_BBGEM_4x_dighit_strip->push_back(strip);
      T->Earm_BBGEM_4x_dighit_adc_0->push_back(pl.GetADC(strip, 0));
      T->Earm_BBGEM_4x_dighit_adc_1->push_back(pl.GetADC(strip, 1));
      T->Earm_BBGEM_4x_dighit_adc_2->push_back(pl.GetADC(strip, 2));
      T->Earm_BBGEM_4x_dighit_adc_3->push_back(pl.GetADC(strip, 3));
      T->Earm_BBGEM_4x_dighit_adc_4->push_back(pl.GetADC(strip, 4));
      T->Earm_BBGEM_4x_dighit_adc_5->push_back(pl.GetADC(strip, 5));
    }else{
      T->Earm_BBGEM_4y_dighit_nstrips++;
      T->Earm_BBGEM_4y_dighit_strip->push_back(strip);
      T->Earm_BBGEM_4y_dighit_adc_0->push_back(pl.GetADC(strip, 0));
      T->Earm_BBGEM_4y_dighit_adc_1->push_back(pl.GetADC(strip, 1));
      T->Earm_BBGEM_4y_dighit_adc_2->push_back(pl.GetADC(strip, 2));
      T->Earm_BBGEM_4y_dighit_adc_3->push_back(pl.GetADC(strip, 3));
      T->Earm_BBGEM_4y_dighit_adc_4->push_back(pl.GetADC(strip, 4));
      T->Earm_BBGEM_4y_dighit_adc_5->push_back(pl.GetADC(strip, 5));
    }
  }else{
    strip = j+pl.GetNStrips()*(pl.Module()-12);
    if(pl.ROangle()==0){
      T->Earm_BBGEM_5x_dighit_nstrips++;
      T->Earm_BBGEM_5x_dighit_strip->push_back(strip);
      T->Earm_BBGEM_5x_dighit_adc_0->push_back(pl.GetADC(strip, 0));
      T->Earm_BBGEM_5x_dighit_adc_1->push_back(pl.GetADC(strip, 1));
      T->Earm_BBGEM_5x_dighit_adc_2->push_back(pl.GetADC(strip, 2));
      T->Earm_BBGEM_5x_dighit_adc_3->push_back(pl.GetADC(strip, 3));
      T->Earm_BBGEM_5x_dighit_adc_4->push_back(pl.GetADC(strip, 4));
      T->Earm_BBGEM_5x_dighit_adc_5->push_back(pl.GetADC(strip, 5));
    }else{
      T->Earm_BBGEM_5y_dighit_nstrips++;
      T->Earm_BBGEM_5y_dighit_strip->push_back(strip);
      T->Earm_BBGEM_5y_dighit_adc_0->push_back(pl.GetADC(strip, 0));
      T->Earm_BBGEM_5y_dighit_adc_1->push_back(pl.GetADC(strip, 1));
      T->Earm_BBGEM_5y_dighit_adc_2->push_back(pl.GetADC(strip, 2));
      T->Earm_BBGEM_5y_dighit_adc_3->push_back(pl.GetADC(strip, 3));
      T->Earm_BBGEM_5y_dighit_adc_4->push_back(pl.GetADC(strip, 4));
      T->Earm_BBGEM_5y_dighit_adc_5->push_back(pl.GetADC(strip, 5));
    }
  }
  
}

/*
void
TGEMSBSDigitizedPlane::Cumulate (const TGEMSBSGEMHit *vv, Short_t type,
			      Short_t clusterID )
{
  // cumulate hits (strips signals)
  if (vv) {
    //if(vv->GetSize()>20) {cout<<vv->GetSize();getchar();}
    for( Int_t j=0; j < vv->GetSize(); j++ ) {
      Double_t tempSumADC=0;
      Int_t idx = vv->GetIdx(j);
      assert( idx >= 0 && idx < fNStrips );
      fType[idx] |= type;
      fTime[idx] = (fTime[idx] < vv->GetTime()) ? fTime[idx] : vv->GetTime();
      fCharge[idx] += vv->GetCharge(j);
      bool was_below = !( fTotADC[idx] > fThreshold );
      for( UInt_t k=0; k<fNSamples; k++ ) {
	Int_t nnn = vv->GetADC(j,k);
	fStripClusterADC[k][idx].push_back(nnn);// new
	//cout << nnn << " ";
	assert( nnn >= 0 );
	if( nnn == 0 ) continue;
	Int_t iadc = idx*fNSamples+k;
	//cout << fStripADC[iadc] << " ";
	fStripADC[iadc] = fStripADC[iadc] + nnn;
  fStripSimADC[iadc] = 0;

	//cout << fStripADC[iadc] << " ";
	fTotADC[idx] += nnn;
	tempSumADC   += nnn;
      }//cout << endl;
      if( was_below && fTotADC[idx] > fThreshold ) {
	assert( fNOT < fNStrips );
	fOverThr[fNOT] = idx;
	++fNOT;
      }
      fStripWeightInCluster[idx].push_back(tempSumADC/vv->GetTotalADC());
      //cout<<((Double_t)fTotADC[idx])/vv->GetTotalADC()<<endl; getchar();
      fStripClusters[idx].push_back(clusterID);

    }
    
    //do cross talk if requested, a big signal along the strips 
    //will induce a smaller signal as the bigger one going to the APV, 
    //the smaller signal will appear on strips that is 
    //about 32 channels away from the big signal
    if (!TGEMSBSSimDigitization::fDoCrossTalk) return;
    Int_t isLeft = fRan.Uniform(1.) < 0.5 ? -1 : 1;
    Double_t factor = TGEMSBSSimDigitization::fCrossFactor +
      fRan.Gaus(0., TGEMSBSSimDigitization::fCrossSigma);
    if (factor <= 0.) return; //no induced signal
    
    for( Int_t j=0; j < vv->GetSize(); j++ ) {
      Int_t idx = vv->GetIdx(j);
      assert( idx >= 0 && idx < fNStrips );
      
      Int_t idxInduce = idx + isLeft*TGEMSBSSimDigitization::fNCStripApart;
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
  fStripSimADC[iadc] = 0;
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

*/

/*
//___________________________________________________________________________________
void
TGEMSBSSimDigitization::Print(Option_t*) const
{
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
  cout << "    Trigger offsets: "; //<< fTriggerOffset 
  for(int i = 0; i<fManager->GetNChamber(); i++)cout << fTriggerOffset[i] << " ";
  cout << endl;
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
}

void
TGEMSBSSimDigitization::PrintCharges() const
{
  cout << " Chb  Pln  Strip  Typ    ADC    Charge      Time\n";
  for (UInt_t ic = 0; ic < fNChambers; ++ic)
    {
      for (UInt_t ip = 0; ip < fNROPlanes[ic]; ++ip)
	for (UInt_t ist = 0; ist < (UInt_t) GetNStrips(ic, ip); ++ist)
	  {
	    if (fDP[ic][ip]->GetCharge (ist) > 0)
	      cout << setw(4) << ic
		   << " " << setw(4) << ip
		   << " " << setw(6) << ist
		   << " " << setw(4) << GetType (ic, ip, ist)
		   << " " << setw(6) << GetTotADC (ic, ip, ist)
		   << fixed << setprecision(1)
		   << " " << setw(9) << GetCharge (ic, ip, ist)
		   << fixed << setprecision(3)
		   << " " << setw(9) << GetTime (ic, ip, ist)
		   << endl;
	  }
    }
}

Double_t
TGEMSBSSimDigitization::CommonMode(UInt_t i_mpd)
{
  if(fCommonModeArray.size() && fDoCommonMode){
    i_mpd = (i_mpd<fCommonModeArray.size() ? i_mpd: 0);
    return fCommonModeArray[i_mpd];
  }else{
    return 0;
  }
}

Double_t
TGEMSBSSimDigitization::ZeroSupThreshold(UInt_t i_mpd)
{
  //i_mpd
  if(fDoZeroSup){
    if(fDoCommonMode){
      return fZeroSup+CommonMode(i_mpd)*fEleSamplingPoints;
    }else{
      return fZeroSup;
    }
  }else{
    return -1000;
  }
}



// Tree methods
void
TGEMSBSSimDigitization::InitTree (const TGEMSBSSpec& spect, const TString& ofile)
{
  // fOFile = new TFile( ofile, "RECREATE");

  // if (fOFile == 0 || fOFile->IsZombie() )
  //   {
  //     cerr << "Error: cannot open output file " << ofile << endl;
  //     delete fOFile; fOFile = 0;
  //     return;
  //   }

  // fOTree = new TTree( treeName, "Tree of digitized values");
  // fOTree -> SetMaxTreeSize(100000000000);
  // //fOTree -> SetMaxTreeSize(1000);


  // // create the tree variables

  // //fOTree->Branch( eventBranchName, "TGEMSBSSimEvent", &fEvent );
}

	    //setting strip sample adc and adding pedestal noise
	    for (UInt_t ss = 0; ss < strip.fNsamp; ++ss){
	      strip.fADC[ss] = GetADC(ich, ip, idx, ss);
	      // cout << strip.fADC[ss] << " ";
	       strip.fADC[ss] += fTrnd.Gaus(0, fPulseNoiseSigma);//allowing negative value, before implementing common mode;
	       if(fDoCommonMode)strip.fADC[ss] += CommonMode(mpd_id);
	      // cout << strip.fADC[ss] << " ";
	       //saturation = 4000;
	      if(strip.fADC[ss]>saturation)strip.fADC[ss]=saturation;
	      // TODO: Shouldn't common mode be added at some point? Specially
	      // before encoding it into unsigned integers which don't
	      // take kindly to negative values?
	      SetSimADC(ich,ip,idx,ss,strip.fADC[ss]);
	      const vector<Int_t>& sclust = GetStripClusterADC(ich, ip, idx, ss);
	      
	      strip.fClusterRatio[ss].Set( sclust.size(), &sclust[0] );
	      //     if(sclust.size()!=0){
	      //	 cout<<idx<<" "<<ss<<" "<<sclust.size()<<" == "<<strip.fClusterRatio[ss].GetSize()<<" == "<<GetStripClusters(ich, ip, idx).size()<<endl; getchar();
	      //    }
	    }//cout << endl;
	    if(GetTotADC(ich, ip, idx)==0)
	      {
		cout<<"Type: "<<GetType(ich, ip, idx)<<" Charge: "<<GetCharge(ich, ip, idx)<<" Time: "
		    <<GetTime(ich, ip, idx)<<" cluster: "<<GetStripClusters(ich, ip, idx).size();
		getchar();
		continue;
	      }
	    
	    strip.fSigType = GetType(ich, ip, idx);
	    strip.fCharge  = GetCharge(ich, ip, idx);
	    strip.fTime1   = GetTime(ich, ip, idx);
	
	    const vector<Short_t>& sc = GetStripClusters(ich, ip, idx);
	    strip.fClusters.Set( sc.size(), &sc[0] );
	    const vector<Double_t>& swc = GetStripWeightInCluster(ich, ip, idx);
	    strip.fStripWeightInCluster.Set( swc.size(), &swc[0] );

	    //fEvent->fGEMStrips.push_back( strip );
	  }
	}
      else{
	//modify this so recording all strips? 
	UInt_t nover = GetNOverThr(ich, ip);
	for (UInt_t iover = 0; iover < nover; iover++) {
	  Short_t idx = GetIdxOverThr(ich, ip, iover);
	  strip.fChan = idx;

	  //setting strip sample adc and adding pedestal noise
	  for (UInt_t ss = 0; ss < strip.fNsamp; ++ss){
	    strip.fADC[ss] = GetADC(ich, ip, idx, ss);
	    // cout << strip.fADC[ss] << " ";
	    strip.fADC[ss] += fTrnd.Gaus(0, fPulseNoiseSigma);//allowing negative value, before implementing common mode;
	    // cout << strip.fADC[ss] << " ";
	    if(strip.fADC[ss]>saturation)strip.fADC[ss]=saturation;

	    const vector<Int_t>& sclust = GetStripClusterADC(ich, ip, idx, ss);
	    strip.fClusterRatio[ss].Set( sclust.size(), &sclust[0] );
	    

	  }//cout << endl;

	  strip.fSigType = GetType(ich, ip, idx);
	  strip.fCharge  = GetCharge(ich, ip, idx);
	  strip.fTime1   = GetTime(ich, ip, idx);
	
	  const vector<Short_t>& sc = GetStripClusters(ich, ip, idx);
	  strip.fClusters.Set( sc.size(), &sc[0] );
	  
	  //fEvent->fGEMStrips.push_back( strip );
	}
      }
    }
  }
  fFilledStrips = true;
}

*/
