#include "SBSDigAuxi.h"
#include "TMath.h"
#define DEBUG 0

using namespace std;

bool UnfoldData(gmn_tree* T, double theta_sbs, double d_hcal, TRandom3* R)
{
  //for the time being
  const double m_e = 511.e-6;
  const double n_lg = 1.68;
  
  bool has_data = false;
  
  int Npe;
  double t;
  
  double x_ref = d_hcal*sin(theta_sbs);
  double z_ref = d_hcal*cos(theta_sbs);
  
  double z_hit, Npe_Edep_ratio, sigma_tgen;

  double z_det, pz, E, beta;
  double sin2thetaC;

  // Process GRINCH data
  if(T->Earm_GRINCH_hit_nhits){
    for(int i = 0; i<T->Earm_GRINCH_hit_nhits; i++){
      //g4sbshitdata *grinchhit = new g4sbshitdata(GRINCH_UNIQUE_DETID, 4);
      //grinchhit->SetData(0, fSource);
      //grinchhit->SetData(1, 
      int(T->Earm_GRINCH_hit_PMT->at(i)/5)-1;//);
      // -1 is to range from 0 to 509 instead of 1 to 510 :/
      //grinchhit->SetData(2, 
      T->Earm_GRINCH_hit_Time_avg->at(i);//);
      //grinchhit->SetData(3, 
      T->Earm_GRINCH_hit_NumPhotoelectrons->at(i);//);
      //fg4sbsHitData.push_back(grinchhit);
    }
    has_data = true;
  }

  // Process hodoscope data
  if(T->Earm_BBHodoScint_hit_nhits){
    for(int i = 0; i<T->Earm_BBHodoScint_hit_nhits; i++){
      for(int j = 0; j<2; j++){//j = 0: close beam PMT, j = 1: far beam PMT
	// Evaluation of number of photoelectrons and time from energy deposit documented at:
	// https://hallaweb.jlab.org/dvcslog/SBS/170711_172759/BB_hodoscope_restudy_update_20170711.pdf
	// TODO: put that stuff in DB...
	Npe = R->Poisson(1.0e7*T->Earm_BBHodoScint_hit_sumedep->at(i)*0.113187*exp(-(0.3+pow(-1, j)*T->Earm_BBHodoScint_hit_xhit->at(i))/1.03533)* 0.24);
	t = T->Earm_BBHodoScint_hit_tavg->at(i)+(0.55+pow(-1, j)*T->Earm_BBHodoScint_hit_xhit->at(i))/0.15;
	//g4sbshitdata *hodopmthit = new g4sbshitdata(HODO_UNIQUE_DETID, 5);
	//hodopmthit->SetData(0, fSource);
	//hodopmthit->SetData(1, 
	T->Earm_BBHodoScint_hit_cell->at(i)*2+j;//);
	//hodopmthit->SetData(2, t);
	//hodopmthit->SetData(3, Npe);
	//hodopmthit->SetData(4, 
	T->Earm_BBHodoScint_hit_sumedep->at(i);//);
	//fg4sbsHitData.push_back(hodopmthit);
      }
    }
    has_data = true;
  }
  
  // Process BB PS data
  if(T->Earm_BBPSTF1_hit_nhits){
    for(int i = 0; i<T->Earm_BBPSTF1_hit_nhits; i++){
      // Evaluation of number of photoelectrons and time from energy deposit documented at:
      if(T->Earm_BBPSTF1_hit_sumedep->at(i)<1.e-4)continue;
      //check probability to generate p.e. yield
      bool genpeyield = true;
      if(T->Earm_BBPSTF1_hit_sumedep->at(i)<1.e-2)
	genpeyield = R->Uniform(0, 1)<=1-exp(0.29-950.*T->Earm_BBPSTF1_hit_sumedep->at(i));
      //if we're go, let's generate
      if(genpeyield){
	beta = sqrt( pow(m_e+T->Earm_BBPSTF1_hit_sumedep->at(i), 2)-m_e*m_e )/(m_e + T->Earm_BBPSTF1_hit_sumedep->at(i));
	sin2thetaC = TMath::Max(1.-1./pow(n_lg*beta, 2), 0.);
	//1500. Used to be 454.: just wrong
	Npe = R->Poisson(300.0*
			   T->Earm_BBPSTF1_hit_sumedep->at(i)*
			   sin2thetaC/(1.-1./(n_lg*n_lg)) 
			   );
	t = T->Earm_BBPSTF1_hit_tavg->at(i)+R->Gaus(3.2-5.805*T->Earm_BBPSTF1_hit_zhit->at(i)-17.77*pow(T->Earm_BBPSTF1_hit_zhit->at(i), 2), 0.5);
	//g4sbshitdata *bbpshit = new g4sbshitdata(BBPS_UNIQUE_DETID, 5);
	//bbpshit->SetData(0, fSource);
	//bbpshit->SetData(1, 
	T->Earm_BBPSTF1_hit_cell->at(i);//);
	//bbpshit->SetData(2, t);
	//bbpshit->SetData(3, Npe);
	//bbpshit->SetData(4, 
	T->Earm_BBPSTF1_hit_sumedep->at(i);//);
	//fg4sbsHitData.push_back(bbpshit);
      }
    }
    has_data = true;
  }
  
  if(T->Earm_BBSHTF1_hit_nhits){
    for(int i = 0; i<T->Earm_BBSHTF1_hit_nhits; i++){
      // Evaluation of number of photoelectrons and time from energy deposit documented at:
      // 
       if(T->Earm_BBSHTF1_hit_sumedep->at(i)<1.e-4)continue;
      //check probability to generate p.e. yield
      bool genpeyield = true;
      if(T->Earm_BBSHTF1_hit_sumedep->at(i)<1.e-2)genpeyield = R->Uniform(0, 1)<=1-exp(0.29-950.*T->Earm_BBSHTF1_hit_sumedep->at(i));
      //if we're go, let's generate
      if(genpeyield){
	beta = sqrt( pow(m_e+T->Earm_BBSHTF1_hit_sumedep->at(i), 2)-m_e*m_e )/(m_e + T->Earm_BBSHTF1_hit_sumedep->at(i));
	sin2thetaC = TMath::Max(1.-1./pow(n_lg*beta, 2), 0.);
	//1800. Used to be 932.: just wrong
	Npe = R->Poisson(360.0*
			   T->Earm_BBSHTF1_hit_sumedep->at(i)*
			   sin2thetaC/(1.-1./(n_lg*n_lg)) 
			   );
	t = T->Earm_BBSHTF1_hit_tavg->at(i)+R->Gaus(2.216-8.601*T->Earm_BBSHTF1_hit_zhit->at(i)-7.469*pow(T->Earm_BBSHTF1_hit_zhit->at(i), 2), 0.8);
	//g4sbshitdata *bbshhit = new g4sbshitdata(BBSH_UNIQUE_DETID, 5);
	//bbshhit->SetData(0, fSource);
	//bbshhit->SetData(1, 
	T->Earm_BBSHTF1_hit_cell->at(i);//);
	//bbshhit->SetData(2, 0);
	//bbshhit->SetData(2, t);
	//bbshhit->SetData(3, Npe);
	//bbshhit->SetData(4, 
	T->Earm_BBSHTF1_hit_sumedep->at(i);//);
	//fg4sbsHitData.push_back(bbshhit);
      }
    }
  }
  
  if(T->Harm_HCalScint_hit_sumedep) {
    for(size_t k = 0; k < T->Harm_HCalScint_hit_sumedep->size(); k++) {
      //g4sbshitdata *hcalscinthit = new g4sbshitdata(HCAL_UNIQUE_DETID,4);
      //hcalscinthit->SetData(0,fSource);
      //hcalscinthit->SetData(1,
      T->Harm_HCalScint_hit_cell->at(k);//);
      //hcalscinthit->SetData(2,1); // data type 1 == sumedep
      //hcalscinthit->SetData(3,T->Harm_HCalScint_hit_sumedep->at(k));
      //fg4sbsHitData.push_back(hcalscinthit);

      z_hit = -(T->Harm_HCalScint_hit_xhitg->at(k)-x_ref)*sin(theta_sbs)+(T->Harm_HCalScint_hit_zhitg->at(k)-z_ref)*cos(theta_sbs);
      
      // Evaluation of number of photoelectrons from energy deposit documented at:
      // https://sbs.jlab.org/DocDB/0000/000043/002/Harm_HCal_Digi_EdepOnly_2.pdf
      // TODO: put that stuff in DB...
      Npe_Edep_ratio = 5.242+11.39*z_hit+10.41*pow(z_hit, 2);
      Npe = R->Poisson(Npe_Edep_ratio*T->Harm_HCalScint_hit_sumedep->at(k)*1.0e3);
      t = R->Gaus(T->Harm_HCalScint_hit_tavg->at(k)+10.11, 1.912);
      
      sigma_tgen = 0.4244+11380/pow(Npe+153.4, 2);
      //for(int l = 0; l<Npe; l++){
      //Generate here,...
      //g4sbshitdata *hcalpmthit = new g4sbshitdata(HCAL_UNIQUE_DETID,4);
      //hcalpmthit->SetData(0,fSource);
      //hcalpmthit->SetData(1,
      T->Harm_HCalScint_hit_cell->at(k);//);
      //hcalpmthit->SetData(2,0); // data type 0 == optical photon
      //hcalpmthit->SetData(3,
      R->Landau(t, sigma_tgen);//);
      //fg4sbsHitData.push_back(hcalpmthit);
      //}
    }
  }
  
  // Now process the GEM data
  /*
  for(std::vector<gem_branch>::iterator it = fTree->GEMs.begin();
       it != fTree->GEMs.end(); it++) {
    if(d_flag>=3)printf("GEM tree %ld\n", std::distance(fTree->GEMs.begin(), it));
    
    TSBSGeant4::GEMData_t &t = it->tree;
    if(t.plane) {
      if(d_flag>=3)printf("%d hits \n", t.nhits);
      for(int k = 0; k < t.nhits; k++) {
        // Don't bother with events that deposited no energy
        if(t.edep->at(k)>0) {
          g4sbshitdata *gemhit = new g4sbshitdata(it->id, 35);
          // For now, just copy the whole tree, we'll trim it later I guess
          // (with slight modifications as prescribed in libsbsgem)
          gemhit->SetData(0,fSource);
          gemhit->SetData(1,t.plane->at(k));
          gemhit->SetData(2,t.strip->at(k));
          gemhit->SetData(3,t.x->at(k));
          gemhit->SetData(4,t.y->at(k));
          gemhit->SetData(5,t.z->at(k));
          gemhit->SetData(6,t.polx->at(k));
          gemhit->SetData(7,t.poly->at(k));
          gemhit->SetData(8,t.polz->at(k));
          gemhit->SetData(9,t.t->at(k));
          gemhit->SetData(10,t.trms->at(k));
          gemhit->SetData(11,t.tmin->at(k));
          gemhit->SetData(12,t.tmax->at(k));
          gemhit->SetData(13,t.tx->at(k));
          gemhit->SetData(14,t.ty->at(k));
          gemhit->SetData(15,t.txp->at(k));
          gemhit->SetData(16,t.typ->at(k));
          gemhit->SetData(17,t.xg->at(k));
          gemhit->SetData(18,t.yg->at(k));
          gemhit->SetData(19,t.zg->at(k));
          gemhit->SetData(20,t.trid->at(k));
          gemhit->SetData(21,t.mid->at(k)+1);
          gemhit->SetData(22,t.pid->at(k));
          gemhit->SetData(23,t.vx->at(k));
          gemhit->SetData(24,t.vy->at(k));
          gemhit->SetData(25,t.vz->at(k));
          gemhit->SetData(26,t.p->at(k));
          gemhit->SetData(27,t.edep->at(k)); // convert to MeV?
          gemhit->SetData(28,t.beta->at(k));
          gemhit->SetData(29,t.xin->at(k));
          gemhit->SetData(30,t.yin->at(k));
          gemhit->SetData(31,t.zin->at(k));
          gemhit->SetData(32,t.xout->at(k));
          gemhit->SetData(33,t.yout->at(k));
          gemhit->SetData(34,t.zout->at(k));
          fg4sbsHitData.push_back(gemhit);
        }
      }
    }
  }  
  */
  return has_data;
}


//
// Class SPEModel
//
/*
SPEModel::SPEModel(const char* detname, 
		   double tau, double sigma, 
		   double t0, double tmin, double tmax)
{
  TF1 fFunc1(Form("fFunc1%s",detname),
  	     TString::Format("TMath::Max(0.,(x/%g)*TMath::Exp(-x/(%g)))", 
  			     tau*tau, tau),
  	     tmin ,tmax);
  TF1 fFunc2(Form("fFunc2%s",detname), "gaus", 
  // 	     TString::Format("%g*TMath::Exp(-((x-%g)**2)/(%g))", 
  // 			     1./TMath::Sqrt(2*TMath::Pi()*sigma), t0, sigma*sigma),
   	     tmin, tmax);
  fFunc2.SetParameters(1.0e0/10.0/(sqrt(2*TMath::Pi())*sigma), 0, sigma);//do it this way for thaat purpose
  
  const int NbinsTotal = int(tmax-tmin)*10;// 10 bins/ns should do... since we will extrapolate after...
  fPulseHisto = new TH1D(Form("fPulseHisto%s",detname), "", NbinsTotal, tmin, tmax);
  double t_i, t_j;
  double ps_i, g_j;
  for(int i = 1; i<=NbinsTotal; i++){
    t_i = fPulseHisto->GetBinCenter(i+1);
    ps_i = fFunc1.Eval(t_i);
    if(sigma>0){
      fFunc2.SetParameter(1, t_i);
      for(int j = 1; j<=NbinsTotal; j++){
	t_j = fPulseHisto->GetBinCenter(j+1);
	g_j = fFunc2.Eval(t_j);
	fPulseHisto->Fill(t_j, ps_i*g_j);
      }
    }else{
      fPulseHisto->Fill(t_i, ps_i);
    }
  }
  
  cout << endl<< detname << " pulse histo built with address << " << fPulseHisto << endl;
}

bool SPEModel::PulseOverThr(double charge, double thr)
{
  if(fPulseHisto->GetMaximum()<thr/charge){
    return false;
  }else{
    return true;
  }
};
 
bool SPEModel::FindLeadTrailTime(double charge, double thr, double &t_lead, double &t_trail)
{
  if(!PulseOverThr(charge, thr)){
    t_lead = 1.0e38;
    t_trail = 1.0e38;
    return false;
  }else{
    double xmax = fPulseHisto->GetBinCenter(fPulseHisto->GetMaximumBin());
    if(fPulseHisto->GetBinContent(1)<thr/charge){
      t_lead = GetHistoX(thr/charge, fPulseHisto->GetBinLowEdge(1), xmax);
    }else{
      t_lead = 1.0e38;
      return false;
    }
    if(fPulseHisto->GetBinContent(fPulseHisto->GetNbinsX())<thr/charge){
      t_trail = GetHistoX(thr/charge, xmax, fPulseHisto->GetBinLowEdge(fPulseHisto->GetNbinsX()+1));
    }else{
      t_trail = 1.0e38;
      return false;
    }
    //why does this function seem to oblitarate fMCHit containers when t_trail is calculated to be 1.e38... GetHistoX?
    if(t_trail>1.e30) cout << "thr/charge" << thr/charge << " >? " <<  fPulseHisto->GetBinContent(fPulseHisto->GetNbinsX()) << " or " << fPulseHisto->GetBinContent(fPulseHisto->GetNbinsX()+1) << endl;
    return true;
  }
}

double SPEModel::GetHistoX(double y, double x1, double x2)
{
  double splineslope;
  for(int k = fPulseHisto->FindBin(x1); k<=fPulseHisto->FindBin(x2); k++){
    if(  ( (fPulseHisto->GetBinContent(k+1)-y)*(fPulseHisto->GetBinContent(k)-y) )<0 ){
      // threshold crossed if diff(TH1::GetBinContent-y) changes sign
      splineslope = (fPulseHisto->GetBinContent(k+1)-fPulseHisto->GetBinContent(k))/(fPulseHisto->GetBinCenter(k+1)-fPulseHisto->GetBinCenter(k));
      
      return fPulseHisto->GetBinCenter(k)+(y-fPulseHisto->GetBinContent(k))/splineslope;
    }
  }
  return 1.0e38;
}


//
// Class PMTSignal
//
PMTSignal::PMTSignal()
  : fSumEdep(0), fNpe(0), fNpeChargeConv(1.0), fADC(0), fEventTime(0), fMCHitSize(0)
{ 
  fLeadTimes.clear();
  fTrailTimes.clear();
  fTDCs.clear();
  
  //R = TRndmManager::GetInstance();
}

PMTSignal::PMTSignal(double npechargeconv)
  : fSumEdep(0), fNpe(0), fNpeChargeConv(npechargeconv), fADC(0), fEventTime(0), fMCHitSize(0)
{ 
  fLeadTimes.clear();
  fTrailTimes.clear();
  fTDCs.clear();
}

void PMTSignal::Fill(SPEModel *model, int npe, double thr, double evttime, int signal)
{
  if(signal==0)fEventTime = evttime;
  fNpe+= npe;
  //if(model->PulseOverThr(fCharge, thr))fNpe//fADC = model->GetCharge()*model->GetADCconversion();
  
  //determine lead and trail times
  double t_lead, t_trail;
  // find the lead and trail time for *this* pulse, not the total pulse
  bool goodtime = model->FindLeadTrailTime(npe*fNpeChargeConv, thr, t_lead, t_trail);
  
  t_lead+=evttime;
  t_trail+=evttime;
  
  if(goodtime){
    //Filter here the lead and trail times
    if(fLeadTimes.size()>0){
      // Check if the lead and trail times are inside an existing lead time- trail time pair
      // we assume here that fLeadTimes and fTrailTimes are same size 
      // *if we do things correctly, that should be the case*
      // we shall keep lead-trail times pair in timing order
      // we neglect pulse overlaps ftm.
      if(fLeadTimes.size()!=fTrailTimes.size()){
	cout << " B - Warning: size of lead times container: " << fLeadTimes.size() 
	     << " != size of trail times container: " << fTrailTimes.size() << endl;
      }
      
      for(size_t i = 0; i<fLeadTimes.size(); i++){
	// possibility of the current pair straddling with others.... :/
	// treat those separately to simplify...
	// tL < tT_i-1 < tL_i < tT
	if(i>0){
	  if(t_lead < fTrailTimes.at(i-1) && fLeadTimes.at(i) < t_trail){
	    //do necessary substitutions first, then erase
	    // tL < tL_i-1 => tL *replaces* tL_i-1
	    if(t_lead < fLeadTimes.at(i-1)){
	      fLeadTimes.erase(fLeadTimes.begin()+i-1);
	      fLeadTimes.insert(fLeadTimes.begin()+i-1, t_lead);
	    }
	    // tT_i < tT => tT *replaces* tT_i
	    if(fTrailTimes.at(i) < t_trail){
	      fTrailTimes.erase(fTrailTimes.begin()+i);
	      fTrailTimes.insert(fTrailTimes.begin()+i, t_trail);
	    }
	    fLeadTimes.erase(fLeadTimes.begin()+i);
	    fTrailTimes.erase(fTrailTimes.begin()+i-1);
	    break;
	  }
	}//end if(i>0)
	// tL < tT_i < tL_i+1 < tT
	if(i<fLeadTimes.size()-1){
	  if(t_lead < fTrailTimes.at(i) && fLeadTimes.at(i+1) < t_trail){
	    //do necessary substitutions first, then erase
	    // tL < tL_i => tL *replaces* tL_i
	    if(t_lead < fLeadTimes.at(i)){
	      fLeadTimes.erase(fLeadTimes.begin()+i);
	      fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	    }
	    // tT_i+1 < tT => tT *replaces* tT_i+1
	    if(fTrailTimes.at(i+1) < t_trail){
	      fTrailTimes.erase(fTrailTimes.begin()+i+1);
	      fTrailTimes.insert(fTrailTimes.begin()+i+1, t_trail);
	    }
	    fLeadTimes.erase(fLeadTimes.begin()+i+1);
	    fTrailTimes.erase(fTrailTimes.begin()+i);
	    break;
	  }
	}//end if(i<fLeadTimes.size()-1)
	
	// if not, 6 cases to consider:
	// tL < tT < tL_i < tT_i => both inserted *before* existing pair 
	if(t_trail < fLeadTimes.at(i)){
	  fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	  fTrailTimes.insert(fTrailTimes.begin()+i, t_trail); 
	  break;
	}
	// tL < tL_i < tT < tT_i => tL *replaces* tL_i
	if(t_lead < fLeadTimes.at(i) && fLeadTimes.at(i) < t_trail && t_trail < fTrailTimes.at(i)){
	  fLeadTimes.erase(fLeadTimes.begin()+i);
	  fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	  break;
	}
	// tL_i < tL < tT < tT_i => tL *replaces* tL_i AND tT *replaces* tT_i
	if(t_lead < fLeadTimes.at(i) && fTrailTimes.at(i) < t_trail){
	  fLeadTimes.erase(fLeadTimes.begin()+i);
	  fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	  fTrailTimes.erase(fTrailTimes.begin()+i);
	  fTrailTimes.insert(fTrailTimes.begin()+i, t_trail);
	  break;
	}
	if(fLeadTimes.at(i) < t_lead && t_trail < fTrailTimes.at(i)){
	  break;
	}
	if(fLeadTimes.at(i) < t_lead && t_lead < fTrailTimes.at(i) && fTrailTimes.at(i) < t_trail){
	  fTrailTimes.erase(fTrailTimes.begin()+i);
	  fTrailTimes.insert(fTrailTimes.begin()+i, t_trail);
	}
	if(fTrailTimes.at(i) < t_lead){
	  fLeadTimes.insert(fLeadTimes.begin()+i+1, t_lead);
	  fTrailTimes.insert(fTrailTimes.begin()+i+1, t_trail);
	  break;
	}
	if(fLeadTimes.size()!=fTrailTimes.size()){
	  cout << " A - Warning: size of lead times container: " << fLeadTimes.size() 
	       << " != size of trail times container: " << fTrailTimes.size() << endl;
	}
	if(i>=fLeadTimes.size())cout << "Warning: i = " << i << " >= size of containers = " << fLeadTimes.size() << endl;
      }
    }else{
      //of course, if initial size was 0, just psuh it back
      //hopefully it will be the case most of the time
      fLeadTimes.push_back(t_lead);
      fTrailTimes.push_back(t_trail);
    }
    if(fLeadTimes.size()!=fTrailTimes.size()){
      cout << " A - Warning: size of lead times container: " << fLeadTimes.size() 
	   << " != size of trail times container: " << fTrailTimes.size() << endl;
    }
  }//end if(t_lead && t_trail<30)
}
*/
/*
void PMTSignal::Digitize(TDigInfo diginfo, int chan)
{
  if(fNpe<=0){
    fADC = R->Gaus(diginfo.Pedestal(chan), diginfo.PedestalNoise(chan));
    return;
  }
    
#if DEBUG>0
  cout << "Charge (C) " << Charge() << " (fC) " << Charge()*1.0e15 << ", ADC conversion (fC/ch) " << diginfo.ADCConversion();
#endif
  
  fADC = TMath::Nint(Charge()*1.0e15/diginfo.ADCConversion()+R->Gaus(diginfo.Pedestal(chan), diginfo.PedestalNoise(chan)));
  //if ADC value bigger than number of ADC bits, ADC saturates
  if( fADC>UInt_t(TMath::Nint( TMath::Power(2, diginfo.ADCBits()) )) ){
    fADC = TMath::Nint( TMath::Power(2, diginfo.ADCBits()) );
  }
  //cout << "PMTSignal::Digitize():  " << fLeadTimes.size() << " - " << fTrailTimes.size() << endl;
  
#if DEBUG>0
  cout << " => ADC = " << fADC << endl;
#endif
  
  UInt_t tdc_value;
  
  // For the sake of going forward, we assume that the signal is the first entry of each vector
  fTDCData.time.clear();
  if(fLeadTimes.size() && fTrailTimes.size()){
    for(size_t i = 0; i<fLeadTimes.size(); i++){
#if DEBUG>0
      cout << " fLeadTimes.at(" << i << ") " << fLeadTimes.at(i) 
	   << " fTrailTimes.at(" << i << ") " << fTrailTimes.at(i) << endl;
#endif
      // trim "all" bits that are above the number of TDC bits - a couple to speed it up
      // (since TDC have a revolving clock, as far as I understand)
      // let's use an arbitrary reference time offset of 1us before the trigger
      tdc_value = TMath::Nint((fLeadTimes.at(i)+1.e3)/diginfo.TDCConversion());
      for(int j = 30; j>=diginfo.TDCBits(); j--){
	tdc_value ^= ( -0 ^ tdc_value) & ( 1 << (j) );
      }
      tdc_value ^= ( -0 ^ tdc_value) & ( 1 << (31) );
      fTDCData.time.push_back(tdc_value);
      //fTDCs.insert(fTDCs.begin()+0, TMath::Nint(fLeadTimes.at(0)*diginfo.TDCConversion()));//bug!!!!
      fTDCs.push_back(tdc_value);//they're already sorted in order, presumably
      // also mark the traling time with setting bin 31 to 1 // need to reconvert then
      tdc_value = TMath::Nint((fTrailTimes.at(i)+1.e3)/diginfo.TDCConversion());
      for(int j = 30; j>=diginfo.TDCBits(); j--){
	tdc_value ^= ( -0 ^ tdc_value) & ( 1 << (j) );
      }
      tdc_value ^= ( -1 ^ tdc_value) & ( 1 << (31) );
      fTDCs.push_back(tdc_value);
      fTDCData.time.push_back(tdc_value);
      
#if DEBUG>0
      cout << " fTDCs.at(0) " << fTDCs.at(0) << " fTDCs.at(1) " << fTDCs.at(1) << endl;
#endif
    }
  }
  
  fSumEdep*=1.0e9;// store in eV.
}

void PMTSignal::Clear(Option_t*)
{
  //cout << " PMTSignal::Clear() " << endl;
  
  fSumEdep = 0;
  fNpe = 0;
  fADC = 0;
  
  fEventTime = 0;
  fLeadTimes.clear();
  fTrailTimes.clear();
  fTDCs.clear();
}

ClassImp(PMTSignal)
ClassImp(SPEModel)
*/

