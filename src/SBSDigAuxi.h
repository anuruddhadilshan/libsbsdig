#ifndef SBSDIGAUXI_H
#define SBSDIGAUXI_H

#include <iostream>
#include <vector>
#include <map>
//#include "g4sbs_types.h"
#include "gmn_tree.h"
#include "TF1.h"
#include "TH1D.h"
#include "TRandom3.h"

bool UnfoldData(gmn_tree* T, double theta_sbs, double d_hcal, TRandom3* R);

/*
//
// classes for signal digitization
//
//________________________________
class SPEModel {
 public:
  SPEModel(const char* detname, double tau, double sigma, double t0 = 0, double tmin = -100, double tmax = +100);
  double Eval(double t){return fPulseHisto->Interpolate(t);};
  bool   PulseOverThr(double charge, double thr);
  bool   FindLeadTrailTime(double charge, double thr, double &t_lead, double &t_trail);
  
 private:
  TH1D *fPulseHisto;//At least we'll have to setup and use one per detector
  double GetHistoX(double y, double x1, double x2);
  //void  BuildHisto(double tau, double sigma);
  
  ClassDef(SPEModel,1);
};
*/
/*
//_________________________________
class PMTSignal {
 public:
  PMTSignal();
  PMTSignal(double npechargeconv);
  void Fill(SPEModel *model, int npe, double thr, double evttime, int signal);
  //void Digitize(TDigInfo diginfo, int chan);
  void Clear(Option_t* opt = "");
  ~PMTSignal(){Clear();};
  
  void AddSumEdep(double edep){
    fSumEdep+= edep;
  };
  void SetNpeChargeConv(double npechargeconv){fNpeChargeConv = npechargeconv;};
    
  double SumEdep(){return fSumEdep;};
  UInt_t Npe(){return fNpe;};
  double Charge(){return fNpe*fNpeChargeConv;};
  UInt_t ADC(){return fADC;};

  double EventTime(){return fEventTime;};
  UInt_t LeadTimesSize(){return fLeadTimes.size();};
  double LeadTime(int i){return fLeadTimes.at(i);};
  UInt_t TrailTimesSize(){return fTrailTimes.size();};
  double TrailTime(int i){return fTrailTimes.at(i);};
  UInt_t TDCSize(){return fTDCs.size();};
  UInt_t TDC(int i){return fTDCs.at(i);};
  //SimEncoder::tdc_data TDCData() { return fTDCData; }
  
  //check vectors size
  bool check_vec_size(bool ignore_edep = false);
  
 private:
  //summing variables for dig...
  double fSumEdep;//Not forced to use it for everything
  UInt_t fNpe;
  double fNpeChargeConv;
  UInt_t fADC;// One unique ADC value ?

  double fEventTime;
  //TDCs: multiple values possible.
  std::vector<double> fLeadTimes;
  std::vector<double> fTrailTimes;
  std::vector<UInt_t> fTDCs;
  //SimEncoder::tdc_data fTDCData;
  //TRndmManager* fRN;
  
  ClassDef(PMTSignal,1);
};
*/

#endif // SBSDIGAUXI_H

