#ifndef _G4SBS_DATA_H
#define _G4SBS_DATA_H

#include <Rtypes.h>
class TTree;

namespace TSBSGeant4 {

  // Purely virtual data structure for detector data in a G4SBS ROOT tree.
  struct VDetData_t {
    VDetData_t() {};
    virtual ~VDetData_t(){};
    // All sub-classes *must* implement a concrete instance of this
    virtual bool SetupBranches(TTree *t, const char* prefix) = 0;
    template<typename T>
    int SetupBranch(TTree* tree,const char* prefix,const char* varname,T &var);
  };

  struct CalData_t : public VDetData_t {
    int nhits;
    std::vector<int>    *row;
    std::vector<int>    *col;
    std::vector<int>    *cell;
    std::vector<int>    *plane;
    //std::vector<int>    *wire;
    std::vector<double> *xcell;
    std::vector<double> *ycell;
    std::vector<double> *zcell;
    std::vector<double> *xcellg;
    std::vector<double> *ycellg;
    std::vector<double> *zcellg;
    std::vector<double> *xhit;
    std::vector<double> *yhit;
    std::vector<double> *zhit;
    std::vector<double> *sumedep;
    std::vector<double> *tavg;
    std::vector<double> *trms;
    std::vector<double> *tmin;
    std::vector<double> *tmax;
    CalData_t() :
      nhits(0), row(0), col(0), cell(0), plane(0), xcell(0), ycell(0), zcell(0),
      xcellg(0), ycellg(0), zcellg(0), xhit(0), yhit(0), zhit(0), sumedep(0),
      tavg(0), trms(0), tmin(0), tmax(0) {}
    virtual ~CalData_t(){};
    virtual bool SetupBranches(TTree *t, const char* prefix);
  };

  // Standard ECal data (on PMT)
  struct ECalData_t : public VDetData_t {
    int nhits;
    std::vector<int>    *PMT;
    std::vector<int>    *row;
    std::vector<int>    *col;
    std::vector<int>    *plane;
    std::vector<double> *xcell;
    std::vector<double> *ycell;
    std::vector<double> *zcell;
    std::vector<double> *xgcell;
    std::vector<double> *ygcell;
    std::vector<double> *zgcell;
    std::vector<int>    *NumPhotoelectrons;
    std::vector<double> *Time_avg;
    std::vector<double> *Time_rms;
    std::vector<double> *Time_min;
    std::vector<double> *Time_max;
    ECalData_t() :
      nhits(0), PMT(0), row(0), col(0), plane(0), xcell(0),
      ycell(0), zcell(0), xgcell(0), ygcell(0), zgcell(0), NumPhotoelectrons(0),
      Time_avg(0), Time_rms(0), Time_min(0), Time_max(0) {}
    virtual ~ECalData_t(){};
    virtual bool SetupBranches(TTree *t, const char* prefix);
  };

  // If individual optical photons are included on the PMT
  //struct ECalPartData_t : public ECalData_t {
  struct ECalPartData_t : public VDetData_t {
    int npart_ECAL;
    std::vector<double> *E;
    std::vector<double> *t;
    std::vector<int>    *part_PMT;
    std::vector<bool>   *detected;
    // Quick default constructor sets all pointers to zero for safety
    ECalPartData_t() : npart_ECAL(0), E(0), t(0), part_PMT(0), detected(0) {};
    virtual ~ECalPartData_t(){};
    virtual bool SetupBranches(TTree *t, const char *prefix);
  };

  // RICHoutput type branches
  struct RICHData_t : public VDetData_t {
    Int_t                  nhits;
    std::vector<int>      *PMT;
    std::vector<int>      *row;
    std::vector<int>      *col;
    std::vector<double>   *xpmt;
    std::vector<double>   *ypmt;
    std::vector<double>   *zpmt;
    std::vector<double>   *xgpmt;
    std::vector<double>   *ygpmt;
    std::vector<double>   *zgpmt;
    std::vector<int>      *NumPhotoelectrons;
    std::vector<double>   *Time_avg;
    std::vector<double>   *Time_rms;
    std::vector<double>   *Time_min;
    std::vector<double>   *Time_max;
    std::vector<int>      *mTrackNo;
    std::vector<double>   *xhit;
    std::vector<double>   *yhit;
    std::vector<double>   *zhit;
    std::vector<double>   *pxhit;
    std::vector<double>   *pyhit;
    std::vector<double>   *pzhit;
    std::vector<double>   *pvx;
    std::vector<double>   *pvy;
    std::vector<double>   *pvz;
    std::vector<double>   *ppx;
    std::vector<double>   *ppy;
    std::vector<double>   *ppz;
    std::vector<int>      *volume_flag;
    RICHData_t() :
      nhits(0), PMT(0), row(0), col(0), xpmt(0), ypmt(0), zpmt(0), xgpmt(0),
      ygpmt(0), zgpmt(0), NumPhotoelectrons(0), Time_avg(0), Time_rms(0),
      Time_min(0), Time_max(0), mTrackNo(0), xhit(0), yhit(0), zhit(0),
      pxhit(0), pyhit(0), pzhit(0), pvx(0), pvy(0), pvz(0), ppx(0), ppy(0),
      ppz(0), volume_flag(0) {};
    virtual ~RICHData_t() {};
    virtual bool SetupBranches(TTree *t, const char *prefix);
  };

  struct GEMData_t : public VDetData_t {
    Int_t                 nhits;
    std::vector<int>      *plane;
    std::vector<int>      *strip;
    std::vector<double>   *x;
    std::vector<double>   *y;
    std::vector<double>   *z;
    std::vector<double>   *polx;
    std::vector<double>   *poly;
    std::vector<double>   *polz;
    std::vector<double>   *t;
    std::vector<double>   *trms;
    std::vector<double>   *tmin;
    std::vector<double>   *tmax;
    std::vector<double>   *tx;
    std::vector<double>   *ty;
    std::vector<double>   *txp;
    std::vector<double>   *typ;
    std::vector<double>   *xg;
    std::vector<double>   *yg;
    std::vector<double>   *zg;
    std::vector<int>      *trid;
    std::vector<int>      *mid;
    std::vector<int>      *pid;
    std::vector<double>   *vx;
    std::vector<double>   *vy;
    std::vector<double>   *vz;
    std::vector<double>   *p;
    std::vector<double>   *edep;
    std::vector<double>   *beta;
    GEMData_t(){};
    virtual ~GEMData_t() {};
    virtual bool SetupBranches(TTree *t, const char *prefix);
  };

  struct TrackerData_t : public VDetData_t {
    Int_t                 ntracks;
    std::vector<int>      *TID;
    std::vector<int>      *PID;
    std::vector<int>      *MID;
    std::vector<int>      *NumHits;
    std::vector<int>      *NumPlanes;
    std::vector<int>      *NDF;
    std::vector<double>   *Chi2fit;
    std::vector<double>   *Chi2true;
    std::vector<double>   *X;
    std::vector<double>   *Y;
    std::vector<double>   *Xp;
    std::vector<double>   *Yp;
    std::vector<double>   *T;
    std::vector<double>   *P;
    std::vector<double>   *Sx;
    std::vector<double>   *Sy;
    std::vector<double>   *Sz;
    std::vector<double>   *Xfit;
    std::vector<double>   *Yfit;
    std::vector<double>   *Xpfit;
    std::vector<double>   *Ypfit;
    TrackerData_t(){};
    virtual ~TrackerData_t() {};
    virtual bool SetupBranches(TTree *t, const char *prefix);
  };

  }

#endif // _G4SBS_DATA_H
