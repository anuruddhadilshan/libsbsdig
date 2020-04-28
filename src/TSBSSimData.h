#ifndef TSBSSIMDATA_H
#define TSBSSIMDATA_H

#include <TVector3.h>
#include "g4sbs_types.h"

#define __DEFAULT_DATA_SIZE 32

////////////////////////////////////////////////////////////////////////////
// Auxilliary class for storing hit data
//
// Stores an arbitrary double data in dynamically allocated
// arrays.  Allows us to add in data as we get it and then check
// to make sure all entries in the array are filled
// // ___________________________________________________________ //
// // hit_data: {PMT row, PMT column, 
// //            N_pe, t,
// //            pdetx, pdety, pdetz, 
// //            Xdetx, Xdety, Xdetz, 
// //            type, 
// //            pprodx, pprody, pprodz, 
// //            Xprodx, Xprody, Xprodz,
// //            OrigVolFlag};
// // the strucutre of the data array is identical to the structure 
// // of the hitdata array defined in TSolEVIOFile class

class g4sbshitdata {
    public:
        //Default constructor. 
  	g4sbshitdata( int detid, unsigned int size = __DEFAULT_DATA_SIZE );
	virtual ~g4sbshitdata();
	
	//Get detector ID
	det_type  GetDetType() const { return fDetType;}
	int      GetDetID() const { return fDetID;}
	int     GetDetUniqueID() const { return ((int)fDetType)*10+fDetID;}

	// Get/set one specific element of the data for this hit
	void  SetData( unsigned int, double );
	double GetData( unsigned int ) const ;
	double *GetData(){ return fData; }//Get all data array 
	
	bool    IsFilled() const ;
	
    protected:
	det_type   fDetType;//detector type (see g4sbs_types.h)
	int         fDetID;//detector ID
	unsigned int fSize;//data array size;
	long long int fFillbits;
	double        *fData;//data array: See in .cxx the sequence of this data array for g4sbs GRINCH/RICH
};

////////////////////////////////////////////////////////////////////////////
// Auxilliary class for storing generated track data
// // ___________________________________________________________ //
// // gendata: {type, 
// //           pprodx, pprody, pprodz, 
// //           Xprodx, Xprody, Xprodz,
// //           OrigVolFlag, weight};
// // the strucutre of this data array is identical to the structure 
// // of the gendata array defined in TSolEVIOFile class

class g4sbsgendata : public g4sbshitdata {
 public:
  g4sbsgendata( int detid, unsigned int size = __DEFAULT_DATA_SIZE );
  ~g4sbsgendata(){;}
  
  int	GetSource() const { return IsFilled()? (int) fData[0] : -1e9; }//G4 particle ID
  int	GetTRID() const { return IsFilled()? (int) fData[1] : -1e9; }//G4 particle ID
  int	GetPID() const { return IsFilled()? (int) fData[2] : -1e9; }//G4 particle ID

  double GetX() const { return IsFilled()? fData[3] : -1e9; }
  double GetY() const { return IsFilled()? fData[4] : -1e9; }
  double GetT() const { return IsFilled()? fData[5] : -1e9; }
  double GetE() const { return IsFilled()? fData[6] : -1e9; }
  
  double  GetWeight() const { return IsFilled()? fData[7] : -1e9; }//cross section
  TVector3 GetP() const { return IsFilled()? TVector3(fData[8], fData[9], fData[10]) : TVector3(-1e9, -1e9, -1e9 ); }//Track momentum 3-vector
  TVector3 GetV() const { return IsFilled()? TVector3(fData[11], fData[12], fData[13]) : TVector3(-1e9, -1e9, -1e9 ); }//Track vtx 3-vector
  
  // This is from libsbsgem/TSBSGeant4File.h
  TVector3 GetMomentumAtTarget() const { return IsFilled()? TVector3(fData[14], fData[15], fData[16]) : TVector3(-1e9, -1e9, -1e9 ); }
  TVector3 GetVertexAtTarget() const { return IsFilled()? TVector3(fData[17], fData[18], fData[19]) : TVector3(-1e9, -1e9, -1e9 ); }
  
};

//
// Output data classes
//
class trackmchit_outdata{
 public: 
  trackmchit_outdata();
  ~trackmchit_outdata();

  UInt_t fNTrackMCHits;
  std::vector<Short_t>   fTrackMCSource;
  std::vector<Short_t>   fTrackMCTRID;
  std::vector<Int_t>     fTrackMCPID;
  std::vector<Double_t>  fTrackMCXhit;
  std::vector<Double_t>  fTrackMCYhit;
  std::vector<Double_t>  fTrackMCThit;
  std::vector<Double_t>  fTrackMCE;
  std::vector<Double_t>  fTrackMCWeight;
  std::vector<Double_t>  fTrackMCtrpx;
  std::vector<Double_t>  fTrackMCtrpy;
  std::vector<Double_t>  fTrackMCtrpz;
  std::vector<Double_t>  fTrackMCtrx;
  std::vector<Double_t>  fTrackMCtry;
  /*
  std::vector<Double_t>  fTrackMCtrz;
  std::vector<Double_t>  fTrackMCtrpx_v;
  std::vector<Double_t>  fTrackMCtrpy_v;
  std::vector<Double_t>  fTrackMCtrpz_v;
  std::vector<Double_t>  fTrackMCtrx_v;
  std::vector<Double_t>  fTrackMCtry_v;
  std::vector<Double_t>  fTrackMCtrz_v;
  */
  
  void Clear();
  bool CheckSize(bool print = false);
};

class simhitmc_outdata{
 public: 
  simhitmc_outdata();
  ~simhitmc_outdata();

  UInt_t fNSimHits;
  std::vector<Short_t>   fSimSource;
  std::vector<Short_t>   fSimTRID;
  std::vector<Int_t>     fSimPID;
  std::vector<Short_t>   fSimChannel;
  std::vector<Double_t>  fSimEdep;
  std::vector<Int_t>     fSimNpe;
  std::vector<Double_t>  fSimTime;
  std::vector<Double_t>  fSimLeadTime;
  std::vector<Double_t>  fSimTrailTime;
  
  void Clear();
  bool CheckSize(bool check_edep = true, 
		 bool check_npe = true, 
		 bool check_times = true, 
		 bool print = false);
};

class simgemhitmc_outdata: public simhitmc_outdata{
 public:
  simgemhitmc_outdata();
  ~simgemhitmc_outdata();
  
  std::vector<Short_t>   fPlane;
  std::vector<Short_t>   fModule;
  std::vector<Short_t>   fSizeX;
  std::vector<Short_t>   fSizeY;
  std::vector<Short_t>   fStartX;
  std::vector<Short_t>   fStartY;
  std::vector<Double_t>  fXpos;
  std::vector<Double_t>  fYpos;
  std::vector<Double_t>  fPX;
  std::vector<Double_t>  fPY;
  std::vector<Double_t>  fPZ;
  
  void Clear();
  bool CheckSize(bool print = false);
};



class simdig_outdata{
 public:
  simdig_outdata();
  ~simdig_outdata();
  
  UInt_t fNHits;
  std::vector<Short_t>  fChannel;
  std::vector<uint32_t> fDataWord;
  std::vector<Int_t>    fADC;
  std::vector<Int_t>    fTDC_L;
  std::vector<Int_t>    fTDC_T;
  
  void Clear();
  bool CheckSize(bool check_adc = true, 
		 bool check_tdc = true, 
		 bool print = false);
};

class simdigsamp_outdata: public simdig_outdata{
 public:
  simdigsamp_outdata();
  ~simdigsamp_outdata();
  
  //it's more complicated that I'd like, but that's the price to pay to have a reasonably efficient data storage for sample data / esp for GEMs -- EF
  // we have to assume elt n = samp n
  std::vector<uint32_t> fNsamps;
  std::vector< std::vector<uint32_t> > fDataWord_samps;
  std::vector< std::vector<Int_t> > fADC_samps; 
 
  void Clear();
  bool CheckSize(bool check_tdc = true, 
		 bool print = false);
};



#endif // TSBSSIMDATA_H
