1libsbsdig library:

The purpose of this library is to digitize the detectors output from G4SBS.
To make a local copy, type: 
```shell
git clone git@github.com:efuchey/libsbsdig
```

Used in scripts such as those in directory "example".
NB: To be able to use the scripts in directory example, you need to link the databases you want to use there.


NB: so far, only Cherenkov digitization has been coded;
GEM digitization exists in another repository, which will have to be merged (TODO_1);
other subsystems have to be written:
Scintillators (TODO_2);
EM calorimeters (TODO_3);

## List of classes

The list of classes and their functions are the following: 

**g4sbs_tree**

Unfolds the G4SBS file data tree; used by class TSBSGeant4File.

**TSBSGeant4File**

Reads the G4SBS output to extract the useful data for the digitization of the selected Cherenkov detector.
This data will be stored by an instance of classes TSBSDetData, 
and will be used by class TSBSSimCherDigitization.cxx.

**TSBSDetData**

Generic library for the library internal data container;

**TSBS{GEM/Cher/Scint/ECal/HCal}Data** inherits or will inherit from TSBSDetData.
They each use a specific data structure.
They should use a standard container however, which is not the case yet (TODO).

The structure of the data stored is the following:

```shell
GEM:
to be merged (see TODO_1);
```

```shell
Cher:
    UInt_t    fDetID;          // Hit detector ID
    UInt_t    fPMTID;          // Hit PMT
    Double_t  fXPMT;           // X coordinate of the PMT in transport coordinates
    Double_t  fYPMT;           // Y coordinate of the PMT in transport coordinates
    Double_t  fNpe;            // Number of photoelectrons
    Double_t  fTime;           // hit time
    Double_t  fTimeRMS;        // hit time RMS: should be small (<0.1 ns)
    UInt_t    fType;           // particle type: 1 = primary, >1 secondaries
    TVector3  fPos_det;        // position of the optical photon at detection
    Int_t     fPID_MCtrack;    // Particle ID of the MC particle which produced the hit (*)
    TVector3  fMom_MCtrack;    // momentum of the MC particle which produced the hit (*)
    TVector3  fVtx_MCtrack;    // vertex of the MC particle which produced the hit (*)
    UInt_t    fOriginVolFlag;  // origin volume of the photon
    //  (*) valid if identified
Note the first 5 parameters are essential for the digitization itself,
the other are just useful complementary MC information.
```

```shell
Scint:
to be written (see TODO_2);
```

```shell
ECal:
to be written (see TODO_3);
```

**TSBSSimDetDigitization**

Generic class to perform the digitization of the subsystem data. 

**TSBSSim{GEM/Cher/Scint/ECal/HCal}Digitization** inherit or will inherit from TSBSSimDetDigitization.
They each use specific digitization database parameters and produce a specific output.

It produces the following output:

```shell
GEM:
to be merged (see TODO_1);
```

```shell
Cher:
  std::vector< std::pair<TArrayI, TArrayI> > fTDCArrays;
  // TDC digital value output arrays (1st array: rise time, 2nd array: fall time)
  std::vector< std::pair<TArrayD, TArrayD> > fTDCtimeArrays;
  // TDC "input" time values (evaluated with the PMT pulse shape)
```

```shell
Scint:
to be written (see TODO_2);
```

```shell
ECal:
to be written (see TODO_3);
```

These values are then fed to an instance of class TSBSSimEvent, 
which holds the data structure to fill the output file.

**TSBSSimEvent**

Holds the data structure to fill the output file. 
The data structure is the following (to be completed as the library progresses).

```shell
  // Event identification
  Int_t     fRunID;               // Run number
  Int_t     fEvtID;               // Event number

  Double_t  fWeight;              // Event weight
  Int_t     fNSignal;             // Number of clusters from trigger track (signal)
  
  // MC tracks
  TClonesArray*   fMCTracks;      //-> Physics tracks
  //std::vector< std::pair< Int_t, std::vector<Short_t> > > fMCClusterHitID;
  // collections of hit IDs coming from the same MC track ID
  // save only clusters with size >= 2

  struct PMTHit {
    // MC hit data
    Short_t  fID;          // Hit number
    Int_t    fSource;      // MC data set source (0 = signal, >0 background)
    Char_t   fType;        // GEANT particle type (1 = primary)
    Int_t    fMCtrackID;   // GEANT track ID
    Int_t    fMCtrackPID;  // GEANT particle ID (defined only if particle is primary)
    Short_t  fOrigVolFlag; // 
    Float_t  fXPMT;        // X coordinate of the PMT in transport coordinates
    Float_t  fYPMT;        // Y coordinate of the PMT in transport coordinates
    Float_t  fNpe;         // Number of photoelectrons
    Double_t fTime;        // Arrival time at electronics
    Double_t fTDCtime[2];  // TDC rising and falling times values wrt trigger time
    // Digitization results for this hit
    Short_t  fDetID;       // Detector ID
    Short_t  fChannel;     // Channel number
    Short_t  fPMTrow;      // Row number: cross reference to Channel number
    Short_t  fPMTcol;      // Column number: cross reference to Channel number
    Char_t   fVETROCID;    // VETROC ID
    uint32_t fTDC[2];      // TDC VETROC words
  };
  std::vector<PMTHit> fPMTHits;  // All MC-generated clusters in the GEMs
```



**TSBSDBManager** (auxilliary class)
 
Useful class to unfold the databases and manage the parameters stored in them:
General parameters for DB (To be completed):

  //variable for data base information
  int fNDetectors;  // number of Cherenkov detectors in arm (usually 1...)
  int fChanPerSlot;  // number of PMTs per VETROC
  int fSlotPerCrate;  // number of VETROC per crate
  
  // Parameters for TSBSGeant4File
  int fg4sbsExpType;// experiment flag. Choices are
  // 1 - GMn;
  // 2 - GEn;
  // 3 - GEp;
  // 4 - SIDIS;
  // 5 - A1n;
  int fg4sbsDetType;// flag to determine which type of GEM should be read. Choices are:
  // 1 - GRINCH
  // 2 - RICH
  
  // Parameters for signal particles
  int fNSigParticle; // number of signal particles
  vector<int>    fSigPID;
  vector<int>    fSigTID;
  
struct GeoInfo{
  int    fZCkovIn;       // Z of the entrance window in the spectrometer central ray;
  int    fNradiator;     // radiator index of refraction;
  int    fLradiator;     // radiator length on central ray;
  //int    fNquartz;       // quartz window index of refraction;
  int    fNPMTs;         // number of PMTs
  int    fNPMTrows;      // number of PMT rows
  int    fNPMTcolsMax;   // max number of PMT columns 
  double fPMTmatrixHext; // horizontal extension, in m, of the PMT matrix (from lower PMT center to higher PMT center)
  double fPMTmatrixVext; // vertical extension, in m, of the PMT matrix (from left PMT center to right PMT center)
  double fPMTdistX;      // projected X distance between the center of 2 PMT tubes in consecutive rows, in m
  double fPMTdistY;      // Y distance between the center of 2 PMT tubes in consecutive columns, in m
  double fX_TCPMT;       // X position of the top close PMT center in the PMT matrix (transport coord)
  double fY_TCPMT;       // Y position of the top close PMT center in the PMT matrix (transport coord)
};


**TSBSSimDecoder** (auxilliary class)

Interprets event buffer from input as TSBSSimEvent objects
(containing digitized simulation data) and unpacks them into
crateslot arrays for low-level decoding by detectors.

**TSBSSpec** (auxilliary class)

Spectrometer class which "carry" the detectors classes. 

**TSBSDet** (auxilliary classes)

Generic detector class; contain the detector geometries. 

**TSBS{GEM/Cher/Scint/ECal/HCal}** inherit or will inherit from TSBSDet.
They will each have different geometry parameters.

```shell
TSBSGEMChamber/Plane:
to be merged (see TODO_1);
```

```shell
TSBSCher:
    int    fNPMTs;         // number of PMTs
    int    fNPMTrows;      // number of PMT rows
    int    fNPMTcolsMax;   // max number of PMT columns 
    double fPMTmatrixHext; // horizontal extension, in m, of the PMT matrix (from lower PMT center to higher PMT center)
    double fPMTmatrixVext; // vertical extension, in m, of the PMT matrix (from left PMT center to right PMT center)
    double fPMTdistX;      // projected X distance between the center of 2 PMT tubes in consecutive rows, in m
    double fPMTdistY;      // Y distance between the center of 2 PMT tubes in consecutive columns, in m
    double fX_TCPMT;       // X position of the top close PMT center in the PMT matrix (transport coord)
    double fY_TCPMT;       // Y position of the top close PMT center in the PMT matrix (transport coord)
```

```shell
TSBSScint:
to be written (see TODO_2);
```

```shell
TSBSECal:
to be written (see TODO_3);
```


