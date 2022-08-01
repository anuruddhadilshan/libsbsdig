libsbsdig program:

The purpose of this prgram is to digitize the detectors output from G4SBS.
To make a local copy, type: 
```shell
git clone git@github.com:JeffersonLab/libsbsdig
```

# The most up-to-date version of the program is now in the master branch 

A complete documentation is available at:
https://redmine.jlab.org/projects/sbs-software/wiki/Documentation_of_libsbsdig

## Use:

```shell
sbsdig db_gmn_conf.dat gmn13.5_elastic_ex.txt 100000
```

First argument is the configuration file;
Second argument is the text file containg the list of root files to digitize;
Third argument is the number of events to digitize.


##List of sources/classes:







## List of classes

The list of classes and their functions are the following: 

**g4sbs_tree**

Unfolds the G4SBS file data tree; used by class TSBSGeant4File.

**TSBSGeant4File**

Reads the G4SBS output to extract the useful data for the digitization of the selected Cherenkov detector.
This data will be stored by an instance of classes TSBSDetData, 
and will be used by class TSBSSimCherDigitization.cxx.


**TSBSSimDigitizer**

Perfoms the digitization on all the list of detectors that have been added to it (via function TSBSSimDigitizer::AddDetector(TSBSSimDetector).

NB: TSBSSimDigitizer need to be added a procedure to superimpose e.g. background simulation on top of signal simulation.

**TSBSSimEvent**

Holds the data structure to fill the output file. 
The data structure is the following (to be completed as the library progresses).

```shell
  // Event identification
  Int_t     fRunID;               // Run number
  Int_t     fEvtID;               // Event number

  Double_t  fWeight;              // Event weight
  Int_t     fNSignal;             // Number of clusters from trigger track (signal)
  
  struct SimDetectorData
  struct DetectorData

```



**TSBSDBManager** (auxilliary class)
 
Useful class to unfold the databases and manage the parameters stored in them:
General parameters for DB (To be completed):

  //variable for data base information
  int fNDetectors;  // number of Cherenkov detectors in arm (usually 1...)
  int fChanPerSlot;  // number of PMTs per VETROC
  int fSlotPerCrate;  // number of VETROC per crate

**TSBSSimDecoder** (auxilliary class)

Interprets event buffer from input as TSBSSimEvent objects
(containing digitized simulation data) and unpacks them into
crateslot arrays for low-level decoding by detectors.

**TSBSSpec** (auxilliary class)

Spectrometer class which "carry" the detectors classes. 

**TSBSSimDetector** (auxilliary classes)

Generic detector class; contain the detector geometries. 

**TSBSSim{GEM/Cher/Scint/ECal/HCal}** inherit or will inherit from TSBSDet.
They will each have different geometry parameters.

```shell
TSBSGEMChamber/Plane:
to be merged (see TODO_1);
```


#### NOTE ABOUT RUNNING REPLAY SCRIPTS WITH ROOT6
May need to include a rootlogon.C in the current directory to get it
to pick up the SBS-offline stuff properly. It requires SBS_ANALYSIS
environmental variable to be set, but all in all it looks something
like this:
```cpp
R__LOAD_LIBRARY($SBS_ANALYSIS/libsbs.so)
R__LOAD_LIBRARY(../libsbsdig.so)

void rootlogon()
{
  // SBS-offline path
  gSystem->AddDynamicPath("${SBS_ANALYSIS}");
  gSystem->AddIncludePath(" -I${SBS_ANALYSIS}");
}
```
