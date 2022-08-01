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


##List of sources/classes (to be completed/updated):

**sbsdig.cxx**

Main program: Reads the config files, builds the list of detectors according to the config files, and processes the digitzation.

**g4sbs_tree**

Unfolds the G4SBS file data tree; uses g4sbs_data.

**g4sbs_data**
Contains the standard branch structure for each detector type. Used by g4sbs_tree.



