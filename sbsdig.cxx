//includes

#include <iostream>
#include <math>

#ifdef __APPLE__
#include "unistd.h"
#endif


//____________________________________________________
int main(int argc, char** argv){
  string inputsigfile, inputbkgdfile = "";
  ULong64_t nentries = -1;
  UShort nbkgd = 0;
  
  if(argc<2){
    cout << "*** Not enough arguments! ***" << endl
	 << " Arguments: list_of_sig_input_files (str, mandatory); " << endl
	 << "           nb_of_sig_evts_to_process (int, def=-1); " << endl
	 << "          list_of_bkgd_input_files (str, def=''); " << endl
	 << "         nb_of_bkgd_files_to_add_to_sig_evt (int, def=0); " << endl;
    return(-1);
  }
  
  inputsigfile = argv[1];
  cout << " Signal input files from: " << inputsigfile << endl;
  if(argc>2)nentries = atoi(argv[2]);
  cout << " Number of (signal) events to process = " << nentries << endl;
  if(argc>4){
    inputbkgdfile = argv[3];
    cout << " Background input files from: " << inputbkgdfile << endl;
    nbkgd = atoi(argv[4]);
    cout << " Number of background files to superimpose to signal = " << nbkgd << endl;
  }
  
  
  
}

