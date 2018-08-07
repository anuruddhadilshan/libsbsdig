#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ defined_in "src/g4sbs_tree.h";
#pragma link C++ defined_in "src/g4sbs_data.h";
#pragma link C++ defined_in "src/TSBSDBManager.h";
#pragma link C++ defined_in "src/TSBSGeant4File.h";
#pragma link C++ defined_in "src/TSBSSimAuxi.h";
#pragma link C++ defined_in "src/TSBSSimDecoder.h";
#pragma link C++ defined_in "src/TSBSSimEvent.h";
#pragma link C++ defined_in "src/TSBSSimFile.h";
#pragma link C++ defined_in "src/TSBSSimData.h";
#pragma link C++ defined_in "src/TSBSSimDigitizer.h";
#pragma link C++ defined_in "src/TSBSSimDetector.h";
#pragma link C++ defined_in "src/TSBSSimCher.h";
#pragma link C++ defined_in "src/TSBSSimHCal.h";
#pragma link C++ defined_in "src/TSBSSimScint.h";
#pragma link C++ defined_in "src/TSBSSpec.h";

//#pragma link C++ defined_in "src/TSBSCher.h";
//#pragma link C++ defined_in "src/TSBSCherData.h";
//#pragma link C++ defined_in "src/TSBSDet.h";
//#pragma link C++ defined_in "src/TSBSDetData.h";
//#pragma link C++ defined_in "src/TSBSSimCherDigitization.h";
//#pragma link C++ defined_in "src/TSBSSimDetDigitization.h";

// Limited stuff in EVIO file.  We don't want to be
// able to call the evio functions in the interpreter

#endif
