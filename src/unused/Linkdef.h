#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ defined_in "src/g4sbs_data.h";
#pragma link C++ defined_in "src/g4sbs_tree.h";
#pragma link C++ defined_in "src/TSBSDBManager.h";
#pragma link C++ defined_in "src/TSBSGeant4File.h";
#pragma link C++ defined_in "src/TSBSSimAuxi.h";
#pragma link C++ defined_in "src/TSBSSimEvent.h";
#pragma link C++ defined_in "src/TSBSSimData.h";
#pragma link C++ defined_in "src/TSBSSimDigitizer.h";
#pragma link C++ defined_in "src/TSBSSimDetector.h";
#pragma link C++ defined_in "src/TSBSSimCher.h";
#pragma link C++ defined_in "src/TSBSSimECal.h";
#pragma link C++ defined_in "src/TSBSSimHCal.h";
#pragma link C++ defined_in "src/TSBSSimGEM.h";
#pragma link C++ defined_in "src/TSBSSimScint.h";
#pragma link C++ defined_in "src/TSBSSpec.h";
//#pragma link C++ defined_in "src/TSBSSimDecoder.h";
//#pragma link C++ defined_in "src/digsim_data.h";
//#pragma link C++ defined_in "src/digsim_tree.h";
//#pragma link C++ defined_in "src/TSBSSimFile.h";
//#pragma link C++ defined_in "src/TSBSSimADC.h";
//#pragma link C++ defined_in "src/TSBSSimTDC.h";
//#pragma link C++ defined_in "src/TSBSSimMPD.h";
//#pragma link C++ defined_in "src/TSBSSimDataEncoder.h";

// Then include all the TGEM stuff
#pragma link C++ defined_in "libsbsgem/TGEMSBSBox.h";
#pragma link C++ defined_in "libsbsgem/TGEMSBSDBManager.h";
#pragma link C++ defined_in "libsbsgem/TGEMSBSGEMChamber.h";
#pragma link C++ defined_in "libsbsgem/TGEMSBSGEMHit.h";
#pragma link C++ defined_in "libsbsgem/TGEMSBSGEMPlane.h";
#pragma link C++ defined_in "libsbsgem/TGEMSBSGEMSimHitData.h";
#pragma link C++ defined_in "libsbsgem/TGEMSBSSimAuxi.h";
#pragma link C++ defined_in "libsbsgem/TGEMSBSSimDigitization.h";
#pragma link C++ defined_in "libsbsgem/TGEMSBSSpec.h";
//#pragma link C++ defined_in "libsbsgem/TGEMSBSGeant4File.h";
//#pragma link C++ defined_in "libsbsgem/TGEMSBSSimDecoder.h";
//#pragma link C++ defined_in "libsbsgem/TGEMSBSSimEvent.h";
//#pragma link C++ defined_in "libsbsgem/TGEMSBSSimFile.h";

// Limited stuff in EVIO file.  We don't want to be
// able to call the evio functions in the interpreter

#endif
