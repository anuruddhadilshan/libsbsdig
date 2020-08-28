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
#include "SBSDigPMTDet.h"
#include "SBSDigGEMDet.h"

bool UnfoldData(gmn_tree* T, double theta_sbs, double d_hcal, TRandom3* R, 
		std::map<int, SBSDigPMTDet*> pmtdets, 
		std::map<int, SBSDigGEMDet*> gemdets, int signal);

/*
*/

#endif // SBSDIGAUXI_H

