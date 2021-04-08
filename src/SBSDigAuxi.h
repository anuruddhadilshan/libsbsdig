#ifndef SBSDIGAUXI_H
#define SBSDIGAUXI_H

#include <iostream>
#include <vector>
#include <map>
#include "g4sbs_tree.h"
#include "TF1.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "SBSDigPMTDet.h"
#include "SBSDigGEMDet.h"

bool UnfoldData(g4sbs_tree* T, double theta_sbs, double d_hcal, TRandom3* R, 
		std::vector<SBSDigPMTDet*> pmtdets, 
		std::vector<int> detmap, 
		std::vector<SBSDigGEMDet*> gemdets, 
		std::vector<int> gemmap, 
		double tzero,
		int signal);


#endif // SBSDIGAUXI_H

