//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Sep 20 15:41:09 2019 by ROOT version 6.08/00
// from TTree digtree/
// found on file: ../digitized/simdig_test.root
//////////////////////////////////////////////////////////

#ifndef digtree_h
#define digtree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class digtree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           RunID;
   Int_t           EvtID;
   Double_t        Weight;
   Int_t           NSignal;
   Int_t           NSimData_sbs_hcal;
   vector<short>   *SimData_sbs_hcal_Chan;
   vector<short>   *SimData_sbs_hcal_Type;
   vector<short>   *SimData_sbs_hcal_Ndata;
   vector<vector<double> > *SimData_sbs_hcal_Data;
   Int_t           NData_sbs_hcal;
   vector<short>   *Data_sbs_hcal_Chan;
   vector<short>   *Data_sbs_hcal_Ndata;
   vector<vector<unsigned int> > *Data_sbs_hcal_Data;
   Int_t           NSimData_sbs_cdet;
   vector<short>   *SimData_sbs_cdet_Chan;
   vector<short>   *SimData_sbs_cdet_Type;
   vector<short>   *SimData_sbs_cdet_Ndata;
   vector<vector<double> > *SimData_sbs_cdet_Data;
   Int_t           NData_sbs_cdet;
   vector<short>   *Data_sbs_cdet_Chan;
   vector<short>   *Data_sbs_cdet_Ndata;
   vector<vector<unsigned int> > *Data_sbs_cdet_Data;
   Int_t           NSimData_bb_sh;
   vector<short>   *SimData_bb_sh_Chan;
   vector<short>   *SimData_bb_sh_Type;
   vector<short>   *SimData_bb_sh_Ndata;
   vector<vector<double> > *SimData_bb_sh_Data;
   Int_t           NData_bb_sh;
   vector<short>   *Data_bb_sh_Chan;
   vector<short>   *Data_bb_sh_Ndata;
   vector<vector<unsigned int> > *Data_bb_sh_Data;
   Int_t           NSimData_bb_ps;
   vector<short>   *SimData_bb_ps_Chan;
   vector<short>   *SimData_bb_ps_Type;
   vector<short>   *SimData_bb_ps_Ndata;
   vector<vector<double> > *SimData_bb_ps_Data;
   Int_t           NData_bb_ps;
   vector<short>   *Data_bb_ps_Chan;
   vector<short>   *Data_bb_ps_Ndata;
   vector<vector<unsigned int> > *Data_bb_ps_Data;
   Int_t           NSimData_bb_hodo;
   vector<short>   *SimData_bb_hodo_Chan;
   vector<short>   *SimData_bb_hodo_Type;
   vector<short>   *SimData_bb_hodo_Ndata;
   vector<vector<double> > *SimData_bb_hodo_Data;
   Int_t           NData_bb_hodo;
   vector<short>   *Data_bb_hodo_Chan;
   vector<short>   *Data_bb_hodo_Ndata;
   vector<vector<unsigned int> > *Data_bb_hodo_Data;
   Int_t           NSimData_bb_grinch;
   vector<short>   *SimData_bb_grinch_Chan;
   vector<short>   *SimData_bb_grinch_Type;
   vector<short>   *SimData_bb_grinch_Ndata;
   vector<vector<double> > *SimData_bb_grinch_Data;
   Int_t           NData_bb_grinch;
   vector<short>   *Data_bb_grinch_Chan;
   vector<short>   *Data_bb_grinch_Ndata;
   vector<vector<unsigned int> > *Data_bb_grinch_Data;

   // List of branches
   TBranch        *b_RunID;   //!
   TBranch        *b_EvtID;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_NSignal;   //!
   TBranch        *b_NSimData_sbs_hcal;   //!
   TBranch        *b_SimData_sbs_hcal_Chan;   //!
   TBranch        *b_SimData_sbs_hcal_Type;   //!
   TBranch        *b_SimData_sbs_hcal_Ndata;   //!
   TBranch        *b_SimData_sbs_hcal_Data;   //!
   TBranch        *b_NData_sbs_hcal;   //!
   TBranch        *b_Data_sbs_hcal_Chan;   //!
   TBranch        *b_Data_sbs_hcal_Ndata;   //!
   TBranch        *b_Data_sbs_hcal_Data;   //!
   TBranch        *b_NSimData_sbs_cdet;   //!
   TBranch        *b_SimData_sbs_cdet_Chan;   //!
   TBranch        *b_SimData_sbs_cdet_Type;   //!
   TBranch        *b_SimData_sbs_cdet_Ndata;   //!
   TBranch        *b_SimData_sbs_cdet_Data;   //!
   TBranch        *b_NData_sbs_cdet;   //!
   TBranch        *b_Data_sbs_cdet_Chan;   //!
   TBranch        *b_Data_sbs_cdet_Ndata;   //!
   TBranch        *b_Data_sbs_cdet_Data;   //!
   TBranch        *b_NSimData_bb_sh;   //!
   TBranch        *b_SimData_bb_sh_Chan;   //!
   TBranch        *b_SimData_bb_sh_Type;   //!
   TBranch        *b_SimData_bb_sh_Ndata;   //!
   TBranch        *b_SimData_bb_sh_Data;   //!
   TBranch        *b_NData_bb_sh;   //!
   TBranch        *b_Data_bb_sh_Chan;   //!
   TBranch        *b_Data_bb_sh_Ndata;   //!
   TBranch        *b_Data_bb_sh_Data;   //!
   TBranch        *b_NSimData_bb_ps;   //!
   TBranch        *b_SimData_bb_ps_Chan;   //!
   TBranch        *b_SimData_bb_ps_Type;   //!
   TBranch        *b_SimData_bb_ps_Ndata;   //!
   TBranch        *b_SimData_bb_ps_Data;   //!
   TBranch        *b_NData_bb_ps;   //!
   TBranch        *b_Data_bb_ps_Chan;   //!
   TBranch        *b_Data_bb_ps_Ndata;   //!
   TBranch        *b_Data_bb_ps_Data;   //!
   TBranch        *b_NSimData_bb_hodo;   //!
   TBranch        *b_SimData_bb_hodo_Chan;   //!
   TBranch        *b_SimData_bb_hodo_Type;   //!
   TBranch        *b_SimData_bb_hodo_Ndata;   //!
   TBranch        *b_SimData_bb_hodo_Data;   //!
   TBranch        *b_NData_bb_hodo;   //!
   TBranch        *b_Data_bb_hodo_Chan;   //!
   TBranch        *b_Data_bb_hodo_Ndata;   //!
   TBranch        *b_Data_bb_hodo_Data;   //!
   TBranch        *b_NSimData_bb_grinch;   //!
   TBranch        *b_SimData_bb_grinch_Chan;   //!
   TBranch        *b_SimData_bb_grinch_Type;   //!
   TBranch        *b_SimData_bb_grinch_Ndata;   //!
   TBranch        *b_SimData_bb_grinch_Data;   //!
   TBranch        *b_NData_bb_grinch;   //!
   TBranch        *b_Data_bb_grinch_Chan;   //!
   TBranch        *b_Data_bb_grinch_Ndata;   //!
   TBranch        *b_Data_bb_grinch_Data;   //!

   digtree(TTree *tree=0);
   virtual ~digtree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef digtree_cxx
digtree::digtree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../digitized/simdig_test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../digitized/simdig_test.root");
      }
      f->GetObject("digtree",tree);

   }
   Init(tree);
}

digtree::~digtree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t digtree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t digtree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void digtree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   SimData_sbs_hcal_Chan = 0;
   SimData_sbs_hcal_Type = 0;
   SimData_sbs_hcal_Ndata = 0;
   SimData_sbs_hcal_Data = 0;
   Data_sbs_hcal_Chan = 0;
   Data_sbs_hcal_Ndata = 0;
   Data_sbs_hcal_Data = 0;
   SimData_sbs_cdet_Chan = 0;
   SimData_sbs_cdet_Type = 0;
   SimData_sbs_cdet_Ndata = 0;
   SimData_sbs_cdet_Data = 0;
   Data_sbs_cdet_Chan = 0;
   Data_sbs_cdet_Ndata = 0;
   Data_sbs_cdet_Data = 0;
   SimData_bb_sh_Chan = 0;
   SimData_bb_sh_Type = 0;
   SimData_bb_sh_Ndata = 0;
   SimData_bb_sh_Data = 0;
   Data_bb_sh_Chan = 0;
   Data_bb_sh_Ndata = 0;
   Data_bb_sh_Data = 0;
   SimData_bb_ps_Chan = 0;
   SimData_bb_ps_Type = 0;
   SimData_bb_ps_Ndata = 0;
   SimData_bb_ps_Data = 0;
   Data_bb_ps_Chan = 0;
   Data_bb_ps_Ndata = 0;
   Data_bb_ps_Data = 0;
   SimData_bb_hodo_Chan = 0;
   SimData_bb_hodo_Type = 0;
   SimData_bb_hodo_Ndata = 0;
   SimData_bb_hodo_Data = 0;
   Data_bb_hodo_Chan = 0;
   Data_bb_hodo_Ndata = 0;
   Data_bb_hodo_Data = 0;
   SimData_bb_grinch_Chan = 0;
   SimData_bb_grinch_Type = 0;
   SimData_bb_grinch_Ndata = 0;
   SimData_bb_grinch_Data = 0;
   Data_bb_grinch_Chan = 0;
   Data_bb_grinch_Ndata = 0;
   Data_bb_grinch_Data = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunID", &RunID, &b_RunID);
   fChain->SetBranchAddress("EvtID", &EvtID, &b_EvtID);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("NSignal", &NSignal, &b_NSignal);
   fChain->SetBranchAddress("NSimData_sbs.hcal", &NSimData_sbs_hcal, &b_NSimData_sbs_hcal);
   fChain->SetBranchAddress("SimData_sbs.hcal_Chan", &SimData_sbs_hcal_Chan, &b_SimData_sbs_hcal_Chan);
   fChain->SetBranchAddress("SimData_sbs.hcal_Type", &SimData_sbs_hcal_Type, &b_SimData_sbs_hcal_Type);
   fChain->SetBranchAddress("SimData_sbs.hcal_Ndata", &SimData_sbs_hcal_Ndata, &b_SimData_sbs_hcal_Ndata);
   fChain->SetBranchAddress("SimData_sbs.hcal_Data", &SimData_sbs_hcal_Data, &b_SimData_sbs_hcal_Data);
   fChain->SetBranchAddress("NData_sbs.hcal", &NData_sbs_hcal, &b_NData_sbs_hcal);
   fChain->SetBranchAddress("Data_sbs.hcal_Chan", &Data_sbs_hcal_Chan, &b_Data_sbs_hcal_Chan);
   fChain->SetBranchAddress("Data_sbs.hcal_Ndata", &Data_sbs_hcal_Ndata, &b_Data_sbs_hcal_Ndata);
   fChain->SetBranchAddress("Data_sbs.hcal_Data", &Data_sbs_hcal_Data, &b_Data_sbs_hcal_Data);
   fChain->SetBranchAddress("NSimData_sbs.cdet", &NSimData_sbs_cdet, &b_NSimData_sbs_cdet);
   fChain->SetBranchAddress("SimData_sbs.cdet_Chan", &SimData_sbs_cdet_Chan, &b_SimData_sbs_cdet_Chan);
   fChain->SetBranchAddress("SimData_sbs.cdet_Type", &SimData_sbs_cdet_Type, &b_SimData_sbs_cdet_Type);
   fChain->SetBranchAddress("SimData_sbs.cdet_Ndata", &SimData_sbs_cdet_Ndata, &b_SimData_sbs_cdet_Ndata);
   fChain->SetBranchAddress("SimData_sbs.cdet_Data", &SimData_sbs_cdet_Data, &b_SimData_sbs_cdet_Data);
   fChain->SetBranchAddress("NData_sbs.cdet", &NData_sbs_cdet, &b_NData_sbs_cdet);
   fChain->SetBranchAddress("Data_sbs.cdet_Chan", &Data_sbs_cdet_Chan, &b_Data_sbs_cdet_Chan);
   fChain->SetBranchAddress("Data_sbs.cdet_Ndata", &Data_sbs_cdet_Ndata, &b_Data_sbs_cdet_Ndata);
   fChain->SetBranchAddress("Data_sbs.cdet_Data", &Data_sbs_cdet_Data, &b_Data_sbs_cdet_Data);
   fChain->SetBranchAddress("NSimData_bb.sh", &NSimData_bb_sh, &b_NSimData_bb_sh);
   fChain->SetBranchAddress("SimData_bb.sh_Chan", &SimData_bb_sh_Chan, &b_SimData_bb_sh_Chan);
   fChain->SetBranchAddress("SimData_bb.sh_Type", &SimData_bb_sh_Type, &b_SimData_bb_sh_Type);
   fChain->SetBranchAddress("SimData_bb.sh_Ndata", &SimData_bb_sh_Ndata, &b_SimData_bb_sh_Ndata);
   fChain->SetBranchAddress("SimData_bb.sh_Data", &SimData_bb_sh_Data, &b_SimData_bb_sh_Data);
   fChain->SetBranchAddress("NData_bb.sh", &NData_bb_sh, &b_NData_bb_sh);
   fChain->SetBranchAddress("Data_bb.sh_Chan", &Data_bb_sh_Chan, &b_Data_bb_sh_Chan);
   fChain->SetBranchAddress("Data_bb.sh_Ndata", &Data_bb_sh_Ndata, &b_Data_bb_sh_Ndata);
   fChain->SetBranchAddress("Data_bb.sh_Data", &Data_bb_sh_Data, &b_Data_bb_sh_Data);
   fChain->SetBranchAddress("NSimData_bb.ps", &NSimData_bb_ps, &b_NSimData_bb_ps);
   fChain->SetBranchAddress("SimData_bb.ps_Chan", &SimData_bb_ps_Chan, &b_SimData_bb_ps_Chan);
   fChain->SetBranchAddress("SimData_bb.ps_Type", &SimData_bb_ps_Type, &b_SimData_bb_ps_Type);
   fChain->SetBranchAddress("SimData_bb.ps_Ndata", &SimData_bb_ps_Ndata, &b_SimData_bb_ps_Ndata);
   fChain->SetBranchAddress("SimData_bb.ps_Data", &SimData_bb_ps_Data, &b_SimData_bb_ps_Data);
   fChain->SetBranchAddress("NData_bb.ps", &NData_bb_ps, &b_NData_bb_ps);
   fChain->SetBranchAddress("Data_bb.ps_Chan", &Data_bb_ps_Chan, &b_Data_bb_ps_Chan);
   fChain->SetBranchAddress("Data_bb.ps_Ndata", &Data_bb_ps_Ndata, &b_Data_bb_ps_Ndata);
   fChain->SetBranchAddress("Data_bb.ps_Data", &Data_bb_ps_Data, &b_Data_bb_ps_Data);
   fChain->SetBranchAddress("NSimData_bb.hodo", &NSimData_bb_hodo, &b_NSimData_bb_hodo);
   fChain->SetBranchAddress("SimData_bb.hodo_Chan", &SimData_bb_hodo_Chan, &b_SimData_bb_hodo_Chan);
   fChain->SetBranchAddress("SimData_bb.hodo_Type", &SimData_bb_hodo_Type, &b_SimData_bb_hodo_Type);
   fChain->SetBranchAddress("SimData_bb.hodo_Ndata", &SimData_bb_hodo_Ndata, &b_SimData_bb_hodo_Ndata);
   fChain->SetBranchAddress("SimData_bb.hodo_Data", &SimData_bb_hodo_Data, &b_SimData_bb_hodo_Data);
   fChain->SetBranchAddress("NData_bb.hodo", &NData_bb_hodo, &b_NData_bb_hodo);
   fChain->SetBranchAddress("Data_bb.hodo_Chan", &Data_bb_hodo_Chan, &b_Data_bb_hodo_Chan);
   fChain->SetBranchAddress("Data_bb.hodo_Ndata", &Data_bb_hodo_Ndata, &b_Data_bb_hodo_Ndata);
   fChain->SetBranchAddress("Data_bb.hodo_Data", &Data_bb_hodo_Data, &b_Data_bb_hodo_Data);
   fChain->SetBranchAddress("NSimData_bb.grinch", &NSimData_bb_grinch, &b_NSimData_bb_grinch);
   fChain->SetBranchAddress("SimData_bb.grinch_Chan", &SimData_bb_grinch_Chan, &b_SimData_bb_grinch_Chan);
   fChain->SetBranchAddress("SimData_bb.grinch_Type", &SimData_bb_grinch_Type, &b_SimData_bb_grinch_Type);
   fChain->SetBranchAddress("SimData_bb.grinch_Ndata", &SimData_bb_grinch_Ndata, &b_SimData_bb_grinch_Ndata);
   fChain->SetBranchAddress("SimData_bb.grinch_Data", &SimData_bb_grinch_Data, &b_SimData_bb_grinch_Data);
   fChain->SetBranchAddress("NData_bb.grinch", &NData_bb_grinch, &b_NData_bb_grinch);
   fChain->SetBranchAddress("Data_bb.grinch_Chan", &Data_bb_grinch_Chan, &b_Data_bb_grinch_Chan);
   fChain->SetBranchAddress("Data_bb.grinch_Ndata", &Data_bb_grinch_Ndata, &b_Data_bb_grinch_Ndata);
   fChain->SetBranchAddress("Data_bb.grinch_Data", &Data_bb_grinch_Data, &b_Data_bb_grinch_Data);
   Notify();
}

Bool_t digtree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void digtree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t digtree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef digtree_cxx
