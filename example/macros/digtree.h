//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov 13 14:40:36 2019 by ROOT version 6.08/00
// from TTree digtree/
// found on file: ../digitized/simdig_NoBkgd_2kevts.root
//////////////////////////////////////////////////////////

#ifndef digtree_h
#define digtree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
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
   Int_t           sbs_hcal_Nsimhits;
   vector<short>   *sbs_hcal_simhit_chan;
   vector<double>  *sbs_hcal_simhit_Edep;
   vector<int>     *sbs_hcal_simhit_npe;
   vector<double>  *sbs_hcal_simhit_time;
   vector<double>  *sbs_hcal_simhit_t_lead;
   vector<double>  *sbs_hcal_simhit_t_trail;
   Int_t           sbs_hcal_Nhits;
   vector<short>   *sbs_hcal_hit_chan;
   vector<unsigned int> *sbs_hcal_hit_dataword;
   Int_t           sbs_cdet_Nsimhits;
   vector<short>   *sbs_cdet_simhit_chan;
   vector<double>  *sbs_cdet_simhit_Edep;
   vector<int>     *sbs_cdet_simhit_npe;
   vector<double>  *sbs_cdet_simhit_time;
   vector<double>  *sbs_cdet_simhit_t_lead;
   vector<double>  *sbs_cdet_simhit_t_trail;
   Int_t           sbs_cdet_Nhits;
   vector<short>   *sbs_cdet_hit_chan;
   vector<unsigned int> *sbs_cdet_hit_dataword;
   Int_t           bb_sh_Nsimhits;
   vector<short>   *bb_sh_simhit_chan;
   vector<double>  *bb_sh_simhit_Edep;
   vector<int>     *bb_sh_simhit_npe;
   vector<double>  *bb_sh_simhit_time;
   Int_t           bb_sh_Nhits;
   vector<short>   *bb_sh_hit_chan;
   vector<unsigned int> *bb_sh_hit_dataword;
   Int_t           bb_ps_Nsimhits;
   vector<short>   *bb_ps_simhit_chan;
   vector<double>  *bb_ps_simhit_Edep;
   vector<int>     *bb_ps_simhit_npe;
   vector<double>  *bb_ps_simhit_time;
   Int_t           bb_ps_Nhits;
   vector<short>   *bb_ps_hit_chan;
   vector<unsigned int> *bb_ps_hit_dataword;
   Int_t           bb_hodo_Nsimhits;
   vector<short>   *bb_hodo_simhit_chan;
   vector<double>  *bb_hodo_simhit_Edep;
   vector<int>     *bb_hodo_simhit_npe;
   vector<double>  *bb_hodo_simhit_time;
   vector<double>  *bb_hodo_simhit_t_lead;
   vector<double>  *bb_hodo_simhit_t_trail;
   Int_t           bb_hodo_Nhits;
   vector<short>   *bb_hodo_hit_chan;
   vector<unsigned int> *bb_hodo_hit_dataword;
   Int_t           bb_grinch_Nsimhits;
   vector<short>   *bb_grinch_simhit_chan;
   vector<int>     *bb_grinch_simhit_npe;
   vector<double>  *bb_grinch_simhit_time;
   vector<double>  *bb_grinch_simhit_t_lead;
   vector<double>  *bb_grinch_simhit_t_trail;
   Int_t           bb_grinch_Nhits;
   vector<short>   *bb_grinch_hit_chan;
   vector<unsigned int> *bb_grinch_hit_dataword;

   // List of branches
   TBranch        *b_RunID;   //!
   TBranch        *b_EvtID;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_NSignal;   //!
   TBranch        *b_sbs_hcal_Nsimhits;   //!
   TBranch        *b_sbs_hcal_simhit_chan;   //!
   TBranch        *b_sbs_hcal_simhit_Edep;   //!
   TBranch        *b_sbs_hcal_simhit_npe;   //!
   TBranch        *b_sbs_hcal_simhit_time;   //!
   TBranch        *b_sbs_hcal_simhit_t_lead;   //!
   TBranch        *b_sbs_hcal_simhit_t_trail;   //!
   TBranch        *b_sbs_hcal_Nhits;   //!
   TBranch        *b_sbs_hcal_hit_chan;   //!
   TBranch        *b_sbs_hcal_hit_dataword;   //!
   TBranch        *b_sbs_cdet_Nsimhits;   //!
   TBranch        *b_sbs_cdet_simhit_chan;   //!
   TBranch        *b_sbs_cdet_simhit_Edep;   //!
   TBranch        *b_sbs_cdet_simhit_npe;   //!
   TBranch        *b_sbs_cdet_simhit_time;   //!
   TBranch        *b_sbs_cdet_simhit_t_lead;   //!
   TBranch        *b_sbs_cdet_simhit_t_trail;   //!
   TBranch        *b_sbs_cdet_Nhits;   //!
   TBranch        *b_sbs_cdet_hit_chan;   //!
   TBranch        *b_sbs_cdet_hit_dataword;   //!
   TBranch        *b_bb_sh_Nsimhits;   //!
   TBranch        *b_bb_sh_simhit_chan;   //!
   TBranch        *b_bb_sh_simhit_Edep;   //!
   TBranch        *b_bb_sh_simhit_npe;   //!
   TBranch        *b_bb_sh_simhit_time;   //!
   TBranch        *b_bb_sh_Nhits;   //!
   TBranch        *b_bb_sh_hit_chan;   //!
   TBranch        *b_bb_sh_hit_dataword;   //!
   TBranch        *b_bb_ps_Nsimhits;   //!
   TBranch        *b_bb_ps_simhit_chan;   //!
   TBranch        *b_bb_ps_simhit_Edep;   //!
   TBranch        *b_bb_ps_simhit_npe;   //!
   TBranch        *b_bb_ps_simhit_time;   //!
   TBranch        *b_bb_ps_Nhits;   //!
   TBranch        *b_bb_ps_hit_chan;   //!
   TBranch        *b_bb_ps_hit_dataword;   //!
   TBranch        *b_bb_hodo_Nsimhits;   //!
   TBranch        *b_bb_hodo_simhit_chan;   //!
   TBranch        *b_bb_hodo_simhit_Edep;   //!
   TBranch        *b_bb_hodo_simhit_npe;   //!
   TBranch        *b_bb_hodo_simhit_time;   //!
   TBranch        *b_bb_hodo_simhit_t_lead;   //!
   TBranch        *b_bb_hodo_simhit_t_trail;   //!
   TBranch        *b_bb_hodo_Nhits;   //!
   TBranch        *b_bb_hodo_hit_chan;   //!
   TBranch        *b_bb_hodo_hit_dataword;   //!
   TBranch        *b_bb_grinch_Nsimhits;   //!
   TBranch        *b_bb_grinch_simhit_chan;   //!
   TBranch        *b_bb_grinch_simhit_npe;   //!
   TBranch        *b_bb_grinch_simhit_time;   //!
   TBranch        *b_bb_grinch_simhit_t_lead;   //!
   TBranch        *b_bb_grinch_simhit_t_trail;   //!
   TBranch        *b_bb_grinch_Nhits;   //!
   TBranch        *b_bb_grinch_hit_chan;   //!
   TBranch        *b_bb_grinch_hit_dataword;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../digitized/simdig_NoBkgd_2kevts.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../digitized/simdig_NoBkgd_2kevts.root");
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
   sbs_hcal_simhit_chan = 0;
   sbs_hcal_simhit_Edep = 0;
   sbs_hcal_simhit_npe = 0;
   sbs_hcal_simhit_time = 0;
   sbs_hcal_simhit_t_lead = 0;
   sbs_hcal_simhit_t_trail = 0;
   sbs_hcal_hit_chan = 0;
   sbs_hcal_hit_dataword = 0;
   sbs_cdet_simhit_chan = 0;
   sbs_cdet_simhit_Edep = 0;
   sbs_cdet_simhit_npe = 0;
   sbs_cdet_simhit_time = 0;
   sbs_cdet_simhit_t_lead = 0;
   sbs_cdet_simhit_t_trail = 0;
   sbs_cdet_hit_chan = 0;
   sbs_cdet_hit_dataword = 0;
   bb_sh_simhit_chan = 0;
   bb_sh_simhit_Edep = 0;
   bb_sh_simhit_npe = 0;
   bb_sh_simhit_time = 0;
   bb_sh_hit_chan = 0;
   bb_sh_hit_dataword = 0;
   bb_ps_simhit_chan = 0;
   bb_ps_simhit_Edep = 0;
   bb_ps_simhit_npe = 0;
   bb_ps_simhit_time = 0;
   bb_ps_hit_chan = 0;
   bb_ps_hit_dataword = 0;
   bb_hodo_simhit_chan = 0;
   bb_hodo_simhit_Edep = 0;
   bb_hodo_simhit_npe = 0;
   bb_hodo_simhit_time = 0;
   bb_hodo_simhit_t_lead = 0;
   bb_hodo_simhit_t_trail = 0;
   bb_hodo_hit_chan = 0;
   bb_hodo_hit_dataword = 0;
   bb_grinch_simhit_chan = 0;
   bb_grinch_simhit_npe = 0;
   bb_grinch_simhit_time = 0;
   bb_grinch_simhit_t_lead = 0;
   bb_grinch_simhit_t_trail = 0;
   bb_grinch_hit_chan = 0;
   bb_grinch_hit_dataword = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunID", &RunID, &b_RunID);
   fChain->SetBranchAddress("EvtID", &EvtID, &b_EvtID);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("NSignal", &NSignal, &b_NSignal);
   fChain->SetBranchAddress("sbs.hcal_Nsimhits", &sbs_hcal_Nsimhits, &b_sbs_hcal_Nsimhits);
   fChain->SetBranchAddress("sbs.hcal_simhit_chan", &sbs_hcal_simhit_chan, &b_sbs_hcal_simhit_chan);
   fChain->SetBranchAddress("sbs.hcal_simhit_Edep", &sbs_hcal_simhit_Edep, &b_sbs_hcal_simhit_Edep);
   fChain->SetBranchAddress("sbs.hcal_simhit_npe", &sbs_hcal_simhit_npe, &b_sbs_hcal_simhit_npe);
   fChain->SetBranchAddress("sbs.hcal_simhit_time", &sbs_hcal_simhit_time, &b_sbs_hcal_simhit_time);
   fChain->SetBranchAddress("sbs.hcal_simhit_t_lead", &sbs_hcal_simhit_t_lead, &b_sbs_hcal_simhit_t_lead);
   fChain->SetBranchAddress("sbs.hcal_simhit_t_trail", &sbs_hcal_simhit_t_trail, &b_sbs_hcal_simhit_t_trail);
   fChain->SetBranchAddress("sbs.hcal_Nhits", &sbs_hcal_Nhits, &b_sbs_hcal_Nhits);
   fChain->SetBranchAddress("sbs.hcal_hit_chan", &sbs_hcal_hit_chan, &b_sbs_hcal_hit_chan);
   fChain->SetBranchAddress("sbs.hcal_hit_dataword", &sbs_hcal_hit_dataword, &b_sbs_hcal_hit_dataword);
   fChain->SetBranchAddress("sbs.cdet_Nsimhits", &sbs_cdet_Nsimhits, &b_sbs_cdet_Nsimhits);
   fChain->SetBranchAddress("sbs.cdet_simhit_chan", &sbs_cdet_simhit_chan, &b_sbs_cdet_simhit_chan);
   fChain->SetBranchAddress("sbs.cdet_simhit_Edep", &sbs_cdet_simhit_Edep, &b_sbs_cdet_simhit_Edep);
   fChain->SetBranchAddress("sbs.cdet_simhit_npe", &sbs_cdet_simhit_npe, &b_sbs_cdet_simhit_npe);
   fChain->SetBranchAddress("sbs.cdet_simhit_time", &sbs_cdet_simhit_time, &b_sbs_cdet_simhit_time);
   fChain->SetBranchAddress("sbs.cdet_simhit_t_lead", &sbs_cdet_simhit_t_lead, &b_sbs_cdet_simhit_t_lead);
   fChain->SetBranchAddress("sbs.cdet_simhit_t_trail", &sbs_cdet_simhit_t_trail, &b_sbs_cdet_simhit_t_trail);
   fChain->SetBranchAddress("sbs.cdet_Nhits", &sbs_cdet_Nhits, &b_sbs_cdet_Nhits);
   fChain->SetBranchAddress("sbs.cdet_hit_chan", &sbs_cdet_hit_chan, &b_sbs_cdet_hit_chan);
   fChain->SetBranchAddress("sbs.cdet_hit_dataword", &sbs_cdet_hit_dataword, &b_sbs_cdet_hit_dataword);
   fChain->SetBranchAddress("bb.sh_Nsimhits", &bb_sh_Nsimhits, &b_bb_sh_Nsimhits);
   fChain->SetBranchAddress("bb.sh_simhit_chan", &bb_sh_simhit_chan, &b_bb_sh_simhit_chan);
   fChain->SetBranchAddress("bb.sh_simhit_Edep", &bb_sh_simhit_Edep, &b_bb_sh_simhit_Edep);
   fChain->SetBranchAddress("bb.sh_simhit_npe", &bb_sh_simhit_npe, &b_bb_sh_simhit_npe);
   fChain->SetBranchAddress("bb.sh_simhit_time", &bb_sh_simhit_time, &b_bb_sh_simhit_time);
   fChain->SetBranchAddress("bb.sh_Nhits", &bb_sh_Nhits, &b_bb_sh_Nhits);
   fChain->SetBranchAddress("bb.sh_hit_chan", &bb_sh_hit_chan, &b_bb_sh_hit_chan);
   fChain->SetBranchAddress("bb.sh_hit_dataword", &bb_sh_hit_dataword, &b_bb_sh_hit_dataword);
   fChain->SetBranchAddress("bb.ps_Nsimhits", &bb_ps_Nsimhits, &b_bb_ps_Nsimhits);
   fChain->SetBranchAddress("bb.ps_simhit_chan", &bb_ps_simhit_chan, &b_bb_ps_simhit_chan);
   fChain->SetBranchAddress("bb.ps_simhit_Edep", &bb_ps_simhit_Edep, &b_bb_ps_simhit_Edep);
   fChain->SetBranchAddress("bb.ps_simhit_npe", &bb_ps_simhit_npe, &b_bb_ps_simhit_npe);
   fChain->SetBranchAddress("bb.ps_simhit_time", &bb_ps_simhit_time, &b_bb_ps_simhit_time);
   fChain->SetBranchAddress("bb.ps_Nhits", &bb_ps_Nhits, &b_bb_ps_Nhits);
   fChain->SetBranchAddress("bb.ps_hit_chan", &bb_ps_hit_chan, &b_bb_ps_hit_chan);
   fChain->SetBranchAddress("bb.ps_hit_dataword", &bb_ps_hit_dataword, &b_bb_ps_hit_dataword);
   fChain->SetBranchAddress("bb.hodo_Nsimhits", &bb_hodo_Nsimhits, &b_bb_hodo_Nsimhits);
   fChain->SetBranchAddress("bb.hodo_simhit_chan", &bb_hodo_simhit_chan, &b_bb_hodo_simhit_chan);
   fChain->SetBranchAddress("bb.hodo_simhit_Edep", &bb_hodo_simhit_Edep, &b_bb_hodo_simhit_Edep);
   fChain->SetBranchAddress("bb.hodo_simhit_npe", &bb_hodo_simhit_npe, &b_bb_hodo_simhit_npe);
   fChain->SetBranchAddress("bb.hodo_simhit_time", &bb_hodo_simhit_time, &b_bb_hodo_simhit_time);
   fChain->SetBranchAddress("bb.hodo_simhit_t_lead", &bb_hodo_simhit_t_lead, &b_bb_hodo_simhit_t_lead);
   fChain->SetBranchAddress("bb.hodo_simhit_t_trail", &bb_hodo_simhit_t_trail, &b_bb_hodo_simhit_t_trail);
   fChain->SetBranchAddress("bb.hodo_Nhits", &bb_hodo_Nhits, &b_bb_hodo_Nhits);
   fChain->SetBranchAddress("bb.hodo_hit_chan", &bb_hodo_hit_chan, &b_bb_hodo_hit_chan);
   fChain->SetBranchAddress("bb.hodo_hit_dataword", &bb_hodo_hit_dataword, &b_bb_hodo_hit_dataword);
   fChain->SetBranchAddress("bb.grinch_Nsimhits", &bb_grinch_Nsimhits, &b_bb_grinch_Nsimhits);
   fChain->SetBranchAddress("bb.grinch_simhit_chan", &bb_grinch_simhit_chan, &b_bb_grinch_simhit_chan);
   fChain->SetBranchAddress("bb.grinch_simhit_npe", &bb_grinch_simhit_npe, &b_bb_grinch_simhit_npe);
   fChain->SetBranchAddress("bb.grinch_simhit_time", &bb_grinch_simhit_time, &b_bb_grinch_simhit_time);
   fChain->SetBranchAddress("bb.grinch_simhit_t_lead", &bb_grinch_simhit_t_lead, &b_bb_grinch_simhit_t_lead);
   fChain->SetBranchAddress("bb.grinch_simhit_t_trail", &bb_grinch_simhit_t_trail, &b_bb_grinch_simhit_t_trail);
   fChain->SetBranchAddress("bb.grinch_Nhits", &bb_grinch_Nhits, &b_bb_grinch_Nhits);
   fChain->SetBranchAddress("bb.grinch_hit_chan", &bb_grinch_hit_chan, &b_bb_grinch_hit_chan);
   fChain->SetBranchAddress("bb.grinch_hit_dataword", &bb_grinch_hit_dataword, &b_bb_grinch_hit_dataword);
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
