//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Dec  7 14:58:30 2014 by ROOT version 5.34/11
// from TTree tree/recbe
// found on file: ../root/run_000100_built.root
//////////////////////////////////////////////////////////

#ifndef simpleTohokuAna_h
#define simpleTohokuAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class simpleTohokuAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           triggerNumber;
   Int_t           adc[66][32];
   Int_t           q[66];
   Int_t           tdcNhit[66];
   Int_t           driftTime[66][32];
   Int_t           clockNumberDriftTime[66][32];

   // List of branches
   TBranch        *b_triggerNumber;   //!
   TBranch        *b_adc;   //!
   TBranch        *b_q;   //!
   TBranch        *b_tdcNhit;   //!
   TBranch        *b_driftTime;   //!
   TBranch        *b_clockNumberDriftTime;   //!

   simpleTohokuAna(TTree *tree=0);
   virtual ~simpleTohokuAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     plotWaveForm(Int_t wireID=0, Int_t firstEventNumber=0);
};

#endif

#ifdef simpleTohokuAna_cxx
simpleTohokuAna::simpleTohokuAna(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../root/run_000251_built.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../root/run_000251_built.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

simpleTohokuAna::~simpleTohokuAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t simpleTohokuAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t simpleTohokuAna::LoadTree(Long64_t entry)
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

void simpleTohokuAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("triggerNumber", &triggerNumber, &b_triggerNumber);
   fChain->SetBranchAddress("adc", adc, &b_adc);
   fChain->SetBranchAddress("q", q, &b_q);
   fChain->SetBranchAddress("tdcNhit", tdcNhit, &b_tdcNhit);
   fChain->SetBranchAddress("driftTime", driftTime, &b_driftTime);
   fChain->SetBranchAddress("clockNumberDriftTime", clockNumberDriftTime, &b_clockNumberDriftTime);
   Notify();
}

Bool_t simpleTohokuAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void simpleTohokuAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t simpleTohokuAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef simpleTohokuAna_cxx
