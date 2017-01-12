//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Dec  6 20:09:02 2014 by ROOT version 5.34/11
// from TTree tree/recbe
// found on file: ../root/run_000079.event.root
//////////////////////////////////////////////////////////

#ifndef hoge_h
#define hoge_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class hoge {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           iev;
   Int_t           nSamples;
   Int_t           t0TU;
   Int_t           t0P2;
   Int_t           t0P3;
   Int_t           t0TD;
   Int_t           triggerNumber;
   Int_t           adc[66][32];
   Int_t           tdc[66][32];
   Int_t           q[66];
   Int_t           tdcNhit[66];
   Int_t           driftTime[66][10];
   Int_t           clockNumberDriftTime[66][10];

   // List of branches
   TBranch        *b_iev;   //!
   TBranch        *b_nSamples;   //!
   TBranch        *b_t0TU;   //!
   TBranch        *b_t0P2;   //!
   TBranch        *b_t0P3;   //!
   TBranch        *b_t0TD;   //!
   TBranch        *b_triggerNumber;   //!
   TBranch        *b_adc;   //!
   TBranch        *b_tdc;   //!
   TBranch        *b_q;   //!
   TBranch        *b_tdcNhit;   //!
   TBranch        *b_driftTime;   //!
   TBranch        *b_clockNumberDriftTime;   //!

   hoge(TTree *tree=0);
   virtual ~hoge();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef hoge_cxx
hoge::hoge(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../root/run_000079.event.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../root/run_000079.event.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

hoge::~hoge()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t hoge::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t hoge::LoadTree(Long64_t entry)
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

void hoge::Init(TTree *tree)
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

   fChain->SetBranchAddress("iev", &iev, &b_iev);
   fChain->SetBranchAddress("nSamples", &nSamples, &b_nSamples);
   fChain->SetBranchAddress("t0TU", &t0TU, &b_t0TU);
   fChain->SetBranchAddress("t0P2", &t0P2, &b_t0P2);
   fChain->SetBranchAddress("t0P3", &t0P3, &b_t0P3);
   fChain->SetBranchAddress("t0TD", &t0TD, &b_t0TD);
   fChain->SetBranchAddress("triggerNumber", &triggerNumber, &b_triggerNumber);
   fChain->SetBranchAddress("adc", adc, &b_adc);
   fChain->SetBranchAddress("tdc", tdc, &b_tdc);
   fChain->SetBranchAddress("q", q, &b_q);
   fChain->SetBranchAddress("tdcNhit", tdcNhit, &b_tdcNhit);
   fChain->SetBranchAddress("driftTime", driftTime, &b_driftTime);
   fChain->SetBranchAddress("clockNumberDriftTime", clockNumberDriftTime, &b_clockNumberDriftTime);
   Notify();
}

Bool_t hoge::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void hoge::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t hoge::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef hoge_cxx
