#define ePlotEvent_cxx
#include "ePlotEvent.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#define NUM_OF_WAVE_SAMPLES 32
#define NUM_OF_WIRES 66
//----------------------------------------------------
// Wire definition : wireNo.(wireNo.if the chamber)
//----------------------------------------------------
// UX 00(1)  01(2)  02(3)
// UY 03(1)  04(2)  05(3)
//
// P3        06(0)  07(1)  08(2)  09(3)  10(4)
// P3 11(5)  12(6)  13(7)  14(8)  15(9)  16(10)
// P3        17(11) 18(12) 19(13) 20(14) 21(15) 
// P3 22(16) 23(17) 24(18) 25(19) 26(20) 27(21)
// P3        28(22) 29(23) 30(24) 31(25) 32(26) 
//
// P2        33(0)  34(1)  35(2)  36(3)  37(4)
// P2 38(5)  39(6)  40(7)  41(8)  42(9)  43(10)
// P2        44(11) 45(12) 46(13) 47(14) 48(15) 
// P2 49(16) 50(17) 51(18) 52(19) 53(20) 54(21)
// P2        55(22) 56(23) 57(24) 58(25) 59(26) 
//
// DY 60(1) 61(2) 62(3)
// DX 63(1) 64(2) 65(3)
//----------------------------------------------------

Int_t padPos[NUM_OF_WIRES] = {
  11, 21, 1,
  12, 22, 2,
      -1, -1, -1, -1, -1,
  -1, -1,  3, 13, 23, -1,
      -1,  4, 14, 24, -1,
  -1, -1,  5, 15, 25, -1,
      -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1,
  -1, -1,  6, 16, 26, -1,
      -1,  7, 17, 27, -1,
  -1, -1,  8, 18, 28, -1,
      -1, -1, -1, -1, -1,
  19, 29, 9,
  20, 30, 10
};

#define MIN_WAVE_CH 200
#define MAX_WAVE_CH 400


Int_t vSample[NUM_OF_WAVE_SAMPLES];

void ePlotEvent::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   //--------------------------------------------------
   // Prepare canvas, histograms,,,
   //--------------------------------------------------
   TCanvas *c1 = new TCanvas("c1","Event Display",1800,800);
   c1->Divide(10,3);
   TH2F *h[NUM_OF_WIRES];
   TGraph *grWave[NUM_OF_WIRES];
   for (Int_t iWire=0;iWire<NUM_OF_WIRES;iWire++) {
     h[iWire] = NULL;
     grWave[iWire]= NULL;
   }
   for (Int_t i=0; i<NUM_OF_WAVE_SAMPLES; i++){
     vSample[i] = i;
   }
   
   //--------------------------------------------------
   //   Event loop 
   //--------------------------------------------------
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // - - - - - - - - - - - - - - - - 
      // Add your event process here
      // - - - - - - - - - - - - - - - -
      for (Int_t iWire=0;iWire<NUM_OF_WIRES;iWire++) {
	if (padPos[iWire]==-1) continue;
	c1->cd(padPos[iWire]);
	grWave[iWire] = new TGraph(NUM_OF_WAVE_SAMPLES,vSample,adc[iWire]);
	grWave[iWire]->GetHistogram()->SetMinimum(MIN_WAVE_CH); 
	grWave[iWire]->GetHistogram()->SetMaximum(MAX_WAVE_CH); 
	grWave[iWire]->Draw("alw");
      }
      printf("click pad\n");
      c1->WaitPrimitive();



   } // End of the event loop
}
