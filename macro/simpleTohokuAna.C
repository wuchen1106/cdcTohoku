#define simpleTohokuAna_cxx
#include "simpleTohokuAna.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void simpleTohokuAna::Loop()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //--------------------------------------------------
  // Prepare canvas, histograms,,,
  //--------------------------------------------------


  //--------------------------------------------------
  //   Event loop 
  //--------------------------------------------------
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // - - - - - - - - - - - - - - - - 
    // Add your event process here
    // - - - - - - - - - - - - - - - -

      
  } // End of the event loop
}

#define NUM_OF_WAVE_SAMPLES 32
#define NUM_OF_WIRES 66
#define MIN_WAVE_CH 190
#define MAX_WAVE_CH 400

//=========================================================================
void simpleTohokuAna::plotWaveForm(Int_t wireID, Int_t firstEventNumber)
//=========================================================================
// Plot waveform of a wire with a charge q and tdcHit informations
//  07-Dec-2014 A.Sato
//=========================================================================
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //--------------------------------------------------
  // Prepare canvas, histograms,,,
  //--------------------------------------------------
  TCanvas *c1 = new TCanvas("c1","Event Display",800,600);
  TGraph *gr_waveForm[NUM_OF_WIRES];
  Int_t vSample[NUM_OF_WAVE_SAMPLES];
  for (Int_t i=0; i<NUM_OF_WAVE_SAMPLES; i++){
    vSample[i] = i;
  }

  //--------------------------------------------------
  //   Event loop 
  //--------------------------------------------------
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=firstEventNumber; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // - - - - - - - - - - - - - - - - 
    // Add your event process here
    // - - - - - - - - - - - - - - - -
    gr_waveForm[wireID] = new TGraph(NUM_OF_WAVE_SAMPLES,vSample,adc[wireID]);
    gr_waveForm[wireID]->GetHistogram()->SetTitle(Form("Waveform: wireID=%d, event=%d",wireID,jentry));
    gr_waveForm[wireID]->GetHistogram()->SetMinimum(MIN_WAVE_CH); 
    gr_waveForm[wireID]->GetHistogram()->SetMaximum(MAX_WAVE_CH);
    gr_waveForm[wireID]->GetHistogram()->SetXTitle("Samplig ID");
    gr_waveForm[wireID]->GetHistogram()->SetYTitle("ADC(ch)");
    gr_waveForm[wireID]->SetMarkerStyle(24);
    gr_waveForm[wireID]->SetMarkerSize(0.8);
    gr_waveForm[wireID]->Draw("apwl");
    TLatex *text1 = new TLatex(NUM_OF_WAVE_SAMPLES/18, 
			       MAX_WAVE_CH-1*(MAX_WAVE_CH-MIN_WAVE_CH)/10,
			       Form("q=%d",q[wireID]));
    text1->Draw();
    TLatex *textTDC[32];
    TMarker *markerTDC[32];
    
    for (Int_t i=0; i<tdcNhit[wireID]; i++) {
      textTDC[i] = new TLatex(clockNumberDriftTime[wireID][i], 
			      MIN_WAVE_CH+10,
			      Form("%d",driftTime[wireID][i]));
      textTDC[i]->SetTextSize(0.02);
      textTDC[i]->Draw();
      markerTDC[i] = new TMarker(clockNumberDriftTime[wireID][i],
				 adc[wireID][clockNumberDriftTime[wireID][i]],
				 20);
      markerTDC[i]->SetMarkerSize(0.7);
      markerTDC[i]->SetMarkerColor(2);
      markerTDC[i]->Draw();
    }
    c1->WaitPrimitive();
  
  } // End of the event loop
}
