#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

int main(int argc, char** argv){
	TChain * c = new TChain("tree","tree");
	c->Add("../root/run_000220_built.root");
	int tdcNhit[66];
	int clockNumberDriftTime[66][32];
	int adc[66][32];
	int driftTime[66][32];
	c->SetBranchAddress("tdcNhit",tdcNhit);
	c->SetBranchAddress("clockNumberDriftTime",clockNumberDriftTime);
	c->SetBranchAddress("adc",adc);

	double height;
	int iwire;
	int adcwire[32];
	TFile * f = new TFile("output.root","RECREATE");
	TTree * t = new TTree("t","t");
	t->Branch("h",height);
	t->Branch("wire",iwire);
	t->Branch("adcwire",adcwire,"adcwire[32]/I");

	for (Long64_t i = 0;i<c->GetEntries(); i++){
		c->GetEntry(i);
		for(iwire = 0; iwire<66; iwire++){
			//		  for(int ihit=0; ihit<tdcNhit[iwire]; ihit++){
			//			 height = -1;
			//			 for(int clk = clockNumberDriftTime[iwire][ihit]; clk<(ihit>=31?32:clockNumberDriftTime[iwire][ihit+1]); clk++){
			//				if (adc[iwire][clk]>height) height=adc[iwire][clk];
			//			 }
			//			 if (height!=-1) t->Fill();
			//		  }
			if (tdcNhit[iwire]<=0) continue;
			height = -1;
			for(int clk = 0; clk<32; clk++){
				if (adc[iwire][clk]>height) height=adc[iwire][clk];
			}
			if (height!=-1){
				adcwire=adc[iwire];
//				std::cout<<"adc["<<iwire<<"] = "<<*(adc[iwire])<<"; adc["<<iwire<<"][0] = "<<adc[iwire][0]<<std::endl;
				t->Fill();
			}
		}
	}
	t->Write();
	return 0;
}
