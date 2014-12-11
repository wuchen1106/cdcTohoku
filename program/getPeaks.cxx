#include <iostream>
#include <stdlib.h>
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

void print_usage(char* prog_name);
int main(int argc, char** argv){
	if (argc!=2){
		print_usage(argv[0]);
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);

	TChain * c = new TChain("tree","tree");
	char inputName[128];
	sprintf(inputName,"../root/run_000%d_built.root",runNo);
	c->Add(inputName);
	int tdcNhit[66];
	int clockNumberDriftTime[66][32];
	int adc[66][32];
	int driftTime[66][32];
	c->SetBranchAddress("tdcNhit",tdcNhit);
	c->SetBranchAddress("clockNumberDriftTime",clockNumberDriftTime);
	c->SetBranchAddress("adc",adc);

	int height[32];
	int iwire;
	int adcwire[32];
	int tdcNhitwire;
	char outputName[128];
	sprintf(outputName,"peak_%d.root",runNo);
	TFile * f = new TFile(outputName,"RECREATE");
	TTree * t = new TTree("t","t");
	t->Branch("h",height,"h[32]/I");
	t->Branch("tdcNhit",&tdcNhitwire);
	t->Branch("wire",&iwire);
	t->Branch("adcwire",adcwire,"adcwire[32]/I");

	Long64_t N = c->GetEntries();
	//N=1;
	for (Long64_t i = 0;i<N; i++){
		if (i%1000==0) printf("%lf%\n",(double)i/N*100);
		c->GetEntry(i);
		for(iwire = 0; iwire<66; iwire++){
			tdcNhitwire = tdcNhit[iwire];
			//printf("wire = %d, Nhit = %d\n",iwire,tdcNhitwire);
			if (tdcNhitwire<=0) continue;
			for(int ihit=0; ihit<tdcNhit[iwire]; ihit++){
			   height[ihit] = -1;
			   //printf("for hit[%d], start from %d and stop at%d\n",ihit,clockNumberDriftTime[iwire][ihit],ihit+1>=tdcNhitwire?31:clockNumberDriftTime[iwire][ihit+1]);
			   for(int clk = clockNumberDriftTime[iwire][ihit]; clk<=(ihit+1>=tdcNhitwire?31:clockNumberDriftTime[iwire][ihit+1]); clk++){
				   if (adc[iwire][clk]>height[ihit]){
				   	   height[ihit]=adc[iwire][clk];
					   //printf("adc[%d][%d]=%d>height[%d]=%d\n",iwire,clk,adc[iwire][clk],ihit,height[ihit]);
				   }
			   }
			}
			for (int j = 0; j < 32 ; j++)
				adcwire[j]=adc[iwire][j];
			t->Fill();
		}
	}
	t->Write();
	return 0;
}

void print_usage(char* prog_name)
{
	fprintf(stderr,"\t%s [runNo]\n",prog_name);
}
