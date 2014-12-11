#include <iostream>
#include <stdlib.h>
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"

void print_usage(char* prog_name);
int main(int argc, char** argv){
	if (argc!=3){
		print_usage(argv[0]);
		return 1;
	}
	int startRunNo = (int)strtol(argv[1],NULL,10);
	int endRunNo = (int)strtol(argv[2],NULL,10);

	TFile * output = new TFile("../root/eff.root","UPDATE");
	for (int runNo = startRunNo; runNo<=endRunNo; runNo++){
		TChain * c = new TChain("tree","tree");
		c->Add(Form("../root/run_%06d_built.root",runNo));
		int tdcNhit[66];
		int ch2x[66] = {0,1,1,
			2,3,3,
			8,7,6,5,4,
			9,8,7,6,5,4,
			8,7,6,5,4,
			9,8,7,6,5,4,
			8,7,6,5,4,
			14,13,12,11,10,
			15,14,13,12,11,10,
			14,13,12,11,10,
			15,14,13,12,11,10,
			14,13,12,11,10,
			16,17,17,
			18,19,19
		};
		int ch2y[66] = {2,1,3,
			2,1,3,
			4,4,4,4,4,
			3,3,3,3,3,3,
			2,2,2,2,2,
			1,1,1,1,1,1,
			0,0,0,0,0,
			4,4,4,4,4,
			3,3,3,3,3,3,
			2,2,2,2,2,
			1,1,1,1,1,1,
			0,0,0,0,0,
			2,1,3,
			2,1,3
		};
		c->SetBranchAddress("tdcNhit",tdcNhit);
		TH2D * h1 = new TH2D(Form("h1_%d",runNo),Form("Run #%d Events With Hit Over All Triggered Events",runNo),20,0,20,5,0,5);
		TH2D * h2 = new TH2D(Form("h2_%d",runNo),Form("Run #%d Events With Hit Over Center Beam Events",runNo),20,0,20,5,0,5);
		Long64_t nEntries = c->GetEntries();
		Long64_t nGoodEvents = 0;
		bool good = false;
		for (Long64_t i = 0;i<nEntries; i++){
			c->GetEntry(i);
			good = false;
			if (tdcNhit[0]&&tdcNhit[3]&&tdcNhit[60]&&tdcNhit[63]){
				good = true;
				nGoodEvents++;
			}
			for(int j = 0; j<66; j++){
				if (good){
					h2->Fill(ch2x[j],ch2y[j],(tdcNhit[j]>0));
				}
				h1->Fill(ch2x[j],ch2y[j],(tdcNhit[j]>0));
			}
		}
		h1->Scale(1./nEntries);
		h2->Scale(1./nGoodEvents);
		gStyle->SetPaintTextFormat("4.3f");
		TCanvas * c1 = new TCanvas("c1","c1",1200,600);
		h1->Draw("COLZTEXT");
		c1->SaveAs(Form("../pdf/%d_eff1.pdf",runNo));
		c1->SaveAs(Form("../pdf/%d_eff1.png",runNo));
		TCanvas * c2 = new TCanvas("c2","c2",1200,600);
		h2->Draw("COLZTEXT");
		c2->SaveAs(Form("../pdf/%d_eff2.pdf",runNo));
		c2->SaveAs(Form("../pdf/%d_eff2.png",runNo));
		printf("Total Events: %d\n",nEntries);
		printf("Total Good Events: %d\n",nGoodEvents);
		printf("Tracking Efficiency: %lf\n",(double)nGoodEvents/nEntries);
		h1->Write();
		h2->Write();
	}
	output->Close();
}

void print_usage(char* prog_name)
{
	fprintf(stderr,"\t%s [startRunNo] [endRunNo]\n",prog_name);
}
