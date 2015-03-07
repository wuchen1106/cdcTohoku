#include "TH1D.h"
#include "TChain.h"
#include <stdlib.h>
#include <iostream>

int main(int argc, char** argv){
	if (argc<2){
		return 1;
	}
	int startNo = (int)strtol(argv[1],NULL,10);

	TChain * ch = new TChain("t","t");
	ch->Add(Form("../root/d_%d.root",startNo));
	int dt[66][32];
	int n[66];
	ch->SetBranchAddress("dt",dt);
	ch->SetBranchAddress("n",n);
	TH1D * h[66];
	for ( int i = 0; i<66; i++){
		h[i] = new TH1D(Form("h%d",i),"h",60,-860,-800);
	}

	Long64_t N = ch->GetEntries();
	for(Long64_t iEvent = 0; iEvent<N; iEvent++){
		if (iEvent%1000==0) std::cout<<(double)iEvent/N*100<<"..."<<std::endl;
		ch->GetEntry(iEvent);
		for ( int i = 0; i<66; i++){
			// FIXME
			//if (n[i]>=0) h[i]->Fill(dt[i][0]);
			if (n[i]>=0&&dt[i][0]!=-1) h[i]->Fill(dt[i][0]);
		}
	}

	for(int i = 0; i <4; i++){
		std::vector<int> indice;
		if (i == 0){for(int j = 0; j<6; j++) indice.push_back(j);}
		if (i == 1){indice.push_back(15);indice.push_back(19);indice.push_back(24);}
		if (i == 2){indice.push_back(42);indice.push_back(46);indice.push_back(51);}
		if (i == 3){for(int j = 60; j<66; j++) indice.push_back(j);}
		double t0 = 0;
		for ( int j = 0; j<indice.size(); j++){
			int ch = indice[j];
			int maxbin = h[ch]->GetMaximumBin();
			int max = h[ch]->GetBinContent(maxbin);
			int min = 0;
			for (int k = 1; k<=5; k++){
				min += h[ch]->GetBinContent(k);
			}
			min /= 5.;
			double th = (max-min)/20.*15+min;
			int ibin = 1;
			for (; ibin<=60; ibin++){
				int height = h[ch]->GetBinContent(ibin);
				if (height>th) break;
			}
			t0 += h[ch]->GetBinCenter(ibin);
			printf("#\t%d %d %lf\n",ch,ibin,h[ch]->GetBinCenter(ibin));
		}
		t0/=(indice.size());
		printf("%d %lf\n",i,t0);
	}
}
