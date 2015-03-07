#include <iostream>
#include <stdlib.h>
#include <string>
#include <math.h>
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

double TIMESHIFT = 0;

void print_usage(char* prog_name);
int main(int argc, char** argv){
	if (argc<2){
		print_usage(argv[0]);
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);
	int nEventMax = 0;
	if (argc>=3) nEventMax = (int)strtol(argv[2],NULL,10);
	std::string suffix = "";
	if (argc>=4){
		suffix  = argv[3];
		suffix=suffix+".";
	}
	if (argc>=5){
		TIMESHIFT = strtof(argv[4],NULL);
	}

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
	c->SetBranchAddress("driftTime",driftTime);
	c->SetBranchAddress("adc",adc);

	char outputName[128];
	int o_iHits[66];
	int o_nHits[66];
	double o_pedestal[66];
	int o_pedestalN[66];
	double o_pedestalChi2[66];
	double o_areaall[66];
	double o_area[66][32];
	int o_time[66][32];
	int o_height[66][32];
	//	sprintf(outputName,"d_%d.fp.garv3t1p1.root",runNo);
	sprintf(outputName,("../root/d_%d."+suffix+"root").c_str(),runNo);
	TFile * f = new TFile(outputName,"RECREATE");
	TTree * t = new TTree("t","t");
	t->Branch("i",o_iHits,"i[66]/I");
	t->Branch("n",o_nHits,"n[66]/I");
	t->Branch("p",o_pedestal,"p[66]/D");
	t->Branch("pn",o_pedestalN,"pn[66]/I");
	t->Branch("pchi2",o_pedestalChi2,"pchi2[66]/D");
	t->Branch("aa",o_areaall,"aa[66]/D");
	t->Branch("a",o_area,"a[66][32]/D");
	t->Branch("dt",o_time,"dt[66][32]/I");
	t->Branch("h",o_height,"h[66][32]/I");
	//FIXME
	t->Branch("clockNumberDriftTime",clockNumberDriftTime,"clockNumberDriftTime[66][32]/I");
	t->Branch("adc",adc,"adc[66][32]/I");

	Long64_t N = c->GetEntries();
	if (nEventMax&&nEventMax<N) N = nEventMax;
	//N=1;
	int tdcNhitwire = -1;
	for (Long64_t i = 0;i<N; i++){
		//FIXME
		if (i%1000==0) std::cout<<(double)i/N*100<<"%..."<<std::endl;
		c->GetEntry(i);
		int cwire = 0;
		for(int iwire = 0; iwire<66; iwire++){
			tdcNhitwire = tdcNhit[iwire];
			o_nHits[cwire] = tdcNhitwire;
			o_iHits[cwire]=-1;
			for(int ihit = 0; ihit<32; ihit++){
				o_height[cwire][ihit] = -1;
				o_time[cwire][ihit] =-1;
			}

			// FIXME
			//for(int ihit=0; ihit<tdcNhitwire; ihit++){
			//   //printf("for hit[%d], start from %d and stop at%d\n",ihit,clockNumberDriftTime[iwire][ihit],ihit+1>=tdcNhitwire?32:clockNumberDriftTime[iwire][ihit+1]);
			//   for(int clk = clockNumberDriftTime[iwire][ihit]; clk<(ihit+1>=tdcNhitwire?32:clockNumberDriftTime[iwire][ihit+1]); clk++){
			//	   if (adc[iwire][clk]>o_height[cwire]){
			//	   	   o_height[cwire]=adc[iwire][clk];
			//	   	   o_time[cwire]=driftTime[iwire][ihit]-TIMESHIFT;
			//	   	   o_iHits[cwire]=ihit;
			//	   	   if( t2x(iwire,o_time[cwire],driftD) ) driftD = -1;
			//	   	   o_dist[cwire] = driftD;
			//		   //printf("adc[%d][%d]=%d>o_height[%d]=%d\n",iwire,clk,adc[iwire][clk],ihit,o_height[ihit]);
			//		   //printf("time[%d]  = driftTime[%d][%d] - TIMESHIFT = %d - %lf = %d\n",cwire,iwire,ihit,driftTime[iwire][ihit],TIMESHIFT,o_time[cwire]);
			//	   }
			//   }
			//}
			o_iHits[cwire]=0;
			o_pedestal[cwire]=0;
			int clk;
			for(clk = 0; clk<(tdcNhitwire<=0?32:clockNumberDriftTime[iwire][0]-1); clk++){
				o_pedestal[cwire]+=adc[iwire][clk];
			}
			o_pedestal[cwire]/=clk;
			o_pedestalN[cwire] = clk;
			o_pedestalChi2[cwire]=0;
			for(clk = 0; clk<(tdcNhitwire<=0?32:clockNumberDriftTime[iwire][0]); clk++){
				o_pedestalChi2[cwire]+=pow(adc[iwire][clk]-o_pedestal[cwire],2);
			}
			o_areaall[cwire] = 0;
			for ( int ihit = 0; ihit<tdcNhitwire; ihit++){
				// FIXME
				//if (tdcNhitwire>0) o_time[cwire][ihit]=driftTime[iwire][ihit]-TIMESHIFT;
				o_time[cwire][ihit]=driftTime[iwire][ihit];
				for(int clk = clockNumberDriftTime[iwire][ihit]; clk<(ihit+1>=tdcNhitwire?32:clockNumberDriftTime[iwire][ihit+1]); clk++){
					if (adc[iwire][clk]<o_height[cwire][ihit]){
						break;
					}
					o_height[cwire][ihit]=adc[iwire][clk];
				}

				o_area[cwire][ihit] = 0;
				for(int clk = clockNumberDriftTime[iwire][ihit]; clk<(ihit+1>=tdcNhitwire?32:clockNumberDriftTime[iwire][ihit+1]); clk++){
					if (adc[iwire][clk]<o_pedestal[cwire]) break;
					o_area[cwire][ihit] += adc[iwire][clk]-o_pedestal[cwire];
				}
				o_areaall[cwire] += o_area[cwire][ihit];
			}

			cwire++;
		}
		t->Fill();
	}
	t->Write();
	return 0;
}

void print_usage(char* prog_name)
{
	fprintf(stderr,"\t%s [runNo] <[nEventMax]>\n",prog_name);
}
