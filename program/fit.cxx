#include <iostream>
#include <iomanip>  
#include <stdlib.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TEllipse.h"
#include "TChain.h"
#include "TLine.h"
#include "TSystem.h"
#include "TCanvas.h"

#define MAX_CHANEL 66 // belle2
#define MAX_CLOCK 32 // belle2
//#define debug

//Origin is at left bottom of 1X tracker                                                                      
double U = 8 ; 
double SEP = 18+600+45*2;
double PI = TMath::Pi();
double T_DELAY = 750+100;
double DRIFT_VEL= 0.03;
double MAX_DRIFTD = U;
int VERBOSE = 0;
int PRINTMODULO = 1000;
long int NEVENTS = 0;
std::vector<std::string> INPUTFILES;

int WID_TX[6]={0,1,2,63,64,65};
int WID_TY[6]={3,4,5,60,61,62};
double POSZ_TX[6]={0,
	               2*U,
	               2*U,
	               SEP+4*U,
	               SEP+6*U,
	               SEP+6*U
	              };
double POSY_TX[6]={0,
	               -U,
	               U,
	               0,
	               -U,
	               U
	              };
double POSZ_TY[6]={4*U,
	               6*U,
	               6*U,
	               SEP,
	               SEP+2*U,
	               SEP+2*U
	              };
double POSX_TY[6]={0,
	               -U,
	               U,
	               0,
	               -U,
	               U
	              };

void print_usage(char* prog_name);
double get_driftD(int*,int);
int get_diffR(const char*,double,double,double,int,int,int,double&,double&,double&);
double GetDist(double ,double , double , double);

//read file
int main(int argc, char* argv[])
{

	//=======================================
	//*************read parameter**********
	int result;
	while((result=getopt(argc,argv,"hd:m:s:v:V:n:p:"))!=-1){
		switch(result){
			/* INPUTS */
			case 'd':
				T_DELAY = atof(optarg);
				printf("Delay time: %lf ns\n",T_DELAY);
				break;
			case 'v':
				DRIFT_VEL = atof(optarg);
				printf("Drift velocity: %lf mm/ns\n",DRIFT_VEL);
				break;
			case 'm':
				MAX_DRIFTD = atof(optarg);
				printf("Maximum drift distance: %lf mm\n",MAX_DRIFTD);
				break;
			case 's':
				SEP = atof(optarg);
				printf("Separation: %lf mm\n",SEP);
				break;
			case 'V':
				VERBOSE = atoi(optarg);
				printf("verbose level: %d\n",VERBOSE);
				break;
			case 'n':
				NEVENTS = atoi(optarg);
				printf("Maximum Number of Events: %ld\n",NEVENTS);
				break;
			case 'p':
				PRINTMODULO = atoi(optarg);
				printf("Print Modulo: %d\n",PRINTMODULO);
				break;
			case '?':
				printf("Wrong option! optopt=%c, optarg=%s\n", optopt, optarg);
				break;
			case 'h':
			default:
				print_usage(argv[0]);
				return 1;
		}
	}

	for (;optind<argc;optind++){
		INPUTFILES.push_back(argv[optind]);
	}
	if (INPUTFILES.size()==0){
		print_usage(argv[0]);
		return 1;
	}

	// Get the files
	TChain * c = new TChain("tree","tree");
	for (int i = 0; i<INPUTFILES.size(); i++){
		c->Add(INPUTFILES[i].c_str());
	}

	int driftTime[MAX_CHANEL][MAX_CLOCK];
	int clockNumberDriftTime[MAX_CHANEL][MAX_CLOCK];
	int adc[MAX_CHANEL][MAX_CLOCK];
	int tdcNhit[MAX_CHANEL];

//	c->SetBranchAddress("adc",&adc);
	c->SetBranchAddress("driftTime",&driftTime);
//	c->SetBranchAddress("clockNumberDriftTime",&clockNumberDriftTime);
	c->SetBranchAddress("tdcNhit",&tdcNhit);

	// Prepare hists and canvas
#ifdef debug
	TCanvas * canvas = new TCanvas("c1","c1",1440,400);
	TH2D * hdis = new TH2D("hdis","Event Display",1024,-10,760,1024,-20,20);
	TEllipse* e1 = 0;
	TEllipse* e2 = 0;
	TEllipse* e3 = 0;
	TLine * l1 = 0;
	TEllipse* e1_best = 0;
	TEllipse* e2_best = 0;
	TEllipse* e3_best = 0;
	TLine * l1_best = 0;
#endif

	TH1D* hres = new TH1D("hres","Residual [mm]",256,-5,5);
	TH1D* hslope = new TH1D("hslope","Slope ",256,-0.01,0.01);
	TH1D* hintersect= new TH1D("hintersect","Intersect [mm]",256,-10,10);

	// Loop in events
	if (!NEVENTS) NEVENTS = c->GetEntries();  
	for (long int iev=0; iev<NEVENTS; iev++) {
		c->GetEntry(iev);

		if (VERBOSE>0&&iev%PRINTMODULO==0) std::cout<<"================================"<<iev<<"=========================="<<std::endl;

		double difRx;
		double R_top;
		double R_bottom;
		double R_choose;
		double slopex;
		double intersectx;

		double difRx_best = 1e9;
		double R_top_best;
		double R_bottom_best;
		double R_choose_best;
		double slopex_best;
		double intersectx_best;
		int it_best,ib_best,ic_best;

#ifdef debug 
		hdis->Draw();
#endif

		//check X hit s
//		for(int it = 0; it<3; it++){
		for(int it = 0; it<=0; it++){
			R_top = get_driftD(driftTime[WID_TX[it]],tdcNhit[WID_TX[it]]);
			if (R_top<0||R_top>MAX_DRIFTD) continue;
//			for(int ib = 3; ib<6; ib++){
			for(int ib = 3; ib<=3; ib++){
				R_bottom = get_driftD(driftTime[WID_TX[ib]],tdcNhit[WID_TX[ib]]);
				if (R_bottom<0||R_bottom>MAX_DRIFTD) continue;
				for(int ic = 0; ic<6; ic++){
					if (ic==it||ic==ib) continue;
					R_choose = get_driftD(driftTime[WID_TX[ic]],tdcNhit[WID_TX[ic]]);
					if (R_choose<0||R_choose>MAX_DRIFTD) continue;
					if (VERBOSE>0&&iev%PRINTMODULO==0)
						std::cout<<"Get ("<<WID_TX[it]<<","<<WID_TX[ib]<<","<<WID_TX[ic]<<"): "<<R_top<<", "<<R_bottom<<", "<<R_choose<<std::endl;
					int result = get_diffR("X",R_top,R_bottom,R_choose,it,ib,ic,difRx,slopex,intersectx);
					if (result) continue;
#ifdef debug 
					canvas->cd();
					if (VERBOSE>0&&iev%PRINTMODULO==0)
						std::cout<<"y="<<slopex<<"*x+"<<intersectx<<": ddhit = "<<R_choose<<", ddfit = "<<R_choose+difRx<<", dddif = "<<difRx<<std::endl;
					l1 = new TLine((-20-intersectx)/slopex,-20, (20-intersectx)/slopex, 20);
					e1 = new TEllipse(POSZ_TX[it],POSY_TX[it],R_top,R_top);
					e2 = new TEllipse(POSZ_TX[ib],POSY_TX[ib],R_bottom,R_bottom);
					e3 = new TEllipse(POSZ_TX[ic],POSY_TX[ic],R_choose,R_choose);
#endif
					if (fabs(difRx_best)>fabs(difRx)){
						difRx_best = difRx;
						slopex_best = slopex;
						intersectx_best = intersectx;
						R_top_best = R_top;
						R_bottom_best =	R_bottom;
						R_choose_best = R_choose;
						it_best = it;
						ib_best = ib;
						ic_best = ic;
#ifdef debug 
						l1_best = l1;
						e1_best = e1;
						e2_best = e2;
						e3_best = e3;
						if (VERBOSE>0&&iev%PRINTMODULO==0)
							std::cout<<"Better!"<<std::endl;
#endif
					}
#ifdef debug 
					e1->SetFillStyle(0);
					e2->SetFillStyle(0);
					e3->SetFillStyle(0);
					e1 -> Draw("SAME");
					e2 -> Draw("SAME");
					e3 -> Draw("SAME");
					l1 -> Draw("SAME");
#endif     
				}
			}
		}
      if (difRx_best!=1e9){
         hres->Fill(difRx_best);
         hslope->Fill(slopex_best);
         hintersect->Fill(intersectx_best);
#ifdef debug 
         l1_best->SetLineColor(kRed);
         e1_best->SetLineColor(kRed);
         e2_best->SetLineColor(kRed);
         e3_best->SetLineColor(kMagenta);
         canvas->Update();
         canvas->SaveAs(Form("%d.png",(int)iev));
#endif     
      }
	}
	TCanvas * canvas2 = new TCanvas();
	hres->Draw();
	canvas2->SaveAs("resolution.png");
	TCanvas * canvas3 = new TCanvas();
	hslope->Draw();
	canvas3->SaveAs("slope.png");
	TCanvas * canvas4 = new TCanvas();
	hintersect->Draw();
	canvas4->SaveAs("intersect.png");
	return 0;
}

// get TDC 
double get_driftD(int* tdc,int n)
{
	double driftD;
	for (int clk=0;clk <n; clk ++ ) {
		driftD = (*(tdc+clk)+T_DELAY)*DRIFT_VEL;
		if(driftD>0&&driftD<MAX_DRIFTD) return driftD; 
	}
	return -1;

}

int get_diffR(const char * opt,double hitRt, double hitRb, double hitRc, int hitIt,int hitIb,int hitIc,double& difR, double&slope, double& intersect)
{

	if (opt!="X") return -1;
	double zc0 = POSZ_TX[hitIc];
	double yc0 = POSY_TX[hitIc];
	double zt0 = POSZ_TX[hitIt];
	double yt0 = POSY_TX[hitIt];
	double zb0 = POSZ_TX[hitIb];
	double yb0 = POSY_TX[hitIb];

	double theta = atan((yb0-yt0)/(zb0-zt0)); // the angle from x axis to top->bottom line
	double D = sqrt(pow(yb0-yt0,2)+pow(zb0-zt0,2)); // the distance between top and bottom wire.

	double zc = (zc0-zt0)*cos(theta)+(yc0-yt0)*sin(theta);
	double yc = -(zc0-zt0)*sin(theta)+(yc0-yt0)*cos(theta);

	double alphaE = asin((hitRt-hitRb)/D);
	double alphaI = asin((hitRt+hitRb)/D);

	double dif[4] = {0};
	dif[0] = GetDist(zc,yc,-tan(alphaE),hitRt/cos(alphaE))-hitRc;
	dif[1] = GetDist(zc,yc,tan(alphaE),-hitRt/cos(alphaE))-hitRc;
	dif[2] = GetDist(zc,yc,-tan(alphaI),hitRt/cos(alphaI))-hitRc;
	dif[3] = GetDist(zc,yc,tan(alphaI),-hitRt/cos(alphaI))-hitRc;
//   std::cout<<"  =>EU: "<<-tan(alphaE)<<" "<<hitRt/cos(alphaE)<<" "<<dif[0]<<std::endl;
//   std::cout<<"  =>ED: "<<tan(alphaE)<<" "<<-hitRt/cos(alphaE)<<" "<<dif[1]<<std::endl;
//   std::cout<<"  =>IU: "<<-tan(alphaI)<<" "<<hitRt/cos(alphaI)<<" "<<dif[2]<<std::endl;
//   std::cout<<"  =>ID: "<<tan(alphaI)<<" "<<-hitRt/cos(alphaI)<<" "<<dif[3]<<std::endl;

	difR = 1e9;
	int index = 0;
	for(int i = 0; i<4; i++){
		if (fabs(difR)>fabs(dif[i])){
			index = i;
			difR = dif[i];
		}
	}

//	std::cout<<"Rhit = "<<hitRc<<", zc = "<<zc<<", yc = "<<yc<<std::endl;
	if (index == 0){
		slope = -tan(alphaE);
		intersect = hitRt/cos(alphaE);
//		std::cout<<"External Up: slope="<<slope<<", intersect = "<<intersect<<", Rfit = "<<dif[0]+hitRc<<", Rdif = "<<dif[0]<<std::endl;
	}
	else if (index == 1){
		slope = tan(alphaE);
		intersect = -hitRt/cos(alphaE);
//		std::cout<<"External Down: slope="<<slope<<", intersect = "<<intersect<<", Rfit = "<<dif[1]+hitRc<<", Rdif = "<<dif[1]<<std::endl;
	}
	else if (index == 2){
		slope = -tan(alphaI);
		intersect = hitRt/cos(alphaI);
//		std::cout<<"Internal Up: slope="<<slope<<", intersect = "<<intersect<<", Rfit = "<<dif[2]+hitRc<<", Rdif = "<<dif[2]<<std::endl;
	}
	else{
		slope = tan(alphaI);
		intersect = -hitRt/cos(alphaI);
//		std::cout<<"Internal Down: slope="<<slope<<", intersect = "<<intersect<<", Rfit = "<<dif[3]+hitRc<<", Rdif = "<<dif[3]<<std::endl;
	}

	double slope_new = (slope*cos(theta)+sin(theta))/(cos(theta)-slope*sin(theta));
	double intersect_new = slope_new*(intersect*sin(theta)-zt0)+intersect*cos(theta)+yt0;
	slope = slope_new;
	intersect = intersect_new;

	return 0;
}

double GetDist(double x ,double y, double a, double b)
	//Dist from (x,y) to y=ax+b
{
	return fabs((a*x+b-y)/(1+a*a));
}

void print_usage(char* prog_name)
{
	fprintf(stderr,"Usage %s [options (args)] [input files]\n",prog_name);
	fprintf(stderr,"[options]\n");
	fprintf(stderr,"\t -d\n");
	fprintf(stderr,"\t\t set delay time [ns]. Default %lf\n",T_DELAY);
	fprintf(stderr,"\t -v\n");
	fprintf(stderr,"\t\t set velocity [mm/ns]. Default %lf\n",DRIFT_VEL);
	fprintf(stderr,"\t -s\n");
	fprintf(stderr,"\t\t set separation [mm]. Default %lf\n",SEP);
	fprintf(stderr,"\t -V\n");
	fprintf(stderr,"\t\t set verbose level. Default %d\n",VERBOSE);
	fprintf(stderr,"\t -m\n");
	fprintf(stderr,"\t\t set maximum drift distance [mm]. Default %lf\n",MAX_DRIFTD);
	fprintf(stderr,"\t -n\n");
	fprintf(stderr,"\t\t Maximum number of events\n");
	fprintf(stderr,"\t -p\n");
	fprintf(stderr,"\t\t print Modulo. Default %d\n",PRINTMODULO);
	fprintf(stderr,"\t -h\n");
	fprintf(stderr,"\t\t Usage message.\n");
	fprintf(stderr,"[example]\n");
	fprintf(stderr,"\t\t%s -d 845 -v 0.0305 -m 8 -n 10 -V 20 -p 1\n",prog_name);
}
