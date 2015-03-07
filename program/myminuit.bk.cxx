#include <iostream>
#include <iomanip>  
#include <stdlib.h>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TEllipse.h"
#include "TChain.h"
#include "TLine.h"
#include "TText.h"
#include "TStyle.h"
#include "TPad.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include "TEllipse.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TGraph.h"

# define PRINT

#define NUM_OF_WAVE_SAMPLES 32
#define NHITS_TO_CHOOSE 4
#define NHITS_TO_CHOOSET 4
#define MIN_WAVE_CH 190
#define MAX_WAVE_CH 400
#define NPAR  2
#define NBINS  256

// ________About Geometry___________
double U = 8 ; 
double SEP = 18+600+45*2;
double PI = TMath::Pi();
Double_t Amin[12]={
	0.5,
	0.5,
	0.5,
	7,
	7,
	7,
	0.5,
	0.5,
	0.5,
	10,
	7,
	0.5
};

// ________About Hit Pattern___________
Int_t patternX = 0;
Int_t patternY = 0;
Int_t indiceX[NHITS_TO_CHOOSE] = {0};
Int_t indiceY[NHITS_TO_CHOOSE] = {0};
Int_t itestX = 0;
Int_t itestY = 0;

// ________About Signal___________
double TIMESHIFT = -835;
double f1 = 1;
double f2 = 1;

// ________About Track___________
Double_t z[12],y[12],errord[12];
Double_t dd[12][NUM_OF_WAVE_SAMPLES];
int dt[12][NUM_OF_WAVE_SAMPLES];
double a[12][NUM_OF_WAVE_SAMPLES];
int h[12][NUM_OF_WAVE_SAMPLES];
Double_t * par = (Double_t *)malloc(NPAR*sizeof(Double_t));
Int_t iHit[12];

// ________About Fitting___________
Int_t fit_channel = 0; // 0: X Tracker; 1: Y Tracker

//______________________________________________________________________________
Double_t func(float z,float y,Double_t *par)
{
	Double_t A = par[0]*z+par[1]-y;
	Double_t value = fabs(A/sqrt(1+par[0]*par[0]));
	return value;
}

//______________________________________________________________________________
void getchi2(Double_t &f, Double_t *par)
{
	Int_t index;

	//calculate chisquare
	Double_t chisq = 0;
	Double_t delta;
	Double_t dfit;
	for (Int_t i=0;i<NHITS_TO_CHOOSE; i++) {
		if (fit_channel == 0)
			index = indiceX[i];
		else if (fit_channel == 1)
			index = indiceY[i];
		dfit = func(z[index],y[index],par);
		delta  = (dd[index][iHit[index]]-dfit)/errord[index];
		chisq += delta*delta;
		//printf("\t%dd: [%dd] %lf-%lf=%lf\n",i,index,dd[index][iHit[index]],dfit,delta);
	}
	//printf("[%lf,%lf]: chi2 = %lf\n",par[0],par[1],chisq);
	f = chisq;
}

//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	getchi2(f,par);
}

//______________________________________________________________________________
void getHits(){
	// ...
	// FIXME
	bool flag[12];
	for (int i = 0; i<12; i++){
		if (fabs(dd[i][iHit[i]]-U/2)<U/2&&a[i][iHit[i]]>Amin[i]) flag[i] = true;
		else flag[i] = false;
	}
	// X tracker
	if (flag[0]&&flag[9]){
		indiceX[0] = 0;
		indiceX[2] = 9;
		if (flag[1]&&!flag[2]&&flag[10]&&!flag[11]){
			patternX = 0;
			indiceX[1] = 1;
			indiceX[3] = 10;
		}
		else if (!flag[1]&&flag[2]&&!flag[10]&&flag[11]){
			patternX = 1;
			indiceX[1] = 2;
			indiceX[3] = 11;
		}
		else if (!flag[1]&&flag[2]&&flag[10]&&!flag[11]){
			patternX = 2;
			indiceX[1]= 2;
			indiceX[3] = 10;
		}
		else if (flag[1]&&!flag[2]&&!flag[10]&&flag[11]){
			patternX = 3;
			indiceX[1]= 1;
			indiceX[3] = 11;
		}
		else{
			itestX = 0;
			patternX = -1;
		}
	}
	else{
		itestX = 0;
		patternX = -1;
	}
	/*
	if (flag[0]&&flag[9]){
		indiceX[0] = 0;
		indiceX[1] = 9;
		if (flag[1]&&!flag[2]&&flag[10]&&!flag[11]){
			patternX = 0;
			itestX = 1;
			indiceX[2] = 10;
		}
		else if (!flag[1]&&flag[2]&&!flag[10]&&flag[11]){
			patternX = 1;
			itestX = 2;
			indiceX[2] = 11;
		}
		else if (!flag[1]&&flag[2]&&flag[10]&&!flag[11]){
			patternX = 2;
			itestX = 2;
			indiceX[2] = 10;
		}
		else if (flag[1]&&!flag[2]&&!flag[10]&&flag[11]){
			patternX = 3;
			itestX = 1;
			indiceX[2] = 11;
		}
		else{
			itestX = 0;
			patternX = -1;
		}
	}
	else{
		itestX = 0;
		patternX = -1;
	}
	*/
	// Y tracker
	if (flag[3]&&flag[6]){
		indiceY[0] = 3;
		indiceY[1] = 6;
		if (flag[4]&&!flag[5]&&flag[7]&&!flag[8]){
			patternY = 0;
			itestY = 4;
			indiceY[2] = 7;
		}
		else if (!flag[4]&&flag[5]&&!flag[7]&&flag[8]){
			patternY = 1;
			itestY = 5;
			indiceY[2] = 8;
		}
		else if (!flag[4]&&flag[5]&&flag[7]&&!flag[8]){
			patternY = 2;
			itestY = 5;
			indiceY[2] = 7;
		}
		else if (flag[4]&&!flag[5]&&!flag[7]&&flag[8]){
			patternY = 3;
			itestY = 4;
			indiceY[2] = 8;
		}
		else{
			itestY = 0;
			patternY = -1;
		}
	}
	else{
		itestY = 0;
		patternY = -1;
	}
	//if (NHITS_TO_CHOOSE)
	//	printf("%dd, %lf,%lf,%lf,%lf\n",pattern,dd[indice[0]][iHit[indice[0]]],dd[itestX][iHit[itestX]],dd[indice[1]][iHit[indice[1]]],dd[indice[2]][iHit[indice[2]]]);
}

//______________________________________________________________________________
void print_usage(char* prog_name)
{
	fprintf(stderr,"\t%s [startNo] <[nFiles] [suffix] [nEventMax]>\n",prog_name);
}

//______________________________________________________________________________
int main(int argc, char** argv){
	if (argc<2){
		print_usage(argv[0]);
		return 1;
	}
	int startNo = (int)strtol(argv[1],NULL,10);
	int nFiles = 1;
	if (argc>=3) nFiles = (int)strtol(argv[2],NULL,10);
	std::string suffix = "";
	if (argc>=4){
		suffix  = argv[3];
		suffix=suffix+".";
	}
	int nEventMax = 0;
	if (argc>=5) nEventMax = (int)strtol(argv[4],NULL,10);
	if (argc>=6){
		TIMESHIFT = strtof(argv[5],NULL);
	}
	double f1 = 1;
	if (argc>=7){
		f1 = strtof(argv[6],NULL);
	}
	double f2 = 1;
	if (argc>=8){
		f2 = strtof(argv[7],NULL);
	}

	//===================Set Geometry============================
	// The z values
	z[0] =0;
	z[1] =2*U;
	z[2] =2*U;
	z[3] =4*U;
	z[4] =6*U;
	z[5] =6*U;
	z[6] =SEP;
	z[7] =SEP+2*U;
	z[8] =SEP+2*U;
	z[9] =SEP+4*U;
	z[10]=SEP+6*U;
	z[11]=SEP+6*U;
	// the y values
	y[0] =0;
	y[1] =-U;
	y[2] =U;
	y[3] =0;
	y[4] =-U;
	y[5] =U;
	y[6] =0;
	y[7] =-U;
	y[8] =U;
	y[9] =0;
	y[10]=-U;
	y[11]=U;
	// The errors on z values
	Float_t error = 0.2;
	errord[0]=error;
	errord[1]=error;
	errord[2]=error;
	errord[3]=error;
	errord[4]=error;
	errord[5]=error;
	errord[6]=error;
	errord[7]=error;
	errord[8]=error;
	errord[9]=error;
	errord[10]=error;
	errord[11]=error;

	//===================Get ROOT Files============================
	TChain * c = new TChain("t","t");
	std::stringstream buf;
	for ( int iFile = startNo; iFile<startNo+nFiles; iFile++){
		buf.str(""); buf.clear();
		buf<<"../root/d_"<<iFile<<".root";
		c->Add(buf.str().c_str());
		printf("Add(%s)\n",buf.str().c_str());
	}
	int adc[12][NUM_OF_WAVE_SAMPLES];
	int clockNumberDriftTime[12][NUM_OF_WAVE_SAMPLES];
	int tdcNhit[12];
	c->SetBranchAddress("dd",dd);
	c->SetBranchAddress("dt",dt);
	c->SetBranchAddress("a",a);
	c->SetBranchAddress("h",h);
	c->SetBranchAddress("adc",adc);
	c->SetBranchAddress("clockNumberDriftTime",clockNumberDriftTime);
	c->SetBranchAddress("i",iHit);
	c->SetBranchAddress("n",tdcNhit);

	//===================Output file============================
	buf.str(""); buf.clear();
	buf<<"../root/fit_"<<startNo<<"_"<<nFiles<<"."<<suffix<<"root";
	TFile * f = new TFile(buf.str().c_str(),"RECREATE"); 
	TTree * t = new TTree("t","t");
	Double_t sliY, iniY;
	Double_t slY, inY;
	Double_t slerY, inerY;
	Double_t chi2Y;
	Double_t chi2iY;
	Double_t resoY;
	Double_t resoiY;
	Double_t XY;
	Double_t TY;
	Double_t dminY;
	Double_t dmaxY;
	Double_t sliX, iniX;
	Double_t slX, inX;
	Double_t slerX, inerX;
	Double_t chi2X;
	Double_t chi2iX;
	Double_t resoX;
	Double_t resoiX;
	Double_t XX;
	Double_t TX;
	Double_t dminX;
	Double_t dmaxX;
	Double_t adctestX;// ADC SUM for test channel (which we use to generate XT and reso)
	Double_t adctestY;
	Double_t adcminX;// Minimum ADC SUM for other channels (which we use to get the track)
	Double_t adcminY;
	Int_t htestX;// Peak height for test channel (which we use to generate XT and reso)
	Int_t htestY;
	Int_t hminX;// Minimum Peak height for other channels (which we use to get the track)
	Int_t hminY;
	Long64_t iEvent;
	t->Branch("patY",&patternY);
	t->Branch("itestY",&itestY);
	t->Branch("chi2Y",&chi2Y);
	t->Branch("chi2iY",&chi2iY);
	t->Branch("resoY",&resoY);
	t->Branch("resoiY",&resoiY);
	t->Branch("slerY",&slerY);
	t->Branch("inerY",&inerY);
	t->Branch("slY",&slY);
	t->Branch("inY",&inY);
	t->Branch("sliY",&sliY);
	t->Branch("iniY",&iniY);
	t->Branch("xY",&XY);
	t->Branch("tY",&TY);
	t->Branch("dminX",&dminX);
	t->Branch("dmaxX",&dmaxX);
	t->Branch("patX",&patternX);
	t->Branch("itestX",&itestX);
	t->Branch("chi2X",&chi2X);
	t->Branch("chi2iX",&chi2iX);
	t->Branch("resoX",&resoX);
	t->Branch("resoiX",&resoiX);
	t->Branch("slerX",&slerX);
	t->Branch("inerX",&inerX);
	t->Branch("slX",&slX);
	t->Branch("inX",&inX);
	t->Branch("sliX",&sliX);
	t->Branch("iniX",&iniX);
	t->Branch("xX",&XX);
	t->Branch("tX",&TX);
	t->Branch("dminX",&dminX);
	t->Branch("dmaxX",&dmaxX);
	t->Branch("aminX",&adcminX);
	t->Branch("aminY",&adcminY);
	t->Branch("atestX",&adctestX);
	t->Branch("atestY",&adctestY);
	t->Branch("hminX",&hminX);
	t->Branch("hminY",&hminY);
	t->Branch("htestX",&htestX);
	t->Branch("htestY",&htestY);
	t->Branch("iEvent",&iEvent);

	//===================Prepare for drawing============================
#ifdef PRINT
	TCanvas * ca  = new TCanvas("ca","ca",1024,768);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	TPad * p1 = new TPad("p1","p1",0,0.5,0.7,1);
	TPad * p2 = new TPad("p2","p2",0,0,0.7,0.5);
	TPad * p3 = new TPad("p3","p3",0.7,0.75,1,1);
	TPad * p4 = new TPad("p4","p4",0.7,0.5,1,0.75);
	TPad * p5 = new TPad("p5","p5",0.7,0.25,1,0.5);
	TPad * p6 = new TPad("p6","p6",0.7,0,1,0.25);
	p1->SetGridx(1);
	p1->SetGridy(1);
	p2->SetGridx(1);
	p2->SetGridy(1);
	p2->SetLogz(1);
	p3->SetGridx(1);
	p3->SetGridy(1);
	p4->SetGridx(1);
	p4->SetGridy(1);
	p5->SetGridx(1);
	p5->SetGridy(1);
	p6->SetGridx(1);
	p6->SetGridy(1);
	p1->Draw();
	p2->Draw();
	p3->Draw();
	p4->Draw();
	p5->Draw();
	p6->Draw();
	TH2D * h1 = new TH2D("h1","h1",NBINS,-0.04,0.04,NBINS,-16,16);
	h1->SetContour(50);
	h1->GetXaxis()->SetTitle("Slope");
	h1->GetYaxis()->SetTitle("Intercept");
	TH2D * hdis = new TH2D("hdis","Event Display",2048,-10,760,2048,-20,20);
	TEllipse* e1 = new TEllipse();
	TEllipse* e2 = new TEllipse();
	TEllipse* e3 = new TEllipse();
	TEllipse* e4 = new TEllipse();
	e1->SetFillStyle(0);
	e2->SetFillStyle(0);
	e2->SetLineColor(kRed);
	e3->SetFillStyle(0);
	e4->SetFillStyle(0);
	TLine * l1 = new TLine();
	TLine * l2 = new TLine();
	l2->SetLineColor(kRed);
	TEllipse * eini = new TEllipse(1,1,0.0001,0.05);
	TEllipse * efit = new TEllipse(1,1,0.0002,0.1);
	eini->SetLineColor(kYellow);
	eini->SetFillColor(kYellow);
	eini->SetFillStyle(1);
	efit->SetLineColor(kRed);
	efit->SetFillStyle(0);
	TText * text1 = new TText(-0.035,13,"");
	TText * text2 = new TText(-0.035,11,"");
	TLatex *text3 = new TLatex(NUM_OF_WAVE_SAMPLES/18, 
			MAX_WAVE_CH-1*(MAX_WAVE_CH-MIN_WAVE_CH)/10,
			"");
	TLatex *textTDC[NHITS_TO_CHOOSET][NUM_OF_WAVE_SAMPLES];
	TMarker *markerTDC[NHITS_TO_CHOOSET][NUM_OF_WAVE_SAMPLES];
	for (Int_t i=0; i<NHITS_TO_CHOOSET; i++) {
		for (Int_t j=0; j<NUM_OF_WAVE_SAMPLES; j++) {
			textTDC[i][j] = new TLatex(0,MIN_WAVE_CH+10,"");
			textTDC[i][j]->SetTextSize(0.04);
			markerTDC[i][j] = new TMarker(0,0,20);
			markerTDC[i][j]->SetMarkerSize(0.7);
		}
	}
	TGraph *gr_waveForm[4];
	Int_t vSample[NUM_OF_WAVE_SAMPLES];
	for (Int_t i=0; i<NUM_OF_WAVE_SAMPLES; i++){
		vSample[i] = i;
	}
#endif

	//===================Prepare for fitting============================
	TMinuit *gMinuit = 0;
	Double_t arglist[10];
	Int_t ierflg = 0;
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;

	//===================Loop in Events============================
	Long64_t N = c->GetEntries();
	if (nEventMax&&nEventMax<N) N = nEventMax;
	printf("Processing... %d Events\n",N);
	for (iEvent = 0;iEvent<N; iEvent++){
		if (iEvent%1000==0) printf("%lf%\n",(double)iEvent/N*100);
		c->GetEntry(iEvent);
		//printf("In Event(%d)\n",iEvent);
		//===================Get Hits============================
		getHits();
		if (patternX<0&&patternY<0) continue;
		dminX = 1e9;
		dmaxX = -1;
		dminY = 1e9;
		dmaxY = -1;
		hminX = 1e9;
		hminY = 1e9;
		adctestX = a[itestX][iHit[itestX]];
		adctestY = a[itestY][iHit[itestY]];
		htestX = h[itestX][iHit[itestX]];
		htestY = h[itestY][iHit[itestY]];
		for ( int i = 0; i<NHITS_TO_CHOOSE; i++){
			double dtempX,dtempY;
			double adctempX,adctempY;
			double htempX,htempY;
			dtempX = dd[indiceX[i]][iHit[indiceX[i]]];
			dtempY = dd[indiceY[i]][iHit[indiceY[i]]];
			adctempX = a[indiceX[i]][iHit[indiceX[i]]];
			adctempY = a[indiceY[i]][iHit[indiceY[i]]];
			htempX = h[indiceX[i]][iHit[indiceX[i]]];
			htempY = h[indiceY[i]][iHit[indiceY[i]]];
			if (dminX>dtempX) dminX=dtempX;
			if (dmaxX<dtempX) dmaxX=dtempX;
			if (dminY>dtempY) dminY=dtempY;
			if (dmaxY<dtempY) dmaxY=dtempY;
			if (adcminX>adctempX) adcminX=adctempX;
			if (adcminY>adctempY) adcminY=adctempY;
			if (hminX>htempX) hminX=htempX;
			if (hminY>htempY) hminY=htempY;
		}

		//===================Set Initial============================
		iniX = dd[indiceX[0]][iHit[indiceX[0]]];
		if (patternX==0||patternX==3) iniX*=-1;
		if (patternX<=1){
			sliX = (dd[indiceX[1]][iHit[indiceX[1]]]-dd[indiceX[0]][iHit[indiceX[0]]])/(z[indiceX[1]]-z[indiceX[0]]);
			if (patternX==0) sliX*=-1;
		}
		else{
			sliX = (dd[indiceX[1]][iHit[indiceX[1]]]+dd[indiceX[0]][iHit[indiceX[0]]])/(z[indiceX[1]]-z[indiceX[0]]);
			if (patternX==2) sliX*=-1;
		}
		iniY = dd[indiceY[0]][iHit[indiceY[0]]];
		if (patternY==0||patternY==3) iniY*=-1;
		if (patternY<=1){
			sliY = (dd[indiceY[1]][iHit[indiceY[1]]]-dd[indiceY[0]][iHit[indiceY[0]]])/(z[indiceY[1]]-z[indiceY[0]]);
			if (patternY==0) sliY*=-1;
		}
		else{
			sliY = (dd[indiceY[1]][iHit[indiceY[1]]]+dd[indiceY[0]][iHit[indiceY[0]]])/(z[indiceY[1]]-z[indiceY[0]]);
			if (patternY==2) sliY*=-1;
		}

		//===================Get Minimum for X Tracker ============================
		fit_channel = 0;
		if(gMinuit) delete gMinuit;
		gMinuit = new TMinuit(NPAR);  //initialize TMinuit with a maximum of 5 params
		gMinuit->SetFCN(fcn);
		arglist[0] = 0;
		arglist[1] = 0;
		gMinuit->SetPrintLevel(-1); // no print
		gMinuit->mnexcm("SET NOW", arglist ,1,ierflg); // no warning 

		arglist[0] = 1;
		gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

		// Set starting values and step sizes for parameters
		gMinuit->mnparm(0, "slope", sliX, 0.00001, sliX-3*fabs(sliX),sliX+3*fabs(sliX),ierflg);
		gMinuit->mnparm(1, "intercept", iniX, 0.0001, iniX-3*fabs(iniX),iniX+fabs(iniX),ierflg);

		// Now ready for minimization step
		arglist[0] = 500.0;
		arglist[1] = 1.0;
		gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

		// Print results
		gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
		//printf("====Rrestul====\n");
		//gMinuit->mnprin(3,amin);
		gMinuit->GetParameter(0, slX, slerX);
		gMinuit->GetParameter(1, inX, inerX);

		//===================Get Minimum for Y Tracker ============================
		fit_channel = 1;
		if(gMinuit) delete gMinuit;
		gMinuit = new TMinuit(NPAR);  //initialize TMinuit with a maximum of 5 params
		gMinuit->SetFCN(fcn);
		arglist[0] = 0;
		arglist[1] = 0;
		gMinuit->SetPrintLevel(-1); // no print
		gMinuit->mnexcm("SET NOW", arglist ,1,ierflg); // no warning 

		arglist[0] = 1;
		gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

		// Set starting values and step sizes for parameters
		gMinuit->mnparm(0, "slope", sliY, 0.00001, sliY-3*fabs(sliY),sliY+3*fabs(sliY),ierflg);
		gMinuit->mnparm(1, "intercept", iniY, 0.0001, iniY-3*fabs(iniY),iniY+fabs(iniY),ierflg);

		// Now ready for minimization step
		arglist[0] = 500.0;
		arglist[1] = 1.0;
		gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

		// Print results
		gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
		//printf("====Rrestul====\n");
		//gMinuit->mnprin(3,amin);
		gMinuit->GetParameter(0, slY, slerY);
		gMinuit->GetParameter(1, inY, inerY);

		//===================Get chi2 and reso, X, T============================
		fit_channel = 0;
		*par = sliX;
		*(par+1) = iniX;
		getchi2(chi2iX,par);
		resoiX = func(z[itestX],y[itestX],par)-dd[itestX][iHit[itestX]];
		*par = slX;
		*(par+1) = inX;
		getchi2(chi2X,par);
		resoX = func(z[itestX],y[itestX],par)-dd[itestX][iHit[itestX]];
		XX=dd[itestX][iHit[itestX]]+resoX;
		TX=dt[itestX][iHit[itestX]];

		fit_channel = 1;
		*par = sliY;
		*(par+1) = iniY;
		getchi2(chi2iY,par);
		resoiY = func(z[itestY],y[itestY],par)-dd[itestY][iHit[itestY]];
		*par = slY;
		*(par+1) = inY;
		getchi2(chi2Y,par);
		resoY = func(z[itestY],y[itestY],par)-dd[itestY][iHit[itestY]];
		XY=dd[itestY][iHit[itestY]]+resoY;
		TY=dt[itestY][iHit[itestY]];

#ifdef PRINT
		//===================Draw X Tracker============================
//		if (patternX>=0){
		if (patternX>=0&&fabs(resoX)>0.6){
			fit_channel = 0;
			p1->cd();
			buf.str("");buf.clear();
			buf<<"Pattern["<<patternX<<"], d_{drift}("<<dd[indiceX[0]][iHit[indiceX[0]]]<<","<<dd[itestX][iHit[itestX]]<<","<<dd[indiceX[1]][iHit[indiceX[1]]]<<","<<dd[indiceX[2]][iHit[indiceX[2]]]
				<<") mm, ini("<<sliX<<","<<iniX
				<<"), fit("<<slX<<","<<inX<<")"
				<<", ADC("<<a[indiceX[0]][iHit[indiceX[0]]]<<","<<a[itestX][iHit[itestX]]<<","<<a[indiceX[1]][iHit[indiceX[1]]]<<","<<a[indiceX[2]][iHit[indiceX[2]]]<<")";
			hdis->SetTitle(buf.str().c_str());
			hdis->Draw();
			l1->SetX1((-20-iniX)/sliX);
			l1->SetY1(-20);
			l1->SetX2((20-iniX)/sliX);
			l1->SetY2(20);
			l2->SetX1((-20-inX)/slX);
			l2->SetY1(-20);
			l2->SetX2((20-inX)/slX);
			l2->SetY2(20);
			e1->SetX1(z[indiceX[0]]);e1->SetY1(y[indiceX[0]]);e1->SetR1(dd[indiceX[0]][iHit[indiceX[0]]]);e1->SetR2(dd[indiceX[0]][iHit[indiceX[0]]]);
			e2->SetX1(z[itestX]);e2->SetY1(y[itestX]);e2->SetR1(dd[itestX][iHit[itestX]]);e2->SetR2(dd[itestX][iHit[itestX]]);
			e3->SetX1(z[indiceX[1]]);e3->SetY1(y[indiceX[1]]);e3->SetR1(dd[indiceX[1]][iHit[indiceX[1]]]);e3->SetR2(dd[indiceX[1]][iHit[indiceX[1]]]);
			e4->SetX1(z[indiceX[2]]);e4->SetY1(y[indiceX[2]]);e4->SetR1(dd[indiceX[2]][iHit[indiceX[2]]]);e4->SetR2(dd[indiceX[2]][iHit[indiceX[2]]]);
			e1->Draw("SAME");
			e2->Draw("SAME");
			e3->Draw("SAME");
			e4->Draw("SAME");
			l1->Draw("SAME");
			l2->Draw("SAME");

			p2->cd();
			//===================Check chi2============================
			buf.str("");buf.clear();
			buf<<"Event #"<<iEvent;
			h1->SetTitle(buf.str().c_str());
			Double_t chisq;
			for(int i = 1; i<=NBINS; i++){
				for(int j = 1; j<=NBINS; j++){
					*par = h1->GetXaxis()->GetBinCenter(i);
					*(par+1) = h1->GetYaxis()->GetBinCenter(j);
					getchi2(chisq,par);
					h1->SetBinContent(i,j,chisq);
					//printf("h[%d,%d]=%lf\n",i,j,chisq);
				}
			}
			h1->Draw("COLZ");
			eini->SetX1(sliX); eini->SetY1(iniX);
			efit->SetX1(slX); efit->SetY1(inX); efit->SetR1(slerX); efit->SetR2(inerX);
			eini->Draw("SAME");
			efit->Draw("SAME");
			buf.str("");buf.clear();
			buf<<"initial: slope = "<<sliX<<", intercept = "<<iniX<<", chi2 = "<<chi2iX<<", reso = "<<resoiX;
			text1->SetText(-0.035,13,buf.str().c_str());
			text1->SetTextSize(0.025);
			text1->Draw("SAME");
			buf.str("");buf.clear();
			buf<<"after fit: slope = "<<slX<<", intercept = "<<inX<<", chi2 = "<<chi2X<<", reso = "<<resoX;
			text2->SetText(-0.035,11,buf.str().c_str());
			text2->SetTextSize(0.025);
			text2->Draw("SAME");

			int index = 0;
			for (int i = 0; i<4; i++){
				if (i==0) p3->cd();
				else if (i==1) p4->cd();
				else if (i==2) p5->cd();
				else if (i==3) p6->cd();
				if (i==3) index = itestX;
				else index = indiceX[i];
				gr_waveForm[i] = new TGraph(NUM_OF_WAVE_SAMPLES,vSample,adc[index]);
				if (i<=1)
					gr_waveForm[i]->GetHistogram()->SetTitle(Form("Waveform: Tracker X up, wireID=%d, event=%d",index%3,iEvent));
				else
					gr_waveForm[i]->GetHistogram()->SetTitle(Form("Waveform: Tracker X down, wireID=%d, event=%d",index%3,iEvent));
				gr_waveForm[i]->GetHistogram()->SetMinimum(MIN_WAVE_CH); 
				gr_waveForm[i]->GetHistogram()->SetMaximum(MAX_WAVE_CH);
				gr_waveForm[i]->GetHistogram()->SetXTitle("Samplig ID");
				gr_waveForm[i]->GetHistogram()->SetYTitle("ADC(ch)");
				gr_waveForm[i]->SetMarkerStyle(24);
				gr_waveForm[i]->SetMarkerSize(0.8);
				gr_waveForm[i]->Draw("apwl");
				//text3->SetText(Form("q=%d",q[index]));
				//text3->Draw();
				for (Int_t j=0; j<tdcNhit[index]; j++) {
					textTDC[i][j]->SetText(clockNumberDriftTime[index][j],MIN_WAVE_CH+10,Form("%d",dt[index][j]));
					textTDC[i][j]->Draw();
					markerTDC[i][j]->SetX(clockNumberDriftTime[index][j]);
					markerTDC[i][j]->SetY(adc[index][clockNumberDriftTime[index][j]]);
					//printf("markerX = clockNumberDriftTime[%d][%d] = %d\n",index,j,clockNumberDriftTime[index][j]);
					//printf("markerY = adc[%d][%d] = %d\n",index,clockNumberDriftTime[index][j],adc[index][clockNumberDriftTime[index][j]]);
					if (a[index][j]>Amin[index]) markerTDC[i][j]->SetMarkerColor(kRed);
					else markerTDC[i][j]->SetMarkerColor(kBlue);
					markerTDC[i][j]->Draw();
				}
			}

			buf.str("");buf.clear();
			buf<<"chi2.X."<<startNo<<"_"<<nFiles<<"."<<iEvent<<"."<<suffix<<"pdf";
			ca->SaveAs(buf.str().c_str());
			buf.str("");buf.clear();
			buf<<"chi2.X."<<startNo<<"_"<<nFiles<<"."<<iEvent<<"."<<suffix<<"png";
			ca->SaveAs(buf.str().c_str());
		}
#endif

		t->Fill();
	}
	t->Write();
	f->Close();

	return 0;
}
