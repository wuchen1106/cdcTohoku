#include <iostream>
#include <iomanip>  
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"
#include "TEllipse.h"
#include "TChain.h"
#include "TLine.h"
#include "TText.h"
#include "TStyle.h"
#include "TPad.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TVector3.h"

//#define PRINTR
//#define EVENTDISPLAY
#define GETADC

#define NUM_OF_WAVE_SAMPLES 32
#define t2xDictNTypes  2
#define t2xDictNPoints 16
#define NUM_OF_WAVE_SAMPLES 32
#define NHITS_TO_CHOOSE 4
#define NHITS_TO_CHOOSET 4
#define MIN_WAVE_CH 200
#define MAX_WAVE_CH 700
#define MAX_WAVE_CH2 700
#define NBINS  256

double PI = TMath::Pi();
int power2_15 = pow(2,15);

// ________About Hit Pattern___________
Int_t patternX = 0;
Int_t patternY = 0;
Int_t indiceX[NHITS_TO_CHOOSE] = {0};
Int_t indiceY[NHITS_TO_CHOOSE] = {0};

// ________About XT___________
TF1 * XTL_func[12];
TF1 * XTR_func[12];

// ________About Cell___________
Double_t U=8;
Double_t Z1[66],Y1[66],X1[66],errord[66];// left hand side
Double_t Z2[66],Y2[66],X2[66]; // right hand side

// ________About Track___________
bool flag[66];
Double_t dt[66];
Double_t dd[66];
int h[66];
int tdcNhit[66];

// ________About Fitting___________
Int_t fit_channel = 0; // 0: X Tracker; 1: Y Tracker
TMinuit *gMinuit = 0;
Double_t arglist[10];
Int_t ierflg = 0;
Double_t amin,edm,errdef;
Int_t nvpar,nparx,icstat;

//______________________________________________________________________________
Double_t get_dist(float x,float y,Double_t sl, Double_t in);
void getchi2(Double_t &f, Double_t sl, Double_t in);
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
Double_t t2x(double t, int iwire, int side);
void do_fit(Double_t sli, Double_t ini);
void getHits();
void print_usage(char* prog_name);
int chg2cht(int i);
int chg2chp2(int i);
int chg2chp3(int i);
int cht2chg(int i);
int chp22chg(int i);
int chp32chg(int i);
int get_bid(int i);
int get_bid_core(int i);
//______________________________________________________________________________
int main(int argc, char** argv){
	if (argc<2){
		print_usage(argv[0]);
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);
	std::string suffix = "";
	if (argc>=3){
		suffix  = argv[2];
		suffix=suffix+".";
	}
	int nEventMax = 0;
	if (argc>=4) nEventMax = (int)strtol(argv[3],NULL,10);
	int iterationNo = 0;
	if (argc>=5) iterationNo = (int)strtol(argv[4],NULL,10);
	int N_MIN = 0;
	if (argc>=6) N_MIN = (int)strtol(argv[5],NULL,10);
	double chi2Xmax,chi2Ymax;
	chi2Xmax=1e9;
	if (argc>=7) chi2Xmax = (int)strtol(argv[6],NULL,10);
	chi2Ymax=chi2Xmax;

	double eff_p2_average = 0;
	double eff_p3_average = 0;

	//===================Get HV============================
	//===================Get npair============================
	TFile * if_run = new TFile("../run_summary/run.eff.root");
	TTree * t_run = (TTree*) if_run->Get("t");
	double i_runNo, HVTXY1,HVTXY2,HVP2,HVP3,THRTXY1,THRTXY2,THRP2,THRP3;
	t_run->SetBranchAddress("runNo",&i_runNo);
	t_run->SetBranchAddress("HVTXY1",&HVTXY1);
	t_run->SetBranchAddress("HVTXY2",&HVTXY2);
	t_run->SetBranchAddress("HVP2",&HVP2);
	t_run->SetBranchAddress("HVP3",&HVP3);
	t_run->SetBranchAddress("THRTXY1",&THRTXY1);
	t_run->SetBranchAddress("THRTXY2",&THRTXY2);
	t_run->SetBranchAddress("THRP2",&THRP2);
	t_run->SetBranchAddress("THRP3",&THRP3);
	for(int i = 0; i<t_run->GetEntries(); i++){
		t_run->GetEntry(i);
		if (i_runNo == runNo) break;
	}
	std::cout<<"runNo#"<<runNo<<": "<<(int)HVTXY1<<" V"<<std::endl;
	TString gastype = "He:CH_{4}(80:20)";
	double npair = 17.96;
	if (runNo>=227||runNo<=147){
		gastype = "He:iC_{4}H_{10}(90:10)";
		npair = 27.96;
	}
	else if (runNo<=194){
		gastype = "He:C_{2}H_{4}(50:50)";
		npair = 56.10;
	}

	//===================Get a2c============================
	TF1 * f_a2c = new TF1("a2c","5.98739+2.6652*x+0.000573394*x*x-5.21769e-05*x*x*x+3.05897e-07*x*x*x*x-7.54057e-10*x*x*x*x*x+8.60252e-13*x*x*x*x*x*x-3.68603e-16*x*x*x*x*x*x*x",-10,800);

	//===================Get XT============================
	std::vector<std::vector<double> > v_xt_dr;
	std::vector<std::vector<double> > v_xt_dl;
	std::vector<std::vector<double> > v_xt_tr;
	std::vector<std::vector<double> > v_xt_tl;
	v_xt_tl.resize(66);
	v_xt_tr.resize(66);
	v_xt_dr.resize(66);
	v_xt_dl.resize(66);
	double xt_sigma_tr = 7.5;
	double xt_sigma_p = 30;
	double xt_sigma_tr2 = 1;
	double xt_sigma_p2 = 1;
	// get xt
	TChain * c_para = new TChain("para","para");
	std::stringstream buf;
	buf.str("");
	buf.clear();
	buf<<"../info/info."<<runNo<<"."<<suffix<<iterationNo<<".root";
	std::cout<<"Adding \""<<buf.str()<<"\""<<std::endl;
	TChain * ch_xt = new TChain("xt","xt");
	ch_xt->Add(buf.str().c_str());
	double xt_t, xt_d, xt_sigma;
	int xt_i, xt_n;
	ch_xt->SetBranchAddress("t",&xt_t);
	ch_xt->SetBranchAddress("d",&xt_d);
	ch_xt->SetBranchAddress("i",&xt_i);
	ch_xt->SetBranchAddress("n",&xt_n);
	ch_xt->SetBranchAddress("sig",&xt_sigma);
	for ( int i = 0; i<ch_xt->GetEntries(); i++){
		ch_xt->GetEntry(i);
		//FIXME
		int cht = chg2cht(xt_i);
		if (cht<0) continue; // only take tracker
		if (xt_n<N_MIN) continue;
		if (xt_sigma>xt_sigma_tr) continue;
		if (xt_d<=0){
			v_xt_tl[xt_i].push_back(xt_t);
			v_xt_dl[xt_i].push_back(xt_d);
		}
		if (xt_d>=0){
			v_xt_tr[xt_i].push_back(xt_t);
			v_xt_dr[xt_i].push_back(xt_d);
		}
	}
	TGraph * gxtl[66];
	TGraph * gxtr[66];
	for ( int i = 0; i<66; i++){
		gxtl[i] = 0;
		gxtr[i] = 0;
	}
	for ( int i = 0 ; i<12; i++){
		int chg = cht2chg(i);
		v_xt_dl[chg].push_back(0); v_xt_tl[chg].push_back(0);
		v_xt_dr[chg].push_back(0); v_xt_tr[chg].push_back(0);
		TGraph * gxtl = new TGraph(v_xt_dl[chg].size(),&v_xt_tl[chg][0],&v_xt_dl[chg][0]);
		XTL_func[i] = new TF1(Form("xtl_%d",i),"pol5",0,250);
		gxtl->Fit(Form("xtl_%d",i),"qN0","");
		gxtl->Print();
		TGraph * gxtr = new TGraph(v_xt_dr[chg].size(),&v_xt_tr[chg][0],&v_xt_dr[chg][0]);
		XTR_func[i] = new TF1(Form("xtr_%d",i),"pol5",0,250);
		gxtr->Fit(Form("xtr_%d",i),"qN0","");
	}

	//===================Set Geometry============================
	// The z values
	// get the inital values
	TFile * if_geom = new TFile("../info/wire-position.root");
	TTree * t_geom = (TTree*) if_geom->Get("t");
	double x1_geom, y1_geom, z1_geom, x2_geom, y2_geom, z2_geom;
	int id_geom;
	t_geom->SetBranchAddress("wire_id",&id_geom);
	t_geom->SetBranchAddress("x1",&x1_geom);
	t_geom->SetBranchAddress("y1",&y1_geom);
	t_geom->SetBranchAddress("z1",&z1_geom);
	t_geom->SetBranchAddress("x2",&x2_geom);
	t_geom->SetBranchAddress("y2",&y2_geom);
	t_geom->SetBranchAddress("z2",&z2_geom);
	for(int i = 0; i<t_geom->GetEntries(); i++){
		t_geom->GetEntry(i);
		if (id_geom<0||id_geom>=66){
			printf("WARNING: Entry[%d], wireID = %d\n",i,id_geom);
		}
		X1[id_geom] = x1_geom;
		Y1[id_geom] = y1_geom;
		Z1[id_geom] = z1_geom+378;
		X2[id_geom] = x2_geom;
		Y2[id_geom] = y2_geom;
		Z2[id_geom] = z2_geom+378;
		//FIXME
		if (id_geom>=33&&id_geom<60){
			Y1[id_geom] += 0.625;
			Y2[id_geom] += 0.625;
		}
		else if (id_geom>=6&&id_geom<32){
			Y1[id_geom] += 0;
			Y2[id_geom] += 0;
		}
	}
	if_geom->Close();
	// The errors on z values
	Float_t error = 0.2;
	for ( int i = 0; i<66; i++){
		errord[i]=error;
	}

	//===================Prepare Histograms============================
	//chi2
	TH1D * hist_chi2X = new TH1D("chi2X","chi2X",128,0,20);
	TH1D * hist_chi2Y = new TH1D("chi2Y","chi2Y",128,0,20);
	TH1D * hist_effX = new TH1D("effX","effX",128,0,20);
	TH1D * hist_effY = new TH1D("effY","effY",128,0,20);
	TH1D * hist_slX = new TH1D("slX","slX",128,-0.02,0.02);
	TH1D * hist_slY = new TH1D("slY","slY",128,-0.02,0.02);

	TH2D * hist_wf0[66];
	TH2D * hist_wf1[66];
	TH2D * hist_wf2[66];
	TH2D* hist_gg[66];
	TGaxis * axis = new TGaxis(-U*sqrt(2),0,-U*sqrt(2),6,1,1e6,50510,"G");
	axis->SetTitle("Gas Gain");
	for ( int i = 0; i<66; i++){
		hist_wf0[i] = new TH2D(Form("wf0_%d",i),Form("wave form in ch%d (good event)",i),NUM_OF_WAVE_SAMPLES*2,-NUM_OF_WAVE_SAMPLES,NUM_OF_WAVE_SAMPLES,600,100,700);
		hist_wf1[i] = new TH2D(Form("wf1_%d",i),Form("wave form in ch%d",i),NUM_OF_WAVE_SAMPLES*2,-NUM_OF_WAVE_SAMPLES,NUM_OF_WAVE_SAMPLES,600,100,700);
		hist_wf2[i] = new TH2D(Form("wf2_%d",i),Form("wave form in ch%d (hit but no pass)",i),NUM_OF_WAVE_SAMPLES*2,-NUM_OF_WAVE_SAMPLES,NUM_OF_WAVE_SAMPLES,600,100,700);
		hist_gg[i] = new TH2D(Form("gg2_%d",i),Form("Gas Gain in ch%d",i),128,-U*sqrt(2),U*sqrt(2),128,0,6);
		hist_gg[i]->GetYaxis()->SetTickLength(0);
		hist_gg[i]->GetYaxis()->SetTitleOffset(3);
		hist_gg[i]->GetYaxis()->SetLabelOffset(3);
	}
	TH2D *hist_wf_p2[8];
	TH2D *hist_wf_p3[8];
	for (int i = 0; i<8; i++){
		hist_wf_p2[i] = new TH2D(Form("wf0p2_%dmm",i),Form("%d~%dmm",i,i+1),NUM_OF_WAVE_SAMPLES*2,-NUM_OF_WAVE_SAMPLES,NUM_OF_WAVE_SAMPLES,600,100,700);
		hist_wf_p3[i] = new TH2D(Form("wf0p3_%dmm",i),Form("%d~%dmm",i,i+1),NUM_OF_WAVE_SAMPLES*2,-NUM_OF_WAVE_SAMPLES,NUM_OF_WAVE_SAMPLES,600,100,700);
	}

	TH1D* adcpeak_p2[8];
	TH1D* adcpeak_p3[8];
	TH1D* adcsum_p2[8];
	TH1D* adcsum_p3[8];
	for (int i = 0; i<8; i++){
		adcpeak_p2[i] = new TH1D(Form("adcpeakp2_%dmm",i),Form("%d~%dmm",i,i+1),MAX_WAVE_CH2/2,0,MAX_WAVE_CH2);
		adcpeak_p3[i] = new TH1D(Form("adcpeakp3_%dmm",i),Form("%d~%dmm",i,i+1),MAX_WAVE_CH2/2,0,MAX_WAVE_CH2);
		adcsum_p2[i] = new TH1D(Form("adcsump2_%dmm",i),Form("%d~%dmm",i,i+1),500,-100,1900);
		adcsum_p3[i] = new TH1D(Form("adcsump3_%dmm",i),Form("%d~%dmm",i,i+1),500,-100,1900);
	}

	TLatex * text_reso1[8];
	TLatex * text_reso2[8];
	TLatex * text_reso3[8];
	for ( int i = 0; i<8; i++){
		text_reso1[i] = new TLatex();
		text_reso1[i]->SetTextSize(0.03);
		text_reso2[i] = new TLatex();
		text_reso2[i]->SetTextSize(0.03);
		text_reso3[i] = new TLatex();
		text_reso3[i]->SetTextSize(0.03);
	}

	TLatex * text = new TLatex();
	text->SetTextSize(0.02);
	TLatex * text2 = new TLatex();
	text2->SetTextSize(0.02);
	text2->SetText(0.1,0.96,gastype+Form(" Tracker(%dV,%d) Proto2(%dV,%d) Proto3(%dV,%d)",(int)HVTXY1,(int)THRTXY1,(int)HVP2,(int)THRP2,(int)HVP3,(int)THRP3));
	TString name;
	TString title;
	TString chamberName;
	TH2D* xt_tr[12];
	TH2D* xt_p2[9];
	TH2D* xt_p3[9];
	TH2D* spot[4];
	TH2D* dt123[16];
	TH1D* reso_p2[8];
	TH1D* reso_p3[8];
	double Neffi_p2[8]={0};
	double Neffi_p3[8]={0};

	int TMAX = 280; // ns
	int TMAX2 = 800; // ns
	int DDSTEP = 32; // steps
	int DDSTEP2 = 21; // steps
	int ICH = 1;
	std::vector<TH1D*> DTHISTS;
	std::vector<TH1D*> DTHISTS2;
	for ( int i  =0;  i< DDSTEP*2; i++){
		for ( int j  =0;  j< 66; j++){
			name  = Form("dt_%d_%d",i,j);
			title = name;
			double t;
			if (j<=5||j>=60) t = TMAX;
			else t = TMAX2;
			TH1D * dthist = new TH1D(name,title,(50+t)/2,-50,t);
			DTHISTS.push_back(dthist);
			name  = Form("dt2_%d_%d",i,j);
			title = name;
			TH1D * dthist2 = new TH1D(name,title,256,-U-2,0);
			DTHISTS2.push_back(dthist2);
		}
	}

	for (int i = 0; i<12; i++){
		name = Form("xt_tr_all_%d",i);
		if (i<3) chamberName = "X Tracker Up";
		else if (i<6) chamberName = "Y Tracker Up";
		else if (i<9) chamberName = "Y Tracker Down";
		else if (i<12) chamberName = "X Tracker Down";
		int wireID = i;
		if (i>=6) wireID+= 54;
		title = Form(chamberName+", wire#%d",wireID);
		xt_tr[i] = new TH2D(name,title,128,-9,9,160,-20,300);
		xt_tr[i]->GetXaxis()->SetTitle("X [mm]");
		xt_tr[i]->GetYaxis()->SetTitle("T [ns]");
	}
	for (int i = 0; i<9; i++){
		name = Form("xt_p3_all_%d",i);
		int wireID;
		if (i<3) wireID = 15-i%3;
		else if (i<6) wireID = 20-i%3;
		else if (i<9) wireID = 26-i%3;
		title = Form("wire#%d",wireID);
		xt_p3[i] = new TH2D(name,title,128,-10,10,420,-40,800);
		xt_p3[i]->GetXaxis()->SetTitle("X [mm]");
		xt_p3[i]->GetYaxis()->SetTitle("T [ns]");
	}
	for (int i = 0; i<9; i++){
		name = Form("xt_p2_all_%d",i);
		int wireID;
		if (i<3) wireID = 42-i%3;
		else if (i<6) wireID = 47-i%3;
		else if (i<9) wireID = 53-i%3;
		title = Form("wire#%d",wireID);
		xt_p2[i] = new TH2D(name,title,128,-10,10,420,-40,800);
		xt_p2[i]->GetXaxis()->SetTitle("X [mm]");
		xt_p2[i]->GetYaxis()->SetTitle("T [ns]");
	}
	for (int i = 0; i<4; i++){
		name = Form("spot_%d",i);
		if (i==0) chamberName = "X Tracker Up";
		else if (i==1) chamberName = "Prototype 3";
		else if (i==2) chamberName = "Prototype 2";
		else if (i==3) chamberName = "X Tracker Down";
		title = "At "+chamberName;
		spot[i] = new TH2D(name,title,64,-10,10,64,-10,10);
		spot[i]->GetXaxis()->SetTitle("X [mm]");
		spot[i]->GetYaxis()->SetTitle("Y [mm]");
	}
	for (int i = 0; i<16; i++){
		name = Form("dt123_%d",i);
		if (i<8)
			dt123[i] = new TH2D(name,name,70,0,280,70,0,280);
		else
			dt123[i] = new TH2D(name,name,64,0,8,64,0,8);
	}
	for (int i = 0; i<8; i++){
		name = Form("reso_p2_%d",i);
		reso_p2[i] = new TH1D(name,name,128,-2,2);
		name = Form("reso_p3_%d",i);
		reso_p3[i] = new TH1D(name,name,128,-2,2);
	}

#ifdef EVENTDISPLAY
	// Event Display
	TCanvas * canvas2  = new TCanvas("canvas2","canvas2",2048,1024);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);

	TPad * pad2[12];
	pad2[0] = new TPad("p1","p1",0,0.5,0.4,1);
	pad2[1] = new TPad("p2","p2",0,0,0.4,0.5);
	pad2[2] = new TPad("p3","p3",0.4,0.75,0.6,1);
	pad2[3] = new TPad("p4","p4",0.4,0.5,0.6,0.75);
	pad2[4] = new TPad("p5","p5",0.4,0.25,0.6,0.5);
	pad2[5] = new TPad("p6","p6",0.4,0,0.6,0.25);
	pad2[6] = new TPad("p7","p7",0.6,0.75,0.8,1);
	pad2[7] = new TPad("p8","p8",0.6,0.5, 0.8,0.75);
	pad2[8] = new TPad("p9","p9",0.6,0.25,0.8,0.5);
	pad2[9] = new TPad("p10","p10",0.6,0,   0.8,0.25);
	pad2[10]= new TPad("p11","p11",0.8,0.5,1,1);
	pad2[11]= new TPad("p12","p12",0.8,0, 1,0.5);
	for ( int i = 0; i<12; i++){
		pad2[i]->SetGridx(1);
		pad2[i]->SetGridy(1);
		pad2[i]->Draw();
	}
	TH2D * hdisX = new TH2D("hdisX","Event Display",2048,-10,760,2048,-20,20);
	TH2D * hdisY = new TH2D("hdisY","Event Display",2048,-10,760,2048,-20,20);
	TEllipse* ellipse[12];
	TText * textEll[12];
	for (int i = 0; i<12; i++){
		ellipse[i] = new TEllipse();
		textEll[i] = new TText(-0.035,11,"");
		if (i==4||i==5||i==10||i==11){
			ellipse[i]->SetLineColor(kGreen);
			textEll[i]->SetTextColor(kBlue);
		}
		else{
			ellipse[i]->SetFillStyle(0);
			textEll[i]->SetTextColor(kRed);
		}
		textEll[i]->SetTextSize(0.03);
	}
	TLine * lineX = new TLine();
	TLine * lineY = new TLine();
	TLatex *textWF[12];
	TLatex *textTDC[12][NUM_OF_WAVE_SAMPLES];
	TMarker *markerTDC[12][NUM_OF_WAVE_SAMPLES];
	TLine *line_dd[12];
	TLine *line_fd[12];
	for (Int_t i=0; i<12; i++) {
		line_fd[i] = new TLine();
		line_dd[i] = new TLine();
		line_fd[i]->SetLineColor(kGreen);
		line_dd[i]->SetLineColor(kRed);
		line_dd[i]->SetLineWidth(0.5);
		line_fd[i]->SetLineWidth(0.5);
		textWF[i] = new TLatex(0,0,"");
		textWF[i]->SetTextSize(0.04);
		textWF[i]->SetTextColor(kRed);
		for (Int_t j=0; j<NUM_OF_WAVE_SAMPLES; j++) {
			textTDC[i][j] = new TLatex(0,MIN_WAVE_CH+10,"");
			textTDC[i][j]->SetTextSize(0.04);
			textTDC[i][j]->SetTextColor(kRed);
			markerTDC[i][j] = new TMarker(0,0,20);
			markerTDC[i][j]->SetMarkerSize(0.55);
		}
	}
	TGraph *gr_waveForm[12];
	Int_t vSample[NUM_OF_WAVE_SAMPLES];
	for (Int_t i=0; i<NUM_OF_WAVE_SAMPLES; i++){
		vSample[i] = i;
	}
#endif


	//===================Get ROOT File============================
	TChain * c = new TChain("t","t");
	buf.str(""); buf.clear();
	buf<<"../root/d_"<<runNo<<".root";
	c->Add(buf.str().c_str());
	std::cout<<"Adding \""<<buf.str()<<"\""<<std::endl;
	double aa[66];
	Int_t ihit[66];
	c->SetBranchAddress("dt",dt);
	c->SetBranchAddress("aa",aa);
	c->SetBranchAddress("h",h);
	c->SetBranchAddress("n",tdcNhit);
	c->SetBranchAddress("i",ihit);

#ifdef GETADC
	TChain * chain2 = new TChain("tree","tree");
	buf.str(""); buf.clear();
	buf<<"../root/run_000"<<runNo<<"_built.root";
	chain2->Add(buf.str().c_str());
	int adc[66][32];
	int tdc[66][32];
	int clockNumberDriftTime[66][32];
	int tdcNhit[66];
	chain2->SetBranchAddress("adc",adc);
	chain2->SetBranchAddress("driftTime",tdc);
	chain2->SetBranchAddress("clockNumberDriftTime",clockNumberDriftTime);
	chain2->SetBranchAddress("tdcNhit",tdcNhit);
#endif

	//===================Output file============================
	buf.str(""); buf.clear();
	buf<<"../root/fit."<<runNo<<"."<<suffix<<iterationNo<<".root";
	TFile * f = new TFile(buf.str().c_str(),"RECREATE"); 
	TTree * t = new TTree("t","t");
	Double_t sliY, iniY;
	Double_t sliX, iniX;
	Double_t chi2iY;
	Double_t chi2iX;

	Long64_t iEvent;
	Double_t slY, inY;
	Double_t slX, inX;
	Double_t slerY, inerY;
	Double_t slerX, inerX;
	Double_t chi2Y;
	Double_t chi2X;

	Double_t o_dist[66];
	Double_t fd[66];

	t->Branch("iEvent",&iEvent);

	t->Branch("dt",&dt,"dt[66]/D");
	t->Branch("h",&h,"h[66]/I");
	t->Branch("aa",&aa,"aa[66]/D");
	t->Branch("i",&ihit,"i[66]/I");

	t->Branch("dd",&o_dist,"dd[66]/D");

	t->Branch("patX",&patternX);
	t->Branch("patY",&patternY);
	t->Branch("chi2Y",&chi2Y);
	t->Branch("chi2X",&chi2X);
	t->Branch("slY",&slY);
	t->Branch("slX",&slX);
	t->Branch("inX",&inX);
	t->Branch("inY",&inY);
	t->Branch("slerY",&slerY);
	t->Branch("slerX",&slerX);
	t->Branch("inerX",&inerX);
	t->Branch("inerY",&inerY);

	int flag0; 
	int flag1;
	int flag2;
	int flag3;
	int flag4;
	int flag5;
	int flag6;
	int flag7;
	int flag8;
	int flag9;
	int flag10;
	int flag11;
	t->Branch("flag0",&flag0);
	t->Branch("flag1",&flag1);
	t->Branch("flag2",&flag2);
	t->Branch("flag3",&flag3);
	t->Branch("flag4",&flag4);
	t->Branch("flag5",&flag5);
	t->Branch("flag6",&flag6);
	t->Branch("flag7",&flag7);
	t->Branch("flag8",&flag8);
	t->Branch("flag9",&flag9);
	t->Branch("flag10",&flag10);
	t->Branch("flag11",&flag11);

	t->Branch("fd",&fd,"fd[66]/D");


	//===================Prepare for fitting============================
	TVector3 vTrackU, vTrackD, vTrack;
	TVector3 vWireL, vWireR, vWire;
	TVector3 vDist;
	TVector3 vD;
	double zU, zD, dist;
	zU = Z1[0];
	zD = Z1[65];

	//===================Loop in Events============================
	Long64_t N = c->GetEntries();
	if (nEventMax&&nEventMax<N) N = nEventMax;
	std::cout<<"Processing "<<N<<" Events ..."<<std::endl;
	for (iEvent = 0;iEvent<N; iEvent++){
		if (iEvent%1000==0) std::cout<<(double)iEvent/N*100<<"..."<<std::endl;
		c->GetEntry(iEvent);
		//===================Get Hits============================
		getHits();
		//FIXME
		if (patternX<0||patternY<0) continue;
		if (patternX>3||patternY>3) continue;
		flag0 = flag[cht2chg(0)]; 
		flag1 = flag[cht2chg(1)];
		flag2 = flag[cht2chg(2)];
		flag3 = flag[cht2chg(3)];
		flag4 = flag[cht2chg(4)];
		flag5 = flag[cht2chg(5)];
		flag6 = flag[cht2chg(6)];
		flag7 = flag[cht2chg(7)];
		flag8 = flag[cht2chg(8)];
		flag9 = flag[cht2chg(9)];
		flag10 = flag[cht2chg(10)];
		flag11 = flag[cht2chg(11)];
		//===================Set Initial============================
		iniX = dd[indiceX[0]];
		if (patternX==0||patternX==3||patternX==9||patternX==11) iniX*=-1;

		if (patternX<=1||(patternX>=8&&patternX<=11)){
			sliX = (dd[indiceX[2]]-dd[indiceX[0]])/(Z1[indiceX[2]]-Z1[indiceX[0]]);
			if (patternX==0||patternX==9||patternX==11) sliX*=-1;
		}
		else{
			sliX = (dd[indiceX[2]]+dd[indiceX[0]])/(Z1[indiceX[2]]-Z1[indiceX[0]]);
			if (patternX==2) sliX*=-1;
		}

		if (patternX==4||patternX==5){
			sliX = (dd[indiceX[2]]-dd[indiceX[1]])/(Z1[indiceX[2]]-Z1[indiceX[1]]);
			if (patternX==5) sliX*=-1;
			iniX = dd[indiceX[0]];
			if (patternX==5) iniX*=-1;
		}
		else if (patternX==6||patternX==7){
			sliX = (dd[indiceX[2]]-dd[indiceX[1]])/(Z1[indiceX[2]]-Z1[indiceX[1]]);
			if (patternX==7) sliX*=-1;
			iniX = dd[indiceX[2]];
			if (patternX==7) iniX*=-1;
		}

		iniY = dd[indiceY[0]];
		if (patternY==0||patternY==3||patternY==9||patternY==11) iniY*=-1;

		if (patternY<=1||(patternY>=8&&patternY<=11)){
			sliY = (dd[indiceY[2]]-dd[indiceY[0]])/(Z1[indiceY[2]]-Z1[indiceY[0]]);
			if (patternY==0||patternY==9||patternY==11) sliY*=-1;
		}
		else{
			sliY = (dd[indiceY[2]]+dd[indiceY[0]])/(Z1[indiceY[2]]-Z1[indiceY[0]]);
			if (patternY==2) sliY*=-1;
		}

		if (patternY==4||patternY==5){
			sliY = (dd[indiceY[2]]-dd[indiceY[1]])/(Z1[indiceY[2]]-Z1[indiceY[1]]);
			if (patternY==5) sliY*=-1;
			iniY = dd[indiceY[0]];
			if (patternY==5) iniY*=-1;
		}
		else if (patternY==6||patternY==7){
			sliY = (dd[indiceY[2]]-dd[indiceY[1]])/(Z1[indiceY[2]]-Z1[indiceY[1]]);
			if (patternY==7) sliY*=-1;
			iniY = dd[indiceY[2]];
			if (patternY==7) iniY*=-1;
		}

		//===================Get The Tracker ============================
		fit_channel = 0;
		do_fit(sliX,iniX);
		gMinuit->GetParameter(0, slX, slerX);
		gMinuit->GetParameter(1, inX, inerX);
		getchi2(chi2X,slX,inX);
		getchi2(chi2iX,sliX,iniX);

		fit_channel = 1;
		do_fit(sliY,iniY);
		gMinuit->GetParameter(0, slY, slerY);
		gMinuit->GetParameter(1, inY, inerY);
		getchi2(chi2Y,slY,inY);
		getchi2(chi2iY,sliY,iniY);

		//===================Get Info for output============================
		vTrackU.SetXYZ(slY*zU+inY,slX*zU+inX,zU);
		vTrackD.SetXYZ(slY*zD+inY,slX*zD+inX,zD);
		vTrack=vTrackD-vTrackU;
		for (int ich = 0; ich<66; ich++){
			vWireL.SetXYZ(X1[ich],Y1[ich],Z1[ich]);
			vWireR.SetXYZ(X2[ich],Y2[ich],Z2[ich]);
			vWire=vWireR-vWireL;
			vDist = vWire.Cross(vTrack);
			vD=vWireL-vTrackU;
			dist = vD.Dot(vDist)/vDist.Mag();
			fd[ich] = dist;
			o_dist[ich] = dd[ich];
		}

		t->Fill();

		// Fill chi2
		hist_chi2X->Fill(chi2X);
		hist_chi2Y->Fill(chi2Y);
		// Fill sl
		hist_slX->Fill(slX);
		hist_slY->Fill(slY);
		// Get xt
		for ( int ch = 0; ch < 66; ch++){
			int cht = chg2cht(ch);
			int chp2 = chg2chp2(ch);
			int chp3 = chg2chp3(ch);
			if (cht>=0){
				xt_tr[cht]->Fill(fd[ch],dt[ch]);
			}
			else if (chp3>=0){
				xt_p3[chp3]->Fill(fd[ch],dt[ch]);
			}
			else if (chp2>=0){
				xt_p2[chp2]->Fill(fd[ch],dt[ch]);
			}
			int IDD = (int)((fd[ch]+U)/(U/DDSTEP));
			int IDD2 = (int)((dt[ch]+40)/((TMAX2+40)/DDSTEP2));
			if (fd[ch]>0) IDD2 += DDSTEP2;
			else IDD2 = DDSTEP2 - IDD2;
			// FIXME:
			if (chi2X<=chi2Xmax&&fabs(fd[ch])<=U){
				DTHISTS[IDD*66+ch]->Fill(dt[ch]);
			}
			if (chi2X<=chi2Xmax&&IDD2>=0&&IDD2<DDSTEP2*2){
				DTHISTS2[IDD2*66+ch]->Fill(-fabs(fd[ch]));
			}
		}
		// Get other plots
		for ( int i = 0; i<4; i++){
			for ( int j = 0; j<4; j++){
				int index = i+j*4;
				int d1 = j%2+1;
				int d2;
				if (i==0) d2 = 0;
				else if (i==1) d2 = 3;
				else if (i==2) d2 = 60;
				else if (i==3) d2 = 63;
				if (j<2)
					dt123[index]->Fill(dt[d2],dt[d2+d1]);
				else
					dt123[index]->Fill(dd[d2],dd[d2+d1]);
			}
		}
		for ( int i = 0; i<4; i++){
			double z = 0;
			if (i==1) z = 278;
			else if (i==2) z = 478;
			else if (i==3) z = 756;
			spot[i]->Fill(slY*z+inY,slX*z+inX);
		}
	}
	t->Write();
	f->Close();

	//===================Analyze time distribution============================
	std::cout<<"# Analyze time distribution"<<std::endl;
	std::vector<double> v_n;
	std::vector<double> v_center;
	std::vector<double> v_dist;
	std::vector<double> v_sigma;
	std::vector<double> v_chi2;
	std::vector<double> v_index;
	TCanvas * ca_temp = 0;
	TF1 * func_temp = 0;
	for ( int k  =0;  k< DDSTEP*2; k++){
		for (int ch = 0; ch<66; ch++){
			TH1D* h = DTHISTS[k*66+ch];
			double dist = (k+0.5)*U/DDSTEP-8;
			double t;
			int cht = chg2cht(ch);
			int chp2 = chg2chp2(ch);
			int chp3 = chg2chp3(ch);
			double sigmai;
			if (cht>=0){
				t = TMAX;
				sigmai = xt_sigma_tr/3.;
			}
			else{
				t = TMAX2;
				sigmai = xt_sigma_p/3.;
			}
			int centeribin = h->GetMaximumBin();
			double centeri = h->GetBinCenter(centeribin);
			double heighti = h->GetBinContent(centeribin);
			if (func_temp) delete func_temp;
			if (centeri<5*sigmai)
				func_temp = new TF1(Form("fDT_%d_%d",k,ch),"landau",-50,t);
			else
				func_temp = new TF1(Form("fDT_%d_%d",k,ch),"gaus",-50,t);
			func_temp->SetParameters(heighti,centeri,sigmai);
			h->Fit(Form("fDT_%d_%d",k,ch),"qN0","",centeri-sigmai*6,centeri+sigmai*6);
			double center = func_temp->GetParameter(1);
			double sigma = func_temp->GetParameter(2);
			h->SetTitle(Form("%lf mm: centeri = %lf, heighti = %lf; center = %lf, sigma = %lf",dist,centeri,heighti,center,sigma));
			double chi2=func_temp->GetChisquare();
			v_n.push_back(h->Integral());
			v_center.push_back(center);
			v_dist.push_back(dist);
			v_sigma.push_back(sigma);
			v_chi2.push_back(chi2);
			v_index.push_back(ch);
#ifdef PRINTR
			if (ch==19||ch==46){
				if (ca_temp) delete ca_temp;
				TCanvas * ca_temp = new TCanvas("ca_temp","ca_temp");
				gPad->SetGridx(1);
				gPad->SetGridy(1);
				h->SetMarkerStyle(20);
				h->SetMarkerSize(0.8);
				h->Draw("P");
				func_temp->Draw("SAME");
				ca_temp->SaveAs(Form("xtbin/%d.%d.new.png",ch,k));
			}
#endif
		}
	}

	std::vector<double> v_n2;
	std::vector<double> v_center2;
	std::vector<double> v_dist2;
	std::vector<double> v_sigma2;
	std::vector<double> v_chi22;
	std::vector<double> v_index2;
	for ( int k  =0;  k< DDSTEP2*2; k++){
		for (int ch = 0; ch<66; ch++){
			TH1D* h = DTHISTS2[k*66+ch];
			double t;
			double dist = fabs(k-DDSTEP2+0.5)*(TMAX2+40)/DDSTEP2-40;
			int cht = chg2cht(ch);
			int chp2 = chg2chp2(ch);
			int chp3 = chg2chp3(ch);
			double sigmai;
			if (cht>=0){
				t = TMAX;
				sigmai = xt_sigma_tr2/3.;
			}
			else{
				t = TMAX2;
				sigmai = xt_sigma_p2/3.;
			}
			int centeribin = h->GetMaximumBin();
			double centeri = h->GetBinCenter(centeribin);
			double heighti = h->GetBinContent(centeribin);
			if (func_temp) delete func_temp;
			if (fabs(fabs(centeri)-U)<1.6*sigmai)
				func_temp = new TF1(Form("fDT2_%d_%d",k,ch),"landau",-U-2,0);
			else
				func_temp = new TF1(Form("fDT2_%d_%d",k,ch),"gaus",-U-2,0);
			func_temp->SetParameters(heighti,centeri,sigmai);
			h->Fit(Form("fDT2_%d_%d",k,ch),"qN0","",centeri-sigmai*6,centeri+sigmai*6);
			double center = func_temp->GetParameter(1);
			double sigma = func_temp->GetParameter(2);
			h->SetTitle(Form("%lf ns: centeri = %lf, heighti = %lf; center = %lf, sigma = %lf",dist,centeri,heighti,center,sigma));
			double chi2=func_temp->GetChisquare();
			v_n2.push_back(h->Integral());
			v_center2.push_back(center);
			v_dist2.push_back(dist);
			v_sigma2.push_back(sigma);
			v_chi22.push_back(chi2);
			v_index2.push_back(ch);
#ifdef PRINTR
			if (ch==19||ch==46){
				if (ca_temp) delete ca_temp;
				TCanvas * ca_temp = new TCanvas("ca_temp","ca_temp");
				gPad->SetGridx(1);
				gPad->SetGridy(1);
				h->SetMarkerStyle(20);
				h->SetMarkerSize(0.8);
				h->Draw("P");
				func_temp->Draw("SAME");
				ca_temp->SaveAs(Form("xtbin/RV.%d.%d.new.png",ch,k));
			}
#endif
		}
	}

	//===================Get new XT============================
	std::cout<<"# Get new XT"<<std::endl;
	TF1 * f1[66];
	TF1 * f2[66];
	TF1 * f1o[66];
	TF1 * f2o[66];
	for ( int ch = 0; ch<66; ch++){
		v_xt_dr[ch].clear();
		v_xt_tr[ch].clear();
		v_xt_dl[ch].clear();
		v_xt_tl[ch].clear();
		//v_xt_tr[ch].push_back(0); v_xt_dr[ch].push_back(0);
		int cht = chg2cht(ch);
		int chp2 = chg2chp2(ch);
		int chp3 = chg2chp3(ch);
		double sigmai;
		if (cht>=0){
			sigmai = xt_sigma_tr/3.;
		}
		else{
			sigmai = xt_sigma_p/3.;
		}
		for ( int k  =0;  k< DDSTEP*2; k++){
			int index = k*66+ch;
			double dist = (k+0.5)*U/DDSTEP-8;
			double center = v_center[index];
			double sigma = v_sigma[index];
			double chi2 = v_chi2[index];
			int n = v_n[index];
			// FIXME
			if (n>N_MIN&&sigma<6*sigmai&&center>-10){
				if (dist>0){v_xt_tr[ch].push_back(center); v_xt_dr[ch].push_back(dist);}
				if (dist<0){v_xt_tl[ch].push_back(center); v_xt_dl[ch].push_back(dist);}
			}
		}
		if (cht>=0){
			sigmai = xt_sigma_tr2/3.;
		}
		else{
			sigmai = xt_sigma_p2/3.;
		}
		std::vector<double> vtempd;
		std::vector<double> vtempt;
		for ( int k  =0;  k< DDSTEP2*2; k++){
			int index = k*66+ch;
			double dist = fabs(k-DDSTEP2+0.5)*(TMAX2+40)/DDSTEP2-40;
			double center = fabs(v_center2[index]);
			double sigma = v_sigma2[index];
			double chi2 = v_chi22[index];
			int n = v_n2[index];
			// FIXME
			double tmax = 0;
			if (k<DDSTEP2){
				if (v_xt_tl[ch].size()>0){
					tmax = v_xt_tl[ch][0];
				}
			}
			else{
				if (v_xt_tr[ch].size()>0){
					tmax = v_xt_tr[ch][v_xt_tr[ch].size()-1];
				}
			}
			if (n>N_MIN&&sigma<3*sigmai&&center<10&&dist>fabs(tmax)){
				if (k<DDSTEP2){vtempt.push_back(dist); vtempd.push_back(-center);}
				else     {v_xt_tr[ch].push_back(dist); v_xt_dr[ch].push_back(center);}
			}
		}
		for (int ipoint = 0; ipoint<v_xt_dl[ch].size(); ipoint++){
			vtempd.push_back(v_xt_dl[ch][ipoint]);
			vtempt.push_back(v_xt_tl[ch][ipoint]);
		}
		v_xt_tl[ch] = vtempt;
		v_xt_dl[ch] = vtempd;
		//v_xt_tl[ch].push_back(0); v_xt_dl[ch].push_back(0);
		if (gxtl[ch]) delete gxtl[ch];
		if (gxtr[ch]) delete gxtr[ch];
		gxtl[ch] = new TGraph(v_xt_tl[ch].size(),&(v_xt_tl[ch][0]),&(v_xt_dl[ch][0]));
		gxtr[ch] = new TGraph(v_xt_tr[ch].size(),&(v_xt_tr[ch][0]),&(v_xt_dr[ch][0]));
		std::cout<<"ch = "<<ch<<std::endl;
		gxtl[ch]->Print();
		gxtr[ch]->Print();
		f1[ch] = new TF1(Form("f1_%d",ch),"pol5",0,TMAX2);
		f2[ch] = new TF1(Form("f2_%d",ch),"pol5",0,TMAX2);
		gxtl[ch]->Fit(Form("f1_%d",ch),"qN0","",0,TMAX2);
		gxtr[ch]->Fit(Form("f2_%d",ch),"qN0","",0,TMAX2);
		if (gxtl[ch]) delete gxtl[ch];
		if (gxtr[ch]) delete gxtr[ch];
		int size;
		double * pointer1,*pointer2;
		if (v_xt_tl[ch].size()>=8) {size = 8; pointer1 = &(v_xt_tl[ch][0]);pointer2 = &(v_xt_dl[ch][0]);}
		else {size = v_xt_tl[ch].size()-1; pointer1 = &(v_xt_tl[ch][0]);pointer2 = &(v_xt_dl[ch][0]);}
		gxtl[ch] = new TGraph(size,pointer1,pointer2);
		if (v_xt_tr[ch].size()>=8) {size = 8; pointer1 = &(v_xt_tr[ch][v_xt_tr[ch].size()-8]);pointer2 = &(v_xt_dr[ch][v_xt_tr[ch].size()-8]);}
		else {size = v_xt_tr[ch].size()-1; pointer1 = &(v_xt_tr[ch][0]);pointer2 = &(v_xt_dr[ch][0]);}
		gxtr[ch] = new TGraph(size,pointer1,pointer2);
		f1o[ch] = new TF1(Form("f1o_%d",ch),"pol1",0,TMAX2);
		f2o[ch] = new TF1(Form("f2o_%d",ch),"pol1",0,TMAX2);
		gxtl[ch]->Fit(Form("f1o_%d",ch),"qN0","",0,TMAX2);
		gxtr[ch]->Fit(Form("f2o_%d",ch),"qN0","",0,TMAX2);
		if (gxtl[ch]) delete gxtl[ch];
		if (gxtr[ch]) delete gxtr[ch];
	}
	double xt_tmin_l[66];
	double xt_tmax_l[66];
	double xt_tmax2_l[66];
	double xt_tmin_r[66];
	double xt_tmax_r[66];
	double xt_tmax2_r[66];
	double xt_dmin_l[66];
	double xt_dmax_l[66];
	double xt_dmin_r[66];
	double xt_dmax_r[66];
	double xt_sl_l[66];
	double xt_sl_r[66];
	TF1 * ftemp = 0;
	for ( int ch = 0; ch < 66; ch++){
		int tmax;
		int cht = chg2cht(ch);
		if (cht>=0){
			tmax = TMAX;
		}
		else{
			tmax = TMAX2;
		}

		double c0,c1,c2,c3,c4,c5;
		c0=f1[ch]->GetParameter(0);
		c1=f1[ch]->GetParameter(1);
		c2=f1[ch]->GetParameter(2);
		c3=f1[ch]->GetParameter(3);
		c4=f1[ch]->GetParameter(4);
		c5=f1[ch]->GetParameter(5);

		if (ftemp) delete ftemp;
		ftemp = new TF1("ftemp","pol4",0,tmax);
		ftemp->SetParameters(c1,
						     c2,
						     c3,
						     c4,
						     c5);
		double tmin_temp = ftemp->GetX(0,0,150);
		if (fabs(ftemp->Eval(tmin_temp))>1e-3){
			xt_tmin_l[ch] = 0;
		}
		else{
			xt_tmin_l[ch] = tmin_temp;
		}
		xt_dmin_l[ch] = f1[ch]->Eval(xt_tmin_l[ch]);

		int v_xt_tl_size = v_xt_tl[ch].size();
		if (v_xt_tl_size==0){
			xt_tmax_l[ch] = 0;
			xt_dmax_l[ch] = 0;
		}
		else{
			xt_tmax_l[ch] = v_xt_tl[ch][0];
			xt_dmax_l[ch] = f1[ch]->Eval(xt_tmax_l[ch]);
		}

		//xt_sl_l[ch] = 5*c5*pow(xt_tmax_l[ch],4)
		//	         +4*c4*pow(xt_tmax_l[ch],3)
		//	         +3*c3*pow(xt_tmax_l[ch],2)
		//	         +2*c2*pow(xt_tmax_l[ch],1)
		//	         +1*c1;
		xt_sl_l[ch] = f1o[ch]->GetParameter(1);
		if (xt_sl_l[ch]==0) xt_sl_l[ch] = -1;
		if (v_xt_tl_size<=1){
			xt_tmax2_l[ch] = xt_tmax_l[ch];
		}
		else {
			xt_tmax2_l[ch] = (-8-xt_dmax_l[ch])/xt_sl_l[ch]+xt_tmax_l[ch];
		}

		c0=f2[ch]->GetParameter(0);
		c1=f2[ch]->GetParameter(1);
		c2=f2[ch]->GetParameter(2);
		c3=f2[ch]->GetParameter(3);
		c4=f2[ch]->GetParameter(4);
		c5=f2[ch]->GetParameter(5);
		ftemp->SetParameters(c1,
						     c2,
						     c3,
						     c4,
						     c5);
		tmin_temp = ftemp->GetX(0,0,150);
		if (fabs(ftemp->Eval(tmin_temp))>1e-3){
			xt_tmin_r[ch] = 0;
		}
		else{
			xt_tmin_r[ch] = tmin_temp;
		}
		xt_dmin_r[ch] = f2[ch]->Eval(xt_tmin_r[ch]);

		int v_xt_tr_size = v_xt_tr[ch].size();
		if (v_xt_tr_size<=2){
			if (v_xt_tr_size<=1){
				xt_tmax_r[ch] = 0;
				xt_dmax_r[ch] = 1;
			}
			else{
				xt_tmax_r[ch] = v_xt_tr[ch][v_xt_tr_size-2];
				xt_dmax_r[ch] = f2[ch]->Eval(xt_tmax_r[ch]);
			}
			xt_tmax2_r[ch] = xt_tmax_r[ch];
		}
		else{
			xt_tmax_r[ch] = v_xt_tr[ch][v_xt_tr_size-2];
		}

		//xt_sl_r[ch] = 5*c5*pow(xt_tmax_r[ch],4)
		//	         +4*c4*pow(xt_tmax_r[ch],3)
		//	         +3*c3*pow(xt_tmax_r[ch],2)
		//	         +2*c2*pow(xt_tmax_r[ch],1)
		//	         +1*c1;
		xt_sl_r[ch] = f2o[ch]->GetParameter(1);
		if (xt_sl_r[ch]==0) xt_sl_r[ch] = 1;
		if (v_xt_tr_size>2){
			xt_dmax_r[ch] = f2[ch]->Eval(xt_tmax_r[ch]);
			xt_tmax2_r[ch] = (8-xt_dmax_r[ch])/xt_sl_r[ch]+xt_tmax_r[ch];
		}

		std::cout<<ch<<":["<<xt_tmax2_l[ch]<<","
			              <<xt_tmax_l[ch]<<","
			              <<xt_tmin_l[ch]<<"]["
			              <<xt_tmin_r[ch]<<","
			              <<xt_tmax_r[ch]<<","
						  <<xt_tmax2_r[ch]<<"]["
						  <<xt_dmax_l[ch]<<","
			              <<xt_dmin_l[ch]<<"]["
			              <<xt_dmin_r[ch]<<","
			              <<xt_dmax_r[ch]<<"]["
			              <<xt_sl_l[ch]<<","
			              <<xt_sl_r[ch]<<"]"<<std::endl;
		std::cout<<ch<<"  f1:"<<f1[ch]->GetParameter(0)<<","
						  <<f1[ch]->GetParameter(1)<<","
						  <<f1[ch]->GetParameter(2)<<","
						  <<f1[ch]->GetParameter(3)<<","
						  <<f1[ch]->GetParameter(4)<<","
						  <<f1[ch]->GetParameter(5)<<std::endl;
		std::cout<<ch<<"  f2:"<<f2[ch]->GetParameter(0)<<","
						  <<f2[ch]->GetParameter(1)<<","
						  <<f2[ch]->GetParameter(2)<<","
						  <<f2[ch]->GetParameter(3)<<","
						  <<f2[ch]->GetParameter(4)<<","
						  <<f2[ch]->GetParameter(5)<<std::endl;
	}

	//===================Draw XT for trackers============================
	std::cout<<"# Draw XT for trackers"<<std::endl;
	TCanvas * canvas;
	canvas = new TCanvas("c","c",1024,768);
	TPad * pad[100];
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	text->SetText(0.1,0.98,"XT Curve. ");;
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<4; i++){
		for (int j = 0; j<3; j++){
			int index = i*3+j;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./4*i,0.95/3*(2-j),1./4*(i+1),0.95/3*(3-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./4*i,1./3*(2-j),1./4*(i+1),1./3*(3-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	std::vector<double> v_temp_x;
	std::vector<double> v_temp_y;
	TGraph * g_temp = 0;
	for (int i = 0; i<3; i++){
		for (int j = 0; j<4; j++){
			int cht = j*3+i;
			int ch = cht;
			if (cht>=6) ch+=54;
			pad[cht]->cd();
			xt_tr[cht]->Draw("COLZ");
			ch_xt->Draw("t:d",Form("i==%d",ch),"SAMELP");
			v_temp_x.clear();
			v_temp_y.clear();
			for( int k = 0; k<=100; k++){
				double time = k/100.*(-xt_tmax2_l[ch]+xt_tmin_l[ch])+xt_tmax2_l[ch];
				v_temp_y.push_back(time);
				if (time>xt_tmax_l[ch]){
					v_temp_x.push_back((time-xt_tmax_l[ch])*xt_sl_l[ch]+xt_dmax_l[ch]);
				}
				else{
					v_temp_x.push_back(f1[ch]->Eval(time));
				}
			}
			for( int k = 0; k<=100; k++){
				double time = k/100.*(xt_tmax2_r[ch]-xt_tmin_r[ch])+xt_tmin_r[ch];
				v_temp_y.push_back(time);
				if (time>xt_tmax_r[ch]){
					v_temp_x.push_back((time-xt_tmax_r[ch])*xt_sl_r[ch]+xt_dmax_r[ch]);
				}
				else{
					v_temp_x.push_back(f2[ch]->Eval(time));
				}
			}
			g_temp = new TGraph(v_temp_x.size(),&(v_temp_x[0]),&(v_temp_y[0]));
			g_temp->SetLineColor(kRed);
			g_temp->SetLineWidth(0.5);
			g_temp->Draw("SAME");
		}
	}

	//===================Get ROOT File============================
	TChain * c2 = new TChain("t","t");
	buf.str(""); buf.clear();
	buf<<"../root/fit."<<runNo<<"."<<suffix<<iterationNo<<".root";
	c2->Add(buf.str().c_str());
	std::cout<<"Adding \""<<buf.str()<<"\""<<std::endl;
	c2->SetBranchAddress("iEvent",&iEvent);
	c2->SetBranchAddress("dt",dt);
	c2->SetBranchAddress("dd",dd);
	c2->SetBranchAddress("fd",fd);
	c2->SetBranchAddress("h",h);
	c2->SetBranchAddress("aa",aa);
	c2->SetBranchAddress("chi2X",&chi2X);
	c2->SetBranchAddress("chi2Y",&chi2Y);
	c2->SetBranchAddress("patX",&patternX);
	c2->SetBranchAddress("patY",&patternY);
	c2->SetBranchAddress("slX",&slX);
	c2->SetBranchAddress("slY",&slY);
	c2->SetBranchAddress("inX",&inX);
	c2->SetBranchAddress("inY",&inY);
	c2->SetBranchAddress("i",&ihit);

	//===================Set Output ROOT File============================
	std::cout<<"# Set Output ROOT File"<<std::endl;
	buf.str(""); buf.clear();
	buf<<"../root/ana."<<runNo<<"."<<suffix<<iterationNo<<".root";
	TFile * o_file2 = new TFile(buf.str().c_str(),"RECREATE");
	TTree * o_tree = new TTree("t","t");
	o_tree->Branch("iEvent",&iEvent);
	o_tree->Branch("dt",&dt,"dt[66]/D");
	o_tree->Branch("dd",&dd,"dd[66]/D");
	o_tree->Branch("fd",&fd,"fd[66]/D");
	o_tree->Branch("h",&h,"h[66]/I");
	o_tree->Branch("aa",&aa,"aa[66]/D");
	o_tree->Branch("i",&ihit,"ihit[66]/I");
	o_tree->Branch("chi2X",&chi2X);
	o_tree->Branch("chi2Y",&chi2Y);
	o_tree->Branch("patX",&patternX);
	o_tree->Branch("patY",&patternY);
	o_tree->Branch("slX",&slX);
	o_tree->Branch("slY",&slY);
	o_tree->Branch("inX",&inX);
	o_tree->Branch("inY",&inY);

	//===================Get reso and output ROOT file============================
	N=c2->GetEntries();
	std::cout<<"Processing "<<N<<" Events ..."<<std::endl;
	for (int iEntry = 0;iEntry<N; iEntry++){
		if (iEntry%1000==0) std::cout<<(double)iEntry/N*100<<"..."<<std::endl;
		c2->GetEntry(iEntry);
		// FIXME:
		if (chi2X>chi2Xmax) continue;
		// FIXME:
		getHits();
		double resop2,resop3;
		resop3 = 1e9;
		resop2 = 1e9;
		double fdp2,fdp3;
		fdp3 = 1e9;
		fdp2 = 1e9;
		for ( int ch = 0; ch < 66; ch++){
			dt[ch] = fmod(dt[ch]+1000,power2_15)-1000;
			int chp = -1;
			if(ch<=32){
				if (ch<=15&&ch>=13) chp = 15-ch;
				else if (ch<=20&&ch>=18) chp = 20-ch+3;
				else if (ch<=26&&ch>=24) chp = 26-ch+6;
			}
			else{
				if (ch<=42&&ch>=40) chp = 42-ch;
				else if (ch<=47&&ch>=45) chp = 47-ch+3;
				else if (ch<=53&&ch>=51) chp = 53-ch+6;
			}
			double dist;
			double tmin, tmax, tmax2;
			double dmin, dmax, dmax2;
			double slope;
			TF1* f;
			if (fd[ch]<0){
				tmin = xt_tmin_l[ch];
				tmax = xt_tmax_l[ch];
				tmax2 = xt_tmax2_l[ch];
				dmin = xt_dmin_l[ch];
				dmax = xt_dmax_l[ch];
				dmax2 = -8;
				slope = xt_sl_l[ch];
				f = f1[ch];
			}
			else{
				tmin = xt_tmin_r[ch];
				tmax = xt_tmax_r[ch];
				tmax2 = xt_tmax2_r[ch];
				dmin = xt_dmin_r[ch];
				dmax = xt_dmax_r[ch];
				dmax2 = 8;
				slope = xt_sl_r[ch];
				f = f2[ch];
			}

			if (dt[ch]>tmax2+10){
				dist=1e9;
			}
			else if (dt[ch]>tmax2){
				dist = dmax2;
			}
			else if (dt[ch]>tmax){
				dist = (dt[ch]-tmax)*slope+dmax;
			}
			else if (dt[ch]<tmin-10){
				dist=-1e9;
			}
			else if (dt[ch]<tmin){
				dist=(dt[ch]<0?0:dt[ch])/tmin*dmin;
			}
			else{
				if (fd[ch]<0){
					dist = f->Eval(dt[ch]);
				}
				else{
					dist = f->Eval(dt[ch]);
				}
				if (fabs(dist)>fabs(dmax)) dist = dmax;
				if (fabs(dist)<fabs(dmin)) dist = dmin;
			}

			dd[ch] = dist;

			double ll = sqrt(2)*16-2*fabs(fd[ch]);
			double gg;
			if (aa[ch]>735.346){
				gg = f_a2c->Eval(735)*1e-15/1.6e-19/npair/ll*10;
			}
			else{
				gg = f_a2c->Eval(aa[ch])*1e-15/1.6e-19/npair/ll*10;
			}
			hist_gg[ch]->Fill(fd[ch],log(gg)/log(10));

			int ir = (int)fabs(fd[ch])/1;
			double reso = fabs(fd[ch])-fabs(dist);
			if (ir<8){
				if (ch==19){
					resop3 = reso;
					fdp3 = fd[ch];
					// FIXME: no need to cut on Y
					//if (chi2X<=chi2Xmax&&chi2Y<=chi2Ymax){
					if (chi2X<=chi2Xmax){
						reso_p3[ir]->Fill(resop3);
						adcpeak_p3[ir]->Fill(h[ch]);
						adcsum_p3[ir]->Fill(aa[ch]);
						Neffi_p3[ir]+=1;
					}
					//printf("%lf-%lf=%lf\n",dist,fabs(fd[ch]),fabs(fd[ch])-dist);
					//printf("%d: p3[%d]: Mean = %lf, RMS = %lf\n",iEntry,ir,reso_p3[ir]->GetMean(),reso_p3[ir]->GetRMS());
				}
				else if (ch==46){
					resop2 = reso;
					fdp2 = fd[ch];
					// FIXME: no need to cut on Y
					//if (chi2X<=chi2Xmax&&chi2Y<=chi2Ymax){
					if (chi2X<=chi2Xmax){
						reso_p2[ir]->Fill(resop2);
						adcpeak_p2[ir]->Fill(h[ch]);
						adcsum_p2[ir]->Fill(aa[ch]);
						Neffi_p2[ir]+=1;
					}
				}
			}
		}
		o_tree->Fill();
#ifdef GETADC
		chain2->GetEntry(iEvent);
		int offset = 0;
		double minfd = 1e9;
		for ( int ich = 0; ich<66; ich++){
			if (fabs(fd[ich])<U&&fabs(fd[ich])<minfd){
				minfd = fabs(fd[ich]);
				offset = clockNumberDriftTime[ich][ihit[ich]];
			}
		}
		for ( int ich = 0; ich<66; ich++){
			int offseti = clockNumberDriftTime[ich][ihit[ich]];
			for ( int iadc = 0; iadc<NUM_OF_WAVE_SAMPLES; iadc++ ){
				if (ihit[ich]>=0){
					if (fabs(fd[ich])<=U){
						hist_wf0[ich]->Fill(iadc-offset,adc[ich][iadc]);
						int ir = (int)fabs(fd[ich])/1;
						if (ich==19){
							hist_wf_p3[ir]->Fill(iadc-offseti,adc[ich][iadc]);
						}
						else if (ich==46){
							hist_wf_p2[ir]->Fill(iadc-offseti,adc[ich][iadc]);
						}
					}
					else{
						hist_wf2[ich]->Fill(iadc-offseti,adc[ich][iadc]);
					}
					if (fabs(fd[ich])>minfd){
						hist_wf1[ich]->Fill(iadc-offseti,adc[ich][iadc]);
					}
				}
			}
		}
#endif
#ifdef EVENTDISPLAY
		// Draw Event Display
		TString displayName = "";
		// FIXME
		if (chi2X<7&&((fabs(resop2)>0.5&&fabs(fdp2)<U)||(fabs(resop3)>0.5&&fabs(fdp3)<U))){
			displayName = "proto";
		}
		else if (chi2X>7||chi2Y>7){
			displayName = "tracker";
		}
		else{
			displayName = "normal";
		}
		if (iEvent>100&&displayName!="proto") continue;
		if (iEvent>1000&&displayName=="proto") continue;
		// X tracker display
		pad2[0]->cd();
		buf.str("");buf.clear();
		buf<<"X Trackers: chi = "<<chi2X<<", in = "<<inX<<", sl = "<<slX;
		hdisX->SetTitle(buf.str().c_str());
		hdisX->GetXaxis()->SetTitle("Z [mm]");
		hdisX->GetYaxis()->SetTitle("Y [mm]");
		hdisX->Draw();
		lineX->SetX1((-20-inX)/slX);
		lineX->SetY1(-20);
		lineX->SetX2((20-inX)/slX);
		lineX->SetY2(20);
		lineX->Draw("SAME");
		for(int i = 0; i<6; i++){
			int ch;
			if (i<4) ch = indiceX[i];
			else if (i==4) ch = 19;
			else if (i==5) ch = 46;
			if(i<4){
				ellipse[i]->SetX1(Z1[ch]);ellipse[i]->SetY1(Y1[ch]);ellipse[i]->SetR1(dd[ch]);ellipse[i]->SetR2(dd[ch]);
			}
			else{
				ellipse[i]->SetX1((Z1[ch]+Z2[ch])/2);ellipse[i]->SetY1((Y1[ch]+Y2[ch])/2);ellipse[i]->SetR1(0.5);ellipse[i]->SetR2(0.1);
			}
			ellipse[i]->Draw("SAME");
			double dist;
			dist = dd[ch];
			buf.str("");buf.clear();
			buf<<"("<<dt[ch]<<","<<dist<<","<<fd[ch]<<")";
			if (i>=4){
				textEll[i]->SetText((Z1[ch]+Z2[ch])/2-50,(Y1[ch]+Y2[ch])/2-15,buf.str().c_str());
			}
			else if (i>=2)
				textEll[i]->SetText(Z1[ch]-200,Y1[ch],buf.str().c_str());
			else
				textEll[i]->SetText(Z1[ch]+50,Y1[ch],buf.str().c_str());
			textEll[i]->Draw("SAME");
		}
		// Y tracker display
		pad2[1]->cd();
		buf.str("");buf.clear();
		buf<<"Y Trackers: chi = "<<chi2Y<<", in = "<<inY<<", sl = "<<slY;
		hdisY->SetTitle(buf.str().c_str());
		hdisY->GetXaxis()->SetTitle("Z [mm]");
		hdisY->GetYaxis()->SetTitle("X [mm]");
		hdisY->Draw();
		lineY->SetX1((-20-inY)/slY);
		lineY->SetY1(-20);
		lineY->SetX2((20-inY)/slY);
		lineY->SetY2(20);
		lineY->Draw("SAME");
		for(int i = 0; i<6; i++){
			int ch;
			if (i<4) ch = indiceY[i];
			else if (i==4) ch = 19;
			else if (i==5) ch = 46;
			if (i<4){
				ellipse[i+6]->SetX1(Z1[ch]);ellipse[i+6]->SetY1(X1[ch]);ellipse[i+6]->SetR1(dd[ch]);ellipse[i+6]->SetR2(dd[ch]);
			}
			else{
				ellipse[i+6]->SetX1((Z1[ch]+Z2[ch])/2);ellipse[i+6]->SetY1((X1[ch]+X2[ch])/2);ellipse[i+6]->SetR1(0.5);ellipse[i+6]->SetR2(0.1);
			}
			ellipse[i+6]->Draw("SAME");
			double dist;
			dist = dd[ch];
			buf.str("");buf.clear();
			buf<<"("<<dt[ch]<<","<<dist<<","<<fd[ch]<<")";
			if (i>=4){
				textEll[i+6]->SetText((Z1[ch]+Z2[ch])/2-50,(X1[ch]+X2[ch])/2-15,buf.str().c_str());
			}
			else if (i>=2)
				textEll[i+6]->SetText(Z1[ch]-200,X1[ch],buf.str().c_str());
			else
				textEll[i+6]->SetText(Z1[ch]+50,X1[ch],buf.str().c_str());
			textEll[i+6]->Draw("SAME");
		}

		for (int i = 2; i<12; i++){
			int ch = -1;
			if (i>=2&&i<6) ch = indiceX[i-2];
			else if (i>=6&&i<10) ch = indiceY[i-6];
			else if (i==10) ch = 19;
			else if (i==11) ch = 46;
			if (ch<0) continue;
			pad2[i]->cd();
			gr_waveForm[i] = new TGraph(NUM_OF_WAVE_SAMPLES,vSample,adc[ch]);
			if (i<4)
				gr_waveForm[i]->GetHistogram()->SetTitle(Form("Waveform: Tracker X up, wireID=%d, event=%d",ch,iEvent));
			else if (i<6)
				gr_waveForm[i]->GetHistogram()->SetTitle(Form("Waveform: Tracker X down, wireID=%d, event=%d",ch,iEvent));
			else if (i<8)
				gr_waveForm[i]->GetHistogram()->SetTitle(Form("Waveform: Tracker Y up, wireID=%d, event=%d",ch,iEvent));
			else if (i<10)
				gr_waveForm[i]->GetHistogram()->SetTitle(Form("Waveform: Tracker Y down, wireID=%d, event=%d",ch,iEvent));
			else if (i==10)
				gr_waveForm[i]->GetHistogram()->SetTitle(Form("Waveform: Prototype 3, wireID=%d, event=%d",ch,iEvent));
			else if (i==11)
				gr_waveForm[i]->GetHistogram()->SetTitle(Form("Waveform: Prototype 2, wireID=%d, event=%d",ch,iEvent));
			gr_waveForm[i]->GetHistogram()->SetMinimum(MIN_WAVE_CH); 
			if (i>=10)
				gr_waveForm[i]->GetHistogram()->SetMaximum(MAX_WAVE_CH2);
			else
				gr_waveForm[i]->GetHistogram()->SetMaximum(MAX_WAVE_CH);
			gr_waveForm[i]->GetHistogram()->SetXTitle("Samplig ID");
			gr_waveForm[i]->GetHistogram()->SetYTitle("ADC(ch)");
			gr_waveForm[i]->SetMarkerStyle(24);
			gr_waveForm[i]->SetMarkerSize(0.8);
			gr_waveForm[i]->Draw("apwl");
			//text3->SetText(Form("q=%d",q[ch]));
			//text3->Draw();
			for (Int_t j=0; j<tdcNhit[ch]; j++) {
				textTDC[i][j]->SetText(clockNumberDriftTime[ch][j],MIN_WAVE_CH+10,Form("%d",(int)(tdc[ch][j])));
				textTDC[i][j]->Draw();
				markerTDC[i][j]->SetX(clockNumberDriftTime[ch][j]);
				markerTDC[i][j]->SetY(adc[ch][clockNumberDriftTime[ch][j]]);
				//printf("markerX = clockNumberDriftTime[%d][%d] = %d\n",ch,j,clockNumberDriftTime[ch][j]);
				//printf("markerY = adc[%d][%d] = %d\n",ch,clockNumberDriftTime[ch][j],adc[ch][clockNumberDriftTime[ch][j]]);
				// FIXME
				markerTDC[i][j]->SetMarkerColor(kRed);
				//if (adc[ch][j]>Hmin[ch]) markerTDC[i][j]->SetMarkerColor(kRed);
				//else markerTDC[i][j]->SetMarkerColor(kBlue);
				markerTDC[i][j]->Draw();
			}
			int max;
			if (i<10) max = MAX_WAVE_CH;
			else max = MAX_WAVE_CH2;
			double dist;
			dist = dd[ch];
			textWF[i]->SetText(5,max-20,Form("%lf-%lf=%lf",fd[ch],dist,fabs(fd[ch])-fabs(dist)));
			textWF[i]->Draw();
			double time_fd, time_dd;
			if (fd[ch]<0){
				time_fd = (f1[ch]->Eval(fd[ch])+1000)*32./1000;
				time_dd = (f1[ch]->Eval(dd[ch])+1000)*32./1000;
			}
			else{
				time_fd = (f2[ch]->Eval(fd[ch])+1000)*32./1000;
				time_dd = (f2[ch]->Eval(dd[ch])+1000)*32./1000;
			}
			line_fd[i]->SetX1(time_fd);
			line_fd[i]->SetY1(MIN_WAVE_CH);
			line_fd[i]->SetX2(time_fd);
			line_fd[i]->SetY2(max-50);
			line_fd[i]->Draw("SAME");
			line_dd[i]->SetX1(time_dd);
			line_dd[i]->SetY1(MIN_WAVE_CH);
			line_dd[i]->SetX2(time_dd);
			line_dd[i]->SetY2(max-50);
			line_dd[i]->Draw("SAME");
		}

		buf.str("");buf.clear();
		buf<<"EventDisplay/"<<displayName<<"/run"<<runNo<<"."<<iEvent<<"."<<suffix<<"png";
		canvas2->SaveAs(buf.str().c_str());
		buf.str("");buf.clear();
		buf<<"EventDisplay/"<<displayName<<"/run"<<runNo<<"."<<iEvent<<"."<<suffix<<"pdf";
		canvas2->SaveAs(buf.str().c_str());
#endif
	}
	o_tree->Write();
	for (int i = 0; i<8; i++){
		reso_p2[i]->Write();
		reso_p3[i]->Write();
		adcpeak_p3[i]->Write();
		adcpeak_p2[i]->Write();
		adcsum_p3[i]->Write();
		adcsum_p2[i]->Write();
	}
	for ( int i = 0; i<66; i++){
		hist_wf0[i]->Write();
		hist_wf1[i]->Write();
		hist_wf2[i]->Write();
		hist_gg[i]->Write();
	}
	for ( int i = 0; i<8; i++){
		hist_wf_p2[i]->Write();
		hist_wf_p3[i]->Write();
	}
	o_file2->Close();

	//===================output new info============================
	//prepare for xt and para
	buf.str("");
	buf.clear();
	buf<<"../info/info."<<runNo<<"."<<suffix<<iterationNo+1<<".root";
	TFile * of_info = new TFile(buf.str().c_str(),"RECREATE");
	TTree * otr_xt = new TTree("xt","xt");
	double o_t,o_d,o_s,o_c;
	int o_i,o_n;
	otr_xt->Branch("t",&o_t);
	otr_xt->Branch("d",&o_d);
	otr_xt->Branch("sig",&o_s);
	otr_xt->Branch("chi2",&o_c);
	otr_xt->Branch("i",&o_i);
	otr_xt->Branch("n",&o_n);

	// Output XT and para
	for ( int ch = 0; ch<66; ch++){
		int cht = -1;
		if (ch<=5||ch>=60){
			cht = ch;
			if (ch>=60) cht-=54;
		}
		for ( int k  =0;  k< DDSTEP*2; k++){
			int index = k*66+ch;
			o_t = v_center[index];
			o_d = v_dist[index];
			o_s = v_sigma[index];
			o_c = v_chi2[index];
			o_i = v_index[index];
			o_n = v_n[index];
			otr_xt->Fill();
		}
	}
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".xt.tracker."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".pdf";
	canvas->SaveAs(buf.str().c_str());
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".xt.tracker."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".png";
	canvas->SaveAs(buf.str().c_str());
	otr_xt->Write();
	of_info->Close();

	//===================Draw others============================
	// Draw resolution of Prototype 2
	double sigmap2[8];
	double effip2[8];
	if (canvas) delete canvas; canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.98,"Spatial Resolution (R_{fit}-R_{hit}) of Prototype 2. ");;
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<4; i++){
		for (int j = 0; j<2; j++){
			int index = j*4+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./4*i,0.95/2*(1-j),1./4*(i+1),0.95/2*(2-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./2*i,0.95/2*(1-j),1./2*(i+1),0.95/2*(2-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<4; i++){
		for (int j = 0; j<2; j++){
			int index = j*4+i;
			pad[index]->cd();
			TF1 * f  = new TF1(Form("f%d_%d",i,j),"gaus",-2,2);
			f->SetParameters(reso_p2[index]->GetBinContent(reso_p2[index]->GetMaximumBin()),0,0.2);
			int binc = reso_p2[index]->GetMaximumBin();
			reso_p2[index]->Fit(Form("f%d_%d",i,j),"qN0","",reso_p2[index]->GetBinCenter(binc-8),reso_p2[index]->GetBinCenter(binc+8));
			//printf("Fit(%lf,%lf)\n",reso_p2[index]->GetBinCenter(binc-8),reso_p2[index]->GetBinCenter(binc+8));
			double center = f->GetParameter(1);
			sigmap2[index] = f->GetParameter(2);
			if (fabs(center)>1) center = 0;
			int binl = reso_p2[index]->FindBin(center-0.6);
			int binr = reso_p2[index]->FindBin(center+0.6);
			if (binl<=0) binl = 1;
			effip2[index] = reso_p2[index]->Integral(binl,binr)/Neffi_p2[index];
			double effi2 = reso_p2[index]->Integral()/Neffi_p2[index];
			eff_p2_average += effip2[index];
			reso_p2[index]->SetTitle(Form("%d~%dmm: %3.1lf%% %3.0lfum",index,index+1,effip2[index]*100,sigmap2[index]*1000));
			reso_p2[index]->Draw();
			f->Draw("SAME");
			int max = reso_p2[index]->GetBinContent(reso_p2[index]->GetMaximumBin());
			text_reso1[index]->SetText(-1,max*0.9,Form("center = %lf, reso = %lf",center, sigmap2[index]));
			text_reso2[index]->SetText(-1,max*0.8,Form("eff_{0.6} = %lf/%lf = %lf",(double)reso_p2[index]->Integral(binl,binr),Neffi_p2[index],effip2[index]));
			text_reso3[index]->SetText(-1,max*0.7,Form("eff_{2} = %lf/%lf = %lf",(double)reso_p2[index]->Integral(),Neffi_p2[index],effi2));
			//text_reso1[index]->Draw("SAME");
			//text_reso2[index]->Draw("SAME");
			//text_reso3[index]->Draw("SAME");
			TLine * line_l = new TLine();
			line_l->SetLineColor(kRed);
			line_l->SetLineWidth(0.5);
			line_l->SetX1(center-0.6);
			line_l->SetY1(0);
			line_l->SetX2(center-0.6);
			line_l->SetY2(max);
			TLine * line_r = new TLine();
			line_r->SetLineColor(kRed);
			line_r->SetLineWidth(0.5);
			line_r->SetX1(center+0.6);
			line_r->SetY1(0);
			line_r->SetX2(center+0.6);
			line_r->SetY2(max);
			line_l->Draw("SAME");
			line_r->Draw("SAME");
		}
	}
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".reso_p2."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".pdf";
	canvas->SaveAs(buf.str().c_str());
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".reso_p2."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".png";
	canvas->SaveAs(buf.str().c_str());

	// Draw resolution of Prototype 3
	double sigmap3[8];
	double effip3[8];
	if (canvas) delete canvas; canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.98,"Spatial Resolution (R_{fit}-R_{hit}) of Prototype 3. ");;
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<4; i++){
		for (int j = 0; j<2; j++){
			int index = j*4+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./4*i,0.95/2*(1-j),1./4*(i+1),0.95/2*(2-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./2*i,0.95/2*(1-j),1./2*(i+1),0.95/2*(2-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<4; i++){
		for (int j = 0; j<2; j++){
			int index = j*4+i;
			pad[index]->cd();
			//reso_p3[index]->Print();
			//printf("[%d,%d], p3[%d]: Mean = %lf, RMS = %lf\n",i,j,index,reso_p3[index]->GetMean(),reso_p3[index]->GetRMS());
			TF1 * f  = new TF1(Form("f%d_%d",i,j),"gaus",-2,2);
			f->SetParameters(reso_p3[index]->GetBinContent(reso_p3[index]->GetMaximumBin()),0,0.2);
			int binc = reso_p3[index]->GetMaximumBin();
			reso_p3[index]->Fit(Form("f%d_%d",i,j),"qN0","",reso_p3[index]->GetBinCenter(binc-8),reso_p3[index]->GetBinCenter(binc+8));
			double center = f->GetParameter(1);
			sigmap3[index] = f->GetParameter(2);
			if (fabs(center)>1) center = 0;
			double binl = reso_p3[index]->FindBin(center-0.6);
			double binr = reso_p3[index]->FindBin(center+0.6);
			effip3[index] = reso_p3[index]->Integral(binl,binr)/Neffi_p3[index];
			double effi2 = reso_p3[index]->Integral()/Neffi_p3[index];
			eff_p3_average += effip3[index];
			// FIXME
			//printf("effi[%d] = Integral(%lf,%lf)/%d = %lf/%d = %lf\n",index,binl,binr,Neffi_p3[index],reso_p3[index]->Integral(binl,binr),Neffi_p3[index],effi);
			reso_p3[index]->SetTitle(Form("%d~%dmm: %3.1lf%% %3.0lfum",index,index+1,effip3[index]*100,sigmap3[index]*1000));
			reso_p3[index]->Draw();
			f->Draw("SAME");
			int max = reso_p3[index]->GetBinContent(reso_p3[index]->GetMaximumBin());
			text_reso1[index]->SetText(-1,max*0.9,Form("center = %lf, reso = %lf",center, sigmap3[index]));
			text_reso2[index]->SetText(-1,max*0.8,Form("eff_{0.6} = %lf/%lf = %lf",(double)reso_p3[index]->Integral(binl,binr),Neffi_p3[index],effip3[index]));
			text_reso3[index]->SetText(-1,max*0.7,Form("eff_{2} = %lf/%lf = %lf",(double)reso_p3[index]->Integral(),Neffi_p3[index],effi2));
			//text_reso1[index]->Draw("SAME");
			//text_reso2[index]->Draw("SAME");
			//text_reso3[index]->Draw("SAME");
			TLine * line_l = new TLine();
			line_l->SetLineColor(kRed);
			line_l->SetLineWidth(0.5);
			line_l->SetX1(center-0.6);
			line_l->SetY1(0);
			line_l->SetX2(center-0.6);
			line_l->SetY2(max);
			TLine * line_r = new TLine();
			line_r->SetLineColor(kRed);
			line_r->SetLineWidth(0.5);
			line_r->SetX1(center+0.6);
			line_r->SetY1(0);
			line_r->SetX2(center+0.6);
			line_r->SetY2(max);
			line_l->Draw("SAME");
			line_r->Draw("SAME");
		}
	}
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".reso_p3."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".pdf";
	canvas->SaveAs(buf.str().c_str());
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".reso_p3."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".png";
	canvas->SaveAs(buf.str().c_str());

	// Draw XT for p3
	if (canvas) delete canvas; canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.98,"XT Curve in Prototype 3. ");;
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<3; i++){
		for (int j = 0; j<3; j++){
			int index = j*3+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./3*i,0.95/3*(2-j),1./3*(i+1),0.95/3*(3-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./3*i,1./3*(2-j),1./3*(i+1),1./3*(3-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<3; i++){
		for (int j = 0; j<3; j++){
			int chp = j*3+i;
			int ch = chp;
			if (chp<3) ch = 15-chp;
			else if (chp<6) ch = 20-chp+3;
			else ch = 26-chp+6;
			pad[chp]->cd();
			xt_p3[chp]->Draw("COLZ");
			v_temp_x.clear();
			v_temp_y.clear();
			for( int k = 0; k<=100; k++){
				double time = k/100.*(-xt_tmax2_l[ch]+xt_tmin_l[ch])+xt_tmax2_l[ch];
				v_temp_y.push_back(time);
				if (time>xt_tmax_l[ch]){
					v_temp_x.push_back((time-xt_tmax_l[ch])*xt_sl_l[ch]+xt_dmax_l[ch]);
				}
				else{
					v_temp_x.push_back(f1[ch]->Eval(time));
				}
			}
			for( int k = 0; k<=100; k++){
				double time = k/100.*(xt_tmax2_r[ch]-xt_tmin_r[ch])+xt_tmin_r[ch];
				v_temp_y.push_back(time);
				if (time>xt_tmax_r[ch]){
					v_temp_x.push_back((time-xt_tmax_r[ch])*xt_sl_r[ch]+xt_dmax_r[ch]);
				}
				else{
					v_temp_x.push_back(f2[ch]->Eval(time));
				}
			}
			g_temp = new TGraph(v_temp_x.size(),&(v_temp_x[0]),&(v_temp_y[0]));
			g_temp->SetLineColor(kRed);
			g_temp->SetLineWidth(0.5);
			g_temp->Draw("SAME");
			g_temp = new TGraph(v_xt_tl[ch].size(),&(v_xt_dl[ch][0]),&(v_xt_tl[ch][0]));
			g_temp->Draw("PSAME");
			g_temp = new TGraph(v_xt_tr[ch].size(),&(v_xt_dr[ch][0]),&(v_xt_tr[ch][0]));
			g_temp->Draw("PSAME");
		}
	}
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".xt.p3."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".pdf";
	canvas->SaveAs(buf.str().c_str());
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".xt.p3."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".png";
	canvas->SaveAs(buf.str().c_str());

	// Draw XT for p2
	if (canvas) delete canvas; canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.98,"XT Curve in Prototype 2. ");;
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<3; i++){
		for (int j = 0; j<3; j++){
			int index = j*3+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./3*i,0.95/3*(2-j),1./3*(i+1),0.95/3*(3-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./3*i,1./3*(2-j),1./3*(i+1),1./3*(3-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<3; i++){
		for (int j = 0; j<3; j++){
			int chp = j*3+i;
			int ch = chp;
			if (chp<3) ch = 42-chp;
			else if (chp<6) ch = 47-chp+3;
			else ch = 53-chp+6;
			pad[chp]->cd();
			xt_p2[chp]->Draw("COLZ");
			v_temp_x.clear();
			v_temp_y.clear();
			for( int k = 0; k<=100; k++){
				double time = k/100.*(-xt_tmax2_l[ch]+xt_tmin_l[ch])+xt_tmax2_l[ch];
				v_temp_y.push_back(time);
				if (time>xt_tmax_l[ch]){
					v_temp_x.push_back((time-xt_tmax_l[ch])*xt_sl_l[ch]+xt_dmax_l[ch]);
				}
				else{
					v_temp_x.push_back(f1[ch]->Eval(time));
				}
			}
			for( int k = 0; k<=100; k++){
				double time = k/100.*(xt_tmax2_r[ch]-xt_tmin_r[ch])+xt_tmin_r[ch];
				v_temp_y.push_back(time);
				if (time>xt_tmax_r[ch]){
					v_temp_x.push_back((time-xt_tmax_r[ch])*xt_sl_r[ch]+xt_dmax_r[ch]);
				}
				else{
					v_temp_x.push_back(f2[ch]->Eval(time));
				}
			}
			g_temp = new TGraph(v_temp_x.size(),&(v_temp_x[0]),&(v_temp_y[0]));
			g_temp->SetLineColor(kRed);
			g_temp->SetLineWidth(0.5);
			g_temp->Draw("SAME");
			g_temp = new TGraph(v_xt_tl[ch].size(),&(v_xt_dl[ch][0]),&(v_xt_tl[ch][0]));
			g_temp->Draw("PSAME");
			g_temp = new TGraph(v_xt_tr[ch].size(),&(v_xt_dr[ch][0]),&(v_xt_tr[ch][0]));
			g_temp->Draw("PSAME");
		}
	}
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".xt.p2."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".pdf";
	canvas->SaveAs(buf.str().c_str());
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".xt.p2."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".png";
	canvas->SaveAs(buf.str().c_str());

	// Draw Wave Form for p2
	if (canvas) delete canvas; canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.98,"Wave form in Prototype 2. ");;
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<3; i++){
		for (int j = 0; j<3; j++){
			int index = j*3+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./3*i,0.95/3*(2-j),1./3*(i+1),0.95/3*(3-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./3*i,1./3*(2-j),1./3*(i+1),1./3*(3-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<3; i++){
		for (int j = 0; j<3; j++){
			int chp = j*3+i;
			int ch = chp;
			if (chp<3) ch = 42-chp;
			else if (chp<6) ch = 47-chp+3;
			else ch = 53-chp+6;
			pad[chp]->cd();
			pad[chp]->SetLogz(1);
			hist_wf1[ch]->Draw("COLZ");
		}
	}
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".wf.p2."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".pdf";
	canvas->SaveAs(buf.str().c_str());
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".wf.p2."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".png";
	canvas->SaveAs(buf.str().c_str());

	// Draw Wave Form for p3
	if (canvas) delete canvas; canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.98,"Wave form in Prototype 3. ");;
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<3; i++){
		for (int j = 0; j<3; j++){
			int index = j*3+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./3*i,0.95/3*(2-j),1./3*(i+1),0.95/3*(3-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./3*i,1./3*(2-j),1./3*(i+1),1./3*(3-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<3; i++){
		for (int j = 0; j<3; j++){
			int chp = j*3+i;
			int ch = chp;
			if (chp<3) ch = 15-chp;
			else if (chp<6) ch = 20-chp+3;
			else ch = 26-chp+6;
			pad[chp]->cd();
			pad[chp]->SetLogz(1);
			hist_wf1[ch]->Draw("COLZ");
		}
	}
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".wf.p3."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".pdf";
	canvas->SaveAs(buf.str().c_str());
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".wf.p3."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".png";
	canvas->SaveAs(buf.str().c_str());

	// Draw wave form (according to fd) of Prototype 2
	if (canvas) delete canvas; canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.98,"Wave From of Prototype 2 center wire. ");;
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<4; i++){
		for (int j = 0; j<2; j++){
			int index = j*4+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./4*i,0.95/2*(1-j),1./4*(i+1),0.95/2*(2-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./2*i,0.95/2*(1-j),1./2*(i+1),0.95/2*(2-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<4; i++){
		for (int j = 0; j<2; j++){
			int index = j*4+i;
			pad[index]->cd();
			pad[index]->SetLogz(1);
			hist_wf_p2[index]->SetTitle(Form("%d~%dmm: %3.1lf%% %3.0lfum",index,index+1,effip2[index]*100,sigmap2[index]*1000));
			hist_wf_p2[index]->Draw("COLZ");
		}
	}
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".wf_p2."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".pdf";
	canvas->SaveAs(buf.str().c_str());
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".wf_p2."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".png";
	canvas->SaveAs(buf.str().c_str());

	// Draw wave form (according to fd) of Prototype 3
	if (canvas) delete canvas; canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.98,"Wave From of Prototype 3 center wire. ");;
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<4; i++){
		for (int j = 0; j<2; j++){
			int index = j*4+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./4*i,0.95/2*(1-j),1./4*(i+1),0.95/2*(2-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./2*i,0.95/2*(1-j),1./2*(i+1),0.95/2*(2-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<4; i++){
		for (int j = 0; j<2; j++){
			int index = j*4+i;
			pad[index]->cd();
			pad[index]->SetLogz(1);
			hist_wf_p3[index]->SetTitle(Form("%d~%dmm: %3.1lf%% %3.0lfum",index,index+1,effip3[index]*100,sigmap3[index]*1000));
			hist_wf_p3[index]->Draw("COLZ");
		}
	}
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".wf_p3."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".pdf";
	canvas->SaveAs(buf.str().c_str());
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".wf_p3."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".png";
	canvas->SaveAs(buf.str().c_str());

	// Draw adc peak of Prototype 2
	if (canvas) delete canvas; canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.98,"ADC Peak Value of Prototype 2 center wire. ");;
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<4; i++){
		for (int j = 0; j<2; j++){
			int index = j*4+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./4*i,0.95/2*(1-j),1./4*(i+1),0.95/2*(2-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./2*i,0.95/2*(1-j),1./2*(i+1),0.95/2*(2-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<4; i++){
		for (int j = 0; j<2; j++){
			int index = j*4+i;
			pad[index]->cd();
			adcpeak_p2[index]->SetTitle(Form("%d~%dmm: %3.1lf%% %3.0lfum",index,index+1,effip2[index]*100,sigmap2[index]*1000));
			adcpeak_p2[index]->Draw();
		}
	}
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".adcpeak_p2."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".pdf";
	canvas->SaveAs(buf.str().c_str());
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".adcpeak_p2."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".png";
	canvas->SaveAs(buf.str().c_str());

	// Draw adc peak of Prototype 3
	if (canvas) delete canvas; canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.98,"ADC Peak Value of Prototype 3 center wire. ");;
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<4; i++){
		for (int j = 0; j<2; j++){
			int index = j*4+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./4*i,0.95/2*(1-j),1./4*(i+1),0.95/2*(2-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./2*i,0.95/2*(1-j),1./2*(i+1),0.95/2*(2-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<4; i++){
		for (int j = 0; j<2; j++){
			int index = j*4+i;
			pad[index]->cd();
			adcpeak_p3[index]->SetTitle(Form("%d~%dmm: %3.1lf%% %3.0lfum",index,index+1,effip3[index]*100,sigmap3[index]*1000));
			adcpeak_p3[index]->Draw();
		}
	}
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".adcpeak_p3."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".pdf";
	canvas->SaveAs(buf.str().c_str());
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".adcpeak_p3."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".png";
	canvas->SaveAs(buf.str().c_str());

	// Draw adc sum of Prototype 2
	if (canvas) delete canvas; canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.98,"ADC Sum of Prototype 2 center wire. ");;
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<4; i++){
		for (int j = 0; j<2; j++){
			int index = j*4+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./4*i,0.95/2*(1-j),1./4*(i+1),0.95/2*(2-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./2*i,0.95/2*(1-j),1./2*(i+1),0.95/2*(2-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<4; i++){
		for (int j = 0; j<2; j++){
			int index = j*4+i;
			pad[index]->cd();
			adcsum_p2[index]->SetTitle(Form("%d~%dmm: %3.1lf%% %3.0lfum",index,index+1,effip2[index]*100,sigmap2[index]*1000));
			adcsum_p2[index]->Draw();
		}
	}
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".adcsum_p2."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".pdf";
	canvas->SaveAs(buf.str().c_str());
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".adcsum_p2."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".png";
	canvas->SaveAs(buf.str().c_str());

	// Draw adc sum of Prototype 3
	if (canvas) delete canvas; canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.98,"ADC Sum of Prototype 3 center wire. ");;
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<4; i++){
		for (int j = 0; j<2; j++){
			int index = j*4+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./4*i,0.95/2*(1-j),1./4*(i+1),0.95/2*(2-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./2*i,0.95/2*(1-j),1./2*(i+1),0.95/2*(2-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<4; i++){
		for (int j = 0; j<2; j++){
			int index = j*4+i;
			pad[index]->cd();
			adcsum_p3[index]->SetTitle(Form("%d~%dmm: %3.1lf%% %3.0lfum",index,index+1,effip3[index]*100,sigmap3[index]*1000));
			adcsum_p3[index]->Draw();
		}
	}
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".adcsum_p3."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".pdf";
	canvas->SaveAs(buf.str().c_str());
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".adcsum_p3."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".png";
	canvas->SaveAs(buf.str().c_str());

	// Draw Gas Gain for p2
	if (canvas) delete canvas; canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.98,"Gas Gain in Prototype 2. ");;
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<3; i++){
		for (int j = 0; j<3; j++){
			int index = j*3+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./3*i,0.95/3*(2-j),1./3*(i+1),0.95/3*(3-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./3*i,1./3*(2-j),1./3*(i+1),1./3*(3-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<3; i++){
		for (int j = 0; j<3; j++){
			int chp = j*3+i;
			int ch = chp;
			if (chp<3) ch = 42-chp;
			else if (chp<6) ch = 47-chp+3;
			else ch = 53-chp+6;
			pad[chp]->cd();
			//pad[chp]->SetLogz(1);
			hist_gg[ch]->Draw("COLZ");
			axis->Draw();
		}
	}
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".gg.p2."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".pdf";
	canvas->SaveAs(buf.str().c_str());
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".gg.p2."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".png";
	canvas->SaveAs(buf.str().c_str());

	// Draw Gas Gain for p3
	if (canvas) delete canvas; canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.98,"Gas Gain in Prototype 3. ");;
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<3; i++){
		for (int j = 0; j<3; j++){
			int index = j*3+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./3*i,0.95/3*(2-j),1./3*(i+1),0.95/3*(3-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./3*i,1./3*(2-j),1./3*(i+1),1./3*(3-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<3; i++){
		for (int j = 0; j<3; j++){
			int chp = j*3+i;
			int ch = chp;
			if (chp<3) ch = 15-chp;
			else if (chp<6) ch = 20-chp+3;
			else ch = 26-chp+6;
			pad[chp]->cd();
			//pad[chp]->SetLogz(1);
			hist_gg[ch]->Draw("COLZ");
			axis->Draw();
		}
	}
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".gg.p3."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".pdf";
	canvas->SaveAs(buf.str().c_str());
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".gg.p3."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".png";
	canvas->SaveAs(buf.str().c_str());

	// Draw beam spot and tracking effi
	for(int i = 0; i<=128; i++){
		hist_effX->SetBinContent(i,(double)hist_chi2X->Integral(0,i)/N);
		hist_effY->SetBinContent(i,(double)hist_chi2Y->Integral(0,i)/N);
	}
	if (canvas) delete canvas; canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.98,"Beam Spot At the Center of Different Chambers. ");;
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<3; i++){
		for (int j = 0; j<2; j++){
			int index = j*3+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./3*i,0.95/2*(1-j),1./3*(i+1),0.95/2*(2-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<3; i++){
		for (int j = 0; j<2; j++){
			int index = j*3+i;
			pad[index]->cd();
			if (i<2){
				int index2 = j*2+i;
				spot[index2]->Draw("COLZ");
			}
			else if (j==0){
				hist_effX->SetLineColor(kRed);
				hist_effY->SetLineColor(kBlue);
				hist_effX->Draw();
				hist_effY->Draw("SAME");
			}
			else if (j==1){
				hist_slX->SetLineColor(kRed);
				hist_slY->SetLineColor(kBlue);
				hist_slX->Draw();
				hist_slY->Draw("SAME");
			}
		}
	}
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".spot."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".pdf";
	canvas->SaveAs(buf.str().c_str());
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".spot."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".png";
	canvas->SaveAs(buf.str().c_str());

	// Draw dt123
	if (canvas) delete canvas; canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.98,"dt123. ");;
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<4; i++){
		for (int j = 0; j<4; j++){
			int index = j*4+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./4*i,0.95/4*(3-j),1./4*(i+1),0.95/4*(4-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./4*i,0.95/4*(3-j),1./4*(i+1),0.95/4*(4-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<4; i++){
		for (int j = 0; j<4; j++){
			int index = j*4+i;
			pad[index]->cd();
			dt123[index]->Draw("COLZ");
		}
	}
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".dt123."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".pdf";
	canvas->SaveAs(buf.str().c_str());
	buf.str("");
	buf.clear();
	buf<<"run"<<runNo<<".dt123."<<suffix<<iterationNo<<".chi2_"<<chi2Xmax<<".png";
	canvas->SaveAs(buf.str().c_str());

	// FIXME
	std::cout<<"#EFF: "<<eff_p2_average/8.<<" "<<eff_p3_average/8.<<std::endl;

	return 0;
}

//______________________________________________________________________________
void getHits(){
	// ...
	int nHitYChs = 0;
	int nHitXChs = 0;
	for (int ch = 0; ch<66; ch++){
		int side = 1;
		if (ch == 2||ch==5||ch==62||ch==65) side = -1;
		dd[ch] = t2x(dt[ch],ch,side);
		if (ch<=5||ch>=60){// tracker hits
			if (fabs(dd[ch]-U/2)<=U/2){
				flag[ch] = true;
				if (ch<=2||ch>=63)
					nHitXChs++;
				else
					nHitYChs++;
			}
			else{
				flag[ch] = false;
			}
		}
		else{ // chamber hits; don't care
			flag[ch] = true;
		}
	}

	// X tracker
	if (nHitXChs==4){
		if (flag[cht2chg(0)]&&flag[cht2chg(9)]){
			indiceX[0] = cht2chg(0);
			indiceX[2] = cht2chg(9);
			if (flag[cht2chg(1)]&&!flag[cht2chg(2)]&&flag[cht2chg(10)]&&!flag[cht2chg(11)]){
				patternX = 1;
				indiceX[1] = cht2chg(1);
				indiceX[3] = cht2chg(10);
			}
			else if (!flag[cht2chg(1)]&&flag[cht2chg(2)]&&!flag[cht2chg(10)]&&flag[cht2chg(11)]){
				patternX = 0;
				indiceX[1] = cht2chg(2);
				indiceX[3] = cht2chg(11);
				dd[cht2chg(0)] = t2x(dt[cht2chg(0)],cht2chg(0),-1);
				dd[cht2chg(9)] = t2x(dt[cht2chg(9)],cht2chg(9),-1);
				if (fabs(dd[cht2chg(0)]-U/2)>U/2){patternX = -4;flag[cht2chg(0)]=false;}
				if (fabs(dd[cht2chg(9)]-U/2)>U/2){patternX = -4;flag[cht2chg(9)]=false;}
			}
			else if (flag[cht2chg(1)]&&!flag[cht2chg(2)]&&!flag[cht2chg(10)]&&flag[cht2chg(11)]){
				patternX = 2;
				indiceX[1]= cht2chg(1);
				indiceX[3] = cht2chg(11);
				dd[cht2chg(9)] = t2x(dt[cht2chg(9)],cht2chg(9),-1);
				if (fabs(dd[cht2chg(9)]-U/2)>U/2){patternX = -4;flag[cht2chg(9)]=false;}
			}
			else if (!flag[cht2chg(1)]&&flag[cht2chg(2)]&&flag[cht2chg(10)]&&!flag[cht2chg(11)]){
				patternX = 3;
				indiceX[1]= cht2chg(2);
				indiceX[3] = cht2chg(10);
				dd[cht2chg(0)] = t2x(dt[cht2chg(0)],cht2chg(0),-1);
				if (fabs(dd[cht2chg(0)]-U/2)>U/2){patternX = -4;flag[cht2chg(0)]=false;}
			}
			else{
				patternX = -4;
			}
		}
	}
	else if (nHitXChs>4){
		patternX = -5;
	}
	else if (nHitXChs==3){
		if (flag[cht2chg(0)]&&flag[cht2chg(1)]&&!flag[cht2chg(2)]&&!flag[cht2chg(9)]&&flag[cht2chg(10)]&&!flag[cht2chg(11)]){
			indiceX[0] = cht2chg(0);
			indiceX[1] = cht2chg(1);
			indiceX[2] = cht2chg(10);
			patternX = 4;
		}
		else if (flag[cht2chg(0)]&&!flag[cht2chg(1)]&&flag[cht2chg(2)]&&!flag[cht2chg(9)]&&!flag[cht2chg(10)]&&flag[cht2chg(11)]){
			indiceX[0] = cht2chg(0);
			indiceX[1] = cht2chg(2);
			indiceX[2] = cht2chg(11);
			patternX = 5;
		}
		else if (!flag[cht2chg(0)]&&flag[cht2chg(1)]&&!flag[cht2chg(2)]&&flag[cht2chg(9)]&&flag[cht2chg(10)]&&!flag[cht2chg(11)]){
			indiceX[0] = cht2chg(1);
			indiceX[1] = cht2chg(10);
			indiceX[2] = cht2chg(9);
			patternX = 6;
		}
		else if (!flag[cht2chg(0)]&&!flag[cht2chg(1)]&&flag[cht2chg(2)]&&flag[cht2chg(9)]&&!flag[cht2chg(10)]&&flag[cht2chg(11)]){
			indiceX[0] = cht2chg(2);
			indiceX[1] = cht2chg(11);
			indiceX[2] = cht2chg(9);
			patternX = 7;
		}
		else if (flag[cht2chg(0)]&&flag[cht2chg(1)]&&!flag[cht2chg(2)]&&flag[cht2chg(9)]&&!flag[cht2chg(10)]&&!flag[cht2chg(11)]){
			indiceX[0] = cht2chg(0);
			indiceX[1] = cht2chg(1);
			indiceX[2] = cht2chg(9);
			patternX = 8;
		}
		else if (flag[cht2chg(0)]&&!flag[cht2chg(1)]&&flag[cht2chg(2)]&&flag[cht2chg(9)]&&!flag[cht2chg(10)]&&!flag[cht2chg(11)]){
			indiceX[0] = cht2chg(0);
			indiceX[1] = cht2chg(2);
			indiceX[2] = cht2chg(9);
			patternX = 9;
		}
		else if (flag[cht2chg(0)]&&!flag[cht2chg(1)]&&!flag[cht2chg(2)]&&flag[cht2chg(9)]&&flag[cht2chg(10)]&&!flag[cht2chg(11)]){
			indiceX[0] = cht2chg(0);
			indiceX[1] = cht2chg(10);
			indiceX[2] = cht2chg(9);
			patternX = 10;
		}
		else if (flag[cht2chg(0)]&&!flag[cht2chg(1)]&&!flag[cht2chg(2)]&&flag[cht2chg(9)]&&!flag[cht2chg(10)]&&flag[cht2chg(11)]){
			indiceX[0] = cht2chg(0);
			indiceX[1] = cht2chg(11);
			indiceX[2] = cht2chg(9);
			patternX = 11;
		}
		else{
			patternX = -3;
		}
	}
	else {
			patternX = -2;
	}

	// Y tracker
	if (nHitYChs==4){
		if (flag[cht2chg(3)]&&flag[cht2chg(6)]){
			indiceY[0] = cht2chg(3);
			indiceY[2] = cht2chg(6);
			if (flag[cht2chg(4)]&&!flag[cht2chg(5)]&&flag[cht2chg(7)]&&!flag[cht2chg(8)]){
				patternY = 0;
				indiceY[1] = cht2chg(4);
				indiceY[3] = cht2chg(7);
			}
			else if (!flag[cht2chg(4)]&&flag[cht2chg(5)]&&!flag[cht2chg(7)]&&flag[cht2chg(8)]){
				patternY = 1;
				indiceY[1] = cht2chg(5);
				indiceY[3] = cht2chg(8);
				dd[cht2chg(3)] = t2x(dt[cht2chg(3)],cht2chg(3),-1);
				dd[cht2chg(6)] = t2x(dt[cht2chg(6)],cht2chg(6),-1);
				if (fabs(dd[cht2chg(3)]-U/2)>U/2){patternY = -4;flag[cht2chg(3)]=false;}
				if (fabs(dd[cht2chg(6)]-U/2)>U/2){patternY = -4;flag[cht2chg(6)]=false;}
			}
			else if (!flag[cht2chg(4)]&&flag[cht2chg(5)]&&flag[cht2chg(7)]&&!flag[cht2chg(8)]){
				patternY = 2;
				indiceY[1] = cht2chg(5);
				indiceY[3] = cht2chg(7);
				dd[cht2chg(3)] = t2x(dt[cht2chg(3)],cht2chg(3),-1);
				if (fabs(dd[cht2chg(3)]-U/2)>U/2){patternY = -4;flag[cht2chg(3)]=false;}
			}
			else if (flag[cht2chg(4)]&&!flag[cht2chg(5)]&&!flag[cht2chg(7)]&&flag[cht2chg(8)]){
				patternY = 3;
				indiceY[1] = cht2chg(4);
				indiceY[3] = cht2chg(8);
				dd[cht2chg(6)] = t2x(dt[cht2chg(6)],cht2chg(6),-1);
				if (fabs(dd[cht2chg(6)]-U/2)>U/2){patternY = -4;flag[cht2chg(6)]=false;}
			}
			else{
				patternY = -4;
			}
		}
	}
	else if (nHitYChs>4){
		patternY = -5;
	}
	else if (nHitYChs==3){
		if (flag[cht2chg(3)]&&flag[cht2chg(5)]&&!flag[cht2chg(4)]&&!flag[cht2chg(6)]&&flag[cht2chg(8)]&&!flag[cht2chg(7)]){
			indiceY[0] = cht2chg(3);
			indiceY[1] = cht2chg(5);
			indiceY[2] = cht2chg(8);
			patternY = 4;
		}
		else if (flag[cht2chg(3)]&&!flag[cht2chg(5)]&&flag[cht2chg(4)]&&!flag[cht2chg(6)]&&!flag[cht2chg(8)]&&flag[cht2chg(7)]){
			indiceY[0] = cht2chg(3);
			indiceY[1] = cht2chg(4);
			indiceY[2] = cht2chg(7);
			patternY = 5;
		}
		else if (!flag[cht2chg(3)]&&flag[cht2chg(5)]&&!flag[cht2chg(4)]&&flag[cht2chg(6)]&&flag[cht2chg(8)]&&!flag[cht2chg(7)]){
			indiceY[0] = cht2chg(5);
			indiceY[1] = cht2chg(8);
			indiceY[2] = cht2chg(6);
			patternY = 6;
		}
		else if (!flag[cht2chg(3)]&&!flag[cht2chg(5)]&&flag[cht2chg(4)]&&flag[cht2chg(6)]&&!flag[cht2chg(8)]&&flag[cht2chg(7)]){
			indiceY[0] = cht2chg(4);
			indiceY[1] = cht2chg(7);
			indiceY[2] = cht2chg(6);
			patternY = 7;
		}
		else if (flag[cht2chg(3)]&&flag[cht2chg(5)]&&!flag[cht2chg(4)]&&flag[cht2chg(6)]&&!flag[cht2chg(8)]&&!flag[cht2chg(7)]){
			indiceY[0] = cht2chg(3);
			indiceY[1] = cht2chg(5);
			indiceY[2] = cht2chg(6);
			patternY = 8;
		}
		else if (flag[cht2chg(3)]&&!flag[cht2chg(5)]&&flag[cht2chg(4)]&&flag[cht2chg(6)]&&!flag[cht2chg(8)]&&!flag[cht2chg(7)]){
			indiceY[0] = cht2chg(3);
			indiceY[1] = cht2chg(4);
			indiceY[2] = cht2chg(6);
			patternY = 9;
		}
		else if (flag[cht2chg(3)]&&!flag[cht2chg(5)]&&!flag[cht2chg(4)]&&flag[cht2chg(6)]&&flag[cht2chg(8)]&&!flag[cht2chg(7)]){
			indiceY[0] = cht2chg(3);
			indiceY[1] = cht2chg(8);
			indiceY[2] = cht2chg(6);
			patternY = 10;
		}
		else if (flag[cht2chg(3)]&&!flag[cht2chg(5)]&&!flag[cht2chg(4)]&&flag[cht2chg(6)]&&!flag[cht2chg(8)]&&flag[cht2chg(7)]){
			indiceY[0] = cht2chg(3);
			indiceY[1] = cht2chg(7);
			indiceY[2] = cht2chg(6);
			patternY = 11;
		}
		else{
			patternY = -3;
		}
	}
	else {
			patternY = -2;
	}

	//std::cout<<"X: "<<nHitXChs<<", "<<patternX<<", ("<<dd[0]<<", "
	//												 <<dd[1]<<", "
	//												 <<dd[2]<<", "
	//												 <<dd[63]<<", "
	//												 <<dd[64]<<", "
	//												 <<dd[65]<<")"
	//                                          <<", ("<<flag[0]<<", "
	//												 <<flag[1]<<", "
	//												 <<flag[2]<<", "
	//												 <<flag[63]<<", "
	//												 <<flag[64]<<", "
	//												 <<flag[65]<<")"
	//												 <<std::endl;
	//std::cout<<"Y: "<<nHitYChs<<", "<<patternY<<", ("<<dd[3]<<", "
	//												 <<dd[4]<<", "
	//												 <<dd[5]<<", "
	//												 <<dd[60]<<", "
	//												 <<dd[61]<<", "
	//												 <<dd[62]<<")"
	//                                          <<", ("<<flag[3]<<", "
	//												 <<flag[4]<<", "
	//												 <<flag[5]<<", "
	//												 <<flag[60]<<", "
	//												 <<flag[61]<<", "
	//												 <<flag[62]<<")"
	//												 <<std::endl;
}

//______________________________________________________________________________
Double_t get_dist(float x,float y,Double_t sl, Double_t in)
{
	Double_t A = sl*x+in-y;
	Double_t value = fabs(A/sqrt(1+sl*sl));
	return value;
}

//______________________________________________________________________________
void getchi2(Double_t &f, Double_t sl, Double_t in)
{
	Int_t index;

	//calculate chisquare
	Double_t chisq = 0;
	Double_t delta;
	Double_t dfit;
	int N = NHITS_TO_CHOOSE;
	if (fit_channel == 0){
		if (patternX<=3){
			N=4;
		}
		else{
			N=3;
		}
	}
	else if (fit_channel == 1){
		if (patternY<=3){
			N=4;
		}
		else{
			N=3;
		}
	}

	for (Int_t i=0;i<NHITS_TO_CHOOSE; i++) {
		if (fit_channel == 0){
			index = indiceX[i];
			dfit = get_dist(Z1[index],Y1[index],sl,in);
		}
		else if (fit_channel == 1){
			index = indiceY[i];
			dfit = get_dist(Z1[index],X1[index],sl,in);
		}
		delta  = (dd[index]-dfit)/errord[index];
		chisq += delta*delta;
	}
	f = chisq;
}

//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	getchi2(f,*par,*(par+1));
}

//______________________________________________________________________________
Double_t t2x(double t, int iwire, int side){
	double x = 0;
	int type = -1;

	if (iwire==0||iwire==3||iwire==60||iwire==63) type = 0;
	else if (iwire<=5||(iwire>=61&&iwire<=65)) type = 1;
	else type = 2;

	if (type<0){
		std::cout<<"iwire = "<<iwire<<"! This wire type is not supported yet!"<<std::endl;
		return -1;
	}

	TF1 *f;
	int cht = iwire;
	double tmax;
	if (iwire>=60) cht -= 54;
	if (type<=1){
		if (side==-1) {
			f = XTL_func[cht];
			tmax = f->GetX(-U);
		}
		else if (side==1){
			f = XTR_func[cht];
			tmax = f->GetX(U);
		}
	}

	if (t<-5){
		return -3;
	}
	else if (t<0){
		return 0;
	}
	else if (t>tmax+5){
		return -4;
	}
	else if (t>tmax){
		t = tmax;
	}

	if (type<=1){
		x = fabs(f->Eval(t));
		//printf("%d:%lf->%lf\n",iwire,t,x);
		return x;
	}
	else if (type==2){
		x=t*8./0.6;
		return x;
	}
}

//______________________________________________________________________________
void do_fit(Double_t sli, Double_t ini){
	if(gMinuit) delete gMinuit;
	gMinuit = new TMinuit(2);  //initialize TMinuit with a maximum of 5 params
	gMinuit->SetFCN(fcn);
	arglist[0] = 0;
	arglist[1] = 0;
	gMinuit->SetPrintLevel(-1); // no print
	gMinuit->mnexcm("SET NOW", arglist ,1,ierflg); // no warning 

	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

	// Set starting values and step sizes for parameters
	gMinuit->mnparm(0, "slope", sli, 0.00001, sli-3*fabs(sli),sli+3*fabs(sli),ierflg);
	gMinuit->mnparm(1, "intercept", ini, 0.0001, ini-3*fabs(ini),ini+fabs(ini),ierflg);

	// Now ready for minimization step
	arglist[0] = 500.0;
	arglist[1] = 1.0;
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

	// Print results
	gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	//printf("====Rrestul====\n");
	//gMinuit->mnprin(3,amin);
}

int chg2cht(int i){
	int cht = -1;
	if (i<6) cht=i;
	else if (i>=60) cht = i-54;
	return cht;
}

int chg2chp3(int i){
	int chp = -1;
	if (i<=15&&i>=13) chp = 15-i;
	else if (i<=20&&i>=18) chp = 20-i+3;
	else if (i<=26&&i>=24) chp = 26-i+6;
	return chp;
}

int chg2chp2(int i){
	int chp = -1;
	if (i<=42&&i>=40) chp = 42-i;
	else if (i<=47&&i>=45) chp = 47-i+3;
	else if (i<=53&&i>=51) chp = 53-i+6;
	return chp;
}

int cht2chg(int i){
	int chg = -1;
	if (i<6) chg=i;
	else if (i<12) chg = i+54;
	return chg;
}

int chp32chg(int i){
	int chg = -1;
	if (i<3)       chg = 15-i;
	else if (i<6)  chg = 20-i+3;
	else if (i<12) chg = 26-i+6;
	return chg;
}

int chp22chg(int i){
	int chg = -1;
	if (i<3)       chg = 42-i;
	else if (i<6)  chg = 47-i+3;
	else if (i<12) chg = 53-i+6;
	return chg;
}

int get_bid(int i){
	int bid = -1;
	if (i<12) bid = 0;
	else if (i<33) bid = 1;
	else if (i<60) bid = 2;
	else bid = 3;
	return bid;
}

int get_bid_core(int i){
	int bid = -1;
	if (i>=0&&i<6) bid = 0;
	else if (i==15||i==19||i==24) bid = 1;
	else if (i==42||i==46||i==51) bid = 2;
	else if (i>=60&&i<66) bid = 3;
	return bid;
}

void print_usage(char* prog_name)
{
	fprintf(stderr,"\t%s [runNo] <[suffix] [nEventMax] [iterationNo] [N_Min]\n",prog_name);
}
