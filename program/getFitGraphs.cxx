#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <math.h>
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TH1D.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TText.h"
#include "TMarker.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TF1.h"
#include "TLatex.h"

#define NUM_OF_WAVE_SAMPLES 32
#define MIN_WAVE_CH 180
#define MAX_WAVE_CH 400
#define MAX_WAVE_CH2 650
#define HMIN 150
#define HMAX 330

int power2_15 = pow(2,15);

//#define EVENTDISPLAY

Double_t Z1[66],Y1[66],X1[66],errord[66];// left hand side
Double_t Z2[66],Y2[66],X2[66]; // right hand side

//______________________________________________________________________________
void print_usage(char* prog_name)
{
	fprintf(stderr,"\t%s [startNo] <[nFiles] [suffix] [nEventMax]>\n",prog_name);
}

int main(int argc, char** argv){
	if (argc<2){
		print_usage(argv[0]);
		return 1;
	}
	int startNo = (int)strtol(argv[1],NULL,10);
	int nFiles = 1;
	if (argc>=3) nFiles = (int)strtol(argv[2],NULL,10);
	TString suffix = "";
	if (argc>=4){
		suffix  = argv[3];
		suffix=suffix+".";
	}
	int iterationNo = 0;
	if (argc>=5){
		iterationNo = strtol(argv[4],NULL,10);
	}
	int nEventMax = 0;
	if (argc>=6) nEventMax = (int)strtol(argv[5],NULL,10);
	double chi2Xmax,chi2Ymax;
	chi2Xmax=1e9;
	if (argc>=7) chi2Xmax = (int)strtol(argv[6],NULL,10);
	chi2Ymax=chi2Xmax;
	std::stringstream buf;

	//===================Prepare Histograms============================
	// Hmin
	TH1D * hist_hmin_p2[9];
	TH1D * hist_hmin_p3[9];
	for (int i = 0; i<9; i++){
		hist_hmin_p2[i] = new TH1D(Form("hhmin_%d_p2",i),"hminp2",HMAX-HMIN,HMIN,HMAX);
		hist_hmin_p3[i] = new TH1D(Form("hhmin_%d_p3",i),"hminp3",HMAX-HMIN,HMIN,HMAX);
	}

	//chi2
	TH1D * hist_chi2X = new TH1D("chi2X","chi2X",128,0,20);
	TH1D * hist_chi2Y = new TH1D("chi2Y","chi2Y",128,0,20);
	TH1D * hist_effX = new TH1D("effX","effX",128,0,20);
	TH1D * hist_effY = new TH1D("effY","effY",128,0,20);

	TLatex * text = new TLatex();
	text->SetTextSize(0.02);
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
	int Neffi_p2[8]={0};
	int Neffi_p3[8]={0};

	int TMAX = 280; // ns
	int TMAX2 = 800; // ns
	double U =8; // mm
	int DDSTEP = 32; // steps
	int ICH = 1;
	std::vector<TH1D*> DTHISTS;
	for ( int i  =0;  i< DDSTEP*2; i++){
		for ( int j  =0;  j< 66; j++){
			name  = Form("dt_%d_%d",i,j);
			title = name;
			double t;
			if (j<=5||j>=60) t = TMAX;
			else t = TMAX2;
			TH1D * dthist = new TH1D(name,title,(50+t)/3,-50,t);
			DTHISTS.push_back(dthist);
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
		xt_tr[i] = new TH2D(name,title,256,-9,9,310,-10,300);
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
		xt_p3[i] = new TH2D(name,title,256,-16,16,850,-50,800);
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
		xt_p2[i] = new TH2D(name,title,256,-16,16,850,-50,800);
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
		spot[i] = new TH2D(name,title,256,-10,10,256,-10,10);
		spot[i]->GetXaxis()->SetTitle("X [mm]");
		spot[i]->GetYaxis()->SetTitle("Y [mm]");
	}
	for (int i = 0; i<16; i++){
		name = Form("dt123_%d",i);
		if (i<8)
			dt123[i] = new TH2D(name,name,280,0,280,280,0,280);
		else
			dt123[i] = new TH2D(name,name,256,0,8,256,0,8);
	}
	for (int i = 0; i<8; i++){
		name = Form("reso_p2_%d",i);
		reso_p2[i] = new TH1D(name,name,128,-2,2);
		name = Form("reso_p3_%d",i);
		reso_p3[i] = new TH1D(name,name,128,-2,2);
	}

#ifdef EVENTDISPLAY
	// Event Display
	int indiceX[4];
	int indiceY[4];
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

	//===================Get para for previous run============================
	// Get previous xt
	TChain * ch_xt = new TChain("xt","xt");
	ch_xt ->Add(Form("../info/info.%d."+suffix+"%d.root",startNo,iterationNo));
	//double i_t,i_d,i_s;
	//int i_i;
	//ch_xt->SetBranchAddress("t",&i_t);
	//ch_xt->SetBranchAddress("d",&i_d);
	//ch_xt->SetBranchAddress("sig",&i_s);
	//ch_xt->SetBranchAddress("i",&i_i);
	ch_xt->SetMarkerStyle(20);
	ch_xt->SetMarkerSize(0.15);
	ch_xt->SetMarkerColor(kBlue);
	ch_xt->SetLineColor(kBlue);
	ch_xt->SetLineWidth(0.5);

	// get previous tuning
	TChain * ch_tune = new TChain("para","para");
	ch_tune ->Add(Form("../info/info.%d."+suffix+"%d.root",startNo,iterationNo));
	double t01, t02, t03, t04;
	double dxu, dyu, dxd, dyd;
	ch_tune->SetBranchAddress("t01",&t01);
	ch_tune->SetBranchAddress("t02",&t02);
	ch_tune->SetBranchAddress("t03",&t03);
	ch_tune->SetBranchAddress("t04",&t04);
	ch_tune->SetBranchAddress("dxu",&dxu);
	ch_tune->SetBranchAddress("dyu",&dyu);
	ch_tune->SetBranchAddress("dxd",&dxd);
	ch_tune->SetBranchAddress("dyd",&dyd);
	ch_tune->GetEntry(0);

	// get Hmin
	double Hmin[66];
	TChain * ch_hmin = new TChain("hmin","hmin");
	ch_hmin->Add(Form("../info/info.%d."+suffix+"%d.root",startNo,iterationNo));
	int idhmin, hhmin;
	ch_hmin->SetBranchAddress("i",&idhmin);
	ch_hmin->SetBranchAddress("h",&hhmin);
	for(int i = 0; i<ch_hmin->GetEntries(); i++){
		ch_hmin->GetEntry(i);
		if (idhmin<0||idhmin>=66){
			printf("WARNING: Entry[%d], wireID = %d\n",i,idhmin);
		}
		Hmin[idhmin] = hhmin;
	}

	// get geometry
	// The z values
	// get the inital values
	TChain * ch_geom = new TChain("t","t");
	ch_geom->Add("../info/wire-position.root");
	double x1_geom, y1_geom, z1_geom, x2_geom, y2_geom, z2_geom;
	int id_geom;
	ch_geom->SetBranchAddress("wire_id",&id_geom);
	ch_geom->SetBranchAddress("x1",&x1_geom);
	ch_geom->SetBranchAddress("y1",&y1_geom);
	ch_geom->SetBranchAddress("z1",&z1_geom);
	ch_geom->SetBranchAddress("x2",&x2_geom);
	ch_geom->SetBranchAddress("y2",&y2_geom);
	ch_geom->SetBranchAddress("z2",&z2_geom);
	for(int i = 0; i<ch_geom->GetEntries(); i++){
		ch_geom->GetEntry(i);
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
		//if (id_geom>5&&id_geom<60){
		//	Y1[id_geom] *= -1;
		//	Y2[id_geom] *= -1;
		//	if (id_geom<33){
		//		Z1[id_geom]+=2*(278-Z1[id_geom]);
		//		Z2[id_geom]+=2*(278-Z2[id_geom]);
		//	}
		//	else{
		//		Y1[id_geom] -= 0.5;
		//		Y2[id_geom] -= 0.5;
		//		Z1[id_geom]+=2*(478-Z1[id_geom]);
		//		Z2[id_geom]+=2*(478-Z2[id_geom]);
		//	}
		//}
		double x0 = 0;
		if (id_geom==1||id_geom==2) x0 = dxu;
		else if (id_geom==4||id_geom==5) x0 = dyu;
		else if (id_geom==61||id_geom==62) x0 = dxd;
		else if (id_geom==64||id_geom==65) x0 = dyd;
		Y1[id_geom] -= x0;
		Y2[id_geom] -= x0;
	}
	// The errors on z values
	Float_t error = 0.2;
	for ( int i = 0; i<66; i++){
		errord[i]=error;
	}

	//===================Get ROOT Files============================
	TChain * c = new TChain("t","t");
	TString filename;
	for ( int iFile = startNo; iFile<startNo+nFiles; iFile++){
		filename = Form("../root/fit_%d_%d."+suffix+"%d.root",iFile,nFiles,iterationNo);
		c->Add(filename);
	}
	double slX,slY,inX,inY;
	double dt[66];
	double dd[66];
	double fd[66];
	double chi2X;
	double chi2Y;
	int h[66];
	double a[66];
	int patX;
	int patY;
	c->SetBranchAddress("dt",dt);
	c->SetBranchAddress("dd",dd);
	c->SetBranchAddress("fd",fd);
	c->SetBranchAddress("h",h);
	c->SetBranchAddress("a",a);
	c->SetBranchAddress("chi2X",&chi2X);
	c->SetBranchAddress("chi2Y",&chi2Y);
	c->SetBranchAddress("patX",&patX);
	c->SetBranchAddress("patY",&patY);
	c->SetBranchAddress("slX",&slX);
	c->SetBranchAddress("inX",&inX);
	c->SetBranchAddress("slY",&slY);
	c->SetBranchAddress("inY",&inY);
#ifdef EVENTDISPLAY
	int adc[66][32];
	int tdc[66][32];
	int clockNumberDriftTime[66][32];
	int tdcNhit[66];
	c->SetBranchAddress("adc",adc);
	c->SetBranchAddress("tdc",tdc);
	c->SetBranchAddress("clockNumberDriftTime",clockNumberDriftTime);
	c->SetBranchAddress("n",tdcNhit);
#endif

	//===================Loop in Events============================
	Long64_t N = c->GetEntries();
	if (nEventMax&&nEventMax<N) N = nEventMax;
	std::cout<<"Processing "<<N<<" Events ..."<<std::endl;
	for (int iEvent = 0;iEvent<N; iEvent++){
		if (iEvent%1000==0) std::cout<<(double)iEvent/N*100<<"..."<<std::endl;
		c->GetEntry(iEvent);
		double t0;
		hist_chi2X->Fill(chi2X);
		hist_chi2Y->Fill(chi2Y);
		// Modify time according to t0
		for ( int ch = 0; ch < 66; ch++){
			dt[ch] = fmod(dt[ch]+1000,power2_15)-1000;
			if (ch<=5) t0 = t01;
			else if (ch<=32) t0 = t02;
			else if (ch<=59) t0 = t03;
			else t0 = t04;
			dt[ch] -= t0;
		}
		// Get xt
		for ( int ch = 0; ch < 66; ch++){
			int cht = -1;
			int chp = -1;
			if (ch<=5||ch>=60){
				cht = ch;
				if (ch>=60) cht-=54;
			}
			else if(ch<=32){
				if (ch<=15&&ch>=13) chp = 15-ch;
				else if (ch<=20&&ch>=18) chp = 20-ch+3;
				else if (ch<=26&&ch>=24) chp = 26-ch+6;
			}
			else{
				if (ch<=42&&ch>=40) chp = 42-ch;
				else if (ch<=47&&ch>=45) chp = 47-ch+3;
				else if (ch<=53&&ch>=51) chp = 53-ch+6;
			}
			if (cht>=0){
				xt_tr[cht]->Fill(fd[ch],dt[ch]);
			}
			else if (chp>=0){
				if(ch<=32){
					xt_p3[chp]->Fill(fd[ch],dt[ch]);
				}
				else{
					xt_p2[chp]->Fill(fd[ch],dt[ch]);
				}
			}
			int IDD = (int)((fd[ch]+U)/(U/DDSTEP));
			// FIXME:
			if (chi2X<=chi2Xmax&&fabs(fd[ch])<=U)
				DTHISTS[IDD*66+ch]->Fill(dt[ch]);
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

	//===================Analyze time distribution============================
	std::cout<<"# Analyze time distribution"<<std::endl;
	std::vector<double> v_n;
	std::vector<double> v_center;
	std::vector<double> v_dist;
	std::vector<double> v_sigma;
	std::vector<double> v_index;
	for ( int k  =0;  k< DDSTEP*2; k++){
		for (int ch = 0; ch<66; ch++){
			TH1D* h = DTHISTS[k*66+ch];
			double dist = (k+0.5)*U/DDSTEP-8;
			double t;
			double sigmai;
			if (ch<=5||ch>=60){
				t = TMAX;
				sigmai = 2.5;
			}
			else{
				t = TMAX2;
				sigmai = 10;
			}
			int centeribin = h->GetMaximumBin();
			double centeri = h->GetBinCenter(centeribin);
			double heighti = h->GetBinContent(centeribin);
			TF1 * f;
			if (centeri<1.5*sigmai)
				f = new TF1(Form("fDT_%d_%d",k,ch),"landau",-50,t);
			else
				f = new TF1(Form("fDT_%d_%d",k,ch),"gaus",-50,t);
			f->SetParameters(heighti,centeri,sigmai);
			h->Fit(Form("fDT_%d_%d",k,ch),"N0","",centeri-sigmai*2,centeri+sigmai*2);
			double center = f->GetParameter(1);
			double sigma = f->GetParameter(2);
			h->SetTitle(Form("%lf mm: centeri = %lf, heighti = %lf; center = %lf, sigma = %lf",dist,centeri,heighti,center,sigma));
			v_n.push_back(h->GetEntries());
			v_center.push_back(center);
			v_dist.push_back(dist);
			v_sigma.push_back(sigma);
			v_index.push_back(ch);
			//FIXME
			//if (ch==0){
			//if (1){
			//	h->Draw();
			//	f->Draw("SAME");
			//	canvas->SaveAs(Form("DTHISTS/%d_%d.png",ch,k));
			//}
			////FIXME
			//printf("%d,%d,%d:v_center.push_back(%lf); v_dist.push_back(%lf)\n",ch,k,v_center.size(),center,dist);
		}
	}

	//===================Get new XT============================
	std::cout<<"# Get new XT"<<std::endl;
	TF1 * f1[66];
	TF1 * f2[66];
	TGraph * gxtl[66];
	TGraph * gxtr[66];
	std::vector<double> vddr;
	std::vector<double> vdtr;
	std::vector<double> vddl;
	std::vector<double> vdtl;
	for ( int ch = 0; ch<66; ch++){
		vddr.clear();
		vdtr.clear();
		vddl.clear();
		vdtl.clear();
		double sigmai;
		if (ch<=5||ch>=60){
			sigmai = 7.5;
		}
		else{
			sigmai = 25;
		}
		for ( int k  =0;  k< DDSTEP*2; k++){
			int index = k*66+ch;
			double dist = (k+0.5)*U/DDSTEP-8;
			double center = v_center[index];
			double sigma = v_sigma[index];
			int n = v_n[index];
			// FIXME
			if (n>100&&sigma<sigmai){
			//if (dist<6){
				if (dist>=0){vdtr.push_back(center); vddr.push_back(dist);}
				if (dist<=0){vdtl.push_back(center); vddl.push_back(dist);}
			}
		}
		gxtl[ch] = new TGraph(vdtl.size(),&(vddl[0]),&(vdtl[0]));
		gxtr[ch] = new TGraph(vdtr.size(),&(vddr[0]),&(vdtr[0]));
		f1[ch] = new TF1(Form("f1_%d",ch),"pol5",-8,0);
		f2[ch] = new TF1(Form("f2_%d",ch),"pol5",0,8);
		gxtl[ch]->Fit(Form("f1_%d",ch),"N0","",-8,0);
		gxtr[ch]->Fit(Form("f2_%d",ch),"N0","",0,8);
	}

	//===================Draw XT for trackers and get t1 x1============================
	std::cout<<"# Draw XT for trackers and get t1 x1"<<std::endl;
	TCanvas * canvas;
	canvas = new TCanvas("c","c",1024,768);
	TPad * pad[100];
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	double t11, t12, t13, t14;
	double dxu1, dyu1, dxd1, dyd1;
	t11  = t01;
	t12  = t02;
	t13  = t03;
	t14  = t04;
	dxu1 = 0;
	dyu1 = 0;
	dxd1 = 0;
	dyd1 = 0;
	text->SetText(0.1,0.975,"XT Curve ");
	text->Draw("SAME");
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
	for (int i = 0; i<3; i++){
		for (int j = 0; j<4; j++){
			int cht = j*3+i;
			int ch = cht;
			if (cht>=6) ch+=54;
			pad[cht]->cd();
			xt_tr[cht]->Draw("COLZ");
			gxtr[ch]->SetMarkerStyle(20);
			gxtr[ch]->SetMarkerSize(0.18);
			gxtr[ch]->SetMarkerColor(kRed);
			gxtr[ch]->SetLineColor(kRed);
			gxtr[ch]->SetLineWidth(0.5);
			gxtr[ch]->Draw("PLSAME");
			gxtl[ch]->SetMarkerStyle(20);
			gxtl[ch]->SetMarkerSize(0.18);
			gxtl[ch]->SetMarkerColor(kRed);
			gxtl[ch]->SetLineColor(kRed);
			gxtl[ch]->SetLineWidth(0.5);
			gxtl[ch]->Draw("PLSAME");
			ch_xt->Draw("t:d",Form("i==%d",ch),"SAMELP");
			f1[ch]->SetLineWidth(0.5);
			f2[ch]->SetLineWidth(0.5);
			f1[ch]->Draw("SAME");
			f2[ch]->Draw("SAME");
			TF1 * f3;
			double t1,x1;
			double t0,x0;
			if (cht<=5) t0 = t01;
			else t0 = t04;
			t1 = t0;
			if (cht==0||cht==3||cht==6||cht==9){
				f3  = new TF1("f3",Form("pow((%lf-%lf)*pow(x,5)+(%lf-%lf)*pow(x,4)+(%lf-%lf)*pow(x,3)+(%lf-%lf)*x*x+(%lf-%lf)*x+%lf-%lf,2)",
							f1[ch]->GetParameter(5),f2[ch]->GetParameter(5),
							f1[ch]->GetParameter(4),f2[ch]->GetParameter(4),
							f1[ch]->GetParameter(3),f2[ch]->GetParameter(3),
							f1[ch]->GetParameter(2),f2[ch]->GetParameter(2),
							f1[ch]->GetParameter(1),f2[ch]->GetParameter(1),
							f1[ch]->GetParameter(0),f2[ch]->GetParameter(0)),-2,2);
				//x1 = f3->GetX(0,-0.5,0.5);
				//t1 += f1[ch]->Eval(x1);
				// FIXME: for now no tuning
				t1 = t0;
				if (cht==0) {t11 = t1;}
				else if (cht==3) {t11 += t1;t11/=2.;}
				else if (cht==6) {t14 = t1;}
				else if (cht==9) {t14 += t1; t14/=2.;}
			}
			else{
				if (cht<6) t1 = t11;
				else t1 = t14;
				if (cht==2||cht==5||cht==8||cht==11){
					f3  = new TF1("f3",Form("abs(pow((%lf)*pow(x,5)+(%lf)*pow(x,4)+(%lf)*pow(x,3)+(%lf)*x*x+(%lf)*x+%lf-%lf,2))",
								f1[ch]->GetParameter(5),
								f1[ch]->GetParameter(4),
								f1[ch]->GetParameter(3),
								f1[ch]->GetParameter(2),
								f1[ch]->GetParameter(1),
								f1[ch]->GetParameter(0),t1-t0),-2,2);
					// FIXME: for now no tuning
					//x1 = f3->GetX(0,-2,0);
				}
				else{
					f3  = new TF1("f3",Form("abs(pow((%lf)*pow(x,5)+(%lf)*pow(x,4)+(%lf)*pow(x,3)+(%lf)*x*x+(%lf)*x+%lf-%lf,2))",
								f2[ch]->GetParameter(5),
								f2[ch]->GetParameter(4),
								f2[ch]->GetParameter(3),
								f2[ch]->GetParameter(2),
								f2[ch]->GetParameter(1),
								f2[ch]->GetParameter(0),t1-t0),-2,2);
					// FIXME: for now no tuning
					//x1 = f3->GetX(0,0,2);
				}
			} 

			// FIXME: for now no tuning
			//if (cht==1||cht==2) {x0 = dxu;dxu1+=(x0+x1)/2.;}
			//else if (cht==4||cht==5) {x0 = dyu;dyu1+=(x0+x1)/2.;}
			//else if (cht==7||cht==8) {x0 = dyd;dyd1+=(x0+x1)/2.;}
			//else if (cht==10||cht==11) {x0 = dxd;dxd1+=(x0+x1)/2.;}
			//x1+=x0;
			if (cht==1||cht==2) {x0 = dxu;;}
			else if (cht==4||cht==5) {x0 = dyu;}
			else if (cht==7||cht==8) {x0 = dyd;}
			else if (cht==10||cht==11) {x0 = dxd;}
			x1=x0;

			f3->SetLineWidth(0.5);
			f3->SetLineColor(kGreen);
			f3->Draw("SAME");
			TLatex * text2 = new TLatex();
			text2->SetTextSize(0.025);
			text2->SetText(-4,236,Form("Before: t0 = %lf, y0 = %lf",t0,x0));
			text2->Draw("SAME");
			TLatex * text3 = new TLatex();
			text3->SetTextSize(0.025);
			text3->SetText(-4,214,Form("After: t0 = %lf, y0 = %lf",t1,x1));
			text3->Draw("SAME");
		}
	}

	//===================Set Output ROOT File============================
	std::cout<<"# Set Output ROOT File"<<std::endl;
	TFile * o_file2 = new TFile(Form("../root/ana.%d."+suffix+"%d.root",startNo,iterationNo),"RECREATE");
	TTree * o_tree = new TTree("t","t");
	o_tree->Branch("dt",&dt,"dt[66]/D");
	o_tree->Branch("fd",&fd,"fd[66]/D");
	o_tree->Branch("dd",&dd,"dd[66]/D");
	o_tree->Branch("h",&h,"h[66]/I");
	o_tree->Branch("a",&a,"a[66]/D");
	o_tree->Branch("chi2X",&chi2X);
	o_tree->Branch("chi2Y",&chi2Y);
	o_tree->Branch("patX",&patX);
	o_tree->Branch("patY",&patY);
	o_tree->Branch("slX",&slX);
	o_tree->Branch("slY",&slY);
	o_tree->Branch("inX",&inX);
	o_tree->Branch("inY",&inY);
	int iEvent = 0;
	o_tree->Branch("iEvent",&iEvent);
#ifdef EVENTDISPLAY
//	o_tree->Branch("adc",&adc,"adc[66][32]/I");
//	o_tree->Branch("tdc",&tdc,"tdc[66][32]/I");
//	o_tree->Branch("clockNumberDriftTime",&clockNumberDriftTime,"clockNumberDriftTime[66][32]/I");
//	o_tree->Branch("n",&tdcNhit,"n[66]]/I");
#endif

	//===================Get reso and output ROOT file============================
	double xt_tmin_l[66];
	double xt_tmax_l[66];
	double xt_tmin_r[66];
	double xt_tmax_r[66];
	double xt_dmin_l[66];
	double xt_dmax_l[66];
	double xt_dmin_r[66];
	double xt_dmax_r[66];
	for ( int ch = 0; ch < 66; ch++){
		xt_tmin_l[ch] = f1[ch]->GetMinimum(-U,0);
		xt_tmax_l[ch] = f1[ch]->GetMaximum(-U,0);
		xt_tmin_r[ch] = f2[ch]->GetMinimum(0,U);
		xt_tmax_r[ch] = f2[ch]->GetMaximum(0,U);
		xt_dmin_l[ch] = f1[ch]->GetMinimumX(-U,0);
		xt_dmax_l[ch] = f1[ch]->GetMaximumX(-U,0);
		xt_dmin_r[ch] = f2[ch]->GetMinimumX(0,U);
		xt_dmax_r[ch] = f2[ch]->GetMaximumX(0,U);
		std::cout<<ch<<":"<<xt_tmax_l[ch]<<","
			              <<xt_tmin_l[ch]<<","
			              <<xt_tmin_r[ch]<<","
						  <<xt_tmax_r[ch]<<","
						  <<xt_dmax_l[ch]<<","
			              <<xt_dmin_l[ch]<<","
			              <<xt_dmin_r[ch]<<","
			              <<xt_dmax_r[ch]<<std::endl;
		std::cout<<ch<<"  f1:"<<f1[ch]->GetParameter(5)<<","
						  <<f1[ch]->GetParameter(4)<<","
						  <<f1[ch]->GetParameter(3)<<","
						  <<f1[ch]->GetParameter(2)<<","
						  <<f1[ch]->GetParameter(1)<<","
						  <<f1[ch]->GetParameter(0)<<std::endl;
		std::cout<<ch<<"  f2:"<<f2[ch]->GetParameter(5)<<","
						  <<f2[ch]->GetParameter(4)<<","
						  <<f2[ch]->GetParameter(3)<<","
						  <<f2[ch]->GetParameter(2)<<","
						  <<f2[ch]->GetParameter(1)<<","
						  <<f2[ch]->GetParameter(0)<<std::endl;
	}
	std::cout<<"Processing "<<N<<" Events ..."<<std::endl;
	for (iEvent = 0;iEvent<N; iEvent++){
		if (iEvent%1000==0) std::cout<<(double)iEvent/N*100<<"..."<<std::endl;
		c->GetEntry(iEvent);
		// FIXME:
		if (chi2X>chi2Xmax) continue;
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
			double t0;
			if (ch<=5) t0 = t01;
			else if (ch<=32) t0 = t02;
			else if (ch<=59) t0 = t03;
			else t0 = t04;
			double dist;
			double tmin, tmax;
			double dmin, dmax;
			TF1* f;
			if (fd[ch]<0){
				tmin = xt_tmin_l[ch]+t0;
				tmax = xt_tmax_l[ch]+t0;
				dmin = xt_dmin_l[ch];
				dmax = xt_dmax_l[ch];
				f = f1[ch];
			}
			else{
				tmin = xt_tmin_r[ch]+t0;
				tmax = xt_tmax_r[ch]+t0;
				dmin = xt_dmin_r[ch];
				dmax = xt_dmax_r[ch];
				f = f2[ch];
			}

			if (dt[ch]>tmax+10){
				dist=1e9;
			}
			else if (dt[ch]>tmax){
				dist = dmax;
			}
			else if (dt[ch]<tmin-10){
				dist=-1e9;
			}
			else if (dt[ch]<tmin){
				dist=(dt[ch]<0?0:dt[ch])/tmin*dmin;
			}
			else{
				if (fd[ch]<0){
					dist = f->GetX(dt[ch]-t0,dmax,dmin);
				}
				else{
					dist = f->GetX(dt[ch]-t0,dmin,dmax);
				}
				if (fabs(dist)>fabs(dmax)) dist = dmax;
				if (fabs(dist)<fabs(dmin)) dist = dmin;
			}

			dd[ch] = dist;

			int ir = (int)fabs(fd[ch])/1;
			double reso = fabs(fd[ch])-fabs(dist);
			if (fabs(reso)>0.5&&chp>=0) {
				if (ch<=32) hist_hmin_p3[chp]->Fill(h[ch]);
				else hist_hmin_p2[chp]->Fill(h[ch]);
			}
			if (ir<8){
				if (ch==19){
					resop3 = reso;
					fdp3 = fd[ch];
					// FIXME: no need to cut on Y
					//if (chi2X<=chi2Xmax&&chi2Y<=chi2Ymax){
					if (chi2X<=chi2Xmax){
						reso_p3[ir]->Fill(resop3);
						Neffi_p3[ir]++;
					}
					//printf("%lf-%lf=%lf\n",dist,fabs(fd[ch]),fabs(fd[ch])-dist);
					//printf("%d: p3[%d]: Mean = %lf, RMS = %lf\n",iEvent,ir,reso_p3[ir]->GetMean(),reso_p3[ir]->GetRMS());
				}
				else if (ch==46){
					resop2 = reso;
					fdp2 = fd[ch];
					// FIXME: no need to cut on Y
					//if (chi2X<=chi2Xmax&&chi2Y<=chi2Ymax){
					if (chi2X<=chi2Xmax){
						reso_p2[ir]->Fill(resop2);
						Neffi_p2[ir]++;
					}
				}
			}
		}
		o_tree->Fill();
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
		indiceX[0] = 0;
		indiceX[2] = 63;
		indiceY[0] = 3;
		indiceY[2] = 60;
		if (patX==1){
			indiceX[1] = 1;
			indiceX[3] = 64;
		}
		else if (patX==0){
			indiceX[1] = 2;
			indiceX[3] = 65;
		}
		else if (patX==3){
			indiceX[1] = 2;
			indiceX[3] = 64;
		}
		else if (patX==2){
			indiceX[1] = 1;
			indiceX[3] = 65;
		}
		if (patY==0){
			indiceY[1] = 4;
			indiceY[3] = 61;
		}
		else if (patY==1){
			indiceY[1] = 5;
			indiceY[3] = 62;
		}
		else if (patY==2){
			indiceY[1] = 5;
			indiceY[3] = 61;
		}
		else if (patY==3){
			indiceY[1] = 4;
			indiceY[3] = 62;
		}
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
			double t0;
			int ch;
			if (i<2) t0 = t01;
			else if (i==4) t0 = t02;
			else if (i==5) t0 = t03;
			else t0 = t04;
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
			buf<<"("<<dt[ch]-t0<<","<<dist<<","<<fd[ch]<<")";
			// FIXME
			//std::cout<<"dt["<<ch<<"] = "<<dt[ch]<<", t0 = "<<t0<<", dist = "<<dist<<", fd = "<<fd[ch]<<std::endl;
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
			double t0;
			int ch;
			if (i<2) t0 = t01;
			else if (i==4) t0 = t02;
			else if (i==5) t0 = t03;
			else t0 = t04;
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
			buf<<"("<<dt[ch]-t0<<","<<dist<<","<<fd[ch]<<")";
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
			double t0;
			if (ch<6) t0 = t01;
			else if (ch<33) t0 = t02;
			else if (ch<60) t0 = t03;
			else t0 = t04;
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
				textTDC[i][j]->SetText(clockNumberDriftTime[ch][j],MIN_WAVE_CH+10,Form("%d",(int)(tdc[ch][j]-t0)));
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
				time_fd = (f1[ch]->Eval(fd[ch])+t0+1000)*32./1000;
				time_dd = (f1[ch]->Eval(dd[ch])+t0+1000)*32./1000;
			}
			else{
				time_fd = (f2[ch]->Eval(fd[ch])+t0+1000)*32./1000;
				time_dd = (f2[ch]->Eval(dd[ch])+t0+1000)*32./1000;
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
		buf<<"EventDisplay/"<<displayName<<"/run"<<startNo<<"_"<<nFiles<<"."<<iEvent<<"."<<suffix<<"png";
		canvas2->SaveAs(buf.str().c_str());
		buf.str("");buf.clear();
		buf<<"EventDisplay/"<<displayName<<"/run"<<startNo<<"_"<<nFiles<<"."<<iEvent<<"."<<suffix<<"pdf";
		canvas2->SaveAs(buf.str().c_str());
#endif
	}
	o_tree->Write();
	o_file2->Close();

	//===================output new info============================
	//prepare for xt and para
	TFile * of_info = new TFile(Form("../info/info.%d."+suffix+"%d.root",startNo,iterationNo+1),"RECREATE");
	TTree * otr_para = new TTree("para","para");
	otr_para->Branch("t01",&t11);
	otr_para->Branch("t02",&t12);
	otr_para->Branch("t03",&t13);
	otr_para->Branch("t04",&t14);
	otr_para->Branch("dxu",&dxu1);
	otr_para->Branch("dyu",&dyu1);
	otr_para->Branch("dxd",&dxd1);
	otr_para->Branch("dyd",&dyd1);
	TTree * otr_xt = new TTree("xt","xt");
	double o_t,o_d,o_s;
	int o_i,o_n;
	otr_xt->Branch("t",&o_t);
	otr_xt->Branch("d",&o_d);
	otr_xt->Branch("sig",&o_s);
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
			double x0 = 0;
			double t0 = 0;
			double x1 = x0;
			double t1 = t0;
			if (cht==1||cht==2) {x0 = dxu;x1 = dxu1;}
			else if (cht==4||cht==5) {x0 = dyu;x1 = dyu1;}
			else if (cht==7||cht==8) {x0 = dyd;x1 = dyd1;}
			else if (cht==10||cht==11) {x0 = dxd;x1 = dxd1;}
			if (cht<=5) {t0=t01;t1=t11;}
			else {t0=t04;t1=t14;}
			o_t = v_center[index]+t0-t1;
			o_d = v_dist[index]+x0-x1;
			o_s = v_sigma[index];
			o_i = v_index[index];
			o_n = v_n[index];
			otr_xt->Fill();
			////FIXME
			//printf("%d,%d,%d:Fill: %lf+%lf-%lf=%lf,%lf+%lf-%lf=%lf)\n",cht,k,index,v_center[index],t0,t1,o_t,v_dist[index],x0,x1,o_d);
		}
	}
	canvas->SaveAs(Form("run%d.xt.tracker."+suffix+"%d.chi2_%d.pdf",startNo,iterationNo,(int)chi2Xmax));
	otr_xt->Write();
	otr_para->Fill();
	otr_para->Write();

	// Output Hmin
	TTree * otr_hmin = new TTree("hmin","hmin");
	otr_hmin->Branch("i",&idhmin);
	otr_hmin->Branch("h",&hhmin);
	//FIXME
	//if (iterationNo==0){
	if (0){
		for (int i = 0; i<18; i++){
			int ch;
			if (i<3) ch = 15-i;
			else if (i<6) ch = 20-i+3;
			else if (i<9)ch = 26-i+6;
			else if (i<12) ch = 42-i+9;
			else if (i<15) ch = 47-i+12;
			else ch = 53-i+15;
			TH1D * h;
			if (i<9) h = hist_hmin_p3[i];
			else h = hist_hmin_p2[i-9];
			int centeribin = h->GetMaximumBin();
			double centeri = h->GetBinCenter(centeribin);
			double heighti = h->GetBinContent(centeribin);
			double sigmai = 4;
			TF1 * f;
			f = new TF1(Form("hmin_%d",ch),"gaus",HMIN,HMAX);
			f->SetParameters(heighti,centeri,sigmai);
			h->Fit(Form("hmin_%d",ch),"N0","",centeri-sigmai*3,centeri+sigmai*3);
			double center = f->GetParameter(1);
			double sigma = f->GetParameter(2);
			Hmin[ch] = center+sigma*3;
		}
	}
	for(idhmin = 0; idhmin<66; idhmin++){
		hhmin = Hmin[idhmin];
		otr_hmin->Fill();
	}
	otr_hmin->Write();
	of_info->Close();

	//===================Draw others============================
	// Draw resolution of Prototype 2
	canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.975,"Spatial Resolution (R_{fit}-R_{hit}) of Prototype 2");
	text->Draw("SAME");
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
			reso_p2[index]->Fit(Form("f%d_%d",i,j),"N0","",reso_p2[index]->GetBinCenter(binc-8),reso_p2[index]->GetBinCenter(binc+8));
			//printf("Fit(%lf,%lf)\n",reso_p2[index]->GetBinCenter(binc-8),reso_p2[index]->GetBinCenter(binc+8));
			double center = f->GetParameter(1);
			double sigma = f->GetParameter(2);
			double binl = reso_p2[index]->FindBin(center-sigma*3);
			double binr = reso_p2[index]->FindBin(center+sigma*3);
			double effi = reso_p2[index]->Integral(binl,binr)/Neffi_p2[index];
			double effi2 = reso_p2[index]->Integral()/Neffi_p2[index];
			printf("effi[%d] = Integral(%lf,%lf)/%d = %lf/%d = %lf\n",index,binl,binr,Neffi_p2[index],reso_p2[index]->Integral(binl,binr),Neffi_p2[index],effi);
			reso_p2[index]->SetTitle(Form("%d-%dmm: center = %lf, reso = %lf, eff = %lf, eff2 = %lf",index,index+1,center,sigma,effi,effi2));
			reso_p2[index]->Draw();
			f->Draw("SAME");
			int max = reso_p3[index]->GetBinContent(reso_p3[index]->GetMaximumBin());
			TLine * line_l = new TLine();
			line_l->SetLineColor(kRed);
			line_l->SetLineWidth(0.5);
			line_l->SetX1(center-sigma*3);
			line_l->SetY1(0);
			line_l->SetX2(center-sigma*3);
			line_l->SetY2(max);
			TLine * line_r = new TLine();
			line_r->SetLineColor(kRed);
			line_r->SetLineWidth(0.5);
			line_r->SetX1(center+sigma*3);
			line_r->SetY1(0);
			line_r->SetX2(center+sigma*3);
			line_r->SetY2(max);
			line_l->Draw("SAME");
			line_r->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.reso_p2."+suffix+"%d.chi2_%d.pdf",startNo,iterationNo,(int)chi2Xmax));

	// Draw resolution of Prototype 3
	canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.975,"Spatial Resolution (R_{fit}-R_{hit}) of Prototype 3");
	text->Draw("SAME");
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
			reso_p3[index]->Fit(Form("f%d_%d",i,j),"N0","",reso_p3[index]->GetBinCenter(binc-8),reso_p3[index]->GetBinCenter(binc+8));
			double center = f->GetParameter(1);
			double sigma = f->GetParameter(2);
			double binl = reso_p3[index]->FindBin(center-sigma*3);
			double binr = reso_p3[index]->FindBin(center+sigma*3);
			double effi = reso_p3[index]->Integral(binl,binr)/Neffi_p3[index];
			double effi2 = reso_p3[index]->Integral()/Neffi_p3[index];
			// FIXME
			//printf("effi[%d] = Integral(%lf,%lf)/%d = %lf/%d = %lf\n",index,binl,binr,Neffi_p3[index],reso_p3[index]->Integral(binl,binr),Neffi_p3[index],effi);
			reso_p3[index]->SetTitle(Form("%d-%dmm: center = %lf, reso = %lf, eff = %lf, eff2 = %lf",index,index+1,center,sigma,effi,effi2));
			reso_p3[index]->Draw();
			f->Draw("SAME");
			int max = reso_p3[index]->GetBinContent(reso_p3[index]->GetMaximumBin());
			TLine * line_l = new TLine();
			line_l->SetLineColor(kRed);
			line_l->SetLineWidth(0.5);
			line_l->SetX1(center-sigma*3);
			line_l->SetY1(0);
			line_l->SetX2(center-sigma*3);
			line_l->SetY2(max);
			TLine * line_r = new TLine();
			line_r->SetLineColor(kRed);
			line_r->SetLineWidth(0.5);
			line_r->SetX1(center+sigma*3);
			line_r->SetY1(0);
			line_r->SetX2(center+sigma*3);
			line_r->SetY2(max);
			line_l->Draw("SAME");
			line_r->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.reso_p3."+suffix+"%d.chi2_%d.pdf",startNo,iterationNo,(int)chi2Xmax));

	// Draw XT for p3
	TF1 * f1_3[9];
	TF1 * f2_3[9];
	for (int i = 0; i<9; i++){
		f1_3[i] = new TF1(Form("f1_3_%d",i),"pol5",-9,2);
		f2_3[i] = new TF1(Form("f2_3_%d",i),"pol5",-2,9);
	}
	canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.975,"XT Curve in Prototype 3");
	text->Draw("SAME");
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
			gxtr[ch]->SetMarkerStyle(20);
			gxtr[ch]->SetMarkerSize(0.18);
			gxtr[ch]->SetMarkerColor(kRed);
			gxtr[ch]->SetLineColor(kRed);
			gxtr[ch]->SetLineWidth(0.5);
			gxtr[ch]->Draw("PLSAME");
			gxtl[ch]->SetMarkerStyle(20);
			gxtl[ch]->SetMarkerSize(0.18);
			gxtl[ch]->SetMarkerColor(kRed);
			gxtl[ch]->SetLineColor(kRed);
			gxtl[ch]->SetLineWidth(0.5);
			gxtl[ch]->Draw("PLSAME");
			f1[ch]->SetLineWidth(0.5);
			f2[ch]->SetLineWidth(0.5);
			f1[ch]->Draw("SAME");
			f2[ch]->Draw("SAME");
			//vddr.clear();
			//vdtr.clear();
			//vddl.clear();
			//vdtl.clear();
			//for ( int k  =0;  k< DDSTEP*2; k++){
			//	double dist = (k+0.5)*U/DDSTEP-8;
			//	TF1 * f  = new TF1(Form("f_%d_%d",k,chp),"gaus");
			//	TH1D * h = DTHISTS[k*66+ch];
			//	int centeribin = h->GetMaximumBin();
			//	double centeri = h->GetBinCenter(centeribin);
			//	double heighti = h->GetBinContent(centeribin);
			//	f->SetParameters(heighti,centeri,5);
			//	h->Fit(Form("f_%d_%d",k,chp),"N0","",centeri-5,centeri+5);
			//	double center = f->GetParameter(1);
			//	double sigma = f->GetParameter(2);
			//	h->SetTitle(Form("%lf mm: center = %lf, sigma = %lf",dist,center,sigma));
			//	//h->Draw();
			//	if (h->GetEntries()>500){
			//		if (dist>=0){vdtr.push_back(center); vddr.push_back(dist);}
			//		if (dist<=0){vdtl.push_back(center); vddl.push_back(dist);}
			//	}
			//}
			//TGraph * gxtl = new TGraph(vdtl.size(),&(vddl[0]),&(vdtl[0]));
			//TGraph * gxtr = new TGraph(vdtr.size(),&(vddr[0]),&(vdtr[0]));
			//gxtl->Fit(Form("f1_3_%d",chp),"N0","",-9,0);
			//gxtr->Fit(Form("f2_3_%d",chp),"N0","",0,9);
		}
	}
	canvas->SaveAs(Form("run%d.xt.p3."+suffix+"%d.chi2_%d.pdf",startNo,iterationNo,(int)chi2Xmax));

	// Draw XT for p2
	TF1 * f1_2[9];
	TF1 * f2_2[9];
	for (int i = 0; i<9; i++){
		f1_2[i] = new TF1(Form("f1_2_%d",i),"pol5",-9,2);
		f2_2[i] = new TF1(Form("f2_2_%d",i),"pol5",-2,9);
	}
	canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.975,"XT Curve in Prototype 2");
	text->Draw("SAME");
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
			gxtr[ch]->SetMarkerStyle(20);
			gxtr[ch]->SetMarkerSize(0.18);
			gxtr[ch]->SetMarkerColor(kRed);
			gxtr[ch]->SetLineColor(kRed);
			gxtr[ch]->SetLineWidth(0.5);
			gxtr[ch]->Draw("PLSAME");
			gxtl[ch]->SetMarkerStyle(20);
			gxtl[ch]->SetMarkerSize(0.18);
			gxtl[ch]->SetMarkerColor(kRed);
			gxtl[ch]->SetLineColor(kRed);
			gxtl[ch]->SetLineWidth(0.5);
			gxtl[ch]->Draw("PLSAME");
			f1[ch]->SetLineWidth(0.5);
			f2[ch]->SetLineWidth(0.5);
			f1[ch]->Draw("SAME");
			f2[ch]->Draw("SAME");
			//vddr.clear();
			//vdtr.clear();
			//vddl.clear();
			//vdtl.clear();
			//for ( int k  =0;  k< DDSTEP*2; k++){
			//	double dist = (k+0.5)*U/DDSTEP-8;
			//	TF1 * f  = new TF1(Form("f_%d_%d",k,chp),"gaus");
			//	TH1D * h = DTHISTS[k*66+ch];
			//	int centeribin = h->GetMaximumBin();
			//	double centeri = h->GetBinCenter(centeribin);
			//	double heighti = h->GetBinContent(centeribin);
			//	f->SetParameters(heighti,centeri,5);
			//	h->Fit(Form("f_%d_%d",k,chp),"N0","",centeri-5,centeri+5);
			//	double center = f->GetParameter(1);
			//	double sigma = f->GetParameter(2);
			//	h->SetTitle(Form("%lf mm: center = %lf, sigma = %lf",dist,center,sigma));
			//	//h->Draw();
			//	if (h->GetEntries()>500){
			//		if (dist>=0){vdtr.push_back(center); vddr.push_back(dist);}
			//		if (dist<=0){vdtl.push_back(center); vddl.push_back(dist);}
			//	}
			//}
			//TGraph * gxtl = new TGraph(vdtl.size(),&(vddl[0]),&(vdtl[0]));
			//TGraph * gxtr = new TGraph(vdtr.size(),&(vddr[0]),&(vdtr[0]));
			//gxtl->Fit(Form("f1_2_%d",chp),"N0","",-9,0);
			//gxtr->Fit(Form("f2_2_%d",chp),"N0","",0,9);
		}
	}
	canvas->SaveAs(Form("run%d.xt.p2."+suffix+"%d.chi2_%d.pdf",startNo,iterationNo,(int)chi2Xmax));

	// Draw beam spot
	canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.975,"Beam Spot At the Center of Different Chambers");
	text->Draw("SAME");
	for (int i = 0; i<2; i++){
		for (int j = 0; j<2; j++){
			int index = j*2+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./2*i,0.95/2*(1-j),1./2*(i+1),0.95/2*(2-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./2*i,0.95/2*(1-j),1./2*(i+1),0.95/2*(2-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<2; i++){
		for (int j = 0; j<2; j++){
			int index = j*2+i;
			pad[index]->cd();
			spot[index]->Draw("COLZ");
		}
	}
	canvas->SaveAs(Form("run%d.spot."+suffix+"%d.chi2_%d.pdf",startNo,iterationNo,(int)chi2Xmax));

	// Draw dt123
	canvas = new TCanvas("c","c",1024,768);
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
	canvas->SaveAs(Form("run%d.dt123."+suffix+"%d.chi2_%d.pdf",startNo,iterationNo,(int)chi2Xmax));

	// Draw tracking effi
	for(int i = 0; i<=128; i++){
		hist_effX->SetBinContent(i,(double)hist_chi2X->Integral(0,i)/N);
		hist_effY->SetBinContent(i,(double)hist_chi2Y->Integral(0,i)/N);
	}
	canvas = new TCanvas("c","c",1024,768);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	hist_effX->SetLineColor(kRed);
	hist_effY->SetLineColor(kBlue);
	hist_effX->Draw();
	hist_effY->Draw("SAME");
	canvas->SaveAs(Form("run%d.eff."+suffix+"%d.chi2_%d.pdf",startNo,iterationNo,(int)chi2Xmax));

	return 0;
}
