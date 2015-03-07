#include <iostream>
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
#include "TH2D.h"
#include "TGraph.h"
#include "TF1.h"
#include "TLatex.h"

#define NUM_OF_WAVE_SAMPLES 32

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
	double A1=0;
	double B1=0;
	double A2=0;
	double B2=0;
	if (argc>=10){
		A1 = strtof(argv[6],NULL);
		B1 = strtof(argv[7],NULL);
		A2 = strtof(argv[8],NULL);
		B2 = strtof(argv[9],NULL);
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
	c->SetBranchAddress("dt",dt);
	c->SetBranchAddress("dd",dd);
	c->SetBranchAddress("fd",fd);
	c->SetBranchAddress("chi2X",&chi2X);
	c->SetBranchAddress("chi2Y",&chi2Y);
	c->SetBranchAddress("slX",&slX);
	c->SetBranchAddress("inX",&inX);
	c->SetBranchAddress("slY",&slY);
	c->SetBranchAddress("inY",&inY);

	TCanvas * canvas;
	TPad * pad[100];
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);

	//===================Get para for previous run============================
	TFile * i_file = new TFile(Form("../info/XT.Tracker.%d."+suffix+"%d.root",startNo,iterationNo));
	if (!i_file) printf("Cannot find \"../info/XT.Tracker.%d."+suffix+"%d.root\"",startNo,iterationNo);
	TTree * i_tree2 = (TTree*) i_file->Get("para");
	double t01, t02, t03, t04;
	double dxu, dyu, dxd, dyd;
	i_tree2->SetBranchAddress("t01",&t01);
	i_tree2->SetBranchAddress("t02",&t02);
	i_tree2->SetBranchAddress("t03",&t03);
	i_tree2->SetBranchAddress("t04",&t04);
	i_tree2->SetBranchAddress("dxu",&dxu);
	i_tree2->SetBranchAddress("dyu",&dyu);
	i_tree2->SetBranchAddress("dxd",&dxd);
	i_tree2->SetBranchAddress("dyd",&dyd);
	i_tree2->GetEntry(0);
	i_file->Close();

	//===================Prepare Histograms============================
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
	double U =8; // mm
	int DDSTEP = 32; // steps
	int ICH = 1;
	std::vector<TH1D*> DTHISTS;
	for ( int i  =0;  i< DDSTEP*2; i++){
		for ( int j  =0;  j< 12; j++){
			name  = Form("dttr_%d_%d",i,j);
			title = name;
			TH1D * dthist = new TH1D(name,title,10+TMAX,-10,TMAX);
			DTHISTS.push_back(dthist);
		}
	}
	std::vector<TH1D*> DTP2HISTS;
	for ( int i  =0;  i< DDSTEP*2; i++){
		for ( int j  =0;  j< 9; j++){
			name  = Form("dtp2_%d_%d",i,j);
			title = name;
			TH1D * dthist = new TH1D(name,title,10+TMAX,-10,TMAX);
			DTP2HISTS.push_back(dthist);
		}
	}
	std::vector<TH1D*> DTP3HISTS;
	for ( int i  =0;  i< DDSTEP*2; i++){
		for ( int j  =0;  j< 9; j++){
			name  = Form("dtp3_%d_%d",i,j);
			title = name;
			TH1D * dthist = new TH1D(name,title,10+TMAX,-10,TMAX);
			DTP3HISTS.push_back(dthist);
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

	//===================Loop in Events============================
	Long64_t N = c->GetEntries();
	if (nEventMax&&nEventMax<N) N = nEventMax;
	printf("Processing... %d Events\n",N);
	for (int iEvent = 0;iEvent<N; iEvent++){
		if (iEvent%1000==0) printf("%lf%\n",(double)iEvent/N*100);
		c->GetEntry(iEvent);
		// FIXME:
		if (chi2X>=10||chi2Y>=10) continue;
		double t0;
		double x0 = 0;
		for ( int ch = 0; ch < 66; ch++){
			if (ch<=5) t0 = t01;
			else if (ch<=32) t0 = t02;
			else if (ch<=59) t0 = t03;
			else t0 = t04;
			dt[ch] -= t0;
			if (ch==1||ch==2) x0 = dxu;
			else if (ch==4||ch==5) x0 = dyu;
			else if (ch==61||ch==62) x0 = dyd;
			else if (ch==64||ch==65) x0 = dxd;
			dd[ch] -= x0;
			fd[ch] -= x0;
		}
		for ( int ch = 0; ch < 66; ch++){
			int cht = -1;
			if (ch<=5||ch>=60){
				cht = ch;
				if (ch>=60) cht-=54;
			}
			else{
				continue;
			}
			xt_tr[cht]->Fill(fd[ch],dt[ch]);
			if (fabs(fd[ch])>U) continue;
			int IDD = (int)((fd[ch]+U)/(U/DDSTEP));
			DTHISTS[IDD*12+cht]->Fill(dt[ch]);
		}
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
			if (ch<=5||ch>=60){
				xt_tr[cht]->Fill(fd[ch],dt[ch]);
			}
			else{
				if (chp>=0){
					if(ch<=32){
						xt_p3[chp]->Fill(fd[ch],dt[ch]);
						if (fabs(fd[ch])>U) continue;
						int IDD = (int)((fd[ch]+U)/(U/DDSTEP));
						DTP3HISTS[IDD*9+chp]->Fill(dt[ch]);
					}
					else{
						xt_p2[chp]->Fill(fd[ch],dt[ch]);
						if (fabs(fd[ch])>U) continue;
						int IDD = (int)((fd[ch]+U)/(U/DDSTEP));
						DTP2HISTS[IDD*9+chp]->Fill(dt[ch]);
					}
				}
			}
			if ((ch==19||ch==46)&&fabs(dt[ch]-498)<500){
				double dist;
				if (ch==19)
					dist = sqrt(A1*dt[ch]+B1*B1)-B1;
				else
					dist = sqrt(A2*dt[ch]+B2*B2)-B2;
				int ir = (int)fabs(fd[ch])/1;
				if (ir<8) {
					if (ch==19){
						reso_p3[ir]->Fill(fabs(fd[ch])-dist);
						Neffi_p3[ir]++;
						//printf("%lf-%lf=%lf\n",dist,fabs(fd[ch]),fabs(fd[ch])-dist);
						//printf("%d: p3[%d]: Mean = %lf, RMS = %lf\n",iEvent,ir,reso_p3[ir]->GetMean(),reso_p3[ir]->GetRMS());
					}
					else if (ch==46){
						reso_p2[ir]->Fill(fabs(fd[ch])-dist);
						Neffi_p2[ir]++;
					}
				}
			}
		}
	}

	//===================Find Peaks============================
	TChain * i_tree = new TChain("xt","xt");
	i_tree ->Add(Form("../info/XT.Tracker.%d."+suffix+"%d.root",startNo,iterationNo));
	//double i_t,i_d,i_s;
	//int i_i;
	//i_tree->SetBranchAddress("t",&i_t);
	//i_tree->SetBranchAddress("d",&i_d);
	//i_tree->SetBranchAddress("sig",&i_s);
	//i_tree->SetBranchAddress("i",&i_i);
	i_tree->SetMarkerStyle(20);
	i_tree->SetMarkerSize(0.15);
	i_tree->SetMarkerColor(kBlue);
	i_tree->SetLineColor(kBlue);
	i_tree->SetLineWidth(0.5);

	TFile * o_file = new TFile(Form("../info/XT.Tracker.%d."+suffix+"%d.root",startNo,iterationNo+1),"RECREATE");
	TTree * o_tree2 = new TTree("para","para");
	double t11, t12, t13, t14;
	double dxu1, dyu1, dxd1, dyd1;
	o_tree2->Branch("t01",&t11);
	o_tree2->Branch("t02",&t12);
	o_tree2->Branch("t03",&t13);
	o_tree2->Branch("t04",&t14);
	o_tree2->Branch("dxu",&dxu1);
	o_tree2->Branch("dyu",&dyu1);
	o_tree2->Branch("dxd",&dxd1);
	o_tree2->Branch("dyd",&dyd1);
	t11  = t01;
	t12  = t02;
	t13  = t03;
	t14  = t04;
	dxu1 = dxu;
	dyu1 = dyu;
	dxd1 = dxd;
	dyd1 = dyd;
	TTree * o_tree = new TTree("xt","xt");
	double o_t,o_d,o_s;
	int o_i,o_n;
	o_tree->Branch("t",&o_t);
	o_tree->Branch("d",&o_d);
	o_tree->Branch("sig",&o_s);
	o_tree->Branch("i",&o_i);
	o_tree->Branch("n",&o_n);
	// Draw XT for trackers
	canvas = new TCanvas("c","c",1024,768);
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
	std::vector<double> vddr;
	std::vector<double> vdtr;
	std::vector<double> vddl;
	std::vector<double> vdtl;
	std::vector<double> v_n;
	std::vector<double> v_center;
	std::vector<double> v_dist;
	std::vector<double> v_sigma;
	std::vector<double> v_index;
	for (int i = 0; i<4; i++){
		for (int j = 0; j<3; j++){
			int index = i*3+j;
			pad[index]->cd();
			vddr.clear();
			vdtr.clear();
			vddl.clear();
			vdtl.clear();
			for ( int k  =0;  k< DDSTEP*2; k++){
				double dist = (k+0.5)*U/DDSTEP-8;
				TF1 * f  = new TF1(Form("f_%d_%d",k,index),"gaus",-10,TMAX);
				int centeribin = DTHISTS[k*12+index]->GetMaximumBin();
				double centeri = DTHISTS[k*12+index]->GetBinCenter(centeribin);
				double heighti = DTHISTS[k*12+index]->GetBinContent(centeribin);
				f->SetParameters(heighti,centeri,2.5);
				DTHISTS[k*12+index]->Fit(Form("f_%d_%d",k,index),"N0","",centeri-5,centeri+5);
				double center = f->GetParameter(1);
				double sigma = f->GetParameter(2);
				DTHISTS[k*12+index]->SetTitle(Form("%lf mm: centeri = %lf, heighti = %lf; center = %lf, sigma = %lf",dist,centeri,heighti,center,sigma));
				v_n.push_back(DTHISTS[k*12+index]->GetEntries());
				v_center.push_back(center);
				v_dist.push_back(dist);
				v_sigma.push_back(sigma);
				v_index.push_back(index);
				////FIXME
				//if (index==0){
				//	DTHISTS[k*12+index]->Draw();
				//	f->Draw("SAME");
				//	canvas->SaveAs(Form("%d.pdf",k));
				//}
				////FIXME
				printf("%d,%d,%d:v_center.push_back(%lf); v_dist.push_back(%lf)\n",index,k,v_center.size(),center,dist);
				if (DTHISTS[k*12+index]->GetEntries()>100&&sigma<8){
					if (dist>=0){vdtr.push_back(center); vddr.push_back(dist);}
					if (dist<=0){vdtl.push_back(center); vddl.push_back(dist);}
				}
			}
			xt_tr[index]->Draw("COLZ");
			TGraph * gxtl = new TGraph(vdtl.size(),&(vddl[0]),&(vdtl[0]));
			TGraph * gxtr = new TGraph(vdtr.size(),&(vddr[0]),&(vdtr[0]));
			gxtr->SetMarkerStyle(20);
			gxtr->SetMarkerSize(0.18);
			gxtr->SetMarkerColor(kRed);
			gxtr->SetLineColor(kRed);
			gxtr->SetLineWidth(0.5);
			gxtr->Draw("PLSAME");
			gxtl->SetMarkerStyle(20);
			gxtl->SetMarkerSize(0.18);
			gxtl->SetMarkerColor(kRed);
			gxtl->SetLineColor(kRed);
			gxtl->SetLineWidth(0.5);
			gxtl->Draw("PLSAME");
			i_tree->Draw("t:d",Form("i==%d",index),"SAMELP");
			TF1 * f1  = new TF1("f1","pol5",-9,2);
			TF1 * f2  = new TF1("f2","pol5",-2,9);
			gxtl->Fit("f1","N0","",-9,0);
			gxtr->Fit("f2","N0","",0,9);
			f1->SetLineWidth(0.5);
			f2->SetLineWidth(0.5);
			f1->Draw("SAME");
			f2->Draw("SAME");
			TF1 * f3;
			double t1,x1;
			double t0,x0;
			if (index<=5) t0 = t01;
			else t0 = t04;
			if (index==0||index==3||index==6||index==9){
				f3  = new TF1("f3",Form("pow((%lf-%lf)*pow(x,5)+(%lf-%lf)*pow(x,4)+(%lf-%lf)*pow(x,3)+(%lf-%lf)*x*x+(%lf-%lf)*x+%lf-%lf,2)",
							f1->GetParameter(5),f2->GetParameter(5),
							f1->GetParameter(4),f2->GetParameter(4),
							f1->GetParameter(3),f2->GetParameter(3),
							f1->GetParameter(2),f2->GetParameter(2),
							f1->GetParameter(1),f2->GetParameter(1),
							f1->GetParameter(0),f2->GetParameter(0)),-2,2);
				x1 = f3->GetX(0,-0.5,0.5);
				t1 = f1->Eval(x1);
				if (index==0) {t1+=t01; t11 = t1;}
				else if (index==3) {t1+=t01;t11 += t1;t11/=2.;}
				else if (index==6) {t1+=t04;t14 = t1;}
				else if (index==9) {t1+=t04;t14 += t1; t14/=2.;}
			}
			else{
				if (index<6) t1 = t11;
				else t1 = t14;
				if (index==2||index==5||index==8||index==11){
					f3  = new TF1("f3",Form("abs(pow((%lf)*pow(x,5)+(%lf)*pow(x,4)+(%lf)*pow(x,3)+(%lf)*x*x+(%lf)*x+%lf-%lf,2))",
								f1->GetParameter(5),
								f1->GetParameter(4),
								f1->GetParameter(3),
								f1->GetParameter(2),
								f1->GetParameter(1),
								f1->GetParameter(0),t1-t0),-2,2);
					x1 = f3->GetX(0,-2,0);
				}
				else{
					f3  = new TF1("f3",Form("abs(pow((%lf)*pow(x,5)+(%lf)*pow(x,4)+(%lf)*pow(x,3)+(%lf)*x*x+(%lf)*x+%lf-%lf,2))",
								f2->GetParameter(5),
								f2->GetParameter(4),
								f2->GetParameter(3),
								f2->GetParameter(2),
								f2->GetParameter(1),
								f2->GetParameter(0),t1-t0),-2,2);
					x1 = f3->GetX(0,0,2);
				}
			} 
			if (index==1||index==2) {x0 = dxu;dxu1+=(x0+x1)/2;}
			else if (index==4||index==5) {x0 = dyu;dyu1+=(x0+x1)/2;}
			else if (index==7||index==8) {x0 = dyd;dyd1+=(x0+x1)/2;}
			else if (index==10||index==11) {x0 = dxd;dxd1+=(x0+x1)/2;}
			x1+=x0;
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
	for ( int cht = 0; cht<12; cht++){
		for ( int k  =0;  k< DDSTEP*2; k++){
			int index = cht*DDSTEP*2+k;
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
			o_tree->Fill();
			////FIXME
			printf("%d,%d,%d:Fill: %lf+%lf-%lf=%lf,%lf+%lf-%lf=%lf)\n",cht,k,index,v_center[index],t0,t1,o_t,v_dist[index],x0,x1,o_d);
		}
	}
	canvas->SaveAs(Form("run%d.xt.tracker."+suffix+"%d.pdf",startNo,iterationNo));
	o_tree->Write();
	o_tree2->Fill();
	o_tree2->Write();
	o_file->Close();

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
			int index = j*3+i;
			pad[index]->cd();
			xt_p3[index]->Draw("COLZ");
			vddr.clear();
			vdtr.clear();
			vddl.clear();
			vdtl.clear();
			for ( int k  =0;  k< DDSTEP*2; k++){
				double dist = (k+0.5)*U/DDSTEP-8;
				TF1 * f  = new TF1(Form("f_%d_%d",k,index),"gaus");
				int centeribin = DTHISTS[k*12+index]->GetMaximumBin();
				double centeri = DTHISTS[k*12+index]->GetBinCenter(centeribin);
				double heighti = DTHISTS[k*12+index]->GetBinContent(centeribin);
				f->SetParameters(heighti,centeri,5);
				DTHISTS[k*12+index]->Fit(Form("f_%d_%d",k,index),"N0","",centeri-5,centeri+5);
				double center = f->GetParameter(1);
				double sigma = f->GetParameter(2);
				DTHISTS[k*12+index]->SetTitle(Form("%lf mm: center = %lf, sigma = %lf",dist,center,sigma));
				//DTHISTS[k*12+index]->Draw();
				if (DTHISTS[k*12+index]->GetEntries()>500){
					if (dist>=0){vdtr.push_back(center); vddr.push_back(dist);}
					if (dist<=0){vdtl.push_back(center); vddl.push_back(dist);}
				}
			}
			TGraph * gxtl = new TGraph(vdtl.size(),&(vddl[0]),&(vdtl[0]));
			TGraph * gxtr = new TGraph(vdtr.size(),&(vddr[0]),&(vdtr[0]));
			gxtl->Fit(Form("f1_3_%d",index),"N0","",-9,0);
			gxtr->Fit(Form("f2_3_%d",index),"N0","",0,9);
		}
	}
	canvas->SaveAs(Form("run%d.xt.p3."+suffix+"%d.pdf",startNo,iterationNo));

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
			int index = j*3+i;
			pad[index]->cd();
			xt_p2[index]->Draw("COLZ");
			vddr.clear();
			vdtr.clear();
			vddl.clear();
			vdtl.clear();
			for ( int k  =0;  k< DDSTEP*2; k++){
				double dist = (k+0.5)*U/DDSTEP-8;
				TF1 * f  = new TF1(Form("f_%d_%d",k,index),"gaus");
				int centeribin = DTHISTS[k*9+index]->GetMaximumBin();
				double centeri = DTHISTS[k*9+index]->GetBinCenter(centeribin);
				double heighti = DTHISTS[k*9+index]->GetBinContent(centeribin);
				f->SetParameters(heighti,centeri,5);
				DTHISTS[k*9+index]->Fit(Form("f_%d_%d",k,index),"N0","",centeri-5,centeri+5);
				double center = f->GetParameter(1);
				double sigma = f->GetParameter(2);
				DTHISTS[k*9+index]->SetTitle(Form("%lf mm: center = %lf, sigma = %lf",dist,center,sigma));
				//DTHISTS[k*+index]->Draw();
				if (DTHISTS[k*9+index]->GetEntries()>500){
					if (dist>=0){vdtr.push_back(center); vddr.push_back(dist);}
					if (dist<=0){vdtl.push_back(center); vddl.push_back(dist);}
				}
			}
			TGraph * gxtl = new TGraph(vdtl.size(),&(vddl[0]),&(vdtl[0]));
			TGraph * gxtr = new TGraph(vdtr.size(),&(vddr[0]),&(vdtr[0]));
			gxtl->Fit(Form("f1_2_%d",index),"N0","",-9,0);
			gxtr->Fit(Form("f2_2_%d",index),"N0","",0,9);
		}
	}
	canvas->SaveAs(Form("run%d.xt.p2."+suffix+"%d.pdf",startNo,iterationNo));

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
	canvas->SaveAs(Form("run%d.spot."+suffix+"%d.pdf",startNo,iterationNo));

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
	canvas->SaveAs(Form("run%d.dt123."+suffix+"%d.pdf",startNo,iterationNo));

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
			TF1 * f  = new TF1(Form("f%d_%d",i,j),"gaus");
			f->SetParameters(reso_p2[index]->GetBinContent(reso_p2[index]->GetMaximumBin()),0,0.2);
			int binc = reso_p2[index]->GetMaximumBin();
			reso_p2[index]->Fit(Form("f%d_%d",i,j),"","",reso_p2[index]->GetBinCenter(binc-8),reso_p2[index]->GetBinCenter(binc+8));
			//printf("Fit(%lf,%lf)\n",reso_p2[index]->GetBinCenter(binc-8),reso_p2[index]->GetBinCenter(binc+8));
			double center = f->GetParameter(1);
			double sigma = f->GetParameter(2);
			double binl = reso_p2[index]->FindBin(center-sigma*3);
			double binr = reso_p2[index]->FindBin(center+sigma*3);
			double effi = reso_p2[index]->Integral(binl,binr)/Neffi_p2[index];
			//printf("effi[%d] = Integral(%lf,%lf)/%d = %lf/%d = %lf\n",index,binl,binr,Neffi_p2[index],reso_p2[index]->Integral(binl,binr),Neffi_p2[index],effi);
			reso_p2[index]->SetTitle(Form("%d-%dmm: center = %lf, reso = %lf, eff = %lf",index*2,index*2+2,center,sigma,effi));
			reso_p2[index]->Draw();
			//f->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.reso_p2."+suffix+"%d.pdf",startNo,iterationNo));

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
			TF1 * f  = new TF1(Form("f%d_%d",i,j),"gaus");
			f->SetParameters(reso_p3[index]->GetBinContent(reso_p3[index]->GetMaximumBin()),0,0.2);
			int binc = reso_p3[index]->GetMaximumBin();
			reso_p3[index]->Fit(Form("f%d_%d",i,j),"","",reso_p3[index]->GetBinCenter(binc-8),reso_p3[index]->GetBinCenter(binc+8));
			double center = f->GetParameter(1);
			double sigma = f->GetParameter(2);
			double binl = reso_p3[index]->FindBin(center-sigma*3);
			double binr = reso_p3[index]->FindBin(center+sigma*3);
			double effi = reso_p3[index]->Integral(binl,binr)/Neffi_p3[index];
			//printf("effi[%d] = Integral(%lf,%lf)/%d = %lf/%d = %lf\n",index,binl,binr,Neffi_p3[index],reso_p3[index]->Integral(binl,binr),Neffi_p3[index],effi);
			reso_p3[index]->SetTitle(Form("%d-%dmm: center = %lf, reso = %lf, eff = %lf",index*2,index*2+2,center,sigma,effi));
			reso_p3[index]->Draw();
			//f->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.reso_p3."+suffix+"%d.pdf",startNo,iterationNo));

	return 0;
}
