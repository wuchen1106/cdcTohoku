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
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"

void print_usage(char* prog_name);
int chg2cht(int i);
int chg2chp2(int i);
int chg2chp3(int i);
int cht2chg(int i);
int chp22chg(int i);
int chp32chg(int i);
int get_bid(int i);
int get_bid_core(int i);
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

	//===================Get input ROOT file============================
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

	//===================Prepare output ROOT file============================
	char outputName[128];
	int o_iHits[66];
	int o_nHits[66];
	int o_pedestalN[66];
	double o_pedestal[66];
	double o_pedestalChi2[66];
	double o_areaall[66];
	double o_area[66];
	double o_time[66];
	int o_height[66];
	sprintf(outputName,("../root/d_%d."+suffix+"root").c_str(),runNo);
	TFile * f = new TFile(outputName,"RECREATE");
	TTree * t = new TTree("t","t");
	t->Branch("i",o_iHits,"i[66]/I");
	t->Branch("n",o_nHits,"n[66]/I");
	t->Branch("p",o_pedestal,"p[66]/D");
	t->Branch("pn",o_pedestalN,"pn[66]/I");
	t->Branch("pchi2",o_pedestalChi2,"pchi2[66]/D");
	t->Branch("aa",o_areaall,"aa[66]/D");
	t->Branch("a",o_area,"a[66]/D");
	t->Branch("dt",o_time,"dt[66]/D");
	t->Branch("h",o_height,"h[66]/I");
	//FIXME better not to record
	//t->Branch("clockNumberDriftTime",clockNumberDriftTime,"clockNumberDriftTime[66][32]/I");
	//t->Branch("adc",adc,"adc[66][32]/I");

	//===================Prepare Histograms============================
	TCanvas * canvas;
	TPad * pad[100];
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);

	TLine *line_ht[12];
	TLatex *text_ht[12];
	for (Int_t i=0; i<12; i++) {
		line_ht[i] = new TLine();
		line_ht[i]->SetLineColor(kRed);
		line_ht[i]->SetLineWidth(0.5);
		text_ht[i] = new TLatex(0,0,"");
		text_ht[i]->SetTextSize(0.04);
		text_ht[i]->SetTextColor(kRed);
	}

	TLatex * text = new TLatex();
	text->SetTextSize(0.02);
	TString name;
	TString title;
	TString chamberName;
	TH2D* hdt_tr[12];
	TH2D* hdt_tr_fp[12];
	TH2D* hdt_p2[9];
	TH2D* hdt_p2_fp[9];
	TH2D* hdt_p3[9];
	TH2D* hdt_p3_fp[9];
	TH1D* n_tr[12];
	TH1D* n_p2[9];
	TH1D* n_p3[9];
	TH1D* n_tr_ac[12];
	TH1D* n_p2_ac[9];
	TH1D* n_p3_ac[9];
	for (int i = 0; i<12; i++){
		name = Form("hdt_tr_all_%d",i);
		if (i<3) chamberName = "X Tracker Up";
		else if (i<6) chamberName = "Y Tracker Up";
		else if (i<9) chamberName = "Y Tracker Down";
		else if (i<12) chamberName = "X Tracker Down";
		int wireID = i;
		if (i>=6) wireID+= 54;
		title = Form(chamberName+", wire#%d",wireID);
		hdt_tr[i] = new TH2D(name,title,250,-100,900,250,200,700);
		name = Form("hdt_tr_fp_%d",i);
		hdt_tr_fp[i] = new TH2D(name,title,250,-100,900,250,200,700);
		hdt_tr_fp[i]->GetXaxis()->SetTitle("Time [ns]");
		hdt_tr_fp[i]->GetYaxis()->SetTitle("Peak Height [ADC ch]");
		hdt_tr[i]->GetXaxis()->SetTitle("Time [ns]");
		hdt_tr[i]->GetYaxis()->SetTitle("Peak Height [ADC ch]");
		name = Form("n_tr_%d",i);
		n_tr[i] = new TH1D(name,title,30,0,15);
		n_tr_ac[i] = new TH1D(name,title,30,0,15);
		n_tr_ac[i]->SetLineColor(kRed);
	}
	for (int i = 0; i<9; i++){
		name = Form("hdt_p3_all_%d",i);
		int wireID;
		if (i<3) wireID = 15-i%3;
		else if (i<6) wireID = 20-i%3;
		else if (i<9) wireID = 26-i%3;
		title = Form("wire#%d",wireID);
		hdt_p3[i] = new TH2D(name,title,250,-100,900,250,200,700);
		name = Form("hdt_p3_fp_%d",i);
		hdt_p3_fp[i] = new TH2D(name,title,250,-100,900,250,200,700);
		hdt_p3_fp[i]->GetXaxis()->SetTitle("Time [ns]");
		hdt_p3_fp[i]->GetYaxis()->SetTitle("Peak Height [ADC ch]");
		hdt_p3[i]->GetXaxis()->SetTitle("Time [ns]");
		hdt_p3[i]->GetYaxis()->SetTitle("Peak Height [ADC ch]");
		name = Form("n_p3_%d",i);
		n_p3[i] = new TH1D(name,title,30,0,15);
		n_p3_ac[i] = new TH1D(name,title,30,0,15);
		n_p3_ac[i]->SetLineColor(kRed);
	}
	for (int i = 0; i<9; i++){
		name = Form("hdt_p2_all_%d",i);
		int wireID;
		if (i<3) wireID = 42-i%3;
		else if (i<6) wireID = 47-i%3;
		else if (i<9) wireID = 53-i%3;
		title = Form("wire#%d",wireID);
		hdt_p2[i] = new TH2D(name,title,250,-100,900,250,200,700);
		name = Form("hdt_p2_fp_%d",i);
		hdt_p2_fp[i] = new TH2D(name,title,250,-100,900,250,200,700);
		hdt_p2_fp[i]->GetXaxis()->SetTitle("Time [ns]");
		hdt_p2_fp[i]->GetYaxis()->SetTitle("Peak Height [ADC ch]");
		hdt_p2[i]->GetXaxis()->SetTitle("Time [ns]");
		hdt_p2[i]->GetYaxis()->SetTitle("Peak Height [ADC ch]");
		name = Form("n_p2_%d",i);
		n_p2[i] = new TH1D(name,title,30,0,15);
		n_p2_ac[i] = new TH1D(name,title,30,0,15);
		n_p2_ac[i]->SetLineColor(kRed);
	}
	// hists
	TH1D * h_dt[66];
	TH1D * h_h[66];
	//TH1D * h_h2[66];
	TH1D * h_h3[66];
	TH1D * h_h31[66];
	TH1D * h_h32[66];
	double pedestal[66]={0};
	for ( int i = 0; i<66; i++){
		h_dt[i] = new TH1D(Form("h_dt%d",i),"h_dt",60,-860,-800);
		h_h[i] = new TH1D(Form("h_h_%d",i),"h_h",500,200,700);
		//h_h2[i] = new TH1D(Form("h_h2_%d",i),"h_h2",500,200,700);
		h_h3[i] = new TH1D(Form("h_h3_%d",i),"h_h3",300,200,500);
		h_h31[i] = new TH1D(Form("h_h31_%d",i),"h_h31",50,300,700);
		h_h32[i] = new TH1D(Form("h_h32_%d",i),"h_h32",100,200,300);
	}

	// Loop in events
	Long64_t N = c->GetEntries();
	if (nEventMax&&nEventMax<N) N = nEventMax;
	std::cout<<"Processing "<<N<<" events..."<<std::endl;
	int tdcNhitwire = -1;
	//N=1;
	for (Long64_t i = 0;i<N; i++){
		if (i%1000==0) std::cout<<(double)i/N*100<<"%..."<<std::endl;
		c->GetEntry(i);
		for(int iwire = 0; iwire<66; iwire++){
			tdcNhitwire = tdcNhit[iwire];
			for ( int ihit = 0; ihit<tdcNhitwire; ihit++){
				h_dt[iwire]->Fill(driftTime[iwire][ihit]);
				int height=0;
				for(int clk = clockNumberDriftTime[iwire][ihit]; clk<(ihit+1>=tdcNhitwire?32:clockNumberDriftTime[iwire][ihit+1]); clk++){
					if (adc[iwire][clk]<height){
						break;
					}
					height=adc[iwire][clk];
				}
				if (driftTime[iwire][ihit]<-850)
					h_h[iwire]->Fill(height);
				//else if (driftTime[iwire][ihit]>-400)
				//	h_h2[iwire]->Fill(height);
				if (driftTime[iwire][ihit]<-725&&driftTime[iwire][ihit]>-825){
					h_h3[iwire]->Fill(height);
					h_h31[iwire]->Fill(height);
					h_h32[iwire]->Fill(height);
				}
			}
			int clk;
			double ped = 0;
			for(clk = 0; clk<(tdcNhitwire<=0?32:clockNumberDriftTime[iwire][0]-1); clk++){
				ped+=adc[iwire][clk];
			}
			if (clk>0) ped/=(double)clk;
			pedestal[iwire] += ped;
		}
	}

	// Get t0 and hmin
	double Hmin[66];
	TF1 *f1 = 0;
	TF1 *f2 = 0;
	TF1 *f3 = 0;
	for ( int i = 0; i<66; i++){
		if(h_h31[i]->GetMaximumBin()>5){
			if (f2) delete f2;
			f2 = new TF1(Form("hmin2_%d",i),"gaus",200,300);
			double center = h_h32[i]->GetBinCenter(h_h32[i]->GetMaximumBin());
			h_h32[i]->Fit(Form("hmin2_%d",i),"qN0","",center-8,center+8);
			Hmin[i] = f2->GetParameter(1)+6*f2->GetParameter(2);
			printf("super large gas gain!\n");
			printf("Hmin[%d] = %lf; %d\n",i,Hmin[i],h_h31[i]->GetMaximumBin());
		}
		else{
			double h_maxbin_300 = h_h32[i]->GetBinCenter(h_h32[i]->GetMaximumBin());
			double h_ratio_300 = ((double)h_h32[i]->GetBinContent(99))/h_h32[i]->GetBinContent(h_h32[i]->GetMaximumBin());
			printf("ratio_300 = %lf, %lf - %lf;\n",h_ratio_300,h_maxbin_300,pedestal[i]/N);
			if (h_ratio_300>0.01&&h_maxbin_300-pedestal[i]/N<8){
				if (f2) delete f2;
				f2 = new TF1(Form("hmin2_%d",i),"gaus",200,300);
				double center = h_h32[i]->GetBinCenter(h_h32[i]->GetMaximumBin());
				f2->SetParameters(h_h32[i]->Integral(),center,10);
				h_h32[i]->Fit(Form("hmin2_%d",i),"qN0","",center-8,center+8);
				if (f2->GetParameter(1)-center>2&&f2->GetParameter(2)>3)
					Hmin[i] = center+2;
				else
					Hmin[i] = f2->GetParameter(1)+5*f2->GetParameter(2);
				printf("Hmin[%d] = %lf; (%lf,%lf,%lf)\n",i,Hmin[i],f2->GetParameter(1),f2->GetParameter(2),center);
			}
			else{
				if (f1) delete f1;
				f1 = new TF1(Form("hmin_%d",i),"gaus",200,700);
				double center = h_h[i]->GetBinCenter(h_h[i]->GetMaximumBin());
				f1->SetParameters(h_h[i]->Integral(),center,10);
				h_h[i]->Fit(Form("hmin_%d",i),"qN0","",center-8,center+8);
				double hmin = f1->GetParameter(1)+3*f1->GetParameter(2);

				if (f3) delete f3;
				f3 = new TF1(Form("hmin3_%d",i),"landau",200,500);
				center = h_h3[i]->GetBinCenter(h_h3[i]->GetMaximumBin());
				f3->SetParameters(h_h3[i]->Integral(),center,10);
				h_h3[i]->Fit(Form("hmin3_%d",i),"qN0","",center-30,center+60);
				double hmin3 = f3->GetParameter(1)-2.5*f3->GetParameter(2);

				if ((double)h_h[i]->GetEntries()>N*8.e-3&&hmin>hmin3)
					Hmin[i] = hmin;
				else
					if (hmin3>600)
						Hmin[i] = f1->GetParameter(1)+f1->GetParameter(2);
					else
						Hmin[i] = hmin3;


				printf("Hmin[%d] = %lf; (%lf,%lf,%lf) (%lf,%lf,%lf)\n",i,Hmin[i],f3->GetParameter(1),f3->GetParameter(2),center,f1->GetParameter(1),f1->GetParameter(2),h_h[i]->GetEntries()/((double)N));
			}
		}
		// FIXME
		if (i>=6&&i<60)
			Hmin[i] = 0;
	}
	double t0[4];
	double t0_nmax[4];
	double t0_temp[66];
	for(int i = 0; i <4; i++){
		std::vector<int> indice;
		if (i == 0){for(int j = 0; j<6; j++) indice.push_back(j);}
		if (i == 1){indice.push_back(15);indice.push_back(19);indice.push_back(24);}
		if (i == 2){indice.push_back(42);indice.push_back(46);indice.push_back(51);}
		if (i == 3){for(int j = 60; j<66; j++) indice.push_back(j);}
		t0[i] = 0;
		for ( int j = 0; j<indice.size(); j++){
			int ch = indice[j];
			int maxbin = h_dt[ch]->GetMaximumBin();
			double max = h_dt[ch]->GetBinContent(maxbin);
			if (t0_nmax[i]<max) t0_nmax[i] = max;
			double min = 0;
			for (int k = 1; k<=5; k++){
				min += h_dt[ch]->GetBinContent(k);
			}
			min /= 5.;
			double th = (max-min)/10+min;
			int ibin = 1;
			for (; ibin<=60; ibin++){
				int height = h_dt[ch]->GetBinContent(ibin);
				if (height>th) break;
			}
			t0_temp[indice[j]] = h_dt[ch]->GetBinCenter(ibin);
			t0[i] += t0_temp[indice[j]];
			printf("#\t%d (%lf,%lf)%5=%lf, %d %lf\n",ch,min,max,th,ibin,h_dt[ch]->GetBinCenter(ibin));
		}
		int N = indice.size();
		double average = t0[i]/(N);
		for ( int j = 0; j<indice.size(); j++){
			if (fabs(average-t0_temp[indice[j]])>5){
				printf("#\t%d %lf-%lf>5\n",indice[j],average,t0_temp[indice[j]]);
				t0[i]-=t0_temp[indice[j]];
				N--;
			}
		}
		t0[i]/=N;
	}

	// Loop in events
	for (Long64_t i = 0;i<N; i++){
		//FIXME
		if (i%1000==0) std::cout<<(double)i/N*100<<"%..."<<std::endl;
		c->GetEntry(i);
		for(int iwire = 0; iwire<66; iwire++){
			// reset
			tdcNhitwire = tdcNhit[iwire];
			o_nHits[iwire] = tdcNhitwire;
			o_iHits[iwire]=-1;
			o_time[iwire] = -1e9;
			o_height[iwire] = -1e9;
			o_area[iwire] = -1e9;

			// Get location
			int cht = chg2cht(iwire);
			int chp2 = chg2chp2(iwire);
			int chp3 = chg2chp3(iwire);
			int bid = get_bid(iwire);

			// Get height and etc
			int o_height_all[66][32];
			double o_time_all[66][32];
			double o_area_all[66][32];
			o_pedestal[iwire]=0;
			int clk;
			for(clk = 0; clk<(tdcNhitwire<=0?32:clockNumberDriftTime[iwire][0]-1); clk++){
				o_pedestal[iwire]+=adc[iwire][clk];
			}
			o_pedestal[iwire]/=clk;
			o_pedestalN[iwire] = clk;
			o_pedestalChi2[iwire]=0;
			for(clk = 0; clk<(tdcNhitwire<=0?32:clockNumberDriftTime[iwire][0]); clk++){
				o_pedestalChi2[iwire]+=pow(adc[iwire][clk]-o_pedestal[iwire],2);
			}
			o_areaall[iwire] = 0;
			for ( int ihit = 0; ihit<tdcNhitwire; ihit++){
				o_time_all[iwire][ihit]=driftTime[iwire][ihit]-t0[bid];
				o_height_all[iwire][ihit] = 0;
				for(int clk = clockNumberDriftTime[iwire][ihit]; clk<(ihit+1>=tdcNhitwire?32:clockNumberDriftTime[iwire][ihit+1]); clk++){
					if (adc[iwire][clk]<o_height_all[iwire][ihit]){
						break;
					}
					o_height_all[iwire][ihit]=adc[iwire][clk];
				}

				o_area_all[iwire][ihit] = 0;
				for(int clk = clockNumberDriftTime[iwire][ihit]; clk<(ihit+1>=tdcNhitwire?32:clockNumberDriftTime[iwire][ihit+1]); clk++){
					if (clk!=clockNumberDriftTime[iwire][ihit]&&adc[iwire][clk]<o_pedestal[iwire]) break;
					o_area_all[iwire][ihit] += adc[iwire][clk]-o_pedestal[iwire];
				}
				o_areaall[iwire] += o_area_all[iwire][ihit];
			}

			// find the peak
			for ( int ihit = 0; ihit<tdcNhit[iwire]; ihit++){
				if (o_height_all[iwire][ihit]>Hmin[iwire]){
					o_iHits[iwire] = ihit;
					o_height[iwire] = o_height_all[iwire][ihit];
					o_area[iwire] = o_area_all[iwire][ihit];
					o_time[iwire] = o_time_all[iwire][ihit];
					break;
				}
			}

			// Fill histograms
			if (cht>=0){
				n_tr[cht]->Fill(tdcNhit[iwire]);
			}
			else if(chp3>=0){
				n_p3[chp3]->Fill(tdcNhit[iwire]);
			}
			else if(chp2>=0){
				n_p2[chp2]->Fill(tdcNhit[iwire]);
			}
			int nGoodPeaks = 0;
			for ( int iHit = 0; iHit<tdcNhit[iwire]; iHit++){
				if (cht>=0){
					if (iHit==0) hdt_tr_fp[cht]->Fill(o_time_all[iwire][iHit],o_height_all[iwire][iHit]);
					hdt_tr[cht]->Fill(o_time_all[iwire][iHit],o_height_all[iwire][iHit]);
				}
				else if (chp3>=0){
					if (iHit==0) hdt_p3_fp[chp3]->Fill(o_time_all[iwire][iHit],o_height_all[iwire][iHit]);
					hdt_p3[chp3]->Fill(o_time_all[iwire][iHit],o_height_all[iwire][iHit]);
				}
				else if (chp2>=0){
					if (iHit==0) hdt_p2_fp[chp2]->Fill(o_time_all[iwire][iHit],o_height_all[iwire][iHit]);
					hdt_p2[chp2]->Fill(o_time_all[iwire][iHit],o_height_all[iwire][iHit]);
				}
				if (o_height_all[iwire][iHit]>Hmin[iwire]) nGoodPeaks++;
			}
			if (cht>=0){
				n_tr_ac[cht]->Fill(nGoodPeaks);
			}
			else if (chp2>=0){
				n_p2_ac[chp2]->Fill(nGoodPeaks);
			}
			else if (chp3>=0){
				n_p3_ac[chp3]->Fill(nGoodPeaks);
			}

		}
		t->Fill();
	}
	t->Write();

	// Draw h VS dt for trackers, all peaks
	canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.975,"Peak Height VS Time, (All Peaks)");
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
	for (int i = 0; i<4; i++){
		for (int j = 0; j<3; j++){
			int index = i*3+j;
			pad[index]->cd();
			hdt_tr[index]->Draw("COLZ");
			double hmin = Hmin[cht2chg(index)];
			line_ht[index]->SetX1(-100);
			line_ht[index]->SetY1(hmin);
			line_ht[index]->SetX2(900);
			line_ht[index]->SetY2(hmin);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(800,hmin+5,Form("%lf",hmin));
			text_ht[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.ht.all.tracker.pdf",runNo));
	canvas->SaveAs(Form("run%d.ht.all.tracker.png",runNo));

	// Draw h VS dt for trackers, first peaks
	canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.975,"Peak Height VS Time, (Only 1st Peak)");
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
	for (int i = 0; i<4; i++){
		for (int j = 0; j<3; j++){
			int index = i*3+j;
			pad[index]->cd();
			hdt_tr_fp[index]->Draw("COLZ");
			line_ht[index]->Draw("SAME");
			text_ht[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.ht.fp.tracker.pdf",runNo));
	canvas->SaveAs(Form("run%d.ht.fp.tracker.png",runNo));

	// Draw n for trackers
	canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.975,"Number of Hits");
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
	for (int i = 0; i<4; i++){
		for (int j = 0; j<3; j++){
			int index = i*3+j;
			pad[index]->cd();
			if(n_tr[index]->GetMaximum()<n_tr_ac[index]->GetMaximum())
				n_tr[index]->GetYaxis()->SetRangeUser(0,n_tr_ac[index]->GetMaximum()*1.2);
			n_tr[index]->Draw();
			n_tr_ac[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.n.tracker.pdf",runNo));
	canvas->SaveAs(Form("run%d.n.tracker.png",runNo));

	// Draw h VS dt for p3, all peaks
	canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.975,"Peak Height VS Time in Prototype 3, (All Peaks)");
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
			hdt_p3[index]->Draw("COLZ");
			double hmin = Hmin[chp32chg(index)];
			line_ht[index]->SetX1(-100);
			line_ht[index]->SetY1(hmin);
			line_ht[index]->SetX2(900);
			line_ht[index]->SetY2(hmin);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(800,hmin+5,Form("%lf",hmin));
			text_ht[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.ht.all.p3.pdf",runNo));
	canvas->SaveAs(Form("run%d.ht.all.p3.png",runNo));

	// Draw h VS dt for p3, first peaks
	canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.975,"Peak Height VS Time in Prototype 3, (Only 1st Peak)");
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
			hdt_p3_fp[index]->Draw("COLZ");
			line_ht[index]->Draw("SAME");
			text_ht[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.ht.fp.p3.pdf",runNo));
	canvas->SaveAs(Form("run%d.ht.fp.p3.png",runNo));

	// Draw h VS dt for p2, all peaks
	canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.975,"Peak Height VS Time in Prototype 2, (All Peaks)");
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
			hdt_p2[index]->Draw("COLZ");
			double hmin = Hmin[chp22chg(index)];
			line_ht[index]->SetX1(-100);
			line_ht[index]->SetY1(hmin);
			line_ht[index]->SetX2(900);
			line_ht[index]->SetY2(hmin);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(800,hmin+5,Form("%lf",hmin));
			text_ht[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.ht.all.p2.pdf",runNo));
	canvas->SaveAs(Form("run%d.ht.all.p2.png",runNo));

	// Draw h VS dt for p2, first peaks
	canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.975,"Peak Height VS Time in Prototype 2, (Only 1st Peak)");
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
			hdt_p2_fp[index]->Draw("COLZ");
			line_ht[index]->Draw("SAME");
			text_ht[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.ht.fp.p2.pdf",runNo));
	canvas->SaveAs(Form("run%d.ht.fp.p2.png",runNo));

	// Draw h VS dt for p3, all peaks
	canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.975,"Number of Hits in Prototype 3");
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
			if(n_p3[index]->GetMaximum()<n_p3_ac[index]->GetMaximum())
				n_p3[index]->GetYaxis()->SetRangeUser(0,n_p3_ac[index]->GetMaximum()*1.2);
			n_p3[index]->Draw();
			n_p3_ac[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.n.p3.pdf",runNo));
	canvas->SaveAs(Form("run%d.n.p3.png",runNo));

	// Draw h VS dt for p2, all peaks
	canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.975,"Number of Hits in Prototype 2");
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
			if(n_p2[index]->GetMaximum()<n_p2_ac[index]->GetMaximum())
				n_p2[index]->GetYaxis()->SetRangeUser(0,n_p2_ac[index]->GetMaximum()*1.2);
			n_p2[index]->Draw();
			n_p2_ac[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.n.p2.pdf",runNo));

	// Draw dt
	canvas = new TCanvas("c","c",1024,768);
	text->SetText(0.1,0.975,"t_{0}");
	text->Draw("SAME");
	for (int i = 0; i<2; i++){
		for (int j = 0; j<2; j++){
			int index = i*2+j;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./2*i,0.95/2*(1-j),1./2*(i+1),0.95/2*(2-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./4*i,1./3*(2-j),1./4*(i+1),1./3*(3-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<2; i++){
		for (int j = 0; j<2; j++){
			int index = i*2+j;
			pad[index]->cd();
			bool drawn = false;
			for ( int ch = 0; ch<66; ch++){
				int bid = get_bid_core(ch);
				if (bid==index){
					if (!drawn) {
						h_dt[ch]->GetYaxis()->SetRangeUser(0,t0_nmax[bid]);
						h_dt[ch]->SetTitle(Form("Board #%d",bid));
						h_dt[ch]->Draw();
						drawn = true;
					}
					else h_dt[ch]->Draw("SAME");
				}
			}
			double t_0 = t0[index];
			line_ht[index]->SetX1(t_0);
			line_ht[index]->SetY1(0);
			line_ht[index]->SetX2(t_0);
			line_ht[index]->SetY2(t0_nmax[index]);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(t_0,0,Form("%lf",t_0));
			text_ht[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.dt.pdf",runNo));
	canvas->SaveAs(Form("run%d.dt.png",runNo));

	for (int i = 0; i<66; i++) h_dt[i]->Write();;
	for (int i = 0; i<66; i++) h_h[i]->Write();;
	//for (int i = 0; i<66; i++) h_h2[i]->Write();;
	for (int i = 0; i<66; i++) h_h3[i]->Write();;
	for (int i = 0; i<66; i++) h_h31[i]->Write();;
	for (int i = 0; i<66; i++) h_h32[i]->Write();;
	for (int i = 0; i<12; i++) hdt_tr[i]->Write();;
	for (int i = 0; i<12; i++) hdt_tr_fp[i]->Write();;
	for (int i = 0; i<9; i++) hdt_p2[i]->Write();;
	for (int i = 0; i<9; i++) hdt_p2_fp[i]->Write();;
	for (int i = 0; i<9; i++) hdt_p3[i]->Write();;
	for (int i = 0; i<9; i++) hdt_p3_fp[i]->Write();;
	for (int i = 0; i<12; i++) n_tr[i]->Write();;
	for (int i = 0; i<9; i++) n_p2[i]->Write();;
	for (int i = 0; i<9; i++) n_p3[i]->Write();;
	for (int i = 0; i<12; i++) n_tr_ac[i]->Write();;
	for (int i = 0; i<9; i++) n_p2_ac[i]->Write();;
	for (int i = 0; i<9; i++) n_p3_ac[i]->Write();;

	return 0;
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
	else if (i>=6) chg = i+54;
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
	fprintf(stderr,"\t%s [runNo] <[nEventMax]>\n",prog_name);
}
