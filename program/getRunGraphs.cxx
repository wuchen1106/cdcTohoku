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

	int Hmin[66];
	TFile * ifhmin = new TFile(Form("../info/info.%d."+suffix+"%d.root",startNo,iterationNo));
	TTree * thmin = (TTree*) ifhmin->Get("hmin");
	int idhmin, hhmin;
	thmin->SetBranchAddress("i",&idhmin);
	thmin->SetBranchAddress("h",&hhmin);
	for(int i = 0; i<thmin->GetEntries(); i++){
		thmin->GetEntry(i);
		if (idhmin<0||idhmin>=66){
			printf("WARNING: Entry[%d], wireID = %d\n",i,idhmin);
		}
		Hmin[idhmin] = hhmin;
	}
	ifhmin->Close();

	//===================Get ROOT Files============================
	TChain * c = new TChain("t","t");
	for ( int iFile = startNo; iFile<startNo+nFiles; iFile++){
		c->Add(Form("../root/d_%d.root",iFile));
	}
	int tdcNhit[66];
	int dt[66][NUM_OF_WAVE_SAMPLES];
	int h[66][NUM_OF_WAVE_SAMPLES];
	c->SetBranchAddress("dt",dt);
	c->SetBranchAddress("h",h);
	c->SetBranchAddress("n",tdcNhit);


	TCanvas * canvas;
	TPad * pad[100];
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);

	//===================Prepare Histograms============================
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
		hdt_tr[i] = new TH2D(name,title,1000,-1000,0,200,200,400);
		name = Form("hdt_tr_fp_%d",i);
		hdt_tr_fp[i] = new TH2D(name,title,1000,-1000,0,200,200,400);
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
		hdt_p3[i] = new TH2D(name,title,1000,-1000,0,200,200,400);
		name = Form("hdt_p3_fp_%d",i);
		hdt_p3_fp[i] = new TH2D(name,title,1000,-1000,0,200,200,400);
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
		hdt_p2[i] = new TH2D(name,title,1000,-1000,0,200,200,400);
		name = Form("hdt_p2_fp_%d",i);
		hdt_p2_fp[i] = new TH2D(name,title,1000,-1000,0,200,200,400);
		hdt_p2_fp[i]->GetXaxis()->SetTitle("Time [ns]");
		hdt_p2_fp[i]->GetYaxis()->SetTitle("Peak Height [ADC ch]");
		hdt_p2[i]->GetXaxis()->SetTitle("Time [ns]");
		hdt_p2[i]->GetYaxis()->SetTitle("Peak Height [ADC ch]");
		name = Form("n_p2_%d",i);
		n_p2[i] = new TH1D(name,title,30,0,15);
		n_p2_ac[i] = new TH1D(name,title,30,0,15);
		n_p2_ac[i]->SetLineColor(kRed);
	}

	//===================Loop in Events============================
	Long64_t N = c->GetEntries();
	if (nEventMax&&nEventMax<N) N = nEventMax;
	printf("Processing... %d Events\n",N);
	for (int iEvent = 0;iEvent<N; iEvent++){
		if (iEvent%1000==0) printf("%lf%\n",(double)iEvent/N*100);
		c->GetEntry(iEvent);
		for ( int ch = 0; ch < 66; ch++){
			int cht = ch;
			if (ch>=60) cht-=54;
			int nGoodPeaks = 0;
			int chp = -1;
			if (ch<=5||ch>=60){
				n_tr[cht]->Fill(tdcNhit[ch]);
			}
			else if(ch<=32){
				if (ch<=15&&ch>=13) chp = 15-ch;
				else if (ch<=20&&ch>=18) chp = 20-ch+3;
				else if (ch<=26&&ch>=24) chp = 26-ch+6;
				if (chp>=0){
					n_p3[chp]->Fill(tdcNhit[ch]);
				}
			}
			else{
				if (ch<=42&&ch>=40) chp = 42-ch;
				else if (ch<=47&&ch>=45) chp = 47-ch+3;
				else if (ch<=53&&ch>=51) chp = 53-ch+6;
				if (chp>=0){
					n_p2[chp]->Fill(tdcNhit[ch]);
				}
			}
			for ( int iHit = 0; iHit<tdcNhit[ch]; iHit++){
				if (ch<=5||ch>=60){
					if (iHit==0) hdt_tr_fp[cht]->Fill(dt[ch][iHit],h[ch][iHit]);
					hdt_tr[cht]->Fill(dt[ch][iHit],h[ch][iHit]);
				}
				else{
					if (chp>=0){
						if(ch<=32){
							if (iHit==0) hdt_p3_fp[chp]->Fill(dt[ch][iHit],h[ch][iHit]);
							hdt_p3[chp]->Fill(dt[ch][iHit],h[ch][iHit]);
						}
						else{
							if (iHit==0) hdt_p2_fp[chp]->Fill(dt[ch][iHit],h[ch][iHit]);
							hdt_p2[chp]->Fill(dt[ch][iHit],h[ch][iHit]);
						}
					}
				}
				if (h[ch][iHit]>Hmin[ch]) nGoodPeaks++;
			}
			if (ch<=5||ch>=60){
				n_tr_ac[cht]->Fill(nGoodPeaks);
			}
			else if (chp>=0){
				if(ch<=32){
					n_p3_ac[chp]->Fill(nGoodPeaks);
				}
				else{
					n_p2_ac[chp]->Fill(nGoodPeaks);
				}
			}
		}
	}

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
		}
	}
	canvas->SaveAs(Form("run%d.ht.all.tracker.pdf",startNo));
	canvas->SaveAs(Form("run%d.ht.all.tracker.png",startNo));

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
		}
	}
	canvas->SaveAs(Form("run%d.ht.fp.tracker.pdf",startNo));
	canvas->SaveAs(Form("run%d.ht.fp.tracker.png",startNo));

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
	canvas->SaveAs(Form("run%d.n.tracker.pdf",startNo));
	canvas->SaveAs(Form("run%d.n.tracker.png",startNo));

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
		}
	}
	canvas->SaveAs(Form("run%d.ht.all.p3.pdf",startNo));
	canvas->SaveAs(Form("run%d.ht.all.p3.png",startNo));

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
		}
	}
	canvas->SaveAs(Form("run%d.ht.fp.p3.pdf",startNo));
	canvas->SaveAs(Form("run%d.ht.fp.p3.png",startNo));

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
		}
	}
	canvas->SaveAs(Form("run%d.ht.all.p2.pdf",startNo));
	canvas->SaveAs(Form("run%d.ht.all.p2.png",startNo));

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
		}
	}
	canvas->SaveAs(Form("run%d.ht.fp.p2.pdf",startNo));
	canvas->SaveAs(Form("run%d.ht.fp.p2.png",startNo));

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
	canvas->SaveAs(Form("run%d.n.p3.pdf",startNo));
	canvas->SaveAs(Form("run%d.n.p3.png",startNo));

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
	canvas->SaveAs(Form("run%d.n.p2.pdf",startNo));
	canvas->SaveAs(Form("run%d.n.p2.png",startNo));

	return 0;
}
