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
	double Neffi_p2[8]={0};
	double Neffi_p3[8]={0};
	TH1D* reso_p2[8];
	TH1D* reso_p3[8];
	TString name;
	for (int i = 0; i<8; i++){
		name = Form("reso_p2_%d",i);
		reso_p2[i] = new TH1D(name,name,128,-2,2);
		name = Form("reso_p3_%d",i);
		reso_p3[i] = new TH1D(name,name,128,-2,2);
	}
	int tr_eff = 0;

	int ntot = ((TTree*)((new TFile(Form("../root/d_%d.root",startNo)))->Get("t")))->GetEntries();
	int ntot2 = ((TTree*)((new TFile(Form("../root/fit_%d_%d."+suffix+"%d.root",startNo,nFiles,iterationNo)))->Get("t")))->GetEntries();

	//===================Get ROOT Files============================
	TChain * c = new TChain("t","t");
	TString filename;
	for ( int iFile = startNo; iFile<startNo+nFiles; iFile++){
		filename = Form("../root/ana.%d."+suffix+"%d.root",iFile,iterationNo);
		c->Add(filename);
	}
	double dd[66];
	double fd[66];
	double chi2X;
	double chi2Y;
	c->SetBranchAddress("dd",dd);
	c->SetBranchAddress("fd",fd);
	c->SetBranchAddress("chi2X",&chi2X);
	c->SetBranchAddress("chi2Y",&chi2Y);

	//===================Loop in Events============================
	Long64_t N = c->GetEntries();
	if (nEventMax&&nEventMax<N) N = nEventMax;
	//std::cout<<"Processing "<<N<<" Events ..."<<std::endl;
	for (int iEvent = 0;iEvent<N; iEvent++){
	//	if (iEvent%1000==0) std::cout<<(double)iEvent/N*100<<"..."<<std::endl;
		c->GetEntry(iEvent);
		if (chi2X>chi2Xmax) continue;
		tr_eff++;
		int ir3 = (int)fabs(fd[19]);
		double reso3 = fabs(fd[19])-fabs(dd[19]);
		int ir2 = (int)fabs(fd[46]);
		double reso2 = fabs(fd[46])-fabs(dd[46]);
		if (ir3<8){
			Neffi_p3[ir3]+=1;
			reso_p3[ir3]->Fill(reso3);
		}
		if (ir2<8){
			Neffi_p2[ir2]+=1;
			reso_p2[ir2]->Fill(reso2);
		}
	}

	double eff0p6,eff2,center,sigma;
	int binl,binr,binc;
	TF1 *f;
	for(int ir = 0; ir<8; ir++){
		if (reso_p2[ir]->Integral()){
			f  = new TF1(Form("f_%d",ir),"gaus",-2,2);
			f->SetParameters(reso_p2[ir]->GetBinContent(reso_p2[ir]->GetMaximumBin()),0,0.2);
			binc = reso_p2[ir]->GetMaximumBin();
			reso_p2[ir]->Fit(Form("f_%d",ir),"qN0","",reso_p2[ir]->GetBinCenter(binc-8),reso_p2[ir]->GetBinCenter(binc+8));
			center = f->GetParameter(1);
			sigma = f->GetParameter(2);
			binl = reso_p2[ir]->FindBin(center-0.6);
			binr = reso_p2[ir]->FindBin(center+0.6);
			eff0p6 = reso_p2[ir]->Integral(binl,binr)/Neffi_p2[ir];
			binl = reso_p2[ir]->FindBin(center-2);
			binr = reso_p2[ir]->FindBin(center+2);
			eff2 = reso_p2[ir]->Integral(binl,binr)/Neffi_p2[ir];
		}
		else{
			sigma = 1e9;
			eff0p6 = 0;
			eff2 = 0;
		}
		std::cout<<sigma<<" "<<eff0p6<<" "<<eff2<<" ";

		if (reso_p3[ir]->Integral()){
			f->SetParameters(reso_p3[ir]->GetBinContent(reso_p3[ir]->GetMaximumBin()),0,0.2);
			binc = reso_p3[ir]->GetMaximumBin();
			reso_p3[ir]->Fit(Form("f_%d",ir),"qN0","",reso_p3[ir]->GetBinCenter(binc-8),reso_p3[ir]->GetBinCenter(binc+8));
			center = f->GetParameter(1);
			sigma = f->GetParameter(2);
			binl = reso_p3[ir]->FindBin(center-0.6);
			binr = reso_p3[ir]->FindBin(center+0.6);
			eff0p6 = reso_p3[ir]->Integral(binl,binr)/Neffi_p3[ir];
			binl = reso_p3[ir]->FindBin(center-2);
			binr = reso_p3[ir]->FindBin(center+2);
			eff2 = reso_p3[ir]->Integral(binl,binr)/Neffi_p3[ir];
		}
		else{
			sigma = 1e9;
			eff0p6 = 0;
			eff2 = 0;
		}
		std::cout<<sigma<<" "<<eff0p6<<" "<<eff2<<" ";
	}
	std::cout<<ntot2/((double)ntot)<<" "<<tr_eff/((double)ntot2)<<std::endl;

	return 0;
}
