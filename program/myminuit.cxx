#include <iostream>
#include <iomanip>  
#include <stdlib.h>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
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
#include "TVector3.h"

//# define PRINT
# define PRINTR

#define t2xDictNTypes  2
#define t2xDictNPoints 16
#define NUM_OF_WAVE_SAMPLES 32
#define NHITS_TO_CHOOSE 4
#define NHITS_TO_CHOOSET 4
#define MIN_WAVE_CH 190
#define MAX_WAVE_CH 400
#define NBINS  256

// ________About Geometry___________
double U = 8 ; 
double SEP = 18+600+45*2;
double PI = TMath::Pi();
double Hmin[66];
int tch[12]={
	0,
	2,
	1,
	3,
	4,
	5,
	60,
	61,
	62,
	63,
	65,
	64
};

// ________About Hit Pattern___________
Int_t patternX = 0;
Int_t patternY = 0;
Int_t indiceX[NHITS_TO_CHOOSE] = {0};
Int_t indiceY[NHITS_TO_CHOOSE] = {0};
Int_t itestX = 0;
Int_t itestY = 0;

// ________About Cells___________
double t01 = 0;
double t02 = 0;
double t03 = 0;
double t04 = 0;
double f1 = 1;
double f2 = 1;
double dxu = 0;
double dxd = 0;
double dyu = 0;
double dyd = 0;

TF1 * XTL_func[12];
TF1 * XTR_func[12];

// ________About Track___________
bool flag[66];
Double_t Z1[66],Y1[66],X1[66],errord[66];// left hand side
Double_t Z2[66],Y2[66],X2[66]; // right hand side
Double_t dd[66][NUM_OF_WAVE_SAMPLES];
int tdcNhit[66];
int dt[66][NUM_OF_WAVE_SAMPLES];
double a[66][NUM_OF_WAVE_SAMPLES];
int h[66][NUM_OF_WAVE_SAMPLES];
Int_t iHit[66];

// ________About Fitting___________
Int_t fit_channel = 0; // 0: X Tracker; 1: Y Tracker
TMinuit *gMinuit = 0;
Double_t arglist[10];
Int_t ierflg = 0;
Double_t amin,edm,errdef;
Int_t nvpar,nparx,icstat;

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
	for (Int_t i=0;i<NHITS_TO_CHOOSE; i++) {
		if (fit_channel == 0){
			index = indiceX[i];
			dfit = get_dist(Z1[index],Y1[index],sl,in);
		}
		else if (fit_channel == 1){
			index = indiceY[i];
			dfit = get_dist(Z1[index],X1[index],sl,in);
		}
		delta  = (dd[index][iHit[index]]-dfit)/errord[index];
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
	double t1,t2,x1,x2,x;
	int type = -1;
	double t0;
	double tori = t;
	if (iwire<=5) t0 = t01;
	else if (iwire<=32) t0 = t02;
	else if (iwire<=59) t0 = t03;
	else t0 = t04;
	t -= t0;

	if (iwire==0||iwire==3||iwire==60||iwire==63) type = 0;
	else if (iwire<=5||(iwire>=61&&iwire<=65)) type = 1;
	else type = 2;

	if (type == 0) t*=f1;
	else if (type == 1) t*=f2;

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
		// FIXME
		/*
		for(int i=0; i<t2xDictNPoints; i++){
			if (t<=t2xDict[type][i]){
				x1=8./t2xDictNPoints*i;
				x2=x1+8./t2xDictNPoints;
				if (i!=0) t1=t2xDict[type][i-1]; else t1=0;
				t2=t2xDict[type][i];
				if (t2!=t1) x=x1+(x2-x1)*(t-t1)/(t2-t1);
				else std::cout<<"t2=t1="<<t1<<"!!!"<<std::endl;
				//if (x<0) std::cout<<"x["<<t<<"] = "<<x<<std::endl;
				if (x==0&&t>0.1) std::cout<<x1<<","<<x2<<","<<t1<<","<<t2<<std::endl;
				return x;
			}
		}
		*/
		x = fabs(f->Eval(t));
		// FIXME
		//printf("%d:%lf->%lf->%lf\n",iwire,tori,t,x);
		return x;
	}
	else if (type==2){
		x=t*8./0.6;
		return x;
	}
}


//______________________________________________________________________________
void getHits(){
	// ...
	int nHitYChs = 0;
	int nHitXChs = 0;
	//for (int ch = 0; ch<66; ch++){
	//	int side = 1;
	//	if (ch == 2||ch==5||ch==62||ch==65) side = -1;
	//	dd[ch][iHit[ch]] = t2x(dt[ch][iHit[ch]],ch,side);
	//	if (ch<=5||ch>=60){
	//		if (fabs(dd[ch][iHit[ch]]-U/2)<U/2&&h[ch][iHit[ch]]>Hmin[ch]){
	//			flag[ch] = true;
	//			if (ch<=2||ch>=63)
	//				nHitXChs++;
	//			else
	//				nHitYChs++;
	//		}
	//		else flag[ch] = false;
	//	}
	//	else{
	//		if (tdcNhit[ch]>0&&h[ch][iHit[ch]]>Hmin[ch]) flag[ch] = true;
	//		else flag[ch] = false;
	//	}
	//}
	for (int ch = 0; ch<66; ch++){
		iHit[ch] = -1;
		for ( int ihit = 0; ihit<tdcNhit[ch]; ihit++){
			if (h[ch][ihit]>Hmin[ch]){
				iHit[ch] = ihit;
				break;
			}
		}
		if (iHit[ch]<0){
			flag[ch] = false;
			continue;
		}
		int side = 1;
		if (ch == 2||ch==5||ch==62||ch==65) side = -1;
		dd[ch][iHit[ch]] = t2x(dt[ch][iHit[ch]],ch,side);
		if (ch<=5||ch>=60){
			if (fabs(dd[ch][iHit[ch]]-U/2)<U/2){
				flag[ch] = true;
				if (ch<=2||ch>=63)
					nHitXChs++;
				else
					nHitYChs++;
			}
			else flag[ch] = false;
		}
		else{
			flag[ch] = true;
		}
	}
	// X tracker
	if (nHitXChs==4){
		if (flag[tch[0]]&&flag[tch[9]]){
			indiceX[0] = tch[0];
			indiceX[2] = tch[9];
			if (flag[tch[1]]&&!flag[tch[2]]&&flag[tch[10]]&&!flag[tch[11]]){
				patternX = 0;
				indiceX[1] = tch[1];
				indiceX[3] = tch[10];
				dd[tch[0]][iHit[tch[0]]] = t2x(dt[tch[0]][iHit[tch[0]]],tch[0],-1);
				dd[tch[9]][iHit[tch[9]]] = t2x(dt[tch[9]][iHit[tch[9]]],tch[9],-1);
				if (fabs(dd[tch[0]][iHit[tch[0]]]-U/2)>U/2)	patternX = -4;
				if (fabs(dd[tch[9]][iHit[tch[9]]]-U/2)>U/2)	patternX = -4;
			}
			else if (!flag[tch[1]]&&flag[tch[2]]&&!flag[tch[10]]&&flag[tch[11]]){
				patternX = 1;
				indiceX[1] = tch[2];
				indiceX[3] = tch[11];
			}
			else if (!flag[tch[1]]&&flag[tch[2]]&&flag[tch[10]]&&!flag[tch[11]]){
				patternX = 2;
				indiceX[1]= tch[2];
				indiceX[3] = tch[10];
				dd[tch[9]][iHit[tch[9]]] = t2x(dt[tch[9]][iHit[tch[9]]],tch[9],-1);
				if (fabs(dd[tch[9]][iHit[tch[9]]]-U/2)>U/2)	patternX = -4;
			}
			else if (flag[tch[1]]&&!flag[tch[2]]&&!flag[tch[10]]&&flag[tch[11]]){
				patternX = 3;
				indiceX[1]= tch[1];
				indiceX[3] = tch[11];
				dd[tch[0]][iHit[tch[0]]] = t2x(dt[tch[0]][iHit[tch[0]]],tch[0],-1);
				if (fabs(dd[tch[0]][iHit[tch[0]]]-U/2)>U/2)	patternX = -4;
			}
			else{
				patternX = -1;
			}
		}
	}
	else if (nHitXChs>4){
		patternX = -2;
	}
	else{
		patternX = -3;
	}
	// Y tracker
	if (nHitYChs==4){
		if (flag[tch[3]]&&flag[tch[6]]){
			indiceY[0] = tch[3];
			indiceY[2] = tch[6];
			if (flag[tch[4]]&&!flag[tch[5]]&&flag[tch[7]]&&!flag[tch[8]]){
				patternY = 0;
				indiceY[1] = tch[4];
				indiceY[3] = tch[7];
			}
			else if (!flag[tch[4]]&&flag[tch[5]]&&!flag[tch[7]]&&flag[tch[8]]){
				patternY = 1;
				indiceY[1] = tch[5];
				indiceY[3] = tch[8];
				dd[tch[3]][iHit[tch[3]]] = t2x(dt[tch[3]][iHit[tch[3]]],tch[3],-1);
				dd[tch[6]][iHit[tch[6]]] = t2x(dt[tch[6]][iHit[tch[6]]],tch[6],-1);
				if (fabs(dd[tch[3]][iHit[tch[3]]]-U/2)>U/2)	patternY = -4;
				if (fabs(dd[tch[6]][iHit[tch[6]]]-U/2)>U/2)	patternY = -4;
			}
			else if (!flag[tch[4]]&&flag[tch[5]]&&flag[tch[7]]&&!flag[tch[8]]){
				patternY = 2;
				indiceY[1] = tch[5];
				indiceY[3] = tch[7];
				dd[tch[3]][iHit[tch[3]]] = t2x(dt[tch[3]][iHit[tch[3]]],tch[3],-1);
				if (fabs(dd[tch[3]][iHit[tch[3]]]-U/2)>U/2)	patternY = -4;
			}
			else if (flag[tch[4]]&&!flag[tch[5]]&&!flag[tch[7]]&&flag[tch[8]]){
				patternY = 3;
				indiceY[1] = tch[4];
				indiceY[3] = tch[8];
				dd[tch[6]][iHit[tch[6]]] = t2x(dt[tch[6]][iHit[tch[6]]],tch[6],-1);
				if (fabs(dd[tch[6]][iHit[tch[6]]]-U/2)>U/2)	patternY = -4;
			}
			else{
				patternY = -1;
			}
		}
	}
	else if (nHitYChs>4){
		patternY = -2;
	}
	else{
		patternY = -3;
	}
	//if (NHITS_TO_CHOOSE)
	//	printf("%dd, %lf,%lf,%lf,%lf\n",pattern,dd[indice[0]][iHit[indice[0]]],dd[itestX][iHit[itestX]],dd[indice[1]][iHit[indice[1]]],dd[indice[2]][iHit[indice[2]]]);
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

//______________________________________________________________________________
void print_usage(char* prog_name)
{
	fprintf(stderr,"\t%s [startNo] <[nFiles] [suffix] [nEventMax] [iterationNo] [N_Min] [f1] [f2]>\n",prog_name);
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
	int iterationNo = 0;
	if (argc>=6) iterationNo = (int)strtol(argv[5],NULL,10);
	int N_MIN = 0;
	if (argc>=7) N_MIN = (int)strtol(argv[6],NULL,10);
	if (argc>=8){
		f1 = strtof(argv[7],NULL);
	}
	if (argc>=9){
		f2 = strtof(argv[8],NULL);
	}

	//===================Set Electronics============================
	//===================Set Cells============================
	// get para
	TChain * c_para = new TChain("para","para");
	std::stringstream buf;
	buf.str("");
	buf.clear();
	buf<<"../info/info."<<startNo<<"."<<suffix<<iterationNo<<".root";
	std::cout<<"Adding \""<<buf.str()<<"\""<<std::endl;
	// get Hmin
	TChain* thmin = new TChain("hmin","hmin");
	thmin->Add(buf.str().c_str());
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
	c_para->Add(buf.str().c_str());
	c_para->SetBranchAddress("t01",&t01);
	c_para->SetBranchAddress("t02",&t02);
	c_para->SetBranchAddress("t03",&t03);
	c_para->SetBranchAddress("t04",&t04);
	c_para->SetBranchAddress("dxu",&dxu);
	c_para->SetBranchAddress("dxd",&dxd);
	c_para->SetBranchAddress("dyu",&dyu);
	c_para->SetBranchAddress("dyd",&dyd);
	c_para->GetEntry(0);
	// get xt
	TChain * c_xt = new TChain("xt","xt");
	c_xt->Add(buf.str().c_str());
	double xt_t, xt_d, xt_sigma;
	int xt_i, xt_n;
	c_xt->SetBranchAddress("t",&xt_t);
	c_xt->SetBranchAddress("d",&xt_d);
	c_xt->SetBranchAddress("i",&xt_i);
	c_xt->SetBranchAddress("n",&xt_n);
	c_xt->SetBranchAddress("sig",&xt_sigma);
	std::vector<std::vector<double> > v_xt_dr;
	std::vector<std::vector<double> > v_xt_dl;
	std::vector<std::vector<double> > v_xt_tr;
	std::vector<std::vector<double> > v_xt_tl;
	v_xt_tl.resize(12);
	v_xt_tr.resize(12);
	v_xt_dr.resize(12);
	v_xt_dl.resize(12);
	for ( int i = 0; i<c_xt->GetEntries(); i++){
		c_xt->GetEntry(i);
		//FIXME
		if (xt_n<N_MIN) continue;
		double sigmai;
		int index = -1;
		if (xt_i<=5||xt_i>=60){
			sigmai = 2.5;
			index = xt_i;
			if (xt_i>=60) index-=54;
		}
		else{
			sigmai = 10;
			continue;
		}
		if (xt_sigma>sigmai*3) continue;
		if (xt_d<=0){
			v_xt_tl[index].push_back(xt_t);
			v_xt_dl[index].push_back(xt_d);
		}
		if (xt_d>=0){
			v_xt_tr[index].push_back(xt_t);
			v_xt_dr[index].push_back(xt_d);
		}
	}
	for ( int i = 0 ; i<12; i++){
		v_xt_dl[i].push_back(0); v_xt_tl[i].push_back(0);
		v_xt_dr[i].push_back(0); v_xt_tr[i].push_back(0);
		TGraph * gxtl = new TGraph(v_xt_dl[i].size(),&v_xt_tl[i][0],&v_xt_dl[i][0]);
		XTL_func[i] = new TF1(Form("xtl_%d",i),"pol5",0,250);
		gxtl->Fit(Form("xtl_%d",i),"N0","");
		gxtl->Print();
		TGraph * gxtr = new TGraph(v_xt_dr[i].size(),&v_xt_tr[i][0],&v_xt_dr[i][0]);
		XTR_func[i] = new TF1(Form("xtr_%d",i),"pol5",0,250);
		gxtr->Fit(Form("xtr_%d",i),"N0","");
#ifdef PRINTR
		TCanvas * ca_temp = new TCanvas("ca_temp","ca_temp");
		gPad->SetGridx(1);
		gPad->SetGridy(1);
		gxtl->SetMarkerStyle(20);
		gxtl->SetMarkerSize(0.8);
		gxtl->SetTitle(Form("tmax = %f",XTL_func[i]->GetX(-U)));
		gxtl->Draw("AP");
		XTL_func[i]->Draw("SAME");
		ca_temp->SaveAs(Form("L%d.pdf",i));
		gxtr->Draw("AP");
		gxtr->SetMarkerStyle(20);
		gxtr->SetMarkerSize(0.8);
		gxtr->SetTitle(Form("tmax = %f",XTR_func[i]->GetX(U)));
		XTR_func[i]->Draw("SAME");
		ca_temp->SaveAs(Form("R%d.pdf",i));
#endif
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
	if_geom->Close();
	// The errors on z values
	Float_t error = 0.2;
	for ( int i = 0; i<66; i++){
		errord[i]=error;
	}

	//===================Get ROOT Files============================
	TChain * c = new TChain("t","t");
	for ( int iFile = startNo; iFile<startNo+nFiles; iFile++){
		buf.str(""); buf.clear();
#ifdef PRINTR
		buf<<"../root/d_"<<iFile<<".root";
#else
		buf<<"../root/d_noadc_"<<iFile<<".root";
#endif
		c->Add(buf.str().c_str());
		std::cout<<"Adding \""<<buf.str()<<"\""<<std::endl;
	}
	int aa[66];
	c->SetBranchAddress("dt",dt);
	c->SetBranchAddress("a",a);
	c->SetBranchAddress("aa",aa);
	c->SetBranchAddress("h",h);
	c->SetBranchAddress("i",iHit);
	c->SetBranchAddress("n",tdcNhit);
#ifdef PRINTR
	int adc[66][NUM_OF_WAVE_SAMPLES];
	int clockNumberDriftTime[66][NUM_OF_WAVE_SAMPLES];
	c->SetBranchAddress("adc",adc);
	c->SetBranchAddress("clockNumberDriftTime",clockNumberDriftTime);
#endif

	//===================Output file============================
	buf.str(""); buf.clear();
	buf<<"../root/fit_"<<startNo<<"_"<<nFiles<<"."<<suffix<<iterationNo<<".root";
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
	Double_t o_fitd[66];
	Double_t o_time[66];
	Int_t o_peak[66];
	Double_t o_sum[66];
	Double_t o_q[66];

	t->Branch("iEvent",&iEvent);
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
	t->Branch("dd",&o_dist,"dd[66]/D");
	t->Branch("fd",&o_fitd,"fd[66]/D");
	t->Branch("dt",&o_time,"dt[66]/D");
	t->Branch("h",&o_peak,"h[66]/I");
	t->Branch("a",&o_sum,"a[66]/D");
	t->Branch("aa",&o_q,"aa[66]/D");
#ifdef PRINTR
	t->Branch("adc",&adc,"adc[66][32]/I");
	t->Branch("tdc",&dt,"tdc[66][32]/I");
	t->Branch("clockNumberDriftTime",&clockNumberDriftTime,"clockNumberDriftTime[66][32]/I");
	t->Branch("n",&tdcNhit,"n[66]/I");
#endif

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
		//printf("In Event(%d)\n",iEvent);
		//===================Get Hits============================
		getHits();
		//FIXME
		//printf("%d %d\n",patternX,patternY);
		if (patternX<0||patternY<0) continue;

		//===================Set Initial============================
		iniX = dd[indiceX[0]][iHit[indiceX[0]]];
		if (patternX==0||patternX==3) iniX*=-1;
		if (patternX<=1){
			sliX = (dd[indiceX[2]][iHit[indiceX[2]]]-dd[indiceX[0]][iHit[indiceX[0]]])/(Z1[indiceX[2]]-Z1[indiceX[0]]);
			if (patternX==0) sliX*=-1;
		}
		else{
			sliX = (dd[indiceX[2]][iHit[indiceX[2]]]+dd[indiceX[0]][iHit[indiceX[0]]])/(Z1[indiceX[2]]-Z1[indiceX[0]]);
			if (patternX==2) sliX*=-1;
		}
		iniY = dd[indiceY[0]][iHit[indiceY[0]]];
		if (patternY==0||patternY==3) iniY*=-1;
		if (patternY<=1){
			sliY = (dd[indiceY[2]][iHit[indiceY[2]]]-dd[indiceY[0]][iHit[indiceY[0]]])/(Z1[indiceY[2]]-Z1[indiceY[0]]);
			if (patternY==0) sliY*=-1;
		}
		else{
			sliY = (dd[indiceY[2]][iHit[indiceY[2]]]+dd[indiceY[0]][iHit[indiceY[0]]])/(Z1[indiceY[2]]-Z1[indiceY[0]]);
			if (patternY==2) sliY*=-1;
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
			o_fitd[ich] = dist;
			o_q[ich] = aa[ich];
			// FIXME
//			if (flag[ich]){
			if (iHit[ich]>=0){
				o_dist[ich] = dd[ich][iHit[ich]];
				o_time[ich] = dt[ich][iHit[ich]];
				o_peak[ich] = h[ich][iHit[ich]];
				o_sum[ich] = a[ich][iHit[ich]];
			}
			else{
				o_dist[ich] = -1e9;
				o_time[ich] = -1e9;
				o_peak[ich] = -1e9;
				o_sum[ich] = -1e9;
			}
		}

#ifdef PRINT
		//===================Draw X Tracker============================
		if (patternX>=0){
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
			e1->SetX1(Z1[indiceX[0]]);e1->SetY1(Y1[indiceX[0]]);e1->SetR1(dd[indiceX[0]][iHit[indiceX[0]]]);e1->SetR2(dd[indiceX[0]][iHit[indiceX[0]]]);
			e2->SetX1(Z1[indiceX[1]]);e2->SetY1(Y1[indiceX[1]]);e2->SetR1(dd[indiceX[1]][iHit[indiceX[1]]]);e2->SetR2(dd[indiceX[1]][iHit[indiceX[1]]]);
			e3->SetX1(Z1[indiceX[2]]);e3->SetY1(Y1[indiceX[2]]);e3->SetR1(dd[indiceX[2]][iHit[indiceX[2]]]);e3->SetR2(dd[indiceX[2]][iHit[indiceX[2]]]);
			e4->SetX1(Z1[indiceX[3]]);e4->SetY1(Y1[indiceX[3]]);e4->SetR1(dd[indiceX[3]][iHit[indiceX[3]]]);e4->SetR2(dd[indiceX[3]][iHit[indiceX[3]]]);
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
					getchi2(chisq,h1->GetXaxis()->GetBinCenter(i),h1->GetYaxis()->GetBinCenter(j));
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
			buf<<"initial: slope = "<<sliX<<", intercept = "<<iniX<<", chi2 = "<<chi2iX;
			text1->SetText(-0.035,13,buf.str().c_str());
			text1->SetTextSize(0.025);
			text1->Draw("SAME");
			buf.str("");buf.clear();
			buf<<"after fit: slope = "<<slX<<", intercept = "<<inX<<", chi2 = "<<chi2X;
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
					if (h[index][j]>Hmin[index]) markerTDC[i][j]->SetMarkerColor(kRed);
					else markerTDC[i][j]->SetMarkerColor(kBlue);
					markerTDC[i][j]->Draw();
				}
			}

			buf.str("");buf.clear();
			buf<<"chi2.X."<<startNo<<"_"<<nFiles<<"."<<iEvent<<"."<<suffix<<"png";
			ca->SaveAs(buf.str().c_str());
			buf.str("");buf.clear();
			buf<<"chi2.X."<<startNo<<"_"<<nFiles<<"."<<iEvent<<"."<<suffix<<"pdf";
			ca->SaveAs(buf.str().c_str());
		}
#endif
		t->Fill();
	}
	t->Write();
	f->Close();

	return 0;
}
