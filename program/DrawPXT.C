draw(TString runNo = "251"){

	TString Name = "X Tracker 1, Run#"+runNo;

	// Get TTree
	TFile * f1 = new TFile("../root/fit_251_1.yzflip.root");
	TFile * f2 = new TFile("../info/wire-position.root");
	TTree * t1 = (TTree*) f1->Get("t");
	TTree * t2 = (TTree*) f2->Get("t");
	double inX,slX;
	t1->SetBranchAddress("inX",&inX);
	t1->SetBranchAddress("slX",&slX);

	new TCanvas("c","c",1024,768);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);

	TPad * pad[9];
	for (int i = 0; i<3; i++){
		for (int j = 0; j<3; j++){
			int index = i*3+j;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./3*i,1./3*(2-j),1./3*(i+1),1./3*(3-j));
			//printf("pad[%d] = new TPad(%s,%lf,%lf,%lf,%lf);\n",index,Form("p%d_%d",i,j),1./4*i,1./3*(2-j),1./4*(i+1),1./3*(3-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	for (int i = 0; i<3; i++){
		for (int j = 0; j<3; j++){
			int index = i*3+j;
			pad[index]->cd();
			int ch;
			if(i==0){
				ch = 15-j;
			}
			else if(i==1){
				ch = 20-j;
			}
			else if(i==2){
				ch = 26-j;
			}
			t1->Draw(Form("dt[%d]:fd[%d]>>h%d(1000,-20,20,1000,-100,900)",ch,ch,index),Form("h[%d]>260",ch),"COLZ");
		}
	}

	/*
	TPad * pad[12];
	for (int i = 0; i<4; i++){
		for (int j = 0; j<3; j++){
			int index = i*3+j;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./4*i,1./3*(2-j),1./4*(i+1),1./3*(3-j));
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
			int ch;
			if(i<=1){
				ch = j+i*3;
			}
			else {
				ch = j+(i-2)*3+60;
			}
			t1->Draw(Form("dt[%d]:fd[%d]>>h%d(300,-10,10,300,0,300)",ch,ch,index),"","COLZ");
		}
	}
	*/

	/*
			// Draw beam spot
			if (i==0){
				double dx = 4*(1-j);
				t2->SetMarkerStyle(20);
				//t2->Draw(Form("(y1*%lf+y2*%lf)/%lf:(z1*%lf+z2*%lf)/%lf+378>>h1(256,218,338,256,-50,50)",75.4+dx,75.4-dx,75.4*2,75.4+dx,75.4-dx,75.4*2),"wire_id>5&&wire_id<33");
				t2->Draw(Form("(y1*%lf+y2*%lf)/%lf:(z1*%lf+z2*%lf)/%lf+378",75.4+dx,75.4-dx,75.4*2,75.4+dx,75.4-dx,75.4*2),"wire_id>5&&wire_id<33");
				//h1->Draw();
				for(int iEntry = 0; iEntry<20; iEntry++){
					t1->GetEntry(iEntry);
					TLine * l = new TLine((-60-inX)/slX,-60,(60-inX)/slX,60);
					l->SetLineColor(kCyan);
					l->Draw("SAME");
				}
			}
	p1->cd();
	//t->Draw("dt[0]:dt[1]>>h1(300,0,300,800,0,800)","abs(dt[0]-400)<=400&&abs(dt[1]-150)<=150","COLZ");
	t->Draw("dt[0]:dt[1]>>h1(300,0,300,300,0,300)","abs(dt[0]-400)<=400&&abs(dt[1]-150)<=150","COLZ");
	h1->SetTitle(Name);
	h1->SetContour(50);
	h1->GetXaxis()->SetTitle("Drift Time of wire 1 [ns]");
	h1->GetYaxis()->SetTitle("Drift Time of wire 0 [ns]");
	h1->Draw("COLZ");

	p2->cd();
	t->Draw("dd[0]:dd[1]>>h2(128,0,8,128,0,8)","abs(dt[0]-400)<=400&&abs(dt[1]-150)<=150","COLZ");
	h2->SetTitle(Name);
	h2->SetContour(50);
	h2->GetXaxis()->SetTitle("Drift Distance of wire 1 [mm]");
	h2->GetYaxis()->SetTitle("Drift Distance of wire 0 [mm]");
	h2->Draw("COLZ");

	p3->cd();
	t->Draw("a[0]>>h31(500,0,500)","abs(dt[0]-400)<=400");
	t->Draw("a[1]>>h32(500,0,500)","abs(dt[1]-150)<=150");
	t->Draw("a[2]>>h33(500,0,500)","abs(dt[2]-150)<=150");
	h33->SetTitle(Name);
	h33->GetXaxis()->SetTitle("ADC Sum[mm]");
	h31->SetLineColor(kRed);
	h32->SetLineColor(kBlue);
	h33->SetLineColor(kGreen);
	h33->Draw("");
	h32->Draw("SAME");
	h31->Draw("SAME");
	TLegend * l1 = new TLegend(0.7,0.7,0.9,0.9);
	l1->AddEntry(h31,"wire 0");
	l1->AddEntry(h32,"wire 1");
	l1->AddEntry(h33,"wire 2");
	l1->Draw("SAME");

	p4->cd();
	t->Draw("dd[0]+dd[1]>>h41(128,0,16)","abs(dd[0]-4)<=4&&abs(dd[1]-4)<=4");
	h41->SetTitle(Name);
	h41->GetXaxis()->SetTitle("The Sum of the Drift Distance of wire 0 and wire 1 [mm]");
	h41->Draw("");
	*/
}
