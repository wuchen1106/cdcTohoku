draw(TString runNo = "251"){

	TString Name = "X Tracker 1, Run#"+runNo;

	new TCanvas("c","c",1440,960);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	TPad * p1 = new TPad("p1","p1",0   ,0.5 ,0.5 ,1   );
	TPad * p2 = new TPad("p2","p2",0.5 ,0.5 ,1   ,1   );
	TPad * p3 = new TPad("p3","p3",0   ,0   ,0.5 ,0.5 );
	TPad * p4 = new TPad("p4","p4",0.5 ,0   ,1   ,0.5 );
	p1->Draw();
	p2->Draw();
	p3->Draw();
	p4->Draw();
	p1->SetGridx(1);
	p1->SetGridy(1);
	p2->SetGridx(1);
	p2->SetGridy(1);
	p3->SetGridx(1);
	p3->SetGridy(1);
	p4->SetGridx(1);
	p4->SetGridy(1);

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
}
