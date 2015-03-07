{
	TChain * c = new TChain("t","t");
	c->Add("../root/d_251.root");
	int dt[66][32];
	int h[66][32];
	int n[66];
	c->SetBranchAddress("dt",dt);
	c->SetBranchAddress("h",h);
	c->SetBranchAddress("n",n);

	TFile * f = new TFile("output.root","RECREATE");
	TTree * ot = new TTree("t","t");
	int o_dt,o_t,o_h;
	ot->Branch("dt",&o_dt);
	ot->Branch("t",&o_t);
	ot->Branch("h",&o_h);

	int ch = 19;
	for ( int i = 0; i<1000000; i++){
		c->GetEntry(i);
		if (n[ch]>=2){
			o_dt = dt[ch][1]-dt[ch][0];
			o_t = dt[ch][1];
			o_h = h[ch][1];
			ot->Fill();
		}
	}
	ot->Write();
	f->Close();
}
