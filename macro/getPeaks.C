{
	int ich = 0;

	TChain * c = new TChain("tree","tree");
	c->Add("../root/run_000220_built.root");
	int tdcNhit[66];
	int clockNumberDriftTime[66][32];
	int adc[66][32];
	c->SetBranchAddress("tdcNhit",tdcNhit);
	c->SetBranchAddress("clockNumberDriftTime",clockNumberDriftTime);
	c->SetBranchAddress("adc",adc);

   double height;
   TFile * f = new TFile("output.root","RECREATE");
   TTree * t = new TTree("t","t");
   t->Branch("h",&height);

   double maxheight = -1;
	for (Long64_t i = 0;i<c->GetEntries(); i++){
      c->GetEntry(i);
      height = -1;
      for(int ihit=0; ihit<tdcNhit[ich]; ihit++){
         maxheight = -1;
         for(int clk = clockNumberDriftTime[ich][ihit]; clk<(ihit>=31?32:clockNumberDriftTime[ich][ihit+1]); clk++){
            if (adc[ich][clk]>maxheight) maxheight=adc[ich][clk];
         }
         if (maxheight!=-1) t->Fill();
      }
   }
   t->Write();
}
