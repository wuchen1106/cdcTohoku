{
	int ich = 0;
	int ch2x[66] = {0,1,1,
		            2,3,3,
		            8,7,6,5,4,
		            4,5,6,7,8,9,
		            8,7,6,5,4,
		            4,5,6,7,8,9,
		            8,7,6,5,4,
	                14,13,12,11,10,
	                10,11,12,13,14,15,
	                14,13,12,11,10,
	                10,11,12,13,14,15,
	                14,13,12,11,10,
	                16,17,17,
	                18,19,19
	};
	int ch2y[66] = {2,1,3,
		            2,1,3,
		            4,4,4,4,4,
		            3,3,3,3,3,3,
		            2,2,2,2,2,
		            1,1,1,1,1,1,
		            0,0,0,0,0,
		            4,4,4,4,4,
		            3,3,3,3,3,3,
		            2,2,2,2,2,
		            1,1,1,1,1,1,
		            0,0,0,0,0,
		            2,1,3,
		            2,1,3
	};

	TChain * c = new TChain("tree","tree");
	c->Add("../root/run_000220_built.root");
	int tdcNhit[66];
	int clockNumberDriftTime[66][32];
	int adc[66][32];
	c->SetBranchAddress("tdcNhit",tdcNhit);
	c->SetBranchAddress("clockNumberDriftTime",clockNumberDriftTime);
	c->SetBranchAddress("adc",adc);
	TCanvas * ca = new TCanvas("c","c",1440,960);
	TPad * p1 = new TPad("p1","p1",0,0,0.5,1);
	TPad * p2 = new TPad("p2","p2",0.5,0,1,1);
	for (Long64_t i = 0;(i+1)*50<c->GetEntries(); i++){
		bool notDrawnYet = true;
		for(int ii = 0; ii<100; ii++){
			c->GetEntry(i*100+ii);
	//		for(int j = 0; j<66; j++){
	//			for(int k = 0; k<32; k++){
	//				printf("%d: tdc[%d][%d] = %d; adc[%d][%d] = %d;\n",i,j,k,tdc[j][k],j,k,adc[j][k]);
	//			}
	//		}
			TH1D * h1 = new TH1D("h1","h1",64,-32,32);
			TH1D * h2 = new TH1D("h2","h2",64,-32,32);
			h1->GetYaxis()->SetRangeUser(200,350);
			h2->GetYaxis()->SetRangeUser(200,350);
			if (tdcNhit[ich]<=0) continue;
			int offset = clockNumberDriftTime[ich][0];
			for(int clk=0; clk<32; clk++){
				if (h>240)
					h1->Fill(clk-offset,adc[ich][clk]);
				else
					h2->Fill(clk-offset,adc[ich][clk]);
			}
			if (notDrawnYet){
				p1->cd();
				h1->Draw("LP");
				p2->cd();
				h2->Draw("LP");
				notDrawnYet=false;
			}
			else{
				p1->cd();
				h1->Draw("LPSAME");
				p2->cd();
				h2->Draw("LPSAME");
			}
		}
		ca->WaitPrimitive();
	}
}
