void adc(int irun, int isub1, int isub2, int boardnumber)
{
  //TFile *f = new TFile(Form("root_2014_11_05/run%04d.root", irun));                               
  // TTree* t = f->Get("t");                                                                        
  TChain* t  = new TChain("tree");

  for (int isub=isub1; isub<=isub2; isub++) {
    t->Add(Form("/Users/SamWong/prototype/root/run_%06d_%03d.root", irun, isub));
  }
  TCanvas* c1 = new TCanvas("c1", "adc", 0, 0, 500, 500);
  c1->Divide(6,8);
  int gr_n=0;
  double gr_x[100];
  double gr_y[100];
  for (int i=0; i<48; i++) {
    c1->cd(i+1);
    t->Draw(Form("adc%d", i));
    TH1F* h1 = new TH1F(Form("h1-%d",i),"", 100, 200, 300);
    t->Draw(Form("adc%d>>h1-%d", i,i), Form("headerver==%d",boardnumber),"goff");
    gr_x[i] = i;
    gr_y[i] = h1->GetRMS();
    gr_n++;
  }
  TCanvas* c3 = new TCanvas("c3", "adc_rms", 0, 600, 500, 400);
  TGraph* gr = new TGraph(gr_n, gr_x, gr_y);
  gr->SetMarkerStyle(20);
  gr->Draw("apl");
}

void tdc(int irun, int isub1, int isub2, int boardnumber)
{
  //TFile *f = new TFile(Form("root_2014_11_05/run%04d.root", irun));                               
  //TTree* t = f->Get("t");                                                                         
  TChain* t  = new TChain("tree");
  for (int isub=isub1; isub<=isub2; isub++) {
    t->Add(Form("/Users/SamWong/prototype/root/run_%06d_%03d.root", irun, isub));
  }
  TCanvas* c1 = new TCanvas("c2", "tdc", 500, 0, 500, 500);
  c1->Divide(8,6);

  //      int chlist[48]                                                                            

  //   long constant = 2^17;                                                                        

  int chlist[48] ;

  for(int j = 0 ; j<48 ; j++)
    {
      chlist[j] = j;
    }
  //   int chlist[6] = {1,3,5,25,29,35};                                                          
#if 1
  for (int i=0; i<48; i++) {
    c1->cd(i+1);
    int ch = chlist[i];
    t->Draw(Form("tdc%d-trigger_time >>htdc%d(50, -1000, 0)", ch, ch),Form("tdc%d!=0 && tdc%d<trig\
ger_time", ch, ch),Form("headerver==%d",boardnumber));


    // t->Draw("tdc%d(%)%d - trigger_time>>htdc%d(50, -2000,100)", ch,2^17, ch),                         
  }
}
#endif

void plot(int irun, int isub1, int isub2, int boardnumber)
{
  adc(irun, isub1, isub2,boardnumber);
  tdc(irun, isub1, isub2,boardnumber);
}


