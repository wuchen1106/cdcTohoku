{
   TFile * f = new TFile("../eff.root");
   for (int i = 220; i<=226; i++ ){
      TString name1 = Form("h1_%d",i);
      TString name2 = Form("h2_%d",i);
      TH1D * h1 = (TH1D*) f->Get(name1);
      TH1D * h2 = (TH1D*) f->Get(name2);
      if (!h1||!h2){
         printf("Cannot find %d\n",i);
         continue;
      }
      double tx1c_1 = h1->GetBinContent(1,3);
      double ty1c_1 = h1->GetBinContent(3,3);
      double p3s_1 = h1->GetBinContent(6,2);
      double p3c_1 = h1->GetBinContent(7,3);
      double p2s_1 = h1->GetBinContent(12,2);
      double p2c_1 = h1->GetBinContent(13,3);
      double ty2c_1 = h1->GetBinContent(17,3);
      double tx2c_1 = h1->GetBinContent(19,3);
      double tx1c_2 = h2->GetBinContent(1,3);
      double ty1c_2 = h2->GetBinContent(3,3);
      double p3s_2 = h2->GetBinContent(6,2);
      double p3c_2 = h2->GetBinContent(7,3);
      double p2s_2 = h2->GetBinContent(12,2);
      double p2c_2 = h2->GetBinContent(13,3);
      double ty2c_2 = h2->GetBinContent(17,3);
      double tx2c_2 = h2->GetBinContent(19,3);
      printf("%d %4.4lf %4.4lf %4.4lf %4.4lf %4.4lf %4.4lf %4.4lf %4.4lf %4.4lf %4.4lf %4.4lf %4.4lf %4.4lf %4.4lf %4.4lf %4.4lf\n",i,tx1c_1,ty1c_1,p3s_1,p3c_1,p2s_1,p2c_1,ty2c_1,tx2c_1,tx1c_2,ty1c_2,p3s_2,p3c_2,p2s_2,p2c_2,ty2c_2,tx2c_2);
   }
}
