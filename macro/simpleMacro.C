{
  gROOT->Reset();
  tree->Draw("q[4]:driftTime[4]");
}
