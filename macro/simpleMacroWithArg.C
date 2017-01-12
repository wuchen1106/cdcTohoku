void simpleMacroWithArg(Int_t wireID) {
  gROOT->Reset();
  tree->Draw(Form("q[%d]:driftTime[%d]",wireID,wireID));
}
