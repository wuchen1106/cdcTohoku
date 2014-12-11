{
  TTree *t = (TTree*)gROOT->FindObject("tree");
  gROOT->LoadMacro("ePlotEvent.C");
  ePlotEvent e(t);
  e.Loop();
}
