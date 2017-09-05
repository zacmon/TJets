void CheckTruthMatch(){

  TFile *f = new TFile("ntuple.root");

  TTree *t = f->Get("JetTree");

  t->Draw("TruthRaw_m>>h0(50,0,250)");
  t->Draw("TruthRaw_m>>h1(50,0,250)","TruthRaw_flavor==0");
  t->Draw("TruthRaw_m>>h2(50,0,250)","TruthRaw_flavor==1");
  t->Draw("TruthRaw_m>>h3(50,0,250)","TruthRaw_flavor==3");

  TH1D *h0 = h0->Clone("h0");
  TH1D *h1 = h1->Clone("h1");
  TH1D *h2 = h2->Clone("h2");
  TH1D *h3 = h3->Clone("h3");

  h0->SetLineColor(1);
  h0->SetLineWidth(2);
  h1->SetLineColor(2);
  h1->SetLineWidth(2);
  h2->SetLineColor(4);
  h2->SetLineWidth(2);
  h3->SetLineColor(7);
  h3->SetLineWidth(2);

  TCanvas *c = new TCanvas("c","c",400,600);
  c->Divide(1,4);
  c->cd(1);
  h0->Draw();
  c->cd(2);
  h1->Draw();
  c->cd(3);
  h2->Draw();
  c->cd(4);
  h3->Draw();






}