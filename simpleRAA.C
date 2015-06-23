void simpleRAA(){

  TFile *ftop = new TFile("inputToRAA/raa_pPb_numerator-Cjet_FirstHalfDataset_ssvhe0_svtxcorrLT2p5_etaCM_v2_bin0_100_eta-2.0To2.0.root");
  TFile *fbot = new TFile("inputToRAA/raa_pp_denomForpA-Cjet_etaCM_v2_bin0_100_eta-2.0To2.0.root");
  TH1D *htop = (TH1D*)ftop->Get("hRawBData");
  TH1D *hbot = (TH1D*)fbot->Get("hRawBMC");

  TH1D *hcln = (TH1D*)htop->Clone("hcln");
  
  htop->Divide(hbot);
  htop->SetYTitle("c Jet R_{pA}^{PYTHIA}");
  htop->SetTitle("");
  htop->Draw();

  TLine *l1 = new TLine(40,1,400,1);
  l1->SetLineStyle(2);
  l1->Draw("same");

  TCanvas *c2 = new TCanvas("c2","",600,600);
  c2->cd();
  hcln->Draw();
  hbot->SetMarkerColor(2);
  hbot->Draw("same");
  
  
}

