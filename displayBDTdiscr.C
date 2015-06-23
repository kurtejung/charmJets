
void displayBDTdiscr(std::string discr="BDTG_BCvL", std::string discr2="BDTG_BvC"){

  TFile *f1 = new TFile("../input/DMesonCJet_QCDJetOnly_pPbMC_ppReco_akPu3PF_convertToJetTree_tightDCuts.root");
  f1->cd();

  double rangeMax = 0.03;
  double rangeMin = -0.03;

  TH2D *hL = new TH2D("hL","",20,rangeMin,rangeMax,20,rangeMin,rangeMax);
  TH2D *hC = new TH2D("hC","",20,rangeMin,rangeMax,20,rangeMin,rangeMax);
  TH2D *hB = new TH2D("hB","",20,rangeMin,rangeMax,20,rangeMin,rangeMax);
  hL->SetTitle("Light Jets");
  hC->SetTitle("Charm Jets");
  hB->SetTitle("Bottom Jets");
  hL->SetXTitle(discr2.c_str());
  hC->SetXTitle(discr2.c_str());
  hB->SetXTitle(discr2.c_str());
  hL->SetYTitle(discr.c_str());
  hC->SetYTitle(discr.c_str());
  hB->SetYTitle(discr.c_str());

  TTree *j = (TTree*)f1->Get("jets");
  
  TCanvas *c1 = new TCanvas("c1","",600,600);
  c1->cd();
  c1->SetLogz();
  j->Draw(Form("%s:%s>>hL",discr.c_str(),discr2.c_str()),"weight*(abs(jteta)<2 && jtpt>80)","colz");
  hL->SetMinimum(hL->GetMaximum()*1e-6);
  TCanvas *c2 = new TCanvas("c2","",600,600);
  c2->cd();
  c2->SetLogz();
  j->Draw(Form("%s:%s>>hB",discr.c_str(),discr2.c_str()),"weight*(abs(jteta)<2 && jtpt>80 && abs(refparton_flavorForB)==5)","colz");
  hB->SetMinimum(hB->GetMaximum()*1e-6);
  TCanvas *c3 = new TCanvas("c3","",600,600);
  c3->cd();
  c3->SetLogz();
  j->Draw(Form("%s:%s>>hC",discr.c_str(),discr2.c_str()),"weight*(abs(jteta)<2 && jtpt>80 && abs(refparton_flavorForB)==4)","colz");
  hC->SetMinimum(hC->GetMaximum()*1e-6);
  hL->SetFillColor(kBlue+2);
  hC->SetFillColor(kGreen+2);
  hB->SetFillColor(kRed+2);

  TH1D *hL_bcl = new TH1D("hL_bcl","",50,rangeMin,rangeMax);
  hL_bcl->SetXTitle("BDTG_BCvL");
  TH1D *hC_bcl = (TH1D*)hL_bcl->Clone("hC_bcl");
  TH1D *hB_bcl = (TH1D*)hL_bcl->Clone("hB_bcl");
  TH1D *hSum_bcl = (TH1D*)hL_bcl->Clone("hSum_bcl");
  TH1D *hL_bc = (TH1D*)hL_bcl->Clone("hL_bc");
  hL_bc->SetXTitle("BDTG_BvC, BDTG_BCvL>0");
  TH1D *hC_bc = (TH1D*)hL_bc->Clone("hC_bc");
  TH1D *hB_bc = (TH1D*)hL_bc->Clone("hB_bc");
  TH1D *hSum_bc = (TH1D*)hL_bc->Clone("hSum_bc");
  TCanvas *c4 = new TCanvas("c4","",600,600);
  c4->cd();
  THStack *hs_BCvL = new THStack("hs_BCvL","");
  j->Project("hSum_bcl","BDTG_BCvL","weight*(abs(jteta)<2 && jtpt>120)");
  j->Project("hL_bcl","BDTG_BCvL","weight*(abs(jteta)<2 && jtpt>120 && abs(refparton_flavorForB)!=5 && abs(refparton_flavorForB)!=4)");
  hL_bcl->SetFillColor(kBlue+2);
  j->Project("hB_bcl","BDTG_BCvL","weight*(abs(jteta)<2 && jtpt>120 && abs(refparton_flavorForB)==5)");
  hB_bcl->SetFillColor(kRed+2);
  j->Project("hC_bcl","BDTG_BCvL","weight*(abs(jteta)<2 && jtpt>120 && abs(refparton_flavorForB)==4)");
  hC_bcl->SetFillColor(kGreen+2);
  hs_BCvL->Add(hB_bcl);
  hs_BCvL->Add(hC_bcl);
  hs_BCvL->Add(hL_bcl);
  hs_BCvL->SetTitle("C vs Light");
  hs_BCvL->Draw("");

  TCanvas *c5 = new TCanvas("c5","",600,600);
  c5->cd();
  THStack *hs_BvC = new THStack("hs_BvC","");
  j->Project("hSum_bc","BDTG_BvC","weight*(abs(jteta)<2 && jtpt>120 && svtxm<2 && BDTG_BCvL>-0.01)");
  j->Project("hL_bc","BDTG_BvC","weight*(abs(jteta)<2 && jtpt>120 && abs(refparton_flavorForB)!=5 && abs(refparton_flavorForB)!=4 && svtxm<2 && BDTG_BCvL>-0.01)");
  hL_bc->SetFillColor(kBlue+2);
  j->Project("hB_bc","BDTG_BvC","weight*(abs(jteta)<2 && jtpt>120 && abs(refparton_flavorForB)==5 && svtxm<2 && BDTG_BCvL>-0.01)");
  hB_bc->SetFillColor(kRed+2);
  j->Project("hC_bc","BDTG_BvC","weight*(abs(jteta)<2 && jtpt>120 && abs(refparton_flavorForB)==4 && svtxm<2 && BDTG_BCvL>-0.01)");
  hC_bc->SetFillColor(kGreen+2);
  hs_BvC->Add(hB_bc);
  hs_BvC->Add(hC_bc);
  hs_BvC->Add(hL_bc);
  hs_BvC->SetTitle("B vs C");
  hs_BvC->Draw("");

  TCanvas *c6 = new TCanvas("c6","",600,600);
  c6->cd();
  TH1D *hCClone = (TH1D*)hC_bcl->Clone("hCClone");
  hCClone->Divide(hSum_bcl);
  hCClone->Draw();

  TCanvas *c7 = new TCanvas("c7","",600,600);
  c7->cd();
  TH1D *hCClone2 = (TH1D*)hC_bc->Clone("hCClone2");
  hCClone2->Divide(hSum_bc);
  hCClone2->Draw();
  
}
