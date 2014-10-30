#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"

void drawBasicDs(){

  TChain *cc = new TChain("ct","");
  cc->AddFile("input/DMesonCJet_NoJetTrg_AllData_RevCuts_FakeJet2Cut_JetAssn-9-17_part1.root");
  cc->AddFile("input/DMesonCJet_NoJetTrg_AllData_RevCuts_FakeJet2Cut_JetAssn-9-17_part2.root");

  TH1D *hMass = new TH1D("hMass","",50,1.775,1.95);
  hMass->Sumw2();
  hMass->SetXTitle("D_{0} Mass (GeV/c^{2})");
  hMass->SetYTitle("Counts");
  TH1D *hMassStar = new TH1D("hMassStar","",50,0.14,0.155);
  hMassStar->Sumw2();
  hMassStar->SetXTitle("D*-D_{0} Mass (GeV/c^{2})");
  hMassStar->SetYTitle("Counts");

  cc->Project("hMass","dCandMass","(dCandType==2 && dCandPt>30)");
  cc->Project("hMassStar","dCandMass-dCandChildMass","dCandType==1 && dCandPt>30");
  TF1* fit_fun = new TF1("fit_fun", "[0]*(1/[2]/sqrt(6.28)*exp(-0.5*pow((x-[1])/[2], 2))) + pol2(3) + [6]*(1/[8]/sqrt(6.28)*exp(-0.5*pow((x-[7])/[8], 2)))", 1.775, 1.95);
  TF1* fit_fun_star = new TF1("fit_fun_star", "[0]*(1/[2]/sqrt(6.28)*exp(-0.5*pow((x-[1])/[2], 2))) + pol2(3) + [6]*(1/[8]/sqrt(6.28)*exp(-0.5*pow((x-[7])/[8], 2)))", 0.14, 0.155);
  
  float zero = 0;
  float tot = hMass->Integral();
  double MesonMass[2] = {0.1455, 1.8648};
  float var_mean = 0.01, var_width = 0.01;
  float p0 = tot, p1 = MesonMass[1], p2 = 0.015;
  fit_fun->SetParameter(0, p0);
  fit_fun->SetParameter(1, p1);
  fit_fun->SetParameter(2, p2);
  fit_fun->SetParLimits(1, TMath::Max(zero, (float)(MesonMass[1]-var_mean)), MesonMass[1]+var_mean);
  fit_fun->SetParLimits(2, TMath::Max(zero, (float)(p2-var_width)),  p2+var_width);
  fit_fun->FixParameter(6, 0);
  fit_fun->FixParameter(7, 0);
  fit_fun->FixParameter(8, 0);

  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->cd();
  hMass->Fit("fit_fun","q", "", 1.775, 1.95);
  
  TF1* fit_fun_bg = (TF1*)fit_fun->Clone("fit_fun_bg");
  fit_fun_bg->SetParameter(0, 0);
  fit_fun_bg->SetParameter(1, 0);
  fit_fun_bg->SetParameter(2, 0);
  fit_fun_bg->SetParameter(6, 0);
  fit_fun_bg->SetParameter(7, 0);
  fit_fun_bg->SetParameter(8, 0);
  
  fit_fun_bg->SetLineColor(8);
  fit_fun_bg->Draw("same");

  TCanvas *c3 = new TCanvas("c3","",800,600);
  c3->cd();
  var_mean = 0.005;
  var_width = 5e-4;
  p2 = 1e-3;
  p0 = hMassStar->Integral(), p1 = MesonMass[0];
  fit_fun_star->SetParameter(0, p0);
  fit_fun_star->SetParameter(1, p1);
  fit_fun_star->SetParameter(2, p2);
  fit_fun_star->SetParLimits(1, TMath::Max(zero, (float)(MesonMass[0]-var_mean)), MesonMass[0]+var_mean);
  fit_fun_star->SetParLimits(2, TMath::Max(zero, (float)(p2-var_width)),  p2+var_width);
  fit_fun_star->SetParameter(6, p0);
  fit_fun_star->SetParameter(7, p1);
  fit_fun_star->SetParameter(8, p2);
  hMassStar->Fit("fit_fun_star","q","",0.14,0.155);

  TF1* fit_fun_bg_star = (TF1*)fit_fun_star->Clone("fit_fun_bg_star");
  fit_fun_bg_star->SetParameter(0, 0);
  fit_fun_bg_star->SetParameter(1, 0);
  fit_fun_bg_star->SetParameter(2, 0);
  fit_fun_bg_star->SetParameter(6, 0);
  fit_fun_bg_star->SetParameter(7, 0);
  fit_fun_bg_star->SetParameter(8, 0);
  
  fit_fun_bg_star->SetLineColor(8);
  fit_fun_bg_star->Draw("same");
  
  TFile *fout = new TFile("basicDHistos.root","recreate");
  fout->cd();
  hMass->Write();
  hMassStar->Write();
  
}
