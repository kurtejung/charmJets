#include <iostream>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TLegend.h"

using namespace std;

void drawDEff(){

  TFile *inf[6];
  inf[0] = new TFile("dCandidateToJetAssns_Dstar_RevCuts_JetAssn-8-11.root");
  inf[1] = new TFile("dCandidateToJetAssns_D02Kpi_RevCuts_JetAssn-8-11.root");
  inf[2] = new TFile("dCandidateToJetAssns_Ds2phipi_RevCuts_JetAssn-8-11.root");
  inf[3] = new TFile("dCandidateToJetAssns_Ds2KstarK_RevCuts_JetAssn-8-11.root");
  inf[4] = new TFile("dCandidateToJetAssns_Dpm_RevCuts_JetAssn-8-11.root");
  inf[5] = new TFile("dCandidateToJetAssns_CJetOnly_RevCuts_JetAssn-8-11.root");
  
  string Dlabels[6] = {"D* ","D_{0} ","D_{s}->#pi#phi ","D_{s}->K*K ","D^{+/-} ","C-Jet Only "};

  TH1D *eff[3][6];
  TH1D *eff_fg[3][6];
  TH1D *eff_bg[3][6];
  TH1D *eff_reject[3][6];
  TH1D *cSpec[6];
  TH1D *bgcSpec[6];

  cout << "checkpoint 0" << endl;
  
  for(int i=0; i<6; i++){
    for(int j=0; j<2; j++){
      if(!inf[i]){
	cout << "file " << i << " not found!! " << endl;
	exit(0);
      }
      if(j<2){
	eff[j][i] = (TH1D*)(inf[i]->Get(Form("cJetEffvsPt_step%d",j+1)))->Clone(Form("cJetEffvsPt_type%d_step%d",i,j));
	eff_fg[j][i] = (TH1D*)(inf[i]->Get(Form("cJetEffvsPt_step%d",j+1)))->Clone(Form("cJetEffvsPt2_type%d_step%d",i,j));
	eff_bg[j][i] = (TH1D*)(inf[i]->Get(Form("cJetEffvsPt_bg_step%d",j+1)))->Clone(Form("cJetEffvsPt_bg_type%d_step%d",i,j));
	eff_reject[j][i] = (TH1D*)(inf[i]->Get(Form("cJetEffvsPt_rej_step%d",j+1)))->Clone(Form("cJetEffvsPt_reject_type%d_step%d",i,j));
      }
    }
    if(i<5){
      eff[2][i] = (TH1D*)(inf[5]->Get(Form("cJetEffvsPt_step3_type%d",i)))->Clone(Form("cJetEffvsPt_type%d_step2",i));
      eff_fg[2][i] = (TH1D*)(inf[5]->Get(Form("cJetEffvsPt_step3_type%d",i)))->Clone(Form("cJetEffvsPt2_type%d_step2",i));
      eff_bg[2][i] = (TH1D*)(inf[5]->Get("cJetEffvsPt_bg_step3"))->Clone(Form("cJetEffvsPt_bg_type%d_step2",i));
      eff_reject[2][i] = (TH1D*)(inf[5]->Get("cJetEffvsPt_rej_step3"))->Clone(Form("cJetEffvsPt_reject_type%d_step2",i));

      cSpec[i] = (TH1D*)(inf[i]->Get(Form("CjetSpec_type%d",i)))->Clone(Form("CjetSpec_type%d",i));
      bgcSpec[i] = (TH1D*)(inf[i]->Get(Form("bgCjetSpec_type%d",i)))->Clone(Form("bgCjetSpec_type%d",i));
    }
    cout << "checkpoint 0.5" << endl;
    if(!eff[2][i]) cout << "oh no! " << 2 << " " << i << endl;
    if(!eff_fg[2][i]) cout << "oh no! " << 2 << " " << i << endl;
    if(!eff_bg[2][i]) cout << "oh no! " << 2 << " " << i << endl;
    if(!eff_reject[2][i]) cout << "oh no! " << 2 << " " << i << endl;
    
  }

  cout << "checkpoint 1" << endl;

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->Divide(2,3);
  for(int i=0; i<5; i++){
    c1->cd(i+1);
    for(int j=0;j<3;j++){
      eff[j][i]->Divide(eff[j][i],eff_bg[j][i],1,1,"B");
      eff[j][i]->SetMinimum(0.7);
    }
    eff[0][i]->Draw();
  }
   cout << "checkpoint 2" << endl;

  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->Divide(2,3);
  for(int i=0; i<6; i++){
    c2->cd(i+1);
    TPad *pad1 = new TPad(Form("pad1_%d",i),"",0.0,0.46,1.0,1.0);
    TPad *pad2 = new TPad(Form("pad2_%d",i),"",0.0,0.0,1.0,0.6);
    pad1->Draw();
    pad1->cd();
    gPad->SetLogy();
    eff_bg[0][i]->SetMinimum(0.1);
    string temp = Dlabels[i]+"+ C-Jet";
    eff_bg[0][i]->SetYTitle(temp.c_str());
    eff_bg[0][i]->Draw();
    //eff_fg[0][i]->Multiply(eff_bg[0][i]);
    eff_fg[0][i]->SetLineColor(2);
    eff_fg[0][i]->SetMarkerColor(2);
    eff_fg[0][i]->Draw("same");
    //eff_reject[0][i]->Add(eff_fg[0][i],-1);
    eff_reject[0][i]->SetMarkerColor(4);
    eff_reject[0][i]->SetLineColor(4);
    eff_reject[0][i]->Draw("same");
    if(i==0){
      TLegend *l1 = new TLegend(0.2,0.5,0.2,0.5);
      l1->AddEntry(eff_bg[0][i],"Total C-Jets");
      l1->AddEntry(eff_fg[0][i],"Matched C-Jets");
      l1->AddEntry(eff_reject[0][i],"Unmatched C-Jets");
      l1->Draw();
    }
    pad1->Update();
    eff_bg[0][i]->GetXaxis()->SetLabelSize(0.1);
    eff_bg[0][i]->GetYaxis()->SetLabelSize(0.1);
    eff_bg[0][i]->GetYaxis()->SetTitleOffset(0.5);
    eff_bg[0][i]->GetYaxis()->SetTitleSize(0.1);
    eff[0][i]->GetXaxis()->SetLabelSize(0.1);
    eff[0][i]->GetYaxis()->SetLabelSize(0.1);
    eff[0][i]->GetYaxis()->SetTitleOffset(0.7);
    eff[0][i]->GetXaxis()->SetTitleOffset(1.05);
    eff[0][i]->GetYaxis()->SetTitleSize(0.09);
    eff[0][i]->GetXaxis()->SetTitleSize(0.09);
    
    eff[0][i]->SetNdivisions(6,"Y");

    c2->cd(i+1);
    pad2->Draw();
    pad2->cd();
    eff[0][i]->SetMinimum(0.11);
    eff[0][i]->SetMaximum(1.09);
    eff[0][i]->SetYTitle("Eff.");
    eff[0][i]->Draw();
    pad2->Update();
  }
  TCanvas *c3 = new TCanvas("c3","",1000,600);
  c3->Divide(2,3);
  for(int i=0; i<6; i++){
    c3->cd(i+1);
    eff[1][i]->SetMinimum(0.11);
    eff[1][i]->SetMaximum(1.09);
    string temp = Dlabels[i]+"reco Efficiency";
    eff[1][i]->SetYTitle(temp.c_str());
    eff[1][i]->SetXTitle("Jet p_{T} (GeV/c)");
    eff[1][i]->Draw();
  }
  
  //World average FF (C->D) from http://arxiv.org/pdf/hep-ex/0408149v1.pdf
  //"D* ","D_{0} ","D_{s}->#pi#phi ","D_{s}->K*K ","D^{+/-} "
  double cbranchingRatios[5] = {0.235,0.549,0.101,0.101,0.232};

  //from July 2012 PDG booklet (D->Hadronic decay probabilities)
  double D2HadronBRs[5] = {0.62*0.0388, 0.0388, 0.0549, 0.054, 0.0913};  

  TCanvas *c4 = new TCanvas("c4","",1000,600);
  c4->Divide(2,3);
  for(int i=0; i<5; i++){
    c4->cd(i+1);
    string temp = Dlabels[i]+" Spectrum";
    cSpec[i]->SetYTitle(temp.c_str());
    //cSpec[i]->Add(bgcSpec[i],-1);
    cSpec[i]->Divide(eff[0][i]);
    cSpec[i]->Divide(eff[1][i]);
    //cSpec[i]->Divide(eff[2][i]);
    cSpec[i]->Scale(1./cbranchingRatios[i]);
    cSpec[i]->Scale(1./D2HadronBRs[i]);
    cSpec[i]->Draw();
  }

}
