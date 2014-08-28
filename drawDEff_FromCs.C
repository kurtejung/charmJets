void drawDEff_FromCs(){

  TH1D *foreground[4][5];
  TH1D *background[4][5];
  
  string canvasTitles[5] = {"D-star","DZero","Ds->Pi+Phi","Ds->K*K","Charged D"};
  string titles[4] = {"dR Selection Eff","D-Reco Eff.","D Branching Ratio","Fragmentation Fraction"};

  TFile *f1 = new TFile("dCandidateToJetAssns_CJetOnly_RevCuts_JetAssn-8-11.root");
  f1->cd();
  for(int i=0; i<4; i++){
    for(int j=0; j<5; j++){
      foreground[i][j] = (TH1D*)(f1->Get(Form("cJetEffvsPt_step%d_type%d",i+1,j)))->Clone(Form("cJetEffvsPt_step%d_type%d",i+1,j));
      background[i][j] = (TH1D*)(f1->Get(Form("cJetEffvsPt_bg_step%d_type%d",i+1,j)))->Clone(Form("cJetEffvsPt_step%d_type%d",i+1,j)); 

      foreground[i][j]->Divide(foreground[i][j],background[i][j],1,1,"B");
    }
  }

  TCanvas *c1[5];
  for(int i=0; i<5; i++){
    c1[i] = new TCanvas(Form("c1Dash_%d",i),"",1200,800);
    c1[i]->SetTitle(canvasTitles[i].c_str());
    c1[i]->Divide(2,2);
    for(int j=0; j<4; j++){
      c1[i]->cd(j+1);
      foreground[j][i]->SetYTitle(titles[j].c_str());
      foreground[j][i]->SetXTitle("Jet p_{T} (GeV/c)");
      foreground[j][i]->Draw();
    }
  }
}
