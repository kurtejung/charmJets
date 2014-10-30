void makeCorrectedCharmSpec(bool usePDGs=1, bool useUnfolded=0, bool useData=1){

  TFile *unfold = new TFile("pPb_cJet_Unfo_inCM_v21_ak3PF_akPu3PF_noGplus_FirstHalfOnly_Converged_usedParameterizedUnfold0_jtpt35_Inc_clo0_chi100_v8_eta_-2.00To2.00_.root");
  TFile *corrs = new TFile("dCandidateToJetAssns_CJetOnly_RevCuts_JetAssn-10-13.root");
  TFile *data = new TFile("dCandidateToJetAssns_halfData_RevCuts_JetAssn-10-3_fixChgSign.root");
  TFile *purFactors;
  if(useData) purFactors = new TFile("purityHistos_data.root");
  else purFactors = new TFile("purityHistos_MC.root");
  //purFactors = new TFile("purityHistos_data.root");

  //purFactors = new TFile("purityXcheck.root");

  //spit out from drawDMesons.C
  double bgScaleFactors[5] = {0.252154, 0.351251, 0.233429, 0.245132, 0.130555};
  TH1F *unfolded[5];
  TH1F *bgunfold[5];
  TH1D *l1corr[5];
  TH1D *l2corr[5];
  TH1D *l3corr[5];
  TH1D *l4corr[5];
  TH1D *backgroundl1[5];
  TH1D *backgroundl2[5];
  TH1D *backgroundl3[5];
  TH1D *backgroundl4[5];
  TH1D *purity[5];

  //D*, D0, Ds, Ds, Dpm
  string labels[5] = {"D*->C-Jet Spectrum","D_{0}->C-Jet Spectrum","D_{s}->C-Jet Spectrum","D_{s}->C-Jet Spectrum","D^{+/-}->C-Jet Spectrum"};
  double brs[5] = {0.0263,0.0388,0.0288,0.0288,0.0913};
  double ffs[5] = {0.235, 0.549, 0.101, 0.101, 0.232};
  double bgMCscaleFactors[5] = {0.4111,0.3981,0.24969,0.2891,0.17444};
  //double bgMCscaleFactors[5] = {0.45833, 0.4165, 0.2758, 0.2758, 0.16904};
  double bgDatascaleFactors[5] = {0.3128,0.409,0.2822,0.2822,0.16972};
  for(int i=0; i<5; i++){
    if(useUnfolded){
      unfolded[i] = (TH1F*)unfold->Get(Form("hMeas_cent%d",i))->Clone(Form("unfolded_%d",i));
    }
    else if(useData){
      unfolded[i] = (TH1F*)data->Get(Form("CjetSpec_type%d",i))->Clone(Form("unfolded_%d",i));
      bgunfold[i] = (TH1F*)data->Get(Form("bgCjetSpec_type%d",i))->Clone(Form("bgunfolded_%d",i));
    }
    else{
       unfolded[i] = (TH1F*)corrs->Get(Form("CjetSpec_type%d",i))->Clone(Form("unfolded_%d",i));
       bgunfold[i] = (TH1F*)corrs->Get(Form("bgCjetSpec_type%d",i))->Clone(Form("bgunfolded_%d",i));
    }
    l1corr[i] = (TH1D*)corrs->Get(Form("cJetEffvsPt_step1_type%d",i))->Clone(Form("l1corr_%d",i));
    l2corr[i] = (TH1D*)corrs->Get(Form("cJetEffvsPt_step2_type%d",i))->Clone(Form("l2corr_%d",i));
    l3corr[i] = (TH1D*)corrs->Get(Form("cJetEffvsPt_step3_type%d",i))->Clone(Form("l3corr_%d",i));
    l4corr[i] = (TH1D*)corrs->Get(Form("cJetEffvsPt_step4_type%d",i))->Clone(Form("l4corr_%d",i));
    backgroundl1[i] = (TH1D*)(corrs->Get(Form("cJetEffvsPt_bg_step1_type%d",i)))->Clone(Form("backgroundl1_%d",i));
    backgroundl2[i] = (TH1D*)(corrs->Get(Form("cJetEffvsPt_bg_step2_type%d",i)))->Clone(Form("backgroundl2_%d",i));
    backgroundl3[i] = (TH1D*)(corrs->Get(Form("cJetEffvsPt_bg_step3_type%d",i)))->Clone(Form("backgroundl3_%d",i));
    backgroundl4[i] = (TH1D*)(corrs->Get(Form("cJetEffvsPt_bg_step4_type%d",i)))->Clone(Form("backgroundl4_%d",i));

    purity[i] = (TH1D*)(purFactors->Get(Form("purity_type%d",i)))->Clone(Form("purity_type%d",i));
    //purity[i] = (TH1D*)(purFactors->Get(Form("purityXchk_type%d",i)))->Clone(Form("purity_type%d",i));
    
    l1corr[i]->Divide(l1corr[i],backgroundl1[i],1,1,"B");
    l2corr[i]->Divide(l2corr[i],backgroundl2[i],1,1,"B");
    l3corr[i]->Divide(l3corr[i],backgroundl3[i],1,1,"B");
    l4corr[i]->Divide(l4corr[i],backgroundl4[i],1,1,"B");

    l2corr[i]->Multiply(l1corr[i]);
    if(usePDGs){
      l2corr[i]->Scale(brs[i]);
      l2corr[i]->Scale(ffs[i]);
    }
    else{
      cout << "l3corr bins: " << l3corr[i]->GetNbinsX() << endl;
      l2corr[i]->Multiply(l3corr[i]);
      cout << "l4corr bins: " << l4corr[i]->GetNbinsX() << endl;
      l2corr[i]->Multiply(l4corr[i]);
    }
    //if(!useData) bgunfold[i]->Scale(bgMCscaleFactors[i]);
    //else bgunfold[i]->Scale(bgDatascaleFactors[i]);
    //unfolded[i]->Add(bgunfold[i],-1);
    cout << "l1corr bins: "<< l2corr[i]->GetNbinsX() << endl;
    unfolded[i]->Divide(l2corr[i]);
    //bgunfold[i]->Divide(l2corr[i]);
    cout << "purity bins: "<< purity[i]->GetNbinsX() << endl;
    cout << "unfolded bin scheme: ";
    for(int j=1; j<11; j++){
      cout << unfolded[i]->GetBinLowEdge(j) << " ";
    }
    cout << "purity bin scheme: ";
    for(int j=1; j<11; j++){
      cout << purity[i]->GetBinLowEdge(j) << " ";
    }
    cout << endl;
    unfolded[i]->Multiply(purity[i]);
    //bgunfold[i]->Multiply(purity[i]);
  }

  TH1D *tmp = (TH1D*)unfolded[1]->Clone("tmp");
  TH1D *results[5];
  TCanvas *c1 = new TCanvas("c1","",1000,700);
  c1->Divide(2,2);
  for(int i=0; i<5; i++){
    if(i==4) c1->cd(4);
    else c1->cd(i+1);
    if(useUnfolded) unfolded[i]->SetXTitle("Unfolded Jet p_{T}");
    else unfolded[i]->SetXTitle("Jet p_{T}");
    unfolded[i]->SetYTitle("dN/dp_{T}");
    unfolded[i]->SetTitle(labels[i].c_str());
    unfolded[i]->SetMarkerColor(1);
    l1corr[i]->SetMarkerColor(2);
    unfolded[i]->GetXaxis()->SetRangeUser(30,400);
    //unfolded[i]->SetMaximum(1E8);
    //unfolded[i]->SetMinimum(1E2);
    unfolded[i]->Divide(tmp);
    unfolded[i]->Draw();
    bgunfold[i]->SetMarkerColor(2);
    //bgunfold[i]->Draw("same");
    results[i] = (TH1D*)unfolded[i]->Clone(Form("results_%d",i));
    results[i]->Add(bgunfold[i],-1);
    results[i]->SetMarkerColor(kGreen+2);
    //results[i]->Draw("same");
    //l1corr[i]->Draw("same");
  }

}
