{
	
	TFile *f1 = new TFile("input/DMesonCJet_pPbData_ppReco_akPu3PF_convertToJetTree_withLHCbVars_looseDCuts.root");
	TFile *f2 = new TFile("input/DMesonCJet_QCDJetOnly_pPbMC_ppReco_akPu3PF_convertToJetTree_looseDCuts.root");
	TFile *f3 = new TFile("input/DMesonCJet_QCDJetOnly_Merged_addLHCbVars_pPbMC_centReweight_ppReco_akPu3PF.root");
	TFile *f4 = new TFile("input/DMesonCJet_MergedHalf_addLHCbVars_NoJetTrgCut_pPbdata_ppReco_akPu3PF.root");

	TTree *jt1 = (TTree*)f1->Get("jets");
	TTree *jt2 = (TTree*)f2->Get("jets");
	TTree *jt3 = (TTree*)f3->Get("ct");
	TTree *jt4 = (TTree*)f4->Get("ct");

	TH1D *hdata_post = new TH1D("hdata_post","",20,0,3); hdata_post->Sumw2();
	TH1D *hmc_l_post = new TH1D("hmc_l_post","",20,0,3); hmc_l_post->Sumw2();
	TH1D *hdata_pre = new TH1D("hdata_pre","",20,0,3); hdata_pre->Sumw2();
	TH1D *hmc_l_pre = new TH1D("hmc_l_pre","",20,0,3); hmc_l_pre->Sumw2();
	TH1D *hmc_l_precut = new TH1D("hmc_l_precut","",20,0,3); hmc_l_precut->Sumw2();
	TH1D *hmc_c_precut = new TH1D("hmc_c_precut","",20,0,3); hmc_c_precut->Sumw2();
	TH1D *hmc_b_precut = new TH1D("hmc_b_precut","",20,0,3); hmc_b_precut->Sumw2();
	TH1D *hdata_precut = new TH1D("hdata_precut","",20,0,3); hdata_precut->Sumw2();

	TH1D *hmc_l_aftercut = new TH1D("hmc_l_aftercut","",20,0,3); hmc_l_aftercut->Sumw2();
	TH1D *hmc_c_aftercut = new TH1D("hmc_c_aftercut","",20,0,3); hmc_c_aftercut->Sumw2();
	TH1D *hmc_b_aftercut = new TH1D("hmc_b_aftercut","",20,0,3); hmc_b_aftercut->Sumw2();
	TH1D *hdata_aftercut = new TH1D("hdata_aftercut","",20,0,3); hdata_aftercut->Sumw2();


	jt1->Draw("discr_prob>>hdata_post","weight*(abs(jteta)<1.2 && rawpt>30 && jtpt>100 && jtpt<250)");
	jt2->Draw("discr_prob>>hmc_l_post","weight*(abs(jteta)<1.2 && refpt>30 && rawpt>30 && jtpt>30 && jtpt<100 && pthat>250)","same");
	jt4->Draw("discr_prob>>hdata_pre","weight*(abs(jteta)<1.2 && rawpt>30 && jtpt>100 && jtpt<250)");
	jt3->Draw("discr_prob>>hmc_l_pre","weight*(abs(jteta)<1.2 && refpt>30 && rawpt>30 && jtpt>30 && jtpt<100 && pthat>250)","same");

	hmc_l_pre->Scale(hdata_pre->Integral()/hmc_l_pre->Integral());
	hmc_l_post->Scale(hdata_post->Integral()/hmc_l_post->Integral());

	TH1D *cln = hdata_pre->Clone("cln");
	cln->Divide(hmc_l_pre);

	TH1D *cln2 = hdata_post->Clone("cln");
	cln2->Divide(hmc_l_post);

	TCanvas *c2 = new TCanvas("c2","",600,600);
	cln->Draw();
	cln2->SetMarkerColor(2);
	cln2->Draw("same");

	jt3->Project("hmc_l_precut","discr_prob","weight*(abs(jteta)<1.2 && jtpt>100 && jtpt<250 && refpt>30 && rawpt>30 && abs(refparton_flavorForB)!=4 && abs(refparton_flavorForB)!=5)");
	jt3->Project("hmc_c_precut","discr_prob","weight*(abs(jteta)<1.2 && jtpt>100 && jtpt<250 && refpt>30 && rawpt>30 && abs(refparton_flavorForB)==4)");
	jt3->Project("hmc_b_precut","discr_prob","weight*(abs(jteta)<1.2 && jtpt>100 && jtpt<250 && refpt>30 && rawpt>30 && abs(refparton_flavorForB)==5)");

	jt4->Project("hdata_precut","discr_prob","weight*(abs(jteta)<1.2 && jtpt>100 && jtpt<250)");

	jt2->Project("hmc_l_aftercut","discr_prob","weight*(abs(jteta)<1.2 && jtpt>100 && jtpt<250 && refpt>30 && rawpt>30 && abs(refparton_flavorForB)!=4 && abs(refparton_flavorForB)!=5)");
	jt2->Project("hmc_c_aftercut","discr_prob","weight*(abs(jteta)<1.2 && jtpt>100 && jtpt<250 && refpt>30 && rawpt>30 && abs(refparton_flavorForB)==4)");
	jt2->Project("hmc_b_aftercut","discr_prob","weight*(abs(jteta)<1.2 && jtpt>100 && jtpt<250 && refpt>30 && rawpt>30 && abs(refparton_flavorForB)==5)");

	jt1->Project("hdata_aftercut","discr_prob","weight*(abs(jteta)<1.2 && jtpt>100 && jtpt<250)");

	hmc_l_aftercut->Scale(hmc_l_precut->Integral()/hmc_l_aftercut->Integral());
	hmc_c_aftercut->Scale(hmc_c_precut->Integral()/hmc_c_aftercut->Integral());
	hmc_b_aftercut->Scale(hmc_b_precut->Integral()/hmc_b_aftercut->Integral());

	//hmc_l_aftercut->Divide(hmc_l_precut);
	//hmc_c_aftercut->Divide(hmc_c_precut);
	//hmc_b_aftercut->Divide(hmc_b_precut);

	TCanvas *c3 = new TCanvas("c3","",600,600);
	hmc_l_aftercut->SetMarkerColor(kBlue+2);
	hmc_c_aftercut->SetMarkerColor(kGreen+2);
	hmc_b_aftercut->SetMarkerColor(kRed+2);

	hmc_l_aftercut->SetTitle("after cuts");
	hmc_l_aftercut->Draw();
	hmc_c_aftercut->Draw("same");
	hmc_b_aftercut->Draw("same");

	hdata_aftercut->Scale(hmc_l_aftercut->Integral()/hdata_aftercut->Integral());
	hdata_aftercut->Draw("same");

	TCanvas *c4 = new TCanvas("c4","",600,600);
	hmc_l_precut->SetMarkerColor(kBlue+2);
	hmc_c_precut->SetMarkerColor(kGreen+2);
	hmc_b_precut->SetMarkerColor(kRed+2);

	hmc_l_precut->SetTitle("before cuts");
	hmc_l_precut->Draw();
	hmc_c_precut->Draw("same");
	hmc_b_precut->Draw("same");

	hdata_precut->Scale(hmc_l_precut->Integral()/hdata_precut->Integral());
	hdata_precut->Draw("same");

}