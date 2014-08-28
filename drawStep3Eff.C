void drawStep3Eff(){

  cJetEffvsPt_step3_type0->Divide(cJetEffvsPt_step3_type0,cJetEffvsPt_bg_step3,1,1,"B");
  cJetEffvsPt_step3_type1->Divide(cJetEffvsPt_step3_type1,cJetEffvsPt_bg_step3,1,1,"B");
  cJetEffvsPt_step3_type2->Divide(cJetEffvsPt_step3_type2,cJetEffvsPt_bg_step3,1,1,"B");
  cJetEffvsPt_step3_type3->Divide(cJetEffvsPt_step3_type3,cJetEffvsPt_bg_step3,1,1,"B");
  cJetEffvsPt_step3_type4->Divide(cJetEffvsPt_step3_type4,cJetEffvsPt_bg_step3,1,1,"B");

  cJetEffvsPt_step3_type1->SetMarkerColor(2);
  cJetEffvsPt_step3_type2->SetMarkerColor(4);
  cJetEffvsPt_step3_type3->SetMarkerColor(6);
  cJetEffvsPt_step3_type4->SetMarkerColor(kOrange+1);

  cJetEffvsPt_step3_type1->SetLineColor(2);
  cJetEffvsPt_step3_type2->SetLineColor(4);
  cJetEffvsPt_step3_type3->SetLineColor(6);
  cJetEffvsPt_step3_type4->SetLineColor(kOrange+1);

  cJetEffvsPt_step3_type0->SetYTitle("Efficiency");
  cJetEffvsPt_step3_type0->SetXTitle("Jet p_{T}");
  cJetEffvsPt_step3_type0->Draw();
  cJetEffvsPt_step3_type1->Draw("same");
  cJetEffvsPt_step3_type2->Draw("same");
  cJetEffvsPt_step3_type4->Draw("same");

  TLegend *l1 = new TLegend(0.1,0.6,0.5,0.9);
  l1->SetFillColor(0);
  l1->AddEntry(cJetEffvsPt_step3_type0,"D* Efficiency");
  l1->AddEntry(cJetEffvsPt_step3_type1,"D_{0} Efficiency");
  l1->AddEntry(cJetEffvsPt_step3_type2,"D_{s} Efficiency");
  l1->AddEntry(cJetEffvsPt_step3_type4,"D^{+/-} Efficiency");
  l1->Draw("same");
  
}
