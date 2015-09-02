#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TPad.h"
#include "TGraphErrors.h"

#include <iostream>

using namespace std;

void drawCJetFrac(){

  bool doSpectrum = 1;
  bool doRatio = 0;

  TFile* fIncl1 = new TFile("pPb_Unfo_inCM_v31_officialMC_ak3PF_akPu3PF_noGplus_FirstHalfOnly_Converged_usedParameterizedUnfold0_jtpt35_Inc_clo0_chi100_v8_eta_-2.00To2.00_.root");
  TFile* fIncl2 = new TFile("pPb_Unfo_inCM_v31_officialMC_Reverse_WithResCorr_ak3PF_akPu3PF_noGplus_SecondHalfOnly_Converged_usedParameterizedUnfold0_jtpt35_Inc_clo0_chi100_v8_eta_-2.00To2.00_.root");

  TH1D *InclSpec = (TH1D*)fIncl1->Get("hRecoRooUnfold0")->Clone("InclSpec");
  TH1D *InclSpecRev = (TH1D*)fIncl2->Get("hRecoRooUnfold0")->Clone("InclSpecRev");
  InclSpec->Add(InclSpecRev);

  cout << "nbins: "<< InclSpec->GetNbinsX();
  cout << "bin scheme:" << endl;
  for(int i=0; i<InclSpec->GetNbinsX(); i++){
    cout << InclSpec->GetBinLowEdge(i+1) << ", ";
  }
  cout << InclSpec->GetBinLowEdge(InclSpec->GetNbinsX())+InclSpec->GetBinWidth(InclSpec->GetNbinsX()) << endl;

  TFile *bj1 = new TFile("pPb_Unfo_inCM_v31_officialMC_ak3PF_akPu3PF_noGplus_FirstHalfOnly_Converged_usedParameterizedUnfold0_jtpt35_cJets_clo0_chi100_v8_eta_-2.00To2.00_.root");
  TFile *bj2 = new TFile("pPb_Unfo_inCM_v31_officialMC_Reverse_WithResCorr_ak3PF_akPu3PF_noGplus_SecondHalfOnly_Converged_usedParameterizedUnfold0_jtpt35_cJets_clo0_chi100_v8_eta_-2.00To2.00_.root");

  TH1D *bSpec = (TH1D*)bj1->Get("hRecoRooUnfold0")->Clone("bSpec");
  TH1D *bSpecRev = (TH1D*)bj2->Get("hRecoRooUnfold0")->Clone("bSpecRev");
  bSpec->Add(bSpecRev);

  double xbinsRebin[12] = {0, 5, 10, 15, 20, 40, 55, 80, 120, 170, 250, 400};
  //double xbinsRebin[12] = {0,5,10,15,20,40,55,70,100,140,250,400};
  bSpec = (TH1D*)bSpec->Rebin(11,bSpec->GetName(),xbinsRebin);
  for(int ibin=0; ibin<11; ibin++){
    //bSpec->SetBinContent(ibin+1,bSpec->GetBinContent(ibin+1)/bSpec->GetBinWidth(ibin+1));
    //bSpec->SetBinError(ibin+1,bSpec->GetBinError(ibin+1)/bSpec->GetBinWidth(ibin+1));
  }
  //bSpec->Scale(35./20.97);

  TCanvas *c2 = new TCanvas("c2","",600,600);
  c2->cd();
  TH1D *bcln = (TH1D*)bSpec->Clone("bSpec");
  TH1D *inclCln = (TH1D*)InclSpec->Clone("inclCln");
  bcln->Draw();
  inclCln->Draw("same");

  TFile *pythia = new TFile("/Users/kjung/bTagTrees/pPb/histos/ppMC_ppReco_akPu3PF_QCDjetTrig_etashift_Fix2Sample_MCWeightFinalWithVz_noTrgSelection_Full.root");
  
  double *tmp = new double[bSpec->GetNbinsX()+1];
  for(int i=1; i<=bSpec->GetNbinsX(); i++){
    tmp[i-1] = bSpec->GetBinLowEdge(i);
  }
  tmp[bSpec->GetNbinsX()] = bSpec->GetBinLowEdge(bSpec->GetNbinsX()) + bSpec->GetBinWidth(bSpec->GetNbinsX());

  TH1D *pBFrac = new TH1D("pBFrac","",bSpec->GetNbinsX(),tmp); pBFrac->Sumw2();
  TH1D *inclForDiv = new TH1D("inclForDiv","",bSpec->GetNbinsX(),tmp); inclForDiv->Sumw2();
  TTree *nt = (TTree*)pythia->Get("nt");
  nt->Draw("refpt>>inclForDiv","weight*(abs(jteta)<2)");
  nt->Draw("refpt>>pBFrac","weight*(abs(jteta)<2 && abs(refparton_flavorForB)==4)");

  //TFile *preUnfpPb = new TFile("output/cJetFitterV3_looseDCuts_fitToJP_officialMC_akPu3PF_FirstHalf_consistentEta_ssvhe20_addlJetCuts_fixBin_bFractionDataTemplate_pPbpp1_gsp0_jetptcut30__nofilter_vsSvtxm_SSVHPat1.7_FixCL0_bin_0_100_eta_-2.00_2.00.root");
  if(!doSpectrum){
    pBFrac->Divide(inclForDiv);
    bSpec->Divide(InclSpec);
    //bSpec = (TH1D*)preUnfpPb->Get("hBFractionData")->Clone("pBFrac");
  }
  else{
    bSpec->Scale(70./(35e9*2.110*0.85*6.9));
  }

  TCanvas *c3= new TCanvas("c3","",600,600);
  c3->cd();
  TPad *p1 = new TPad("p1_1","",0,0.3,1.0,1.0);
  p1->SetBottomMargin(0.01);
  TPad *p2 = new TPad("p2_1","",0,0,1.0,0.3);
  p2->SetTopMargin(0.01);
  p2->SetBottomMargin(0.4);
  if(doSpectrum) p1->Draw();
  if(doSpectrum) p2->Draw();
  //c3->SetTopMargin(0.06993);

  //pBFrac->SetLineStyle(7);
  pBFrac->SetLineWidth(2);
  pBFrac->SetLineColor(4);
  pBFrac->SetMinimum(0);
  pBFrac->SetXTitle("c jet p_{T} [GeV/c]");
  if(!doSpectrum) pBFrac->SetYTitle("c jet Fraction");
  else pBFrac->SetYTitle("#frac{1}{T_{pA}} #frac{d#sigma}{dp_{T}}");
  pBFrac->GetXaxis()->SetRangeUser(0,400);
  if(!doSpectrum) pBFrac->GetYaxis()->SetRangeUser(0,0.2);
  else if(doRatio) bSpec->GetYaxis()->SetRangeUser(0,2.5);
  else{
    pBFrac->GetYaxis()->SetRangeUser(5e-9,1e-3);
  }
  pBFrac->GetXaxis()->CenterTitle(1);
  for(int i=1; i<=pBFrac->GetNbinsX(); i++){
    cout << " bin center " << pBFrac->GetBinCenter(i) << endl;
    if(pBFrac->GetBinCenter(i)<55){
      pBFrac->SetBinContent(i,-1);
      bSpec->SetBinContent(i,-1);
    } 
  }

  cout << "cp0 "<< endl;

  if(doSpectrum) p1->cd();

  if(doRatio)bSpec->Divide(pBFrac);
  if(!doRatio) pBFrac->Draw("");

  cout << "cp1" << endl;

  const int nBins = 11;
  double xbins[nBins], ybins[nBins], xerr[nBins], yerr[nBins];
  double yerrTot[nBins] = {0,0,0,0,0,0.20337, 0.16156, 0.1361, 0.1636, 0.1626, 0.1501};
  double jecUnc[nBins] = {0,0,0,0,0,0.106701, 0.104898, 0.100017, 0.0972753, 0.0955836, 0.099491};
  for(int jj=0; jj<nBins; jj++){
    if(doSpectrum) yerrTot[jj] = sqrt(pow(yerrTot[jj],2)+pow(jecUnc[jj],2));
  }
  for(int i=1; i<=bSpec->GetNbinsX(); i++){
    xbins[i-1] = bSpec->GetBinLowEdge(i)+bSpec->GetBinWidth(i)/2.;
    ybins[i-1] = bSpec->GetBinContent(i);
    xerr[i-1] = bSpec->GetBinWidth(i)*0.495;
    yerr[i-1] = yerrTot[i-1]*ybins[i-1];
  }

  TGraphErrors *systErr = new TGraphErrors(nBins,xbins,ybins,xerr,yerr);
  systErr->SetFillColor(kYellow);
  if(doRatio){
    bSpec->SetXTitle("c jet p_{T}");
    bSpec->SetYTitle("R_{pA}^{PYTHIA}");
    bSpec->Draw("");
  }
  systErr->Draw("2,5,same");
  bSpec->SetMarkerStyle(20);
  bSpec->SetMarkerColor(1);
  pBFrac->SetMarkerColor(1);
  pBFrac->SetMarkerStyle(25);
  pBFrac->SetFillColor(4);//kGreen+3);
  //pBFrac->SetFillStyle(3244);
  pBFrac->Draw("E2,same");
  TH1D *pBFracCln = (TH1D*)pBFrac->Clone("pBFracCln");
  pBFracCln->SetFillStyle(0);
  pBFracCln->Draw("same");
  bSpec->Draw("same");

  TLegend *l1 = new TLegend(0.42,0.7377,0.889,0.8898);
  l1->SetFillColor(0);
  l1->AddEntry(systErr,"pPb data, -2 < #eta_{lab} < 2","fp");
  l1->AddEntry(pBFrac,"PYTHIA Z2 (5.02 TeV)","lp");
  //l1->AddEntry(ppBfrac,"pp data, -2 < #eta < 2","lp");
  //l1->AddEntry(pBFrac_pp,"PYTHIA Z2 (2.76 TeV), -2 < #eta < 2","lp");
  l1->Draw("same");

  TLatex *cmsP = new TLatex(15,0.177,"CMS ");
  cmsP->SetTextFont(62);
  cmsP->SetTextSize(0.0558);
  cmsP->Draw("same");
  TLatex *cmsPre = new TLatex(71,0.177,"Preliminary ");
  cmsPre->SetTextFont(52);
  cmsPre->SetTextSize(0.0558);
  cmsPre->Draw("same");
  TLatex *l5 = new TLatex(239.5,0.204,"35 nb^{-1} (5.02 TeV)");//; PbPb L = 150 #mub^{-1}");
  l5->SetTextFont(43);
  l5->SetTextSize(25);
  l5->Draw("same");

  if(doSpectrum) p2->cd();

  cout << " checkpoint0 " << endl;

  TH1D *hpARatio = (TH1D*)bSpec->Clone("hpARatio");
  hpARatio->SetTitle("");
  hpARatio->SetYTitle("Data/PYTHIA");
  hpARatio->SetXTitle("c jet p_{T} [GeV/c]");
  hpARatio->SetNdivisions(505,"Y");
  hpARatio->SetLabelSize(0.12,"Y");
  hpARatio->SetLabelSize(0.14,"X");
  hpARatio->SetTitleOffset(0.95,"X");
  hpARatio->SetTitleOffset(0.4,"Y");
  hpARatio->SetTitleSize(0.18,"X");
  hpARatio->SetTitleSize(0.15,"Y");
  hpARatio->SetTickLength(0.07,"X");
  hpARatio->SetMaximum(2.2);
  hpARatio->SetMinimum(-0.2);
  hpARatio->Divide(pBFrac);

  cout << "checkpoint1 " << endl;
  double yerr3[nBins], ybins3[nBins];
  for(int i=1; i<=hpARatio->GetNbinsX(); i++){
    if(hpARatio->GetBinLowEdge(i)<50){
      hpARatio->SetBinContent(i,-999);
    } 
    else{
      yerr3[i-1] = hpARatio->GetBinContent(i)*yerrTot[i-1];
      ybins3[i-1] = hpARatio->GetBinContent(i);
    }
  }
  if(doSpectrum) hpARatio->Draw();
  TLine *line1 = new TLine(0,1,400,1);
  line1->SetLineStyle(2);
  if(doSpectrum) line1->Draw("same");
  TGraphErrors *systErrRatiopA = new TGraphErrors(nBins,xbins,ybins3,xerr,yerr3);
  systErrRatiopA->SetFillColor(kYellow);
  if(doSpectrum) systErrRatiopA->Draw("2,5,same");
  if(doSpectrum) hpARatio->Draw("same");

  cout <<"checkpoint2 "<< endl;

  /***************************************************************************/

  //now plot the pp stuff...
  TCanvas *c4 = new TCanvas("c4","",600,600);
  c4->cd();
  TPad *p12 = new TPad("p1_2","",0,0.3,1.0,1.0);
  p12->SetBottomMargin(0.01);
  TPad *p22 = new TPad("p2_2","",0,0,1.0,0.3);
  p22->SetTopMargin(0.01);
  p22->SetBottomMargin(0.4);
  if(doSpectrum) p12->Draw();
  if(doSpectrum) p22->Draw();
  TFile *ppData = new TFile("output/pp2p76_NewFormatV11_ak3PF_fixBin_bFractionMCTemplate_pPbpp1_gsp0_jetptcut30__nofilter_vsSvtxm_SSVHPat1.7FixCL0_bin_0_100_eta_-2_2.root");
  TH1D *ppBfrac, *ppEff;
  if(!doSpectrum){
    ppBfrac = (TH1D*)bj1->Get("hRecoRooUnfold1")->Clone("ppBfrac");
    TH1D *ppInc = (TH1D*)fIncl1->Get("hRecoRooUnfold1")->Clone("ppInc");
    ppBfrac->Divide(ppInc);
  }
  else{
    //ppBfrac = (TH1D*)ppData->Get("hRawBData")->Clone("ppBfrac");
    //ppEff = (TH1D*)ppData->Get("hBEfficiencyMC")->Clone("ppEff");
    //ppBfrac->Divide(ppEff);
    ppBfrac = (TH1D*)bj1->Get("hRecoRooUnfold1")->Clone("ppBfrac");
    ppBfrac->Scale(1e3/(4.81e12*0.85));
    for(int i=1; i<=ppBfrac->GetNbinsX(); i++){
      //ppBfrac->SetBinError(i,ppBfrac->GetBinError(i)/ppBfrac->GetBinWidth(i));
      //ppBfrac->SetBinContent(i,ppBfrac->GetBinContent(i)/ppBfrac->GetBinWidth(i));
    }
  }

  if(doSpectrum) p12->cd();

  //setup axes
  TH1D *haxis = new TH1D("haxis","",1,0,250);
  if(!doSpectrum) haxis->GetYaxis()->SetRangeUser(0,0.2);
  else haxis->GetYaxis()->SetRangeUser(7e-8,1e-3);
  if(!doSpectrum) haxis->SetYTitle("c jet Fraction");
  else haxis->SetYTitle("#frac{d#sigma}{dp_{T}}");
  haxis->SetXTitle("c jet p_{T} [GeV/c]");
  haxis->GetXaxis()->CenterTitle(1);
  haxis->SetLabelSize(0.06,"Y");
  haxis->Draw();

  ppBfrac->SetMarkerColor(1);
  ppBfrac->SetXTitle("c jet p_{T} [GeV/c]");
  if(!doSpectrum) ppBfrac->SetYTitle("c jet Fraction");
  else ppBfrac->SetYTitle("#frac{d#sigma}{dp_{T}}");
  ppBfrac->GetXaxis()->SetRangeUser(0,250);
  if(!doSpectrum) ppBfrac->GetYaxis()->SetRangeUser(0,0.2);
  else ppBfrac->GetYaxis()->SetRangeUser(1e-8,1e-3);
  ppBfrac->GetXaxis()->CenterTitle(1);
  //ppBfrac->Draw("");

  double xbins2[nBins], ybins2[nBins], xerr2[nBins], yerr2[nBins];
  double yerrTot2[nBins] = {0.223221, 0.173437, 0.205388, 0.221439, 0.348013, 1.11841};
  double jecUnc2[nBins] = {0.106701, 0.104898, 0.100017, 0.0972753, 0.0955836, 0.099491};
  double dMesonErr = 0.055;
  for(int jj=0; jj<nBins; jj++){
    //yerrTot2[jj] = sqrt(pow(yerrTot2[jj],2)+pow(dMesonErr,2));
    if(doSpectrum) yerrTot2[jj] = sqrt(pow(yerrTot2[jj],2)+pow(jecUnc2[jj],2));
  }
  int startpoint = 6;
  if(!doSpectrum) startpoint=6;
  for(int i=startpoint; i<=ppBfrac->GetNbinsX(); i++){
    xbins2[i-startpoint] = ppBfrac->GetBinLowEdge(i)+ppBfrac->GetBinWidth(i)/2.;
    ybins2[i-startpoint] = ppBfrac->GetBinContent(i);
    xerr2[i-startpoint] = ppBfrac->GetBinWidth(i)*0.495;
    yerr2[i-startpoint] = yerrTot2[i-6]*ybins2[i-6];
  }
  TGraphErrors *systErr2 = new TGraphErrors(nBins,xbins2,ybins2,xerr2,yerr2);
  systErr2->SetFillColor(kRed-9);
  systErr2->Draw("2,5,same");
  ppBfrac->Draw("same");

  TFile *ppMC = new TFile("input/DMesonCJet_QCDJetOnly_ppMC_centReweight_ppReco_akPu3PF_2100.root");
  TH1D *pBFrac_pp = new TH1D("pBFrac_pp","",bSpec->GetNbinsX(),tmp); pBFrac_pp->Sumw2();
  TH1D *inclForDiv_pp = new TH1D("inclForDiv_pp","",bSpec->GetNbinsX(),tmp); inclForDiv_pp->Sumw2();
  TTree *ntpp = (TTree*)ppMC->Get("ct");
  ntpp->Project("inclForDiv_pp","refpt","weight*(abs(jteta)<2)");
  ntpp->Project("pBFrac_pp","refpt","weight*(abs(jteta)<2 && abs(refparton_flavorForB)==4)");
  if(!doSpectrum) pBFrac_pp->Divide(inclForDiv_pp);
  else{
    for(int i=1; i<=pBFrac_pp->GetNbinsX(); i++){
      if(pBFrac_pp->GetBinCenter(i)<40) pBFrac_pp->SetBinContent(i,-1);
      //pBFrac_pp->SetBinError(i,pBFrac_pp->GetBinError(i)/pBFrac_pp->GetBinWidth(i));
      //pBFrac_pp->SetBinContent(i,pBFrac_pp->GetBinContent(i)/pBFrac_pp->GetBinWidth(i));
    }
  }
  pBFrac_pp->SetMarkerStyle(25);
  pBFrac_pp->SetMarkerColor(4);
  pBFrac_pp->Draw("same");

  for(int i=1; i<=pBFrac->GetNbinsX(); i++){
    cout << "bin " << i << " : " << pBFrac->GetBinContent(i) << " " << bSpec->GetBinContent(i) << endl;
    if(bSpec->GetBinContent(i)>0) cout << "err: "<< 1-pBFrac->GetBinContent(i)/bSpec->GetBinContent(i) << endl;
    if(pBFrac->GetBinCenter(i)<40){
      pBFrac->SetBinContent(i,-1);
      bSpec->SetBinContent(i,-1);
      pBFrac_pp->SetBinContent(i,-1);
      ppBfrac->SetBinContent(i,-1);
    }
  }

  TLegend *l2 = new TLegend(0.42,0.7377,0.889,0.8898);
  l2->SetFillColor(0);
  //l2->AddEntry(systErr,"pPb data, -2 < #eta_{CM} < 2","fp");
  //l2->AddEntry(pBFrac,"PYTHIA Z2 + HIJING (5.02 TeV)","lp");
  l2->AddEntry(systErr2,"pp data, -2 < #eta < 2","pf");
  l2->AddEntry(pBFrac_pp,"PYTHIA Z2 (2.76 TeV)","lp");
  l2->Draw("same");

  //50.85,0.0926233
  TLatex *cmsP2 = new TLatex(15,0.177,"CMS ");
  cmsP2->SetTextFont(62);
  cmsP2->SetTextSize(0.0558);
  cmsP2->Draw("same");
  TLatex *cmsPre2 = new TLatex(71,0.177,"Preliminary ");
  cmsPre2->SetTextFont(52);
  cmsPre2->SetTextSize(0.0558);
  cmsPre2->Draw("same");
  TLatex *l4 = new TLatex(146,0.204,"4.8 pb^{-1} (2.76 TeV)");//; PbPb L = 150 #mub^{-1}");
  l4->SetTextFont(43);
  l4->SetTextSize(25);
  l4->Draw("same");

  gPad->RedrawAxis();

  if(doSpectrum) p22->cd();
  TH1D *hppRatio = (TH1D*)ppBfrac->Clone("hppRatio");
  hppRatio->SetTitle("");
  hppRatio->SetYTitle("Data/PYTHIA");
  hppRatio->SetNdivisions(505,"Y");
  hppRatio->SetLabelSize(0.12,"Y");
  hppRatio->SetLabelSize(0.14,"X");
  hppRatio->SetTitleOffset(0.95,"X");
  hppRatio->SetTitleOffset(0.4,"Y");
  hppRatio->SetTitleSize(0.18,"X");
  hppRatio->SetTitleSize(0.15,"Y");
  hppRatio->SetTickLength(0.07,"X");
  hppRatio->SetMaximum(2.2);
  hppRatio->SetMinimum(-0.2);
  hppRatio->Divide(pBFrac_pp);
  double yerr4[nBins], ybins4[nBins];
  for(int i=1; i<=hppRatio->GetNbinsX(); i++){
    if(hppRatio->GetBinLowEdge(i)<40){
      hppRatio->SetBinContent(i,-999);
    } 
    else{
      yerr4[i-6] = hppRatio->GetBinContent(i)*yerrTot2[i-6];
      ybins4[i-6] = hppRatio->GetBinContent(i);
    }
  }
  if(doSpectrum) hppRatio->Draw();
  TLine *line2 = new TLine(0,1,250,1);
  line2->SetLineStyle(2);
  if(doSpectrum) line2->Draw("same");
  TGraphErrors *systErrRatio = new TGraphErrors(nBins,xbins2,ybins4,xerr2,yerr4);
  systErrRatio->SetFillColor(kRed-9);
  if(doSpectrum) systErrRatio->Draw("2,5,same");
  if(doSpectrum) hppRatio->Draw("same");


}
