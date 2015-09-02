
#include <iostream>
#include <fstream>
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TH2F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TCanvas.h"

using namespace std;


TH1F *transformReco2Gen(TH1F *h, TH2F *hXform){
  
  if(h->GetNbinsX() != hXform->GetNbinsX()){ cout << "Warning! Matrix binning and histogram binning non-compatible!" << endl; 
    cout << "histo bins: " << h->GetNbinsX() << " " << "matrix bins: "<< hXform->GetNbinsX() << endl;  
    exit(0);
  }

  TH1F *htemp = (TH1F*)h->Clone(h->GetName());
  for(int i=0;i<hXform->GetNbinsX();i++){
    float genMC = 0;
    float genErrMC = 0;
    for(int j=0;j<hXform->GetNbinsY();j++){
      
      float coeff = hXform->GetBinContent(i+1,j+1);
      float recoMC = h->GetBinContent(j+1);
      float recoErrMC = h->GetBinError(j+1);
      genMC += coeff * recoMC;
      genErrMC +=  coeff * recoErrMC * coeff * recoErrMC;
    }
    genErrMC = sqrt(genErrMC);
    htemp->SetBinContent(i+1, genMC);
    //htemp->SetBinError(i+1, genErrMC);
  }
  return htemp;
}

TH1F *zeroErrors(TH1F *h){
  for(int i=0;i<h->GetNbinsX();i++) h->SetBinError(i+1,0);

  return h;
}

void plotSystematicComponents(bool doTransform=0, bool ppPbPb=1, bool plotSymmetrized=1)
{
  gStyle->SetErrorX(0);
  gStyle->SetHistLineWidth(2);
  gROOT->ForceStyle(1);

  TFile *fMatrix = NULL;
  TH2F *xNorm=NULL;

  if(doTransform){
    if(ppPbPb){
      fMatrix= new TFile("output/reco2GenMatrix_pA.root");
      xNorm = (TH2F*)(fMatrix->Get("hRecoVsGenNorm"))->Clone("xNorm");
    }
    else{
      fMatrix = new TFile("output/reco2GenMatrix_pA.root");
      xNorm = (TH2F*)(fMatrix->Get("hRecoVsGenNorm"))->Clone("xNorm");
    }
  }

  TFile *f0, *f1a, *f1b, *f2, *f3, *f4a, *f4b;

  if(ppPbPb){
    // Default parameters and LTJP systematic
    f0 = new TFile("output/cJetFitterV3_looseDCuts_fitToJP_officialMC_akPu3PF_FirstHalf_consistentEta_ssvhe20_addlJetCuts_fixBin_bFractionMCTemplate_pPbpp1_gsp0_jetptcut30__nofilter_vsSvtxm_SSVHPat1.7_FixCL0_bin_0_100_eta_-2.00_2.00.root");
    // Working point variation
    f1a = new TFile("output/cJetFitterV3_looseDCuts_fitToJP_officialMC_akPu3PF_FirstHalf_consistentEta_ssvhe20_addlJetCuts_fixBin_bFractionMCTemplate_pPbpp1_gsp0_jetptcut30__nofilter_vsSvtxm_SSVHPat1.4_FixCL0_bin_0_100_eta_-2.00_2.00.root");
    f1b = new TFile("output/cJetFitterV3_looseDCuts_fitToJP_officialMC_akPu3PF_FirstHalf_consistentEta_ssvhe20_addlJetCuts_fixBin_bFractionMCTemplate_pPbpp1_gsp0_jetptcut30__nofilter_vsSvtxm_SSVHPat2.5_FixCL0_bin_0_100_eta_-2.00_2.00.root");
    // Data-driven template
    f2 = new TFile("output/cJetFitterV3_looseDCuts_fitToJP_officialMC_akPu3PF_FirstHalf_consistentEta_ssvhe20_addlJetCuts_fixBin_bFractionMCTemplate_pPbpp1_gsp0_jetptcut30__nofilter_vsSvtxm_SSVHPat1.7_FixCL0_bin_0_100_eta_-2.00_2.00.root");
    // fix charm
    f3 = new TFile("output/cJetFitterV3_looseDCuts_fitToJP_officialMC_akPu3PF_FirstHalf_consistentEta_ssvhe20_addlJetCuts_fixBin_bFractionMCTemplate_pPbpp1_gsp0_jetptcut30__nofilter_vsSvtxm_SSVHPat1.7_FixCL0_bin_0_100_eta_-2.00_2.00.root");
    //gspUpDown
    f4a = new TFile("output/cJetFitterV3_looseDCuts_fitToJP_officialMC_akPu3PF_FirstHalf_consistentEta_ssvhe20_addlJetCuts_fixBin_bFractionMCTemplate_pPbpp1_gsp2_jetptcut30__nofilter_vsSvtxm_SSVHP_at1.7_FixCL0_bin_0_100_eta_-2.00_2.00.root");
    f4b = new TFile("output/cJetFitterV3_looseDCuts_fitToJP_officialMC_akPu3PF_FirstHalf_consistentEta_ssvhe20_addlJetCuts_fixBin_bFractionMCTemplate_pPbpp1_gsp3_jetptcut30__nofilter_vsSvtxm_SSVHP_at1.7_FixCL0_bin_0_100_eta_-2.00_2.00.root");
  }
  else{
    // Default parameters and LTJP systematic
    f0 = new TFile("output/pp2p76_NewFormatV11_ak3PF_fixBin_bFractionMCTemplate_pPbpp1_gsp0_jetptcut30__nofilter_vsSvtxm_SSVHPat1.7FixCL0_bin_0_100_eta_-2_2.root");
    // Working point variation
    f1a = new TFile("output/pp2p76_NewFormatV11_ak3PF_fixBin_bFractionMCTemplate_pPbpp1_jetptcut30__nofilter_vsSvtxm_SSVHPat1.4FixCL0_bin_0_100_eta_-2_2.root");
    f1b = new TFile("output/pp2p76_NewFormatV11_ak3PF_fixBin_bFractionMCTemplate_pPbpp1_jetptcut30__nofilter_vsSvtxm_SSVHPat2.5FixCL0_bin_0_100_eta_-2_2.root");
    // Data-driven template
    f2 = new TFile("output/pp2p76_NewFormatV11_ak3PF_fixBin_bFractionMCTemplate_pPbpp1_gsp0_jetptcut30__nofilter_vsSvtxm_SSVHPat1.7FixCL0_bin_0_100_eta_-2_2.root");
    // fix charm
    f3 = new TFile("output/pp2p76_NewFormatV11_ak3PF_fixBin_bFractionMCTemplate_pPbpp1_gsp0_jetptcut30__nofilter_vsSvtxm_SSVHPat1.7FixCL0_bin_0_100_eta_-2_2.root");
    //gsp up and down
    f4a = new TFile("output/pp2p76_NewFormatV11_ak3PF_fixBin_bFractionMCTemplate_pPbpp1_gsp2_jetptcut30__nofilter_vsSvtxm_SSVHPat1.7FixCL0_bin_0_100_eta_-2_2.root");
    f4b = new TFile("output/pp2p76_NewFormatV11_ak3PF_fixBin_bFractionMCTemplate_pPbpp1_gsp3_jetptcut30__nofilter_vsSvtxm_SSVHPat1.7FixCL0_bin_0_100_eta_-2_2.root");
    //Eric's Gplus PU filter
    //f5 = new TFile("output/NewFormatV5_bFractionMCTemplate_pppp1_gPlusFilter_SSVHEat2.0FixCL0_bin_0_40_eta_0_2.root");
  }

  double workPointpA[6] = {0.126768,0.060995,0.0985326,0.072506,0.0867271,0.0860764};
  double workPointpp[6] = {0.0616832, 0.0963795, 0.070848, 0.116304, 0.107444, 0.12458};


  // Divide by efficiency

  // LTJP and default
  TH1F *hRawSpec0;
  if(ppPbPb) hRawSpec0 =  (TH1F*)f0->Get("hRawBData");
  else hRawSpec0 =  (TH1F*)f0->Get("hRawBData");
  TH1F *hEff0 =  (TH1F*)f0->Get("hBEfficiencyMC");
  TH1F *hEffLTJP =  (TH1F*)f0->Get("hBEfficiencyDataLTJP");

  TH1F *hDefault = (TH1F*)hRawSpec0->Clone("hDefault");
  TH1F *hLTJP = (TH1F*)hRawSpec0->Clone("hLTJP");
  
  hDefault->Divide(hEff0);
  if(ppPbPb) hLTJP->Divide(hEffLTJP);
  else hLTJP->Divide(hEffLTJP);
  double jpSyst[8] = {1.11,1.05,1.05,1.1,1.1,1.1,1.2,1.2};
  //double ppjpSyst[6] = {1.1,1.01,1.01,1.13,1.06,1.07};
  TF1 *fp2 = new TF1("fp2","[0]+x*[1]+x*x*[2]",40,400);
  fp2->SetParameter(0,0.545537);
  fp2->SetParameter(1,0.00652);
  fp2->SetParameter(2,-2.032e-5);
  for(int i=1; i<=hLTJP->GetNbinsX(); i++){
    if(!ppPbPb) hLTJP->SetBinContent(i,fp2->Eval(hLTJP->GetBinCenter(i))*hDefault->GetBinContent(i));
  }
  
  if(doTransform){
    hDefault = transformReco2Gen(hDefault, xNorm);
    hLTJP = transformReco2Gen(hLTJP, xNorm);
  }


  // Working point variation
  TH1F *hRawSpec1a;
  if(ppPbPb) hRawSpec1a = (TH1F*)f1a->Get("hRawBData");
  else hRawSpec1a =  (TH1F*)f1a->Get("hRawBData");
  TH1F *hEff1a =  (TH1F*)f1a->Get("hBEfficiencyMC");
  
  TH1F *hRawSpec1b;
  if(ppPbPb) hRawSpec1b = (TH1F*)f1b->Get("hRawBData");
  else hRawSpec1b =  (TH1F*)f1b->Get("hRawBData");
  TH1F *hEff1b =  (TH1F*)f1b->Get("hBEfficiencyMC");
  
  TH1F *hWorkingPointUp = (TH1F*)hRawSpec1a->Clone("hWorkingPointUp");
  TH1F *hWorkingPointDown = (TH1F*)hRawSpec1b->Clone("hWorkingPointDown");
  
  hWorkingPointUp->Divide(hEff1a);
  hWorkingPointDown->Divide(hEff1b);

 if(doTransform){
  hWorkingPointUp = transformReco2Gen(hWorkingPointUp, xNorm);
  hWorkingPointDown = transformReco2Gen(hWorkingPointDown, xNorm);
 }

 //if(ppPbPb){
   // Data-driven template
   TH1F *hRawSpec2 =  (TH1F*)f2->Get("hRawBData");
   TH1F *hEff2 =  (TH1F*)f2->Get("hBEfficiencyMC");
   
   TH1F *hDataDriven = (TH1F*)hRawSpec2->Clone("hDataDriven");
   hDataDriven->Divide(hEff2);

   //Changed to D-Meson template uncertainty
  for(int i=1; i<=hDataDriven->GetNbinsX(); i++){
    hDataDriven->SetBinContent(i,1.055*hDefault->GetBinContent(i));
  }


   //double xbins[6] = {40,60,95,120,170,400};
   //hDataDriven = (TH1F*)hDataDriven->Rebin(5,hDataDriven->GetName(),xbins);
   //hDataDriven->Scale(hDefault->Integral()/hDataDriven->Integral());
   
  // if(doTransform) hDataDriven = transformReco2Gen(hDataDriven, xNorm);
   
   // Fix charm - changed to MC template Unc. (toyMC)
   TH1F *hRawSpec3 =  (TH1F*)f3->Get("hRawBData");
   TH1F *hEff3 =  (TH1F*)f3->Get("hBEfficiencyMC");
   
   TH1F *hCharm = (TH1F*)hRawSpec3->Clone("hCharm");
   hCharm->Divide(hEff3);
   double templateUnc[6] = {1.103,1.115,1.054,1.092,1.079,1.10};
   double pptemplateUnc[6] = {1.08, 1.079, 1.184, 1.164, 1.325, 2.0};
  
   if(doTransform) hCharm = transformReco2Gen(hCharm, xNorm);
   for(int ibin=1; ibin<=hCharm->GetNbinsX(); ibin++){
     if(ppPbPb){
       hCharm->SetBinContent(ibin,hCharm->GetBinContent(ibin)*templateUnc[ibin-1]);
       hCharm->SetBinError(ibin,hCharm->GetBinError(ibin)*templateUnc[ibin-1]);
     }
     else{
      hCharm->SetBinContent(ibin,hCharm->GetBinContent(ibin)*pptemplateUnc[ibin-1]);
      hCharm->SetBinError(ibin,hCharm->GetBinError(ibin)*pptemplateUnc[ibin-1]);
    }
  }
 //}
 // gluon systematics
 TH1F *hRawSpec4a;
 if(ppPbPb) hRawSpec4a = (TH1F*)f4a->Get("hRawBData");
 else hRawSpec4a = (TH1F*)f4a->Get("hRawBData");
 TH1F *hEff4a =  (TH1F*)f4a->Get("hBEfficiencyMC");
 
 TH1F *hGlueUp = (TH1F*)hRawSpec4a->Clone("hGlueUp");
 hGlueUp->Divide(hEff4a);
 //double xbins[6] = {40,60,95,120,170,400};
 //hGlueUp = (TH1F*)hGlueUp->Rebin(5,hGlueUp->GetName(),xbins);

 TH1F *hRawSpec4b;
 if(ppPbPb) hRawSpec4b = (TH1F*)f4b->Get("hRawBData");
 else hRawSpec4b = (TH1F*)f4b->Get("hRawBData");
 TH1F *hEff4b =  (TH1F*)f4b->Get("hBEfficiencyMC");
 
 TH1F *hGlueDown = (TH1F*)hRawSpec4b->Clone("hGlueDown");
 hGlueDown->Divide(hEff4b);
 //hGlueDown = (TH1F*)hGlueDown->Rebin(5,hGlueDown->GetName(),xbins);
 
 if(doTransform) hGlueUp = transformReco2Gen(hGlueUp, xNorm);
 if(doTransform) hGlueDown = transformReco2Gen(hGlueDown, xNorm);

 // pu systematics
 /* TH1F *hRawSpec5 =  (TH1F*)f5->Get("hRawBData");
 TH1F *hEff5 =  (TH1F*)f5->Get("hBEfficiencyMC");
 
 TH1F *hPU = hRawSpec5->Clone("hPU");
 hPU->Scale(1./7.54082720175149213e-01); //adjustment in luminosity for PU rejection
 hPU->Divide(hEff5);
 
 if(doTransform) hPU = transformReco2Gen(hPU, xNorm);*/

 //plotting style stuff

 hLTJP->SetMarkerColor(kRed+2);
 hWorkingPointUp->SetMarkerColor(kSpring+2);
 hWorkingPointDown->SetMarkerColor(kGreen+2);
 hDataDriven->SetMarkerColor(kMagenta+2);
 hCharm->SetMarkerColor(kCyan+2);
 hGlueUp->SetMarkerColor(kYellow+2);
 hGlueDown->SetMarkerColor(kOrange+2);
 //hPU->SetMarkerColor(kBlue-1);
 
 hLTJP->SetLineColor(kRed+2);
 hWorkingPointUp->SetLineColor(kSpring+2);
 hWorkingPointDown->SetLineColor(kGreen+2);
 hDataDriven->SetLineColor(kMagenta+2);
 hCharm->SetLineColor(kCyan+2);
 hGlueUp->SetLineColor(kYellow+2);
 hGlueDown->SetLineColor(kOrange+2);
 // hPU->SetLineColor(kBlue-1);

 if(!plotSymmetrized){
   hLTJP->SetMarkerStyle(24);
   hWorkingPointUp->SetMarkerStyle(26);
   hWorkingPointDown->SetMarkerStyle(32);
   hDataDriven->SetMarkerStyle(3);
   hCharm->SetMarkerStyle(4);
   hGlueUp->SetMarkerStyle(2);
   hGlueDown->SetMarkerStyle(5);
   //hPU->SetMarkerStyle(28);

  hLTJP->SetLineStyle(7);
  hWorkingPointUp->SetLineStyle(6);
  hWorkingPointDown->SetLineStyle(8);
  hDataDriven->SetLineStyle(3);
  hCharm->SetLineStyle(4);
  hGlueUp->SetLineStyle(2);
  hGlueDown->SetLineStyle(5);
  //hPU->SetLineStyle(10);
   
  hLTJP->SetLineWidth(3);
  hWorkingPointUp->SetLineWidth(3);
  hWorkingPointDown->SetLineWidth(3);
  hDataDriven->SetLineWidth(3);
  hCharm->SetLineWidth(3);
  hGlueUp->SetLineWidth(3);
  hGlueDown->SetLineWidth(3);
  //hPU->SetLineWidth(3);
  
 }



 hDefault->Draw();
 hLTJP->Draw("same");
 
 hWorkingPointUp->Draw("same");
 hWorkingPointDown->Draw("same");
 hDataDriven->Draw("same");
 hCharm->Draw("same");
 hGlueUp->Draw("same");
 hGlueDown->Draw("same");
 //hPU->Draw("same");
 

  // Now take fractional systematics

  TCanvas *c2=new TCanvas("c2","c2",600,600);

  TH1F *hFracTotal = (TH1F*)hDefault->Clone("hFracTotal");
  TH1F *hFracLTJP = (TH1F*)hLTJP->Clone("hFracLTJP");
  TH1F *hFracWorkingPointUp = (TH1F*)hWorkingPointUp->Clone("hFracWorkingPointUp");
  TH1F *hFracWorkingPointDown = (TH1F*)hWorkingPointDown->Clone("hFracWorkingPointDown");
  TH1F *hFracDataDriven = (TH1F*)hDataDriven->Clone("hFracDataDriven");
  TH1F *hFracCharm = (TH1F*)hCharm->Clone("hFracCharm");
  TH1F *hFracGlueUp = (TH1F*)hGlueUp->Clone("hFracGlueUp");
  TH1F *hFracGlueDown = (TH1F*)hGlueDown->Clone("hFracGlueDown");
  //TH1F *hFracPU = hPU->Clone("hFracPU");
  
  hFracTotal->Reset();
  hFracLTJP->Reset();
  hFracWorkingPointUp->Reset();
  hFracWorkingPointDown->Reset();
  hFracDataDriven->Reset();
  hFracCharm->Reset();
  hFracGlueUp->Reset();
  hFracGlueDown->Reset();
  //hFracPU->Reset();

  TH1F *hFracWorkingPoint = (TH1F*)hWorkingPointUp->Clone("hFracWorkingPoint");
  TH1F *hFracGlue = (TH1F*)hGlueUp->Clone("hFracGlue");
  
  TH1F *hFracTotalInv = (TH1F*)hFracTotal->Clone("hFracTotalInv");
  TH1F *hFracLTJPInv = (TH1F*)hFracLTJP->Clone("hFracLTJPInv");
  TH1F *hFracWorkingPointInv = (TH1F*)hFracWorkingPoint->Clone("hFracWorkingPointInv");
  TH1F *hFracDataDrivenInv = (TH1F*)hFracDataDriven->Clone("hFracDataDrivenInv");
  TH1F *hFracCharmInv = (TH1F*)hFracCharm->Clone("hFracCharmInv");
  TH1F *hFracGlueInv = (TH1F*)hFracGlue->Clone("hFracGlueInv");
  //TH1F *hFracPUInv = (TH1F*)hFracPU->Clone("hFracPUInv");

  for(int i =0;i < hDefault->GetNbinsX(); i++){

    float valDefault = hDefault->GetBinContent(i+1);

    float fracLTJP = hLTJP->GetBinContent(i+1)/valDefault-1.;
    float fracWorkingPointUp = hWorkingPointUp->GetBinContent(i+1)/valDefault-1.;
    float fracWorkingPointDown = hWorkingPointDown->GetBinContent(i+1)/valDefault-1.;
    float fracDataDriven = hDataDriven->GetBinContent(i+1)/valDefault-1.;
    float fracCharm = hCharm->GetBinContent(i+1)/valDefault-1.;
    float fracGlueUp = hGlueUp->GetBinContent(i+1)/valDefault-1.;
    float fracGlueDown = hGlueDown->GetBinContent(i+1)/valDefault-1.;
    //float fracPU = hPU->GetBinContent(i+1)/valDefault-1.;
    
    hFracLTJP->SetBinContent(i+1,fracLTJP);
    hFracWorkingPointUp->SetBinContent(i+1,fracWorkingPointUp);
    hFracWorkingPointDown->SetBinContent(i+1,fracWorkingPointDown);
    hFracDataDriven->SetBinContent(i+1,fracDataDriven);
    hFracCharm->SetBinContent(i+1,fracCharm);
    hFracGlueUp->SetBinContent(i+1,fracGlueUp);
    hFracGlueDown->SetBinContent(i+1,fracGlueDown);
    //hFracPU->SetBinContent(i+1, fracPU);

    float fracWorkingPoint = TMath::Max(fabs(fracWorkingPointUp),fabs(fracWorkingPointDown));
    float fracGlue = TMath::Max(fabs(fracGlueUp),fabs(fracGlueDown));

    //Redefine working point variation to come from values obtained in calcTaggerVariance.C
    if(ppPbPb) fracWorkingPoint = workPointpA[i];
    else fracWorkingPoint = workPointpp[i];

    hFracWorkingPoint->SetBinContent(i+1,fracWorkingPoint);
    hFracGlue->SetBinContent(i+1,fracGlue);

    hFracLTJPInv->SetBinContent(i+1,-1.*fracLTJP);
    hFracWorkingPointInv->SetBinContent(i+1,-1.*fracWorkingPoint);
    hFracDataDrivenInv->SetBinContent(i+1,-1.*fracDataDriven);
    hFracCharmInv->SetBinContent(i+1,-1.*fracCharm);
    hFracGlueInv->SetBinContent(i+1,-1.*fracGlue);
    //hFracPUInv->SetBinContent(i+1, -1.*fracPU);
        

    float fracTotal;
    if(ppPbPb) fracTotal = sqrt(fracLTJP*fracLTJP+
			   fracWorkingPoint*fracWorkingPoint+
			   fracDataDriven*fracDataDriven+
			   fracCharm*fracCharm
				//fracPU*fracPU+
				+fracGlue*fracGlue
			   );
    else fracTotal = sqrt(fracLTJP*fracLTJP+
			   fracWorkingPoint*fracWorkingPoint+
         fracDataDriven*fracDataDriven+
         fracCharm*fracCharm+
			   //fracPU*fracPU+
			   fracGlue*fracGlue
			   );
    
    hFracTotal->SetBinContent(i+1,fracTotal);
    cout << "bin " << i+1 << " sysErr: "<< fracTotal << endl;
    hFracTotalInv->SetBinContent(i+1,-1.*fracTotal);

    // Don't need errors
    /*
    float effDefault = hDefault->GetBinError(i+1);
    float relErrDefault = effDefault/valDefault;

    float relErrLTJP = hLTJP->GetBinError(i+1)/hLTJP->GetBinContent(i+1);
    float relErrWorkingPointUp = hWorkingPointUp->GetBinError(i+1)/hWorkingPointUp->GetBinContent(i+1);
    float relErrWorkingPointDown = hWorkingPointDown->GetBinError(i+1)/hWorkingPointDown->GetBinContent(i+1);
    float relErrDataDriven = hDataDriven->GetBinError(i+1)/hDataDriven->GetBinContent(i+1);
    float relErrCharm = hCharm->GetBinError(i+1)/hCharm->GetBinContent(i+1);


    float errFracLTJP = fracLTJP*sqrt(relErrDefault*relErrDefault+relErrLTJP*relErrLTJP);
    float errFracWorkingPointUp = fracWorkingPointUp*sqrt(relErrDefault*relErrDefault+relErrWorkingPointUp*relErrWorkingPointUp);
    float errFracWorkingPointDown = fracWorkingPointDown*sqrt(relErrDefault*relErrDefault+relErrWorkingPointDown*relErrWorkingPointDown);
    float errFracDataDriven = fracDataDriven*sqrt(relErrDefault*relErrDefault+relErrDataDriven*relErrDataDriven);
    float errFracCharm = fracCharm*sqrt(relErrDefault*relErrDefault+relErrCharm*relErrCharm);

    hFracLTJP->SetBinError(i+1,errFracLTJP);
    hFracWorkingPointUp->SetBinError(i+1,errFracWorkingPointUp);
    hFracWorkingPointDown->SetBinError(i+1,errFracWorkingPointDown);
    hFracDataDriven->SetBinError(i+1,errFracDataDriven);
    hFracCharm->SetBinError(i+1,errFracCharm);
    */
    


  }

  hFracTotal->GetYaxis()->SetTitle("Relative Error");
  hFracTotal->SetTitle("");
  hFracTotal->GetYaxis()->SetTitleOffset(1.3);
  hFracTotal->GetXaxis()->SetNdivisions(505);
  
  hFracTotal->SetMaximum(0.6);
  hFracTotal->SetMinimum(-0.4);

  hFracTotal->SetFillStyle(1001);
  hFracTotalInv->SetFillStyle(1001);
  hFracTotal->SetFillColor(kGray);
  hFracTotalInv->SetFillColor(kGray);
  
  hFracTotal->Draw("h");
  hFracTotalInv->Draw("h,same");

  // hFracLTJP->SetFillColor(hFracLTJP->GetMarkerColor()-2);
  // if(plotSymmetrized)hFracLTJP->SetFillStyle(3003);
  hFracLTJP->Draw("same");

  hFracLTJPInv->SetMaximum(0.5);
  hFracLTJPInv->SetMinimum(-0.5);
  //hFracLTJPInv->SetFillColor(hFracLTJP->GetMarkerColor()-2);
  //if(plotSymmetrized)hFracLTJPInv->SetFillStyle(3003);
  if(plotSymmetrized)hFracLTJPInv->Draw("same");


  // hFracWorkingPointUp->SetFillColor(hFracWorkingPointUp->GetMarkerColor()-2);
  // if(plotSymmetrized)hFracWorkingPointUp->SetFillStyle(3006);
  if(!plotSymmetrized)hFracWorkingPointUp->Draw("same");
  //hFracWorkingPointDown->SetFillColor(hFracWorkingPointDown->GetMarkerColor()-2);
  //if(plotSymmetrized)hFracWorkingPointDown->SetFillStyle(3006);
  if(!plotSymmetrized)hFracWorkingPointDown->Draw("same");

  hFracWorkingPoint = zeroErrors(hFracWorkingPoint);
  hFracWorkingPointInv = zeroErrors(hFracWorkingPointInv);

  hFracWorkingPoint = zeroErrors(hFracWorkingPoint);
  // hFracWorkingPoint->SetFillColor(hFracWorkingPoint->GetMarkerColor()-2);
  // if(plotSymmetrized)hFracWorkingPoint->SetFillStyle(3007);
  if(plotSymmetrized)hFracWorkingPoint->Draw("same");
  // hFracWorkingPointInv->SetFillColor(hFracWorkingPoint->GetMarkerColor()-2);
  // if(plotSymmetrized)hFracWorkingPointInv->SetFillStyle(3007);
  if(plotSymmetrized)hFracWorkingPointInv->Draw("same");



  //hFracDataDriven->SetFillColor(hFracDataDriven->GetMarkerColor()-2);
  // if(plotSymmetrized)hFracDataDriven->SetFillStyle(3004);
  hFracDataDriven->Draw("same");
  // hFracDataDrivenInv->SetFillColor(hFracDataDriven->GetMarkerColor()-2);
  // if(plotSymmetrized)hFracDataDrivenInv->SetFillStyle(3004);
  if(plotSymmetrized)hFracDataDrivenInv->Draw("same");


  // hFracCharm->SetFillColor(hFracCharm->GetMarkerColor()-2);
  // if(plotSymmetrized)hFracCharm->SetFillStyle(3005);
  hFracCharm->Draw("same");
  // hFracCharmInv->SetFillColor(hFracCharm->GetMarkerColor()-2);
  // if(plotSymmetrized)hFracCharmInv->SetFillStyle(3005);
  if(plotSymmetrized)hFracCharmInv->Draw("same");
  
  hFracGlue = zeroErrors(hFracGlue);
  hFracGlueInv = zeroErrors(hFracGlueInv);

  // hFracGlueUp->SetFillColor(hFracGlueUp->GetMarkerColor()-2);
  //if(plotSymmetrized)hFracGlueUp->SetFillStyle(3007);
  if(!plotSymmetrized)hFracGlueUp->Draw("same");
  // hFracGlueDown->SetFillColor(hFracGlue->GetMarkerColor()-2);
  // if(plotSymmetrized)hFracGlueDown->SetFillStyle(3007);
  if(!plotSymmetrized)hFracGlueDown->Draw("same");

  //hFracGlue->SetFillColor(hFracGlue->GetMarkerColor()-2);
  // if(plotSymmetrized)hFracGlue->SetFillStyle(3007);
  if(plotSymmetrized)hFracGlue->Draw("same");
  //hFracGlueInv->SetFillColor(hFracGlue->GetMarkerColor()-2);
  // if(plotSymmetrized)hFracGlueInv->SetFillStyle(3007);
  if(plotSymmetrized)hFracGlueInv->Draw("same");

  // hFracPU->SetFillColor(hFracPU->GetMarkerColor()-2);
  // if(plotSymmetrized)hFracPU->SetFillStyle(3008);
  //MARKED hFracPU->Draw("same");
  // hFracPUInv->SetFillColor(hFracPU->GetMarkerColor()-2);
  // if(plotSymmetrized)hFracPUInv->SetFillStyle(3008);
  //MARKED if(plotSymmetrized)hFracPUInv->Draw("same");

  string legStyle = "l";
  //if(plotSymmetrized) legStyle="f";

  TLegend *leg=new TLegend(0.225,0.6555,0.624,0.905);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hFracTotal,"Total","lf");
  leg->AddEntry(hFracLTJP,"Reference Tagger",legStyle.c_str());
  if(plotSymmetrized)leg->AddEntry(hFracWorkingPoint,"Working point",legStyle.c_str());
  else{
    leg->AddEntry(hFracWorkingPointDown,"Working point down",legStyle.c_str());
    leg->AddEntry(hFracWorkingPointUp,"Working point up",legStyle.c_str());
  }

  //if(ppPbPb) leg->AddEntry(hFracDataDriven,"Data-driven template",legStyle.c_str());
  leg->AddEntry(hFracCharm,"MC template",legStyle.c_str());
  leg->AddEntry(hFracDataDriven, "D-Meson Decay",legStyle.c_str());
  if(plotSymmetrized)leg->AddEntry(hFracGlue,"Gluon Splitting",legStyle.c_str());
  // else{
  //  leg->AddEntry(hFracGlueUp,"Glue up",legStyle.c_str());
  //  leg->AddEntry(hFracGlueDown,"Glue down",legStyle.c_str());
  // }
  //leg->AddEntry(hFracPU,"GPlus PU Filter",legStyle.c_str());
  leg->Draw();

  TLatex *lt1 = new TLatex(183.5,0.455,"#sqrt{s} = 5.02 TeV");
  if(!ppPbPb) lt1->SetText(183.5,0.455,"#sqrt{s} = 2.76 TeV");
  lt1->SetTextSize(0.042);
  lt1->Draw();
  TLatex *lt2 = new TLatex(217,0.52,"CMS Preliminary");
  lt2->SetTextSize(0.05);
  lt2->Draw();

  TFile *fout=NULL;
  if(ppPbPb)fout = new TFile("systematics_components_pPb.root","recreate");
  else fout = new TFile("systematics_components_pp.root","recreate");
  hFracLTJP->Write();
  hFracWorkingPointUp->Write();
  hFracWorkingPointDown->Write();
  hFracWorkingPoint->Write();
  if(ppPbPb) hFracDataDriven->Write();
  if(ppPbPb) hFracCharm->Write();
  hFracGlue->Write();
  hFracGlueUp->Write();
  hFracGlueDown->Write();
  //hFracPU->Write();

  hFracLTJPInv->Write();
  hFracWorkingPointInv->Write();
  if(ppPbPb) hFracDataDrivenInv->Write();
  if(ppPbPb) hFracCharmInv->Write();
  //hFracPUInv->Write();

  fout->Close();
  
}
