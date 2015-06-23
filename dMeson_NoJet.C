
#include <iostream>
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"

using std::endl;
using std::cout; 
using std::vector;

#ifdef __CINT__
#pragma link C++ class vector<vector<int> >+;
#endif

const int cTypes = 5;
const double tightSignalCutLow[cTypes] = {0.144, 1.84, 1.93, 1.93, 1.85};
const double tightSignalCutHi[cTypes] = {0.147, 1.885, 1.99, 1.99, 1.885};
const double MesonMass[cTypes] = {0.1455,1.8648,1.9685,1.9685,1.8696};

double draw_fit(int ich, TH1D* hfitter, int nrb, float rlow, float rhigh, float& sig, float& err, float &massSigma, float &purity, float &purityErr)
{
  hfitter->Rebin(nrb);
  // h->SetMarkerSize(0.8);
  //h->SetLineColor(2);
  //h->SetMarkerColor(2);
  //h->SetMarkerStyle(20);
  //h->GetXaxis()->SetNdivisions(505);
  
  hfitter->GetXaxis()->SetRangeUser(rlow, rhigh);
  //.. fit with a Gaussian and pol or expo
  TF1* fit_fun = new TF1("fit_fun", "[0]*(1/[2]/sqrt(6.28)*exp(-0.5*pow((x-[1])/[2], 2))) + pol2(3) + [6]*(1/[8]/sqrt(6.28)*exp(-0.5*pow((x-[7])/[8], 2)))", rlow, rhigh);
  //TF1* fit_fun = new TF1("fit_fun", "[0]*(1/[2]/sqrt(6.28)*exp(-0.5*pow((x-[1])/[2], 2))) + expo(3) + [6]*(1/[8]/sqrt(6.28)*exp(-0.5*pow((x-[7])/[8], 2)))", rlow, rhigh);
  float tot = hfitter->Integral();

  TF1* bg_fit = new TF1("bg_fit","pol2(0)");
  
  float var_mean = 0.01, var_width = 0.01;
  float p0 = tot, p1 = MesonMass[ich], p2 = 0.015;
  if(ich==0) {
    var_mean = 0.005;
    var_width = 5e-4;
    p2 = 1e-3;
  }
  int pass = 1;
  float p6 = tot, p7 = MesonMass[4], p8 = 0.015; //..2nd peak is D+/-
  float zero = 0;
  while(pass) {
    //.. initial value 
    fit_fun->SetParameter(0, p0);
    fit_fun->SetParameter(1, p1);
    fit_fun->SetParameter(2, p2);
    if(ich==2) {//.. 2nd gaussian for Ds->phi+pi
      fit_fun->SetParameter(6, p6);
      fit_fun->SetParameter(7, p7);
      fit_fun->SetParameter(8, p8);
    }
    
    //.. fit constraint ..
    //fit_fun->SetParLimits(0, 0, 1e20);
    fit_fun->SetParLimits(1, TMath::Max(zero, (float)(MesonMass[ich]-var_mean)), MesonMass[ich]+var_mean);
    fit_fun->SetParLimits(2, TMath::Max(zero, (float)(p2-var_width)),  p2+var_width);
    
    if(ich==2) {
      //fit_fun->SetParLimits(6, 0, 1e20);
      fit_fun->SetParLimits(7, TMath::Max(zero, (float)(MesonMass[4]-var_mean)), MesonMass[4]+var_mean);
      fit_fun->SetParLimits(8, TMath::Max(zero, (float)(p8-var_width)),  p8+var_width);
    } else {//.. all D except Ds->phi+pi have only one peak. remove the 2nd peak
      fit_fun->FixParameter(6, 0);
      fit_fun->FixParameter(7, 0);
      fit_fun->FixParameter(8, 0);
    }
    
    //.. using likelyhood. Chi2 bias total counts
    for(int ii = 0; ii<10; ii++)
      hfitter->Fit(fit_fun,  "q", "", rlow, rhigh);
    //!!!-- I can't get the fit error from likelihood fit right 
    //!!!-- the mean is usable for systematic check
    //h->Fit(fit_fun,  "WL", "", rlow, rhigh);
    bg_fit->SetParameter(0,fit_fun->GetParameter(3));
    bg_fit->SetParameter(1,fit_fun->GetParameter(4));
    bg_fit->SetParameter(2,fit_fun->GetParameter(5));
    bg_fit->SetParError(0,fit_fun->GetParError(3));
    bg_fit->SetParError(1,fit_fun->GetParError(4));
    bg_fit->SetParError(2,fit_fun->GetParError(5));
    
    cout << "fit_fun par0: " << fit_fun->GetParameter(0) << " par1: "<< fit_fun->GetParameter(1) << " par2: "<< fit_fun->GetParameter(2) << endl;
    
    //.. draw foreground and background ..
    cout<<" hist: "<<hfitter->GetName()<<endl;
    hfitter->Draw();
    
    TF1* fit_fun_bg = (TF1*)fit_fun->Clone("fit_fun_bg");
    fit_fun_bg->SetParameter(0, 0);
    fit_fun_bg->SetParameter(1, 0);
    fit_fun_bg->SetParameter(2, 0);
    fit_fun_bg->SetParameter(6, 0);
    fit_fun_bg->SetParameter(7, 0);
    fit_fun_bg->SetParameter(8, 0);
    
    fit_fun_bg->SetLineColor(8);
    fit_fun_bg->Draw("same");
    gPad->Update();
    
    pass = 0;
    //.. check if need to fit again ...
    /* cout<<" good fit? (0: no refit, 1: refit w/o range change, 2: refit w/ range change, 3: change both range and fit pars"<<endl;
    scanf("%d", &pass);
    
    if(pass==1) {//.. change the par range only
    cout<<" var(mean), var(width) ?"<<endl;
    scanf("%f%f", &var_mean, &var_width);
    } else if(pass==2) {//.. change the fitting range only
    cout<<" rlow, rhigh? "<<endl;
    scanf("%f%f", &rlow, &rhigh);
    } else if(pass==3) {// change both par and fitting range
    cout<<" var(mean), var(width) ?"<<endl;
    scanf("%f%f", &var_mean, &var_width);
    cout<<" rlow, rhigh? "<<endl;
    scanf("%f%f", &rlow, &rhigh);
    }*/
  }
  
  massSigma = fit_fun->GetParameter(2);

  double sigInt = fit_fun->Integral(tightSignalCutLow[ich], tightSignalCutHi[ich]);
  double sigIntErr = fit_fun->IntegralError(tightSignalCutLow[ich], tightSignalCutHi[ich]);
  double bgInt = bg_fit->Integral(tightSignalCutLow[ich], tightSignalCutHi[ich]);
  double bgIntErr = bg_fit->IntegralError(tightSignalCutLow[ich], tightSignalCutHi[ich]);
  purity = sigInt/(sigInt+bgInt);
  purityErr = sqrt(pow(sigIntErr/sigInt,2)+pow(bgIntErr/bgInt,2));
  purityErr = sqrt(pow(purityErr/(sigInt+bgInt),2)+pow(sigIntErr/sigInt,2));
  cout << "sigInt: "<< sigInt << " sigIntErr: " << sigIntErr << endl;
  cout << "bgInt: "<< bgInt << " bgIntErr: " << bgIntErr << endl;
  cout << "purity: "<< purity << " purErr: "<< purityErr << endl;
  
  // correct mass bin width
  sig = fit_fun->GetParameter(0)/hfitter->GetBinWidth(1);
  err = fit_fun->GetParError(0)/hfitter->GetBinWidth(1);
  
  cout<<"total number of D: "<<sig<<"+/-"<<err<<endl;
  return sig;
}

void dMeson_NoJet(){

  TFile *fin = new TFile("input/DMesonCJet_MergedHalf_addLHCbVars_NoJetTrgCut_pPbdata_ppReco_akPu3PF.root");
  TTree *mt = (TTree*)fin->Get("ct");

  TFile *finMC = new TFile("input/DMesonCJet_QCDJetOnly_Merged_addLHCbVars_pPbMC_centReweight_ppReco_akPu3PF.root");
  TTree *mt2 = (TTree*)finMC->Get("ct");

  if(!mt){ cout << "file or tree not found!!" << endl; exit(0); }
  vector<double> *dCandPt=0, *dCandMass=0, *dCandEta=0, *dCandPhi=0, *dCandChildMass=0, *dCandType=0, *dCandCharge1=0, *dCandCharge2=0;
  vector<double> *ffls3d=0, *falpha0=0, *fprob=0;
  double weight;
  mt->SetBranchAddress("dCandPt",&dCandPt);
  mt->SetBranchAddress("dCandEta",&dCandEta);
  mt->SetBranchAddress("dCandPhi",&dCandPhi);
  mt->SetBranchAddress("dCandMass",&dCandMass);
  mt->SetBranchAddress("dCandType",&dCandType);
  mt->SetBranchAddress("dCandChildMass",&dCandChildMass);
  mt->SetBranchAddress("dCandCharge1",&dCandCharge1);
  mt->SetBranchAddress("dCandCharge2",&dCandCharge2);
  mt->SetBranchAddress("ffls3d",&ffls3d);
  mt->SetBranchAddress("falpha0",&falpha0);
  mt->SetBranchAddress("fprob",&fprob);
  mt->SetBranchAddress("weight",&weight);

  mt2->SetBranchAddress("dCandPt",&dCandPt);
  mt2->SetBranchAddress("dCandEta",&dCandEta);
  mt2->SetBranchAddress("dCandPhi",&dCandPhi);
  mt2->SetBranchAddress("dCandMass",&dCandMass);
  mt2->SetBranchAddress("dCandType",&dCandType);
  mt2->SetBranchAddress("dCandChildMass",&dCandChildMass);
  mt2->SetBranchAddress("dCandCharge1",&dCandCharge1);
  mt2->SetBranchAddress("dCandCharge2",&dCandCharge2);
  mt2->SetBranchAddress("ffls3d",&ffls3d);
  mt2->SetBranchAddress("falpha0",&falpha0);
  mt2->SetBranchAddress("fprob",&fprob);
  mt2->SetBranchAddress("weight",&weight);

  const double ptBins[9]={0,8,12,16,20,25,35,60,500};
  TH1D *mass[8], *massMC[8];
  TH1D *purity = new TH1D("purity","",8,ptBins); purity->Sumw2();
  TH1D *purityMC = new TH1D("purityMC","",8,ptBins); purityMC->Sumw2();
  for(int i=0; i<8; i++){
    mass[i] = new TH1D(Form("mass_pt%i",i),"",30,1.7,2.0); mass[i]->Sumw2();
    massMC[i] = new TH1D(Form("massMC_pt%i",i),"",30,1.7,2.0); massMC[i]->Sumw2();
  }

  cout << "nentries: "<< mt->GetEntries() << endl;
  int nentries = mt->GetEntries();
  for(int ientry=0; ientry<nentries; ientry++){
    mt->GetEntry(ientry);
    for(unsigned int icand=0; icand<dCandPt->size(); icand++){
      if(dCandCharge1->at(icand)!=dCandCharge2->at(icand) && dCandType->at(icand)==2){ //D0's only for right now
        if(ffls3d->at(icand)>3.5 && falpha0->at(icand)<0.06 && fprob->at(icand)>0.25){

         double pt = dCandPt->at(icand);
         int jj=0;
         while(pt>ptBins[jj] && pt<150) jj++;
         if(pt>150 || jj>8) jj = 8;
	   //cout << "pt " << pt << endl;
	   //cout << "filling ptbin: "<< jj << " with mass "<< dCandMass->at(icand) << " weight " << weight << endl;
         mass[jj-1]->Fill(dCandMass->at(icand),weight);
       }
     }
   }
 }

 cout << "nentries mc: "<< mt2->GetEntries() << endl;
 nentries = mt2->GetEntries();
 for(int ientry=0; ientry<nentries; ientry++){
  mt2->GetEntry(ientry);
  for(unsigned int icand=0; icand<dCandPt->size(); icand++){
      if(dCandCharge1->at(icand)!=dCandCharge2->at(icand) && dCandType->at(icand)==2){ //D0's only for right now
        if(ffls3d->at(icand)>3.5 && falpha0->at(icand)<0.06 && fprob->at(icand)>0.25){

         double pt = dCandPt->at(icand);
         int jj=0;
         while(pt>ptBins[jj] && pt<150) jj++;
         if(pt>150 || jj>8) jj = 8;
     //cout << "pt " << pt << endl;
     //cout << "filling ptbin: "<< jj << " with mass "<< dCandMass->at(icand) << " weight " << weight << endl;
         massMC[jj-1]->Fill(dCandMass->at(icand),weight);
       }
     }
   }
 }

  TCanvas *call = new TCanvas("call","",800,600);
  call->Divide(2,2);
  for(int i=0; i<4; i++){
    call->cd(i+1);
    mass[i]->Draw();
    
    float signal=0, err=0, massSigma=0, purityVal=0, purityErr=0;
    draw_fit(1, mass[i], 1, 1.7, 2.0, signal, err, massSigma, purityVal, purityErr);
    purity->SetBinContent(i+1,purityVal);
    purity->SetBinError(i+1,purityErr);
  }

  TCanvas *cmc = new TCanvas("cmc","",800,600);
  cmc->Divide(2,2);
  for(int i=0; i<4; i++){
    cmc->cd(i+1);
    massMC[i]->SetMarkerColor(2);
    massMC[i]->Draw();
    float signal=0, err=0, massSigma=0, purityVal=0, purityErr=0;
    draw_fit(1, massMC[i], 1, 1.7, 2.0, signal, err, massSigma, purityVal, purityErr);
    purityMC->SetBinContent(i+1,purityVal);
    purityMC->SetBinError(i+1,purityErr);
  }
  
  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->cd();
  purity->Draw();
  purityMC->SetMarkerColor(2);
  purityMC->Draw("same");

  TFile *fout = new TFile("bgScaleFactors.root","recreate");
  fout->cd();
  purity->Write();
  purityMC->Write();

}
