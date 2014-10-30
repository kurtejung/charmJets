#include <iostream>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"

#ifdef __CINT__
#pragma link C++ class vector<vector<int> >+;
#endif

using namespace std;

const double lowMassBin[cTypes] = {0.14, 1.775, 1.8, 1.8, 1.75};
const double highMassBin[cTypes] = {0.155, 1.95, 2.1, 2.1, 2.0};
const double daughterLowMassBin[cTypes] = {1.8, -1, 1.00, 1.00, 0.5};
const double daughterHighMassBin[cTypes] = {1.92, 1, 1.04, 1.04, 2.05};

void draw_fit(int ich, TH1D* hfitter, int nrb, float rlow, float rhigh, float& sig, float& err)
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
    
    //.. check if need to fit again ...
    cout<<" good fit? (0: no refit, 1: refit w/o range change, 2: refit w/ range change, 3: change both range and fit pars"<<endl;
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
    }
  }
  
  // correct mass bin width
  sig = fit_fun->GetParameter(0)/hfitter->GetBinWidth(1);
  err = fit_fun->GetParError(0)/hfitter->GetBinWidth(1);
  
  cout<<"total number of D: "<<sig<<"+/-"<<err<<endl;
}

//rejects signal region when fitting a pol2 
Double_t poly2reject(Double_t *x, Double_t *par){

  if(x[0] > tightSignalCutLow[globalCharmType] && x[0] < tightSignalCutHi[globalCharmType]){
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0] + par[2]*pow(x[0],2);
}

bool extraDCuts(double mass, double daughterMass, int charge1, int charge2, int type, bool isFG){

  bool dsel = false;
  int type = dCand.type;

  if(dCand.charge1 != dCand.charge2){
    if(isFG){
      if(dCand.type==0){ //dstar 
	if(dCand.mass-dCand.daughterMass>lowMassBin[type] && dCand.mass-dCand.daughterMass < highMassBin[type]){
	  if(dCand.daughterMass>daughterLowMassBin[type] && dCand.daughterMass<daughterHighMassBin[type]){
	    dsel = true;
	  }
	}
      }
      else if(dCand.type==1 || dCand.type==5){ //d0 has no daughters, no selections on D+/-
	if(dCand.mass>lowMassBin[type] && dCand.mass < highMassBin[type]) 
	  dsel = true;
      }
      else{
	if(dCand.mass>lowMassBin[type] && dCand.mass < highMassBin[type] &&
	   dCand.daughterMass>daughterLowMassBin[type] && dCand.daughterMass<daughterHighMassBin[type]) 
	  dsel = true;
      }
    }
    else{ //if using background cuts...
      if(dCand.type==0){
	if(dCand.mass-dCand.daughterMass>lowMassBin[type] && dCand.mass-dCand.daughterMass < highMassBin[type]){
	  dsel = true;
	}
      }
      else if(dCand.type>0){
	if(dCand.mass>lowMassBin[type] && dCand.mass < highMassBin[type]){
	  dsel = true;
	}
      }
    }
  }

  return dsel;
}


void drawDoubleCounting(bool isMC=1){

  const int jetPtBins = 6;
  double ptBoundaries[jetPtBins+1] = {20,40,60,80,100,200,400};
  
  double weight;
  vector<double> *dCandPt=0, *dCandMass=0, *dCandEta=0, *dCandPhi=0, *jtpt=0, *jteta=0, *jtphi=0, *refpt=0, *rawpt=0;
  vector<double> *dCandType=0, *dCandChildMass=0, *dCandCharge1=0, *dCandCharge2=0;
  vector<double> *dCandMatchGenPdg_dau2=0, *dCandMatchGenPdg_index1=0, *dCandMatchGenPdg_index2=0;
  vector<double> *refparton_flavorForB=0;
  vector<vector<int> > *dCandParentPartonIdx=0;
  vector<int> *genMatch=0, *dCandGenMatchPdg=0, *dCandGenMatchDoubleCountedPdg=0, *dCandGenMatchIdx=0, *dCandGenMatchDoubleCountedIdx=0;
  TBranch *bdcandpt=0, *bdcandmass=0, *bdcandeta=0, *bdcandphi=0, *bdcandtype=0, *bdcandchildmass, *bjtpt=0, *bjteta=0, *bjtphi=0, *bdcandcharge1=0, *bdcandcharge2=0;
  vector<vector<int> > *hasGenD=0, *hasGenDwithKpi=0;

  TH1D *matchedD[5];
  TH1D *unmatchedD[5];
  TH1D *jetMatchedD[5];
  TH1D *jetUnmatchedD[5];
  for(int i=0; i<5; i++){
    matchedD[i] = new TH1D(Form("matchedD_%d",i),"",jetPtBins,ptBoundaries);
    unmatchedD[i] = new TH1D(Form("unmatchedD_%d",i),"",jetPtBins,ptBoundaries);
    jetMatchedD[i] = new TH1D(Form("jetMatchedD_%d",i),"",jetPtBins,ptBoundaries);
    jetUnmatchedD[i] = new TH1D(Form("jetUnmatchedD_%d",i),"",jetPtBins,ptBoundaries);
  }

  //d*, d0, ds->phipi, ds->k*k, dpm
  TFile* fin = new TFile("input/DMesonCJet_CJetOnly_pPbMC_ppReco_akPu3PF_448_take4.root","OLD");
  if(!fin){
    cout << "input file not found! " <<endl;
    exit(0);
  }
  fin->cd();
  TTree *ct = (TTree*)fin->Get("ct");

  ct->SetBranchAddress("dCandPt",&dCandPt, &bdcandpt);
  ct->SetBranchAddress("dCandEta",&dCandEta, &bdcandeta);
  ct->SetBranchAddress("dCandPhi",&dCandPhi, &bdcandphi);
  ct->SetBranchAddress("dCandMass",&dCandMass, &bdcandmass);
  ct->SetBranchAddress("dCandType",&dCandType, &bdcandtype);
  ct->SetBranchAddress("dCandChildMass",&dCandChildMass, &bdcandchildmass);
  ct->SetBranchAddress("dCandCharge1",&dCandCharge1, &bdcandcharge1);
  ct->SetBranchAddress("dCandCharge2",&dCandCharge2, &bdcandcharge2);
  if(isMC){
    ct->SetBranchAddress("dCandMatchGenPdg_dau2",&dCandMatchGenPdg_dau2);
    ct->SetBranchAddress("dCandMatchGenPdg_index1",&dCandMatchGenPdg_index1);
    ct->SetBranchAddress("dCandMatchGenPdg_index2",&dCandMatchGenPdg_index2);
    ct->SetBranchAddress("dCandGenMatchPdg",&dCandGenMatchPdg);
    ct->SetBranchAddress("dCandGenMatchIdx",&dCandGenMatchIdx);
    ct->SetBranchAddress("dCandGenMatchDoubleCountedPdg",&dCandGenMatchDoubleCountedPdg);
    ct->SetBranchAddress("dCandGenMatchDoubleCountedIdx",&dCandGenMatchDoubleCountedIdx);
    // ct->SetBranchAddress("dCandGenMatchPdgNoKpi",&dCandGenMatchPdgNoKpi);
  }

  ct->SetBranchAddress("jtpt",&jtpt, &bjtpt);
  ct->SetBranchAddress("jteta",&jteta, &bjteta);
  ct->SetBranchAddress("jtphi",&jtphi, &bjtphi);
  ct->SetBranchAddress("rawpt",&rawpt);
  if(isMC){
    ct->SetBranchAddress("refpt",&refpt);
    ct->SetBranchAddress("refparton_flavorForB",&refparton_flavorForB);
    ct->SetBranchAddress("dCandParentPartonIdx",&dCandParentPartonIdx);
    ct->SetBranchAddress("genMatch",&genMatch);
    ct->SetBranchAddress("hasGenD",&hasGenD);
    ct->SetBranchAddress("hasGenDwithKpi",&hasGenDwithKpi);
  }
  ct->SetBranchAddress("weight",&weight);

  TTree *gt;
  vector<double> *genPt=0, *genEta=0, *genPhi=0;
  vector<int> *genIdx=0, *genPdg=0;
  if(isMC){
    gt = (TTree*)fin->Get("gTree");
    gt->SetBranchAddress("genPt",&genPt);
    gt->SetBranchAddress("genEta",&genEta);
    gt->SetBranchAddress("genPhi",&genPhi);
    gt->SetBranchAddress("genIdx",&genIdx); //event-level index of particle
    gt->SetBranchAddress("genPdg",&genPdg);
  }

  int nentries = ct->GetEntries();
  for(int ientry=0; ientry<nentries; ientry++){
    ct->GetEntry(ientry);
    if(isMC) gt->GetEntry(ientry);
    if (ientry%1000000==0) cout<<" i = "<<ientry<<" out of "<<nentries<<" ("<<(int)(100*(float)ientry/(float)nentries)<<"%)"<<endl;

    if(isMC) weight = 1.;

    for(unsigned int ijet=0; ijet<jtpt->size(); ijet++){
      for(int ij=0; ij<5; ij++){
	bool foundFlag = false;
	if(hasGenDwithKpi->at(ijet)[ij] && genPt->size()>0){
	  for(unsigned int igen=0; igen<genPt->size(); igen++){
	    for(unsigned int icand=0; icand<dCandPt->size(); icand++){
	      if(extraDCuts(dCandMass->at(icand), dCandChildMass->at(icand), dCandType->at(icand), 1)){
		if(dCandGenMatchIdx->at(icand)==genIdx->at(igen) && genPdg->at(igen)==pdgs[ij]){
		  jetMatchedD[ij]->Fill(refpt->at(ijet));
		  foundFlag = true;
		}
	      }
	    }
	  }
	}
	if(!foundFlag){
	  jetUnmatchedD[ij]->Fill(refpt->at(ijet));
	}
      }
    }

    for(int icand=0; icand<dCandPt->size(); icand++){
      bool foundFlag = false;
      for(int igen=0; igen<genPt->size(); igen++){
	for(int ij=0; ij<5; ij++){
	  if(genPdg->at(igen) == pdgs[ij]){
	    if(genIdx->at(igen) == dCandGenMatchIdx->at(icand)){
	      matchedD[ij]->Fill(dCandMass->at(icand));
	      foundFlag = true;
	    }
	  }
	}
      }
      if(!foundFlag){
	unmatchedD[ij]->Fill(dCandMass->at(icand));
      }
    }
    
    const int cTypes = 5;
    TCanvas *cc[cTypes];
    for(int k=0; k<cTypes; k++){
      cc[k] = new TCanvas(Form("cc_type%d",k),Form("cc_type%d",k),1000,800);
      cc[k]->Divide(2,3);
      for(int i=0; i<jetPtBins; i++){
	draw_fit(k,spec[i][j][k],1,(float)range1,(float)range4,dmesonSig[i][k],dmesonSigErr[i][k]);
    }
  }
}
