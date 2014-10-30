
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <math.h>

#include <vector>

#ifdef __CINT__
#pragma link C++ class vector<vector<int> >+;
#endif

using namespace std;

void getGenCJetSpectra(bool useWA=0){
  
  const int jetBins = 10;
  double xjetBins[jetBins+1] = {0,5,10,15,20,40,60,80,100,200,400};
  double brs[5] = {0.0263,0.0388,0.0288,0.0288,0.0913};
  double ffs[5] = {0.235, 0.549, 0.101, 0.101, 0.232};
  double bgMCscaleFactors[5] = {0.4111,0.3981,0.24969,0.2891,0.17444};
  int pdgs[5] = {413,421,431,431,411};
  string labels[5] = {"D*->C-Jet Spectrum","D_{0}->C-Jet Spectrum","D_{s}->C-Jet Spectrum","D_{s}->C-Jet Spectrum","D^{+/-}->C-Jet Spectrum"};
  TH1D *spectra[5], *bgspectra[5], *mctruth[5], *l4corr[5], *backgroundl4[5];
  TH1D *l3corr[5], *backgroundl3[5], *l2corr[5], *backgroundl2[5], *l1corr[5], *backgroundl1[5];
  TH1D *purity[5], *purityXchk[5], *purityData[5];
  TFile *f1 = new TFile("input/DMesonCJet_CJetOnly_pPbMC_ppReco_akPu3PF_512.root");
  TFile *corrs = new TFile("dCandidateToJetAssns_CJetOnly_RevCuts_JetAssn-10-30.root");
  TFile *purFactors = new TFile("purityHistos_MC.root");
  TFile *purFactors2 = new TFile("purityHistos_data.root");
  TTree *ct = (TTree*)f1->Get("ct");

  TF1 *feffline[5];
  for(int i=0; i<5; i++){
    feffline[i] = (TF1*)corrs->Get(Form("feffline_%d",i));
  }
  
  vector<vector<int> > *hasGenD=0;
  vector<vector<int> > *hasGenDwithKpi=0, *hasGenDwithKpiIdx=0;
  vector<double> *refpt=0, *jteta=0, *jtphi=0;
  vector<double> *dCandPt=0, *dCandEta=0, *dCandPhi=0, *dCandMass=0, *dCandChildMass=0;
  vector<vector<int> > *dCandParentPartonIdx=0;
  vector<int> *genMatch=0, *dCandType=0, *dCandCharge1=0, *dCandCharge2=0, *dCandGenMatchPdg=0, *dCandGenMatchIdx=0, *dCandDoubleCountedGenMatchPdg=0;
  double weight;
  ct->SetBranchAddress("refpt",&refpt);
  ct->SetBranchAddress("jteta",&jteta);
  ct->SetBranchAddress("jtphi",&jtphi);
  ct->SetBranchAddress("hasGenD",&hasGenD);
  ct->SetBranchAddress("hasGenDwithKpi",&hasGenDwithKpi);
  ct->SetBranchAddress("hasGenDwithKpiIdx",&hasGenDwithKpiIdx);
  ct->SetBranchAddress("genMatch",&genMatch);
  ct->SetBranchAddress("weight",&weight);

  ct->SetBranchAddress("dCandParentPartonIdx",&dCandParentPartonIdx);
  ct->SetBranchAddress("dCandPt",&dCandPt);
  ct->SetBranchAddress("dCandEta",&dCandEta);
  ct->SetBranchAddress("dCandPhi",&dCandPhi);
  ct->SetBranchAddress("dCandType",&dCandType);
  ct->SetBranchAddress("dCandMass",&dCandMass);
  ct->SetBranchAddress("dCandChildMass",&dCandChildMass);
  ct->SetBranchAddress("dCandCharge1",&dCandCharge1);
  ct->SetBranchAddress("dCandCharge2",&dCandCharge2);
  ct->SetBranchAddress("dCandGenMatchPdg",&dCandGenMatchPdg);
  ct->SetBranchAddress("dCandGenMatchIdx",&dCandGenMatchIdx);
  ct->SetBranchAddress("dCandGenMatchDoubleCountedPdg",&dCandDoubleCountedGenMatchPdg);

  for(int i=0; i<5; i++){
    spectra[i] = new TH1D(Form("spectra_%d",i),"",jetBins,xjetBins);
    spectra[i]->Sumw2();
    bgspectra[i] = new TH1D(Form("bgspectra_%d",i),"",jetBins,xjetBins);
    bgspectra[i]->Sumw2();
    mctruth[i] = new TH1D(Form("mctruth_%d",i),"",jetBins,xjetBins);
    mctruth[i]->Sumw2();
  }

  const double lowMassBin[5] = {0.14, 1.775, 1.8, 1.8, 1.75};
  const double highMassBin[5] = {0.155, 1.95, 2.1, 2.1, 2.0};
  const double daughterLowMassBin[5] = {1.8, -1, 1.00, 1.00, 0.5};
  const double daughterHighMassBin[5] = {1.92, 1, 1.04, 1.04, 2.05};

  const double tightSignalCutLow[5] = {0.144, 1.84, 1.93, 1.93, 1.85};
  const double tightSignalCutHi[5] = {0.147, 1.885, 1.99, 1.99, 1.885};
 
  for(int ientry=0; ientry<ct->GetEntries(); ientry++){
    ct->GetEntry(ientry);
    
    for(unsigned int ijet=0; ijet<refpt->size(); ijet++){
      for(int ij=0; ij<5; ij++){
	if(hasGenDwithKpi->at(ijet)[ij]){
	  mctruth[ij]->Fill(refpt->at(ijet),weight);
	}
      }
      double assocCandPt[5]={0,0,0,0,0};
      double deta, dphi;
      double fgdjetR = 999;
      double bgdjetR = 999;
      bool dsel[5]={0,0,0,0,0};
      double candDR[5]={999,999,999,999,999};
      for(unsigned int icand=0; icand<dCandPt->size(); icand++){

	int type = dCandType->at(icand) - 1;
	//make extra cuts
	//if(fabs(dCandGenMatchPdg->at(icand))==pdgs[type] && !dsel[type] && dCandGenMatchIdx->at(icand)==hasGenDwithKpiIdx->at(ijet)[type]){ //make spectra 100% pure for a minute...
	if(!dsel[type]){
	  if(dCandCharge1->at(icand) != dCandCharge2->at(icand)){
	    if(dCandType->at(icand)==1){ //dstar 
	      if(dCandMass->at(icand)-dCandChildMass->at(icand)>lowMassBin[type] && dCandMass->at(icand)-dCandChildMass->at(icand) < highMassBin[type]){
		if(dCandChildMass->at(icand)>daughterLowMassBin[type] && dCandChildMass->at(icand)<daughterHighMassBin[type]){
		  if(dCandMass->at(icand)-dCandChildMass->at(icand)>tightSignalCutLow[type] && dCandMass->at(icand)-dCandChildMass->at(icand)<tightSignalCutHi[type]){
		    dsel[type]=true;
		    //spectra[type]->Fill(refpt->at(ijet),weight);
		  }
		}
	      }
	    }
	    else if(dCandType->at(icand)==2 || dCandType->at(icand)==5){ //d0 has no daughters, no selections on D+/-
	      if(dCandMass->at(icand)>lowMassBin[type] && dCandMass->at(icand) < highMassBin[type]) {
		if(dCandMass->at(icand)>tightSignalCutLow[type] && dCandMass->at(icand)<tightSignalCutHi[type]){
		  dsel[type]=true;
		  //spectra[type]->Fill(refpt->at(ijet),weight);
		}
	      }
	    }
	    else{
	      if(dCandMass->at(icand)>lowMassBin[type] && dCandMass->at(icand) < highMassBin[type] &&
		 dCandChildMass->at(icand)>daughterLowMassBin[type] && dCandChildMass->at(icand)<daughterHighMassBin[type]){
		if(dCandMass->at(icand)>tightSignalCutLow[type] && dCandMass->at(icand)<tightSignalCutHi[type]){
		  dsel[type] = true;
		  //spectra[type]->Fill(refpt->at(ijet),weight);
		}
	      }
	    }
	  }
	}
	if(dsel[type]){
	  deta = fabs(dCandEta->at(icand)-jteta->at(ijet));
	  dphi = fabs(dCandPhi->at(icand)-jtphi->at(ijet));
	  assocCandPt[type] = dCandPt->at(icand);
	  if(candDR[type]>sqrt(pow(deta,2)+pow(dphi,2))) candDR[type] = sqrt(pow(deta,2)+pow(dphi,2));
	}
      }
      for(int type=0; type<5; type++){ 
	if(dsel[type]){
	  if(candDR[type]<0.3){
	    spectra[type]->Fill(refpt->at(ijet),weight);//*feffline[type]->Eval(assocCandPt[type]));
	  }
	}
      }
    }
  }
  
  for(int i=0; i<5; i++){
    
    l4corr[i] = (TH1D*)corrs->Get(Form("cJetEffvsPt_step4_type%d",i))->Clone(Form("l4corr_%d",i));
    backgroundl4[i] = (TH1D*)(corrs->Get(Form("cJetEffvsPt_bg_step4_type%d",i)))->Clone(Form("backgroundl4_%d",i));
    l4corr[i]->Divide(l4corr[i],backgroundl4[i],1,1,"B");

    l3corr[i] = (TH1D*)corrs->Get(Form("cJetEffvsPt_step3_type%d",i))->Clone(Form("l3corr_%d",i));
    backgroundl3[i] = (TH1D*)(corrs->Get(Form("cJetEffvsPt_bg_step3_type%d",i)))->Clone(Form("backgroundl3_%d",i));
    l3corr[i]->Divide(l3corr[i],backgroundl3[i],1,1,"B");

    l2corr[i] = (TH1D*)corrs->Get(Form("cJetEffvsPt_step2_type%d",i))->Clone(Form("l2corr_%d",i));
    backgroundl2[i] = (TH1D*)(corrs->Get(Form("cJetEffvsPt_bg_step2_type%d",i)))->Clone(Form("backgroundl2_%d",i));
    l2corr[i]->Divide(l2corr[i],backgroundl2[i],1,1,"B");

    l1corr[i] = (TH1D*)corrs->Get(Form("cJetEffvsPt_step1_type%d",i))->Clone(Form("l1corr_%d",i));
    backgroundl1[i] = (TH1D*)(corrs->Get(Form("cJetEffvsPt_bg_step1_type%d",i)))->Clone(Form("backgroundl1_%d",i));
    l1corr[i]->Divide(l1corr[i],backgroundl1[i],1,1,"B");

    purity[i] = (TH1D*)(purFactors->Get(Form("purity_type%d",i)))->Clone(Form("purity_type%d",i));
    purityData[i] = (TH1D*)(purFactors2->Get(Form("purity_type%d",i)))->Clone(Form("purityd_type%d",i));

    bgspectra[i]->Scale(bgMCscaleFactors[i]);
    //spectra[i]->Add(bgspectra[i],-1);

    if(useWA){
      spectra[i]->Scale(1./ffs[i]);
      spectra[i]->Scale(1./brs[i]);
      spectra[i]->Multiply(purity[i]);
    }
    else{
      spectra[i]->Divide(l4corr[i]);
      spectra[i]->Divide(l3corr[i]);
      spectra[i]->Divide(l2corr[i]);
      spectra[i]->Divide(l1corr[i]);

      mctruth[i]->Divide(l4corr[i]);
      mctruth[i]->Divide(l3corr[i]);
      
      purityXchk[i] = (TH1D*)mctruth[i]->Clone(Form("purityXchk_type%d",i));
      purityXchk[i]->Divide(spectra[i]);
      
      spectra[i]->Multiply(purity[i]);
      //spectra[i]->Multiply(purityXchk[i]);
    }

    
  }

  TH1D *tmp = (TH1D*)spectra[1]->Clone("tmp");
  TCanvas *c2 = new TCanvas("c2","",1000,700);
  c2->Divide(2,2);
  for(int i=0; i<5; i++){
    if(i==4) c2->cd(4);
    else c2->cd(i+1);
    spectra[i]->SetXTitle("Ref Jet p_{T}");
    spectra[i]->SetYTitle("Ratio to C(D0) spectrum");
    spectra[i]->SetTitle(labels[i].c_str());
    spectra[i]->GetXaxis()->SetRangeUser(40,400);
    //spectra[i]->SetMaximum(2);
    //spectra[i]->SetMinimum(0);
    //spectra[i]->Divide(tmp);
    //spectra[i]->Divide(mctruth[i]);
    spectra[i]->Draw();
    //purity[i]->Draw();
    purityXchk[i]->SetMarkerColor(kGreen+2);
    //purityXchk[i]->Draw("same");
    purityData[i]->SetMarkerColor(2);
    //purityData[i]->Draw("same");
    mctruth[i]->SetMarkerColor(kGreen+2);
    //mctruth[i]->Divide(tmp);
    mctruth[i]->Draw("same");
    bgspectra[i]->SetMarkerColor(2);
    //bgspectra[i]->Draw("same");
  }

  TFile *fmc = new TFile("purityXcheck.root","recreate");
  fmc->cd();
  for(int i=0; i<5; i++){
    purityXchk[i]->Write();
  }

}
