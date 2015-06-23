#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLatex.h"

using namespace std;

void drawDCandStoB(){

  TFile *f = new TFile("input/DMesonCJet_MergedHalf_addLHCbVars_NoJetTrgCut_pPbdata_ppReco_akPu3PF.root");
  TTree *ct = (TTree*)f->Get("ct");

  double ptCuts[11] = {3.5, 4.5, 5.5, 7.0, 9.0, 11.0, 13.0, 16.0, 20., 28., 40.};
  double d0cuts[10] = {3.81, 3.81, 3.91, 3.92, 3.77, 3.57, 3.37, 2.88, 2.71, 2.63};
  double alphaCuts[10] = {0.056, 0.056, 0.048, 0.060, 0.063, 0.057, 0.053, 0.059, 0.059, 0.042};
  double vtxProbCuts[10] = {0.267, 0.267, 0.209, 0.113, 0.093, 0.057, 0.75, 0.025, 0.68, 0.018};

  vector<double> *dCandPt=0, *dCandEta=0, *dCandMass=0, *dCandType=0, *falpha0=0, *fprob0=0, *ffls3d=0;
  double weight;
  ct->SetBranchAddress("dCandPt",&dCandPt);
  ct->SetBranchAddress("dCandEta",&dCandEta);
  ct->SetBranchAddress("dCandMass",&dCandMass);
  ct->SetBranchAddress("dCandType",&dCandType);
  ct->SetBranchAddress("falpha0",&falpha0);
  ct->SetBranchAddress("fprob",&fprob0);
  ct->SetBranchAddress("ffls3d",&ffls3d);
  ct->SetBranchAddress("weight",&weight);

  TH1D *candMass[10];
  for(int i=0; i<10; i++){
    candMass[i] = new TH1D(Form("candMass_%d",i),"",25,1.75,2.05);
    candMass[i]->Sumw2();
    candMass[i]->SetXTitle("D-Candidate Mass");
    candMass[i]->SetYTitle("Counts");
  }
  
  for(int ientries=0; ientries<ct->GetEntries(); ientries++){
    ct->GetEntry(ientries);

    for(unsigned int icand=0; icand<dCandPt->size(); icand++){
      int j=0;
      if(dCandPt->at(icand)<3.5) continue;
      while((dCandPt->at(icand)>ptCuts[j]) && (j<10)) j++;
      j--;
      //cout << "dCandPt: "<< dCandPt->at(icand) << " bin: "<< ptCuts[j] << " to " << ptCuts[j+1] << endl;
      
      if(fabs(dCandEta->at(icand))<1.0 && dCandType->at(icand)==2){
	if(falpha0->at(icand)<alphaCuts[j] && fprob0->at(icand)>vtxProbCuts[j] && ffls3d->at(icand)>d0cuts[j]){
	  candMass[j]->Fill(dCandMass->at(icand));
	}
      }						       						       
    }
  }

  TLatex *l1[10];
  for(int i=0; i<9; i++){
    l1[i] = new TLatex(1.9,500,Form("%.1f < D pT < %.1f",ptCuts[i],ptCuts[i+1]));
  }
  
  TCanvas *cc = new TCanvas("cc","",1000,1000);
  cc->Divide(3,3);
  for(int i=0; i<9; i++){
    cc->cd(i+1);
    candMass[i]->Draw();
    l1[i]->SetY(candMass[i]->GetMaximum());
    l1[i]->Draw("same");
  }
  
}
