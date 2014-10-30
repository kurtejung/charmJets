
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TChain.h"

using namespace std;

void drawJetComb(){
  
  double weight;
  int HLT_PAZeroBiasPixel_SingleTrack_v1;
  int HLT_Jet20_NoJetID_v1, HLT_Jet40_NoJetID_v1, HLT_Jet60_NoJetID_v1, HLT_Jet80_NoJetID_v1, HLT_Jet100_NoJetID_v1;
  int HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl, L1_ZeroBias_Prescl;
  int HLT_Jet20_NoJetID_v1_Prescl, HLT_Jet40_NoJetID_v1_Prescl, HLT_Jet60_NoJetID_v1_Prescl, HLT_Jet80_NoJetID_v1_Prescl, HLT_Jet100_NoJetID_v1_Prescl;
  vector<double> *jtpt=0;
  vector<double> *rawpt=0;

  TH1D *combination = new TH1D("combination","",50,0,400);
  combination->Sumw2();
  TH1D *component[6];
  for(int i=0; i<6; i++){
    component[i] = new TH1D(Form("component_%d",i),"",50,0,400); component[i]->Sumw2();
  }

  TFile *fin1 = new TFile("input/DMesonCJet_NoJetTrgCut_pPbdata_MinBias_ppReco_akPu3PF.root");
  TTree *mb = (TTree*)fin1->Get("ct");
  
  mb->SetBranchAddress("jtpt",&jtpt);
  mb->SetBranchAddress("rawpt",&rawpt);
  mb->SetBranchAddress("weight",&weight);
  mb->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&HLT_PAZeroBiasPixel_SingleTrack_v1);
  mb->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl",&HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl);
  mb->SetBranchAddress("L1_ZeroBias_Prescl",&L1_ZeroBias_Prescl);
  for(int ientry=0; ientry<mb->GetEntries(); ientry++){
    mb->GetEntry(ientry);
    if(HLT_PAZeroBiasPixel_SingleTrack_v1){
      for(unsigned int ijet=0; ijet<jtpt->size(); ijet++){
	if(rawpt->at(ijet)>23){
	  combination->Fill(jtpt->at(ijet),weight);
	  component[0]->Fill(jtpt->at(ijet),HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl*L1_ZeroBias_Prescl);
	}
      }
    }
  }
  cout << "MinBias complete..." << endl;

  TChain *ct = new TChain("ct","");
  ct->Add("input/DMesonCJet_NoJetTrg_AllData_RevCuts_FakeJet2Cut_JetAssn-9-17_part1.root");
  ct->Add("input/DMesonCJet_NoJetTrg_AllData_RevCuts_FakeJet2Cut_JetAssn-9-17_part2.root");
  //TTree *ct = (TTree*)fin2->GetTree("ct");

  ct->SetBranchAddress("jtpt",&jtpt);
  ct->SetBranchAddress("rawpt",&rawpt);
  ct->SetBranchAddress("weight",&weight);
  ct->SetBranchAddress("HLT_Jet20_NoJetID_v1",&HLT_Jet20_NoJetID_v1);
  ct->SetBranchAddress("HLT_Jet40_NoJetID_v1",&HLT_Jet40_NoJetID_v1);
  ct->SetBranchAddress("HLT_Jet60_NoJetID_v1",&HLT_Jet60_NoJetID_v1);
  ct->SetBranchAddress("HLT_Jet80_NoJetID_v1",&HLT_Jet80_NoJetID_v1);
  ct->SetBranchAddress("HLT_Jet100_NoJetID_v1",&HLT_Jet100_NoJetID_v1);
  ct->SetBranchAddress("HLT_Jet20_NoJetID_v1_Prescl",&HLT_Jet20_NoJetID_v1_Prescl);
  ct->SetBranchAddress("HLT_Jet40_NoJetID_v1_Prescl",&HLT_Jet40_NoJetID_v1_Prescl);
  ct->SetBranchAddress("HLT_Jet60_NoJetID_v1_Prescl",&HLT_Jet60_NoJetID_v1_Prescl);
  ct->SetBranchAddress("HLT_Jet80_NoJetID_v1_Prescl",&HLT_Jet80_NoJetID_v1_Prescl);
  ct->SetBranchAddress("HLT_Jet100_NoJetID_v1_Prescl",&HLT_Jet100_NoJetID_v1_Prescl);
  
  for(int ientry=0; ientry<ct->GetEntries(); ientry++){
    ct->GetEntry(ientry);
    
    for(unsigned int ijet=0; ijet<jtpt->size(); ijet++){
      if(rawpt->at(ijet)>23){
	combination->Fill(jtpt->at(ijet),weight);
	if(HLT_Jet20_NoJetID_v1){
	  component[1]->Fill(jtpt->at(ijet),HLT_Jet20_NoJetID_v1_Prescl);
	}
	if(HLT_Jet40_NoJetID_v1){
	  component[2]->Fill(jtpt->at(ijet),HLT_Jet40_NoJetID_v1_Prescl);
	}
	if(HLT_Jet60_NoJetID_v1){
	  component[3]->Fill(jtpt->at(ijet),HLT_Jet60_NoJetID_v1_Prescl);
	}
	if(HLT_Jet80_NoJetID_v1){
	  component[4]->Fill(jtpt->at(ijet),HLT_Jet80_NoJetID_v1_Prescl);
	}
	if(HLT_Jet100_NoJetID_v1){
	  component[5]->Fill(jtpt->at(ijet),HLT_Jet100_NoJetID_v1_Prescl);
	}
      }
    }
  }
  int colors[6] = {2,4,5,6,8,20};
  combination->Draw();
  for(int i=0; i<6; i++){
    component[i]->SetLineColor(colors[i]);
    component[i]->SetMarkerColor(colors[i]);
    component[i]->Draw("same");
  }

  TFile *fout = new TFile("output_jet.root","recreate");
  fout->cd();
  for(int i=0; i<6; i++){
    component[i]->Write();
    combination->Write();
  }
}
