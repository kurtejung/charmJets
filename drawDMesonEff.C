
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#ifdef __CINT__
#pragma link C++ class vector<vector<int> >+;
#endif

using namespace std;

void drawDMesonEff(){

  TFile *f1 = new TFile("DMesonCJet_CJet_PbPbMC_PbPbReco_akPu3PF_993.root");
  TTree *ct = (TTree*)f1->Get("ct");
  
  vector<double> *dCandPt=0, *dCandMass=0, *dCandEta=0, *dCandPhi=0, *jtpt=0, *refpt=0, *rawpt=0;
  vector<double> *dCandType=0, *dCandChildMass=0, *dCandCharge1=0, *dCandCharge2=0;
  vector<double> *refparton_flavorForB=0;
  vector<int> *dCandGenMatchPdg=0, *dCandGenMatchDoubleCountedPdg=0, *dCandGenMatchIdx=0, *dCandGenMatchDoubleCountedIdx=0;
  vector<vector<int> > *hasGenD=0, *hasGenDwithKpi=0, *hasGenDwithKpiIdx=0, *hasGenDIdx=0;;
  ct->SetBranchAddress("dCandPt",&dCandPt);
  ct->SetBranchAddress("jtpt",&jtpt);
  ct->SetBranchAddress("refpt",&refpt);
  ct->SetBranchAddress("rawpt",&rawpt);
  ct->SetBranchAddress("hasGenDwithKpiIdx",&hasGenDwithKpiIdx);
  ct->SetBranchAddress("hasGenDIdx",&hasGenDIdx);
  ct->SetBranchAddress("dCandEta",&dCandEta);
  ct->SetBranchAddress("dCandPhi",&dCandPhi);
  ct->SetBranchAddress("dCandMass",&dCandMass);
  ct->SetBranchAddress("dCandType",&dCandType);
  ct->SetBranchAddress("dCandChildMass",&dCandChildMass);
  ct->SetBranchAddress("dCandCharge1",&dCandCharge1);
  ct->SetBranchAddress("dCandCharge2",&dCandCharge2);
  ct->SetBranchAddress("dCandGenMatchPdg",&dCandGenMatchPdg);
  ct->SetBranchAddress("dCandGenMatchIdx",&dCandGenMatchIdx);
  ct->SetBranchAddress("dCandGenMatchDoubleCountedPdg",&dCandGenMatchDoubleCountedPdg);
  ct->SetBranchAddress("dCandGenMatchDoubleCountedIdx",&dCandGenMatchDoubleCountedIdx);
  ct->SetBranchAddress("refparton_flavorForB",&refparton_flavorForB);
    
  vector<double> *genPt=0, *genEta=0, *genPhi=0;
  vector<int> *genIdx=0, *genPdg=0, *genPtIdx2=0;
  vector<bool> *hasHadronicDecay=0;
  TTree *gt = (TTree*)f1->Get("gTree");
  gt->SetBranchAddress("genPt",&genPt);
  gt->SetBranchAddress("genEta",&genEta);
  gt->SetBranchAddress("genPhi",&genPhi);
  gt->SetBranchAddress("genIdx",&genIdx); //event-level index of particle
  gt->SetBranchAddress("genPdg",&genPdg);
  gt->SetBranchAddress("genPtIdx2",&genPtIdx2);
  gt->SetBranchAddress("hasHadronicDecay",&hasHadronicDecay);

  TH1D *Deff_fg[5];
  TH1D *Deff_bg[5];
  TH1D *Deff_jet_fg[5];
  TH1D *Deff_jet_bg[5];
  for(int i=0; i<5; i++){
    Deff_bg[i] = new TH1D(Form("Deff_bg_%d",i),"",10,0,100); Deff_bg[i]->Sumw2();
    Deff_fg[i] = new TH1D(Form("Deff_fg_%d",i),"",10,0,100); Deff_fg[i]->Sumw2();
    Deff_jet_bg[i] = new TH1D(Form("Deff_jet_bg_%d",i),"",20,0,200); Deff_jet_bg[i]->Sumw2();
    Deff_jet_fg[i] = new TH1D(Form("Deff_jet_fg_%d",i),"",20,0,200); Deff_jet_fg[i]->Sumw2();
  }

  int pdgs[5] = {413, 421, 431, 431, 411};

  for(int ientry=0; ientry<ct->GetEntries(); ientry++){
    ct->GetEntry(ientry);
    gt->GetEntry(ientry);

    //for all gen-level d-mesons, how many are reconstructed?
    //do reconstruction efficiency starting from the gen-level d-meson jet association
    for(unsigned int ijet=0; ijet<jtpt->size(); ijet++){
      if(refpt->at(ijet)<20) continue;
      if(rawpt->at(ijet)<23) continue;
      if(jtpt->at(ijet)<80 || jtpt->at(ijet)>100) continue;
      //if(TMath::Abs(refparton_flavorForB->at(ijet))==4){
	for(int ij=0; ij<5; ij++){
	  bool FF2=false;
	  bool foundFlag=false;
	  double jetptToFill=0;
	  double mesonptToFill=0;
	  for(unsigned int igen=0; igen<genPt->size(); igen++){
	    //if(genPt->at(igen)<10) continue;
	    if(hasGenDwithKpiIdx->at(ijet)[ij] && abs(genPdg->at(igen))==pdgs[ij] /*&& hasHadronicDecay->at(igen)*/){
	      FF2 = true;
	      jetptToFill = jtpt->at(ijet);
	      mesonptToFill = genPt->at(igen);
	      for(unsigned int icand=0; icand<dCandPt->size(); icand++){
		if((dCandGenMatchIdx->at(icand)==genIdx->at(igen)) && !foundFlag && dCandType->at(icand)-1 == ij ){
		  foundFlag=true;
		}
	      }
	    }
	  }
	  if(FF2){
	    Deff_bg[ij]->Fill(mesonptToFill);
	    Deff_jet_bg[ij]->Fill(jetptToFill);
	  }
	  if(foundFlag){
	    Deff_fg[ij]->Fill(mesonptToFill);
	    Deff_jet_fg[ij]->Fill(jetptToFill);
	  }
	}
	//   }
    }
  }

  for(int i=0; i<5; i++){
    //Deff_fg[i]->Divide(Deff_bg[i]);
    //Deff_jet_fg[i]->Divide(Deff_jet_bg[i]);
  }

  Deff_jet_fg[1]->SetXTitle("charm object p_{T}");
  Deff_jet_fg[1]->SetYTitle("Associated Counts");
  Deff_jet_fg[1]->Draw();
  Deff_jet_fg[1]->SetMarkerColor(2);
  Deff_fg[1]->Draw("same");

}
