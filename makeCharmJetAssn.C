#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLatex.h"

#ifdef __CINT__
#pragma link C++ class vector<vector<int> >+;
#endif

class candidate{
public:
  double mass;
  double pt;
  double daughterMass;
  double eta;
  double phi;
  int type;
  int charge1;
  int charge2;
  int genPdg_dau2;
  int genPdg_idx1;
  int genPdg_idx2;

  candidate(){
  }
  ~candidate(){
  }
  candidate(double m, double p, double dm, double e, double ph, int t, int c1, int c2, int gpdg_d2, int gpdg_idx1, int gpdg_idx2){
    mass = m;
    pt = p;
    daughterMass = dm;
    eta = e;
    phi = ph;
    type = t;
    charge1 = c1;
    charge2 = c2;
    genPdg_dau2 = gpdg_d2;
    genPdg_idx1 = gpdg_idx1;
    genPdg_idx2 = gpdg_idx2;
  }
  void setStuff(double m, double p, double dm, double e, double ph, int t, int c1, int c2, int gpdg_d2, int gpdg_idx1, int gpdg_idx2){
    mass = m;
    pt = p;
    daughterMass = dm;
    eta = e;
    phi = ph;
    type = t;
    charge1 = c1;
    charge2 = c2;
    genPdg_dau2 = gpdg_d2;
    genPdg_idx1 = gpdg_idx1;
    genPdg_idx2 = gpdg_idx2;
  }
};

//Mass selections to filter on resonances: D*, D0, Ds->PhiPi, Ds->KK*, D+-
const int cTypes = 5;
const double lowMassBin[cTypes] = {0.14, 1.775, 1.8, 1.8, 1.75};
const double highMassBin[cTypes] = {0.155, 1.95, 2.1, 2.1, 2.0};
const double daughterLowMassBin[cTypes] = {1.8, -1, 1.00, 1.00, 0.5};
const double daughterHighMassBin[cTypes] = {1.92, 1, 1.04, 1.04, 2.05}; 

double findClosestJet(vector<double> *dCandEta, vector<double> *dCandPhi, double jetEta, double jetPhi){

  //check to ensure the vector sizes are the same
  //Can be commented out without issue if assert isn't playing nice with CMSSW
  assert(dCandEta->size() == dCandPhi->size());
  //
  unsigned int vecSize = dCandEta->size();
  int candno = -1;
  double closestApproach=9999;
  for(unsigned int icand=0; icand<vecSize; icand++){
    double deta = fabs(dCandEta->at(icand) - jetEta);
    double dphi = fabs(dCandPhi->at(icand) - jetPhi);
    if(sqrt(pow(deta,2)+pow(dphi,2)) < closestApproach){
      closestApproach = sqrt(pow(deta,2)+pow(dphi,2));
      candno = icand;
    }
  }

  return candno;
}

double snapToBin(TH1D *input, double cut, bool ishi){
  int bin = input->FindBin(cut);
  double snapped=0.;

  if(ishi) snapped = input->GetBinLowEdge(bin)+input->GetBinWidth(bin);
  else snapped = input->GetBinLowEdge(bin);

  return snapped;  
}

double findClosestJet(double dCandEta, double dCandPhi, vector<double> *jetEta, vector<double> *jetPhi){

  //check to ensure the vector sizes are the same
  //Can be commented out without issue if assert isn't playing nice with CMSSW
  assert(jetEta->size() == jetPhi->size());
  //
  unsigned int vecSize = jetEta->size();
  int candno = -1;
  double closestApproach=9999;
  for(unsigned int icand=0; icand<vecSize; icand++){
    double deta = fabs(jetEta->at(icand) - dCandEta);
    double dphi = fabs(jetPhi->at(icand) - dCandPhi);
    if(sqrt(pow(deta,2)+pow(dphi,2)) < closestApproach){
      closestApproach = sqrt(pow(deta,2)+pow(dphi,2));
      candno = icand;
    }
  }

  return candno;
}

bool extraDCuts(candidate dCand, bool isFG){

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
      else if(dCand.type==1 || dCand.type==4){ //d0 has no daughters, no selections on D+/-
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

void makeCharmJetAssn(int isMC=0, bool doDraw = 0){

  const int ptBins = 6;
  const int jetPtBins = 6;
  const int djetBins = 1;

  const double djetR_testCut = 0.3;
  const bool doJet80Only = false;

  bool doMassDR = false; //turn on if you wanna apply dR selections to the meson mass plots

  int ptBoundaries[jetPtBins+1] = {20,40,60,80,100,200,400};
  int dCandptBoundaries[ptBins+1] = {0,5,10,20,40,100,400};
  double djetBoundaries[djetBins+1] = {0,9999};
  string Dlabels[6] = {"D^{+/-} ","D_{s}->K*K ","D_{s}->#pi#phi ","D* ","D_{0} ","C-Jet Only "};

  //Ds is not implemented yet - should be pdg=431
  int pdgs[6] = {411, -999, -999, 413, 421, -999}; 

  TFile *fout;
  if(!isMC) fout = new TFile("dCandidateToJetAssns_halfData_RevCuts_JetAssn-10-3_fixChgSign.root","recreate");
  else if(isMC==1) fout = new TFile("dCandidateToJetAssns_Dpm_RevCuts_JetAssn-8-11.root","recreate");
  else if(isMC==2) fout = new TFile("dCandidateToJetAssns_Ds2KstarK_RevCuts_JetAssn-8-11.root","recreate");
  else if(isMC==3) fout = new TFile("dCandidateToJetAssns_Ds2phipi_RevCuts_JetAssn-8-11.root","recreate");
  else if(isMC==4) fout = new TFile("dCandidateToJetAssns_Dstar_RevCuts_JetAssn-8-11.root","recreate");
  else if(isMC==5) fout = new TFile("dCandidateToJetAssns_D02Kpi_RevCuts_JetAssn-8-11.root","recreate");
  else if(isMC==6) fout = new TFile("dCandidateToJetAssns_CJetOnly_RevCuts_JetAssn-8-11.root","recreate");
  else{
    cout << "isMC cannot be more than 5! Is this a new D-decay type??" << endl;
    exit(0);
  }
 
  //d*, d0, ds->phipi, ds->k*k, dpm

  const double tightSignalCutLow[cTypes] = {0.144, 1.84, 1.93, 1.93, 1.85};
  const double tightSignalCutHi[cTypes] = {0.147, 1.885, 1.99, 1.99, 1.885};

  const int jetBins = 10;
  double xjetBins[jetBins+1] = {0,5,10,15,20,40,60,80,100,200,400};
  //double xjetBins[jetBins+1] = {20,40,60,80,100,200,400};

  const int ndrbins = 20;
  double xdrbins[ndrbins+1] = {-10,0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.3,0.4,0.6,1.,1.5,2,2.5,3,3.5,4,20};

  TH1D *allJetSpec = new TH1D("allJetSpec","",jetBins,xjetBins); allJetSpec->Sumw2();
  TH1D *CjetSpec[cTypes];
  TH1D *bgCjetSpec[cTypes];

  TH1D *jetFF[ptBins];
  for(int i=0; i<jetPtBins; i++){
    jetFF[i] = new TH1D(Form("jetFF_%d",i),"",30,0,1.5); jetFF[i]->Sumw2();
  }

  TH1D *proximity[jetPtBins][cTypes];
  TH1D *proximityEta[jetPtBins][cTypes];
  TH1D *proximityPhi[jetPtBins][cTypes];
  TH1D *bgProximity[jetPtBins][cTypes];
  TH1D *bgProximityEta[jetPtBins][cTypes];
  TH1D *bgProximityPhi[jetPtBins][cTypes];
  TH1D *assns[jetPtBins][djetBins][cTypes];
  TH1D *assnsUnweight[jetPtBins][djetBins][cTypes];
  TH1D *bgAssns[jetPtBins][djetBins][cTypes];
  TH1D *bgAssnsUnweight[jetPtBins][djetBins][cTypes];

  TH1D *cJetEffvsPt = new TH1D("cJetEffvsPt","",jetBins,xjetBins); cJetEffvsPt->Sumw2();
  TH1D *cJetEffvsdR = new TH1D("cJetEffvsdR","",ndrbins,xdrbins); cJetEffvsdR->Sumw2();
  TH1D *cJetEffvsPt_bg = new TH1D("cJetEffvsPt_bg","",jetBins,xjetBins); cJetEffvsPt_bg->Sumw2();
  TH1D *cJetEffvsdR_bg = new TH1D("cJetEffvsdR_bg","",ndrbins,xdrbins); cJetEffvsdR_bg->Sumw2();

  TH1D *cJetPurvsPt = new TH1D("cJetPurvsPt","",jetBins,xjetBins); cJetPurvsPt->Sumw2();
  TH1D *cJetEffvsPt_rej = new TH1D("cJetEffvsPt_rej","",jetBins,xjetBins); cJetEffvsPt_rej->Sumw2();
  TH1D *cJetPurvsdR = new TH1D("cJetPurvsdR","",ndrbins,xdrbins); cJetPurvsdR->Sumw2();
  TH1D *cJetPurvsPt_bg = new TH1D("cJetPurvsP_bg","",jetBins,xjetBins); cJetPurvsPt_bg->Sumw2();
  TH1D *cJetPurvsdR_bg = new TH1D("cJetPurvsdR_bg","",ndrbins,xdrbins); cJetPurvsdR_bg->Sumw2();

  TH1D *jetAssnBot = new TH1D("jetAssnBot","",jetBins,xjetBins); jetAssnBot->Sumw2();
  TH1D *jetAssnTop = new TH1D("jetAssnTop","",jetBins,xjetBins); jetAssnTop->Sumw2();

  TH1D *dRTestHist = new TH1D("dRTestHist","",30,0,3.14); dRTestHist->Sumw2();
  TH1D *unmatchedCJets = new TH1D("unmatchedCJets","",jetBins,xjetBins); unmatchedCJets->Sumw2();
  TH1D *cJetEffvsPt_bg_step1 = new TH1D("cJetEffvsPt_bg_step1","",jetBins,xjetBins); cJetEffvsPt_bg_step1->Sumw2();
  TH1D *cJetEffvsPt_step1 = new TH1D("cJetEffvsPt_step1","",jetBins,xjetBins); cJetEffvsPt_step1->Sumw2();
  TH1D *cJetEffvsPt_rej_step1 = new TH1D("cJetEffvsPt_rej_step1","",jetBins,xjetBins); cJetEffvsPt_rej_step1->Sumw2();

  TH1D *cJetEffvsPt_bg_step2 = new TH1D("cJetEffvsPt_bg_step2","",jetBins,xjetBins); cJetEffvsPt_bg_step2->Sumw2();
  TH1D *cJetEffvsPt_step2 = new TH1D("cJetEffvsPt_step2","",jetBins,xjetBins); cJetEffvsPt_step2->Sumw2();
  TH1D *cJetEffvsPt_rej_step2 = new TH1D("cJetEffvsPt_rej_step2","",jetBins,xjetBins); cJetEffvsPt_rej_step2->Sumw2();

  TH1D *cJetEffvsPt_step3[5]; //one for each d-meson -> hadron branching ratio
  for(int i=0; i<5; i++){
    cJetEffvsPt_step3[i] = new TH1D(Form("cJetEffvsPt_step3_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_step3[i]->Sumw2();
  }
  TH1D *cJetEffvsPt_bg_step3 = new TH1D("cJetEffvsPt_bg_step3","",jetBins,xjetBins); cJetEffvsPt_bg_step3->Sumw2();
  TH1D *cJetEffvsPt_rej_step3 = new TH1D("cJetEffvsPt_rej_step3","",jetBins,xjetBins); cJetEffvsPt_rej_step3->Sumw2();

  for(int i=0; i<jetPtBins; i++){
    for(int j=0; j<djetBins; j++){
      for(int k=0; k<cTypes; k++){
	assns[i][j][k] = new TH1D(Form("assns_pt%d_djet%d_type%d",i,j,k),"",35,lowMassBin[k],highMassBin[k]);
	assns[i][j][k]->Sumw2();
	assns[i][j][k]->SetXTitle("Candidate Mass [GeV/c^{2}]");
	assns[i][j][k]->SetYTitle("Weighted Counts");
	
	assnsUnweight[i][j][k] = new TH1D(Form("assns_unweight_pt%d_djet%d_type%d",i,j,k),"",35,lowMassBin[k],highMassBin[k]);
	assnsUnweight[i][j][k]->Sumw2();
	assnsUnweight[i][j][k]->SetXTitle("Candidate Mass [GeV/c^{2}]");
	assnsUnweight[i][j][k]->SetYTitle("Unweighted Counts");
	
	bgAssns[i][j][k] = new TH1D(Form("bgAssns_pt%d_djet%d_type%d",i,j,k),"",35,lowMassBin[k],highMassBin[k]);
	//bgAssnsUnweight[i][j][k]->Sumw2();
	bgAssns[i][j][k]->SetXTitle("Candidate Mass [GeV/c^{2}]");
	bgAssns[i][j][k]->SetYTitle("Unweighted Counts");

	bgAssnsUnweight[i][j][k] = new TH1D(Form("bgAssns_unweight_pt%d_djet%d_type%d",i,j,k),"",35,lowMassBin[k],highMassBin[k]);
	//bgAssnsUnweight[i][j][k]->Sumw2();
	bgAssnsUnweight[i][j][k]->SetXTitle("Candidate Mass [GeV/c^{2}]");
	bgAssnsUnweight[i][j][k]->SetYTitle("Unweighted Counts");

      }
    }
  }
  
 
  for(int i=0; i<jetPtBins; i++){
    for(int k=0; k<cTypes; k++){
      proximity[i][k] = new TH1D(Form("proximity_pt%d_type%d",i,k),"",ndrbins,xdrbins);
      proximityEta[i][k] = new TH1D(Form("proximityEta_pt%d_type%d",i,k),"",ndrbins,xdrbins);
      proximityPhi[i][k] = new TH1D(Form("proximityPhi_pt%d_type%d",i,k),"",ndrbins,xdrbins);
      bgProximity[i][k] = new TH1D(Form("bgProximity_pt%d_type%d",i,k),"",ndrbins,xdrbins);
      bgProximityEta[i][k] = new TH1D(Form("bgProximityEta_pt%d_type%d",i,k),"",ndrbins,xdrbins);
      bgProximityPhi[i][k] = new TH1D(Form("bgProximityPhi_pt%d_type%d",i,k),"",ndrbins,xdrbins);
    }
  }
  
  for(int k=0; k<cTypes; k++){
    CjetSpec[k] = new TH1D(Form("CjetSpec_type%d",k),"",jetBins,xjetBins); CjetSpec[k]->Sumw2();
    bgCjetSpec[k] = new TH1D(Form("bgCjetSpec_type%d",k),"",jetBins,xjetBins); bgCjetSpec[k]->Sumw2();
  }

  double weight;
  int HLT_Jet20_NoJetID_v1, HLT_Jet40_NoJetID_v1, HLT_Jet60_NoJetID_v1, HLT_Jet80_NoJetID_v1, HLT_Jet100_NoJetID_v1;
  vector<double> *dCandPt=0, *dCandMass=0, *dCandEta=0, *dCandPhi=0, *jtpt=0, *jteta=0, *jtphi=0, *rawpt=0;
  vector<double> *dCandType=0, *dCandChildMass=0, *dCandCharge1=0, *dCandCharge2=0;
  vector<double> *dCandMatchGenPdg_dau2=0, *dCandMatchGenPdg_index1=0, *dCandMatchGenPdg_index2=0;
  vector<double> *refparton_flavorForB=0;
  vector<vector<int> > *dCandParentPartonIdx=0;
  vector<int> *genMatch=0, *hasGenD=0, *dCandGenMatchPdg=0;
  TBranch *bdcandpt=0, *bdcandmass=0, *bdcandeta=0, *bdcandphi=0, *bdcandtype=0, *bdcandchildmass, *bjtpt=0, *bjteta=0, *bjtphi=0, *bdcandcharge1=0, *bdcandcharge2=0;

  TFile *fin;
  if(!isMC) fin = new TFile("input/DMesonCJet_NoJetTrgCut_pPbdata_HighPt_JianSubset_ppReco_akPu3PF.root","OLD");
  if(isMC==1) fin = new TFile("DMesonCJet_DpmEmbed_pPbMC_ppReco_akPu3PF_1.root","OLD");
  if(isMC==2) fin = new TFile("DMesonCJet_Ds2KstarKEmbed_pPbMC_ppReco_akPu3PF_1.root","OLD");
  if(isMC==3) fin = new TFile("DMesonCJet_Ds2PhiPi_pPbMC_ppReco_akPu3PF_1.root","OLD");
  if(isMC==4) fin = new TFile("DMesonCJet_Dstar_pPbMC_ppReco_akPu3PF_1.root","OLD");
  if(isMC==5) fin = new TFile("DMesonCJet_Dzero_pPbMC_ppReco_akPu3PF_1.root","OLD");
  if(isMC==6) fin = new TFile("DMesonCJet_CJetOnly_pPbMC_ppReco_akPu3PF_9.root","OLD");
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
  }

  ct->SetBranchAddress("jtpt",&jtpt, &bjtpt);
  ct->SetBranchAddress("jteta",&jteta, &bjteta);
  ct->SetBranchAddress("jtphi",&jtphi, &bjtphi);
  ct->SetBranchAddress("rawpt",&rawpt);
  if(isMC){
    ct->SetBranchAddress("refparton_flavorForB",&refparton_flavorForB);
    ct->SetBranchAddress("dCandParentPartonIdx",&dCandParentPartonIdx);
    ct->SetBranchAddress("genMatch",&genMatch);
    ct->SetBranchAddress("hasGenD",&hasGenD);
  }
  ct->SetBranchAddress("weight",&weight);
  ct->SetBranchAddress("HLT_Jet20_NoJetID_v1",&HLT_Jet20_NoJetID_v1);
  ct->SetBranchAddress("HLT_Jet40_NoJetID_v1",&HLT_Jet40_NoJetID_v1);
  ct->SetBranchAddress("HLT_Jet60_NoJetID_v1",&HLT_Jet60_NoJetID_v1);
  ct->SetBranchAddress("HLT_Jet80_NoJetID_v1",&HLT_Jet80_NoJetID_v1);
  ct->SetBranchAddress("HLT_Jet100_NoJetID_v1",&HLT_Jet100_NoJetID_v1);
  
  if(doJet80Only) cout << "WARNING! Applying Jet 80 trigger selection!!" << endl;
  int nentries = ct->GetEntries();
  for(int ientry=0; ientry<nentries; ientry++){
    ct->GetEntry(ientry);
    if (ientry%1000000==0) cout<<" i = "<<ientry<<" out of "<<nentries<<" ("<<(int)(100*(float)ientry/(float)nentries)<<"%)"<<endl;

    if(isMC){
      if(dCandPt->size() != dCandMass->size() || dCandPt->size() != dCandEta->size() || dCandPt->size() != dCandPhi->size() || dCandPt->size() != dCandParentPartonIdx->size()){
	cout << "Warning! dCandidate sizes are different in entry: "<< ientry << "!" << endl;
	cout << "dCandPt size: "<< dCandPt->size() << " dcandparentSize: "<< dCandParentPartonIdx->size() << endl;
	cout << "Skipping event..." << endl;
	continue;
      }
    }
    else if(dCandPt->size() != dCandMass->size() || dCandPt->size() != dCandEta->size() || dCandPt->size() != dCandPhi->size()){
      cout << "Warning! dCandidate sizes are different in entry: "<< ientry << "!" << endl;
      cout << "Skipping event..." << endl;
      continue;
    }

    double maxjetpt=-1;
    //for all events with d-mesons, how many have an associated jet?
    if(jtpt->size()>0){
      for(unsigned int icand=0; icand<jtpt->size(); icand++){
	if(jtpt->at(icand) > maxjetpt) maxjetpt = jtpt->at(icand);
      }
      jetAssnBot->Fill(maxjetpt);
      if(dCandPt->size()>0){
	jetAssnTop->Fill(maxjetpt);
      }
    }
    
    if(doJet80Only){ if(!HLT_Jet80_NoJetID_v1){ continue; }}
    else if(!HLT_Jet20_NoJetID_v1 && !HLT_Jet40_NoJetID_v1 && !HLT_Jet60_NoJetID_v1 && !HLT_Jet80_NoJetID_v1 && !HLT_Jet100_NoJetID_v1) continue;

    for(unsigned int ijet=0; ijet<jtpt->size(); ijet++){
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
	  if(candDR[type]>sqrt(pow(deta,2)+pow(dphi,2))) candDR[type] = sqrt(pow(deta,2)+pow(dphi,2));
	}
      }
      for(int type=0; type<5; type++){ 
	if(dsel[type]){
	  if(candDR[type]<0.3){
	    CjetSpec[type]->Fill(rawpt->at(ijet),weight);
	  }
	}
      }
    }


    for(unsigned int ijet=0; ijet<jtpt->size(); ijet++){
      allJetSpec->Fill(jtpt->at(ijet));
      bool dRFlag = false;
      bool dRFlag2 = false;
      bool anotherTestFlag=false;
      double minDR=999;
      bool checkCspec[5] = {0,0,0,0,0};
      bool bgcheckCspec[5] = {0,0,0,0,0};
      if(rawpt->at(ijet)<23) continue;
      int ipt=0;
      while(jtpt->at(ijet) > ptBoundaries[ipt+1] && ipt<(jetPtBins-1)) ipt++;
      for(unsigned int icand=0; icand<dCandPt->size(); icand++){

	//if(dCandPt->at(icand) > maxjetpt) continue;
	int iCandpt=0;
	//int idjet=0;
	
	//double jetprox = findClosestJet(dCandEta, dCandPhi, jteta->at(ijet), jtphi->at(ijet));
	//if(jetprox==-1) continue;
	double deta = fabs(dCandEta->at(icand)-jteta->at(ijet));
	double dphi = fabs(dCandPhi->at(icand)-jtphi->at(ijet));
	double djetR = sqrt(pow(deta,2)+pow(dphi,2));

	//while(djetR>djetBoundaries[idjet+1] && idjet<(djetBins-1)) idjet++;     
	
	while(dCandPt->at(icand) > dCandptBoundaries[iCandpt+1] && iCandpt<(ptBins-1)) iCandpt++;
	int type = dCandType->at(icand)-1;
	
	//ipt = iCandpt;
	candidate tempCand;
        tempCand.setStuff(dCandMass->at(icand), dCandPt->at(icand), dCandChildMass->at(icand), dCandEta->at(icand), dCandPhi->at(icand), type, dCandCharge1->at(icand), dCandCharge2->at(icand),-999,-999,-999);

	if(extraDCuts(tempCand,1)){
	  if(dCandType->at(icand)==1){ //start with D* meson since its different
	    if(djetR<djetR_testCut || !doMassDR){
	      assns[ipt][0][type]->Fill(dCandMass->at(icand) - dCandChildMass->at(icand), weight);
	      assnsUnweight[ipt][0][type]->Fill(dCandMass->at(icand)-dCandChildMass->at(icand));
	    }
	    double lowball = snapToBin(assns[ipt][0][type],tightSignalCutLow[type],0);
	    double hiball = snapToBin(assns[ipt][0][type],tightSignalCutHi[type],1);
	    if((dCandMass->at(icand)-dCandChildMass->at(icand)) > lowball && (dCandMass->at(icand)-dCandChildMass->at(icand)) < hiball){
	      proximity[ipt][type]->Fill(djetR);
	      proximityEta[ipt][type]->Fill(deta);
	      proximityPhi[ipt][type]->Fill(dphi);
	      if(djetR<djetR_testCut){
		checkCspec[type] = 1; //CjetSpec[type]->Fill(jtpt->at(ijet));
		jetFF[ipt]->Fill(dCandPt->at(icand)/jtpt->at(ijet));
	      }
	    }
	    else{
	      if(djetR<djetR_testCut){
		bgcheckCspec[type] = 1;//->Fill(jtpt->at(ijet));
	      }
	      bgProximity[ipt][type]->Fill(djetR);
	      bgProximityEta[ipt][type]->Fill(deta);
	      bgProximityPhi[ipt][type]->Fill(dphi);
	    }
	  }
	  else{ //do all other d-mesons now
	    if(djetR<djetR_testCut || !doMassDR){
	      assns[ipt][0][type]->Fill(dCandMass->at(icand), weight);
	      assnsUnweight[ipt][0][type]->Fill(dCandMass->at(icand));
	    }
	    double lowball = snapToBin(assns[ipt][0][type],tightSignalCutLow[type],0);
	    double hiball = snapToBin(assns[ipt][0][type],tightSignalCutHi[type],1);
	    if(dCandMass->at(icand) > lowball && dCandMass->at(icand) < hiball){
	      proximity[ipt][type]->Fill(djetR);
	      proximityEta[ipt][type]->Fill(deta);
	      proximityPhi[ipt][type]->Fill(dphi);
	      if(djetR<djetR_testCut){
		checkCspec[type] = 1; //CjetSpec[type]->Fill(jtpt->at(ijet));
		jetFF[ipt]->Fill(dCandPt->at(icand)/jtpt->at(ijet));
	      }
	    }
	    else{
	      if(djetR<djetR_testCut){
		bgcheckCspec[type] = 1; //bgCjetSpec[type]->Fill(jtpt->at(ijet));
	      }
	      bgProximity[ipt][type]->Fill(djetR);
	      bgProximityEta[ipt][type]->Fill(deta);
	      bgProximityPhi[ipt][type]->Fill(dphi);
	    }
	  }
	}
	else if(extraDCuts(tempCand,0)){
	  if(djetR<djetR_testCut || !doMassDR){
	    if(dCandType->at(icand)==1){
	      bgAssnsUnweight[ipt][0][type]->Fill(dCandMass->at(icand) - dCandChildMass->at(icand));
	      bgAssns[ipt][0][type]->Fill((dCandMass->at(icand) - dCandChildMass->at(icand)),weight);
	    }
	    else{
	      bgAssnsUnweight[ipt][0][type]->Fill(dCandMass->at(icand));
	      bgAssns[ipt][0][type]->Fill(dCandMass->at(icand),weight);
	    }
	  }
	}
      } //close candidate loop
      for(int ii=0; ii<5; ii++){
	if(checkCspec[ii]){
	  //CjetSpec[ii]->Fill(jtpt->at(ijet),weight);
	}
	else if(bgcheckCspec[ii]){
	  //bgCjetSpec[ii]->Fill(jtpt->at(ijet),weight);
	}
      }
      
      if(isMC && !anotherTestFlag && TMath::Abs(refparton_flavorForB->at(ijet))==4){
	unmatchedCJets->Fill(jtpt->at(ijet));
      }
    } //close jet loop
  }
  ct->ResetBranchAddresses();

  //TLatex *twrite[ptBins][djetBins][cTypes];
  TLatex *twrite2[cTypes];
  string dtypes[cTypes] = {"Dstar","D0","Ds->PhiPi","Ds->K*K","D+/-"};

  if(doDraw){
    TCanvas *cc[djetBins][cTypes];
    for(int k=0; k<cTypes; k++){
      for(int j=0; j<djetBins; j++){
	cc[j][k] = new TCanvas(Form("cc_djet%d_type%d",j,k),Form("cc_djet%d_type%d",j,k),1000,800);
	cc[j][k]->Divide(3,3);
	for(int i=0; i<ptBins; i++){
	  cc[j][k]->cd(i+1);
	  assnsUnweight[i][j][k]->Draw();
	  cc[j][k]->Update();
	  if(i==0) twrite2[i] = new TLatex(1.75,gPad->GetUymax()*0.85,Form("jet p_{T}: %d-%d, djet: %.1f-%.1f dType: %s",dCandptBoundaries[i],dCandptBoundaries[i+1], djetBoundaries[j],djetBoundaries[j+1], dtypes[k].c_str()));
	  if(i==0) twrite2[i]->Draw("same");
	  //twrite[i][j][k] = new TLatex(3,gPad->GetUymax()*0.75, Form("#Deltajet: %.1f-%.1f",djetBoundaries[j],djetBoundaries[j+1]));
	  //twrite[i][j][k]->Draw("same");
	}
      }
    }
  }
  
  fout->cd();
  for(int i=0; i<jetPtBins; i++){
    for(int j=0; j<djetBins; j++){
      for(int k=0; k<cTypes; k++){
      assns[i][j][k]->Write();
      assnsUnweight[i][j][k]->Write();
      bgAssns[i][j][k]->Write();
      bgAssnsUnweight[i][j][k]->Write();
      }
    }
  }

  for(int i=0; i<jetPtBins; i++){
    for(int k=0; k<cTypes; k++){
      proximity[i][k]->Write();
      proximityEta[i][k]->Write();
      proximityPhi[i][k]->Write();
      bgProximity[i][k]->Write();
      bgProximityEta[i][k]->Write();
      bgProximityPhi[i][k]->Write();
    }
  }
  allJetSpec->Write();
  for(int k=0; k<cTypes; k++){
    CjetSpec[k]->Write();
    bgCjetSpec[k]->Write();
  }

  //cJetEffvsPt->Divide(cJetEffvsPt,cJetEffvsPt_bg,1,1,"B");
  //cJetPurvsPt->Divide(cJetPurvsPt,cJetPurvsPt_bg,1,1,"B");
  cJetEffvsPt->SetYTitle( (Dlabels[isMC-1]+(string)Form("dR<%g C-Jet Efficiency",djetR_testCut)).c_str());
  cJetEffvsPt->SetXTitle("c-Jet p_{T}");
  cJetPurvsPt->SetYTitle( (Dlabels[isMC-1]+(string)Form("dR<%g C-Jet Purity",djetR_testCut)).c_str());
  cJetPurvsPt->SetXTitle("c-Jet p_{T}");
  cJetEffvsPt->Write();
  cJetEffvsPt_bg->Write();
  cJetEffvsPt_rej->Write();
  cJetPurvsPt->Write();

  //cJetEffvsPt_take2->Divide(cJetEffvsPt_take2,cJetEffvsPt_bg_take2,1,1,"B");
  cJetEffvsPt_step1->Write();
  cJetEffvsPt_bg_step1->Write();
  cJetEffvsPt_rej_step1->Write();
  cJetEffvsPt_step2->Write();
  cJetEffvsPt_bg_step2->Write();
  cJetEffvsPt_rej_step2->Write();
  cJetEffvsPt_step3[3] = (TH1D*)cJetEffvsPt_step3[2]->Clone("cJetEffvsPt_step3_type3");
  for(int i=0; i<5; i++){
    cJetEffvsPt_step3[i]->Write();
  }
  cJetEffvsPt_bg_step3->Write();
  cJetEffvsPt_rej_step3->Write();

  jetAssnTop->Write();
  jetAssnBot->Write();

  dRTestHist->Write();
  unmatchedCJets->Write();

  Color_t colors[ptBins] = {kBlack,kRed,kOrange+1,kGreen+2,kBlue+2,kMagenta+2};
  for(int i=0; i<ptBins; i++){
    jetFF[i]->Scale(1./jetFF[i]->Integral());
    jetFF[i]->SetMarkerColor(colors[i]);
    jetFF[i]->Write();
  }

  fout->Close();
}
