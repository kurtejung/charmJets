#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>
#include <set>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TF1.h"

#ifdef __CINT__
#pragma link C++ class vector<vector<int> >+;
#endif

using namespace std;

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

void makeCharmJetAssn_CJetSample(int isMC=6, bool doDraw = 0){

  const int ptBins = 6;
  const int jetPtBins = 6;
  const int djetBins = 1;

  const bool fillDMassWithoutJets = false;
  const bool findRecoEffVsDMesonPt = true;

  const double djetR_testCut = 0.3;

  int ptBoundaries[jetPtBins+1] = {20,40,60,80,100,200,400};
  int dCandptBoundaries[ptBins+1] = {0,10,15,20,40,60,100};
  double djetBoundaries[djetBins+1] = {0,9999};
  string Dlabels[6] = {"D* ","D_{0} ","D_{s}->#pi#phi ","D_{s}->K*K ","D^{+/-} ","C-Jet Only "};

  // type { 1,2,3,4,5,?}
  //d*, d0, ds->phipi, ds->k*k, dpm
  int pdgs[6] = {413, 421, 431, 431, 411, -999};

  //TFile *fpur = new TFile("purityFactors_MC_redux.root","recreate");
  TFile *fout;
  if(!isMC) fout = new TFile("dCandidateToJetAssns_NoJetTrg_AllData_RevCuts_JetAssn-9-12.root","recreate");
  else if(isMC==1) fout = new TFile("dCandidateToJetAssns_Dstar_RevCuts_JetAssn-8-11.root","recreate");
  else if(isMC==2) fout = new TFile("dCandidateToJetAssns_D02Kpi_RevCuts_JetAssn-8-11.root","recreate");
  else if(isMC==3) fout = new TFile("dCandidateToJetAssns_Ds2phipi_RevCuts_JetAssn-8-11.root","recreate");
  else if(isMC==4) fout = new TFile("dCandidateToJetAssns_Ds2KstarK_RevCuts_JetAssn-8-11.root","recreate");
  else if(isMC==5) fout = new TFile("dCandidateToJetAssns_Dpm_RevCuts_JetAssn-8-11.root","recreate");
  else if(isMC==6) fout = new TFile("dCandidateToJetAssns_CJetOnly_RevCuts_JetAssn-10-30.root","recreate");
  else{
    cout << "isMC cannot be more than 6! Is this a new D-decay type??" << endl;
    exit(0);
  }

  const double tightSignalCutLow[cTypes] = {0.144, 1.84, 1.93, 1.93, 1.85};
  const double tightSignalCutHi[cTypes] = {0.147, 1.885, 1.99, 1.99, 1.885};

  const int jetBins = 10;
  double xjetBins[jetBins+1] = {0,5,10,15,20,40,60,80,100,200,400};

  const int dmesonBins = 12;
  double xmesonBins[dmesonBins+1] = {0,4,8,12,15,18,20,25,30,40,60,80,100};

  const int ndrbins = 20;
  double xdrbins[ndrbins+1] = {-10,0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.3,0.4,0.6,1.,1.5,2,2.5,3,3.5,4,20};

  TH1D *allJetSpec = new TH1D("allJetSpec","",jetBins,xjetBins); allJetSpec->Sumw2();
  TH2D *CjetSpec[cTypes];
  TH2D *bgCjetSpec[cTypes];
  TH1D *purity_bg[cTypes];
  TH1D *purity[cTypes];

  TH1D *proximity[ptBins][cTypes];
  TH1D *proximityEta[ptBins][cTypes];
  TH1D *proximityPhi[ptBins][cTypes];
  TH1D *bgProximity[ptBins][cTypes];
  TH1D *bgProximityEta[ptBins][cTypes];
  TH1D *bgProximityPhi[ptBins][cTypes];
  TH1D *assns[ptBins][djetBins][cTypes];
  TH1D *assnsUnweight[ptBins][djetBins][cTypes];
  TH1D *bgAssns[ptBins][djetBins][cTypes];
  TH1D *bgAssnsUnweight[ptBins][djetBins][cTypes];
  TH1D *dcAssns[ptBins][djetBins][cTypes];
  TH1D *nondcAssns[ptBins][djetBins][cTypes];

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

  TH1D *jetFF[ptBins];
  for(int i=0; i<ptBins; i++){
    jetFF[i] = new TH1D(Form("jetFF_%d",i),"",30,0,1.5); jetFF[i]->Sumw2();
  }

  TH1D *jetDRMatched = new TH1D("jetDRMatched","",30,0,0.3); jetDRMatched->Sumw2();
  TH1D *jetDR = new TH1D("jetDR","",30,0,0.3); jetDR->Sumw2();
  TH2D *jetCorr = new TH2D("jetCorr","",30,0,200,30,0,400); jetCorr->Sumw2();

  TH2D *cMesonRecoEff[5];
  TH2D *cMesonRecoEff_bg[5];
  for(int i=0; i<5; i++){
    cMesonRecoEff[i] = new TH2D(Form("cMesonRecoEff_%d",i),"",20,0,50,jetBins,xjetBins); cMesonRecoEff[i]->Sumw2();
    cMesonRecoEff_bg[i] = new TH2D(Form("cMesonRecoEff_bg_%d",i),"",20,0,50,jetBins,xjetBins); cMesonRecoEff_bg[i]->Sumw2();
  }
  TH1D *cJetEffvsPt_step1[5];
  TH1D *cJetEffvsPt_step2[5];
  TH1D *cJetEffvsPt_step3[5]; //one for each d-meson type
  TH1D *cJetEffvsPt_step4[5];

  TH1D *cJetEffvsPt_bg_step1[5];
  TH1D *cJetEffvsPt_bg_step2[5];
  TH1D *cJetEffvsPt_bg_step3[5];
  TH1D *cJetEffvsPt_bg_step4[5];

  TH1D *cJetEffvsPt_rej_step1[5];
  TH1D *cJetEffvsPt_rej_step2[5];
  TH1D *cJetEffvsPt_rej_step3[5];
  TH1D *cJetEffvsPt_rej_step4[5];

  TH1D *cMesonRecoEfftest[5];
  TH1D *cMesonRecoEfftest_bg[5];

  for(int i=0; i<5; i++){
    cMesonRecoEfftest[i] = new TH1D(Form("cMesonRecoEfftest_%d",i),"",dmesonBins,xmesonBins); cMesonRecoEfftest[i]->Sumw2();
    cMesonRecoEfftest_bg[i] = new TH1D(Form("cMesonRecoEfftest_bg_%d",i),"",dmesonBins,xmesonBins); cMesonRecoEfftest_bg[i]->Sumw2();
  }

  for(int i=0; i<5; i++){
    cJetEffvsPt_step1[i] = new TH1D(Form("cJetEffvsPt_step1_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_step1[i]->Sumw2();
    cJetEffvsPt_step2[i] = new TH1D(Form("cJetEffvsPt_step2_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_step2[i]->Sumw2();
    cJetEffvsPt_step3[i] = new TH1D(Form("cJetEffvsPt_step3_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_step3[i]->Sumw2();
    cJetEffvsPt_step4[i] = new TH1D(Form("cJetEffvsPt_step4_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_step4[i]->Sumw2();

    cJetEffvsPt_bg_step1[i] = new TH1D(Form("cJetEffvsPt_bg_step1_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_bg_step1[i]->Sumw2();
    cJetEffvsPt_bg_step2[i] = new TH1D(Form("cJetEffvsPt_bg_step2_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_bg_step2[i]->Sumw2();
    cJetEffvsPt_bg_step3[i] = new TH1D(Form("cJetEffvsPt_bg_step3_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_bg_step3[i]->Sumw2();
    cJetEffvsPt_bg_step4[i] = new TH1D(Form("cJetEffvsPt_bg_step4_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_bg_step4[i]->Sumw2();

    cJetEffvsPt_rej_step1[i] = new TH1D(Form("cJetEffvsPt_rej_step1_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_rej_step1[i]->Sumw2();
    cJetEffvsPt_rej_step2[i] = new TH1D(Form("cJetEffvsPt_rej_step2_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_rej_step2[i]->Sumw2();
    cJetEffvsPt_rej_step3[i] = new TH1D(Form("cJetEffvsPt_rej_step3_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_rej_step3[i]->Sumw2();
    cJetEffvsPt_rej_step4[i] = new TH1D(Form("cJetEffvsPt_rej_step4_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_rej_step4[i]->Sumw2();
  }

  for(int i=0; i<ptBins; i++){
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
	//bgAssns[i][j][k]->Sumw2();
	bgAssns[i][j][k]->SetXTitle("Candidate Mass [GeV/c^{2}]");
	bgAssns[i][j][k]->SetYTitle("Weighted Counts");

	bgAssnsUnweight[i][j][k] = new TH1D(Form("bgAssns_unweight_pt%d_djet%d_type%d",i,j,k),"",35,lowMassBin[k],highMassBin[k]);
	//bgAssnsUnweight[i][j][k]->Sumw2();
	bgAssnsUnweight[i][j][k]->SetXTitle("Candidate Mass [GeV/c^{2}]");
	bgAssnsUnweight[i][j][k]->SetYTitle("Unweighted Counts");

	dcAssns[i][j][k] = new TH1D(Form("DC_assns_pt%d_djet%d_type%d",i,j,k),"",35,lowMassBin[k],highMassBin[k]);
	dcAssns[i][j][k]->Sumw2();
	dcAssns[i][j][k]->SetXTitle("Candidate Mass [GeV/c^{2}]");
	dcAssns[i][j][k]->SetYTitle("Unweighted Counts");

	nondcAssns[i][j][k] = new TH1D(Form("non_DC_assns_pt%d_djet%d_type%d",i,j,k),"",35,lowMassBin[k],highMassBin[k]);
	nondcAssns[i][j][k]->Sumw2();
	nondcAssns[i][j][k]->SetXTitle("Candidate Mass [GeV/c^{2}]");
	nondcAssns[i][j][k]->SetYTitle("Unweighted Counts");
      }
    }
  }
  
 
  for(int i=0; i<ptBins; i++){
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
    CjetSpec[k] = new TH2D(Form("CjetSpec_type%d",k),"",jetBins,xjetBins,dmesonBins,xmesonBins); CjetSpec[k]->Sumw2();
    bgCjetSpec[k] = new TH2D(Form("bgCjetSpec_type%d",k),"",jetBins,xjetBins,dmesonBins,xmesonBins); bgCjetSpec[k]->Sumw2();
    purity[k] = new TH1D(Form("purity_type%d",k),"",jetBins,xjetBins); purity[k]->Sumw2();
    purity_bg[k] = new TH1D(Form("purity_bg_type%d",k),"",jetBins,xjetBins); purity_bg[k]->Sumw2();
  }

  double weight;
  int HLT_Jet20_NoJetID_v1, HLT_Jet40_NoJetID_v1, HLT_Jet60_NoJetID_v1, HLT_Jet80_NoJetID_v1, HLT_Jet100_NoJetID_v1;
  vector<double> *dCandPt=0, *dCandMass=0, *dCandEta=0, *dCandPhi=0, *jtpt=0, *jteta=0, *jtphi=0, *refpt=0, *rawpt=0;
  vector<double> *dCandType=0, *dCandChildMass=0, *dCandCharge1=0, *dCandCharge2=0;
  vector<double> *dCandMatchGenPdg_dau2=0, *dCandMatchGenPdg_index1=0, *dCandMatchGenPdg_index2=0;
  vector<double> *refparton_flavorForB=0;
  vector<vector<int> > *dCandParentPartonIdx=0;
  vector<int> *genMatch=0, *dCandGenMatchPdg=0, *dCandGenMatchDoubleCountedPdg=0, *dCandGenMatchIdx=0, *dCandGenMatchDoubleCountedIdx=0;
  TBranch *bdcandpt=0, *bdcandmass=0, *bdcandeta=0, *bdcandphi=0, *bdcandtype=0, *bdcandchildmass, *bjtpt=0, *bjteta=0, *bjtphi=0, *bdcandcharge1=0, *bdcandcharge2=0;
  vector<vector<int> > *hasGenD=0, *hasGenDwithKpi=0, *hasGenDwithKpiIdx=0;

  //d*, d0, ds->phipi, ds->k*k, dpm
  TFile *fin;
  if(!isMC) fin = new TFile("input/DMesonCJet_RevCuts_NoJetTrgCut_pPbdata_ppReco_akPu3PF.root","OLD");
  if(isMC==1) fin = new TFile("DMesonCJet_Dstar_pPbMC_ppReco_akPu3PF_1.root","OLD");
  if(isMC==2) fin = new TFile("DMesonCJet_Dzero_pPbMC_ppReco_akPu3PF_1.root","OLD");
  if(isMC==3) fin = new TFile("DMesonCJet_Ds2PhiPi_pPbMC_ppReco_akPu3PF_76.root","OLD");
  if(isMC==4) fin = new TFile("DMesonCJet_Ds2KstarKEmbed_pPbMC_ppReco_akPu3PF_1.root","OLD");
  if(isMC==5) fin = new TFile("DMesonCJet_DpmEmbed_pPbMC_ppReco_akPu3PF_1.root","OLD");
  if(isMC==6) fin = new TFile("input/DMesonCJet_CJetOnly_pPbMC_ppReco_akPu3PF_512.root","OLD");
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
    ct->SetBranchAddress("hasGenDwithKpiIdx",&hasGenDwithKpiIdx);
  }
  ct->SetBranchAddress("weight",&weight);
  ct->SetBranchAddress("HLT_Jet20_NoJetID_v1",&HLT_Jet20_NoJetID_v1);
  ct->SetBranchAddress("HLT_Jet40_NoJetID_v1",&HLT_Jet40_NoJetID_v1);
  ct->SetBranchAddress("HLT_Jet60_NoJetID_v1",&HLT_Jet60_NoJetID_v1);
  ct->SetBranchAddress("HLT_Jet80_NoJetID_v1",&HLT_Jet80_NoJetID_v1);
  ct->SetBranchAddress("HLT_Jet100_NoJetID_v1",&HLT_Jet100_NoJetID_v1);

  TTree *gt;
  vector<double> *genPt=0, *genEta=0, *genPhi=0;
  vector<int> *genIdx=0, *genPdg=0, *genPtIdx2=0;
  vector<bool> *hasHadronicDecay=0;
  if(isMC){
    gt = (TTree*)fin->Get("gTree");
    gt->SetBranchAddress("genPt",&genPt);
    gt->SetBranchAddress("genEta",&genEta);
    gt->SetBranchAddress("genPhi",&genPhi);
    gt->SetBranchAddress("genIdx",&genIdx); //event-level index of particle
    gt->SetBranchAddress("genPdg",&genPdg);
    gt->SetBranchAddress("genPtIdx2",&genPtIdx2);
    gt->SetBranchAddress("hasHadronicDecay",&hasHadronicDecay);
  }

  TF1 *feffline[5];
  for(int ij=0; ij<5; ij++){
    feffline[ij] = new TF1(Form("feffline_%d",ij),"[0]*exp(-[1]/x)",0,100);
    feffline[ij]->SetParameter(0,0.5);
    feffline[ij]->SetParameter(1,2);
  }

  int nentries = ct->GetEntries();
  if(findRecoEffVsDMesonPt){
    for(int ientry=0; ientry<nentries; ientry++){
      ct->GetEntry(ientry);
      if(isMC) gt->GetEntry(ientry);
      //for all gen-level d-mesons, how many are reconstructed?
      //do reconstruction efficiency starting from the gen-level d-meson jet association
      for(unsigned int ijet=0; ijet<jtpt->size(); ijet++){
	if(refpt->at(ijet)<20) continue;
	if(rawpt->at(ijet)<23) continue;
	if(TMath::Abs(refparton_flavorForB->at(ijet))==4){
	  for(int ij=0; ij<5; ij++){
	    bool FF2=false;
	    bool foundFlag=false;
	    double genPtToFill=0;
	    for(unsigned int igen=0; igen<genPt->size(); igen++){
	      if(hasGenDwithKpiIdx->at(ijet)[ij]==genIdx->at(igen) /*&& hasHadronicDecay->at(igen)*/ && abs(genPdg->at(igen))==pdgs[ij]){
		FF2 = true;
		genPtToFill = genPt->at(igen);
		for(unsigned int icand=0; icand<dCandPt->size(); icand++){
		  if((dCandGenMatchIdx->at(icand)==genIdx->at(igen) || dCandGenMatchDoubleCountedIdx->at(icand)==genIdx->at(igen)) && !foundFlag && dCandType->at(icand)-1 == ij ){
		    foundFlag=true;
		  }
		}
	      }
	    }
	    if(FF2) cJetEffvsPt_bg_step2[ij]->Fill(genPtToFill);
	    if(foundFlag){
	      cJetEffvsPt_step2[ij]->Fill(genPtToFill);
	    }
	    else{
	      cJetEffvsPt_rej_step2[ij]->Fill(genPtToFill);
	    }
	  }
	}
      }
    }
    for(int ij=0; ij<5; ij++){
      cJetEffvsPt_step2[ij]->Divide(cJetEffvsPt_step2[ij],cJetEffvsPt_bg_step2[ij],1,1,"B");
      cJetEffvsPt_step2[ij]->Fit(feffline[ij],"q0","",0,100);
    }
  }

  for(int ientry=0; ientry<nentries; ientry++){
    ct->GetEntry(ientry);
    if(isMC) gt->GetEntry(ientry);
    if (ientry%1000000==0) cout<<" i = "<<ientry<<" out of "<<nentries<<" ("<<(int)(100*(float)ientry/(float)nentries)<<"%)"<<endl;

    //if(isMC) weight = 1.;
    if(dCandPt->size() != dCandMass->size() || dCandPt->size() != dCandEta->size() || dCandPt->size() != dCandPhi->size() || dCandPt->size() != dCandParentPartonIdx->size()){
      cout << "Warning! dCandidate sizes are different in entry: "<< ientry << "!" << endl;
      cout << "dCandPt size: "<< dCandPt->size() << " dcandparentSize: "<< dCandParentPartonIdx->size() << endl;
      cout << "Skipping event..." << endl;
      continue;
    }
    
    if(!HLT_Jet20_NoJetID_v1 && !HLT_Jet40_NoJetID_v1 && !HLT_Jet60_NoJetID_v1 && !HLT_Jet80_NoJetID_v1 && !HLT_Jet100_NoJetID_v1) continue;

    if(dCandPt->size() != dCandParentPartonIdx->size()) cout << "warning! candidate parent parton != dcandPt!" << endl;
    //if(dCandPt->size() != hasGenD->size()) cout << "warning! d-match flag != dCandPt!" << endl;
    //cout << "dcandpt size: "<< dCandPt->size() << " hasgend size: "<< hasGenD->size() << endl;

    for(unsigned int igen=0; igen<genPt->size(); igen++){
      for(unsigned int ijet=0; ijet<jtpt->size(); ijet++){
	if(refpt->at(ijet)<20) continue;
	if(rawpt->at(ijet)<23) continue;
	//if(genMatch->at(ijet)==genIdx->at(igen)){
	double deta = fabs(genEta->at(igen)-jteta->at(ijet));
	double dphi = fabs(genPhi->at(igen)-jtphi->at(ijet));
	double djetR = sqrt(pow(deta,2)+pow(dphi,2));
	if(djetR<0.3){
	  for(int ij=0; ij<5; ij++){
	    if(abs(genPdg->at(igen))==pdgs[ij]){
	      cMesonRecoEff_bg[ij]->Fill(genPt->at(igen),refpt->at(ijet));
	      //cJetEffvsPt_bg_step2[ij]->Fill(refpt->at(ijet)); //fill directly from gen particles (doesnt work for Ds or D*)
	    }
	  }
	}
      }
    }

    for(unsigned int icand=0; icand<dCandPt->size(); icand++){
      for(unsigned int ijet=0; ijet<jtpt->size(); ijet++){
	if(refpt->at(ijet)<20) continue;
	if(rawpt->at(ijet)<23) continue;
	for(unsigned int icandParton=0; icandParton<dCandParentPartonIdx->at(icand).size(); icandParton++){
	  if(genMatch->at(ijet) == dCandParentPartonIdx->at(icand).at(icandParton)){
	    for(int ij=0; ij<5; ij++){
	      if(abs(dCandGenMatchPdg->at(icand))==pdgs[ij]){
		cMesonRecoEff[ij]->Fill(dCandPt->at(icand),refpt->at(ijet));
	      }
	    }
	  }
	}
      }
    }

    //do efficiency for dR matching (step 1)
    for(unsigned int ijet=0; ijet<jtpt->size(); ijet++){
      if(refpt->at(ijet)<20) continue;
      if(rawpt->at(ijet)<23) continue;
      for(int ij=0; ij<5; ij++){
	bool FF1=false;
	bool FF2=false;
	for(unsigned int igen=0; igen<genPt->size(); igen++){
	  for(unsigned int icand=0; icand<dCandPt->size(); icand++){
	    if(dCandGenMatchIdx->at(icand)==genIdx->at(igen) && dCandType->at(icand)-1 == ij && hasGenDwithKpiIdx->at(ijet)[ij]==genIdx->at(igen)){
	      FF1 = true;
	      double djetR = sqrt(pow((dCandEta->at(icand)-jteta->at(ijet)),2)+pow((dCandPhi->at(icand)-jtphi->at(ijet)),2));
	      if(djetR<djetR_testCut){
		FF2=true;
	      }
	    }
	  }
	}
	if(FF1) cJetEffvsPt_bg_step1[ij]->Fill(refpt->at(ijet));
	if(FF2) cJetEffvsPt_step1[ij]->Fill(refpt->at(ijet));
	else cJetEffvsPt_rej_step1[ij]->Fill(refpt->at(ijet));	      
      }
    }


    //for all gen-level d-mesons, how many are reconstructed?
    //do reconstruction efficiency starting from the gen-level d-meson jet association
    for(unsigned int ijet=0; ijet<jtpt->size(); ijet++){
      if(refpt->at(ijet)<20) continue;
      if(rawpt->at(ijet)<23) continue;
      if(TMath::Abs(refparton_flavorForB->at(ijet))==4){
	for(int ij=0; ij<5; ij++){
	  bool FF2=false;
	  bool foundFlag=false;
	  for(unsigned int igen=0; igen<genPt->size(); igen++){
	    if(hasGenDwithKpiIdx->at(ijet)[ij]==genIdx->at(igen) /*&& hasHadronicDecay->at(igen)*/ && abs(genPdg->at(igen))==pdgs[ij]){
	      FF2 = true;
	      for(unsigned int icand=0; icand<dCandPt->size(); icand++){
		if((dCandGenMatchIdx->at(icand)==genIdx->at(igen) || dCandGenMatchDoubleCountedIdx->at(icand)==genIdx->at(igen)) && !foundFlag && dCandType->at(icand)-1 == ij ){
		  foundFlag=true;
		}
	      }
	    }
	    if(!findRecoEffVsDMesonPt){
	      if(FF2) cJetEffvsPt_bg_step2[ij]->Fill(refpt->at(ijet));
	      if(foundFlag){
		cJetEffvsPt_step2[ij]->Fill(refpt->at(ijet));
	      }
	      else{
		cJetEffvsPt_rej_step2[ij]->Fill(refpt->at(ijet));
	      }
	    }
	  }
	}
      }
    }
      
    //check d-meson reco eff as a function of d-meson pt
    for(int ij=0; ij<5; ij++){
      for(unsigned int igen=0; igen<genPt->size(); igen++){
	bool FF2=false;
	bool foundFlag=false;
	for(unsigned int ijet=0; ijet<jtpt->size(); ijet++){
	  if(hasGenDwithKpiIdx->at(ijet)[ij]==genIdx->at(igen) && abs(genPdg->at(igen))==pdgs[ij]){ //require D->Kpi decay!!
	    FF2=true;
	    for(unsigned int icand=0; icand<dCandPt->size(); icand++){
	      if(dCandGenMatchIdx->at(icand)==genIdx->at(igen)){
		foundFlag=true;
	      }
	    }
	  }
	}
	if(FF2) cMesonRecoEfftest_bg[ij]->Fill(genPtIdx2->at(igen));
	if(foundFlag) cMesonRecoEfftest[ij]->Fill(genPtIdx2->at(igen));
      }
    }
    /*if(FF2) cMesonRecoEfftest_bg[ij]->Fill(dCandPtTemp);
      if(foundFlag){
      cMesonRecoEfftest[ij]->Fill(dCandPtTemp);
      }*/
	    
    
    //for all events with d-mesons, how many have an associated jet?
    double maxcandpt=0;
    if(jtpt->size()>0){
      for(unsigned int icand=0; icand<jtpt->size(); icand++){
	if(jtpt->at(icand) > maxcandpt) maxcandpt = jtpt->at(icand);
      }
      jetAssnBot->Fill(maxcandpt);
      if(dCandPt->size()>0){
	jetAssnTop->Fill(maxcandpt);
      }
    }
    //check out d-mesons without jet correspondance...
    if(fillDMassWithoutJets){
      for(unsigned int icand=0; icand<dCandPt->size(); icand++){
	candidate tempCand;
	int type = dCandType->at(icand)-1;
        tempCand.setStuff(dCandMass->at(icand), dCandPt->at(icand), dCandChildMass->at(icand), dCandEta->at(icand), dCandPhi->at(icand), type, dCandCharge1->at(icand), dCandCharge2->at(icand),-999,-999,-999);
	int iCandpt=0;
	while(dCandPt->at(icand) > dCandptBoundaries[iCandpt+1] && iCandpt<(ptBins-1)) iCandpt++;
	//if(abs(dCandGenMatchPdg->at(icand))==pdgs[type]){
	if(extraDCuts(tempCand,1)){
	  if(dCandType->at(icand)==1){ //start with D* meson since its different
	    assns[iCandpt][0][type]->Fill(dCandMass->at(icand) - dCandChildMass->at(icand), weight);
	    assnsUnweight[iCandpt][0][type]->Fill(dCandMass->at(icand)-dCandChildMass->at(icand));
	    bool foundFlag = false;
	    for(unsigned int igen=0; igen<genPt->size(); igen++){
	      if(genIdx->at(igen) == dCandGenMatchIdx->at(icand) && abs(genPdg->at(igen))==pdgs[type]){
		dcAssns[iCandpt][0][type]->Fill(dCandMass->at(icand)-dCandChildMass->at(icand));
		foundFlag = true;
	      }
	    }
	    if(!foundFlag){
	      nondcAssns[iCandpt][0][type]->Fill(dCandMass->at(icand)-dCandChildMass->at(icand));
	    }
	  }
	  else{ //do all other d-mesons now
	    //if(abs(dCandGenMatchPdg->at(icand)) == pdgs[type]){
	    assns[iCandpt][0][type]->Fill(dCandMass->at(icand), weight);
	    assnsUnweight[iCandpt][0][type]->Fill(dCandMass->at(icand));
	    //}
	    bool foundFlag = false;
	    for(unsigned int igen=0; igen<genPt->size(); igen++){
	      if(genIdx->at(igen) == dCandGenMatchIdx->at(icand) && abs(genPdg->at(igen))==pdgs[type]){
		dcAssns[iCandpt][0][type]->Fill(dCandMass->at(icand));
		foundFlag = true;
	      }
	    }
	    if(!foundFlag){
	      nondcAssns[iCandpt][0][type]->Fill(dCandMass->at(icand));
	    }
	  }
	}
	else if(extraDCuts(tempCand,0)){
	  if(dCandType->at(icand)==1){
	    bgAssnsUnweight[iCandpt][0][type]->Fill((dCandMass->at(icand) - dCandChildMass->at(icand)));
	    bgAssns[iCandpt][0][type]->Fill(dCandMass->at(icand) - dCandChildMass->at(icand),weight);
	  }
	  else{
	    bgAssnsUnweight[iCandpt][0][type]->Fill((dCandMass->at(icand)));
	    bgAssns[iCandpt][0][type]->Fill(dCandMass->at(icand),weight);
	  }
	}
      }
    }
    
    for(unsigned int ijet=0; ijet<jtpt->size(); ijet++){
      if(refpt->at(ijet)<20) continue;
      if(rawpt->at(ijet)<23) continue;
      //if(ijet>1) continue;
      allJetSpec->Fill(jtpt->at(ijet));
      bool dRFlag = false;
      bool dRFlag2 = false;
      bool anotherTestFlag=false;
      bool checkDCjets[5] = {0,0,0,0,0};
      bool checkbgDCjets[5] = {0,0,0,0,0};
      double minDR=999;
      int candJetMatchIdx=-1;
      int candJetPartonMatchIdx=-1;
      for(int ij=0; ij<5; ij++){
	if((hasGenDwithKpi->at(ijet)[ij]) && TMath::Abs(refparton_flavorForB->at(ijet))==4){
	  //cJetEffvsPt_bg_step2[ij]->Fill(jtpt->at(ijet));
	}
      }

      for(unsigned int icand=0; icand<dCandPt->size(); icand++){
	double deta = fabs(dCandEta->at(icand)-jteta->at(ijet));
	double dphi = fabs(dCandPhi->at(icand)-jtphi->at(ijet));
	double djetR = sqrt(pow(deta,2)+pow(dphi,2));
	int type = dCandType->at(icand)-1;
	if(djetR<djetR_testCut){
	  purity_bg[type]->Fill(jtpt->at(ijet));
	  if(TMath::Abs(refparton_flavorForB->at(ijet))==4){
	    purity[type]->Fill(jtpt->at(ijet));
	  }
	}
      }

      double maxcandpt2 = -1;
      int maxcandno = -1;
      double candDR[5]={999,999,999,999,999};
      double bgcandDR[5]={999,999,999,999,999};
      double assocPt[5]={0,0,0,0,0};
      for(unsigned int icand=0; icand<dCandPt->size(); icand++){
	if(dCandPt->at(icand)>maxcandpt2){
	  maxcandpt2 = dCandPt->at(icand);
	  maxcandno = icand;
	}
      }
      //if(maxcandno>0){
      for(unsigned int icand=0; icand<dCandPt->size(); icand++){
	//int icand = maxcandno; //------------------------flip this switch on to just look at the highest pt d-meson
	//if(dCandPt->at(icand) > maxcandpt) continue;
	int ipt=0;
	int iCandpt=0;
	//int idjet=0;
	
	//double jetprox = findClosestJet(dCandEta, dCandPhi, jteta->at(ijet), jtphi->at(ijet));
	//if(jetprox==-1) continue;
	double deta = fabs(dCandEta->at(icand)-jteta->at(ijet));
	double dphi = fabs(dCandPhi->at(icand)-jtphi->at(ijet));
	double djetR = sqrt(pow(deta,2)+pow(dphi,2));
	
	if(isMC){

	  for(unsigned int icandParton=0; icandParton<dCandParentPartonIdx->at(icand).size(); icandParton++){
	    if(genMatch->at(ijet) == dCandParentPartonIdx->at(icand).at(icandParton)) anotherTestFlag=true;
	  }
	  
	  //check all partons to make sure *one of them* is the D-Meson in the jet cone
	  //if(abs(dCandGenMatchPdg->at(icand))==pdgs[isMC-1]){
	    for(unsigned int icandParton=0; icandParton<dCandParentPartonIdx->at(icand).size(); icandParton++){
	      if(genMatch->at(ijet) == dCandParentPartonIdx->at(icand).at(icandParton)){// && TMath::Abs(refparton_flavorForB->at(ijet))==4){
		if(djetR<minDR){
		  minDR = djetR;
		  candJetMatchIdx=icand;
		  candJetPartonMatchIdx=icandParton;
		}
	      }
	    }
	    //}
	  for(unsigned int icandParton=0; icandParton<dCandParentPartonIdx->at(icand).size(); icandParton++){
	    if(genMatch->at(ijet) == dCandParentPartonIdx->at(icand).at(icandParton) && dRFlag==false){
	      dRTestHist->Fill(djetR);
	      cJetEffvsPt_bg->Fill(refpt->at(ijet));
	      if(djetR<djetR_testCut){ 
		cJetEffvsPt->Fill(refpt->at(ijet));
		cJetPurvsPt->Fill(refpt->at(ijet));
	      }
	      else{
		cJetEffvsPt_rej->Fill(refpt->at(ijet));
	      }
	      dRFlag=true;
	    }
	    if(djetR<djetR_testCut && dRFlag2==false){
	      cJetPurvsPt_bg->Fill(refpt->at(ijet));
	      dRFlag2=true;
	    }
	  }
	}
	//while(djetR>djetBoundaries[idjet+1] && idjet<(djetBins-1)) idjet++;     
	while(jtpt->at(ijet) > ptBoundaries[ipt+1] && ipt<(jetPtBins-1)) ipt++;
	while(dCandPt->at(icand) > dCandptBoundaries[iCandpt+1] && iCandpt<(ptBins-1)) iCandpt++;
	if(iCandpt > jetPtBins-1) iCandpt=jetPtBins-1;
	int type = dCandType->at(icand)-1;
	//	if(ipt>5 || iCandpt>4 || type>4) cout << "ipt: "<< ipt << " icandpt: "<< iCandpt << " type: "<< type << endl;
	//ipt = iCandpt;
	candidate tempCand;
	//if(!isMC) 
        tempCand.setStuff(dCandMass->at(icand), dCandPt->at(icand), dCandChildMass->at(icand), dCandEta->at(icand), dCandPhi->at(icand), type, dCandCharge1->at(icand), dCandCharge2->at(icand),-999,-999,-999);
	//else tempCand.setStuff(dCandMass->at(icand), dCandPt->at(icand), dCandChildMass->at(icand), dCandEta->at(icand), dCandPhi->at(icand), type, dCandCharge1->at(icand), dCandCharge2->at(icand),dCandMatchGenPdg_dau2->at(icand),dCandMatchGenPdg_index1->at(icand),dCandMatchGenPdg_index2->at(icand));
	if(extraDCuts(tempCand,1)){
	  if(dCandType->at(icand)==1){ //start with D* meson since its different
	    if(!fillDMassWithoutJets){
	      if(djetR<djetR_testCut){
		assns[ipt][0][type]->Fill(dCandMass->at(icand) - dCandChildMass->at(icand), weight);
		assnsUnweight[ipt][0][type]->Fill(dCandMass->at(icand)-dCandChildMass->at(icand));
		if(abs(dCandGenMatchDoubleCountedPdg->at(icand)) == pdgs[0]){
		  dcAssns[ipt][0][type]->Fill(dCandMass->at(icand)-dCandChildMass->at(icand));
		}
	      }
	    }
	    double lowball = snapToBin(assns[ipt][0][type],tightSignalCutLow[type],0);
	    double hiball = snapToBin(assns[ipt][0][type],tightSignalCutHi[type],1);
	    if((dCandMass->at(icand)-dCandChildMass->at(icand)) > lowball && (dCandMass->at(icand)-dCandChildMass->at(icand)) < hiball){
	      checkDCjets[type]=1;
	      proximity[ipt][type]->Fill(djetR);
	      proximityEta[ipt][type]->Fill(deta);
	      proximityPhi[ipt][type]->Fill(dphi);
	      jetDR->Fill(djetR);
	      if(djetR<djetR_testCut){
		//CjetSpec[type]->Fill(jtpt->at(ijet));
		for(unsigned int igen=0; igen<genPt->size(); igen++){
		  for(int ij=0; ij<5; ij++){
		    if(hasGenDwithKpiIdx->at(ijet)[ij]==abs(genIdx->at(igen))){
		      if(dCandGenMatchIdx->at(icand)==genIdx->at(igen) && abs(genPdg->at(igen))==pdgs[ij]){
			jetFF[ipt]->Fill(dCandPt->at(icand)/jtpt->at(ijet),weight);
			jetCorr->Fill(dCandPt->at(icand),jtpt->at(ijet),weight);
			jetDRMatched->Fill(djetR);
		      }
		    }
		  }
		}
	      }
	    }
	    else{
	      checkbgDCjets[type]=1;
	      //bgCjetSpec[type]->Fill(jtpt->at(ijet));
	      bgProximity[ipt][type]->Fill(djetR);
	      bgProximityEta[ipt][type]->Fill(deta);
	      bgProximityPhi[ipt][type]->Fill(dphi);
	    }
	  }
	  else{ //do all other d-mesons now
	    //if(abs(dCandGenMatchPdg->at(icand)) == pdgs[type]){
	    if(!fillDMassWithoutJets){
	      if(djetR<djetR_testCut){
		assns[ipt][0][type]->Fill(dCandMass->at(icand), weight);
		assnsUnweight[ipt][0][type]->Fill(dCandMass->at(icand));
		//}
		if(abs(dCandGenMatchDoubleCountedPdg->at(icand)) == pdgs[type]){
		  dcAssns[ipt][0][type]->Fill(dCandMass->at(icand));
		}
	      }
	    }
	    double lowball = snapToBin(assns[ipt][0][type],tightSignalCutLow[type],0);
	    double hiball = snapToBin(assns[ipt][0][type],tightSignalCutHi[type],1);
	    if(dCandMass->at(icand) > lowball && dCandMass->at(icand) < hiball){
	      proximity[ipt][type]->Fill(djetR);
	      proximityEta[ipt][type]->Fill(deta);
	      proximityPhi[ipt][type]->Fill(dphi);
	      checkDCjets[type]=1;
	      if(djetR<djetR_testCut){
		jetDR->Fill(djetR);
		for(unsigned int igen=0; igen<genPt->size(); igen++){
		  for(int ij=0; ij<5; ij++){
		    if(hasGenDwithKpiIdx->at(ijet)[ij]==abs(genIdx->at(igen))){
		      if(dCandGenMatchIdx->at(icand)==genIdx->at(igen) && abs(genPdg->at(igen))==pdgs[ij]){
			jetFF[ipt]->Fill(dCandPt->at(icand)/jtpt->at(ijet),weight);
			jetCorr->Fill(dCandPt->at(icand),jtpt->at(ijet),weight);
			jetDRMatched->Fill(djetR);
		      }
		    }
		  }
		}
	      }
	    }
	    else{
	      checkbgDCjets[type]=1;
	      bgProximity[ipt][type]->Fill(djetR);
	      bgProximityEta[ipt][type]->Fill(deta);
	      bgProximityPhi[ipt][type]->Fill(dphi);
	    }
	  }
	}
	else if(extraDCuts(tempCand,0)){
	  if(!fillDMassWithoutJets){
	    if(dCandType->at(icand)==1){
	      bgAssns[ipt][0][type]->Fill(dCandMass->at(icand) - dCandChildMass->at(icand),weight);
	      bgAssnsUnweight[ipt][0][type]->Fill((dCandMass->at(icand) - dCandChildMass->at(icand)));
	    }
	    else{
	      bgAssns[ipt][0][type]->Fill(dCandMass->at(icand),weight);
	      bgAssnsUnweight[ipt][0][type]->Fill((dCandMass->at(icand)));
	    }
	  }
	}
		// int pdgs[6] = {413, 421, 431, 431, 411, -999}; 
	// Dstar, D0, Ds, Ds, Dpm
	// type { 1,2,3,4,5}
	//do efficiency for D reco efficiency (step 2)
	for(unsigned int icandParton=0; icandParton<dCandParentPartonIdx->at(icand).size(); icandParton++){
	  if(genMatch->at(ijet) == dCandParentPartonIdx->at(icand).at(icandParton)){// && TMath::Abs(refparton_flavorForB->at(ijet))==4){
	    for(int ij=0; ij<5; ij++){
	      if(hasGenDwithKpi->at(ijet)[ij]){
		//cJetEffvsPt_bg_step2[ij]->Fill(jtpt->at(ijet)); //filled upstairs before dcand loop
		if(abs(dCandGenMatchPdg->at(icand)) == pdgs[ij]){
		  //step2cJetEffvsPt_step2[ij]->Fill(jtpt->at(ijet));
		}
		else{
		  //cJetEffvsPt_rej_step2[ij]->Fill(jtpt->at(ijet));
		}
	      }
	    }
	    break;
	  }
	}
	
	if(checkDCjets[type]){
	  deta = fabs(dCandEta->at(icand)-jteta->at(ijet));
	  dphi = fabs(dCandPhi->at(icand)-jtphi->at(ijet));
	  if(candDR[type]>sqrt(pow(deta,2)+pow(dphi,2))){
	    candDR[type] = sqrt(pow(deta,2)+pow(dphi,2));
	    assocPt[type] = dCandPt->at(icand);
	  }
	}
	if(checkbgDCjets[type]){
	  deta = fabs(dCandEta->at(icand)-jteta->at(ijet));
	  dphi = fabs(dCandPhi->at(icand)-jtphi->at(ijet));
	  if(bgcandDR[type]>sqrt(pow(deta,2)+pow(dphi,2))){
	    bgcandDR[type] = sqrt(pow(deta,2)+pow(dphi,2));
	    assocPt[type] = dCandPt->at(icand);
	  }
	}
      } //close candidate loop

      //fill jet spectra based on whats fg and bg
      for(int itype=0; itype<5; itype++){
	if(candDR[itype]<0.3){
	  CjetSpec[itype]->Fill(refpt->at(ijet),assocPt[itype],weight*feffline[itype]->Eval(assocPt[itype]));
	}
	else if(bgcandDR[itype]<0.3){
	  bgCjetSpec[itype]->Fill(refpt->at(ijet),assocPt[itype],weight*feffline[itype]->Eval(assocPt[itype]));
	}
      }

      //do efficiency for D branching ratio (step 3)
      if(TMath::Abs(refparton_flavorForB->at(ijet))==4){
	  //for(unsigned int icandParton=0; icandParton<dCandParentPartonIdx->at(icand).size(); icandParton++){
	    for(int ij=0; ij<5; ij++){
	      //if(genMatch->at(ijet) == dCandParentPartonIdx->at(icand).at(icandParton)){
	      if(hasGenD->at(ijet)[ij]){
		cJetEffvsPt_bg_step3[ij]->Fill(refpt->at(ijet));
		if(hasGenDwithKpi->at(ijet)[ij]){
		  cJetEffvsPt_step3[ij]->Fill(refpt->at(ijet));
		}
		else{
		  cJetEffvsPt_rej_step3[ij]->Fill(refpt->at(ijet));
		}
	      }
	    }
	    //  }
      }
      
      //do efficiency for c->D fragmentation fraction
	 if(TMath::Abs(refparton_flavorForB->at(ijet))==4){
	  for(int ij=0; ij<5; ij++){
	    cJetEffvsPt_bg_step4[ij]->Fill(refpt->at(ijet));
	    // for(unsigned int icandParton=0; icandParton<dCandParentPartonIdx->at(icand).size(); icandParton++){
	      //if(genMatch->at(ijet) == dCandParentPartonIdx->at(icand).at(icandParton)){
		if(hasGenD->at(ijet)[ij]){
		  cJetEffvsPt_step4[ij]->Fill(refpt->at(ijet));
		}
		//   }
	      else cJetEffvsPt_rej_step4[ij]->Fill(refpt->at(ijet));
		// }
	  }
	}
    
      if(!anotherTestFlag && TMath::Abs(refparton_flavorForB->at(ijet))==4){
	unmatchedCJets->Fill(refpt->at(ijet));
      }
      
      if(candJetMatchIdx>0){
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
  for(int i=0; i<ptBins; i++){
    for(int j=0; j<djetBins; j++){
      for(int k=0; k<cTypes; k++){
      assns[i][j][k]->Write();
      assnsUnweight[i][j][k]->Write();
      bgAssnsUnweight[i][j][k]->Write();
      bgAssns[i][j][k]->Write();
      dcAssns[i][j][k]->Write();
      nondcAssns[i][j][k]->Write();
      }
    }
  }

  for(int i=0; i<ptBins; i++){
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
  for(int i=0; i<5; i++){
    cJetEffvsPt_step1[i]->Write();
    cJetEffvsPt_bg_step1[i]->Write();
    cJetEffvsPt_rej_step1[i]->Write();
    cJetEffvsPt_step2[i]->Write();
    cJetEffvsPt_bg_step2[i]->Write();
    cJetEffvsPt_rej_step2[i]->Write();
    cJetEffvsPt_step3[i]->Write();
    cJetEffvsPt_bg_step3[i]->Write();
    cJetEffvsPt_rej_step3[i]->Write();
    cJetEffvsPt_step4[i]->Write();
    cJetEffvsPt_bg_step4[i]->Write();
    cJetEffvsPt_rej_step4[i]->Write();
  }

  jetAssnTop->Write();
  jetAssnBot->Write();

  dRTestHist->Write();
  unmatchedCJets->Write();

  for(int i=0; i<5; i++){
    cMesonRecoEff[i]->Write();
    cMesonRecoEff_bg[i]->Write();
  }

  for(int i=0; i<cTypes; i++){
    purity[i]->Divide(purity[i],purity_bg[i],1,1,"B");
    purity[i]->Write();
  }

  Color_t colors[ptBins] = {kBlack,kRed,kOrange+1,kGreen+2,kBlue+2,kMagenta+2};
  for(int i=0; i<ptBins; i++){
    jetFF[i]->Scale(1./jetFF[i]->Integral());
    jetFF[i]->SetMarkerColor(colors[i]);
    jetFF[i]->Write();
  }
  
  for(int i=0; i<5; i++){
    cMesonRecoEfftest[i]->Divide(cMesonRecoEfftest_bg[i]);
    cMesonRecoEfftest[i]->Write();
    cMesonRecoEfftest_bg[i]->Write();

    feffline[i]->Write();
  }

  
  
  jetCorr->SetMinimum(1E-8);
  jetCorr->Write();

  jetDR->Write();
  jetDRMatched->Write();
  fout->Close();

  //fpur->cd();

  //fpur->Close();
}
