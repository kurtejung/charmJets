#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>
#include <set>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TF1.h"
#include "TMath.h"
#include "TVector3.h"

#ifdef __CINT__
#pragma link C++ class vector<vector<int> >+;
#endif

using namespace std;

const double PI=3.1415926;

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
const double lowMassBin[cTypes] = {0.14, 1.68, 1.8, 1.8, 1.68};
const double highMassBin[cTypes] = {0.155, 1.95, 2.1, 2.1, 2.05};
const double daughterLowMassBin[cTypes] = {1.8, -1, 1.00, 1.00, 0.5};
const double daughterHighMassBin[cTypes] = {1.92, 1, 1.04, 1.04, 2.05};
const double sigmaM[cTypes] = {6.1e-4, 1.32e-2, 1.55e-2, 1.67e-2, 1.31e-2};
const double MesonMass[cTypes] = {0.1455,1.8648,1.9685,1.9685,1.8696};

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

double rotateMass(double mD, double kpt, double keta, double kphi, double pipt, double pieta, double piphi, int method = 2){

  //m^2 = E^2-p^2
  //p = pT*cosh(eta)

  //Method 1:
  //Rotation of Kaon by R_z(180) [rotation about z axis by 180 deg] =
  // -> M_D - M_D(rot) = 2px(pi)*px(k) + 2py(pi)*py(k)
  // so M_D(rot) = M_D - 2pt(pi)*pt(k)[ cos(phi)_pi*cos(phi)_k + sin(phi)_pi*sin(phi)_k ]
		
  double diff = 4*pipt*kpt;
  //cout << "piphi: " << piphi << " kphi: "<< kphi << endl;
  double angles = cos(piphi-kphi); //cos(piphi)*cos(kphi) + sin(piphi)*sin(kphi); [[cos(A-B) = cos(A)cos(B)+sin(A)sin(B)]]
  diff *=angles;
  double rotMass = sqrt(abs(pow(mD,2) - diff));

  //Method 2 (direct recalculation using pT and phi):
  double newkphi = kphi+PI;
  TVector3 kmomentum(kpt*cos(newkphi),kpt*sin(newkphi),kpt*cosh(keta));
  TVector3 pimomentum(pipt*cos(piphi),pipt*sin(piphi),pipt*cosh(pieta));
  double mkaon = 0.493667; //gev
  double mpion = 0.13957018; //gev

  double kenergy2 = kmomentum.Mag2() + pow(mkaon,2);
  double pienergy2 = pimomentum.Mag2() + pow(mpion,2);
  double newmass2 = pow(mkaon,2) + pow(mpion,2) + 2*sqrt(kenergy2)*sqrt(pienergy2) - 2*kmomentum.Dot(pimomentum);

  //cout << "method 1: " << rotMass << " method 2: "<< sqrt(newmass2) << endl;
  if(rotMass<2 || sqrt(newmass2)<2) cout << "real mass found!! " << rotMass << endl; 
  
  if(method==1) return rotMass;
  else return sqrt(newmass2);
}

void makeCharmJetAssn_CJetSample(int isMC=6, bool doDraw = 0){

  const int ptBins = 9;
  const int jetPtBins = 6;
  const int djetBins = 1;

  const double sigbg = 4.0; //change the distance away from the mass peak to start bg estimation in # of sigmas
  //this needs to be at least 2 sigma!!

  const bool fillDMassWithoutJets = false;
  const bool findRecoEffVsDMesonPt = false; //can't do this anymore - d reco eff goes to 0 so some jets are just lost

  const double djetR_testCut = 0.3;
  double drCut = djetR_testCut;

  int ptBoundaries[jetPtBins+1] = {20,40,60,80,100,200,400};
  int dCandptBoundaries[ptBins+1] = {0,4,6,8,12,15,20,40,60,100};
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
  else if(isMC==6) fout = new TFile("dCandidateToJetAssns_CJetOnly_RevCuts_JetAssn-12-12.root","recreate");
  else if(isMC==7) fout = new TFile("dCandidateToJetAssns_officialQCDJetOnly_RevCuts_JetAssn-11-11.root","recreate");
  else{
    cout << "isMC cannot be more than 6! Is this a new D-decay type??" << endl;
    exit(0);
  }

  //mass-window selections
  const double tightSignalCutLow[cTypes] = {0.144, 1.84, 1.93, 1.93, 1.85};
  const double tightSignalCutHi[cTypes] = {0.147, 1.885, 1.99, 1.99, 1.885};
  
  const int jetBins = 10;
  double xjetBins[jetBins+1] = {0,5,10,15,20,40,60,80,100,200,400};

  const int dmesonBins = 19;
  double xmesonBins[dmesonBins+1] = {0,4,6,8,10,12,14,16,18,20,25,30,40,60,80,100,125,150,175,200};

  const int ndrbins = 18;
  double xdrbins[ndrbins+1] = {0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.3,0.4,0.6,1.,1.5,2,2.5,3,3.5,4};

  TH1D *allJetSpec = new TH1D("allJetSpec","",jetBins,xjetBins); allJetSpec->Sumw2();
  TH1D *CjetSpec[cTypes];
  TH1D *bgCjetSpec[cTypes];
  TH1D *purity_bg[cTypes];
  TH1D *purity[cTypes];
  TH1D *DMesonSpec[cTypes];

  TH1D *proximity[ptBins][cTypes];
  TH1D *Jetproximity[ptBins][cTypes];
  TH1D *proximityTruth[ptBins][cTypes];
  TH1D *proximityComb[ptBins][cTypes];
  TH1D *proximityEta[ptBins][cTypes];
  TH1D *proximityPhi[ptBins][cTypes];
  TH1D *bgProximity[ptBins][cTypes];
  TH1D *bgRotate[ptBins][djetBins][cTypes];
  TH1D *bgProximityTruth[ptBins][cTypes];
  TH1D *bgJetproximity[ptBins][cTypes];
  TH1D *bgProximityEta[ptBins][cTypes];
  TH1D *bgProximityPhi[ptBins][cTypes];
  TH1D *assns[ptBins][djetBins][cTypes];
  TH1D *truthAssns[ptBins][djetBins][cTypes];
  TH1D *combAssns[ptBins][djetBins][cTypes];
  TH1D *assnsUnweight[ptBins][djetBins][cTypes];
  TH1D *bgAssns[ptBins][djetBins][cTypes];
  TH1D *bgAssnsUnweight[ptBins][djetBins][cTypes];
  TH1D *dcAssns[ptBins][djetBins][cTypes];
  TH1D *nondcAssns[ptBins][djetBins][cTypes];
  TH1D *bgProximityLeftSide[ptBins][cTypes];
  TH1D *bgProximityRightSide[ptBins][cTypes];
  TH1D *bgProximityLeftSideTruth[ptBins][cTypes];
  TH1D *bgProximityRightSideTruth[ptBins][cTypes];

  TH1D *mesonMults[3]; //find how many mesons are inside the jet cone for each jet type (light,c,b)
  TH1D *mesonMultsNoDstar[3];
  for(int i=0; i<3; i++){
    mesonMults[i] = new TH1D(Form("mesonMults_%d",i),"",10,0,10);
    mesonMultsNoDstar[i] = new TH1D(Form("mesonMultsNoDstar_%d",i),"",10,0,10);
  }

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
    if(findRecoEffVsDMesonPt) cJetEffvsPt_step2[i] = new TH1D(Form("cJetEffvsPt_step2_type%d",i),"",dmesonBins,xmesonBins);
    else cJetEffvsPt_step2[i] = new TH1D(Form("cJetEffvsPt_step2_type%d",i),"",jetBins,xjetBins);
    cJetEffvsPt_step2[i]->Sumw2();
    cJetEffvsPt_step2[i]->SetDirectory(0);
    cJetEffvsPt_step3[i] = new TH1D(Form("cJetEffvsPt_step3_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_step3[i]->Sumw2();
    cJetEffvsPt_step4[i] = new TH1D(Form("cJetEffvsPt_step4_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_step4[i]->Sumw2();

    cJetEffvsPt_bg_step1[i] = new TH1D(Form("cJetEffvsPt_bg_step1_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_bg_step1[i]->Sumw2();
    if(findRecoEffVsDMesonPt) cJetEffvsPt_bg_step2[i] = new TH1D(Form("cJetEffvsPt_bg_step2_type%d",i),"",dmesonBins,xmesonBins);
    else cJetEffvsPt_bg_step2[i] = new TH1D(Form("cJetEffvsPt_bg_step2_type%d",i),"",jetBins,xjetBins);
    cJetEffvsPt_bg_step2[i]->Sumw2();
    cJetEffvsPt_bg_step2[i]->SetDirectory(0);
    cJetEffvsPt_bg_step3[i] = new TH1D(Form("cJetEffvsPt_bg_step3_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_bg_step3[i]->Sumw2();
    cJetEffvsPt_bg_step4[i] = new TH1D(Form("cJetEffvsPt_bg_step4_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_bg_step4[i]->Sumw2();

    cJetEffvsPt_rej_step1[i] = new TH1D(Form("cJetEffvsPt_rej_step1_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_rej_step1[i]->Sumw2();
    if(findRecoEffVsDMesonPt) cJetEffvsPt_rej_step2[i] = new TH1D(Form("cJetEffvsPt_rej_step2_type%d",i),"",dmesonBins,xmesonBins);
    else cJetEffvsPt_rej_step2[i] = new TH1D(Form("cJetEffvsPt_rej_step2_type%d",i),"",jetBins,xjetBins);
    cJetEffvsPt_rej_step2[i]->Sumw2();
    cJetEffvsPt_rej_step2[i]->SetDirectory(0);
    cJetEffvsPt_rej_step3[i] = new TH1D(Form("cJetEffvsPt_rej_step3_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_rej_step3[i]->Sumw2();
    cJetEffvsPt_rej_step4[i] = new TH1D(Form("cJetEffvsPt_rej_step4_type%d",i),"",jetBins,xjetBins); cJetEffvsPt_rej_step4[i]->Sumw2();
  }

  for(int i=0; i<ptBins; i++){
    for(int j=0; j<djetBins; j++){
      for(int k=0; k<cTypes; k++){
	assns[i][j][k] = new TH1D(Form("assns_pt%d_djet%d_type%d",i,j,k),"",50,lowMassBin[k],highMassBin[k]);
	assns[i][j][k]->Sumw2();
	assns[i][j][k]->SetXTitle("Candidate Mass [GeV/c^{2}]");
	assns[i][j][k]->SetYTitle("Weighted Counts");

	combAssns[i][j][k] = new TH1D(Form("combAssns_pt%d_djet%d_type%d",i,j,k),"",50,lowMassBin[k],highMassBin[k]);
	combAssns[i][j][k]->Sumw2();
	combAssns[i][j][k]->SetXTitle("Candidate Mass [GeV/c^{2}]");
	combAssns[i][j][k]->SetYTitle("Weighted Counts");

	truthAssns[i][j][k] = new TH1D(Form("truthAssns_pt%d_djet%d_type%d",i,j,k),"",50,lowMassBin[k],highMassBin[k]);
	truthAssns[i][j][k]->Sumw2();
	truthAssns[i][j][k]->SetXTitle("Candidate Mass [GeV/c^{2}]");
	truthAssns[i][j][k]->SetYTitle("Weighted Counts");
	
	assnsUnweight[i][j][k] = new TH1D(Form("assns_unweight_pt%d_djet%d_type%d",i,j,k),"",50,lowMassBin[k],highMassBin[k]);
	assnsUnweight[i][j][k]->Sumw2();
	assnsUnweight[i][j][k]->SetXTitle("Candidate Mass [GeV/c^{2}]");
	assnsUnweight[i][j][k]->SetYTitle("Unweighted Counts");

	bgAssns[i][j][k] = new TH1D(Form("bgAssns_pt%d_djet%d_type%d",i,j,k),"",50,lowMassBin[k],highMassBin[k]);
	//bgAssns[i][j][k]->Sumw2();
	bgAssns[i][j][k]->SetXTitle("Candidate Mass [GeV/c^{2}]");
	bgAssns[i][j][k]->SetYTitle("Weighted Counts");

	bgAssnsUnweight[i][j][k] = new TH1D(Form("bgAssns_unweight_pt%d_djet%d_type%d",i,j,k),"",50,lowMassBin[k],highMassBin[k]);
	//bgAssnsUnweight[i][j][k]->Sumw2();
	bgAssnsUnweight[i][j][k]->SetXTitle("Candidate Mass [GeV/c^{2}]");
	bgAssnsUnweight[i][j][k]->SetYTitle("Unweighted Counts");

	dcAssns[i][j][k] = new TH1D(Form("DC_assns_pt%d_djet%d_type%d",i,j,k),"",50,lowMassBin[k],highMassBin[k]);
	dcAssns[i][j][k]->Sumw2();
	dcAssns[i][j][k]->SetXTitle("Candidate Mass [GeV/c^{2}]");
	dcAssns[i][j][k]->SetYTitle("Unweighted Counts");

	nondcAssns[i][j][k] = new TH1D(Form("non_DC_assns_pt%d_djet%d_type%d",i,j,k),"",50,lowMassBin[k],highMassBin[k]);
	nondcAssns[i][j][k]->Sumw2();
	nondcAssns[i][j][k]->SetXTitle("Candidate Mass [GeV/c^{2}]");
	nondcAssns[i][j][k]->SetYTitle("Unweighted Counts");

	//bgRotate[i][j][k] = new TH1D(Form("bgRotate_pt%d_djet%d_type%d",i,j,k),"",50,lowMassBin[k],highMassBin[k]);
	bgRotate[i][j][k] = new TH1D(Form("bgRotate_pt%d_djet%d_type%d",i,j,k),"",50,0,100);
	bgRotate[i][j][k]->Sumw2();
	bgRotate[i][j][k]->SetXTitle("Candidate Mass [GeV/c^{2}]");
	bgRotate[i][j][k]->SetYTitle("Unweighted Counts");
      }
    }
  }
  
 
  for(int i=0; i<ptBins; i++){
    for(int k=0; k<cTypes; k++){
      proximity[i][k] = new TH1D(Form("proximity_pt%d_type%d",i,k),"",ndrbins,xdrbins);
      Jetproximity[i][k] = new TH1D(Form("Jetproximity_pt%d_type%d",i,k),"",ndrbins,xdrbins);
      proximityTruth[i][k] = new TH1D(Form("proximityTruth_pt%d_type%d",i,k),"",ndrbins,xdrbins);
      proximityComb[i][k] = new TH1D(Form("proximityComb_pt%d_type%d",i,k),"",ndrbins,xdrbins);
      proximityEta[i][k] = new TH1D(Form("proximityEta_pt%d_type%d",i,k),"",ndrbins,xdrbins);
      proximityPhi[i][k] = new TH1D(Form("proximityPhi_pt%d_type%d",i,k),"",ndrbins,xdrbins);
      bgProximity[i][k] = new TH1D(Form("bgProximity_pt%d_type%d",i,k),"",ndrbins,xdrbins);
      bgProximityTruth[i][k] = new TH1D(Form("bgProximityTruth_pt%d_type%d",i,k),"",ndrbins,xdrbins); bgProximityTruth[i][k]->Sumw2();
      bgJetproximity[i][k] = new TH1D(Form("bgJetproximity_pt%d_type%d",i,k),"",ndrbins,xdrbins);
      bgProximityEta[i][k] = new TH1D(Form("bgProximityEta_pt%d_type%d",i,k),"",ndrbins,xdrbins);
      bgProximityPhi[i][k] = new TH1D(Form("bgProximityPhi_pt%d_type%d",i,k),"",ndrbins,xdrbins);
      bgProximityLeftSide[i][k] = new TH1D(Form("bgProximityLeftSide_pt%d_type%d",i,k),"",ndrbins,xdrbins);
      bgProximityRightSide[i][k] = new TH1D(Form("bgProximityRightSide_pt%d_type%d",i,k),"",ndrbins,xdrbins);
      bgProximityLeftSideTruth[i][k] = new TH1D(Form("bgProximityLeftSideTruth_pt%d_type%d",i,k),"",ndrbins,xdrbins);
      bgProximityRightSideTruth[i][k] = new TH1D(Form("bgProximityRightSideTruth_pt%d_type%d",i,k),"",ndrbins,xdrbins);
    }
  }
  
  for(int k=0; k<cTypes; k++){
    CjetSpec[k] = new TH1D(Form("CjetSpec_type%d",k),"",jetBins,xjetBins); CjetSpec[k]->Sumw2();
    bgCjetSpec[k] = new TH1D(Form("bgCjetSpec_type%d",k),"",jetBins,xjetBins); bgCjetSpec[k]->Sumw2();
    purity[k] = new TH1D(Form("purity_type%d",k),"",jetBins,xjetBins); purity[k]->Sumw2();
    purity_bg[k] = new TH1D(Form("purity_bg_type%d",k),"",jetBins,xjetBins); purity_bg[k]->Sumw2();
    DMesonSpec[k] = new TH1D(Form("DMesonSpec_%d",k),"",20,0,100); DMesonSpec[k]->Sumw2();
  }

  double weight;
  int HLT_Jet20_NoJetID_v1, HLT_Jet40_NoJetID_v1, HLT_Jet60_NoJetID_v1, HLT_Jet80_NoJetID_v1, HLT_Jet100_NoJetID_v1;
  vector<double> *dCandPt=0, *dCandMass=0, *dCandEta=0, *dCandPhi=0, *jtpt=0, *jteta=0, *jtphi=0, *refpt=0, *rawpt=0;
  vector<double> *dCandType=0, *dCandChildMass=0, *dCandCharge1=0, *dCandCharge2=0;
  vector<double> *dCandKPt=0, *dCandKEta=0, *dCandKPhi=0, *dCandPiPt=0, *dCandPiEta=0, *dCandPiPhi=0;
  vector<int> *dCandMatchGenPdg_dau2=0, *dCandMatchGenPdg_index1=0, *dCandMatchGenPdg_index2=0;
  vector<int> *refparton_flavorForB=0;
  vector<vector<int> > *dCandParentPartonIdx=0;
  vector<int> *genMatch=0, *dCandGenMatchPdg=0, *dCandGenMatchDoubleCountedPdg=0, *dCandGenMatchIdx=0, *dCandGenMatchDoubleCountedIdx=0;
  TBranch *bdcandpt=0, *bdcandmass=0, *bdcandeta=0, *bdcandphi=0, *bdcandtype=0, *bdcandchildmass, *bjtpt=0, *bjteta=0, *bjtphi=0, *bdcandcharge1=0, *bdcandcharge2=0;
  vector<vector<int> > *hasGenD=0, *hasGenDIdx=0, *hasGenDwithKpi=0, *hasGenDwithKpiIdx=0;

  //d*, d0, ds->phipi, ds->k*k, dpm
  TFile *fin;
  if(!isMC) fin = new TFile("input/DMesonCJet_RevCuts_NoJetTrgCut_pPbdata_ppReco_akPu3PF.root","OLD");
  if(isMC==1) fin = new TFile("DMesonCJet_Dstar_pPbMC_ppReco_akPu3PF_1.root","OLD");
  if(isMC==2) fin = new TFile("DMesonCJet_Dzero_pPbMC_ppReco_akPu3PF_1.root","OLD");
  if(isMC==3) fin = new TFile("DMesonCJet_Ds2PhiPi_pPbMC_ppReco_akPu3PF_76.root","OLD");
  if(isMC==4) fin = new TFile("DMesonCJet_Ds2KstarKEmbed_pPbMC_ppReco_akPu3PF_1.root","OLD");
  if(isMC==5) fin = new TFile("DMesonCJet_DpmEmbed_pPbMC_ppReco_akPu3PF_1.root","OLD");
  if(isMC==6) fin = new TFile("input/DMesonCJet_CJetOnly_pPbMC_ppReco_akPu3PF_1202.root","OLD");
  if(isMC==7) fin = new TFile("input/DMesonCJet_officialQCDJetOnly_pPbMC_ppReco_akPu3PF_1341.root","OLD");
  if(!fin){
    cout << "input file not found! " <<endl;
    exit(0);
  }

  TTree *ct = (TTree*)fin->Get("ct");

  ct->SetBranchAddress("dCandPt",&dCandPt, &bdcandpt);
  ct->SetBranchAddress("dCandEta",&dCandEta, &bdcandeta);
  ct->SetBranchAddress("dCandPhi",&dCandPhi, &bdcandphi);
  ct->SetBranchAddress("dCandMass",&dCandMass, &bdcandmass);
  ct->SetBranchAddress("dCandType",&dCandType, &bdcandtype);
  ct->SetBranchAddress("dCandChildMass",&dCandChildMass, &bdcandchildmass);
  ct->SetBranchAddress("dCandCharge1",&dCandCharge1, &bdcandcharge1);
  ct->SetBranchAddress("dCandCharge2",&dCandCharge2, &bdcandcharge2);
  ct->SetBranchAddress("dCandKPt",&dCandKPt);
  ct->SetBranchAddress("dCandKEta",&dCandKEta);
  ct->SetBranchAddress("dCandKPhi",&dCandKPhi);
  ct->SetBranchAddress("dCandPiPt",&dCandPiPt);
  ct->SetBranchAddress("dCandPiEta",&dCandPiEta);
  ct->SetBranchAddress("dCandPiPhi",&dCandPiPhi);
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
    ct->SetBranchAddress("hasGenDIdx",&hasGenDIdx);
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
  vector<int> *genIdx=0, *genPdg=0;
  vector<int> *hasHadronicDecay=0;
  if(isMC){
    gt = (TTree*)fin->Get("gTree");
    gt->SetBranchAddress("genPt",&genPt);
    gt->SetBranchAddress("genEta",&genEta);
    gt->SetBranchAddress("genPhi",&genPhi);
    gt->SetBranchAddress("genIdx",&genIdx); //event-level index of particle
    gt->SetBranchAddress("genPdg",&genPdg);
    //gt->SetBranchAddress("genPtIdx2",&genPtIdx2);
    gt->SetBranchAddress("hasHadronicDecay",&hasHadronicDecay);
  }

  TF1 *feffline[5];
  for(int ij=0; ij<5; ij++){
    feffline[ij] = new TF1(Form("feffline_%d",ij),"[0]*exp(-[1]/x)+[2]",0,400);
    feffline[ij]->SetParameter(0,0.5);
    feffline[ij]->SetParameter(1,10);
    feffline[ij]->SetParameter(2,5);
  }

  int nentries = ct->GetEntries();
  if(findRecoEffVsDMesonPt){
    for(int ientry=0; ientry<nentries; ientry++){
      ct->GetEntry(ientry);
      if(isMC) gt->GetEntry(ientry);

      if(!HLT_Jet20_NoJetID_v1 && !HLT_Jet40_NoJetID_v1 && !HLT_Jet60_NoJetID_v1 && !HLT_Jet80_NoJetID_v1 && !HLT_Jet100_NoJetID_v1){
      //cout << "no triggers fired! skipping event..." << endl;
      continue;
    }
      if(jtpt->size()==0) continue;
      bool cFlag=false;
      //for all gen-level d-mesons, how many are reconstructed?
      //do reconstruction efficiency starting from the gen-level d-meson jet association
      for(unsigned int ijet=0; ijet<jtpt->size(); ijet++){
	if(refpt->at(ijet)<20) continue;
	if(rawpt->at(ijet)<23) continue;
	if(TMath::Abs(refparton_flavorForB->at(ijet))!=4) continue;
	//if(cFlag) continue;
	cFlag=true; //check to make sure multiple jets don't loop through the same mesons twice
     	  for(int ij=0; ij<5; ij++){
	    for(unsigned int igen=0; igen<genPt->size(); igen++){
	      double genPtToFill=0;
	      bool FF2=false;
	      bool foundFlag=false;
	      if(hasGenDwithKpiIdx->at(ijet)[ij]==genIdx->at(igen) && abs(genPdg->at(igen))==pdgs[ij] && hasHadronicDecay->at(igen)){
		FF2 = true;
		genPtToFill = genPt->at(igen);
		for(unsigned int icand=0; icand<dCandPt->size(); icand++){
		  if(dCandGenMatchIdx->at(icand)==genIdx->at(igen) && dCandGenMatchIdx->at(icand) != -999){
		    foundFlag=true;
		  }
		}
	      }
	    
	    
	      if(FF2) cJetEffvsPt_bg_step2[ij]->Fill(genPtToFill);
	      if(foundFlag && FF2){
		cJetEffvsPt_step2[ij]->Fill(genPtToFill);
	      }
	      else{
		cJetEffvsPt_rej_step2[ij]->Fill(genPtToFill);
	      }
	      if(ij==1 && FF2) cout << "event: "<< ientry << " jet pt: "<< refpt->at(ijet) << " meson pt: "<< genPtToFill << " bg filling, type: " << ij << endl; 
	      if(ij==1 && foundFlag && FF2) cout << "event: "<< ientry << " pt: "<< genPtToFill << " type: "<< ij << endl;
	    }
	  }    
      }
    }
  

    for(int ij=0; ij<5; ij++){
      cJetEffvsPt_step2[ij]->Divide(cJetEffvsPt_step2[ij],cJetEffvsPt_bg_step2[ij],1,1,"B");
      cJetEffvsPt_step2[ij]->Fit(feffline[ij],"qN","",5,100);
    }
  }

  //nentries = 500;
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
    
    if(!HLT_Jet20_NoJetID_v1 && !HLT_Jet40_NoJetID_v1 && !HLT_Jet60_NoJetID_v1 && !HLT_Jet80_NoJetID_v1 && !HLT_Jet100_NoJetID_v1){
      //cout << "no triggers fired! skipping event..." << endl;
      continue;
    }
    
    if(dCandPt->size() != dCandParentPartonIdx->size()) cout << "warning! candidate parent parton != dcandPt!" << endl;
    //if(dCandPt->size() != hasGenD->size()) cout << "warning! d-match flag != dCandPt!" << endl;
    //cout << "dcandpt size: "<< dCandPt->size() << " hasgend size: "<< hasGenD->size() << endl;

    for(unsigned int ijet=0; ijet<jtpt->size(); ijet++){
      if(refpt->at(ijet)<20) continue;
      if(rawpt->at(ijet)<23) continue;
      int jtype=-1;
      int dcounter=0;
      int dcounterNoDstar=0;
      int jettype = abs(refparton_flavorForB->at(ijet));
      if(jettype==5) jtype=2;
      else if(jettype==4) jtype=1;
      else jtype=0;
      for(unsigned int igen=0; igen<genPt->size(); igen++){
	double deta = fabs(genEta->at(igen)-jteta->at(ijet));
	double dphi = fabs(genPhi->at(igen)-jtphi->at(ijet));
	double djetR = sqrt(pow(deta,2)+pow(dphi,2));
	if(djetR<0.3){
	  dcounter++;
	  if(abs(genPdg->at(igen))==421){
	    dcounterNoDstar++;
	  }
	}
      }
      mesonMults[jtype]->Fill(dcounter);
      mesonMultsNoDstar[jtype]->Fill(dcounterNoDstar);
    }
    

    for(unsigned int igen=0; igen<genPt->size(); igen++){
      for(unsigned int ijet=0; ijet<jtpt->size(); ijet++){
	if(refpt->at(ijet)<20) continue;
	if(rawpt->at(ijet)<23) continue;
	//if(genMatch->at(ijet)==genIdx->at(igen)){
	double deta = fabs(genEta->at(igen)-jteta->at(ijet));
	double dphi = fabs(genPhi->at(igen)-jtphi->at(ijet));
	double djetR = sqrt(pow(deta,2)+pow(dphi,2));
	if(djetR<drCut){
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
      if(TMath::Abs(refparton_flavorForB->at(ijet))!=4) continue;
      for(int ij=0; ij<5; ij++){
	bool FF1=false;
	bool FF2=false;
	for(unsigned int icand=0; icand<dCandPt->size(); icand++){
	  for(unsigned int igen=0; igen<genPt->size(); igen++){
	    if(hasGenDwithKpiIdx->at(ijet)[ij]==genIdx->at(igen) && abs(dCandGenMatchIdx->at(icand))==genIdx->at(igen)){
	      FF1 = true;
	      double djetR = sqrt(pow((dCandEta->at(icand)-jteta->at(ijet)),2)+pow((dCandPhi->at(icand)-jtphi->at(ijet)),2));
	      if(djetR<djetR_testCut){
		FF2=true;
	      }
	    }
	  }
	}
	if(FF1){
	  cJetEffvsPt_bg_step1[ij]->Fill(refpt->at(ijet),1.);
	  if(FF2) cJetEffvsPt_step1[ij]->Fill(refpt->at(ijet),1.);
	  else cJetEffvsPt_rej_step1[ij]->Fill(refpt->at(ijet),1.);
	}
      }
    }


    //for all gen-level d-mesons, how many are reconstructed?
    //do reconstruction efficiency starting from the gen-level d-meson jet association
    if(!findRecoEffVsDMesonPt){
      for(unsigned int ijet=0; ijet<jtpt->size(); ijet++){
	if(refpt->at(ijet)<20) continue;
	if(rawpt->at(ijet)<23) continue;
	if(TMath::Abs(refparton_flavorForB->at(ijet))==4){
	  for(int ij=0; ij<5; ij++){
	    bool FF2=false;
	    bool foundFlag=false;
	    for(unsigned int igen=0; igen<genPt->size(); igen++){
	      if(hasGenDwithKpiIdx->at(ijet)[ij]==genIdx->at(igen) && abs(genPdg->at(igen))==pdgs[ij]){
		FF2 = true;
		for(unsigned int icand=0; icand<dCandPt->size(); icand++){
		  //wtf ok so using double-counted idx works for D0 but regular works for D+/-...
		  if(dCandGenMatchIdx->at(icand)==genIdx->at(igen) && dCandGenMatchIdx->at(icand) != -999 && dCandType->at(icand)-1==ij){
		    foundFlag=true;
		  }
		}
	      }
	    }
	    if(FF2){
	      cJetEffvsPt_bg_step2[ij]->Fill(refpt->at(ijet),1.);
	      if(foundFlag){
		cJetEffvsPt_step2[ij]->Fill(refpt->at(ijet),1.);
	      }
	      else{
		cJetEffvsPt_rej_step2[ij]->Fill(refpt->at(ijet),1.);
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
	if(FF2) cMesonRecoEfftest_bg[ij]->Fill(genPt->at(igen));
	if(foundFlag) cMesonRecoEfftest[ij]->Fill(genPt->at(igen));
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
	  double rotation = rotateMass(dCandMass->at(icand), dCandKPt->at(icand), dCandKEta->at(icand), dCandKPhi->at(icand), dCandPiPt->at(icand), dCandPiEta->at(icand), dCandPiPhi->at(icand));
	  if(dCandType->at(icand)==1){
	    bgRotate[iCandpt][0][type]->Fill(rotation - dCandChildMass->at(icand));
	    bgAssnsUnweight[iCandpt][0][type]->Fill((dCandMass->at(icand) - dCandChildMass->at(icand)));
	    bgAssns[iCandpt][0][type]->Fill(dCandMass->at(icand) - dCandChildMass->at(icand),weight);
	  }
	  else{
	    bgRotate[iCandpt][0][type]->Fill(rotation);
	    bgAssnsUnweight[iCandpt][0][type]->Fill((dCandMass->at(icand)));
	    bgAssns[iCandpt][0][type]->Fill(dCandMass->at(icand),weight);
	  }
	}
      }
    }
    
    for(unsigned int ijet=0; ijet<jtpt->size(); ijet++){
      if(refpt->at(ijet)<20) continue;
      if(rawpt->at(ijet)<23) continue;
      if(TMath::Abs(refparton_flavorForB->at(ijet))!=4) continue;
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
      int closestDRidx[5] = {-1,-1,-1,-1,-1};
      double closestDR[5] = {999,999,999,999,999};
      
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
      int ipt=0;
      while(refpt->at(ijet) > ptBoundaries[ipt+1] && ipt<(jetPtBins-1)) ipt++;
      
      //if(maxcandno>0){
      for(unsigned int icand=0; icand<dCandPt->size(); icand++){
	//if(dCandPt->at(icand) < 10) continue;
	//int icand = maxcandno; //------------------------flip this switch on to just look at the highest pt d-meson
	//if(dCandPt->at(icand) > maxcandpt) continue;
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
	while(dCandPt->at(icand) > dCandptBoundaries[iCandpt+1] && iCandpt<(ptBins-1)) iCandpt++;
	if(iCandpt > jetPtBins-1) iCandpt=jetPtBins-1;
	int type = dCandType->at(icand)-1;

	for(int ij=0; ij<5; ij++){
	  // if(abs(dCandGenMatchPdg->at(icand))==pdgs[ij] && dCandGenMatchIdx->at(icand)==hasGenDwithKpiIdx->at(ijet)[ij] && dCandGenMatchIdx->at(icand) != -999 && dCandType->at(icand)-1==ij) proximityTruth[ipt][ij]->Fill(djetR);
	}
	
	double rotation = rotateMass(dCandMass->at(icand), dCandKPt->at(icand), dCandKEta->at(icand), dCandKPhi->at(icand), dCandPiPt->at(icand), dCandPiEta->at(icand), dCandPiPhi->at(icand));
	if(type==0) bgRotate[ipt][0][type]->Fill(rotation - dCandChildMass->at(icand));
	else bgRotate[ipt][0][type]->Fill(rotation);
	
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
	      
	      if(abs(dCandGenMatchPdg->at(icand))==pdgs[type] && dCandGenMatchIdx->at(icand)==hasGenDwithKpiIdx->at(ijet)[type]){
		truthAssns[ipt][0][type]->Fill(dCandMass->at(icand)-dCandChildMass->at(icand));
	      }
	      if(abs(dCandGenMatchPdg->at(icand))==pdgs[type] && dCandGenMatchIdx->at(icand)!=hasGenDwithKpiIdx->at(ijet)[type]){
		combAssns[ipt][0][type]->Fill(dCandMass->at(icand)-dCandChildMass->at(icand));
	      }
	      }
	    }
	    if(djetR<closestDR[type]){ closestDR[type]=djetR; closestDRidx[type]=icand; }
	    double lowball = snapToBin(assns[ipt][0][type],tightSignalCutLow[type],0);
	    double hiball = snapToBin(assns[ipt][0][type],tightSignalCutHi[type],1);
	    double lowballbg = snapToBin(assns[ipt][0][type], MesonMass[type]-(sigmaM[type]*sigbg),0);
	    double hiballbg = snapToBin(assns[ipt][0][type], MesonMass[type]+(sigmaM[type]*sigbg),1);
	    if((dCandMass->at(icand)-dCandChildMass->at(icand)) > lowball && (dCandMass->at(icand)-dCandChildMass->at(icand)) < hiball){
	      checkDCjets[type]=1;
	      proximity[ipt][type]->Fill(djetR);
	      //if(abs(dCandGenMatchPdg->at(icand))==pdgs[type] && hasGenDwithKpiIdx->at(ijet)[type]!=dCandGenMatchIdx->at(icand)) proximityComb[ipt][type]->Fill(djetR);
	      proximityEta[ipt][type]->Fill(deta);
	      proximityPhi[ipt][type]->Fill(dphi);
	      if(abs(dCandGenMatchPdg->at(icand)) != pdgs[type] && abs(dCandGenMatchDoubleCountedPdg->at(icand)) != pdgs[type]) bgProximityTruth[ipt][type]->Fill(djetR);
	      if(djetR<0) cout << "oh no! djetR is negative!! " << djetR << endl;
	      // if(abs(dCandGenMatchPdg->at(icand))==pdgs[type] && hasGenDwithKpiIdx->at(ijet)[type]==dCandGenMatchIdx->at(icand)) proximityTruth[ipt][type]->Fill(djetR);
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
	    else if(((dCandMass->at(icand)-dCandChildMass->at(icand)) < lowballbg && dCandMass->at(icand)-dCandChildMass->at(icand) > lowMassBin[type]) || (dCandMass->at(icand)-dCandChildMass->at(icand)) > hiballbg){
	      checkbgDCjets[type]=1;
	      //bgCjetSpec[type]->Fill(jtpt->at(ijet));
	      bgProximity[ipt][type]->Fill(djetR);
	      bgProximityEta[ipt][type]->Fill(deta);
	      bgProximityPhi[ipt][type]->Fill(dphi);
	    }
	    if((dCandMass->at(icand)-dCandChildMass->at(icand)) < lowballbg){
	      if(abs(dCandGenMatchPdg->at(icand)) != pdgs[type]){
		bgProximityLeftSide[ipt][type]->Fill(djetR);
	      }
	    }
	    if((dCandMass->at(icand)-dCandChildMass->at(icand)) > hiballbg){
	      if(abs(dCandGenMatchPdg->at(icand)) != pdgs[type]){
		bgProximityRightSide[ipt][type]->Fill(djetR);
	      }
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
		if(abs(dCandGenMatchPdg->at(icand))==pdgs[type] && dCandGenMatchIdx->at(icand)==hasGenDwithKpiIdx->at(ijet)[type]){
		  truthAssns[ipt][0][type]->Fill(dCandMass->at(icand));
		}
		if(abs(dCandGenMatchPdg->at(icand))==pdgs[type] && dCandGenMatchIdx->at(icand)!=hasGenDwithKpiIdx->at(ijet)[type]){
		  combAssns[ipt][0][type]->Fill(dCandMass->at(icand));
		}
	      }
	    }
	    if(djetR<closestDR[type]){ closestDR[type]=djetR; closestDRidx[type]=icand; }
	    double lowball = snapToBin(assns[ipt][0][type],tightSignalCutLow[type],0);
	    double hiball = snapToBin(assns[ipt][0][type],tightSignalCutHi[type],1);
	    double lowballbg = snapToBin(assns[ipt][0][type], MesonMass[type]-(sigmaM[type]*sigbg),0);
	    double hiballbg = snapToBin(assns[ipt][0][type], MesonMass[type]+(sigmaM[type]*sigbg),1);
	    if(dCandMass->at(icand) > lowball && dCandMass->at(icand) < hiball){
	      proximity[ipt][type]->Fill(djetR);
	      //if(abs(dCandGenMatchPdg->at(icand))==pdgs[type] && hasGenDwithKpiIdx->at(ijet)[type]!=dCandGenMatchIdx->at(icand)) proximityComb[ipt][type]->Fill(djetR);
	      proximityEta[ipt][type]->Fill(deta);
	      proximityPhi[ipt][type]->Fill(dphi);
	      if(abs(dCandGenMatchPdg->at(icand)) != pdgs[type]  && abs(dCandGenMatchDoubleCountedPdg->at(icand)) != pdgs[type]) bgProximityTruth[ipt][type]->Fill(djetR);
	      checkDCjets[type]=1;
	      // if(abs(dCandGenMatchPdg->at(icand))==pdgs[type] && hasGenDwithKpiIdx->at(ijet)[type]==dCandGenMatchIdx->at(icand)) proximityTruth[ipt][type]->Fill(djetR);
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
	    else if((dCandMass->at(icand) < lowballbg && dCandMass->at(icand) > lowMassBin[type])  || dCandMass->at(icand) > hiballbg){ //only grab bg distributions from the region outside 5-sigma!
	      checkbgDCjets[type]=1;
	      bgProximity[ipt][type]->Fill(djetR);
	      bgProximityEta[ipt][type]->Fill(deta);
	      bgProximityPhi[ipt][type]->Fill(dphi);
	    }

	    if((dCandMass->at(icand)) < lowballbg){
	      bgProximityLeftSide[ipt][type]->Fill(djetR);
	      if(abs(dCandGenMatchPdg->at(icand)) != pdgs[type]){
		bgProximityLeftSideTruth[ipt][type]->Fill(djetR);
	      }
	    }
	    if((dCandMass->at(icand)) > hiballbg){
	      bgProximityRightSide[ipt][type]->Fill(djetR);
	      if(abs(dCandGenMatchPdg->at(icand)) != pdgs[type]){
		bgProximityRightSideTruth[ipt][type]->Fill(djetR);
	      }
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
	
	if(checkDCjets[type]){
	  deta = fabs(dCandEta->at(icand)-jteta->at(ijet));
	  dphi = fabs(dCandPhi->at(icand)-jtphi->at(ijet));
	  if(candDR[type]>sqrt(pow(deta,2)+pow(dphi,2))){
	    candDR[type] = sqrt(pow(deta,2)+pow(dphi,2));
	    for(unsigned int igen=0; igen<genPt->size(); igen++){
	      if(genIdx->at(igen) == dCandGenMatchIdx->at(icand) && abs(genPdg->at(igen))==pdgs[type]){
		assocPt[type] = genPt->at(igen);
	      }
	    }
	  }
	}
	if(checkbgDCjets[type]){
	  deta = fabs(dCandEta->at(icand)-jteta->at(ijet));
	  dphi = fabs(dCandPhi->at(icand)-jtphi->at(ijet));
	  if(bgcandDR[type]>sqrt(pow(deta,2)+pow(dphi,2))){
	    bgcandDR[type] = sqrt(pow(deta,2)+pow(dphi,2));
	    for(unsigned int igen=0; igen<genPt->size(); igen++){
	      if(genIdx->at(igen) == dCandGenMatchIdx->at(icand) && abs(genPdg->at(igen))==pdgs[type]){
		assocPt[type] = genPt->at(igen);
	      }
	    }
	  }
	}
      } //close candidate loop
      
      for(int itype=0; itype<5; itype++){
	if(closestDRidx[itype]>=0){
	  if(itype==0){
	    double lowball = snapToBin(assns[ipt][0][itype],tightSignalCutLow[itype],0);
	    double hiball = snapToBin(assns[ipt][0][itype],tightSignalCutHi[itype],1);
	    if((dCandMass->at(closestDRidx[itype])-dCandChildMass->at(closestDRidx[itype])) > lowball && (dCandMass->at(closestDRidx[itype])-dCandChildMass->at(closestDRidx[itype])) < hiball){
	      Jetproximity[ipt][itype]->Fill(closestDR[itype]);
	      if(abs(dCandGenMatchPdg->at(closestDRidx[itype]))==pdgs[itype] && hasGenDwithKpiIdx->at(ijet)[itype]==dCandGenMatchIdx->at(closestDRidx[itype])) proximityTruth[ipt][itype]->Fill(closestDR[itype]);
	      else if(abs(dCandGenMatchPdg->at(closestDRidx[itype]))==pdgs[itype]) proximityComb[ipt][itype]->Fill(closestDR[itype]);
	    }
	    else{
	      bgJetproximity[ipt][itype]->Fill(closestDR[itype]);
	    }
	  }
	  else{
	    double lowball = snapToBin(assns[ipt][0][itype],tightSignalCutLow[itype],0);
	    double hiball = snapToBin(assns[ipt][0][itype],tightSignalCutHi[itype],1);
	    if(dCandMass->at(closestDRidx[itype]) > lowball && dCandMass->at(closestDRidx[itype]) < hiball){
	      Jetproximity[ipt][itype]->Fill(closestDR[itype]);
	      if(abs(dCandGenMatchPdg->at(closestDRidx[itype]))==pdgs[itype] && hasGenDwithKpiIdx->at(ijet)[itype]==dCandGenMatchIdx->at(closestDRidx[itype])) proximityTruth[ipt][itype]->Fill(closestDR[itype]);
	      else if(abs(dCandGenMatchPdg->at(closestDRidx[itype]))==pdgs[itype]) proximityComb[ipt][itype]->Fill(closestDR[itype]);
	    }
	    else{
	      bgJetproximity[ipt][itype]->Fill(closestDR[itype]);
	    }
	  }
	}
      }
	    
      //fill jet spectra based on whats fg and bg
      for(int itype=0; itype<5; itype++){
	if(candDR[itype]<drCut){
	  //CjetSpec[itype]->Fill(refpt->at(ijet),assocPt[itype],weight/feffline[itype]->Eval(assocPt[itype]));
	  double efficiency = cJetEffvsPt_step2[itype]->GetBinContent(cJetEffvsPt_step2[itype]->FindBin(assocPt[itype]));
	  CjetSpec[itype]->Fill(refpt->at(ijet),1.);///efficiency);
	  DMesonSpec[itype]->Fill(assocPt[itype],weight);///feffline[itype]->Eval(assocPt[itype]));
	}
	else if(bgcandDR[itype]<drCut){
	  bgCjetSpec[itype]->Fill(refpt->at(ijet),1.);///feffline[itype]->Eval(assocPt[itype]));
	}
      }

      //do efficiency for D branching ratio (step 3)
      if(TMath::Abs(refparton_flavorForB->at(ijet))==4){
	for(int ij=0; ij<5; ij++){
	  bool step3fg=false;
	  bool step3bg=false;
	  for(unsigned int igen=0; igen<genPt->size(); igen++){
	    if(hasGenDIdx->at(ijet)[ij]==genIdx->at(igen)){
	      step3bg=true;
	    }
	    if(hasGenDwithKpiIdx->at(ijet)[ij]==genIdx->at(igen)){
	      step3fg=true;
	    }
	  }
	  if(step3bg){
	    cJetEffvsPt_bg_step3[ij]->Fill(refpt->at(ijet),1.);
	    if(step3fg){
	      cJetEffvsPt_step3[ij]->Fill(refpt->at(ijet),1.);
	    }
	    else{
	      cJetEffvsPt_rej_step3[ij]->Fill(refpt->at(ijet),1.);
	    }
	  }
	}
      }    
      
      //do efficiency for c->D fragmentation fraction
      if(TMath::Abs(refparton_flavorForB->at(ijet))==4){
	for(int ij=0; ij<5; ij++){
	  cJetEffvsPt_bg_step4[ij]->Fill(refpt->at(ijet),1.);
	  if(hasGenD->at(ijet)[ij]){
	    cJetEffvsPt_step4[ij]->Fill(refpt->at(ijet),1.);
	  }
	  else cJetEffvsPt_rej_step4[ij]->Fill(refpt->at(ijet),1.);
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
  
  TLatex *twrite3[cTypes];
  TH1D *qcd_jet[cTypes];
  //TFile *fqcd = new TFile("dCandidateToJetAssns_officialQCDJetOnly_RevCuts_JetAssn-11-11.root");
  TCanvas *ctemp = new TCanvas("ctemp","",800,800);
  ctemp->Divide(3,2);
  for(int ij=0; ij<5; ij++){
    //qcd_jet[ij] = (TH1D*)fqcd->Get(Form("cJetEffvsPt_step2_type%d",ij))->Clone(Form("qcd_jet_%d",ij));
    ctemp->cd(ij+1);
    cJetEffvsPt_step2[ij]->SetXTitle("gen D pT");
    cJetEffvsPt_step2[ij]->SetYTitle("Reco Efficiency");
    cJetEffvsPt_step2[ij]->SetMaximum(0.9);
    cJetEffvsPt_step2[ij]->Draw();
    //qcd_jet[ij]->SetMarkerColor(kBlue+2);
    //qcd_jet[ij]->Draw("same");
    cJetEffvsPt_bg_step2[ij]->SetMarkerColor(2);
    //cJetEffvsPt_bg_step2[ij]->Draw("same");
    feffline[ij]->Draw("same");
    twrite3[ij] = new TLatex(20,gPad->GetUymax()*0.65,Form("Meson Type: %s",dtypes[ij].c_str()));
    twrite3[ij]->Draw("same");
  }
  
  fout->cd();
  for(int i=0; i<ptBins; i++){
    for(int j=0; j<djetBins; j++){
      for(int k=0; k<cTypes; k++){
      assns[i][j][k]->Write();
      truthAssns[i][j][k]->Write();
      combAssns[i][j][k]->Write();
      assnsUnweight[i][j][k]->Write();
      bgAssnsUnweight[i][j][k]->Write();
      bgAssns[i][j][k]->Write();
      dcAssns[i][j][k]->Write();
      nondcAssns[i][j][k]->Write();
      bgRotate[i][j][k]->Write();
      }
    }
  }

  for(int i=0; i<ptBins; i++){
    for(int k=0; k<cTypes; k++){
      proximity[i][k]->Write();
      Jetproximity[i][k]->Write();
      proximityComb[i][k]->Write();
      proximityTruth[i][k]->Write();
      proximityEta[i][k]->Write();
      proximityPhi[i][k]->Write();
      bgProximity[i][k]->Write();
      bgProximityTruth[i][k]->Write();
      bgJetproximity[i][k]->Write();
      bgProximityEta[i][k]->Write();
      bgProximityPhi[i][k]->Write();
      bgProximityLeftSide[i][k]->Write();
      bgProximityRightSide[i][k]->Write();
      bgProximityLeftSideTruth[i][k]->Write();
      bgProximityRightSideTruth[i][k]->Write();
    }
  }
  for(int i=0; i<3; i++){
    mesonMults[i]->Write();
    mesonMultsNoDstar[i]->Write();
  }
  
  allJetSpec->Write();
  for(int k=0; k<cTypes; k++){
    CjetSpec[k]->Write();
    bgCjetSpec[k]->Write();
    DMesonSpec[k]->Write();
  }

  //cJetEffvsPt->Divide(cJetEffvsPt,cJetEffvsPt_bg,1,1,"B");
  //cJetPurvsPt->Divide(cJetPurvsPt,cJetPurvsPt_bg,1,1,"B");
  cJetEffvsPt->SetYTitle( (Dlabels[TMath::Max(isMC-1,5)]+(string)Form("dR<%g C-Jet Efficiency",djetR_testCut)).c_str());
  cJetEffvsPt->SetXTitle("c-Jet p_{T}");
  cJetPurvsPt->SetYTitle( (Dlabels[TMath::Max(isMC-1,5)]+(string)Form("dR<%g C-Jet Purity",djetR_testCut)).c_str());
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
