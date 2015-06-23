#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TMath.h"
#include <vector>
#include <iostream>
#include <algorithm>

//to get the trained discriminator information from TMVA:
#include "/Users/kjung/charmJets/pPb/tmvaTesting/weights/TMVAClassification_BDTG_CvL_medCuts.class.C"
//#include "/Users/kjung/charmJets/pPb/tmvaTesting/weights/TMVAClassification_BDTG_CvL_noSvtx.class.C"
#include "/Users/kjung/charmJets/pPb/tmvaTesting/weights/TMVAClassification_BDTG_BvC_medCuts.class.C"
//#include "/Users/kjung/charmJets/pPb/tmvaTesting/weights/TMVAClassification_BDTG_BvC_noSvtx.class.C"
//#include "/Users/kjung/charmJets/pPb/tmvaTesting/weights/TMVAClassification_RuleFit.class.C"
//end

#ifdef __CINT__
#pragma link C++ class vector<vector<double> >;
#else
template class std::vector<std::vector<double> >;
#endif

using std::vector;
using std::endl;
using std::cout;

int isMC=false;

double findDR(int isMC, vector<double> *dCandPt, vector<double> *dCandEta, vector<double> *dCandPhi, vector<double> *dCandMass, 
  vector<double> *dCandChildMass, vector<double> *dCandCharge1, vector<double> *dCandCharge2, vector<double> *dCandType, 
  vector<int> *dCandGenMatchedPdg, double jteta, double jtphi, double &dMass, double &dChildMass, double &dType, int &dPDG){

  double drRet = 999;
  bool sigCheck = false;

  //only check dR against signal D*'s and D0's for now
  for(unsigned int i=0; i<dCandPt->size(); i++){
    if(dCandCharge1->at(i) != dCandCharge2->at(i)){
      if(dCandType->at(i)==1 && dCandMass->at(i)-dCandChildMass->at(i)>0.14 && dCandMass->at(i)-dCandChildMass->at(i)<0.155 &&
         dCandChildMass->at(i)>1.8 && dCandChildMass->at(i)<1.92){
	//dMass = dCandMass->at(i);
	//dChildMass = dCandChildMass->at(i);
	//dType = dCandType->at(i);
         sigCheck = true;
       }
       else if(dCandType->at(i)==2 && dCandMass->at(i)>1.8 && dCandMass->at(i)<1.92){
	//dMass = dCandMass->at(i);
	//dType = dCandType->at(i);
        sigCheck = true;
      }
    }

    if(sigCheck){
      double drTest = sqrt(pow(dCandEta->at(i)-jteta,2)+pow(dCandPhi->at(i)-jtphi,2));
      if(drTest<drRet){
       drRet = drTest;
       dMass = dCandMass->at(i);
       dChildMass = dCandChildMass->at(i);
       dType = dCandType->at(i);
       if(isMC) dPDG = dCandGenMatchedPdg->at(i);
     }
   }

 }
 if(!sigCheck){
   dMass = -999;
   dChildMass = -999;
   dType = -999;
   dPDG = -999;
  }
  return drRet;
}

double eventFilter(double svtxmcorr, double svtxdl, double svtxpt, vector<double> *dCandMass, vector<double> *dCandType, 
  vector<double> *dCandChildMass, double jteta, double jtphi, vector<double> *dCandPt, vector<double> *dCandEta, 
  vector<double> *dCandPhi, vector<double> *ffls3d, vector<double> *falpha0, vector<double> *fprob0, 
  vector<int> *dCandGenMatchPdg, TH1D *bgAdjuster){

  //cout << "isMC: " << isMC << endl;

  bool doDR = true;
  bool checkDMeson=false;
  bool hasRealD=false;
  double weightToRet=0;
  double maxMesonPt=0;

  double ptCuts[11] = {3.5, 4.5, 5.5, 7.0, 9.0, 11.0, 13.0, 16.0, 20., 28., 40.};
  double d0cuts[10] = {3.81, 3.81, 3.91, 3.92, 3.77, 3.57, 3.37, 2.88, 2.71, 2.63};
  double alphaCuts[10] = {0.056, 0.056, 0.048, 0.060, 0.063, 0.057, 0.053, 0.059, 0.059, 0.042};
  double vtxProbCuts[10] = {0.267, 0.267, 0.209, 0.113, 0.093, 0.057, 0.75, 0.025, 0.68, 0.018};

  //tight cuts = [1.80-1.92]
  //med cuts = [1.76-1.96]
  //loose cuts = [1.7-2.05]

  double d0MassLow = 1.80;
  double d0MassHigh = 1.92;

  TRandom3 *tt = new TRandom3();
  double random = tt->Rndm();

  //if(svtxmcorr > 0.4 && (svtxdl/svtxpt)<0.15){
  for(unsigned int icand=0; icand<dCandMass->size(); icand++){
    if(dCandPt->at(icand)<3.5) continue;
    int j=0;
    while((dCandPt->at(icand)>ptCuts[j]) && (j<10)) j++;
    j--;
      //cout << "dCandPt: "<< dCandPt->at(icand) << " bin: "<< ptCuts[j] << " to " << ptCuts[j+1] << endl;

      //if(fabs(dCandEta->at(icand))<1.0){
    if(falpha0->at(icand)<alphaCuts[j] && fprob0->at(icand)>vtxProbCuts[j] && ffls3d->at(icand)>d0cuts[j]){

      double dMassDiff = dCandMass->at(icand)-dCandChildMass->at(icand);
      if(dCandMass->at(icand) > 1.7 && (dCandType->at(icand)==1 && dMassDiff>0.144 && dMassDiff<0.147 && dCandChildMass->at(icand)>d0MassLow && dCandChildMass->at(icand)<d0MassHigh)){
        if(!doDR || sqrt(pow(dCandEta->at(icand)-jteta,2)+pow(dCandPhi->at(icand)-jtphi,2))<0.5){
          checkDMeson=true;
          if(isMC){
            if(abs(dCandGenMatchPdg->at(icand))==413) hasRealD=true;
            if(dCandPt->at(icand)>maxMesonPt){ 
              maxMesonPt = dCandPt->at(icand);
              weightToRet = bgAdjuster->GetBinContent(bgAdjuster->FindBin(dCandPt->at(icand)));
            }
          }
        }
      }
      else if(dCandMass->at(icand)>d0MassLow && dCandType->at(icand)==2 && dCandMass->at(icand)<d0MassHigh){
        if(!doDR || sqrt(pow(dCandEta->at(icand)-jteta,2)+pow(dCandPhi->at(icand)-jtphi,2))<0.5){
          checkDMeson=true;
          if(isMC){
            if(abs(dCandGenMatchPdg->at(icand))==421) hasRealD=true;
            if(dCandPt->at(icand)>maxMesonPt){ 
              maxMesonPt = dCandPt->at(icand);
              weightToRet = bgAdjuster->GetBinContent(bgAdjuster->FindBin(dCandPt->at(icand)));
            }
          }
        }
      }
    }
  }

 if(!isMC && checkDMeson) return 1;
 else if(!isMC) return 0;

 else if(hasRealD) return 1;
 else if(checkDMeson){ return 1; }
 // weightToRet = bgAdjuster->GetBinContent(bgAdjuster->FindBin(maxMesonPt));
 // if(random<weightToRet) return weightToRet;
 // else return 0;
 // } 
 else return 0;
  
}

void convertVectorToDoubleTree(int MC=0){

  isMC = MC;

  bool refillTMVAvars = true;
  bool dopthat80Only = false;

  TFile *DMesonBGWeight = new TFile("../bgScaleFactors.root");
  TH1D *corrFact = (TH1D*)DMesonBGWeight->Get("purity");
  TH1D *purityMC = (TH1D*)DMesonBGWeight->Get("purityMC");
  corrFact->Divide(purityMC);

  TFile *fout;
  if(dopthat80Only && isMC==1) fout = new TFile("../input/DMesonCJet_QCDJetOnly_pthat80Only_pPbMC_ppReco_akPu3PF_convertToJetTree.root","RECREATE");
  else if(dopthat80Only && isMC==2) fout = new TFile("../input/DMesonCJet_CJetOnly_pthat80Only_pPbMC_ppReco_akPu3PF_convertToJetTree_temp.root","RECREATE");
  else if(isMC==1) fout = new TFile("../input/DMesonCJet_QCDJetOnly_pPbMC_ppReco_akPu3PF_convertToJetTree_tightDCutsWithDR.root","RECREATE");
  else if(isMC==2) fout = new TFile("../input/DMesonCJet_CJetOnly_pPbMC_ppReco_akPu3PF_convertToJetTree_temp.root","RECREATE");
  else fout = new TFile("../input/DMesonCJet_pPbData_ppReco_akPu3PF_convertToJetTree_withLHCbVars_tightDCutsWithDR_part1.root","RECREATE");

  TTree *jets = new TTree("jets","");
  TTree *jetsNoSvtx = new TTree("jetsNoSvtx","");
  TTree *dMesons = new TTree("dMesons","");

  double a_dCandPt, a_dCandEta, a_dCandPhi, a_dCandMass, a_dCandType, a_dCandChildMass, a_dCandCharge1, a_dCandCharge2;
  double a_jtpt, a_jteta, a_jtphi, a_refpt, a_rawpt, a_djetR, a_jtm;
  double a_refparton_flavorForB, a_svtxm, a_svtxmcorr, a_svtxdl, a_svtxpt, a_svtxdls, a_sv2Trkdl, a_sv2Trkdls, a_svtxTrkSumChi2, a_svtxTrkNetCharge, a_svtxNtrkInCone, a_svtxptFrac;
  double a_closestDMass, a_closestDChildMass, a_closestDType, a_ffls3d, a_falpha0, a_fprob;
  int a_closestDPdg;
  int a_genMatch, a_bin, a_run, a_subid, a_evtSelection;
  double a_weight, a_vz, a_bdt, a_bdt_bvc, a_bdt_noVtx, a_bdt_bvc_noVtx, a_ruleFit;
  double a_discr_csvSimple, a_discr_ssvHighEff, a_discr_ssvHighPur, a_discr_cJetHighEff, a_discr_cJetHighPur, a_closestDPt, a_discr_prob, a_discr_tcHighEff, a_discr_tcHighPur;

  int a_nIP, a_svtxntrk, a_jtntrks;
  double a_neutralMax, a_neutralSum, a_chargedMax, a_chargedSum, a_photonMax, a_photonSum, a_muSum, a_eSum;
  double a_svJetDeltaR, a_trackSumJetDeltaR, a_svtxmEnergyFrac;
  double a_trackSip2dSigAboveCharm, a_trackSip2dValAboveCharm, a_trackSip3dValAboveCharm, a_trackSip3dSigAboveCharm;
  double a_trackIP2dSig[4], a_trackIP3dSig[4], a_trackIP2d[4], a_trackIP3d[4], a_ipProb0[4], a_ipPt[4], a_trackPtRel[4], a_trackPPar[4], a_trackPParRatio[4], a_ipDist2Jet[4], a_ipClosest2Jet[4], a_trackDeltaR[4], a_trackPtRatio[4];
  double pthat;

  dMesons->Branch("dCandPt",&a_dCandPt);
  dMesons->Branch("dCandEta",&a_dCandEta);
  dMesons->Branch("dCandPhi",&a_dCandPhi);
  dMesons->Branch("dCandMass",&a_dCandMass);
  dMesons->Branch("dCandType",&a_dCandType);
  dMesons->Branch("dCandChildMass",&a_dCandChildMass);
  dMesons->Branch("dCandCharge1",&a_dCandCharge1);
  dMesons->Branch("dCandCharge2",&a_dCandCharge2);
  dMesons->Branch("weight",&a_weight);
  dMesons->Branch("fprob",&a_fprob);
  dMesons->Branch("falpha0",&a_falpha0);
  dMesons->Branch("ffls3d",&a_ffls3d);
  
  jets->Branch("jtpt",&a_jtpt);
  jets->Branch("jteta",&a_jteta);
  jets->Branch("jtphi",&a_jtphi);
  jets->Branch("rawpt",&a_rawpt);
  jets->Branch("jtm",&a_jtm);
  jets->Branch("jtntrks",&a_jtntrks);
  jets->Branch("svtxmEnergyFrac",&a_svtxmEnergyFrac);
  if(isMC){
    jets->Branch("refpt",&a_refpt);
    jets->Branch("refparton_flavorForB",&a_refparton_flavorForB);
    jets->Branch("genMatch",&a_genMatch);
    jets->Branch("subid",&a_subid);
  }
  jets->Branch("djetR",&a_djetR);
  jets->Branch("discr_csvSimple",&a_discr_csvSimple);
  jets->Branch("discr_ssvHighEff",&a_discr_ssvHighEff);
  jets->Branch("discr_ssvHighPur",&a_discr_ssvHighPur);
  if(isMC>1) jets->Branch("discr_cJetHighEff",&a_discr_cJetHighEff);
  if(isMC>1) jets->Branch("discr_cJetHighPur",&a_discr_cJetHighPur);
  jets->Branch("discr_tcHighEff",&a_discr_tcHighEff);
  jets->Branch("discr_tcHighPur",&a_discr_tcHighPur);
  jets->Branch("discr_prob",&a_discr_prob);
  jets->Branch("closestDPt",&a_closestDPt);
  jets->Branch("closestDMass",&a_closestDMass);
  jets->Branch("closestDChildMass",&a_closestDChildMass);
  jets->Branch("closestDType",&a_closestDType);
  if(isMC) jets->Branch("closestDPdg",&a_closestDPdg);
  jets->Branch("weight",&a_weight);
  jets->Branch("run",&a_run);
  jets->Branch("evtSelection",&a_evtSelection);
  jets->Branch("vz",&a_vz);
  jets->Branch("bin",&a_bin);
  jets->Branch("pthat",&pthat);
  jets->Branch("svtxm",&a_svtxm);
  jets->Branch("svtxmcorr",&a_svtxmcorr);
  jets->Branch("svtxdl",&a_svtxdl);
  jets->Branch("svtxdls",&a_svtxdls);
  jets->Branch("svtxntrk",&a_svtxntrk);
  jets->Branch("svtxpt",&a_svtxpt);
  jets->Branch("svtxptFrac",&a_svtxptFrac);
  jets->Branch("sv2Trkdl",&a_sv2Trkdl);
  jets->Branch("sv2Trkdls",&a_sv2Trkdls);
  jets->Branch("svtxTrkSumChi2",&a_svtxTrkSumChi2);
  jets->Branch("svtxTrkNetCharge",&a_svtxTrkNetCharge);
  jets->Branch("svtxNtrkInCone",&a_svtxNtrkInCone);
  jets->Branch("nIP",&a_nIP);
  for(int i=0; i<4; i++){
    jets->Branch(Form("ipPt_%d",i),&a_ipPt[i]);
    jets->Branch(Form("trackIP2dSig_%d",i),&a_trackIP2dSig[i]);
    jets->Branch(Form("trackIP3dSig_%d",i),&a_trackIP3dSig[i]);
    jets->Branch(Form("trackIP2d_%d",i),&a_trackIP2d[i]);
    jets->Branch(Form("trackIP3d_%d",i),&a_trackIP3d[i]);
    jets->Branch(Form("ipProb0_%d",i),&a_ipProb0[i]);
    jets->Branch(Form("trackPtRel_%d",i),&a_trackPtRel[i]);
    jets->Branch(Form("trackPtRatio_%d",i),&a_trackPtRatio[i]);
    jets->Branch(Form("trackPPar_%d",i),&a_trackPPar[i]);
    jets->Branch(Form("trackPParRatio_%d",i),&a_trackPParRatio[i]);
    jets->Branch(Form("trackJetDist_%d",i),&a_ipDist2Jet[i]);
    jets->Branch(Form("trackDecayLenVal_%d",i),&a_ipClosest2Jet[i]);
    jets->Branch(Form("trackDeltaR_%d",i),&a_trackDeltaR[i]);
  }
  jets->Branch("trackSip2dSigAboveCharm",&a_trackSip2dSigAboveCharm);
  jets->Branch("trackSip3dSigAboveCharm",&a_trackSip3dSigAboveCharm);
  jets->Branch("trackSip2dValAboveCharm",&a_trackSip2dValAboveCharm);
  jets->Branch("trackSip3dValAboveCharm",&a_trackSip3dValAboveCharm);
  jets->Branch("svJetDeltaR",&a_svJetDeltaR);
  jets->Branch("trackSumJetDeltaR",&a_trackSumJetDeltaR);
  
  jets->Branch("neutralMax",&a_neutralMax);
  jets->Branch("neutralSum",&a_neutralSum);
  jets->Branch("chargedMax",&a_chargedMax);
  jets->Branch("chargedSum",&a_chargedSum);
  jets->Branch("photonMax",&a_photonMax);
  jets->Branch("photonSum",&a_photonSum);
  jets->Branch("muSum",&a_muSum);
  jets->Branch("eSum",&a_eSum);
  
  if(refillTMVAvars) jets->Branch("BDTG_BCvL",&a_bdt);
  if(refillTMVAvars) jets->Branch("BDTG_BvC",&a_bdt_bvc);
  //if(refillTMVAvars) jets->Branch("ruleFit",&a_ruleFit);

//****** START NO SVTX TREE ******//
  jetsNoSvtx->Branch("jtpt",&a_jtpt);
  jetsNoSvtx->Branch("jteta",&a_jteta);
  jetsNoSvtx->Branch("jtphi",&a_jtphi);
  jetsNoSvtx->Branch("rawpt",&a_rawpt);
  jetsNoSvtx->Branch("jtm",&a_jtm);
  jetsNoSvtx->Branch("jtntrks",&a_jtntrks);
  jetsNoSvtx->Branch("svtxmEnergyFrac",&a_svtxmEnergyFrac);
  if(isMC){
    jetsNoSvtx->Branch("refpt",&a_refpt);
    jetsNoSvtx->Branch("refparton_flavorForB",&a_refparton_flavorForB);
    jetsNoSvtx->Branch("genMatch",&a_genMatch);
    jetsNoSvtx->Branch("subid",&a_subid);
  }
  jetsNoSvtx->Branch("djetR",&a_djetR);
  jetsNoSvtx->Branch("discr_csvSimple",&a_discr_csvSimple);
  jetsNoSvtx->Branch("discr_ssvHighEff",&a_discr_ssvHighEff);
  jetsNoSvtx->Branch("discr_ssvHighPur",&a_discr_ssvHighPur);
  if(isMC>1) jetsNoSvtx->Branch("discr_cJetHighEff",&a_discr_cJetHighEff);
  if(isMC>1) jetsNoSvtx->Branch("discr_cJetHighPur",&a_discr_cJetHighPur);
  jetsNoSvtx->Branch("discr_tcHighEff",&a_discr_tcHighEff);
  jetsNoSvtx->Branch("discr_tcHighPur",&a_discr_tcHighPur);
  jetsNoSvtx->Branch("discr_prob",&a_discr_prob);
  jetsNoSvtx->Branch("closestDPt",&a_closestDPt);
  jetsNoSvtx->Branch("closestDMass",&a_closestDMass);
  jetsNoSvtx->Branch("closestDChildMass",&a_closestDChildMass);
  jetsNoSvtx->Branch("closestDType",&a_closestDType);
  if(isMC) jetsNoSvtx->Branch("closestDPdg",&a_closestDPdg);
  jetsNoSvtx->Branch("weight",&a_weight);
  jetsNoSvtx->Branch("run",&a_run);
  jetsNoSvtx->Branch("evtSelection",&a_evtSelection);
  jetsNoSvtx->Branch("vz",&a_vz);
  jetsNoSvtx->Branch("bin",&a_bin);
  jetsNoSvtx->Branch("pthat",&pthat);
  jetsNoSvtx->Branch("svtxm",&a_svtxm);
  jetsNoSvtx->Branch("svtxmcorr",&a_svtxmcorr);
  jetsNoSvtx->Branch("svtxdl",&a_svtxdl);
  jetsNoSvtx->Branch("svtxdls",&a_svtxdls);
  jetsNoSvtx->Branch("svtxntrk",&a_svtxntrk);
  jetsNoSvtx->Branch("svtxpt",&a_svtxpt);
  jetsNoSvtx->Branch("svtxptFrac",&a_svtxptFrac);
  jetsNoSvtx->Branch("sv2Trkdl",&a_sv2Trkdl);
  jetsNoSvtx->Branch("sv2Trkdls",&a_sv2Trkdls);
  jetsNoSvtx->Branch("svtxTrkSumChi2",&a_svtxTrkSumChi2);
  jetsNoSvtx->Branch("svtxTrkNetCharge",&a_svtxTrkNetCharge);
  jetsNoSvtx->Branch("svtxNtrkInCone",&a_svtxNtrkInCone);
  jetsNoSvtx->Branch("nIP",&a_nIP);
  for(int i=0; i<4; i++){
    jetsNoSvtx->Branch(Form("ipPt_%d",i),&a_ipPt[i]);
    jetsNoSvtx->Branch(Form("trackIP2dSig_%d",i),&a_trackIP2dSig[i]);
    jetsNoSvtx->Branch(Form("trackIP3dSig_%d",i),&a_trackIP3dSig[i]);
    jetsNoSvtx->Branch(Form("trackIP2d_%d",i),&a_trackIP2d[i]);
    jetsNoSvtx->Branch(Form("trackIP3d_%d",i),&a_trackIP3d[i]);
    jetsNoSvtx->Branch(Form("ipProb0_%d",i),&a_ipProb0[i]);
    jetsNoSvtx->Branch(Form("trackPtRel_%d",i),&a_trackPtRel[i]);
    jetsNoSvtx->Branch(Form("trackPtRatio_%d",i),&a_trackPtRatio[i]);
    jetsNoSvtx->Branch(Form("trackPPar_%d",i),&a_trackPPar[i]);
    jetsNoSvtx->Branch(Form("trackPParRatio_%d",i),&a_trackPParRatio[i]);
    jetsNoSvtx->Branch(Form("trackJetDist_%d",i),&a_ipDist2Jet[i]);
    jetsNoSvtx->Branch(Form("trackDecayLenVal_%d",i),&a_ipClosest2Jet[i]);
    jetsNoSvtx->Branch(Form("trackDeltaR_%d",i),&a_trackDeltaR[i]);
  }
  jetsNoSvtx->Branch("trackSip2dSigAboveCharm",&a_trackSip2dSigAboveCharm);
  jetsNoSvtx->Branch("trackSip3dSigAboveCharm",&a_trackSip3dSigAboveCharm);
  jetsNoSvtx->Branch("trackSip2dValAboveCharm",&a_trackSip2dValAboveCharm);
  jetsNoSvtx->Branch("trackSip3dValAboveCharm",&a_trackSip3dValAboveCharm);
  jetsNoSvtx->Branch("svJetDeltaR",&a_svJetDeltaR);
  jetsNoSvtx->Branch("trackSumJetDeltaR",&a_trackSumJetDeltaR);
  
  jetsNoSvtx->Branch("neutralMax",&a_neutralMax);
  jetsNoSvtx->Branch("neutralSum",&a_neutralSum);
  jetsNoSvtx->Branch("chargedMax",&a_chargedMax);
  jetsNoSvtx->Branch("chargedSum",&a_chargedSum);
  jetsNoSvtx->Branch("photonMax",&a_photonMax);
  jetsNoSvtx->Branch("photonSum",&a_photonSum);
  jetsNoSvtx->Branch("muSum",&a_muSum);
  jetsNoSvtx->Branch("eSum",&a_eSum);
  
  if(refillTMVAvars) jetsNoSvtx->Branch("BDTG_BCvL",&a_bdt_noVtx);
  if(refillTMVAvars) jetsNoSvtx->Branch("BDTG_BvC",&a_bdt_bvc_noVtx);


  TFile *fin;
  TChain *t;
  TChain *finData;
  if(isMC==1) fin = new TFile("../input/DMesonCJet_QCDJetOnly_Merged_addLHCbVars_pPbMC_centReweight_ppReco_akPu3PF.root");
  else if(isMC==2) fin = new TFile("../input/DMesonCJet_CJetOnly_pPbMC_ppReco_akPu3PF_1059.root");
  else{
    finData = new TChain("ct");
    finData->Add("../input/DMesonCJet_FullMerge_addLHCbVars_NoJetTrgCut_pPbdata_ppReco_akPu3PF_part1.root");
    finData->Add("../input/DMesonCJet_FullMerge_addLHCbVars_NoJetTrgCut_pPbdata_ppReco_akPu3PF_part2.root");
    finData->Add("../input/DMesonCJet_FullMerge_addLHCbVars_NoJetTrgCut_pPbdata_ppReco_akPu3PF_part3.root");
    cout << "chain complete!" << endl;
  }
  if(isMC) t = (TChain*)fin->Get("ct");
  else t = finData;
  if(!isMC) cout << "tree obtained from chain!" << endl;
  
  double weight;
  float vz;
  int bin, run, subid;
  vector<double> *dCandPt=0, *dCandMass=0, *dCandEta=0, *dCandPhi=0, *jtpt=0, *jteta=0, *jtphi=0, *refpt=0, *rawpt=0, *djetR=0, *jtm=0;
  vector<double> *dCandType=0, *dCandChildMass=0, *dCandCharge1=0, *dCandCharge2=0, *ffls3d=0, *falpha0=0, *fprob=0;
  vector<double> *closestDMass=0, *closestDType=0, *closestDChildMass=0;
  vector<int> *refparton_flavorForB=0;
  vector<int> *genMatch=0, *dCandGenMatchPdg=0, *jtntrks=0, *svtxTrkNetCharge=0, *svtxNtrkInCone=0;
  vector<double> *discr_csvSimple=0, *discr_ssvHighEff=0, *discr_ssvHighPur=0, *discr_cJetHighEff=0, *discr_cJetHighPur=0, *discr_tcHighEff=0, *discr_tcHighPur=0, *closestDPt=0, *discr_prob=0, *svtxm=0, *svtxdl=0, *svtxmcorr=0, *sv2Trkdl=0, *sv2Trkdls=0, *svtxTrkSumChi2=0;
  vector<bool> *isGSP=0;
  vector<double> *svtxdls=0, *svtxpt=0, *neutralMax=0, *neutralSum=0, *chargedMax=0, *chargedSum=0, *photonMax=0, *photonSum=0, *muSum=0, *eSum=0, *trackSip2dSigAboveCharm=0, *trackSip3dSigAboveCharm=0, *trackSip2dValAboveCharm=0, *trackSip3dValAboveCharm=0, *trackSumJetDeltaR=0, *svJetDeltaR=0;
  vector<int> *nIP=0, *svtxntrk=0;
  vector<double> *trackIP2dSig=0, *trackIP3dSig=0, *trackIP2d=0, *trackIP3d=0, *ipProb0=0, *ipPt=0, *trackPtRel=0, *trackPtRatio=0, *trackPPar=0, *trackDeltaR=0, *trackPParRatio=0, *ipDist2Jet=0, *ipClosest2Jet=0;

  //if(!isMC) dCandGenMatchPdg = new vector<int>;

  //string inputArr_withVtx[] =  {"djetR", "nIP", "svtxm", "svtxdl", "svtxdls", "svtxntrk", "jteta", "ipProb0_0", "ipPt_0", "trackIP2dSig_0", "trackIP3dSig_0", "trackIP2d_0", "trackIP3d_0", "trackPtRel_0", "trackPPar_0", "trackPParRatio_0", "trackJetDist_0", "trackDecayLenVal_0", "trackDeltaR_0", "trackPtRatio_0", "ipProb0_1", "ipPt_1", "trackIP2dSig_1", "trackIP3dSig_1", "trackIP2d_1", "trackIP3d_1", "trackPtRel_1", "trackPPar_1", "trackPParRatio_1", "trackJetDist_1", "trackDecayLenVal_1", "trackDeltaR_1", "trackPtRatio_1", "ipProb0_2", "ipPt_2", "trackIP2dSig_2", "trackIP3dSig_2", "trackIP2d_2", "trackIP3d_2", "trackPtRel_2", "trackPPar_2", "trackPParRatio_2", "trackJetDist_2", "trackDecayLenVal_2", "trackDeltaR_2", "trackPtRatio_2", "ipProb0_3", "ipPt_3", "trackIP2dSig_3", "trackIP3dSig_3", "trackIP2d_3", "trackIP3d_3", "trackPtRel_3", "trackPPar_3", "trackPParRatio_3", "trackJetDist_3", "trackDecayLenVal_3", "trackDeltaR_3", "trackPtRatio_3", "trackSip2dValAboveCharm", "trackSip3dValAboveCharm", "svJetDeltaR", "trackSumJetDeltaR" };
  string inputArr_withVtx[] = { "svtxptFrac", "nIP", "svtxm", "svtxmEnergyFrac", "svtxpt", "svtxmcorr", "svtxdl", "svtxdls", "svtxntrk", "sv2Trkdl", "svtxTrkSumChi2", "svtxTrkNetCharge", "svtxNtrkInCone", "closestDMass", "closestDType", "closestDPt", "trackIP2dSig_0", "trackIP3dSig_0", "trackIP2d_0", "trackIP3d_0", "trackIP2dSig_1", "trackIP3dSig_1", "trackIP2d_1", "trackIP3d_1", "trackIP2dSig_2", "trackIP3dSig_2", "trackIP2d_2", "trackIP3d_2", "trackJetDist_0", "trackDecayLenVal_0", "trackSumJetDeltaR" };
  //string inputArr_noVtx[] = { "trackIP2dSig_0", "trackIP3dSig_0", "trackIP2d_0", "trackIP3d_0", "trackIP2dSig_1", "trackIP3dSig_1", "trackIP2d_1", "trackIP3d_1", "trackIP2dSig_2", "trackIP3dSig_2", "trackIP2d_2", "trackIP3d_2", "trackJetDist_0", "trackDecayLenVal_0", "svJetDeltaR", "trackSumJetDeltaR" };
  vector<string> inputVars_withVtx;
  for(int i=0; i<31; i++) inputVars_withVtx.push_back(inputArr_withVtx[i]);
  //vector<string> inputVars_noVtx;
  //for(int i=0; i<16; i++) inputVars_noVtx.push_back(inputArr_noVtx[i]);

  ReadBDTG_CvL_withSvtx testBDTG_CvL_withSvtx(inputVars_withVtx);
  ReadBDTG_BvC_withSvtx testBDTG_BvC_withSvtx(inputVars_withVtx);
  //ReadBDTG_CvL_noSvtx testBDTG_CvL_noSvtx(inputVars_noVtx);
  //ReadBDTG_BvC_noSvtx testBDTG_BvC_noSvtx(inputVars_noVtx);

  cout << " checkpoint 1" << endl;
  
  t->SetBranchAddress("dCandPt",&dCandPt);
  t->SetBranchAddress("dCandEta",&dCandEta);
  t->SetBranchAddress("dCandPhi",&dCandPhi);
  t->SetBranchAddress("dCandMass",&dCandMass);
  t->SetBranchAddress("dCandType",&dCandType);
  t->SetBranchAddress("dCandChildMass",&dCandChildMass);
  t->SetBranchAddress("dCandCharge1",&dCandCharge1);
  t->SetBranchAddress("dCandCharge2",&dCandCharge2);
  t->SetBranchAddress("ffls3d",&ffls3d);
  t->SetBranchAddress("falpha0",&falpha0);
  t->SetBranchAddress("fprob",&fprob);
  if(isMC) t->SetBranchAddress("dCandGenMatchPdg",&dCandGenMatchPdg);

  t->SetBranchAddress("jtpt",&jtpt);
  t->SetBranchAddress("jteta",&jteta);
  t->SetBranchAddress("jtphi",&jtphi);
  t->SetBranchAddress("jtm",&jtm);
  t->SetBranchAddress("jtntrks",&jtntrks);
  t->SetBranchAddress("rawpt",&rawpt);
  t->SetBranchAddress("djetR",&djetR);
  t->SetBranchAddress("svtxm",&svtxm);
  t->SetBranchAddress("svtxmcorr",&svtxmcorr);
  t->SetBranchAddress("svtxdl",&svtxdl);
  t->SetBranchAddress("svtxntrk",&svtxntrk);
  t->SetBranchAddress("sv2Trkdl",&sv2Trkdl);
  t->SetBranchAddress("sv2Trkdls",&sv2Trkdls);
  t->SetBranchAddress("svtxTrkSumChi2",&svtxTrkSumChi2);
  t->SetBranchAddress("svtxTrkNetCharge",&svtxTrkNetCharge);
  t->SetBranchAddress("svtxNtrkInCone",&svtxNtrkInCone);
  if(isMC){
    t->SetBranchAddress("refpt",&refpt);
    t->SetBranchAddress("refparton_flavorForB",&refparton_flavorForB);
    t->SetBranchAddress("genMatch",&genMatch);
    t->SetBranchAddress("subid",&subid);
    t->SetBranchAddress("pthat",&pthat);
  }
  t->SetBranchAddress("closestDMass",&closestDMass);
  t->SetBranchAddress("closestDChildMass",&closestDChildMass);
  t->SetBranchAddress("closestDType",&closestDType);
  t->SetBranchAddress("closestDPt",&closestDPt);
  t->SetBranchAddress("discr_csvSimple",&discr_csvSimple);
  t->SetBranchAddress("discr_prob",&discr_prob);
  t->SetBranchAddress("discr_ssvHighEff",&discr_ssvHighEff);
  t->SetBranchAddress("discr_ssvHighPur",&discr_ssvHighPur);
  if(isMC>1)t->SetBranchAddress("discr_cJetHighEff",&discr_cJetHighEff);
  if(isMC>1)t->SetBranchAddress("discr_cJetHighPur",&discr_cJetHighPur);
  t->SetBranchAddress("discr_tcHighEff",&discr_tcHighEff);
  t->SetBranchAddress("discr_tcHighPur",&discr_tcHighPur);

  t->SetBranchAddress("svtxdls",&svtxdls);
  t->SetBranchAddress("svtxntrk",&svtxntrk);
  t->SetBranchAddress("svtxpt",&svtxpt);
  t->SetBranchAddress("nIP",&nIP);
  t->SetBranchAddress("ipPt",&ipPt);
  t->SetBranchAddress("trackIP2dSig",&trackIP2dSig);
  t->SetBranchAddress("trackIP3dSig",&trackIP3dSig);
  t->SetBranchAddress("trackIP2d",&trackIP2d);
  t->SetBranchAddress("trackIP3d",&trackIP3d);
  t->SetBranchAddress("ipProb0",&ipProb0);
  t->SetBranchAddress("neutralMax",&neutralMax);
  t->SetBranchAddress("neutralSum",&neutralSum);
  t->SetBranchAddress("chargedMax",&chargedMax);
  t->SetBranchAddress("chargedSum",&chargedSum);
  t->SetBranchAddress("photonMax",&photonMax);
  t->SetBranchAddress("photonSum",&photonSum);
  t->SetBranchAddress("muSum",&muSum);
  t->SetBranchAddress("eSum",&eSum);

  t->SetBranchAddress("svJetDeltaR",&svJetDeltaR);

  t->SetBranchAddress("ipDist2Jet",&ipDist2Jet);
  t->SetBranchAddress("ipClosest2Jet",&ipClosest2Jet);
  t->SetBranchAddress("trackPtRel",&trackPtRel);
  t->SetBranchAddress("trackPPar",&trackPPar);
  t->SetBranchAddress("trackPParRatio",&trackPParRatio);
  t->SetBranchAddress("trackSip2dSigAboveCharm",&trackSip2dSigAboveCharm);
  t->SetBranchAddress("trackSip3dSigAboveCharm",&trackSip3dSigAboveCharm);
  t->SetBranchAddress("trackSip2dValAboveCharm",&trackSip2dValAboveCharm);
  t->SetBranchAddress("trackSip3dValAboveCharm",&trackSip3dValAboveCharm);
  t->SetBranchAddress("trackSumJetDeltaR",&trackSumJetDeltaR);
  t->SetBranchAddress("trackPtRatio",&trackPtRatio);
  t->SetBranchAddress("trackDeltaR",&trackDeltaR);
  
  t->SetBranchAddress("weight",&weight);
  t->SetBranchAddress("bin",&bin);
  t->SetBranchAddress("vz",&vz);
  if(!isMC) t->SetBranchAddress("run",&run);

  cout << " addresses set" << endl;
  int nPass[3]={0,0,0};
  int ptEntries[3] = {0,0,0};
  cout << "nEntries: " << t->GetEntries() << endl;
  for(int ievt=0; ievt<t->GetEntries()/2; ievt++){
    //if(!(ievt%100000)) cout << "entry: " << ievt << endl;
    bool flagPass[3] = {0,0,0};
    t->GetEntry(ievt);
    a_evtSelection=1;

    if(abs(vz)>15) continue;
    if(isMC && dopthat80Only && (pthat<80)) continue;
    a_vz = vz;
    if(!isMC) a_run = run;
    else a_run = 1;
    a_bin = bin;
    a_weight = weight;
    a_subid = 0;
    for(unsigned int ii=0; ii<dCandPt->size(); ii++){
      double bgweight;
      if(isMC) bgweight = eventFilter(1, 1, 2, dCandMass, dCandType, dCandChildMass,0, 0, dCandPt, dCandEta, dCandPhi, ffls3d, falpha0, fprob, dCandGenMatchPdg,corrFact);
      else bgweight = eventFilter(1, 1, 2, dCandMass, dCandType, dCandChildMass,0, 0, dCandPt, dCandEta, dCandPhi, ffls3d, falpha0, fprob, 0,corrFact);
      if(!bgweight) continue;
      a_weight*=bgweight;
      a_dCandPt = dCandPt->at(ii);
      a_dCandEta = dCandEta->at(ii);
      a_dCandPhi = dCandPhi->at(ii);
      a_dCandMass = dCandMass->at(ii);
      a_dCandType = dCandType->at(ii);
      a_dCandChildMass = dCandChildMass->at(ii);
      a_dCandCharge1 = dCandCharge1->at(ii);
      a_dCandCharge2 = dCandCharge2->at(ii);
      a_ffls3d = ffls3d->at(ii);
      a_falpha0 = falpha0->at(ii);
      a_fprob = fprob->at(ii);
      dMesons->Fill();
    }
    for(unsigned int jj=0; jj<jtpt->size(); jj++){
      if(isMC && refpt->at(jj)<20) continue;

      if(a_jtpt>40 && a_jtpt<80) ptEntries[0]++;
      else if(a_jtpt>80 && a_jtpt<150) ptEntries[1]++;
      else if(a_jtpt>150 && a_jtpt<400) ptEntries[2]++;

      double bgweight;
      if(isMC) bgweight = eventFilter(svtxmcorr->at(jj), svtxdl->at(jj), svtxpt->at(jj), dCandMass, dCandType, dCandChildMass,jteta->at(jj), jtphi->at(jj), dCandPt, dCandEta, dCandPhi, ffls3d, falpha0, fprob, dCandGenMatchPdg,corrFact);
      else bgweight = eventFilter(svtxmcorr->at(jj), svtxdl->at(jj), svtxpt->at(jj), dCandMass, dCandType, dCandChildMass,jteta->at(jj), jtphi->at(jj), dCandPt, dCandEta, dCandPhi, ffls3d, falpha0, fprob, 0,corrFact);
      if(!bgweight) continue;

      //if(discr_ssvHighEff->at(jj)<0.5) continue;

      //flagPass=true;
      a_jtpt = jtpt->at(jj);

      if(a_jtpt>40 && a_jtpt<80) flagPass[0] = true;
      else if(a_jtpt>80 && a_jtpt<150) flagPass[1] = true;
      else if(a_jtpt>150 && a_jtpt<400) flagPass[2] = true;

      a_jteta = jteta->at(jj);
      a_jtphi = jtphi->at(jj);
      a_jtm = jtm->at(jj);
      a_jtntrks = jtntrks->at(jj);
      if(isMC) a_refpt = refpt->at(jj);
      a_rawpt = rawpt->at(jj);
      a_discr_ssvHighEff = discr_ssvHighEff->at(jj);
      a_discr_ssvHighPur = discr_ssvHighPur->at(jj);
      if(isMC>1)a_discr_cJetHighEff = discr_cJetHighEff->at(jj);
      if(isMC>1)a_discr_cJetHighPur = discr_cJetHighPur->at(jj);
      a_discr_csvSimple = discr_csvSimple->at(jj);
      a_discr_tcHighEff = discr_tcHighEff->at(jj);
      a_discr_tcHighPur = discr_tcHighPur->at(jj);
      a_discr_prob = discr_prob->at(jj);
      if(svtxm->at(jj)>0) a_svtxm = svtxm->at(jj);
      else a_svtxm = -1;
      if(svtxmcorr->at(jj)>0) a_svtxmcorr = svtxmcorr->at(jj);
      else a_svtxmcorr = -1;
      if(svtxmcorr->at(jj)>10) a_svtxmcorr = 10.;
      a_svtxmEnergyFrac = TMath::Sqrt(TMath::Power(svtxm->at(jj),2)+TMath::Power(svtxpt->at(jj),2))/TMath::Sqrt(TMath::Power(jtm->at(jj),2)+TMath::Power(jtpt->at(jj),2));
      a_svtxdl = svtxdl->at(jj);
      a_svtxpt = svtxpt->at(jj);
      a_svtxptFrac = svtxpt->at(jj)/rawpt->at(jj);
      a_svtxdls = svtxdls->at(jj);
      a_svtxntrk = svtxntrk->at(jj);
      a_sv2Trkdl = sv2Trkdl->at(jj);
      a_sv2Trkdls = sv2Trkdls->at(jj);
      a_svtxTrkSumChi2 = svtxTrkSumChi2->at(jj);
      a_svtxTrkNetCharge = svtxTrkNetCharge->at(jj);
      a_svtxNtrkInCone = svtxNtrkInCone->at(jj);
      //cout << "ip2dsig size:" << trackIP2dSig->size() << endl;
      if(trackIP2dSig->size()>0){
       a_nIP = trackIP2dSig->size();
	//track-level variables are filled to each jet.  IPsigs, etc, are a collection of tracks per EVENT, not per jet.
	//first sort all the vectors
       std::sort(trackIP2dSig->begin(), trackIP2dSig->end(), std::greater<double>());
       std::sort(trackIP3dSig->begin(), trackIP3dSig->end(), std::greater<double>());
       std::sort(trackIP2d->begin(), trackIP2d->end(), std::greater<double>());
       std::sort(trackIP3d->begin(), trackIP3d->end(), std::greater<double>());
       std::sort(ipProb0->begin(), ipProb0->end(), std::greater<double>());
       std::sort(ipPt->begin(), ipPt->end(), std::greater<double>());
       std::sort(trackPtRel->begin(), trackPtRel->end(), std::greater<double>());
       std::sort(trackPPar->begin(), trackPPar->end(), std::greater<double>());
       std::sort(trackPParRatio->begin(), trackPParRatio->end(), std::greater<double>());
       std::sort(ipDist2Jet->begin(), ipDist2Jet->end(), std::greater<double>());
       std::sort(ipClosest2Jet->begin(), ipClosest2Jet->end(), std::greater<double>());
       std::sort(trackDeltaR->begin(), trackDeltaR->end(), std::less<double>());
       std::sort(trackPtRatio->begin(), trackPtRatio->end(), std::greater<double>());

       int nip2dsig=0, nip3dsig=0, nip2d=0, nip3d=0, nipProb0=0, nipPt=0;
       a_nIP = (a_nIP > 4 ? 4:a_nIP);
       for(int iip=0; iip<a_nIP; iip++){
         a_trackIP2dSig[iip] = trackIP2dSig->at(iip);
         a_trackIP3dSig[iip] = trackIP3dSig->at(iip);
         a_trackIP2d[iip] = trackIP2d->at(iip);
         a_trackIP3d[iip] = trackIP3d->at(iip);
         a_ipProb0[iip] = ipProb0->at(iip);
         a_ipPt[iip] = ipPt->at(iip);
         a_trackPtRel[iip] = trackPtRel->at(iip);
         a_trackPPar[iip] = trackPPar->at(iip);
         a_trackPParRatio[iip] = trackPParRatio->at(iip);
         a_ipDist2Jet[iip] = ipDist2Jet->at(iip);
         a_ipClosest2Jet[iip] = ipClosest2Jet->at(iip);
         a_trackDeltaR[iip] = trackDeltaR->at(iip);
         a_trackPtRatio[iip] = trackPtRatio->at(iip);
       }

       if(a_jtpt>40 && a_jtpt<80){
         if(flagPass){
          nPass[0]++;
        }
      } 
      else if(a_jtpt>80 && a_jtpt<150){
        if(flagPass){
         nPass[1]++;
       }
     }
     else if(a_jtpt>150 && a_jtpt<400){
       if(flagPass){
        nPass[2]++;
      }
    }
  }
   a_trackSip2dSigAboveCharm = trackSip2dSigAboveCharm->at(jj);
     a_trackSip3dSigAboveCharm = trackSip3dSigAboveCharm->at(jj);
     a_trackSip2dValAboveCharm = trackSip2dValAboveCharm->at(jj);
     a_trackSip3dValAboveCharm = trackSip3dValAboveCharm->at(jj);
     a_svJetDeltaR = svJetDeltaR->at(jj);
     a_trackSumJetDeltaR = trackSumJetDeltaR->at(jj);

     a_neutralMax = neutralMax->at(jj);
     a_neutralSum = neutralSum->at(jj);
     a_chargedMax = chargedMax->at(jj);
     a_chargedSum = chargedSum->at(jj);
     a_photonMax = photonMax->at(jj);
     a_photonSum = photonSum->at(jj);
     a_muSum = muSum->at(jj);
     a_eSum = eSum->at(jj);
     a_closestDPt = closestDPt->at(jj);
     a_closestDMass = closestDMass->at(jj);
     a_closestDChildMass = closestDChildMass->at(jj);
     a_closestDType = closestDType->at(jj);
     if(a_closestDPt==0) a_closestDPt=-999;
     if(isMC){
       a_refparton_flavorForB = refparton_flavorForB->at(jj);
       a_genMatch = genMatch->at(jj);
     }
     a_djetR = djetR->at(jj);

      //a_djetR = findDR(isMC, dCandPt, dCandEta, dCandPhi, dCandMass, dCandChildMass, dCandCharge1, dCandCharge2, dCandType, dCandGenMatchPdg, jteta->at(jj), jtphi->at(jj), a_closestDMass, a_closestDChildMass, a_closestDType, a_closestDPdg);

      //double arrBDTvals[] = { a_djetR, static_cast<double>(a_nIP), a_svtxm, a_svtxdl, a_svtxdls, static_cast<double>(a_svtxntrk), a_jteta, a_ipProb0[0], a_ipPt[0], a_trackIP2dSig[0], a_trackIP3dSig[0], a_trackIP2d[0], a_trackIP3d[0], a_trackPtRel[0], a_trackPPar[0], a_trackPParRatio[0], a_ipDist2Jet[0], a_ipClosest2Jet[0], a_trackDeltaR[0], a_trackPtRatio[0], a_ipProb0[1], a_ipPt[1], a_trackIP2dSig[1], a_trackIP3dSig[1], a_trackIP2d[1], a_trackIP3d[1], a_trackPtRel[1], a_trackPPar[1], a_trackPParRatio[1], a_ipDist2Jet[1], a_ipClosest2Jet[1], a_trackDeltaR[1], a_trackPtRatio[1], a_ipProb0[2], a_ipPt[2], a_trackIP2dSig[2], a_trackIP3dSig[2], a_trackIP2d[2], a_trackIP3d[2], a_trackPtRel[2], a_trackPPar[2], a_trackPParRatio[2], a_ipDist2Jet[2], a_ipClosest2Jet[2], a_trackDeltaR[2], a_trackPtRatio[2], a_ipProb0[3], a_ipPt[3], a_trackIP2dSig[3], a_trackIP3dSig[3], a_trackIP2d[3], a_trackIP3d[3], a_trackPtRel[3], a_trackPPar[3], a_trackPParRatio[3], a_ipDist2Jet[3], a_ipClosest2Jet[3], a_trackDeltaR[3], a_trackPtRatio[3], a_trackSip2dValAboveCharm, a_trackSip3dValAboveCharm, a_svJetDeltaR, a_trackSumJetDeltaR };

     double arrBDTvals[] = { a_svtxptFrac, static_cast<double>(a_nIP), a_svtxm, a_svtxmEnergyFrac, a_svtxpt, a_svtxmcorr, a_svtxdl, a_svtxdls, static_cast<double>(a_svtxntrk), a_svtxTrkSumChi2, a_svtxTrkNetCharge, a_svtxNtrkInCone,  a_closestDMass, a_closestDType, a_closestDPt, a_trackIP2dSig[0], a_trackIP3dSig[0], a_trackIP2d[0], a_trackIP3d[0], a_trackIP2dSig[1], a_trackIP3dSig[1], a_trackIP2d[1], a_trackIP3d[1], a_trackIP2dSig[2], a_trackIP3dSig[2], a_trackIP2d[2], a_trackIP3d[2], a_ipDist2Jet[0], a_ipClosest2Jet[0], a_svJetDeltaR, a_trackSumJetDeltaR};

     //double arrBDTvals2[] = { a_trackIP2dSig[0], a_trackIP3dSig[0], a_trackIP2d[0], a_trackIP3d[0], a_trackIP2dSig[1], a_trackIP3dSig[1], a_trackIP2d[1], a_trackIP3d[1], a_trackIP2dSig[2], a_trackIP3dSig[2], a_trackIP2d[2], a_trackIP3d[2], a_ipDist2Jet[0], a_ipClosest2Jet[0], a_svJetDeltaR, a_trackSumJetDeltaR };

     vector<double> tempBDTvals;
     for(int i=0; i<31; i++) tempBDTvals.push_back(arrBDTvals[i]);
      //vector<double> tempBDTvals2;
    //for(int i=0; i<16; i++) tempBDTvals2.push_back(arrBDTvals2[i]);

      //vector<double> tempBDTvals {a_djetR, a_discr_ssvHighEff, a_discr_csvSimple, a_discr_prob, a_closestDPt, a_svtxm, a_svtxdl};
      if(refillTMVAvars){
        a_bdt = testBDTG_CvL_withSvtx.GetMvaValue(tempBDTvals);
        a_bdt_bvc = testBDTG_BvC_withSvtx.GetMvaValue(tempBDTvals);
        //a_bdt_noVtx = testBDTG_CvL_noSvtx.GetMvaValue(tempBDTvals2);
        //a_bdt_bvc_noVtx = testBDTG_BvC_noSvtx.GetMvaValue(tempBDTvals2);
      }
      //if(a_svtxm>0) jets->Fill();
      //else jetsNoSvtx->Fill();
      jets->Fill();

      //clear out stuff
      for(int iip=0; iip<4; iip++){
       a_trackIP2dSig[iip] = 0;
       a_trackIP3dSig[iip] = 0;
       a_trackIP2d[iip] = 0;
       a_trackIP3d[iip] = 0;
       a_ipProb0[iip] = 0;
       a_ipPt[iip] = 0;
       a_trackPtRel[iip] = 0;
       a_trackPPar[iip] = 0;
       a_trackPParRatio[iip] = 0;
       a_ipDist2Jet[iip] = 0;
       a_ipClosest2Jet[iip] = 0;
       a_trackDeltaR[iip] = 0;
       a_trackPtRatio[iip] = 0;
     }
   }
 }

 for(int i=0; i<3; i++){
   cout << "Cut Efficiency ptBin " << i << ": " << nPass[i] << " / " << ptEntries[i] << " ( " << 100*(double)nPass[i]/(double)ptEntries[i]<< "% )" << endl;
 }

 fout->cd();
 dMesons->Write();
 jets->Write();
 //jetsNoSvtx->Write();
 fout->Close();
}
