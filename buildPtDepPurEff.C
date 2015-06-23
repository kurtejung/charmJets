#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TMath.h"
#include <vector>
#include <iostream>

#ifdef __CINT__
#pragma link C++ class vector<vector<double> >;
#else
template class std::vector<std::vector<double> >;
#endif

using std::vector;
using std::endl;
using std::cout;

	
	double eventFilter(double svtxmcorr, double svtxdl, double svtxpt, vector<double> *dCandMass, vector<double> *dCandType, 
  vector<double> *dCandChildMass, double jteta, double jtphi, vector<double> *dCandPt, vector<double> *dCandEta, 
  vector<double> *dCandPhi, vector<double> *ffls3d, vector<double> *falpha0, vector<double> *fprob0, 
  vector<int> *dCandGenMatchPdg, bool doDR){

  //cout << "isMC: " << isMC << endl;

  int isMC=1;

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
        if(!doDR || sqrt(pow(dCandEta->at(icand)-jteta,2)+pow(dCandPhi->at(icand)-jtphi,2))<0.3){
          checkDMeson=true;
          if(isMC){
            if(abs(dCandGenMatchPdg->at(icand))==413) hasRealD=true;
            if(dCandPt->at(icand)>maxMesonPt){ 
              maxMesonPt = dCandPt->at(icand);
              //weightToRet = bgAdjuster->GetBinContent(bgAdjuster->FindBin(dCandPt->at(icand)));
            }
          }
        }
      }
      else if(dCandMass->at(icand)>d0MassLow && dCandType->at(icand)==2 && dCandMass->at(icand)<d0MassHigh){
        if(!doDR || sqrt(pow(dCandEta->at(icand)-jteta,2)+pow(dCandPhi->at(icand)-jtphi,2))<0.3){
          checkDMeson=true;
          if(isMC){
            if(abs(dCandGenMatchPdg->at(icand))==421) hasRealD=true;
            if(dCandPt->at(icand)>maxMesonPt){ 
              maxMesonPt = dCandPt->at(icand);
              //weightToRet = bgAdjuster->GetBinContent(bgAdjuster->FindBin(dCandPt->at(icand)));
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

void buildPtDepPurEff(){

	TFile *fout = new TFile("dRandCandPurEff_light.root","recreate");

	TH1D *dRPurity = new TH1D("dRPurity","",10,40,400);
	TH1D *dREff = new TH1D("dREff","",10,40,400);

	TH1D *dCandFiltPur = new TH1D("dCandFiltPur","",10,40,400);
	TH1D *dCandFiltEff = new TH1D("dCandFiltEff","",10,40,400);

	TFile *fin = new TFile("input/DMesonCJet_QCDJetOnly_Merged_addLHCbVars_pPbMC_centReweight_ppReco_akPu3PF.root");
	TTree *t = (TTree*)fin->Get("ct");

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
  t->SetBranchAddress("dCandGenMatchPdg",&dCandGenMatchPdg);

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

    t->SetBranchAddress("refpt",&refpt);
    t->SetBranchAddress("refparton_flavorForB",&refparton_flavorForB);
    t->SetBranchAddress("genMatch",&genMatch);
    t->SetBranchAddress("subid",&subid);
    
  t->SetBranchAddress("closestDMass",&closestDMass);
  t->SetBranchAddress("closestDChildMass",&closestDChildMass);
  t->SetBranchAddress("closestDType",&closestDType);
  t->SetBranchAddress("closestDPt",&closestDPt);
  t->SetBranchAddress("discr_csvSimple",&discr_csvSimple);
  t->SetBranchAddress("discr_prob",&discr_prob);
  t->SetBranchAddress("discr_ssvHighEff",&discr_ssvHighEff);
  t->SetBranchAddress("discr_ssvHighPur",&discr_ssvHighPur);
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

  cout << " addresses set" << endl;
  int totalJetPassDR[10] = {0,0,0,0,0,0,0,0,0,0};
  int totalcJet[10] = {0,0,0,0,0,0,0,0,0,0};
  int cJetPassDR[10] = {0,0,0,0,0,0,0,0,0,0};
  int totalJetPassCand[10] = {0,0,0,0,0,0,0,0,0,0};
  int cJetPassCand[10] = {0,0,0,0,0,0,0,0,0,0};
  cout << "nEntries: " << t->GetEntries() << endl;
  for(int ievt=0; ievt<t->GetEntries(); ievt++){
    t->GetEntry(ievt);

  	for(unsigned int jj=0; jj<jtpt->size(); jj++){

  		if(jtpt->at(jj)<40) continue;
  		int jtptBin = (jtpt->at(jj)-40)/36;

  		int dRFilter = eventFilter(svtxmcorr->at(jj), svtxdl->at(jj), svtxpt->at(jj), dCandMass, dCandType, dCandChildMass,jteta->at(jj), jtphi->at(jj), dCandPt, dCandEta, dCandPhi, ffls3d, falpha0, fprob, dCandGenMatchPdg, 1);
  		int candFilter = eventFilter(svtxmcorr->at(jj), svtxdl->at(jj), svtxpt->at(jj), dCandMass, dCandType, dCandChildMass,jteta->at(jj), jtphi->at(jj), dCandPt, dCandEta, dCandPhi, ffls3d, falpha0, fprob, dCandGenMatchPdg, 0);

  		if(abs(refparton_flavorForB->at(jj))!=4 && abs(refparton_flavorForB->at(jj))!=5) {
  			if(dRFilter) cJetPassDR[jtptBin]++;
  			if(candFilter) cJetPassCand[jtptBin]++;
  			totalcJet[jtptBin]++;
  		}
     if(dRFilter) totalJetPassDR[jtptBin]++;
     if(candFilter) totalJetPassCand[jtptBin]++;
  	}
  }

  for(int i=0; i<10; i++){
  	dRPurity->SetBinContent(i, (double)cJetPassDR[i]/(double)totalJetPassDR[i]);
  	dREff->SetBinContent(i, (double)cJetPassDR[i]/(double)totalcJet[i]);

  	dCandFiltPur->SetBinContent(i, (double)cJetPassCand[i]/(double)totalJetPassCand[i]);
  	dCandFiltEff->SetBinContent(i, (double)cJetPassCand[i]/(double)totalcJet[i]);
  }

  fout->cd();
  dRPurity->SetXTitle("Jet p_{T}");
  dRPurity->SetYTitle("Light-Jet Purity");
  dRPurity->Write();
  dREff->SetXTitle("Jet p_{T}");
  dREff->SetYTitle("Light-Jet Efficiency");
  dREff->Write();
  dCandFiltEff->SetMarkerColor(2);
  dCandFiltEff->SetLineColor(2);
  dCandFiltEff->Write();
  dCandFiltPur->SetMarkerColor(2);
  dCandFiltPur->SetLineColor(2);
  dCandFiltPur->Write();
  fout->Close();

}