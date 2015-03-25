#include "TTree.h"
#include <vector>
#include "TFile.h"

//to get the trained discriminator information from TMVA:
#include "/Users/kjung/charmJets/pPb/tmvaTesting/weights/TMVAClassification_BDT.class.C"
#include "/Users/kjung/charmJets/pPb/tmvaTesting/weights/TMVAClassification_RuleFit.class.C"
//end

#ifdef __CINT__
#pragma link C++ class vector<vector<double> >;
#else
template class std::vector<std::vector<double> >;
#endif

double findDR(int isMC, vector<double> *dCandPt, vector<double> *dCandEta, vector<double> *dCandPhi, vector<double> *dCandMass, vector<double> *dCandChildMass, vector<double> *dCandCharge1, vector<double> *dCandCharge2, vector<double> *dCandType, vector<int> *dCandGenMatchedPdg, double jteta, double jtphi, double &dMass, double &dChildMass, double &dType, int &dPDG){

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

void convertVectorToDoubleTree(int isMC=1){

  bool refillTMVAvars = false;

  TFile *fout;
  if(isMC==1) fout = new TFile("../input/DMesonCJet_QCDJetOnly_pPbMC_ppReco_akPu3PF_convertToJetTree.root","RECREATE");
  else if(isMC==2) fout = new TFile("../input/DMesonCJet_CJetOnly_pPbMC_ppReco_akPu3PF_convertToJetTree_temp.root","RECREATE");
  else fout = new TFile("../input/DMesonCJet_pPbData_ppReco_akPu3PF_convertToJetTree.root","RECREATE");

  TTree *jets = new TTree("jets","");
  TTree *dMesons = new TTree("dMesons","");

  double a_dCandPt, a_dCandEta, a_dCandPhi, a_dCandMass, a_dCandType, a_dCandChildMass, a_dCandCharge1, a_dCandCharge2;
  double a_jtpt, a_jteta, a_jtphi, a_refpt, a_rawpt, a_djetR, a_jtm;
  double a_refparton_flavorForB, a_svtxm, a_svtxdl, a_svtxpt, a_svtxdls;
  double a_closestDMass, a_closestDChildMass, a_closestDType;
  int a_closestDPdg;
  int a_genMatch, a_bin, a_run, a_subid, a_evtSelection;
  double a_weight, a_vz, a_bdt, a_ruleFit;
  double a_discr_csvSimple, a_discr_ssvHighEff, a_discr_ssvHighPur, a_discr_cJetHighEff, a_discr_cJetHighPur, a_closestDPt, a_discr_prob, a_discr_tcHighEff, a_discr_tcHighPur;

  int a_nIP, a_svtxntrk, a_jtntrks;
  double a_neutralMax, a_neutralSum, a_chargedMax, a_chargedSum, a_photonMax, a_photonSum, a_muSum, a_eSum; 
  double a_trackIP2dSig, a_trackIP3dSig, a_trackIP2d, a_trackIP3d, a_ipProb0, a_ipPt;
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
  
  jets->Branch("jtpt",&a_jtpt);
  jets->Branch("jteta",&a_jteta);
  jets->Branch("jtphi",&a_jtphi);
  jets->Branch("rawpt",&a_rawpt);
  jets->Branch("jtm",&a_jtm);
  jets->Branch("jtntrks",&a_jtntrks);
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
  jets->Branch("svtxdl",&a_svtxdl);
  jets->Branch("svtxdls",&a_svtxdls);
  jets->Branch("svtxntrk",&a_svtxntrk);
  jets->Branch("svtxpt",&a_svtxpt);
  jets->Branch("nIP",&a_nIP);
  jets->Branch("ipPt",&a_ipPt);
  jets->Branch("trackIP2dSig",&a_trackIP2dSig);
  jets->Branch("trackIP3dSig",&a_trackIP3dSig);
  jets->Branch("trackIP2d",&a_trackIP2d);
  jets->Branch("trackIP3d",&a_trackIP3d);
  jets->Branch("ipProb0",&a_ipProb0);
  jets->Branch("neutralMax",&a_neutralMax);
  jets->Branch("neutralSum",&a_neutralSum);
  jets->Branch("chargedMax",&a_chargedMax);
  jets->Branch("chargedSum",&a_chargedSum);
  jets->Branch("photonMax",&a_photonMax);
  jets->Branch("photonSum",&a_photonSum);
  jets->Branch("muSum",&a_muSum);
  jets->Branch("eSum",&a_eSum);
  
  if(refillTMVAvars) jets->Branch("bdt",&a_bdt);
  if(refillTMVAvars) jets->Branch("ruleFit",&a_ruleFit);

  TFile *fin;
  if(isMC==1) fin = new TFile("../input/DMesonCJet_QCDJetOnly_pPbMC_centReweight_IVFVtx_ppReco_akPu3PF_pthat80Only_458.root");
  else if(isMC==2) fin = new TFile("../input/DMesonCJet_CJetOnly_pPbMC_ppReco_akPu3PF_1190.root");
  else fin = new TFile("../input/DMesonCJet_JetTrg_pPbdata_Full_ppReco_akPu3PF_2-19.root");
  TTree *t = (TTree*)fin->Get("ct");
  
  double weight;
  float vz;
  int bin, run, subid;
  vector<double> *dCandPt=0, *dCandMass=0, *dCandEta=0, *dCandPhi=0, *jtpt=0, *jteta=0, *jtphi=0, *refpt=0, *rawpt=0, *djetR=0, *jtm=0;
  vector<double> *dCandType=0, *dCandChildMass=0, *dCandCharge1=0, *dCandCharge2=0;
  vector<double> *closestDMass=0, *closestDType=0, *closestDChildMass=0;
  vector<int> *refparton_flavorForB=0;
  vector<int> *genMatch=0, *dCandGenMatchPdg=0, *jtntrks=0;
  vector<double> *discr_csvSimple=0, *discr_ssvHighEff=0, *discr_ssvHighPur=0, *discr_cJetHighEff=0, *discr_cJetHighPur=0, *discr_tcHighEff=0, *discr_tcHighPur=0, *closestDPt=0, *discr_prob=0, *svtxm=0, *svtxdl=0;
  vector<bool> *isGSP=0;
  vector<double> *svtxdls=0, *svtxpt=0, *neutralMax=0, *neutralSum=0, *chargedMax=0, *chargedSum=0, *photonMax=0, *photonSum=0, *muSum=0, *eSum=0;
  vector<int> *nIP=0, *svtxntrk=0;
  vector<vector<double> > *trackIP2dSig=0, *trackIP3dSig=0, *trackIP2d=0, *trackIP3d=0, *ipProb0=0, *ipPt=0;

  //if(!isMC) dCandGenMatchPdg = new vector<int>;

  vector<string> inputVars { "djetR", "discr_ssvHighEff", "discr_csvSimple", "discr_prob", "closestDPt", "svtxm", "svtxdl" };
  ReadBDT testBDT(inputVars);
  ReadRuleFit testRuleFit(inputVars);
  
  t->SetBranchAddress("dCandPt",&dCandPt);
  t->SetBranchAddress("dCandEta",&dCandEta);
  t->SetBranchAddress("dCandPhi",&dCandPhi);
  t->SetBranchAddress("dCandMass",&dCandMass);
  t->SetBranchAddress("dCandType",&dCandType);
  t->SetBranchAddress("dCandChildMass",&dCandChildMass);
  t->SetBranchAddress("dCandCharge1",&dCandCharge1);
  t->SetBranchAddress("dCandCharge2",&dCandCharge2);
  if(isMC) t->SetBranchAddress("dCandGenMatchPdg",&dCandGenMatchPdg);

  t->SetBranchAddress("jtpt",&jtpt);
  t->SetBranchAddress("jteta",&jteta);
  t->SetBranchAddress("jtphi",&jtphi);
  t->SetBranchAddress("jtm",&jtm);
  t->SetBranchAddress("jtntrks",&jtntrks);
  t->SetBranchAddress("rawpt",&rawpt);
  t->SetBranchAddress("djetR",&djetR);
  t->SetBranchAddress("svtxm",&svtxm);
  t->SetBranchAddress("svtxdl",&svtxdl);
  t->SetBranchAddress("svtxntrk",&svtxntrk);
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
  
  t->SetBranchAddress("weight",&weight);
  t->SetBranchAddress("bin",&bin);
  t->SetBranchAddress("vz",&vz);
  if(!isMC) t->SetBranchAddress("run",&run);

  for(int ievt=0; ievt<t->GetEntries(); ievt++){
    t->GetEntry(ievt);
    a_evtSelection=1;
    if(abs(vz)>15) continue;
    a_vz = vz;
    if(!isMC) a_run = run;
    else a_run = 1;
    a_bin = bin;
    a_weight = weight;
    a_subid = 0;
    for(unsigned int ii=0; ii<dCandPt->size(); ii++){
      a_dCandPt = dCandPt->at(ii);
      a_dCandEta = dCandEta->at(ii);
      a_dCandPhi = dCandPhi->at(ii);
      a_dCandMass = dCandMass->at(ii);
      a_dCandType = dCandType->at(ii);
      a_dCandChildMass = dCandChildMass->at(ii);
      a_dCandCharge1 = dCandCharge1->at(ii);
      a_dCandCharge2 = dCandCharge2->at(ii);
      dMesons->Fill();
    }
    for(unsigned int jj=0; jj<jtpt->size(); jj++){
      if(isMC && refpt->at(jj)<20) continue;
      a_jtpt = jtpt->at(jj);
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
      a_svtxdl = svtxdl->at(jj);
      a_svtxpt = svtxpt->at(jj);
      a_svtxdls = svtxdls->at(jj);
      a_svtxntrk = svtxntrk->at(jj);
      //cout << "ip2dsig size:" << trackIP2dSig->size() << endl;
      if(trackIP2dSig->size()>0){
	a_nIP = trackIP2dSig->at(0).size();
	//track-level variables are filled to each jet.  IPsigs, etc, are a collection of tracks per EVENT, not per jet.
	int nip2dsig=0, nip3dsig=0, nip2d=0, nip3d=0, nipProb0=0, nipPt=0;
	for(int iip=0; iip<a_nIP; iip++){
	  if(trackIP2dSig->at(0).at(iip)>0){ a_trackIP2dSig += trackIP2dSig->at(0).at(iip); nip2dsig++; }
	  if(trackIP3dSig->at(0).at(iip)>0){ a_trackIP3dSig += trackIP3dSig->at(0).at(iip); nip3dsig++; }
	  if(trackIP2d->at(0).at(iip)>0){ a_trackIP2d += trackIP2d->at(0).at(iip); nip2d++; }
	  if(trackIP3d->at(0).at(iip)>0){ a_trackIP3d += trackIP3d->at(0).at(iip); nip3d++; }
	  if(ipProb0->at(0).at(iip)>0){ a_ipProb0 += ipProb0->at(0).at(iip); nipProb0++; }
	  if(ipPt->at(0).at(iip)>0){ a_ipPt += ipPt->at(0).at(iip); nipPt++; }
	}
	if(nip2dsig) a_trackIP2dSig/=(double)nip2dsig;
	if(nip3dsig) a_trackIP3dSig/=(double)nip3dsig;
	if(nip2d) a_trackIP2d/=(double)nip2d;
	if(nip3d) a_trackIP3d/=(double)nip3d;
	if(nipProb0) a_ipProb0/=(double)nipProb0;
	if(nipPt) a_ipPt/=(double)nipPt;
      }
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
      //a_djetR = djetR->at(jj);
      
      a_djetR = findDR(isMC, dCandPt, dCandEta, dCandPhi, dCandMass, dCandChildMass, dCandCharge1, dCandCharge2, dCandType, dCandGenMatchPdg, jteta->at(jj), jtphi->at(jj), a_closestDMass, a_closestDChildMass, a_closestDType, a_closestDPdg);
      
      vector<double> tempBDTvals {a_djetR, a_discr_ssvHighEff, a_discr_csvSimple, a_discr_prob, a_closestDPt, a_svtxm, a_svtxdl};
      if(refillTMVAvars) a_bdt = testBDT.GetMvaValue(tempBDTvals);
      if(refillTMVAvars) a_ruleFit = testRuleFit.GetMvaValue(tempBDTvals);
      
      jets->Fill();
      //clear out stuff
      a_trackIP2dSig=0;
      a_trackIP3dSig=0;
      a_trackIP2d=0;
      a_trackIP3d=0;
      a_ipProb0=0;
      a_ipPt=0;
    }
  }
  fout->cd();
  dMesons->Write();
  jets->Write();
  fout->Close();
}
