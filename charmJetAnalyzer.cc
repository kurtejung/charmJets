#include "TTree.h"
#include "TH1D.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>

#include "UserCode/OpenHF/interface/hfcand_v1.hh"
#include "UserCode/OpenHF/interface/snglhfcand_v1.hh"

using namespace std;

// ** constants ******
const int MAXJETS = 100;
const int MAXPARTICLES = 5000;
//******************

snglhfcand* convertGenToHFCand(float pt, float eta, float phi, float pdg, int gindex1, int gindex2, int gindex_dau2){

  snglhfcand* cc = new snglhfcand_v1();
  cc->set_fpt(pt);
  cc->set_feta(eta);
  cc->set_fphi(phi);
  //push back all permutations of daughter particles since it doesn't matter how they're organized for MC truth
  vector<int> vgindex1, vgindex2, vgindex_dau2;
  vgindex1.push_back(gindex1);
  vgindex1.push_back(gindex2);
  if(gindex_dau2>0) vgindex1.push_back(gindex_dau2);
  vgindex2.push_back(gindex2);
  vgindex2.push_back(gindex1);
  if(gindex_dau2>0) vgindex2.push_back(gindex_dau2);  
  if(gindex_dau2>0) vgindex_dau2.push_back(gindex_dau2);
  if(gindex_dau2>0) vgindex_dau2.push_back(gindex1);
  if(gindex_dau2>0) vgindex_dau2.push_back(gindex2);
  cc->set_gIndex1(vgindex1);
  cc->set_gIndex2(vgindex2);
  cc->set_gIndex_dau2(vgindex_dau2);
  if(abs(pdg)==413) cc->set_type(1);
  if(abs(pdg)==421) cc->set_type(2);
  if(abs(pdg)==431) cc->set_type(3);
  if(abs(pdg)==411) cc->set_type(5);
  return cc;
}

int findGenDstrange(snglhfcand* cand, int* pdg, int* nMothers, int (*motherIdx)[200], int* nDaughters, bool checkDoubleCounting, bool retIdx){

  if(cand->get_type() != 3 && cand->get_type() != 4 ){
    cout << "how did you get here??" << endl;
    return -999;
  }
  
  vector<int>  gIndex_dau2 = cand->get_gIndex_dau2();
  vector<int>  gIndex1 = cand->get_gIndex1();
  vector<int>  gIndex2 = cand->get_gIndex2();
  
  
  for(unsigned int n1 = 0; n1 < gIndex1.size(); n1++ )
    {
      for(unsigned int n2 = 0; n2 < gIndex2.size(); n2++ )
	{
	  for (unsigned int n3 = 0; n3 < gIndex_dau2.size(); n3++ )
	    {
	      
	      int dau1 = gIndex1[n1];
	      int dau2 = gIndex2[n2];
	      int dau3 = gIndex_dau2[n3];
	      
	      if ( dau3 == -999 || dau1 == -999 || dau2 == -999 )     continue;
	      
	      //mathed to three different gen particles
	      if ( dau1 == dau2 || dau1 == dau3 || dau2 == dau3)     continue;
	      	     
	      //only one mother
	      if (!( nMothers[dau1] == 1 && nMothers[dau2] == 1 && nMothers[dau3] == 1 ))   continue;
	      
	      //mother index the same and not -999
	      if ( (motherIdx[dau1][0] != motherIdx[dau2][0]) || (motherIdx[dau1][0] == -999) ) continue;
	      if(cand->get_type() == 3){
		//Ds->pi phi -> pi K K
		if(!checkDoubleCounting && !( ( abs(pdg[dau1]) == 321 && abs(pdg[dau2]) == 321 ) ) )
		  continue;
		if ( abs(pdg[dau3]) != 211 )   continue;
		int phiIdx = motherIdx[dau1][0];
		//cout << "checkpoint1 pdg: " << pdg[phiIdx] << " ndaught: "<< nDaughters[phiIdx] << endl;
		if( ! (abs (pdg[phiIdx]) == 333)) continue;
		if(nDaughters[phiIdx] != 2)  continue;
		if( nMothers[phiIdx] != 1 || motherIdx[phiIdx][0] != motherIdx[dau3][0] || motherIdx[dau3][0] == -999 )    continue;
		
		int allmotherindex = motherIdx[phiIdx][0];
	       
		if( ! (abs(pdg[allmotherindex]) == 431)) continue;
		if( nDaughters[allmotherindex] != 2)  continue;
		//cout << "found Ds->PiPhi!" << endl;
		return pdg[allmotherindex];
	      }
	      if(cand->get_type() == 4){
		//Ds -> K*K -> pi K K
		if(!checkDoubleCounting && !( ( abs(pdg[dau1]) == 211 && abs(pdg[dau2]) == 321 ) ) )
		  continue;
		if(checkDoubleCounting &&  !( ( abs(pdg[dau1]) == 321 && abs(pdg[dau2]) == 211 ) ) ) continue;
		if ( abs(pdg[dau3]) != 321 )   continue;
		int kstarIdx = motherIdx[dau1][0];
		if( ! (abs (pdg[kstarIdx]) == 313)) continue;
		if(nDaughters[kstarIdx] != 2)  continue;
		if( nMothers[kstarIdx] != 1 || motherIdx[kstarIdx][0] != motherIdx[dau3][0] || motherIdx[dau3][0] == -999 )    continue;
		
		int allmotherindex = motherIdx[kstarIdx][0];
	       
		if( ! (abs(pdg[allmotherindex]) == 431)) continue;
		if( nDaughters[allmotherindex] != 2)  continue;
		if(!retIdx){
		  return pdg[allmotherindex];
		}
		else{
		  return allmotherindex;
		}
	      }
	    }
	}
    }
  return -999;
}

int findGenDstar(snglhfcand* cand, int* pdg, int* nMothers, int (*motherIdx)[200], int* nDaughters, bool checkDoubleCounting, bool retIdx){

  if(cand->get_type() != 1){
    int retval = findGenDstrange(cand,pdg,nMothers,motherIdx,nDaughters,checkDoubleCounting,retIdx);
    return retval;
  }

  vector<int>  gIndex_dau2 = cand->get_gIndex_dau2();
  vector<int>  gIndex1 = cand->get_gIndex1();
  vector<int>  gIndex2 = cand->get_gIndex2();
  
  
  for(unsigned int n1 = 0; n1 < gIndex1.size(); n1++ )
    {
      for(unsigned int n2 = 0; n2 < gIndex2.size(); n2++ )
	{
	  for (unsigned int n3 = 0; n3 < gIndex_dau2.size(); n3++ )
	    {
	      
	      int dau1 = gIndex1[n1];
	      int dau2 = gIndex2[n2];
	      int dau3 = gIndex_dau2[n3];
	      
	      if ( dau3 == -999 || dau1 == -999 || dau2 == -999 )     continue;
	      
	      //mathed to three different gen particles
	      if ( dau1 == dau2 || dau1 == dau3 || dau2 == dau3)     continue;
	      	      
	      // one kaon, two pion
	      if(!checkDoubleCounting && !( ( abs(pdg[dau1]) == 211 && abs(pdg[dau2]) == 321 ) ) ) continue;
              if(checkDoubleCounting && !( ( abs(pdg[dau1]) == 321 && abs(pdg[dau2]) == 211 ) ) ) continue;
              if ( abs(pdg[dau3]) != 211 )   continue;
	      
	      //only one mother
	      if (!( nMothers[dau1] == 1 && nMothers[dau2] == 1 && nMothers[dau3] == 1 ))   continue;

	      //mother index the same and not -999
	      if ( (motherIdx[dau1][0] != motherIdx[dau2][0]) || (motherIdx[dau1][0] == -999) ) continue;
	      int d0index = motherIdx[dau1][0];
	      //is D charged and just three daughters
	      if( ! (abs (pdg[d0index]) == 421)) continue;
	      if(nDaughters[d0index] != 2)  continue;
	      if( nMothers[d0index] != 1 || motherIdx[d0index][0] != motherIdx[dau3][0] || motherIdx[dau3][0] == -999 )    continue;
	      
	      int allmotherindex = motherIdx[d0index][0];
	      if(! (abs(pdg[allmotherindex]) == 413)) continue;
	      if( nDaughters[allmotherindex] != 2)  continue;
	      if(!retIdx){
		return pdg[allmotherindex];
	      }
	      else{
		return allmotherindex;
	      }	      
	    }
	}
    }
  return -999;
}

int findGenDcharged(snglhfcand* cand, int* pdg, int* nMothers, int (*motherIdx)[200], int* nDaughters, bool checkDoubleCounting, bool retIdx){

  if(cand->get_type() != 5){
    int returnVal = findGenDstar(cand, pdg, nMothers, motherIdx, nDaughters, checkDoubleCounting,retIdx);
    // if(returnVal != -999) cout << "returning 2 " << returnVal << endl;
    return returnVal;
  }

  vector<int>  gIndex_dau2 = cand->get_gIndex_dau2();
  vector<int>  gIndex1 = cand->get_gIndex1();
  vector<int>  gIndex2 = cand->get_gIndex2();
  
  for(unsigned int n1 = 0; n1 < gIndex1.size(); n1++ )
    {
      for(unsigned int n2 = 0; n2 < gIndex2.size(); n2++ )
	{
	  for (unsigned int n3 = 0; n3 < gIndex_dau2.size(); n3++ )
	    {
	      
	      int dau1 = gIndex1[n1];
	      int dau2 = gIndex2[n2];
	      int dau3 = gIndex_dau2[n3];
	      
	      if ( dau3 == -999 || dau1 == -999 || dau2 == -999 )     continue;

	      //mathed to three different gen particles
	      if ( dau1 == dau2 || dau1 == dau3 || dau2 == dau3 )     continue;

	      // one kaon, two pion
	      if(!checkDoubleCounting && !(( abs(pdg[dau1]) == 211 && abs(pdg[dau2]) == 321 && abs(pdg[dau3]) == 211 )) )  continue;
	      if(checkDoubleCounting && !(( abs(pdg[dau1]) == 321 && abs(pdg[dau2]) == 211 && abs(pdg[dau3]) == 211 )) )  continue;
	      
	      //only one mother
	      if (!( nMothers[dau1] == 1 && nMothers[dau2] == 1 && nMothers[dau3] == 1 ))   continue;
	      
	      //mother index the same and not -999
	      if ( (motherIdx[dau1][0] != motherIdx[dau2][0]) || (motherIdx[dau1][0] != motherIdx[dau3][0]) || (motherIdx[dau1][0] == -999) ) continue;
	      int allmotherindex = motherIdx[dau1][0];
	      //is D charged and just three daughters
	      if(! (abs (pdg[allmotherindex]) == 411)) continue;
	      if(nDaughters[allmotherindex] != 3)  continue;
	      
	      if(!retIdx){
		return pdg[allmotherindex];
	      }
	      else{
		return allmotherindex;
	      }
	    }
	}
    }
  return -999;
}

int findGenD0(snglhfcand* cand, int* pdg, int* nMothers, int (*motherIdx)[200], int* nDaughters, bool checkDoubleCounting, bool retIdx)
{
  if(cand->get_type() != 2){
    int returnVal = findGenDcharged(cand, pdg, nMothers, motherIdx, nDaughters, checkDoubleCounting,retIdx);
    //  if(returnVal!=-999) cout << "returning 1 " << returnVal << endl;
    return returnVal;
  }
  
  vector<int>  gIndex1 = cand->get_gIndex1();
  vector<int>  gIndex2 = cand->get_gIndex2();

  for(unsigned int n1 = 0; n1 < gIndex1.size(); n1++ )
    {
      for(unsigned int n2 = 0; n2 < gIndex2.size(); n2++ )
	{
	  int dau1 = gIndex1[n1];
	  int dau2 = gIndex2[n2];
	  
	  int debug=0;
	  if(debug){
	    if(nMothers[dau1]>0 && nMothers[dau2]>0){
	      cout << "nMothers: "<< nMothers[dau1] << " "<< nMothers[dau2] << endl;
	      cout << "pdgs: "<< pdg[dau1] << " "<< pdg[dau2] << endl;
	      cout << "motherIdx: "<< motherIdx[dau1][0] << " " << motherIdx[dau2][0] << endl;
	      cout << "motherPdgs: " << pdg[motherIdx[dau1][0]] << " " << pdg[motherIdx[dau2][0]] << endl;
	    }
	  }
	  if ( dau1 == -999 || dau2 == -999 ) continue;
	  //mathed to two different gen particles
	  if ( dau1 == dau2 ) continue;
	  // one kaon, one pion //321,211
	  if(!checkDoubleCounting && ! ((abs(pdg[dau1]) == 211 && abs(pdg[dau2]) == 321))  )  continue;
	  if(checkDoubleCounting && ! ((abs(pdg[dau1]) == 321 && abs(pdg[dau2]) == 211))  )  continue;
	  //only one mother
	  if (!( nMothers[dau1] == 1 && nMothers[dau2] == 1 ))   continue;
	  //mother index the same and not -999
	  //but only if we require a hadronic decay - else we have to iterate through until we find a d-meson
	  if ( (motherIdx[dau1][0] != motherIdx[dau2][0]) || motherIdx[dau1][0] == -999 ) continue;
	  
	  int allmotherindex = motherIdx[dau1][0];
	  //is D charged and just two daughters
	  if(! (abs (pdg[allmotherindex]) == 421)) continue;
	  if(nDaughters[allmotherindex] != 2)  continue;
	  if(debug){
	    cout << "passed!" << endl;
	  }
	  //cout << "candidate matched to gen" << endl;
	  if(!retIdx){
	    return pdg[allmotherindex];
	  }
	  else{
	    return allmotherindex;
	  }
	}
    }
  return -999;
}

vector<int> searchCQuarks(int idx1, vector<int> parentCandidates, int* t_nMother, int* t_genPdg, int (*t_motherIdx)[200]){
  //if the c-quark isn't in the parentage info, push it back
  bool checkFlag=false;
  if(abs(t_genPdg[idx1])==4){
    for(unsigned int ivec=0; ivec<parentCandidates.size(); ivec++){
      if(idx1 == parentCandidates.at(ivec)) checkFlag=true;
    }
    if(!checkFlag){
      parentCandidates.push_back(idx1);
    }
  }
  //else recursively iterate through the parentage
  for(int imom=0; imom<t_nMother[idx1]; imom++){	
    if(idx1>2 && t_motherIdx[idx1][imom] != idx1){
      int idx2 = t_motherIdx[idx1][imom];
      parentCandidates = searchCQuarks(idx2, parentCandidates, t_nMother, t_genPdg, t_motherIdx);
    }
  }
  return parentCandidates;
}

vector<int> hasDDecay(int idx1, vector<int> DMesons, vector<double> duplicates, vector<int> &indexes, int* t_nDaughters, int* t_genPdg, int (*t_daughterIdx)[200], float* t_genPt, bool hashadronicDecay, TH1D *withKpi, TH1D* withoutKpi){


  if(abs(t_genPdg[idx1])==413 || abs(t_genPdg[idx1])==421 || abs(t_genPdg[idx1])==431 || abs(t_genPdg[idx1])==411) {
    //check to ensure i dont put in the same meson twice
    bool checkFlag=false;
    for(unsigned int isize=0; isize<duplicates.size(); isize++){
      if(t_genPt[idx1] == duplicates.at(isize)) checkFlag=true;
    }
    if(!checkFlag){
      if(!hashadronicDecay){
	DMesons.push_back(abs(t_genPdg[idx1]));
	indexes.push_back(idx1);
	duplicates.push_back(t_genPt[idx1]);
	if(t_genPdg[idx1]==421) withoutKpi->Fill(t_genPt[idx1]);
      }
      else{
	if(abs(t_genPdg[idx1])==413){ //start with D*
	  //cout << "D* Ndaughters: "<< t_nDaughters[idx1] << " daughter indexes: "<< t_genPdg[t_daughterIdx[idx1][0]] << " " << t_genPdg[t_daughterIdx[idx1][1]] << endl;
	  if(t_nDaughters[idx1]==2){
	    int idx2 = t_daughterIdx[idx1][0];
	    if(t_nDaughters[idx2]==2 && ((abs(t_genPdg[t_daughterIdx[idx2][0]])==321 && abs(t_genPdg[t_daughterIdx[idx2][1]])==211))){
	      //  cout << "D* found!" << endl;
	      DMesons.push_back(abs(t_genPdg[idx1]));
	      indexes.push_back(idx1);
	      duplicates.push_back(t_genPt[idx1]);
	    }
	  }
	}
	if(abs(t_genPdg[idx1])==421){ //do D0
	  // cout << "D0 Ndaughters: "<< t_nDaughters[idx1] << " daughter indexes: "<< t_genPdg[t_daughterIdx[idx1][0]] << " " << t_genPdg[t_daughterIdx[idx1][1]] << endl;
	  if(t_nDaughters[idx1]==2 && (abs(t_genPdg[t_daughterIdx[idx1][0]])==321 && abs(t_genPdg[t_daughterIdx[idx1][1]])==211)){
	    //  cout << "D0 found!" << endl;
	    DMesons.push_back(abs(t_genPdg[idx1]));
	    indexes.push_back(idx1);
	    duplicates.push_back(t_genPt[idx1]);
	    withKpi->Fill(t_genPt[idx1]);
	  }
	}
	if(abs(t_genPdg[idx1])==411){ //do Dpm
	  // cout << "D+/- Ndaughters: "<< t_nDaughters[idx1] << " daughter indexes: "<< t_genPdg[t_daughterIdx[idx1][0]] << " " << t_genPdg[t_daughterIdx[idx1][1]] << " " << t_genPdg[t_daughterIdx[idx1][2]] << endl;
	  if(t_nDaughters[idx1]==3 && ((abs(t_genPdg[t_daughterIdx[idx1][0]])==321 && abs(t_genPdg[t_daughterIdx[idx1][1]])==211 && abs(t_genPdg[t_daughterIdx[idx1][2]])==211) || (abs(t_genPdg[t_daughterIdx[idx1][0]])==211 && abs(t_genPdg[t_daughterIdx[idx1][1]])==321 && abs(t_genPdg[t_daughterIdx[idx1][2]])==211))){
	    // cout << "Dpm found!" << endl;
	    DMesons.push_back(abs(t_genPdg[idx1]));
	    indexes.push_back(idx1);
	    duplicates.push_back(t_genPt[idx1]);
	  }
	}
	if(abs(t_genPdg[idx1])==431){ //do Ds
	  // cout << "Ds Ndaughters: "<< t_nDaughters[idx1] << " daughter indexes: "<< t_genPdg[t_daughterIdx[idx1][0]] << " " << t_genPdg[t_daughterIdx[idx1][1]] << endl;
	  if(t_nDaughters[idx1]==2 && ((abs(t_genPdg[t_daughterIdx[idx1][0]])==321 && abs(t_genPdg[t_daughterIdx[idx1][1]])==313) || (abs(t_genPdg[t_daughterIdx[idx1][0]])==333 && abs(t_genPdg[t_daughterIdx[idx1][1]])==211) || (abs(t_genPdg[t_daughterIdx[idx1][0]])==313 && abs(t_genPdg[t_daughterIdx[idx1][1]])==321) || (abs(t_genPdg[t_daughterIdx[idx1][0]])==211 && abs(t_genPdg[t_daughterIdx[idx1][1]])==333)) ){
	    //cout << "Ds found!" << endl;
	    DMesons.push_back(abs(t_genPdg[idx1]));
	    indexes.push_back(idx1);
	    duplicates.push_back(t_genPt[idx1]);
	  }
	}
      }
    } 
  }
  //cout << "start pdg: "<< t_genPdg[idx1] << endl;
  //cout << "nDaughters: " << t_nDaughters[idx1] << endl;
  for(int idaught=0; idaught<t_nDaughters[idx1]; idaught++){
    //cout << "daughter " << idaught << " pdg: " << t_genPdg[t_daughterIdx[idx1][idaught]] << endl;
    if(t_daughterIdx[idx1][idaught] != idx1 && t_daughterIdx[idx1][idaught]>0){
      int idx2 = t_daughterIdx[idx1][idaught];
      
      DMesons = hasDDecay(idx2, DMesons, duplicates, indexes, t_nDaughters, t_genPdg, t_daughterIdx, t_genPt, hashadronicDecay, withKpi, withoutKpi);
    }
    //cout << "check on return idx: "<< DMesons << endl;
  }
  return DMesons;
}

bool passDSelections(snglhfcand* cand){

  bool pass = false;

  //alpha0 changed from 0.5 to 0.2, and chi2 changed from 2.0 to 3.0
  if(cand->get_ffls3d() > 2 && //secondary vertex significance
     cand->get_falpha0() < 0.2 && //opening angle between svtx displacement and dmeson momentum
     cand->get_fprob() > 0.05 && //svtx probability
     cand->get_fdr() < 0.25 && //displacement ?
     cand->get_fchi2() < 3.0 && // chi2
     cand->get_fm()>1.65 && cand->get_fm()<2.2){ //&& // mass selections so ntuples arent huge
    //cand->get_fq1() != cand->get_fpt2()){ //ensure decay particles have opposite sign
    pass = true;
  }

  return pass;

}

double trigComb(bool *trg, int *pscl, double pt){

  double weight=0;

  if(trg[4] && pt>=100) weight = pscl[4];
  if(trg[3] && pt>=80 && pt<100) weight = pscl[3];
  if(trg[2] && pt>=60 && pt<80) weight = pscl[2];
  if(trg[1] && pt>=40 && pt<60) weight = pscl[1];
  if(trg[0] && pt>=20 && pt<40) weight = pscl[0];
  
  return weight;
  
}


void charmJetAnalyzer(std::string filelist, int startfile, int endfile, int isMC=0, int usePUsub=1){

  double minJetPt = 20;
  double minRawJetPt = 23;
  double maxJetEta = 2.5;
  
  int recalculatedEntries[2] = {0,0};
  int pthatbins[3] = {30,80,9999};
  double xsecs[3] = {3.378E-02, 4.412E-04, 0};
  int pthatEntries[2] = {251382,102823};
  bool recalculateEntries = true;
  if(recalculateEntries && isMC){
    cout << "recalculating entries in MC files" << endl;
    TChain *tc = new TChain("akPu3PFJetAnalyzer/t","");
    std::ifstream instrre(filelist.c_str(), std::ifstream::in);
    std::string filename;
    
    int ifile=0;
    while(ifile<startfile){ instrre>>filename; ifile++; }
    
    cout << "running from " << ifile << " to " << endfile << endl;
    while(ifile<endfile){
      
      ifile++;
      instrre>>filename;
      tc->Add(filename.c_str());
    }

    for(int i=0; i<2; i++){
      recalculatedEntries[i] = tc->GetEntries(Form("pthat>%d && pthat<%d",pthatbins[i],pthatbins[i+1]));
      cout << "entries between pthat " << pthatbins[i] << " and " << pthatbins[i+1] << ": " << recalculatedEntries[i] << endl;
    }
  }

  TFile *fout = NULL;
  if(isMC==1) fout=new TFile(Form("DMesonCJet_DpmEmbed_pPbMC_ppReco_akPu3PF_%d.root",endfile),"recreate");
  else if(isMC==2) fout=new TFile(Form("DMesonCJet_Ds2KstarKEmbed_pPbMC_ppReco_akPu3PF_%d.root",endfile),"recreate");
  else if(isMC==3) fout=new TFile(Form("DMesonCJet_Ds2PhiPi_pPbMC_ppReco_akPu3PF_%d.root",endfile),"recreate");
  else if(isMC==4) fout=new TFile(Form("DMesonCJet_Dstar_pPbMC_ppReco_akPu3PF_%d.root",endfile),"recreate");
  else if(isMC==5) fout=new TFile(Form("DMesonCJet_Dzero_pPbMC_ppReco_akPu3PF_%d.root",endfile),"recreate");
  else if(isMC==6) fout=new TFile(Form("DMesonCJet_CJetOnly_pPbMC_ppReco_akPu3PF_%d.root",endfile),"recreate");
  else fout=new TFile(Form("DMesonCJet_NoJetTrgCut_pPbdata_ppReco_akPu3PF_%d.root",endfile),"recreate");

  if(isMC==1) cout << "Assuming we're using the D plus/minus MC!!" << endl;
  else if(isMC==2) cout << "Assuming we're using the Ds -> K*K MC!! "<< endl;
  else if(isMC==3) cout << "Assuming we're using the Ds -> Pi Phi MC!! "<< endl;
  else if(isMC==4) cout << "Assuming we're using the D* MC!!" << endl;
  else if(isMC==5) cout << "Assuming we're using the D0 MC!!" << endl;
  else if(isMC==6) cout << "Assuming you have a basic C-Jet MC" << endl;
  else cout << "Assuming we're using the High-pT Dataset!" << endl;

  TH1D *D0withKpi = new TH1D("D0withKpi","",40,0,150); D0withKpi->Sumw2();
  TH1D *D0withoutKpi = new TH1D("D0withoutKpi","",40,0,150); D0withoutKpi->Sumw2();

  float t_jtpt[MAXJETS], t_jteta[MAXJETS], t_jtphi[MAXJETS], t_rawpt[MAXJETS], t_discr_prob[MAXJETS], t_discr_ssvHighEff[MAXJETS], t_discr_ssvHighPur[MAXJETS], t_discr_csvSimple[MAXJETS], t_svtxm[MAXJETS], t_svtxdl[MAXJETS], t_chargedMax[MAXJETS], t_chargedSum[MAXJETS], t_neutralMax[MAXJETS], t_neutralSum[MAXJETS], t_photonMax[MAXJETS], t_photonSum[MAXJETS], t_eSum[MAXJETS], t_muSum[MAXJETS];
  int t_nsvtx[MAXJETS], t_svtxntrk[MAXJETS];
  int t_nref, t_mult;
  int t_genMatch[MAXJETS];
  
  int t_genPdg[MAXPARTICLES], t_motherIdx[MAXPARTICLES][200], t_daughterIdx[MAXPARTICLES][200], t_nMothers[MAXPARTICLES], t_nDaughters[MAXPARTICLES], t_gensubid[MAXPARTICLES];
  float t_genpt[MAXPARTICLES], t_geneta[MAXPARTICLES], t_genphi[MAXPARTICLES];

  int nref, bin;
  vector<double> jtpt, jteta, jtphi, refpt, rawpt, refparton_flavorForB, discr_prob, discr_ssvHighEff, discr_ssvHighPur, discr_csvSimple, svtxm, nsvtx, svtxntrk, svtxdl, subid, chargedMax, chargedSum, neutralMax, neutralSum, photonMax, photonSum, eSum, muSum;
  vector<int> genMatch;
  int HLT_Jet20_NoJetID_v1, HLT_Jet40_NoJetID_v1, HLT_Jet60_NoJetID_v1, HLT_Jet80_NoJetID_v1, HLT_Jet100_NoJetID_v1;
  int HLT_Jet20_NoJetID_v1_Prescl, HLT_Jet40_NoJetID_v1_Prescl, HLT_Jet60_NoJetID_v1_Prescl, HLT_Jet80_NoJetID_v1_Prescl, HLT_Jet100_NoJetID_v1_Prescl;
  float HLT_JetObjectPt, t_pthat;
  double triggerPt, pthat, resCorr, weight;
  int pHBHENoiseFilter, pprimaryvertexFilter, pcollisionEventSelection;

  int t_refparton_flavorForB[100], t_subid[100];
  float t_refpt[100];
  hfcand_v1* hfcandidate = new hfcand_v1;

  vector<double> dGenPt, dCandPt, dCandMass, dCandEta, dCandPhi, dCandDr, dCandChildMass, dCandType, dCandCharge1, dCandCharge2;
  vector<double> dCandMatchGenPt_dau2, dCandMatchGenEta_dau2, dCandMatchGenPhi_dau2;
  vector<int> dCandMatchGenPdg_dau2, dCandMatchGenSubid_dau2;
  vector<double> dCandMatchGenPt_index1, dCandMatchGenEta_index1, dCandMatchGenPhi_index1;
  vector<int> dCandMatchGenPdg_index1, dCandMatchGenSubid_index1;
  vector<double> dCandMatchGenPt_index2, dCandMatchGenEta_index2, dCandMatchGenPhi_index2; 
  vector<int> dCandMatchGenPdg_index2, dCandMatchGenSubid_index2;
  
  vector<vector<int> > dCandMatchParentPdg;
  vector<vector<int> > dCandParentPartonIdx;

  vector<vector<int> > hasGenD, hasGenDwithKpi, hasGenDIdx, hasGenDwithKpiIdx;
  vector<int> dCandGenMatchPdg, dCandGenMatchDoubleCountedPdg;
  vector<int> dCandGenMatchIdx, dCandGenMatchDoubleCountedIdx;

  //Initialize the tree to be saved
  TTree *ct = new TTree("ct","");
  ct->Branch("nref",&nref);
  ct->Branch("jtpt",&jtpt);
  ct->Branch("jteta",&jteta);
  ct->Branch("jtphi",&jtphi);
  if(isMC) ct->Branch("refpt",&refpt);
  ct->Branch("rawpt",&rawpt);
  if(isMC) ct->Branch("refparton_flavorForB",&refparton_flavorForB);
  ct->Branch("discr_prob",&discr_prob);
  ct->Branch("discr_ssvHighEff",&discr_ssvHighEff);
  ct->Branch("discr_ssvHighPur",&discr_ssvHighPur);
  ct->Branch("discr_csvSimple",&discr_csvSimple);
  if(isMC){ 
    ct->Branch("genMatch",&genMatch);
    ct->Branch("hasGenD",&hasGenD);
    ct->Branch("hasGenDIdx",&hasGenDIdx);
    ct->Branch("hasGenDwithKpi",&hasGenDwithKpi);
    ct->Branch("hasGenDwithKpiIdx",&hasGenDwithKpiIdx);
  }
  ct->Branch("svtxm",&svtxm);
  ct->Branch("svtxdl",&svtxdl);
  ct->Branch("svtxntrk",&svtxntrk);
  ct->Branch("bin",&bin,"bin/I");
  ct->Branch("HLT_Jet20_NoJetID_v1",&HLT_Jet20_NoJetID_v1,"HLT_Jet20_NoJetID_v1/I");
  ct->Branch("HLT_Jet40_NoJetID_v1",&HLT_Jet40_NoJetID_v1,"HLT_Jet40_NoJetID_v1/I");
  ct->Branch("HLT_Jet60_NoJetID_v1",&HLT_Jet60_NoJetID_v1,"HLT_Jet60_NoJetID_v1/I");
  ct->Branch("HLT_Jet80_NoJetID_v1",&HLT_Jet80_NoJetID_v1,"HLT_Jet80_NoJetID_v1/I");
  ct->Branch("HLT_Jet100_NoJetID_v1",&HLT_Jet100_NoJetID_v1,"HLT_Jet100_NoJetID_v1/I");
  ct->Branch("HLT_Jet20_NoJetID_v1_Prescl",&HLT_Jet20_NoJetID_v1_Prescl,"HLT_Jet20_NoJetID_v1_Prescl/I");
  ct->Branch("HLT_Jet40_NoJetID_v1_Prescl",&HLT_Jet40_NoJetID_v1_Prescl,"HLT_Jet40_NoJetID_v1_Prescl/I");
  ct->Branch("HLT_Jet60_NoJetID_v1_Prescl",&HLT_Jet60_NoJetID_v1_Prescl,"HLT_Jet60_NoJetID_v1_Prescl/I");
  ct->Branch("HLT_Jet80_NoJetID_v1_Prescl",&HLT_Jet80_NoJetID_v1_Prescl,"HLT_Jet80_NoJetID_v1_Prescl/I");
  ct->Branch("HLT_Jet100_NoJetID_v1_Prescl",&HLT_Jet100_NoJetID_v1_Prescl,"HLT_Jet100_NoJetID_v1_Prescl/I");
  ct->Branch("triggerPt",&triggerPt,"triggerPt/D");
  if(isMC){
    ct->Branch("pthat",&pthat,"pthat/D");
    ct->Branch("subid",&subid,"subid/I");
    ct->Branch("resCorr",&resCorr,"resCorr/D");
  }
  ct->Branch("weight",&weight,"weight/D");
  /*ct->Branch("neutralMax",&neutralMax);
  ct->Branch("neutralSum",&neutralSum);
  ct->Branch("chargedMax",&chargedMax);
  ct->Branch("chargedSum",&chargedSum);
  ct->Branch("photonMax",&photonMax);
  ct->Branch("photonSum",&photonSum);
  ct->Branch("muSum",&muSum);
  ct->Branch("eSum",&eSum);*/

  ct->Branch("dGenPt",&dGenPt);
  ct->Branch("dCandPt",&dCandPt);
  ct->Branch("dCandMass",&dCandMass);
  ct->Branch("dCandEta",&dCandEta);
  ct->Branch("dCandPhi",&dCandPhi);
  ct->Branch("dCandDr",&dCandDr);
  ct->Branch("dCandChildMass",&dCandChildMass);
  ct->Branch("dCandType",&dCandType);
  ct->Branch("dCandCharge1",&dCandCharge1);
  ct->Branch("dCandCharge2",&dCandCharge2);

  if(isMC){
    ct->Branch("dCandMatchGenPt_dau2",&dCandMatchGenPt_dau2);
    ct->Branch("dCandMatchGenEta_dau2",&dCandMatchGenEta_dau2);
    ct->Branch("dCandMatchGenPhi_dau2",&dCandMatchGenPhi_dau2);
    ct->Branch("dCandMatchGenPdg_dau2",&dCandMatchGenPdg_dau2);
    ct->Branch("dCandMatchGenSubid_dau2",&dCandMatchGenSubid_dau2);

    ct->Branch("dCandMatchGenPt_index1",&dCandMatchGenPt_index1);
    ct->Branch("dCandMatchGenEta_index1",&dCandMatchGenEta_index1);
    ct->Branch("dCandMatchGenPhi_index1",&dCandMatchGenPhi_index1);
    ct->Branch("dCandMatchGenPdg_index1",&dCandMatchGenPdg_index1);
    ct->Branch("dCandMatchGenSubid_index1",&dCandMatchGenSubid_index1);

    ct->Branch("dCandMatchGenPt_index2",&dCandMatchGenPt_index2);
    ct->Branch("dCandMatchGenEta_index2",&dCandMatchGenEta_index2);
    ct->Branch("dCandMatchGenPhi_index2",&dCandMatchGenPhi_index2);
    ct->Branch("dCandMatchGenPdg_index2",&dCandMatchGenPdg_index2);
    ct->Branch("dCandMatchGenSubid_index2",&dCandMatchGenSubid_index2);

    ct->Branch("dCandMatchParentPdg",&dCandMatchParentPdg);

    ct->Branch("dCandParentPartonIdx",&dCandParentPartonIdx);
    ct->Branch("dCandGenMatchPdg",&dCandGenMatchPdg);
    ct->Branch("dCandGenMatchIdx",&dCandGenMatchIdx);
    ct->Branch("dCandGenMatchDoubleCountedPdg",&dCandGenMatchDoubleCountedPdg);
    ct->Branch("dCandGenMatchDoubleCountedIdx",&dCandGenMatchDoubleCountedIdx);
  }
  
  vector<double> genPt, genEta, genPhi;
  vector<int> genIdx, genPdg, genDauIdx1, genDauIdx2, genDauDau2;
  TTree *gTree = new TTree("gTree","");
  gTree->Branch("genPt",&genPt);
  gTree->Branch("genEta",&genEta);
  gTree->Branch("genPhi",&genPhi);
  gTree->Branch("genIdx",&genIdx); //event-level index of particle
  gTree->Branch("genPdg",&genPdg);
  gTree->Branch("genDauIdx1",&genDauIdx1);
  gTree->Branch("genDauIdx2",&genDauIdx2);
  gTree->Branch("genDauDau2",&genDauDau2);

  //**** Initialize reader tree variables ***** /
  /*  trigO *HLT_Jet_NoJetID_v1_trigObject[6];
  for(int i=0; i<6; i++){
    HLT_Jet_NoJetID_v1_trigObject[i] = new trigO;
    }*/
  //********************************************* /


  std::ifstream instr(filelist.c_str(), std::ifstream::in);
  std::string filename;
  
  TFile *fin;
  int ifile=0;
  while(ifile<startfile){ instr>>filename; ifile++; }

  cout << "running from " << ifile << " to " << endfile << endl;
  while(ifile<endfile){

    ifile++;
    instr>>filename;
    
    std::cout << "File: " << filename << std::endl;
    fin = TFile::Open(filename.c_str());
    if(!fin){
      cout << "Warning! File not available! Skipping..." << endl;
      continue;
    }

    TTree *t;
    TTree *hftree;
    TTree *HltTree;
    //TTree *HltRerun;
    TTree *evtTree;
    TTree *skimTree;
    TTree *genTree=0;
    TNtuple *HltObject;

    if(usePUsub) t = (TTree*) fin->Get("akPu3PFJetAnalyzer/t");
    else t = (TTree*) fin->Get("ak3PFJetAnalyzer/t");
  
    if(isMC){
      hftree = (TTree*)fin->Get("HFtree/hftree");
      genTree = (TTree*)fin->Get("HiGenParticleAna/hi");
    }
    else{
      hftree = (TTree*)fin->Get("HFtree/hftree");
    }
    HltTree = (TTree*)fin->Get("hltanalysis/HltTree");
    HltObject = (TNtuple*)fin->Get("hltobject/jetObjTree");
    //HltRerun = (TTree*)fin->Get("hltReRun/hltTree");
    evtTree = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
    skimTree = (TTree*)fin->Get("skimanalysis/HltTree");

    if(t->GetEntries() != hftree->GetEntries() || t->GetEntries() != HltTree->GetEntries() || t->GetEntries() != evtTree->GetEntries() || t->GetEntries() != skimTree->GetEntries()){
      cout << "WARNING! TREES HAVE DIFFERENT NENTRIES! DESYNCED??" << endl;
      exit(0);
    }
    
    HltTree->SetBranchAddress("HLT_PAJet20_NoJetID_v1",&HLT_Jet20_NoJetID_v1);
    HltTree->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&HLT_Jet40_NoJetID_v1);
    HltTree->SetBranchAddress("HLT_PAJet60_NoJetID_v1",&HLT_Jet60_NoJetID_v1);
    HltTree->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&HLT_Jet80_NoJetID_v1);
    HltTree->SetBranchAddress("HLT_PAJet100_NoJetID_v1",&HLT_Jet100_NoJetID_v1);
    HltTree->SetBranchAddress("HLT_PAJet20_NoJetID_v1_Prescl",&HLT_Jet20_NoJetID_v1_Prescl);
    HltTree->SetBranchAddress("HLT_PAJet40_NoJetID_v1_Prescl",&HLT_Jet40_NoJetID_v1_Prescl);
    HltTree->SetBranchAddress("HLT_PAJet60_NoJetID_v1_Prescl",&HLT_Jet60_NoJetID_v1_Prescl);
    HltTree->SetBranchAddress("HLT_PAJet80_NoJetID_v1_Prescl",&HLT_Jet80_NoJetID_v1_Prescl);
    HltTree->SetBranchAddress("HLT_PAJet100_NoJetID_v1_Prescl",&HLT_Jet100_NoJetID_v1_Prescl);
    HltObject->SetBranchAddress("pt",&HLT_JetObjectPt);
    if(!isMC) skimTree->SetBranchAddress("pPAcollisionEventSelectionPA",&pcollisionEventSelection);
    skimTree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
    //skimTree->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
    evtTree->SetBranchAddress("hiBin",&bin);
    if(isMC){
      t->SetBranchAddress("refparton_flavorForB",t_refparton_flavorForB);
      t->SetBranchAddress("refpt",t_refpt);
      t->SetBranchAddress("subid",t_subid);
      t->SetBranchAddress("matchedGenID",t_genMatch);
      t->SetBranchAddress("pthat",&t_pthat);
    }
    hftree->SetBranchAddress("hfcandidate",&hfcandidate);

    t->SetBranchAddress("nref",&t_nref);
    t->SetBranchAddress("rawpt",t_rawpt);
    t->SetBranchAddress("jtpt",t_jtpt);
    t->SetBranchAddress("jteta",t_jteta);
    t->SetBranchAddress("jtphi",t_jtphi);
    t->SetBranchAddress("discr_ssvHighEff",t_discr_ssvHighEff);
    t->SetBranchAddress("discr_ssvHighPur",t_discr_ssvHighPur);
    t->SetBranchAddress("discr_csvSimple",t_discr_csvSimple);
    t->SetBranchAddress("discr_prob",t_discr_prob);
    t->SetBranchAddress("nsvtx",t_nsvtx);
    t->SetBranchAddress("svtxntrk",t_svtxntrk);
    t->SetBranchAddress("svtxdl",t_svtxdl);
    t->SetBranchAddress("svtxm",t_svtxm);
    t->SetBranchAddress("chargedMax",t_chargedMax);
    t->SetBranchAddress("photonMax",t_photonMax);
    t->SetBranchAddress("neutralMax",t_neutralMax);
    t->SetBranchAddress("chargedSum",t_chargedSum);
    t->SetBranchAddress("photonSum",t_photonSum);
    t->SetBranchAddress("neutralSum",t_neutralSum);
    t->SetBranchAddress("muSum",t_muSum);
    t->SetBranchAddress("eSum",t_eSum);

    if(isMC){
      genTree->SetBranchAddress("mult",&t_mult);
      genTree->SetBranchAddress("pdg",t_genPdg);
      genTree->SetBranchAddress("pt",t_genpt);
      genTree->SetBranchAddress("eta",t_geneta);
      genTree->SetBranchAddress("phi",t_genphi);
      genTree->SetBranchAddress("sube",t_gensubid);
      genTree->SetBranchAddress("motherIdx",t_motherIdx);
      genTree->SetBranchAddress("daughterIdx",t_daughterIdx);
      genTree->SetBranchAddress("nMothers",t_nMothers);
      genTree->SetBranchAddress("nDaughters",t_nDaughters);
    }
    
    int nentries = t->GetEntries();
    for(int ientry=0; ientry<nentries; ientry++){
      if (ientry%100000==0) cout<<" i = "<<ientry<<" out of "<<nentries<<" ("<<(int)(100*(float)ientry/(float)nentries)<<"%)"<<endl;
      
      hftree->GetEntry(ientry);
      t->GetEntry(ientry);
      evtTree->GetEntry(ientry);
      HltTree->GetEntry(ientry);
      //HltRerun->GetEntry(ientry);
      skimTree->GetEntry(ientry);
      if(isMC) genTree->GetEntry(ientry);

      if(!isMC){
       	if(!pHBHENoiseFilter || !pcollisionEventSelection) continue;
      }
      if(!HLT_Jet20_NoJetID_v1 && !HLT_Jet40_NoJetID_v1 && !HLT_Jet60_NoJetID_v1 && !HLT_Jet80_NoJetID_v1 && !HLT_Jet100_NoJetID_v1) continue;

      //Fill all the gen tree stuff!
      if(isMC){
	for(int igen=0; igen<t_mult; igen++){
	  if(t_nDaughters[igen]<2 || t_nDaughters[igen]>3) continue;
	  if(abs(t_genPdg[igen])==411 || abs(t_genPdg[igen])==413 || abs(t_genPdg[igen])==421 || abs(t_genPdg[igen])==431){
	    int gindex1 = t_daughterIdx[igen][0];
	    int gindex2 = t_daughterIdx[igen][1];
	    int gindex_dau2 = -1;
	    if(t_nDaughters[igen]==3){ 
	      gindex_dau2 = t_daughterIdx[igen][2];
	    }
	    //snglhfcand *temp = convertGenToHFCand(t_genpt[igen], t_geneta[igen], t_genphi[igen], t_genPdg[igen], gindex1, gindex2, gindex_dau2);
	    //int pb = findGenD0(temp, t_genPdg, t_nMothers, t_motherIdx, t_nDaughters,1,0);
	    // int pb2 = findGenD0(temp, t_genPdg, t_nMothers, t_motherIdx, t_nDaughters,0,0);
	    //if((abs(pb)>400 && abs(pb)<500)){// || (abs(pb2)>400 && abs(pb2)<500)){
	    genPt.push_back(t_genpt[igen]);
	    genPhi.push_back(t_genphi[igen]);
	    genEta.push_back(t_geneta[igen]);
	    genIdx.push_back(igen);
	    genPdg.push_back(t_genPdg[igen]);
	    genDauIdx1.push_back(gindex1);
	    genDauIdx2.push_back(gindex2);
	    genDauDau2.push_back(gindex_dau2);
	    // }
	  }
	}
      }
      
      bool isNoise=false;
      for(int ij=0; ij<t_nref; ij++){
	if(t_jtpt[ij]>minJetPt && TMath::Abs(t_jteta[ij])<maxJetEta){
	  //if((t_chargedSum[ij]+t_photonSum[ij]+t_neutralSum[ij]+t_muSum[ij]+t_eSum[ij])/t_jtpt[ij]>1.01){  ///OPTION 1
	  if(t_neutralMax[ij] / TMath::Max(t_chargedSum[ij],t_neutralSum[ij]) < 0.97){ ///OPTION 2
	    
	    isNoise=true;
	  }
	}
      }
      if(!isMC && isNoise) continue;

      //double trigPt[5][100];

      bool trgDec[5] = {(bool)HLT_Jet20_NoJetID_v1, (bool)HLT_Jet40_NoJetID_v1, (bool)HLT_Jet60_NoJetID_v1, (bool)HLT_Jet80_NoJetID_v1, (bool)HLT_Jet100_NoJetID_v1};
      int treePrescl[5] = {HLT_Jet20_NoJetID_v1_Prescl, HLT_Jet40_NoJetID_v1_Prescl, HLT_Jet60_NoJetID_v1_Prescl, HLT_Jet80_NoJetID_v1_Prescl, HLT_Jet100_NoJetID_v1_Prescl};
      /*int trgObjSize[5];
      for(int ii=0; ii<5; ii++){ trgObjSize[ii] = HLT_Jet_NoJetID_v1_trigObject[ii]->size();}
      //Fill the trigger Pt/Eta/Phi from the TriggerObjects
      for(int ii=0; ii<5; ii++){
	for(int iObj=0; iObj<trgObjSize[ii]; iObj++){
	  trigPt[ii][iObj] = HLT_Jet_NoJetID_v1_trigObject[ii]->at(iObj).pt();
	}
	}*/

      int maxtrg= -1;
      for(int ii=4; ii>=0; ii--){
	if(trgDec[ii]==1){
	  maxtrg=ii;
	  break;
	}
      }
      /* float maximumTrgPt=0;
      for(int ii=0; ii<5; ii++){
	if(trgDec[ii]){
	  for(int iObj=0; iObj<trgObjSize[ii]; iObj++){
	    if(trigPt[ii][iObj]>maximumTrgPt) maximumTrgPt = trigPt[ii][iObj];
	  } 
	}
	}*/

      if(recalculateEntries){
	for(int i=0; i<2; i++){
	  pthatEntries[i] = recalculatedEntries[i];
	}
      }
      if(!isMC) weight = trigComb(trgDec, treePrescl, HLT_JetObjectPt);
      else{
	if(t_pthat>80) weight = xsecs[1]/pthatEntries[1];
	else weight = (xsecs[0]-xsecs[1])/pthatEntries[0];
      }

      nref = t_nref;

      for(int ijet=0; ijet<nref; ijet++){

	if(t_jtpt[ijet]>minJetPt && abs(t_jteta[ijet])<maxJetEta && t_rawpt[ijet]>minRawJetPt){
	  jtpt.push_back(t_jtpt[ijet]);
	  jteta.push_back(t_jteta[ijet]);
	  jtphi.push_back(t_jtphi[ijet]);
	  rawpt.push_back(t_rawpt[ijet]);
	  if(isMC) refpt.push_back(t_refpt[ijet]);
	  if(isMC) refparton_flavorForB.push_back(t_refparton_flavorForB[ijet]);
	  discr_prob.push_back(t_discr_prob[ijet]);
	  discr_ssvHighEff.push_back(t_discr_ssvHighEff[ijet]);
	  discr_ssvHighPur.push_back(t_discr_ssvHighPur[ijet]);
	  discr_csvSimple.push_back(t_discr_csvSimple[ijet]);
	  svtxm.push_back(t_svtxm[ijet]);
	  nsvtx.push_back(t_nsvtx[ijet]);
	  svtxntrk.push_back(t_svtxntrk[ijet]);
	  svtxdl.push_back(t_svtxdl[ijet]);
	  if(isMC) subid.push_back(t_subid[ijet]);
	  if(isMC) genMatch.push_back(t_genMatch[ijet]);
	  if(isMC){
	    vector<int> t_hasGenD;
	    vector<double> genPts, dummy2;
	    vector<int> t_hasGenDwithKpi;
	    vector<int> idxes, idxes_withKpi;
	    int tempD[5] = {0,0,0,0,0};
	    int tempDwKpi[5] = {0,0,0,0,0};
	    int idxD[5] = {0,0,0,0,0};
	    int idxDwithKpi[5] = {0,0,0,0,0};
	    t_hasGenD = hasDDecay(t_genMatch[ijet], t_hasGenD, genPts, idxes, t_nDaughters, t_genPdg, t_daughterIdx, t_genpt, 0, D0withKpi, D0withoutKpi);
	    t_hasGenDwithKpi = hasDDecay(t_genMatch[ijet], t_hasGenDwithKpi, dummy2, idxes_withKpi, t_nDaughters, t_genPdg, t_daughterIdx, t_genpt, 1, D0withKpi, D0withoutKpi);
	    dGenPt = genPts;
	    // cout << "hasGenDwithKpi" << endl;
	    for(unsigned iii=0; iii<t_hasGenDwithKpi.size(); iii++){
	      //cout << t_hasGenDwithKpi.at(iii) << " ";
	    }
	    //cout << endl;

	    //fun with bitwise operations
	    //5 bits, 1 or 0 corresponding to whether each d-meson is in the jet
	    //d*, d0, ds->phipi, ds->k*k, dpm
	    //pdgs = {413, 421, 431, 431, 411, -999};
	    // 1 0 0 1 0 = 18 means, e.g. that there's a D* and Ds->k*k in the jet cone.
	    if(t_hasGenD.size() != idxes.size()){ cout << "size mismatch!" << endl; exit(0); }
	    for(unsigned int ires=0; ires<t_hasGenD.size(); ires++){
	      if(t_hasGenD.at(ires)==413){ tempD[0] = 1; idxD[0] = idxes.at(ires); }
	      if(t_hasGenD.at(ires)==421){ tempD[1] = 1; idxD[1] = idxes.at(ires); }
	      if(t_hasGenD.at(ires)==431){ tempD[2] = 1; idxD[2] = idxes.at(ires); }
	      if(t_hasGenD.at(ires)==431){ tempD[3] = 1; idxD[3] = idxes.at(ires); }
	      if(t_hasGenD.at(ires)==411){ tempD[4] = 1; idxD[4] = idxes.at(ires); }
	    }
	    for(unsigned int ires=0; ires<t_hasGenDwithKpi.size(); ires++){
	      if(t_hasGenDwithKpi.at(ires)==413){ tempDwKpi[0] = 1; idxDwithKpi[0] = idxes_withKpi.at(ires); } 
	      if(t_hasGenDwithKpi.at(ires)==421){ tempDwKpi[1] = 1; idxDwithKpi[1] = idxes_withKpi.at(ires); } 
	      if(t_hasGenDwithKpi.at(ires)==431){ tempDwKpi[2] = 1; idxDwithKpi[2] = idxes_withKpi.at(ires); } 
	      if(t_hasGenDwithKpi.at(ires)==431){ tempDwKpi[3] = 1; idxDwithKpi[3] = idxes_withKpi.at(ires); } 
	      if(t_hasGenDwithKpi.at(ires)==411){ tempDwKpi[4] = 1; idxDwithKpi[4] = idxes_withKpi.at(ires); } 
	    }
	    vector<int> Dcands, DcandsWithKpi, DcandsIdx, DcandsWithKpiIdx;
	    for(int ijk=0; ijk<5; ijk++){
	      Dcands.push_back(tempD[ijk]);
	      DcandsWithKpi.push_back(tempDwKpi[ijk]);
	      DcandsIdx.push_back(idxD[ijk]);
	      DcandsWithKpiIdx.push_back(idxDwithKpi[ijk]);
	    }
	    hasGenD.push_back(Dcands);
	    hasGenDwithKpi.push_back(DcandsWithKpi);
	    hasGenDIdx.push_back(DcandsIdx);
	    hasGenDwithKpiIdx.push_back(DcandsWithKpiIdx);
	    
	    //cout << "DcandsWithKpi: "<< DcandsWithKpi << endl;
	  }
	  /* chargedMax.push_back(t_chargedMax[ijet]);
	  chargedSum.push_back(t_chargedSum[ijet]);
	  neutralMax.push_back(t_neutralMax[ijet]);
	  neutralSum.push_back(t_neutralSum[ijet]);
	  photonMax.push_back(t_photonMax[ijet]);
	  photonSum.push_back(t_photonSum[ijet]);
	  eSum.push_back(t_eSum[ijet]);
	  muSum.push_back(t_muSum[ijet]);*/
	}
      }
      //cout << "event: "<< ientry << " ncandidates: "<< hfcandidate->get_nhfcand() << endl;
      for(int icand=0; icand<hfcandidate->get_nhfcand(); icand++){

	snglhfcand* cand = hfcandidate->get_hfcand(icand);
	
	if(passDSelections(cand)){
	  // cout << "candidate number: "<< icand << ", pt: " << cand->get_fpt() << endl;
	  dCandPt.push_back(cand->get_fpt());
	  dCandMass.push_back(cand->get_fm());
	  dCandEta.push_back(cand->get_feta());
	  dCandPhi.push_back(cand->get_fphi());
	  dCandDr.push_back(cand->get_fdr());
	  dCandChildMass.push_back(cand->get_fmdau1());
	  dCandType.push_back(cand->get_type());
	  dCandCharge1.push_back(cand->get_fq1());
	  dCandCharge2.push_back(cand->get_fq2()); //bug in producer code!!
	  if(isMC){
	    int pb = findGenD0(cand, t_genPdg, t_nMothers, t_motherIdx, t_nDaughters,0,0);
	    //if(pb !=-999){ cout << "candidate instance: "<< icand << " pushing back " << pb << endl;
	    //   cout << "candidate pt: "<< cand->get_fpt() << endl;
	    // }
	    int pbidx = findGenD0(cand, t_genPdg, t_nMothers, t_motherIdx, t_nDaughters,0,1);
	    dCandGenMatchPdg.push_back(pb);
	    dCandGenMatchIdx.push_back(pbidx);
	    int pbDC = findGenD0(cand, t_genPdg, t_nMothers, t_motherIdx, t_nDaughters,1,0);
	    int pbDCidx = findGenD0(cand, t_genPdg, t_nMothers, t_motherIdx, t_nDaughters,1,1);
	    dCandGenMatchDoubleCountedPdg.push_back(pbDC);
	    dCandGenMatchDoubleCountedIdx.push_back(pbDCidx);
	    
	    int dau2_mom=0, idx1_mom=0, idx2_mom=0;
	    
	    //do soft particle daughter
	    int gIndex_dau2, gIndex1, gIndex2;
	    if(cand->get_gIndex_dau2().size()>0){
	      gIndex_dau2 = cand->get_gIndex_dau2().at(0);
	      dCandMatchGenPdg_dau2.push_back(t_genPdg[gIndex_dau2]);
	      dCandMatchGenPt_dau2.push_back(t_genpt[gIndex_dau2]);
	      dCandMatchGenEta_dau2.push_back(t_geneta[gIndex_dau2]);
	      dCandMatchGenPhi_dau2.push_back(t_genphi[gIndex_dau2]);
	      dCandMatchGenSubid_dau2.push_back(t_gensubid[gIndex_dau2]);
	      if(t_nMothers[gIndex_dau2]==0) dau2_mom=-999;
	      else if(t_nMothers[gIndex_dau2]==1) dau2_mom = t_genPdg[t_motherIdx[gIndex_dau2][0]];
	      else{
		cout << "warning! Final state candidate particle (dau2) has more than one mom!" << endl;
		cout << "debug: pdg: " << t_genPdg[gIndex_dau2] << " nMothers: "<< t_nMothers[gIndex_dau2] << endl;
		dau2_mom=-999;
	      }
	    }
	    else{
	      dCandMatchGenPdg_dau2.push_back(-999);
	      dCandMatchGenPt_dau2.push_back(-999);
	      dCandMatchGenEta_dau2.push_back(-999);
	      dCandMatchGenPhi_dau2.push_back(-999);
	      dCandMatchGenSubid_dau2.push_back(-999);
	    }
	    //do first daughter
	    if(cand->get_gIndex1().size()>0){
	      gIndex1 = cand->get_gIndex1().at(0);
	      dCandMatchGenPdg_index1.push_back(t_genPdg[gIndex1]);
	      dCandMatchGenPt_index1.push_back(t_genpt[gIndex1]);
	      dCandMatchGenEta_index1.push_back(t_geneta[gIndex1]);
	      dCandMatchGenPhi_index1.push_back(t_genphi[gIndex1]);
	      dCandMatchGenSubid_index1.push_back(t_gensubid[gIndex1]);
	      if(t_nMothers[gIndex1]==0) idx1_mom=-999;
	      else if(t_nMothers[gIndex1]==1) idx1_mom = t_genPdg[t_motherIdx[gIndex1][0]];
	      else{
		cout << "warning! Final state candidate particle (idx1) has more than one mom!" << endl;
		cout << "debug: pdg: " << t_genPdg[gIndex1] << " nMothers: "<< t_nMothers[gIndex1] << endl;
		idx1_mom=-999;
	      }
	    }
	    else{
	      dCandMatchGenPdg_index1.push_back(-999);
	      dCandMatchGenPt_index1.push_back(-999);
	      dCandMatchGenEta_index1.push_back(-999);
	      dCandMatchGenPhi_index1.push_back(-999);
	      dCandMatchGenSubid_index1.push_back(-999);
	    }
	    //do second daughter
	    if(cand->get_gIndex2().size()>0){
	      gIndex2 = cand->get_gIndex2().at(0);
	      dCandMatchGenPdg_index2.push_back(t_genPdg[gIndex2]);
	      dCandMatchGenPt_index2.push_back(t_genpt[gIndex2]);
	      dCandMatchGenEta_index2.push_back(t_geneta[gIndex2]);
	      dCandMatchGenPhi_index2.push_back(t_genphi[gIndex2]);
	      dCandMatchGenSubid_index2.push_back(t_gensubid[gIndex2]);
	      if(t_nMothers[gIndex2]==0) idx2_mom=-999;
	      else if(t_nMothers[gIndex2]==1) idx2_mom = t_genPdg[t_motherIdx[gIndex2][0]];
	      else{
		cout << "warning! Final state candidate particle (idx2) has more than one mom!" << endl;
		cout << "debug: pdg: " << t_genPdg[gIndex2] << " nMothers: "<< t_nMothers[gIndex2] << endl;
		idx2_mom=-999;
	      }
	    }
	    else{
	      dCandMatchGenPdg_index2.push_back(-999);
	      dCandMatchGenPt_index2.push_back(-999);
	      dCandMatchGenEta_index2.push_back(-999);
	      dCandMatchGenPhi_index2.push_back(-999);
	      dCandMatchGenSubid_index2.push_back(-999);
	    }
	    
	    vector<int> tempArr;
	    tempArr.push_back(dau2_mom);
	    tempArr.push_back(idx1_mom);
	    tempArr.push_back(idx2_mom);
	    dCandMatchParentPdg.push_back(tempArr);
	    tempArr.clear();

	    //recursivelty iterate through all possible parentage combinations to get all c-quark parents.
	    int startingIndex=-999;
	    if(cand->get_gIndex1().size()>0) startingIndex = gIndex1;
	    else if(cand->get_gIndex2().size()>0) startingIndex = gIndex2;
	    else if(cand->get_gIndex_dau2().size()>0) startingIndex = gIndex_dau2;
	    vector<int> parentCandidates;
	    if(startingIndex>0) parentCandidates = searchCQuarks(startingIndex, parentCandidates, t_nMothers, t_genPdg, t_motherIdx);
	    else parentCandidates.push_back(-999);
	    dCandParentPartonIdx.push_back(parentCandidates);
	  }
	}
      }
      if(dCandPt.size()>0 || jtpt.size()>0){
	ct->Fill();
	gTree->Fill();
      }

      //clear out all stuff from vector collections!
      genPt.clear();
      genPhi.clear();
      genIdx.clear();
      genPdg.clear();
      genEta.clear();
      genDauIdx1.clear();
      genDauIdx2.clear();
      genDauDau2.clear();

      jtpt.clear();
      jteta.clear();
      jtphi.clear();
      rawpt.clear();
      refpt.clear();
      refparton_flavorForB.clear();
      discr_prob.clear();
      discr_ssvHighEff.clear();
      discr_ssvHighPur.clear();
      discr_csvSimple.clear();
      svtxm.clear();
      nsvtx.clear();
      svtxntrk.clear();
      svtxdl.clear();
      subid.clear();
      genMatch.clear();
      hasGenD.clear();
      hasGenDIdx.clear();
      hasGenDwithKpi.clear();
      hasGenDwithKpiIdx.clear();

      dCandPt.clear();
      dCandMass.clear();
      dCandEta.clear();
      dCandPhi.clear();
      dCandDr.clear();
      dCandChildMass.clear();
      dCandType.clear();
      dCandCharge1.clear();
      dCandCharge2.clear();
      if(isMC){
	dCandMatchGenPdg_dau2.clear();
	dCandMatchGenPdg_index1.clear();
	dCandMatchGenPdg_index2.clear();
	dCandMatchGenPt_dau2.clear();
	dCandMatchGenPt_index1.clear();
	dCandMatchGenPt_index2.clear();
	dCandMatchGenEta_dau2.clear();
	dCandMatchGenEta_index1.clear();
	dCandMatchGenEta_index2.clear();
	dCandMatchGenPhi_dau2.clear();
	dCandMatchGenPhi_index1.clear();
	dCandMatchGenPhi_index2.clear();
	dCandMatchGenSubid_dau2.clear();
	dCandMatchGenSubid_index1.clear();
	dCandMatchGenSubid_index2.clear();
	dCandMatchParentPdg.clear();
	dCandParentPartonIdx.clear();
	dCandGenMatchPdg.clear();
	dCandGenMatchIdx.clear();
	dCandGenMatchDoubleCountedPdg.clear();
	dCandGenMatchDoubleCountedIdx.clear();
      }
      
      hfcandidate->Reset();
      
    }

    fin->Close();
  }
  fout->cd();
  ct->Write();
  gTree->Write();
  D0withKpi->Write();
  D0withoutKpi->Write();
  fout->Close();
}
