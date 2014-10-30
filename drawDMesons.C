#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TF1.h"

using namespace std;

//annoying global stuff
const int ptBins = 5;
const int jetPtBins = 6;
const int djetBins = 4;
const int cTypes = 5;
const double tightSignalCutLow[cTypes] = {0.144, 1.84, 1.93, 1.93, 1.85};
const double tightSignalCutHi[cTypes] = {0.147, 1.885, 1.99, 1.99, 1.885};
int globalCharmType = -1;
const double branchingRatio[cTypes] = {0.024,0.0388,0.0228,0.0263,0.0913};
const double MesonMass[cTypes] = {0.1455,1.8648,1.9685,1.9685,1.8696};

//draw fitting for a poly2 + gaussian
double draw_fit(int ich, TH1D* hfitter, int nrb, float rlow, float rhigh, float& sig, float& err)
{
  hfitter->Rebin(nrb);
  // h->SetMarkerSize(0.8);
  //h->SetLineColor(2);
  //h->SetMarkerColor(2);
  //h->SetMarkerStyle(20);
  //h->GetXaxis()->SetNdivisions(505);

  hfitter->GetXaxis()->SetRangeUser(rlow, rhigh);
  //.. fit with a Gaussian and pol or expo
  TF1* fit_fun = new TF1("fit_fun", "[0]*(1/[2]/sqrt(6.28)*exp(-0.5*pow((x-[1])/[2], 2))) + pol2(3) + [6]*(1/[8]/sqrt(6.28)*exp(-0.5*pow((x-[7])/[8], 2)))", rlow, rhigh);
  //TF1* fit_fun = new TF1("fit_fun", "[0]*(1/[2]/sqrt(6.28)*exp(-0.5*pow((x-[1])/[2], 2))) + expo(3) + [6]*(1/[8]/sqrt(6.28)*exp(-0.5*pow((x-[7])/[8], 2)))", rlow, rhigh);
  float tot = hfitter->Integral();
  
  float var_mean = 0.01, var_width = 0.01;
  float p0 = tot, p1 = MesonMass[ich], p2 = 0.015;
  if(ich==0) {
    var_mean = 0.005;
    var_width = 5e-4;
    p2 = 1e-3;
  }
  int pass = 1;
  float p6 = tot, p7 = MesonMass[4], p8 = 0.015; //..2nd peak is D+/-
  float zero = 0;
  while(pass) {
    //.. initial value 
    fit_fun->SetParameter(0, p0);
    fit_fun->SetParameter(1, p1);
    fit_fun->SetParameter(2, p2);
    if(ich==2) {//.. 2nd gaussian for Ds->phi+pi
      fit_fun->SetParameter(6, p6);
      fit_fun->SetParameter(7, p7);
      fit_fun->SetParameter(8, p8);
    }
    
    //.. fit constraint ..
    //fit_fun->SetParLimits(0, 0, 1e20);
    fit_fun->SetParLimits(1, TMath::Max(zero, (float)(MesonMass[ich]-var_mean)), MesonMass[ich]+var_mean);
    fit_fun->SetParLimits(2, TMath::Max(zero, (float)(p2-var_width)),  p2+var_width);
    
    if(ich==2) {
      //fit_fun->SetParLimits(6, 0, 1e20);
      fit_fun->SetParLimits(7, TMath::Max(zero, (float)(MesonMass[4]-var_mean)), MesonMass[4]+var_mean);
      fit_fun->SetParLimits(8, TMath::Max(zero, (float)(p8-var_width)),  p8+var_width);
    } else {//.. all D except Ds->phi+pi have only one peak. remove the 2nd peak
      fit_fun->FixParameter(6, 0);
      fit_fun->FixParameter(7, 0);
      fit_fun->FixParameter(8, 0);
    }
    
    //.. using likelyhood. Chi2 bias total counts
    for(int ii = 0; ii<10; ii++)
      hfitter->Fit(fit_fun,  "q", "", rlow, rhigh);
    //!!!-- I can't get the fit error from likelihood fit right 
    //!!!-- the mean is usable for systematic check
    //h->Fit(fit_fun,  "WL", "", rlow, rhigh);
    
    //.. draw foreground and background ..
    cout<<" hist: "<<hfitter->GetName()<<endl;
    hfitter->Draw();
    
    TF1* fit_fun_bg = (TF1*)fit_fun->Clone("fit_fun_bg");
    fit_fun_bg->SetParameter(0, 0);
    fit_fun_bg->SetParameter(1, 0);
    fit_fun_bg->SetParameter(2, 0);
    fit_fun_bg->SetParameter(6, 0);
    fit_fun_bg->SetParameter(7, 0);
    fit_fun_bg->SetParameter(8, 0);
    
    fit_fun_bg->SetLineColor(8);
    fit_fun_bg->Draw("same");
    gPad->Update();
    
    pass = 0;
    //.. check if need to fit again ...
    /* cout<<" good fit? (0: no refit, 1: refit w/o range change, 2: refit w/ range change, 3: change both range and fit pars"<<endl;
    scanf("%d", &pass);
    
    if(pass==1) {//.. change the par range only
      cout<<" var(mean), var(width) ?"<<endl;
      scanf("%f%f", &var_mean, &var_width);
    } else if(pass==2) {//.. change the fitting range only
      cout<<" rlow, rhigh? "<<endl;
      scanf("%f%f", &rlow, &rhigh);
    } else if(pass==3) {// change both par and fitting range
      cout<<" var(mean), var(width) ?"<<endl;
      scanf("%f%f", &var_mean, &var_width);
      cout<<" rlow, rhigh? "<<endl;
      scanf("%f%f", &rlow, &rhigh);
      }*/
  }
  
  // correct mass bin width
  sig = fit_fun->GetParameter(0)/hfitter->GetBinWidth(1);
  err = fit_fun->GetParError(0)/hfitter->GetBinWidth(1);
  
  cout<<"total number of D: "<<sig<<"+/-"<<err<<endl;
  return sig;
}

//rejects signal region when fitting a pol2 
Double_t poly2reject(Double_t *x, Double_t *par){

  if(x[0] > tightSignalCutLow[globalCharmType] && x[0] < tightSignalCutHi[globalCharmType]){
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0] + par[2]*pow(x[0],2);
}

TH1D *divideWidths(TH1D* input){

  for(int i=1; i<=input->GetNbinsX(); i++){
    input->SetBinError(i, input->GetBinError(i)/input->GetBinWidth(i));
    input->SetBinContent(i, input->GetBinContent(i)/input->GetBinWidth(i));
  }
  return input;
}

void drawDMesons(bool drawMesonMass = 0, bool doRebin=0, int rebinVal=2, int isMC=1, bool doDCStudy=false){

  int jetPtBoundaries[jetPtBins+1] = {20,40,60,80,100,200,400};
  //int jetPtBoundaries[jetPtBins+1] = {20,400)};
  //int mesonPtBoundaries[jetPtBins+1+3] = {0,10,15,20,40,400};
  // int ptBoundaries[ptBins+1] = {0,10,15,20,40,400};
  // double djetBoundaries[djetBins+1] = {0,0.1,0.2,0.5,1.};

  if(!isMC) doDCStudy = false;

  const int purityBins = 10;
  double xpurbins[purityBins+1] = {0,5,10,15,20,40,60,80,100,200,400};

  string dtypes[cTypes] = {"D*","D^{0}","D_{s} -> #phi+#pi","D_{s} -> K*K","D^{+/-}"};
  string xTitle[cTypes] = {"D* - D^{0} Mass", "D^{0} Mass", "D_{s} Mass", "D_{s} Mass", "D^{+/-} Mass"}; 

  TH1D *spec[jetPtBins][djetBins][cTypes];
  TH1D *bgspec[jetPtBins][djetBins][cTypes];
  TH1D *sigspec[jetPtBins][djetBins][cTypes];
  TH1D *dcspec[jetPtBins][djetBins][cTypes];
  TH1D *nondcspec[jetPtBins][djetBins][cTypes];

  TH1D *dr[jetPtBins][cTypes];
  TH1D *bgdr[jetPtBins][cTypes];

  TF1 *fits[jetPtBins][cTypes];
  TF1 *fitDraw[jetPtBins][cTypes];

  TH1D *purity[cTypes];
  for(int i=0; i<cTypes; i++){
    purity[i] = new TH1D(Form("purity_type%d",i),"",purityBins,xpurbins);
    purity[i]->Sumw2();
  }

  TFile *fin;
  if(!isMC) fin = new TFile("dCandidateToJetAssns_halfData_RevCuts_JetAssn-10-3_fixChgSign.root");
  if(isMC==1) fin = new TFile("dCandidateToJetAssns_CJetOnly_RevCuts_JetAssn-10-30.root"); //dCandidateToJetAssns_CJetOnly_RevCuts_JetAssn-10-13.root");
  char* histoname = new char[100];
  for(int i=0; i<jetPtBins; i++){
    int j=0; //for(int j=0; j<djetBins; j++){
      for(int k=0; k<cTypes; k++){
	histoname = Form("assns_unweight_pt%d_djet%d_type%d",i,j,k);
	spec[i][j][k] = (TH1D*)fin->Get(histoname)->Clone(histoname);
	spec[i][j][k]->SetXTitle(xTitle[k].c_str());
	spec[i][j][k]->SetYTitle("Counts");
	if(doRebin) spec[i][j][k]->Rebin(rebinVal);

	histoname = Form("bgAssns_pt%d_djet%d_type%d",i,j,k);
	bgspec[i][j][k] = (TH1D*)fin->Get(histoname)->Clone(histoname);
	bgspec[i][j][k]->SetXTitle(xTitle[k].c_str());
	bgspec[i][j][k]->SetYTitle("Counts");
	if(doRebin) bgspec[i][j][k]->Rebin(rebinVal);

	histoname = Form("proximity_pt%d_type%d",i,k);
	dr[i][k] = (TH1D*)fin->Get(histoname)->Clone(histoname);
	dr[i][k]->SetXTitle("Jet to Meson dR");
	dr[i][k]->SetYTitle("dN/d#DeltaR");
	if(doRebin) dr[i][k]->Rebin(rebinVal);

	histoname = Form("bgProximity_pt%d_type%d",i,k);
	bgdr[i][k] = (TH1D*)fin->Get(histoname)->Clone(histoname);
	bgdr[i][k]->SetXTitle("Jet to Meson dR");
	bgdr[i][k]->SetYTitle("dN/d#DeltaR");
	if(doRebin) bgdr[i][k]->Rebin(rebinVal);

	if(doDCStudy){
	  histoname = Form("DC_assns_pt%d_djet%d_type%d",i,j,k);
	  dcspec[i][j][k] = (TH1D*)fin->Get(histoname)->Clone(histoname);
	  dcspec[i][j][k]->SetXTitle(xTitle[k].c_str());
	  dcspec[i][j][k]->SetYTitle("Counts");
	  dcspec[i][j][k]->SetMarkerColor(kGreen+2);
	  dcspec[i][j][k]->SetLineColor(kGreen+2);
	  if(doRebin) dcspec[i][j][k]->Rebin(rebinVal);

	  histoname = Form("non_DC_assns_pt%d_djet%d_type%d",i,j,k);
	  nondcspec[i][j][k] = (TH1D*)fin->Get(histoname)->Clone(histoname);
	  nondcspec[i][j][k]->SetXTitle(xTitle[k].c_str());
	  nondcspec[i][j][k]->SetYTitle("Counts");
	  nondcspec[i][j][k]->SetMarkerColor(kOrange+2);
	  nondcspec[i][j][k]->SetLineColor(kOrange+2);
	  if(doRebin) nondcspec[i][j][k]->Rebin(rebinVal);
	}
	fits[i][k] = new TF1(Form("fitfg_pt%d_type%d",i,k),poly2reject,spec[i][j][k]->GetBinCenter(1),spec[i][j][k]->GetBinCenter(spec[i][j][k]->GetNbinsX()),3);
	fitDraw[i][k] = new TF1(Form("fitfgDraw_pt%d_type%d",i,k),"pol2",spec[i][j][k]->GetBinCenter(1),spec[i][j][k]->GetBinCenter(spec[i][j][k]->GetNbinsX()));
      }
      //}
  }
  
  double specScale[jetPtBins][cTypes];
  int nJets[jetPtBins][cTypes];
  float dmesonSig[jetPtBins][cTypes];
  float dmesonSigErr[jetPtBins][cTypes];

  double doubleCountingtop[jetPtBins][cTypes];
  double doubleCountingbot[jetPtBins][cTypes];
  
  TLatex *twrite[jetPtBins][djetBins][cTypes];
  TLatex *twrite2[cTypes][djetBins];
  TCanvas *cc[djetBins][cTypes];
  TCanvas *cc2[djetBins][cTypes];
  for(int k=0; k<cTypes; k++){
    globalCharmType = k;
    int j=0; //for(int j=0; j<djetBins; j++){
      cc[j][k] = new TCanvas(Form("cc_djet%d_type%d",j,k),Form("cc_djet%d_type%d",j,k),1000,800);
      cc2[j][k] = new TCanvas(Form("cc2_djet%d_type%d",j,k),Form("cc2_djet%d_type%d",j,k),1000,800);
      if(drawMesonMass) cc[j][k]->Divide(2,3);
      else cc[j][k]->Divide(2,3);
      cc2[j][k]->Divide(2,3);
      for(int i=0; i<jetPtBins; i++){
	if(dr[i][k]->Integral()==0) continue;
	//if(!drawMesonMass && i>=6) continue;
	if(i>6) continue;
	cc[j][k]->cd(i+1);
	//if(!drawMesonMass) gPad->SetLogy();
	//gStyle->SetErrorX(spec[i][j][k]->GetBinWidth(3)*40);
	int lbin = bgspec[i][j][k]->FindBin(0.149);
	int hbin = bgspec[i][j][k]->FindBin(0.154);
	float scaler = spec[i][j][k]->Integral(lbin,hbin)/bgspec[i][j][k]->Integral(lbin,hbin);
	bgspec[i][j][k]->Scale(scaler);	
	bgspec[i][j][k]->SetLineColor(4);

	spec[i][j][k]->Fit(fits[i][k],"qN0");
	fitDraw[i][k]->SetParameters(fits[i][k]->GetParameter(0),fits[i][k]->GetParameter(1),fits[i][k]->GetParameter(2));

	double range1 = spec[i][j][k]->GetBinLowEdge(1);
	double range2 = spec[i][j][k]->GetBinLowEdge(spec[i][j][k]->FindBin(tightSignalCutLow[k]));
	double range3 = spec[i][j][k]->GetBinLowEdge(spec[i][j][k]->FindBin(tightSignalCutHi[k])+1);
	double range4 = spec[i][j][k]->GetBinLowEdge(spec[i][j][k]->GetNbinsX())+spec[i][j][k]->GetBinWidth(spec[i][j][k]->GetNbinsX());
	//Scale the dR plots based on the signal region bg integral over the sideband region bg integral
	double region1Int = fitDraw[i][k]->Integral(range1, range2);
	double region2Int = fitDraw[i][k]->Integral(range3, range4);
	double region3Int = fitDraw[i][k]->Integral(range2, range3);
	double scaler3 = region3Int/(region1Int+region2Int);
	//cout << "type " << k << " pt: "<< i << " " << region1Int/spec[i][j][k]->GetBinWidth(1) << " " << region2Int/spec[i][j][k]->GetBinWidth(1) << " " << region3Int/spec[i][j][k]->GetBinWidth(1) << endl;
	//cout << "norm: "<< scaler3 << endl;
	//cout << "ensure reg1+2 = bgdr: "<< bgdr[i][k]->Integral() << " " << (region1Int+region2Int)/spec[i][j][k]->GetBinWidth(1) << endl;
	
	//original mass-based bg scaling method:
	bgdr[i][k]->Scale(scaler3);

	//dr-based bg scaling:
	int lowbin = dr[i][k]->FindBin(0.65);
	int highbin = dr[i][k]->FindBin(1.99);
	//scaler3 = dr[i][k]->Integral(lowbin,highbin) / bgdr[i][k]->Integral(lowbin,highbin);
	//bgdr[i][k]->Scale(scaler3);
	
	//try a lee-yang zeroes-esque method
	double maxdiff=999;
	for(int ibin=1; ibin<=dr[i][k]->GetNbinsX(); ibin++){
	  if(dr[i][k]->GetBinContent(ibin)/bgdr[i][k]->GetBinContent(ibin) < maxdiff) maxdiff = dr[i][k]->GetBinContent(ibin)/bgdr[i][k]->GetBinContent(ibin);
	}
	//scaler3 = maxdiff;
	//bgdr[i][k]->Scale(scaler3);
	bgdr[i][k]->SetLineColor(4);

	//REWORKED FIT WITH GAUSSIAN ON POLY2
	double realfitted = draw_fit(k,spec[i][j][k],1,(float)range1,(float)range4,dmesonSig[i][k],dmesonSigErr[i][k]);
	if(doDCStudy){
	  cc2[j][k]->cd(i+1);
	  double dcbot = draw_fit(k,nondcspec[i][j][k],1,(float)range1,(float)range4,dmesonSig[i][k],dmesonSigErr[i][k]);
	  dcspec[i][j][k]->Draw("same");
	  
	  doubleCountingtop[i][k] = dcspec[i][j][k]->Integral();
	  doubleCountingbot[i][k] = dcbot;
	}	
	cc[j][k]->cd(i+1);

	specScale[i][k] = scaler3;
	nJets[i][k] = region3Int;

	int tlbin = spec[i][j][k]->FindBin(tightSignalCutLow[k]);
	int thbin = spec[i][j][k]->FindBin(tightSignalCutHi[k]);
	cout << "tlbin: "<< tlbin << " thbin: "<< thbin << endl;
	cout << "spec: "<< spec[i][j][k]->Integral(tlbin,thbin) << endl;
	cout << "type: "<< k << " pt: "<< i << " fg integral: "<< spec[i][j][k]->Integral(tlbin,thbin,"") << " bg: "<< region3Int/spec[i][j][k]->GetBinWidth(1) << " signal: " << spec[i][j][k]->Integral(tlbin,thbin,"")-region3Int/spec[i][j][k]->GetBinWidth(1) << " signal to bg: " << (spec[i][j][k]->Integral(tlbin,thbin,"") - region3Int/spec[i][j][k]->GetBinWidth(1)) / (region3Int/spec[i][j][k]->GetBinWidth(1)) << endl;

	if(drawMesonMass){
	  sigspec[i][j][k] = (TH1D*)spec[i][j][k]->Clone(Form("signal_pt%d_type%d",i,k));
	  spec[i][j][k]->Draw();
	  if(doDCStudy) dcspec[i][j][k]->Draw("same");
	  fitDraw[i][k]->Draw("same");
	  if(k==0) bgspec[i][j][k]->Draw("same");
	  sigspec[i][j][k]->Add(bgspec[i][j][k],-1);
	  sigspec[i][j][k]->SetLineColor(2);
	  if(k==0) sigspec[i][j][k]->Draw("hist,same");
	}
	else{
	  dr[i][k] = divideWidths(dr[i][k]);
	  bgdr[i][k] = divideWidths(bgdr[i][k]);
	  sigspec[i][j][k] = (TH1D*)dr[i][k]->Clone(Form("signal_pt%d_type%d",i,k));
	  dr[i][k]->SetMinimum(1);
	  dr[i][k]->GetXaxis()->SetRangeUser(0.,4.);
	  dr[i][k]->Draw();
	  bgdr[i][k]->Draw("same,hist");
	  sigspec[i][j][k]->Add(bgdr[i][k],-1);
	  sigspec[i][j][k]->SetLineColor(2);
	  sigspec[i][j][k]->Draw("hist,same");
	  
	}
	
	if(!drawMesonMass){
	  double purDRcut = 0.299;
	  //now calculate purity functions vs jet pt
	  double purContentTopErr = 0;
	  int bincc = sigspec[i][j][k]->FindBin(purDRcut);
	  cout << "bin cross-check" << sigspec[i][j][k]->GetBinLowEdge(bincc) << " " << sigspec[i][j][k]->GetBinLowEdge(bincc)+sigspec[i][j][k]->GetBinWidth(bincc) << endl;
	  double purContentTop = sigspec[i][j][k]->IntegralAndError(1,sigspec[i][j][k]->FindBin(purDRcut),purContentTopErr,""); //you have to use "width" option here if sigspec and dr have different binning schemes - we're luckily ok here

	  double purContentBotErr = 0;
	  double purContentBot = dr[i][k]->IntegralAndError(1,sigspec[i][j][k]->FindBin(purDRcut),purContentBotErr,"");
	  //cout << "consistency check: " << purity[k]->GetBinLowEdge(i+5) << " " << jetPtBoundaries[i] << endl;
	  if(purContentBot!=0){
	    purity[k]->SetBinContent(i+5,purContentTop/purContentBot);
	    purity[k]->SetBinError(i+5, sqrt(pow(purContentTopErr/purContentTop,2)+pow(purContentBotErr/purContentBot,2)));
	  }
	  else{
	    purity[k]->SetBinContent(i+5,0);
	    purity[k]->SetBinError(i+5,0);
	  }
	}
	else{
	  double purContentTopErr = 0;
	  double purContentTop;
	  //purContentTop = sigspec[i][j][k]->IntegralAndError(tlbin,thbin,purContentTopErr,"width");
	  purContentTop = region3Int;
	  //purContentTop = fitDraw[i][k]->Integral(range1,range4);
	  double purContentBotErr = 0;
	  double purContentBot = spec[i][j][k]->IntegralAndError(tlbin,thbin,purContentBotErr,"width");
	  //double purContentBot = spec[i][j][k]->IntegralAndError(1,spec[i][j][k]->GetNbinsX(),purContentBotErr,"width");
	  if(purContentBot!=0){
	    //if(k==0) purity[k]->SetBinContent(i+5,purContentTop/purContentBot);
	    purity[k]->SetBinContent(i+5,(purContentBot-purContentTop)/purContentBot);
	    purity[k]->SetBinError(i+5, sqrt(pow(purContentTopErr/purContentTop,2)+pow(purContentBotErr/purContentBot,2)));
	  }
	  else{
	    purity[k]->SetBinContent(i+5,0);
	    purity[k]->SetBinError(i+5,0);
	  }
	}

	//cout << "S/B cross-check, fg: "<< dr[i][k]->Integral() << " bg: "<< bgdr[i][k]->Integral() << " s/b: "<< (dr[i][k]->Integral("") - bgdr[i][k]->Integral("")) / bgdr[i][k]->Integral("") << endl;
	//cout << "crosscheck2 " << bgdr[i][k]->Integral() << " " << region3Int/spec[i][j][k]->GetBinWidth(1) << endl;

	cc[j][k]->SetLogy();
	cc[j][k]->Update();

	if(i==0) {
	  if(!drawMesonMass) twrite2[i][j] = new TLatex(0.532,gPad->GetUymax()*0.85,Form("Jet p_{T}: %d-%d, dType: %s",jetPtBoundaries[i],jetPtBoundaries[i+1], dtypes[k].c_str()));
	  else if(k==0) twrite2[i][j] = new TLatex(0.142,gPad->GetUymax()*0.85,Form("Jet p_{T}: %d-%d, dType: %s",jetPtBoundaries[i],jetPtBoundaries[i+1], dtypes[k].c_str()));
	  else twrite2[i][j] = new TLatex(1.8,gPad->GetUymax()*0.85,Form("Jet p_{T}: %d-%d, dType: %s",jetPtBoundaries[i],jetPtBoundaries[i+1], dtypes[k].c_str()));
	  twrite2[i][j]->SetTextSize(0.066);
	  twrite2[i][j]->Draw("same");
	}
	else{
	  if(!drawMesonMass) twrite[i][j][k] = new TLatex(0.532,gPad->GetUymax()*0.75, Form("p_{T}: %d-%d",jetPtBoundaries[i],jetPtBoundaries[i+1]));
	  else if(k==0) twrite[i][j][k] = new TLatex(0.142,gPad->GetUymax()*0.75, Form("p_{T}: %d-%d",jetPtBoundaries[i],jetPtBoundaries[i+1]));
	  else twrite[i][j][k] = new TLatex(1.8,gPad->GetUymax()*0.85, Form("p_{T}: %d-%d",jetPtBoundaries[i],jetPtBoundaries[i+1]));
	  twrite[i][j][k]->SetTextSize(0.07);
	  twrite[i][j][k]->Draw("same");
	}
      }
      cout << "type " << k << " purity avg: " << purity[k]->Integral("width")/(purity[k]->GetBinLowEdge(purity[k]->GetNbinsX())+purity[k]->GetBinWidth(purity[k]->GetNbinsX())-purity[k]->GetBinLowEdge(1)) << endl;
      //}
  }
  TH1D *allSpec = (TH1D*)fin->Get("allJetSpec");
  
  TH1D *cSpec[cTypes];
  TH1D *bgcSpec[cTypes];
  for(int k=0; k<cTypes; k++){
    cSpec[k] = (TH1D*)fin->Get(Form("CjetSpec_type%d",k));
    bgcSpec[k] = (TH1D*)fin->Get(Form("bgCjetSpec_type%d",k));
  }

  
  //Take weighted average of scale factors for background of spectra plots
  for(int k=0; k<cTypes; k++){
    /*  specScale[0][k]*=nJets[0][k];
    for(int i=1; i<ptBins; i++){
      specScale[0][k] += specScale[i][k]*nJets[i][k];
      nJets[0][k] += nJets[i][k];
    }
    specScale[0][k]/=nJets[0][k];*/
    cout << "type " << k << " scale factor " << specScale[0][k] << endl;
  }

  allSpec = divideWidths(allSpec);
  allSpec->SetMarkerStyle(20);
  for(int k=0; k<cTypes; k++){
    cSpec[k] = divideWidths(cSpec[k]);
    bgcSpec[k] = divideWidths(bgcSpec[k]);
    bgcSpec[k]->Scale(specScale[0][k]);
    cSpec[k]->SetMarkerStyle(20);
  }
  
  TCanvas *ccan = new TCanvas("ccan","",800,600);
  ccan->Divide(1,2);
  ccan->cd(1);
  allSpec->SetYTitle("dN/dp_{T}");
  allSpec->Draw();
  //cSpec[0]->Draw();
  for(int k=0; k<cTypes; k++){
    cSpec[k]->SetMarkerColor(k+2);
    cSpec[k]->SetLineColor(k+2);
    bgcSpec[k]->SetMarkerColor(k+2);
    cSpec[k]->Add(bgcSpec[k],-1);
    //cSpec[k]->Draw("same");
    //bgcSpec[k]->Draw("same");
  }
  TH1D *cspeccln = (TH1D*)cSpec[0]->Clone("cspeccln");
  //cspeccln->Scale(1./branchingRatio[0]);
  for(int k=1; k<cTypes; k++){
    //cSpec[k]->Scale(1./branchingRatio[k]);
    cspeccln->Add(cSpec[k]);
  }
  cspeccln->Draw("same");
  ccan->cd(2);
  TH1D *cspeccln2 = (TH1D*)cspeccln->Clone("cspeccln2");
  cspeccln2->Divide(allSpec);
  cspeccln2->SetYTitle("c-Jet Fraction");
  cspeccln2->SetXTitle("Jet p_{T}");
  cspeccln2->Draw();

  TCanvas *cpur = new TCanvas("cpur","",800,600);
  cpur->Divide(2,2);
  for(int i=0; i<5; i++){
    if(i==4) cpur->cd(4);
    else cpur->cd(i+1);
    purity[i]->SetXTitle("jet p_{T}");
    purity[i]->SetYTitle(Form("Purity, %s",dtypes[i].c_str()));
    purity[i]->Draw();
  }

  TFile *foutput;
  if(isMC) foutput = new TFile("purityHistos_MC.root","recreate");
  else foutput = new TFile("purityHistos_data.root","recreate");
  foutput->cd();
  for(int i=0; i<5; i++){
    purity[i]->GetYaxis()->SetRangeUser(0,1);
    purity[i]->Write();
  }

  cout << "** Double counting statistics **" << endl;
  for(int k=0; k<5; k++){
    if(k!=2 && k!=3){
      for(int i=1; i<jetPtBins; i++){
	doubleCountingtop[0][k]+=doubleCountingtop[i][k];
	doubleCountingbot[0][k]+=doubleCountingbot[i][k];
      }
      cout << "D-type " << k << " " << (doubleCountingbot[0][k]/doubleCountingtop[0][k])*100 << "%" << endl;
    }
  }
}
