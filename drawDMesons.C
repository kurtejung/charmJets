
//annoying global stuff
const int ptBins = 9;
const int jetPtBins = 6;
const int djetBins = 4;
const int cTypes = 5;
const double tightSignalCutLow[cTypes] = {0.144, 1.84, 1.93, 1.93, 1.85};
const double tightSignalCutHi[cTypes] = {0.147, 1.885, 1.99, 1.99, 1.885};
int globalCharmType = -1;
const double branchingRatio[cTypes] = {0.024,0.0388,0.0228,0.0263,0.0913};

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

void drawDMesons(bool drawMesonMass = 1, bool doRebin=0, int rebinVal=2, isMC==1){

  int jetPtBoundaries[jetPtBins+1] = {20,40,60,100,140,200,400};
  int ptBoundaries[ptBins+1] = {0,10,20,30,40,50,60,70,80,90};
  double djetBoundaries[djetBins+1] = {0,0.1,0.2,0.5,1};

  string dtypes[cTypes] = {"D*","D^{0}","D_{s} -> #phi+#pi","D_{s} -> K*K","D^{+/-}"};
  string xTitle[cTypes] = {"D* - D^{0} Mass", "D^{0} Mass", "D_{s} Mass", "D_{s} Mass", "D^{+/-} Mass"}; 

  TH1D *spec[ptBins][djetBins][cTypes];
  TH1D *bgspec[ptBins][djetBins][cTypes];
  TH1D *sigspec[ptBins][djetBins][cTypes];

  TH1D *dr[ptBins][cTypes];
  TH1D *bgdr[ptBins][cTypes];

  TF1 *fits[ptBins][cTypes];
  TF1 *fitDraw[ptBins][cTypes];

  TFile *fin;
  if(!isMC) fin = new TFile("output/dCandidateToJetAssns_NoJetTrg_AllData_RevCuts_JetAssn_dR8-7-2.root");
  if(isMC==1) fin = new TFile("output/dCandidateToJetAssns_DpmEmbedMC_RevCuts_JetAssn-8-5.root");
  char* histoname = new char[100];
  for(int i=0; i<ptBins; i++){
    int j=0; //for(int j=0; j<djetBins; j++){
      for(int k=0; k<cTypes; k++){
	histoname = Form("assns_unweight_pt%d_djet%d_type%d",i,j,k);
	spec[i][j][k] = (TH1D*)fin->Get(histoname)->Clone(histoname);
	spec[i][j][k]->SetXTitle(xTitle[k].c_str());
	spec[i][j][k]->SetYTitle("Counts");
	if(doRebin) spec[i][j][k]->Rebin(rebinVal);

	histoname = Form("bgAssns_unweight_pt%d_djet%d_type%d",i,j,k);
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

	fits[i][k] = new TF1(Form("fitfg_pt%d_type%d",i,k),poly2reject,spec[i][j][k]->GetBinCenter(1),spec[i][j][k]->GetBinCenter(spec[i][j][k]->GetNbinsX()),3);
	fitDraw[i][k] = new TF1(Form("fitfgDraw_pt%d_type%d",i,k),"pol2",spec[i][j][k]->GetBinCenter(1),spec[i][j][k]->GetBinCenter(spec[i][j][k]->GetNbinsX()));
      }
      //}
  }
  
  double specScale[ptBins][cTypes];
  int nJets[ptBins][cTypes];

  TLatex *twrite[ptBins][djetBins][cTypes];
  TLatex *twrite2[cTypes][djetBins];
  TCanvas *cc[djetBins][cTypes];
  for(int k=0; k<cTypes; k++){
    globalCharmType = k;
    int j=0; //for(int j=0; j<djetBins; j++){
      cc[j][k] = new TCanvas(Form("cc_djet%d_type%d",j,k),Form("cc_djet%d_type%d",j,k),1000,800);
      if(drawMesonMass) cc[j][k]->Divide(2,3);
      else cc[j][k]->Divide(2,3);
      for(int i=0; i<ptBins; i++){
	//if(!drawMesonMass && i>=6) continue;
	if(i>=6) continue;
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
	double range2 = spec[i][j][k]->GetBinLowEdge(spec[i][j][k]->FindBin(tightSignalCutLow[k])+1);
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
	bgdr[i][k]->Scale(scaler3);
	bgdr[i][k]->SetLineColor(4);

	specScale[i][k] = scaler3;
	nJets[i][k] = region3Int;

	int tlbin = spec[i][j][k]->FindBin(tightSignalCutLow[k]);
	int thbin = spec[i][j][k]->FindBin(tightSignalCutHi[k]);
	//cout << "tlbin: "<< tlbin << " thbin: "<< thbin << endl;
	//cout << "spec: "<< spec[i][j][k]->Integral(tlbin,thbin) << endl;
	//cout << "type: "<< k << " pt: "<< i << " fg integral: "<< spec[i][j][k]->Integral(tlbin,thbin,"") << " bg: "<< region3Int/spec[i][j][k]->GetBinWidth(1) <<" signal to bg: " << (spec[i][j][k]->Integral(tlbin,thbin,"") - region3Int/spec[i][j][k]->GetBinWidth(1)) / (region3Int/spec[i][j][k]->GetBinWidth(1)) << endl;

	if(drawMesonMass){
	  sigspec[i][j][k] = (TH1D*)spec[i][j][k]->Clone(Form("signal_pt%d_type%d",i,k));
	  spec[i][j][k]->Draw();
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
	  dr[i][k]->GetXaxis()->SetRangeUser(0.,1.);
	  dr[i][k]->Draw();
	  bgdr[i][k]->Draw("same,hist");
	  sigspec[i][j][k]->Add(bgdr[i][k],-1);
	  sigspec[i][j][k]->SetLineColor(2);
	  sigspec[i][j][k]->Draw("hist,same");
	  
	}

	//cout << "S/B cross-check, fg: "<< dr[i][k]->Integral() << " bg: "<< bgdr[i][k]->Integral() << " s/b: "<< (dr[i][k]->Integral("") - bgdr[i][k]->Integral("")) / bgdr[i][k]->Integral("") << endl;
	//cout << "crosscheck2 " << bgdr[i][k]->Integral() << " " << region3Int/spec[i][j][k]->GetBinWidth(1) << endl;

	cc[j][k]->SetLogy();
	cc[j][k]->Update();

	if(i==0) {
	  if(!drawMesonMass) twrite2[i][j] = new TLatex(0.532,gPad->GetUymax()*0.85,Form("jet p_{T}: %d-%d, dType: %s",jetPtBoundaries[i],jetPtBoundaries[i+1], dtypes[k].c_str()));
	  else if(k==0) twrite2[i][j] = new TLatex(0.142,gPad->GetUymax()*0.85,Form("jet p_{T}: %d-%d, dType: %s",jetPtBoundaries[i],jetPtBoundaries[i+1], dtypes[k].c_str()));
	  else twrite2[i][j] = new TLatex(1.8,gPad->GetUymax()*0.85,Form("jet p_{T}: %d-%d, dType: %s",jetPtBoundaries[i],jetPtBoundaries[i+1], dtypes[k].c_str()));
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
  
}
