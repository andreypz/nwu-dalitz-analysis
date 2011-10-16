#define nC 20
#define F0 6
#define F1 6
#define F2 7
#define F3 4

struct optimalCuts {
  Float_t SB;
  Float_t SrootB;
  Float_t SrootSB;
  Int_t bin1;
  Int_t bin2;
};

void makePlot(Int_t sel=1, TString hPath="00")
{
  gROOT->ProcessLine(".L ./utils.C");
  gROOT->ProcessLine(".L ../data/tdrstyle.C");
  setTDRStyle();

  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLegendBorderSize(0);
  cout.precision(3); cout.setf(ios::fixed, ios::floatfield);
  TH1::SetDefaultSumw2(kTRUE);

  Float_t intLumi = 215.1 + 927.6 + 370.9 + 663.0 ; //double mu	

  Bool_t doTest = 0, doSB = 0, doEBEE=0;
  Bool_t doPhotons = 0, makeZjetsQt = 0;
  Bool_t doOverview= 0;

  //Types of met: met - pfMet, met1 - type1 corrected, met2 - pfMet passed Noise filters, 
  //met3 - projMet, met4 - puProj corrected met (those two are passed Noise filters) 
  TString metType("met2");   TString mtType("mt2"); 
  TString ssel("none");

  if (sel==1)  {ssel = "muon";     TString imgpath("~/afs/public_html/higgs/overview/muon/"); }
  if (sel==2)  {ssel = "electron"; TString imgpath("~/afs/public_html/higgs/overview/electron/"); }

  TString histoPath = Form("%s/%s", hPath.Data(), ssel.Data());
  cout<<"histoPath:  "<<histoPath.Data()<<"  int Lumi: "<<intLumi<<endl;

  if(sel==1){
    TFile* fda_2011A_DoubleMu_May10  = new TFile(Form("./%s/hhhh_DoubleMu_May10.root",histoPath.Data()));
    TFile* fda_2011A_DoubleMu_PromptV4  = new TFile(Form("./%s/hhhh_DoubleMu_PromptV4.root",histoPath.Data()));
  }
  TFile* fda_Data  = new TFile(Form("./%s/m_Data_%i.root", hPath.Data(),  sel));  //Merged Data
  TFile  *fData = (TFile*)fda_Data;
  
  //if (sel==1) TFile  *fData = (TFile*)fda_2011A_DoubleMu_May10;
  //if (sel==1) TFile  *fData = (TFile*)fda_2011A_DoubleMu_PromptV4;
  
  // TFile* fmc_ZllG      = new TFile(Form("./%s/hhhh_WZ.root", histoPath.Data() ));
  // TFile* fmc_Wjets     = new TFile(Form("./%s/hhhh_Wjets.root",histoPath.Data() ));

  TFile* fmc_ZZ        = new TFile(Form("./%s/hhhh_ZZ.root",histoPath.Data() ));
  TFile* fmc_WW        = new TFile(Form("./%s/hhhh_WW.root",histoPath.Data() ));
  TFile* fmc_WZ        = new TFile(Form("./%s/hhhh_WZ.root",histoPath.Data() ));

  TFile* fmc_ggH200    = new TFile(Form("./%s/hhhh_ggHZZ200.root",histoPath.Data() ));
  TFile* fmc_ggH250    = new TFile(Form("./%s/hhhh_ggHZZ250.root",histoPath.Data() ));
  TFile* fmc_ggH300    = new TFile(Form("./%s/hhhh_ggHZZ300.root",histoPath.Data() ));
  TFile* fmc_ggH350    = new TFile(Form("./%s/hhhh_ggHZZ350.root",histoPath.Data() ));
  TFile* fmc_ggH400    = new TFile(Form("./%s/hhhh_ggHZZ400.root",histoPath.Data() ));
  TFile* fmc_ggH450    = new TFile(Form("./%s/hhhh_ggHZZ450.root",histoPath.Data() ));
  TFile* fmc_ggH500    = new TFile(Form("./%s/hhhh_ggHZZ500.root",histoPath.Data() ));
  TFile* fmc_ggH550    = new TFile(Form("./%s/hhhh_ggHZZ550.root",histoPath.Data() ));

  //TFile* fmc_ggH600    = new TFile(Form("./%s/hhhh_SignalM600_HToZZ.root",histoPath.Data() ));

  if(doPhotons)  TFile* fmc_Zjets  = new TFile(Form("./%s/m_DataPh_%i.root", hPath.Data(), sel ));
  else           TFile* fmc_Zjets  = new TFile(Form("./%s/m_Zjets_%i.root", hPath.Data(), sel ));
  TFile* fmc_tt     = new TFile(Form("./%s/m_ttbar_%i.root", hPath.Data(), sel ));
  TFile* fmc_Top    = new TFile(Form("./%s/m_Top_%i.root", hPath.Data(), sel ));
 
  //List of background samples to Stack
  list_bg = new TList();
  list_bg->Add(fmc_Top);
  list_bg->Add(fmc_tt);
  list_bg->Add(fmc_WZ);
  list_bg->Add(fmc_WW);
  list_bg->Add(fmc_ZZ);
  list_bg->Add(fmc_Zjets);
  //list_bg->Add(fmc_Wjets);

  list_bg2 = new TList();
  list_bg2->Add(fmc_Top);
  list_bg2->Add(fmc_tt);
  list_bg2->Add(fmc_WZ);
  list_bg2->Add(fmc_WW);
  list_bg2->Add(fmc_ZZ);
  //list_bg2->Add(fmc_Wjets);

  list_overlay = new TList();
  list_overlay->Add(fData); 
  list_overlay->Add(fmc_ggH250);
  //list_overlay->Add(fmc_ggH300);
  list_overlay->Add(fmc_ggH400);

  list_signal = new TList();
  //list_signal->Add(fmc_ggH200);
  list_signal->Add(fmc_ggH250);
  list_signal->Add(fmc_ggH300);
  list_signal->Add(fmc_ggH350);
  list_signal->Add(fmc_ggH400);
  list_signal->Add(fmc_ggH450);
  list_signal->Add(fmc_ggH500);
  list_signal->Add(fmc_ggH550);

  PrintYields(list_bg, list_signal, fData, intLumi, histoPath, "twiki");
  PrintYields(list_bg, list_signal, fData, intLumi, histoPath, "tex");


  THStack *hs_met_et[nC], *hs_met2_et[nC], *hs_met_over_qt[nC], *hs_met2_over_qt[nC], *hs_di_qt[nC], *hs_met_et_ovQt[nC], *hs_mt[nC], *hs_mtZ[nC];
  THStack *hs_met0_et[nC], *hs_met1_et[nC], *hs_met2_et[nC], *hs_met3_et[nC], *hs_met4_et[nC], *hs_met5_et[nC], *hs_met6_et[nC], *hs_met7_et[nC], *hs_met8_et[nC], *hs_met9_et[nC], *hs_met10_et[nC];
  THStack *hs_jet_N[nC], *hs_jet_dRlep1[nC], *hs_jet_dRlep2[nC];
  THStack *hs_jet_b_N[nC], *hs_jet_b_Nssv[nC], *hs_jet_b_N25[nC], *hs_jet_b_N30[nC], hs_jet_b_pt[nC];
  THStack *hs_di_mass[nC], *hs_di_mass_EB[nC], *hs_di_mass_EE[nC], *hs_di_mass_EX[nC];
  THStack *hs_met_dPhiLeadJet1[nC], *hs_met_dPhiLeadJet2[nC], *hs_met_dPhiClosJet1[nC], *hs_met_dPhiClosJet2[nC];
  THStack *hs_vtx_nPV_raw[nC], *hs_vtx_nPV_weight[nC];
  //THStack *hs_vtx_ndof_1[nC], *hs_vtx_ndof_2[nC];

  for(Int_t n = 0; n<nC; n++)
    {
      hs_mt[n]           = makeStack(list_bg, Form("%s_%i", mtType.Data(), n), intLumi);

      hs_met0_et[n]      = makeStack(list_bg, Form("met0_et_%i", n), intLumi); //pfMet, no noise filters
      hs_met1_et[n]      = makeStack(list_bg, Form("met1_et_%i", n), intLumi); // pfMet type1
      hs_met2_et[n]      = makeStack(list_bg, Form("met2_et_%i", n), intLumi); //pfMet <--- default , no corrs
      hs_met3_et[n]      = makeStack(list_bg, Form("met3_et_%i", n), intLumi); //puCorr
      hs_met4_et[n]      = makeStack(list_bg, Form("met4_et_%i", n), intLumi); //proj
      hs_met5_et[n]      = makeStack(list_bg, Form("met5_et_%i", n), intLumi); //Zproj
      hs_met6_et[n]      = makeStack(list_bg, Form("met6_et_%i", n), intLumi); //red1
      hs_met7_et[n]      = makeStack(list_bg, Form("met7_et_%i", n), intLumi); //red2
      hs_met8_et[n]      = makeStack(list_bg, Form("met8_et_%i", n), intLumi); //comp (track, vtx asso)
      hs_met9_et[n]      = makeStack(list_bg, Form("met9_et_%i", n), intLumi); //comp (track, vtx asso)
      hs_met10_et[n]     = makeStack(list_bg, Form("met10_et_%i", n), intLumi); //comp (track, vtx asso)
      
      hs_di_qt[n]      = makeStack(list_bg, Form("di_qt_%i",n), intLumi);
      hs_di_mass[n]    = makeStack(list_bg, Form("di_mass_%i",n), intLumi);
      hs_di_mass_EB[n] = makeStack(list_bg, Form("di_mass_EB_%i",n), intLumi);
      hs_di_mass_EE[n] = makeStack(list_bg, Form("di_mass_EE_%i",n), intLumi);
      hs_di_mass_EX[n] = makeStack(list_bg, Form("di_mass_EX_%i",n), intLumi);
      hs_jet_N[n]      = makeStack(list_bg, Form("jet_N_%i",n), intLumi);
      hs_jet_b_N[n]    = makeStack(list_bg, Form("jet_b_N_%i",n), intLumi); 
      hs_jet_b_Nssv[n] = makeStack(list_bg, Form("jet_b_Nssv_%i",n), intLumi); 
      hs_jet_b_N25[n]  = makeStack(list_bg, Form("jet_b_N25_%i",n), intLumi); 
      hs_jet_b_N30[n]  = makeStack(list_bg, Form("jet_b_N30_%i",n), intLumi); 

      hs_vtx_nPV_raw[n]    = makeStack(list_bg, Form("vtx_nPV_raw_%i",n), intLumi); 
      hs_vtx_nPV_weight[n] = makeStack(list_bg, Form("vtx_nPV_weight_%i",n), intLumi); 
      //hs_vtx_ndof_1        = makeStack(list_bg, Form("vtx_ndof_1_%i",n), intLumi);
      //hs_vtx_ndof_2        = makeStack(list_bg, Form("vtx_ndof_2_%i",n), intLumi);

      hs_met_dPhiClosJet1[n] = makeStack(list_bg, Form("met2_dPhiClosJet1_%i",n), intLumi);
      hs_met_dPhiClosJet2[n] = makeStack(list_bg, Form("met2_dPhiClosJet2_%i",n), intLumi);

      /*hs_mtZ[n]          = makeStack(list_bg, Form("mtZ_%i", n), intLumi);

      hs_met_over_qt[n]  = makeStack(list_bg, Form("%s_over_qt_%i", metType.Data(), n), intLumi);
      hs_met2_over_qt[n] = makeStack(list_bg, Form("met2_over_qt_%i", n), intLumi);
      hs_met_et_ovQt[n]  = makeStack(list_bg, Form("%s_et_ovQt_%i", metType.Data(), n), intLumi);

      hs_jet_dRlep1[n]   = makeStack(list_bg, Form("jet_dRlep1_%i",n), intLumi);
      hs_jet_dRlep2[n]   = makeStack(list_bg, Form("jet_dRlep2_%i",n), intLumi);

      hs_met_dPhiLeadJet1[n] = makeStack(list_bg, Form("met2_dPhiLeadJet1_%i",n), intLumi);
      hs_met_dPhiLeadJet2[n] = makeStack(list_bg, Form("met2_dPhiLeadJet2_%i",n), intLumi);
      */
    }
    

  TH1 * forLegend[21];
  list_legend = new TList();
  list_legend -> AddAll(list_bg);
  list_legend -> AddAll(list_overlay);

  TFile *ff[10]; 
  Int_t size = list_legend->GetSize();
  TH1 *hh[10]; 
  if(size>10) {cout<<"To many plots to overlay"<<endl; return 0;}
  for(Int_t n=0; n<size; n++)
    {
      ff[n] = (TFile*)list_legend->At(n);
      hh[n] = (TH1*) ff[n]->Get("Andrey/met0_et_0")->Clone();
      cout<<n<<"   "<<ff[n] -> GetName()<<endl;
      forLegend[n] = (TH1F*)hh[n];
    }
 
  leg01 = new TLegend(0.53,0.7,0.95,0.95);
  leg01 -> SetNColumns(2);
  leg01 -> SetTextSize(0.04);
  leg01->AddEntry(forLegend[6], "Data","epl");
  leg01->AddEntry(forLegend[5],  "Z + jets","f");
  leg01->AddEntry(forLegend[4],  "ZZ","f");
  //leg01->AddEntry(forLegend[10],  "W + jets","f");
  leg01->AddEntry(forLegend[3],  "WW","f");
  leg01->AddEntry(forLegend[0],  "tW","f");
  leg01->AddEntry(forLegend[2],  "WZ","f");
  leg01->AddEntry(forLegend[1], "ttbar","f");  
  leg01->AddEntry(forLegend[7], "10xH250","f");
  //leg01->AddEntry(forLegend[8], "10xH400","f");
  //  leg01->AddEntry(forLegend[9], "10xH400","f");
  
  leg01->SetFillColor(kWhite);
  
  hs_vtx_nPV_raw[F0] -> Draw("hist");
  c1 -> SaveAs(imgpath+"ov01.png");
  TCanvas *c2 = new TCanvas("c2","example",600,700);
  

  if(doSB){
    cout<<"doing S/B"<<endl;

    TString testpath("~/afs/public_html/test/");  

    Double_t x_dPhi[200], x_met[200], x_mQt[200], x_err[200];
    Double_t StoB[200], StoB_err[200];
    Double_t sqrtStoB[200], sqrtStoB_err[200];
    
    gROOT->ProcessLine(".L ./xSecAndColors.C");

    Float_t nEv, scaleFactor1, scaleFactor2, cs;
    cs  = getXsecOrColors("DYmumu", 3);
    nEv = getXsecOrColors("DYmumu", 4);
    scaleFactor1      = intLumi*cs/nEv;
    cout<<cs<<"  "<<nEv<<"  "<<scaleFactor1<<endl;

    cs  = getXsecOrColors("ggHZZ400", 3);
    nEv = getXsecOrColors("ggHZZ400", 4);
    scaleFactor2      = intLumi*cs/nEv;
    TFile * sigFile = fmc_ggH400;
    Float_t  y1 = 4, y2 = 8;
    Int_t nDots    = 100;
    Int_t startpoint= 70;
    Float_t phiCut = 0.;

    TCanvas *c3 = new TCanvas("c3","for s/b plots",600,500);
    c3 -> cd();
    for(Int_t i =0; i<nDots; i++) x_err[i]=0;
    
    StoBPlot("pfMet",phiCut, fmc_Zjets,scaleFactor1,  sigFile,scaleFactor2,  StoB, StoB_err, sqrtStoB, sqrtStoB_err, x_met, startpoint);
    gr1 = new TGraphErrors(nDots, x_met, StoB, x_err, StoB_err);
    gr2 = new TGraphErrors(nDots, x_met, sqrtStoB, x_err, sqrtStoB_err);
    gr1 -> SetTitle(";pfMet cut; S/B");   
    twoScales(gr1, gr2, startpoint,startpoint+100, 0,y1,  y2, c3);
    c3 -> SaveAs(testpath+"p01.png");

    StoBPlot("pfMet1",phiCut, fmc_Zjets,scaleFactor1,  sigFile,scaleFactor2,  StoB, StoB_err, sqrtStoB, sqrtStoB_err, x_met, startpoint);
    gr1 = new TGraphErrors(nDots, x_met, StoB, x_err, StoB_err);
    gr2 = new TGraphErrors(nDots, x_met, sqrtStoB, x_err, sqrtStoB_err);
    gr1 -> SetTitle(";type 1 pfMet cut; S/B");   
    twoScales(gr1, gr2, startpoint,startpoint+100, 0,y1,  y2, c3);
    c3 -> SaveAs(testpath+"p02.png");
    
    StoBPlot("redMet2",phiCut, fmc_Zjets,scaleFactor1,  sigFile,scaleFactor2,  StoB, StoB_err, sqrtStoB, sqrtStoB_err, x_met, startpoint);
    gr1 = new TGraphErrors(nDots, x_met, StoB, x_err, StoB_err);
    gr2 = new TGraphErrors(nDots, x_met, sqrtStoB, x_err, sqrtStoB_err);
    gr1 -> SetTitle(";redMet2 cut; S/B");   
    twoScales(gr1, gr2, startpoint,startpoint+100, 0,y1,  y2, c3);
    c3 -> SaveAs(testpath+"p03.png");
    c3 -> SaveAs(testpath+"p04.png");

    /*
    StoBPlot("projMet",phiCut, fmc_Zjets,scaleFactor1,  sigFile,scaleFactor2,  StoB, StoB_err, sqrtStoB, sqrtStoB_err, x_met);
    gr1 = new TGraphErrors(nDots, x_met, StoB, x_err, StoB_err);
    gr2 = new TGraphErrors(nDots, x_met, sqrtStoB, x_err, sqrtStoB_err);
    gr1 -> SetTitle(";projMet cut; S/B");   
    twoScales(gr1, gr2, 0,200, 0,y1,  y2, c3);
    c3 -> SaveAs(testpath+"p03.png");

    StoBPlot("ZprojMet",phiCut, fmc_Zjets,scaleFactor1,  sigFile,scaleFactor2,  StoB, StoB_err, sqrtStoB, sqrtStoB_err, x_met);
    gr1 = new TGraphErrors(nDots, x_met, StoB, x_err, StoB_err);
    gr2 = new TGraphErrors(nDots, x_met, sqrtStoB, x_err, sqrtStoB_err);
    gr1 -> SetTitle(";ZprojMet cut; S/B");   
    twoScales(gr1, gr2, 0,200, 0,y1,  y2, c3);
    c3 -> SaveAs(testpath+"p04.png");
    

    StoBPlot("redMet1",phiCut, fmc_Zjets,scaleFactor1,  sigFile,scaleFactor2,  StoB, StoB_err, sqrtStoB, sqrtStoB_err, x_met);
    gr1 = new TGraphErrors(nDots, x_met, StoB, x_err, StoB_err);
    gr2 = new TGraphErrors(nDots, x_met, sqrtStoB, x_err, sqrtStoB_err);
    gr1 -> SetTitle(";redMet1 cut; S/B");   
    twoScales(gr1, gr2, 0,200, 0,y1,  y2, c3);
    c3 -> SaveAs(testpath+"p05.png");


    StoBPlot("puCorrMet",phiCut, fmc_Zjets,scaleFactor1,  sigFile,scaleFactor2,  StoB, StoB_err, sqrtStoB, sqrtStoB_err, x_met);
    gr1 = new TGraphErrors(nDots, x_met, StoB, x_err, StoB_err);
    gr2 = new TGraphErrors(nDots, x_met, sqrtStoB, x_err, sqrtStoB_err);
    gr1 -> SetTitle(";pu corr Met cut; S/B");   
    twoScales(gr1, gr2, 0,200, 0,y1,  y2, c3);
    c3 -> SaveAs(testpath+"p07.png");
    

    StoBPlot("compMet",phiCut, fmc_Zjets,scaleFactor1,  sigFile,scaleFactor2,  StoB, StoB_err, sqrtStoB, sqrtStoB_err, x_met);
    gr1 = new TGraphErrors(nDots, x_met, StoB, x_err, StoB_err);
    gr2 = new TGraphErrors(nDots, x_met, sqrtStoB, x_err, sqrtStoB_err);
    gr1 -> SetTitle(";compMet (vert assoc) cut; S/B");   
    twoScales(gr1, gr2, 0,200, 0, 1., 0.4, c3);
    c3 -> SaveAs(testpath+"p08.png");
    */

  /*
    fmc_Zjets->cd("Andrey");
    cutTree -> Draw("ct_compMet>>h1(80,0,400)","ct_evtWeight*(ct_dPhiMetJet>0.28 &&  ct_nJets>=2)","hist");
    sigFile->cd("Andrey");
    h1 -> Scale(scaleFactor1);
    cutTree -> Draw("ct_compMet>>h2(80,0,400)","ct_evtWeight*(ct_dPhiMetJet>0.28 &&  ct_nJets>=2)","hist same");
    h1 -> SetTitle("nJet >= 2; compMet (vert assoc); Events; ");
    h2 -> SetLineColor(kBlue);
    h2 -> Scale(scaleFactor2);
    c3 -> SetLogy();
    c3 -> SaveAs(testpath+"p05.png");
    c3 -> SetLogy(0);
    */
  }
  
  
  if(doTest){
    TString testpath("~/afs/public_html/test/");  
    

    //drawMuliPlot("#Delta#phi(MET, closest jet), p_{T}>20, |#eta|<2.4", 1, 0.001, 1000000, 0,5, hs_met_dPhiClosJet2[F0], c2, leg01, list_overlay, intLumi);
    //c2 -> SaveAs(testpath+"p03.png");

    drawMuliPlot("pfMET", 1, 0.001, 1000000, 0,3, hs_met2_et[6], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p01.png");
    drawMuliPlot("pfMET type1", 1, 0.001, 1000000, 0,3, hs_met1_et[6], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p02.png");
    drawMuliPlot("projMET", 1, 0.001, 1000000, 0,3, hs_met4_et[6], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p03.png");
    drawMuliPlot("ZprojMET", 1, 0.001, 1000000, 0,3, hs_met5_et[6], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p04.png");
    drawMuliPlot("redMET1", 1, 0.001, 1000000, 0,3, hs_met6_et[6], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p05.png");
    drawMuliPlot("redMET2", 1, 0.001, 1000000, 0,3, hs_met7_et[6], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p06.png");
    drawMuliPlot("puCorrMET", 1, 0.001, 1000000, 0,3, hs_met3_et[6], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p07.png");
    drawMuliPlot("compMET", 1, 0.001, 1000000, 0,3, hs_met8_et[6], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p08.png");
    //drawMuliPlot("pfMET", 1, 0.001, 1000000, 0,3, hs_met9_et[6], c2, leg01, list_overlay, intLumi); 
    //c2 -> SaveAs(testpath+"p09.png");
    //drawMuliPlot("pfMET", 1, 0.001, 1000000, 0,3, hs_met10_et[6], c2, leg01, list_overlay, intLumi); 
    //c2 -> SaveAs(testpath+"p10.png");

    /*
    drawMuliPlot("MT", 1, 0.001, 1000000, 0,3, hs_mt[6], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p02.png");

    cout<<"dbg"<<endl;

    drawMuliPlot("MT 250", 0, 0.001, 20, 0,5, hs_mt[7], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p03.png");
    drawMuliPlot("MT 300", 0, 0.001, 20, 0,5, hs_mt[8], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p04.png");

    drawMuliPlot("N b-jets", 0, 0.001, 30, 0.,3, hs_jet_b_N[15], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(testpath+"p05.png");
    drawMuliPlot("N b-jets ssv", 0, 0.001, 30, 0.,3, hs_jet_b_Nssv[15], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(testpath+"p06.png");

    drawMuliPlot("N b-jets 25", 0, 0.001, 30, 0.,3, hs_jet_b_N25[15], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(testpath+"p07.png");
    drawMuliPlot("N b-jets 30", 0, 0.001, 30, 0.,3, hs_jet_b_N30[15], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(testpath+"p08.png");
    //drawMuliPlot("N jets", 0, 0.001, 30, 0.,3, hs_jet_N[15], c2, leg01, list_overlay, intLumi);
    //c2 -> SaveAs(testpath+"p06.png");

    fData -> cd("Andrey");
    Float_t a = 0, aer=0;
    for(Int_t i=1; i<=4; i++){
      a = jet_b_N_15 -> GetBinContent(i);
      cout<<a<<endl;
    }
    fmc_tt -> cd("Andrey");
    jet_b_N_15 -> Scale(intLumi/1000);
    for(Int_t i=1; i<=4; i++){
      a   = jet_b_N_15 -> GetBinContent(i);
      aer = jet_b_N_15 -> GetBinError(i);
      cout<<a<<" +/- "<<aer<<endl;
    }

    fmc_Top -> cd("Andrey");
    jet_b_N_15 -> Scale(intLumi/1000);
    for(Int_t i=1; i<=4; i++){
      a   = jet_b_N_15 -> GetBinContent(i);
      aer = jet_b_N_15 -> GetBinError(i);
      cout<<a<<" +/- "<<aer<<endl;
    }
*/
    
    /*
    drawMuliPlot("pfMET", 1, 0.001, 1e6, 0,3, hs_met_et[6], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p07.png");
    drawMuliPlot("pfMET type1", 1, 0.001, 1e6, 0,3, hs_met1_et[6], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p08.png");

    //drawMuliPlot("pu corrMET", 1, 0.001, 1000000, 0,5, hs_met3_et[6], c2, leg01, list_overlay, intLumi); 
    // c2 -> SaveAs(testpath+"p08.png");

    drawMuliPlot("pfMET", 1, 0.001, 1000000, 0,5, hs_met_et[9], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p09.png");
    drawMuliPlot("pu corrMET", 1, 0.001, 1000000, 0,5, hs_met3_et[9], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p10.png");


    drawMuliPlot("nVtx raw", 1, 0.001, 1000000, 0,2, hs_vtx_nPV_raw[6], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p11.png");
    drawMuliPlot("nVtx reweighted", 1, 0.001, 1000000, 0,2, hs_vtx_nPV_weight[6], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p12.png");

    drawMuliPlot("nVtx raw", 0, 0.001, 30, 0,5, hs_vtx_nPV_raw[8], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p13.png");
    drawMuliPlot("nVtx reweighted", 0, 0.001, 30, 0,5, hs_vtx_nPV_weight[8], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p14.png");
    */

    //drawMuliPlot("N b-jets pt>25", 0, 0.0, 30, 0.,5, hs_jet_b_N25[7], c2, leg01, list_overlay, intLumi);
    //c2 -> SaveAs(testpath+"p13.png");
    //drawMuliPlot("N b-jets pt>30", 0, 0.0, 30, 0.,5, hs_jet_b_N30[7], c2, leg01, list_overlay, intLumi);
    //c2 -> SaveAs(testpath+"p14.png");

    // Mll in electons, EE vs EB
    if(doEBEE){
      drawMuliPlot("Leptons in Barrel, |#eta|<1.444,  M(ll)", 1, 0.001, 1000000, 0,2, hs_di_mass_EB[F0], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(testpath+"p11.png");
      drawMuliPlot("Leptons in Endcap, |#eta|>1.566, M(ll)", 1, 0.001, 1000000, 0,2, hs_di_mass_EE[F0], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(testpath+"p12.png");
      //drawMuliPlot("Leptons in EB/EE, mixed, M(ll)", 1, 0.001, 1000000, 0,4, hs_di_mass_EX[F0], c2, leg01, list_overlay, intLumi);
      //c2 -> SaveAs(testpath+"p03.png");
    }


    /*

    drawMuliPlot("#Delta#phi(MET, closest jet), p_{T}>20, |#eta|<5?", 1, 0.001, 1000000, 0,5, hs_met_dPhiClosJet1[6], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(testpath+"p05.png");
    drawMuliPlot("#Delta#phi(MET, closest jet), p_{T}>20, |#eta|<2.4", 1, 0.001, 1000000, 0,5, hs_met_dPhiClosJet2[6], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(testpath+"p06.png");
    
    //drawMuliPlot("projMET/q_{T}", 1, 0.001, 1000000, 0,3, hs_met_over_qt[5], c2, leg01, list_overlay, intLumi);
    //c2 -> SaveAs(testpath+"p03.png");
    
    c1 -> cd();
    TH2 *hh2 = (TH2*)fmc_ggH400 -> Get("met3_et_ovQt_6") ;
    hh2-> SetFillColor(kCyan);
    hh2 -> SetTitle(";projMET; projMET/q_{T}");
    hh2-> Draw("box1");
    hs_met_et_ovQt[6]->Sum()-> Clone("h_bkg1");
    h_bkg1->Draw("box1 same");
    c1 -> SaveAs(testpath+"p07.png");
    hh2->Delete();
    h_bkg1 -> Delete();

    c1 -> SetLogz(0);
    TH2 *hh2 = (TH2*)fmc_ggH400 -> Get("met3_et_ovQt_5") ;
    hh2-> SetFillColor(kCyan);
    hh2 -> SetTitle(";projMET; projMET/q_{T}");
    hh2-> Draw("box1");
    hs_met_et_ovQt[5]->Sum()-> Clone("h_bkg");
    h_bkg->SetLineColor(kOrange);
    h_bkg->Draw("box1 same");

    leg02 = new TLegend(0.75,0.8,0.9,0.89);
    leg02 -> SetTextSize(0.04);
    leg02->AddEntry(hh2, "ggH400", "f");
    leg02->AddEntry(h_bkg, "All bkg", "f");
    leg02->SetFillColor(kWhite);
    leg02 -> Draw();
    c1 -> SaveAs(testpath+"p09.png");



    TCutG *mg1 = new TCutG("mg1",5);
    mg1-> SetPoint(0,110,0.5);
    mg1-> SetPoint(1,110,0.9);
    mg1-> SetPoint(2,400,2.5);
    mg1-> SetPoint(3,400,1.9);
    mg1-> SetPoint(4,110,0.5);


    Double_t mg1_iBkg = mg1->IntegralHist(h_bkg);     
    Double_t mg1_iSig = mg1->IntegralHist(hh2);     

    Double_t iBkg = h_bkg -> Integral(12,41, 0,41);
    Double_t iSig = hh2   -> Integral(12,41, 0,41);
    cout<<endl;
    cout<<"at 110: S/B = "<<0.1*iSig<<"/"<<iBkg<<" = "<<0.1*iSig/iBkg<<endl;
    cout<<"    2D, S/B = "<<0.1*mg1_iSig<<"/"<<mg1_iBkg<<" = "<<0.1*mg1_iSig/mg1_iBkg<<endl;
    cout<<endl;

    hh2->Draw(Form("[%s] box1", mg1->GetName()));
    h_bkg->Draw(Form("[%s] box1 same", mg1->GetName()));
   
    leg02-> Draw();
    c1 -> SaveAs(testpath+"p10.png");


    TCutG *mg2 = new TCutG("mg2",5);
    mg2-> SetPoint(0,110,0);
    mg2-> SetPoint(1,110,4);
    mg2-> SetPoint(2,400,4);
    mg2-> SetPoint(3,400,0);
    mg2-> SetPoint(4,100,0);
    hh2->Draw(Form("[%s] box1", mg2->GetName()));
    h_bkg->Draw(Form("[%s] box1 same", mg2->GetName()));

    leg02-> Draw();
    c1 -> SaveAs(testpath+"p11.png");
    */

    /*
    TF1 * f1 = new TF1 ("f1","[0] + [1]*x +[2]*x*x",0,400);
    f1->SetParameter(0, 1.108);
    f1->SetParameter(1, 0.002429);
    f1->SetParameter(2, -1.655e-06);
    f1 -> Draw();
    f1 -> SetMinimum(0);
    f1 -> SetMaximum(2.8);
    
    c1 -> SaveAs(testpath+"p06.png");
    */
    
    //    drawMuliPlot("N jets", 1, 0.001, 1000000, 0,5, hs_jet_N[F0], c2, leg01, list_overlay, intLumi);
    // c2 -> SaveAs(testpath+"p09.png");
    
    /*
    if(hPath=="v18"){
      drawMuliPlot("#DeltaR(jet, lepton1)", 1, 0.001, 1000000, 0,5, hs_jet_dRlep1[4], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(testpath+"p09.png");
      drawMuliPlot("#DeltaR(jet, lepton2)", 1, 0.001, 1000000, 0,5, hs_jet_dRlep2[4], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(testpath+"p10.png");
    }
    drawMuliPlot("projMET/q_{T}", 1, 0.001, 1000000, 0,3, ph_met_over_qt[5], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(testpath+"p13.png");

    drawMuliPlot("projMET/q_{T}", 1, 0.001, 1000000, 0,3, hs_met_over_qt[5], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(testpath+"p14.png");
    

    drawMuliPlot("MT", 1, 0.001, 1000000, 0,3, hs_mt[7], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(testpath+"p15.png");
    drawMuliPlot("MT", 1, 0.001, 1000000, 0,3, ph_mt[7], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(testpath+"p16.png");

    drawMuliPlot("MT", 0, 0.001, 10, 0,3, hs_mt[6], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(testpath+"p17.png");
    drawMuliPlot("MT", 0, 0.001, 10, 0,3, ph_mt[6], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(testpath+"p18.png");

    */
  }
  
  if(doOverview){

      drawMuliPlot("pfMET", 1, 0.001, 1000000, 0,3, hs_met2_et[F0], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov01.png");
      drawMuliPlot("projMET", 1, 0.001, 1000000, 0,3, hs_met4_et[F0], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov02.png");

      drawMuliPlot("MT", 1, 0.001, 1000000, 0,2, hs_mt[F0], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov03.png");
    
      drawMuliPlot("q_{T} (di-lepton p_{T})", 1, 0.001, 1000000, 0,3, hs_di_qt[F0], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov04.png");
      
      drawMuliPlot("M(ll)", 1, 0.001, 1000000, 0,2, hs_di_mass[F0], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov05.png");
      drawMuliPlot("N jets", 1, 0.001, 1000000, 0,5, hs_jet_N[F0], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov06.png");

      drawMuliPlot("N b-jets", 1, 0.001, 1000000, 0,5, hs_jet_b_N[F0], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov07.png");
      drawMuliPlot("N b-jets (pt>30)", 1, 0.001, 1000000, 0,5, hs_jet_b_N30[F0], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov08.png");
      
      drawMuliPlot("#Delta#phi(MET, closest jet), p_{T}>20, |#eta|<2.4", 1, 0.001, 10000000, 0,2, hs_met_dPhiClosJet1[F0], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov09.png");

    /*
    drawMuliPlot("projMET/q_{T}", 1, 0.001, 1000000, 0,3, hs_met_over_qt[F0], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(imgpath+"ov02.png");
    
    drawMuliPlot("projMET", 1, 0.001, 1000000, 0,5, hs_met_et[F1], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(imgpath+"ov03.png");
    drawMuliPlot("projMET/q_{T}", 1, 0.001, 1000000, 0,5, hs_met_over_qt[F1], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(imgpath+"ov04.png");
  

    drawMuliPlot("projMET", 1, 0.001, 1000000, 0,3, hs_met_et[F2], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(imgpath+"ov05.png");
    drawMuliPlot("projMET/q_{T}", 1, 0.001, 1000000, 0,3, hs_met_over_qt[F2], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(imgpath+"ov06.png");
    */
    
      /*    
      drawMuliPlot("M(ll)", 1, 0.001, 10000000, 0,2, hs_di_mass[F2], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov11.png");
      drawMuliPlot("MT", 1, 0.001, 1000000, 0,2, hs_mt[F2], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov12.png");
      
      
      drawMuliPlot("N jets", 1, 0.001, 1000000, 0,5, hs_jet_N[F1], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov15.png");
      drawMuliPlot("N b-jets", 1, 0.001, 1000000, 0,5, hs_jet_b_N[F1], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov16.png");
      
      
      drawMuliPlot("N jets", 1, 0.001, 1000000, 0,5, hs_jet_N[F3], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov17.png");
      drawMuliPlot("N b-jets", 1, 0.001, 1000000, 0,5, hs_jet_b_N[F3], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov18.png");
      
      
      // drawMuliPlot("N jets", 1, 0.001, 1000000, 0,5, hs_jet_N[F2], c2, leg01, list_overlay, intLumi);
      // c2 -> SaveAs(imgpath+"ov17.png");
      // drawMuliPlot("N b-jets", 1, 0.001, 1000000, 0,5, hs_jet_b_N[F2], c2, leg01, list_overlay, intLumi);
      //   c2 -> SaveAs(imgpath+"ov18.png");
      
      
      drawMuliPlot("N b-jets", 1, 0.001, 1000000, 0,5, hs_jet_b_N[7], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov19.png");
      drawMuliPlot("N b-jets", 1, 0.001, 1000000, 0,5, hs_jet_b_N[8], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov20.png");
      
      drawMuliPlot("#Delta#phi(MET, lead jet), p_{T}>20, |#eta|<2.4", 1, 0.001, 10000000, 0,2, hs_met_dPhiLeadJet1[F0], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov21.png");
      
      drawMuliPlot("#Delta#phi(MET, lead jet), p_{T}>20, |#eta|<2.4", 1, 0.001, 1000000, 0,5, hs_met_dPhiLeadJet1[F1], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov23.png");
      drawMuliPlot("#Delta#phi(MET, closest jet), p_{T}>20, |#eta|<2.4", 1, 0.001, 1000000, 0,5, hs_met_dPhiClosJet1[F1], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov24.png");
      
      drawMuliPlot("#Delta#phi(MET, lead jet), p_{T}>20, |#eta|<2.4", 1, 0.001, 10000000, 0,2, hs_met_dPhiLeadJet1[F2], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov25.png");
      drawMuliPlot("#Delta#phi(MET, closest jet), p_{T}>20, |#eta|<2.4", 1, 0.001, 10000000, 0,2, hs_met_dPhiClosJet1[F2], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov26.png");
      
      
      drawMuliPlot("q_{T} (di-lepton p_{T})", 1, 0.001, 1000000, 0,5, hs_di_qt[F0], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov27.png");
      drawMuliPlot("q_{T} (di_lepton p_{T})", 1, 0.001, 1000000, 0,5, hs_di_qt[F1], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov28.png");
      
      drawMuliPlot("pfMET", 1, 0.001, 1000000, 0,5, hs_met2_et[F0], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov29.png");
      drawMuliPlot("MT (using Z mass)", 1, 0.001, 1000000, 0,5, hs_mtZ[F0], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov30.png");
      
      drawMuliPlot("pfMET", 1, 0.001, 1000000, 0,5, hs_met2_et[F1], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov31.png");
      drawMuliPlot("MT (using Z mass)", 1, 0.001, 1000000, 0,5, hs_mtZ[F1], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov32.png");
  */
  
  }

  cout<<"\n\n   -- end of job  --"<<endl;
  
}
/*
THStack* makeStack(TList *sourcelist, TString name, Float_t lumi)
{
  //cout<<"making stack!"<<endl;
  THStack * hs = new THStack(Form("hs_11_%i",1),"Stacked hist");
  //sourcelist ->Print();

  TFile *first_source = (TFile*)sourcelist->First();
  first_source->cd("Andrey");
  TDirectory *current_sourcedir = gDirectory;

  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key, *oldkey=0;
  while ( (key = (TKey*)nextkey())) 
    {
      if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;
      //first_source->cd();
 
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1")) continue;
	
	if(strncmp(Form("%s",name.Data()),key->GetName(), name.Length()) ==0)
	  {
	    //cout<<"first "<<key->GetName()<<endl;
	    //hs -> Add((TH1*)key->ReadObj());
	    TH1 *hh1 = (TH1*)key->ReadObj()->Clone();
	    hh1 -> Scale(lumi/1000);
	    hs  -> Add(hh1);
	    //hh1 -> Print();
	    //hs -> Add((TH1*)key->ReadObj()->Clone());
	    //hs -> Add( ((TH1*)key->ReadObj()->Clone())->Scale(lumi/1000));
	  }
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      while ( nextsource )
	{
	  nextsource->cd("Andrey");
	  TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(key->GetName());
	  
	  if (key2) 
	    {
	      //TH1 *h2 = (TH1*)key2->ReadObj();
	      if(strncmp(Form("%s",name.Data()),key2->GetName(), name.Length()) ==0)
		{
		  //cout<<"nextsource "<<key2->GetName()<<endl;
		  TH1 *hh2 = (TH1*)key2->ReadObj()->Clone();
		  hh2 -> Scale(lumi/1000);
		  hs -> Add(hh2);
		  //hs -> Add((TH1*)key2->ReadObj()->Clone());
		  //hs -> Add( ((TH1*)key2->ReadObj()->Clone())->Scale(lumi/1000));
		}
	      //delete h2;
	    }
	  
	  nextsource = (TFile*)sourcelist->After( nextsource );
	}
      
    }



  return hs;
}
*/
 /*
void drawMuliPlot(TString xtitle, Int_t isLog, Float_t y1min, Float_t y1max, Float_t y2min, Float_t y2max, THStack *hs, TCanvas *cc, TLegend *leg, TList *list, Float_t lumi)
{
  //Find the samm objects as in hs stack
  TString name = Form("Andrey/%s", hs->GetHists()->First()->GetName());
  //data -> ls();
  cout<<"drawMultiPlot:: " <<name<<endl;

  TFile *ff[5];   Int_t size = list->GetSize();
  TH1 *hh[5]; 
  if(size>5) {cout<<"To many plots to overlay"<<endl; return 0;}
  for(Int_t n=0; n<size; n++)
    {
      ff[n] = (TFile*)list->At(n);
      hh[n] = (TH1*)ff[n]->Get( name.Data() )->Clone();
      //cout<<n<<" file  "<<ff[n] -> GetName()<<"   histoname: "<<hh[n]->GetName()<<endl;
    }


  for(Int_t n=1; n<size; n++)
    hh[n] ->Scale(lumi/1000./10);

  TH1 *h_data = hh[0]->Clone();

  //h_data -> Print();

   cc ->cd();
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  cc -> cd();
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1 -> cd();
  pad1 -> SetLogy(isLog);

  hs -> Draw("hist");
  hs -> SetTitle(";; Events");
  hs -> SetMinimum(y1min);
  hs -> SetMaximum(y1max);
  h_data  -> Draw("same e1pl"); //DATA
 for(Int_t n=1; n<size; n++)
    hh[n] -> Draw("same hist");
  leg -> Draw();

  cc ->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetBottomMargin(0.25);
  pad2->SetTopMargin(0);
  pad2->Draw();
  pad2 -> cd();
  
  TH1F* h1 = (TH1F*) h_data->Clone("");
  TH1F* h2 = (TH1F*) hs->Sum();
  h1 -> Divide(h2);
  h1 -> SetTitle(Form(";%s; Data/MC",xtitle.Data()));
  h1 -> SetMaximum(y2max);
  h1 -> SetMinimum(y2min);
  h1->GetYaxis()->SetNdivisions(206);
  h1->GetYaxis()->SetTitleOffset(0.4);
  h1->SetTitleSize(0.1,"XYZ");
  h1 ->SetLabelSize(0.1,"XY");
  h1 -> Draw("ep");

  pad1->cd();
  TLatex *prelim;
  prelim = new TLatex(0.30,0.95, Form("CMS Preliminary       #it{L_{int}} = %5.f pb^{-1}",lumi));
  prelim -> SetNDC(); prelim->SetTextSize(0.03); 
  prelim -> Draw();

}
 */
  /*
TH1 *THStack::Sum()
{
  //cout<<"doing sum  " <<endl;
  TList * mylist = (TList*)this->GetHists();
  TIter next(mylist);
  TH1 *hh  = (TH1*) mylist -> First() ->Clone();
  hh -> SetLineColor(kBlack);
  hh -> SetFillStyle(0);
  TObject *obj; 
  while ((obj = next()))
    {
      // cout<<obj->GetName()<<endl;
      //skip first object since it's used by creating the histogram
      if(obj == mylist->First()) continue;
      hh -> Add((TH1*)obj);     
    }
  //cout<<"end of sum"<<endl;
  return hh;
}
*/

/*
void PrintYields(TList *bgList, TList *sigList, TFile *dataFile, Float_t lumi, TString path, string option="tex")
{
  Int_t bins = 0;
  ofstream oo;
  oo.open(Form("%s/yields_%s.txt",path.Data(), option.c_str()), ofstream::out);
  oo.precision(1); oo.setf(ios::fixed, ios::floatfield);

  TString beginLine("");
  TString endLine("");
  TString separator("");
  TString title1("");
  TString title2("");
  TString pmSign("");
  if(option=="tex"){
    title1 = " sel & WW & WZ   \\\\ \\hline";
    title2 = " sel & WW & WZ   \\\\ \\hline";
    beginLine = " ";
    endLine   = "\\\\ \\hline";
    separator = "\t &";
    pmSign = " $\pm$ ";
  }
  if(option=="twiki"){
    title1 = "| *cut*        | *top*  | *ttbar*  | *WZ*  | *WW* | *ZZ*  | *Zjets*  | *Data*  | *Total bg* | | | |";
    title2 = "| *Higgs mass* | *top*  | *ttbar*  | *WZ*  | *WW* | *ZZ*  | *Zjets*  | *Data*  | *Total bg* | *higgs* | *S/B* |";
    beginLine = "| ";
    endLine   = "\t |";
    separator = "\t |";
    pmSign = " &plusmn; ";
  }

  string cutNames[18] = {
    "0. Total",
    "1. trigger",
    "2. Vtx, cosmic, 2 lept",
    "3. Z mass",
    "4. soft 3d muon veto",
    "5. qT > 25",
    "6. b-veto",
    "7.  H200",
    "8.  H250",
    "9.  H300",
    "10. H350",
    "11. H400",
    "12. H450",
    "13. H500",
    "14. H550",
    "15. H300 filt",
    "16.",
    "17.",};

  oo<<title1<<endl;

  for(Int_t j = 0; j<=14; j++)
    {
      if(j<=6)  oo.precision(0);
      else  oo.precision(2);
      if(j==7) continue; //skip 200 mass
      if(j==8) oo<<title2<<endl;

      Int_t size = bgList->GetSize();
      TFile *ff[10];    TH1 *hh[10];
      if(size>10) {cout<<"Yields: to many MC samples to print-out"<<endl; return 0;}

      oo<<beginLine<<cutNames[j]<<separator;

      Float_t total_bg = 0, total_bgError =0;

      for(Int_t n=0; n<size; n++)
	{
	  ff[n] = (TFile*)bgList->At(n);
	  if(j<6) hh[n] = (TH1*)ff[n]->Get( Form("Andrey/met0_et_%i",j) )->Clone();
	  else hh[n] = (TH1*)ff[n]->Get( Form("Andrey/met0_et_%i",j) )->Clone();
	  //cout<<n<<"Yields:: file  "<<ff[n] -> GetName()<<"   histoname: "<<hh[n]->GetName()<<endl;
	
	  hh[n] -> Scale(lumi/1000);
	  bins = hh[n] -> GetNbinsX();
	  //Float_t val = hh[n]->Integral(0,bins+1);
	  Double_t err;
	  Double_t val = hh[n]->TH1::IntegralAndError(0,bins+1, err);
	  if(j<=6) oo<<val<<separator;
	  else oo<<val<<pmSign<<err<<separator;
	  total_bg += val;
	  total_bgError += err;
	}      


      if(j<=7){
	TH1* data = (TH1*)dataFile->Get(  Form("Andrey/met0_et_%i",j) )->Clone();
	Int_t nBins = data->GetNbinsX();
	Float_t iBkg = total_bg;
	oo.precision(0);
	oo<<data->Integral(0,nBins+1);

	oo<<separator<<iBkg<<separator<<separator<<endLine<<endl;
      }
      else{
	TH1* data = (TH1*)dataFile->Get(  Form("Andrey/met0_et_%i",j) )->Clone();

	TFile* sigFile  = (TFile*)sigList->At(j-8); //higgs 
	if(j==15)  
	  TFile* sigFile  = (TFile*)sigList->At(1);
       
	TH1* sig  = (TH1*)sigFile ->Get(  Form("Andrey/met0_et_%i",j))->Clone();
	sig -> Scale(lumi/1000);
	Int_t sig_nBins = sig->GetNbinsX();
	Int_t bkg_nBins = hh[0]->GetNbinsX();
	Int_t nBins=0;
	if (sig_nBins!=bkg_nBins) cout<<" IN yields \n WARNING:  different number of bis"<<endl;
	else nBins = sig_nBins;
	
	//Float_t iSig = sig -> Integral();  //Count overflows
	Float_t iSig = sig -> Integral(0,nBins+1);  //Count overflows
	Float_t iBkg = total_bg;
	Float_t iBkg_err = total_bgError;
	Float_t SB   = 0.1*iSig/iBkg;        //Higgs signal is 10*real in the histograms
	//Float_t SrootB  = 0.1*iSig/sqrt(iBkg);
	//Float_t SrootSB = 0.1*iSig/sqrt(iBkg + 0.1*iSig);

	oo.precision(0);
	oo<<data->Integral(0,nBins+1);

	oo.precision(2);
	oo<<separator<<iBkg<<pmSign<<iBkg_err<<separator<<0.1*iSig<<separator;
	oo.precision(3);
	oo<<SB<<endLine<<endl;  
      }      
  }

  oo<<"end of printing yields"<<endl;
  oo.close();
  
}
*/


 
void AddPhoton(TH1 *ph, THStack *st, Float_t lumiToScale, Float_t lumiUsed){
  ph -> SetLineColor(kRed);
  ph -> SetLineWidth(2);
  ph -> Scale(lumiToScale/lumiUsed);
  st  -> Add(ph);
}

optimalCuts calculateOptimalCuts(Int_t opt = 1, TH1 *bkg, TH1 *sig) {
  //Option values opt can be 1,2,3
  // opt = 1 - optimize S/B
  // opt = 2 - optimize S/sqrt(B)
  // opt = 3 - optimize S/sqrt(S+B)

  cout<<  bkg -> GetName()<<endl;

  Int_t sig_nBins = sig->GetNbinsX();
  Int_t bkg_nBins = bkg->GetNbinsX();
  Int_t nBins=0;
  if (sig_nBins!=bkg_nBins) cout<<" IN calculateOptimalCuts \n WARNING:  different number of bis"<<endl;
  else nBins = sig_nBins;
  
  Float_t iSig = sig -> Integral(0,nBins+1);  //Count overflows
  Float_t iBkg = bkg -> Integral(0,nBins+1);
  Float_t SB      = 0.1*iSig/iBkg;        //Higgs signal is 10*real in the histograms
  Float_t SrootB  = 0.1*iSig/sqrt(iBkg);
  Float_t SrootSB = 0.1*iSig/sqrt(iBkg + 0.1*iSig);

  optimalCuts answer;
  answer.SB      = SB;
  answer.SrootB  = SrootB;
  answer.SrootSB = SrootSB;
  answer.bin1 = 0;
  answer.bin2 = nBins+1;
  cout<<"Initial values of S/B = "<<answer.SB<<"   S/sqrt(B) = "<<SrootB<<"   S/sqrt(B+S) = "<<SrootSB<<endl;

  for(Int_t a=0; a<nBins;a++)
    for(Int_t b=nBins+1; b>a+5; b--)
      {      
	iSig = sig -> Integral(a,b);
	iBkg = bkg -> Integral(a,b);
	SB      = 0.1*iSig/iBkg;        //Higgs signal is 10*real in the histograms
	SrootB  = 0.1*iSig/sqrt(iBkg);
	SrootSB = 0.1*iSig/sqrt(iBkg + 0.1*iSig);
	if(opt==1)
	  if (SB>answer.SB && SB<100){ //<100 to prevent infinity (zero BG)
	    answer.SB = SB;
	    answer.SrootB = SrootB;
	    answer.SrootSB = SrootSB;
	    answer.bin1 = a;
	    answer.bin2 = b;
	  }
	if(opt==2)
	  if (SrootB>answer.SrootB && SB<100){
	    answer.SB = SB;
	    answer.SrootB = SrootB;
	    answer.SrootSB = SrootSB;
	    answer.bin1 = a;
	    answer.bin2 = b;
	  }
	if(opt==3)
	  if (SrootSB>answer.SrootSB && SB<100){
	    answer.SB = SB;
	    answer.SrootB = SrootB;
	    answer.SrootSB = SrootSB;
	    answer.bin1 = a;
	    answer.bin2 = b;
	  }
	

      }
  
  return answer;
}




void twoScales(TGraph *g1, TGraph *g2, Float_t x1, Float_t x2, Float_t y1, Float_t y2, Float_t y3, TCanvas *c1) {
  c1 -> Clear();

    g1->SetFillColor(kBlue);
    g1->SetFillStyle(3010);
    g2->SetFillColor(kRed);
    g2->SetFillStyle(3010);


  //compute the pad range with suitable margins
  Double_t ymin = y1;
  Double_t ymax = y2;
  Double_t dy = (ymax-ymin)/0.76; //10 per cent margins top and bottom
  Double_t xmin = x1;
  Double_t xmax = x2;
  Double_t dx = (xmax-xmin)/0.76; //10 per cent margins left and right

  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  TPad *pad2 = new TPad("pad2","",0,0,1,1);
  pad2->SetFillStyle(4000); //will be transparent
  
  pad1->Draw();
  pad1->cd();
  // Margins:
  pad1->SetTopMargin(0.1);
  pad1->SetBottomMargin(0.14);
  pad1->SetLeftMargin(0.14);
  pad1->SetRightMargin(0.1);
  pad1->SetTicky(0);

  g1->Draw("AL3");
  g1-> SetMaximum(y3);
  g1-> SetMinimum(0);

  g1->GetHistogram()->SetAxisRange(xmin,xmax);
  g1 -> GetYaxis() -> SetLabelColor(kBlue);
  g1 -> GetYaxis() -> SetTitle("S/B");
  g1 -> GetYaxis() -> SetTitleOffset(1.0);
  g1 -> GetYaxis() -> SetTitleSize(0.05);
  pad1->Modified();
  c1->cd();

  pad2->Range(xmin-0.14*dx,ymin-0.14*dy,xmax+0.1*dx,ymax+0.1*dy);  //correspond to Margins in gstyle
  pad2->Draw();
  pad2->cd();
  
  g2->Draw("same L3");
   
  // draw axis on the right side of the pad
  TGaxis *axis1 = new TGaxis(x2,y1,x2,y2,y1,y2,510,"+L");
  axis1->SetLabelColor(kRed);
  axis1 -> SetTitle("S/#sqrt{S+B}");
  axis1->Draw();

  //TGaxis *axis2 = new TGaxis(x1,y1,x1,y2,y1,y2,510,"+L");
  //axis2->SetLabelColor(kRed);
  //axis2 -> SetTitle("S/#sqrt{S+B}");
  //axis2->Draw();

  pad2->Update();

  pad1 -> cd();
  leg03 = new TLegend(0.2,0.73,0.4,0.89);
  leg03 -> SetTextSize(0.04);
  leg03->AddEntry(g1, "S/B", "f");
  leg03->AddEntry(g2, "#frac{S}{#sqrt{S+B}}", "f");
  leg03->SetFillColor(kWhite);
  leg03 -> Draw();
  	

}

StoBPlot(string var= "pfMet", Float_t phiCut = 0, TFile *bgFile, Float_t sc1, TFile *sigFile, Float_t sc2,  Double_t StoB[], Double_t StoB_err[], Double_t sqrtStoB[],Double_t  sqrtStoB_err[], Double_t x_met[], Int_t startpoint){

  /*  bgFile -> cd("Andrey");
  TString ct_sample("a");

  TBranch *branch  = cutTree->GetBranch("ct_sample");
  branch->SetAddress(&ct_sample);
  cutTree -> GetEntry(0);
  cout<<"test:  "<<ct_sample<<endl;
  */

  Float_t dPhiCut = phiCut; //400 etc
  
  for(Int_t i=0; i<100; i++)
    {
      x_met[i]  = startpoint+ 1*i;

      //cout<<i<<"  "<<intLumi<<endl;
      bgFile -> cd("Andrey");
      //cutTree -> Draw(Form("ct_%s>>h1",var.c_str()),Form("ct_evtWeight*(ct_dPhiMetJet>%f && ct_MT>216 &&ct_MT<272 && ct_%s>%f)", dPhiCut, var.c_str(), x_met[i]),"hist");
      cutTree -> Draw(Form("ct_%s>>h1",var.c_str()),Form("ct_evtWeight*(ct_dPhiMetJet>%f && ct_%s>%f)", dPhiCut, var.c_str(), x_met[i]),"hist");
      //cutTree -> Draw(Form("ct_%s>>h1",var.c_str()),Form("ct_evtWeight*(ct_dPhiMetJet>%f && ct_%s>%f && ct_nJets>=2)", dPhiCut, var.c_str(), x_met[i]),"hist");
      
      TH1F * temp = (TH1F*)h1;
      temp->Scale(sc1);
      Double_t bg_err =0 ;
      //temp -> Print();
      Double_t bg =  temp->TH1::IntegralAndError(-1,2000, bg_err);
      
      sigFile -> cd("Andrey");
      //cutTree -> Draw(Form("ct_%s>>h2",var.c_str()),Form("ct_evtWeight*(ct_dPhiMetJet>%f && ct_MT>216 &&ct_MT<272 && ct_%s>%f)", dPhiCut, var.c_str(), x_met[i]),"hist");
      cutTree -> Draw(Form("ct_%s>>h2",var.c_str()),Form("ct_evtWeight*(ct_dPhiMetJet>%f && ct_%s>%f)", dPhiCut, var.c_str(), x_met[i]),"hist");
//cutTree -> Draw(Form("ct_%s>>h2",var.c_str()),Form("ct_evtWeight*(ct_dPhiMetJet>%f && ct_%s>%f && ct_nJets>=2)", dPhiCut, var.c_str(), x_met[i]),"hist");

            
      TH1F * temp = (TH1F*)h2;
      temp->Scale(sc2);
      Double_t sig_err = 0;
      //Double_t sig = 	temp -> Integral();
      Double_t sig = temp->TH1::IntegralAndError(-1,2000, sig_err);
            
      if(bg!=0) {
	StoB[i]    = sig/bg;
	StoB_err[i] = StoB[i]*sqrt( pow(sig_err/sig,2) + pow(bg_err/bg,2));
      }
      else {StoB[i]=0; StoB_err[i]=0;}
      //x_err[i] =0;
      
      if((sig+bg)  != 0) {
	sqrtStoB[i]     = sig/sqrt(sig+bg);
	sqrtStoB_err[i] =  sqrtStoB[i]*sqrt( pow(sig_err/sig,2) + 0.25* (pow(sig_err,2) + pow(bg_err,2))/pow(sig+bg,2)  );
      }
      else {sqrtStoB[i] = 0; sqrtStoB_err[i] =0;}
      
      //is this correct error propagation?	
      
      cout<<i<<"  "<<x_met[i]<<"    B="<<bg<<" S="<<sig<<"  S/B="<<StoB[i]<<"   S/sq(S+B)="<<sqrtStoB[i]<<endl;
      
    }
}

