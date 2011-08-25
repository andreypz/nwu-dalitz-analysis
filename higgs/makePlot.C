#define nC 20
#define F0 5
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
  gROOT->ProcessLine(".L ../data/tdrstyle.C");
  setTDRStyle();

  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLegendBorderSize(0);
  cout.precision(3); cout.setf(ios::fixed, ios::floatfield);
  TH1::SetDefaultSumw2(kTRUE);

  //Float_t intLumi = 191.; //Note: in v11 and below the histograms are already normylized to 191, after that - to 1000
  // Float_t intLumi = 191. + 963.6 *0.94; //-6% of lumi 191 =0.94* 204
  Float_t intLumi = 928.2;
  //Float_t intLumi = 201.2 + 928.2;
 

  //Types of met: met - pfMet, met1 - type1 corrected, met2 - pfMet passed Noise filters, 
  //met3 - projMet, met4 - puProj corrected met (those two are passed Noise filters) 
  TString metType("met2");   TString mtType("mt2"); 
  TString ssel("none");

  if (sel==1)  {ssel = "muon"; TString imgpath("~/afs/public_html/higgs/overview/muon/"); }
  if (sel==2)  {ssel = "electron"; TString imgpath("~/afs/public_html/higgs/overview/electron/"); }

  TString histoPath =  hPath.Data();
  cout<<"histoPath:  "<<histoPath.Data()<<"  int Lumi: "<<intLumi<<endl;

  Bool_t doPhotons = 0, makeZjetsQt=0, doEBEE=0, doOverview=0;
  Bool_t doSB = 0, doTest=1;

  TFile* fda_2011A_DoubleMu_May10  = new TFile(Form("./%s/hhhh_DoubleMu_May10.root",histoPath.Data()));
  TFile* fda_2011A_DoubleEl_May10  = new TFile(Form("./%s/hhhh_DoubleMu_May10.root",histoPath.Data()));

  TFile* fda_2011A_DoubleMu_PromptV4  = new TFile(Form("./%s/hhhh_DoubleMu_PromptV4.root",histoPath.Data()));
  TFile* fda_2011A_DoubleEl_PromptV4  = new TFile(Form("./%s/hhhh_DoubleMu_PromptV4.root",histoPath.Data()));

  if (sel==1) TFile  *fData = (TFile*)fda_2011A_DoubleMu_PromptV4;
  if (sel==2) TFile  *fData = (TFile*)fda_2011A_DoubleEl_PromptV4;
  
  TFile* fda_Data  = new TFile(Form("./m_Data_%i.root", sel));  //Merged Data
  //TFile  *fData = (TFile*)fda_Data;
 
  // TFile* fmc_ZllG      = new TFile(Form("./%s/hhhh_WZ.root", histoPath.Data() ));
  // TFile* fmc_Wjets     = new TFile(Form("./%s/hhhh_Wjets.root",histoPath.Data() ));

  TFile* fmc_ZZ        = new TFile(Form("./%s/hhhh_ZZ.root",histoPath.Data() ));
  //TFile* fmc_tt        = new TFile(Form("./%s/hhhh_ttbar_1.root",histoPath.Data() ));
  TFile* fmc_WW        = new TFile(Form("./%s/hhhh_WW.root",histoPath.Data() ));
  TFile* fmc_WZ        = new TFile(Form("./%s/hhhh_WZ.root",histoPath.Data() ));
  //TFile* fmc_Top       = new TFile(Form("./%s/hhhh_tW.root",histoPath.Data() ));

  TFile* fmc_ggH200    = new TFile(Form("./%s/hhhh_ggHZZ200.root",histoPath.Data() ));
  //TFile* fmc_ggH250    = new TFile(Form("./%s/hhhh_ggHZZ250.root",histoPath.Data() ));
  TFile* fmc_ggH300    = new TFile(Form("./%s/hhhh_ggHZZ300.root",histoPath.Data() ));
  //TFile* fmc_ggH350    = new TFile(Form("./%s/hhhh_ggHZZ350.root",histoPath.Data() ));
  TFile* fmc_ggH400    = new TFile(Form("./%s/hhhh_ggHZZ400.root",histoPath.Data() ));



  //TFile* fmc_ggH600    = new TFile(Form("./%s/hhhh_SignalM600_HToZZ.root",histoPath.Data() ));
  //TFile* fmc_  = new TFile(Form("./%s/dir_%i__/hhhh.root",histoPath.Data() ));

  TFile* fmc_tt   = new TFile(Form("./m_ttbar_%i.root",sel ));
  TFile* fmc_Zjets  = new TFile(Form("./m_Zjets_%i.root",sel ));
  TFile* fmc_Top    = new TFile(Form("./m_Top_%i.root",sel ));
 
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
  list_overlay->Add(fmc_ggH200);
  list_overlay->Add(fmc_ggH300);
  //  list_overlay->Add(fmc_ggH400);

  list_signal = new TList();
  //list_signal->Add(fmc_ggH200);
  list_signal->Add(fmc_ggH300);
  list_signal->Add(fmc_ggH400);

  PrintYields(list_bg, list_signal, fData, sel, histoPath, "twiki");


  THStack *hs_met_et[nC], *hs_met2_et[nC], *hs_met_over_qt[nC], *hs_met2_over_qt[nC], *hs_di_qt[nC], *hs_met_et_ovQt[nC], *hs_mt[nC], *hs_mtZ[nC];
  THStack *hs_jet_N[nC], *hs_jet_dRlep1[nC], *hs_jet_dRlep2[nC], *hs_jet_b_N[nC], hs_jet_b_pt[nC];
  THStack *hs_di_mass[nC], *hs_di_mass_EB[nC], *hs_di_mass_EE[nC], *hs_di_mass_EX[nC];
  THStack *hs_met_dPhiLeadJet1[nC], *hs_met_dPhiLeadJet2[nC], *hs_met_dPhiClosJet1[nC], *hs_met_dPhiClosJet2[nC];
  THStack *hs_vtx_nPV_raw[nC], *hs_vtx_nPV_weight[nC];

  for(Int_t n = 0; n<nC; n++)
    {
      hs_met_et[n]       = makeStack(list_bg, Form("%s_et_%i", metType.Data(), n), intLumi);
      hs_mt[n]           = makeStack(list_bg, Form("%s_%i", mtType.Data(), n), intLumi);

      //hs_met2_et[n]     = makeStack(list_bg, Form("met2_et_%i", n), intLumi);
      hs_di_qt[n]      = makeStack(list_bg, Form("di_qt_%i",n), intLumi);
      hs_di_mass[n]    = makeStack(list_bg, Form("di_mass_%i",n), intLumi);
      hs_di_mass_EB[n] = makeStack(list_bg, Form("di_mass_EB_%i",n), intLumi);
      hs_di_mass_EE[n] = makeStack(list_bg, Form("di_mass_EE_%i",n), intLumi);
      hs_di_mass_EX[n] = makeStack(list_bg, Form("di_mass_EX_%i",n), intLumi);
      hs_jet_N[n]      = makeStack(list_bg, Form("jet_N_%i",n), intLumi);
      hs_jet_b_N[n]    = makeStack(list_bg, Form("jet_b_N_%i",n), intLumi); 

      hs_vtx_nPV_raw[n]    = makeStack(list_bg, Form("vtx_nPV_raw_%i",n), intLumi); 
      hs_vtx_nPV_weight[n] = makeStack(list_bg, Form("vtx_nPV_weight_%i",n), intLumi); 

      /*      hs_mtZ[n]          = makeStack(list_bg, Form("mtZ_%i", n), intLumi);

      hs_met_over_qt[n]  = makeStack(list_bg, Form("%s_over_qt_%i", metType.Data(), n), intLumi);
      hs_met2_over_qt[n] = makeStack(list_bg, Form("met2_over_qt_%i", n), intLumi);
      hs_met_et_ovQt[n]  = makeStack(list_bg, Form("%s_et_ovQt_%i", metType.Data(), n), intLumi);

      hs_jet_dRlep1[n]   = makeStack(list_bg, Form("jet_dRlep1_%i",n), intLumi);
      hs_jet_dRlep2[n]   = makeStack(list_bg, Form("jet_dRlep2_%i",n), intLumi);

      hs_met_dPhiLeadJet1[n] = makeStack(list_bg, Form("met2_dPhiLeadJet1_%i",n), intLumi);
      hs_met_dPhiLeadJet2[n] = makeStack(list_bg, Form("met2_dPhiLeadJet2_%i",n), intLumi);
      hs_met_dPhiClosJet1[n] = makeStack(list_bg, Form("met2_dPhiClosJet1_%i",n), intLumi);
      hs_met_dPhiClosJet2[n] = makeStack(list_bg, Form("met2_dPhiClosJet2_%i",n), intLumi);
      */
    }
  
  /*
  THStack *ph_met_et[nC], *ph_met2_et[nC], *ph_met_over_qt[nC], *ph_met2_over_qt[nC], *ph_di_qt[nC], *ph_met_et_ovQt[nC], *ph_mt[nC], *ph_mtZ[nC];
  THStack *ph_jet_N[nC], *ph_jet_dRlep1[nC], *ph_jet_dRlep2[nC], *ph_jet_b_N[nC], ph_jet_b_pt[nC];
  THStack *ph_di_mass[nC], *ph_di_mass_EB[nC], *ph_di_mass_EE[nC], *ph_di_mass_EX[nC];
  THStack *ph_met_dPhiLeadJet1[nC], *ph_met_dPhiLeadJet2[nC], *ph_met_dPhiClosJet1[nC], *ph_met_dPhiClosJet2[nC];
    
  for(Int_t n = 0; n<nC; n++)
    {
      ph_met_et[n]       = makeStack(list_bg2, Form("%s_et_%i", metType.Data(), n), intLumi);
      ph_met2_et[n]      = makeStack(list_bg2, Form("met2_et_%i", n), intLumi);
      ph_mt[n]           = makeStack(list_bg2, Form("%s_%i", mtType.Data(), n), intLumi);
      ph_mtZ[n]          = makeStack(list_bg2, Form("mtZ_%i", n), intLumi);

      ph_met_over_qt[n]  = makeStack(list_bg2, Form("%s_over_qt_%i", metType.Data(), n), intLumi);
      ph_met2_over_qt[n] = makeStack(list_bg2, Form("met2_over_qt_%i", n), intLumi);
      ph_met_et_ovQt[n]  = makeStack(list_bg2, Form("%s_et_ovQt_%i", metType.Data(), n), intLumi);

      ph_di_qt[n]      = makeStack(list_bg2, Form("di_qt_%i",n), intLumi);
      ph_di_mass[n]    = makeStack(list_bg2, Form("di_mass_%i",n), intLumi);
      ph_di_mass_EB[n] = makeStack(list_bg2, Form("di_mass_EB_%i",n), intLumi);
      ph_di_mass_EE[n] = makeStack(list_bg2, Form("di_mass_EE_%i",n), intLumi);
      ph_di_mass_EX[n] = makeStack(list_bg2, Form("di_mass_EX_%i",n), intLumi);
      ph_jet_N[n]      = makeStack(list_bg2, Form("jet_N_%i",n), intLumi);
      ph_jet_b_N[n]    = makeStack(list_bg2, Form("jet_b_N_%i",n), intLumi); 

      ph_jet_dRlep1[n]   = makeStack(list_bg2, Form("jet_dRlep1_%i",n), intLumi);
      ph_jet_dRlep2[n]   = makeStack(list_bg2, Form("jet_dRlep2_%i",n), intLumi);


      ph_met_dPhiLeadJet1[n] = makeStack(list_bg2, Form("met2_dPhiLeadJet1_%i",n), intLumi);
      ph_met_dPhiLeadJet2[n] = makeStack(list_bg2, Form("met2_dPhiLeadJet2_%i",n), intLumi);
      ph_met_dPhiClosJet1[n] = makeStack(list_bg2, Form("met2_dPhiClosJet1_%i",n), intLumi);
      ph_met_dPhiClosJet2[n] = makeStack(list_bg2, Form("met2_dPhiClosJet2_%i",n), intLumi);
    }
  */  
  

  TH1 * forLegend[21];
  list_legend = new TList();
  list_legend -> AddAll(list_bg);
  list_legend -> AddAll(list_overlay);
  //  TIter next(list_legend);

  TFile *ff[10]; 
  Int_t size = list_legend->GetSize();
  TH1 *hh[10]; 
  if(size>10) {cout<<"To many plots to overlay"<<endl; return 0;}
  for(Int_t n=0; n<size; n++)
    {
      ff[n] = (TFile*)list_legend->At(n);
      hh[n] = (TH1*) ff[n]->Get("Andrey/met2_et_6")->Clone();
      //cout<<n<<"   "<<ff[n] -> GetName()<<endl;
      forLegend[n] = (TH1F*)hh[n];
    }
 
  leg01 = new TLegend(0.53,0.7,0.95,0.95);
  leg01 -> SetNColumns(2);
  leg01 -> SetTextSize(0.04);
  
  leg01->AddEntry(forLegend[6], "Data","epl");
  leg01->AddEntry(forLegend[5],  "Z + jets","f");
  leg01->AddEntry(forLegend[4],  "ZZ","f");
  //leg01->AddEntry(forLegend[5],  "W + jets","f");
  leg01->AddEntry(forLegend[2],  "WW","f");
  leg01->AddEntry(forLegend[0],  "tW","f");
  leg01->AddEntry(forLegend[3],  "WZ","f");
  leg01->AddEntry(forLegend[1], "ttbar","f");  
  leg01->AddEntry(forLegend[7], "10xH200","f");
  leg01->AddEntry(forLegend[8], "10xH300","f");
  //  leg01->AddEntry(forLegend[9], "10xH400","f");
  
  leg01->SetFillColor(kWhite);
  
  hs_met_et[F0] -> Draw("hist");
  c1 -> SaveAs(imgpath+"ov01.png");
  TCanvas *c2 = new TCanvas("c2","example",600,700);
  
  

  if(doSB){
    cout<<"doing S/B"<<endl;
    optimalCuts result;

    /*
    result =  calculateOptimalCuts(1, hs_met_dPhiClosJet2[5]->Sum(), (TH1*)fmc_ggH400->Get("met2_dPhiClosJet2_5"));
    cout<<"from "<<(TMath::Pi()/40)*(result.bin1-1)<<" to "<< (TMath::Pi()/40)*result.bin2<<"\n S/B = "<<result.SB<<"  S/sqrt(B) = "<<result.SrootB<<"  S/sqrt(S+B) = "<<result.SrootSB<<endl;
    result =  calculateOptimalCuts(2, hs_met_dPhiClosJet2[5]->Sum(), (TH1*)fmc_ggH400->Get("met2_dPhiClosJet2_5"));
    cout<<"from "<<(TMath::Pi()/40)*(result.bin1-1)<<" to "<< (TMath::Pi()/40)*result.bin2<<"\n S/B = "<<result.SB<<"  S/sqrt(B) = "<<result.SrootB<<"  S/sqrt(S+B) = "<<result.SrootSB<<endl;
    result =  calculateOptimalCuts(3, hs_met_dPhiClosJet2[5]->Sum(), (TH1*)fmc_ggH400->Get("met2_dPhiClosJet2_5"));
    cout<<"from "<<(TMath::Pi()/40)*(result.bin1-1)<<" to "<< (TMath::Pi()/40)*result.bin2<<"\n S/B = "<<result.SB<<"  S/sqrt(B) = "<<result.SrootB<<"  S/sqrt(S+B) = "<<result.SrootSB<<endl;
    */

    for(Int_t ii=1; ii<=3; ii++){
      result =  calculateOptimalCuts(ii, hs_met_over_qt[6]->Sum(), (TH1*)fmc_ggH400->Get("met3_over_qt_6"));
      cout<<"from "<<0.1*(result.bin1-1)<<" to "<< 0.1*result.bin2<<"\n S/B = "<<result.SB<<"  S/sqrt(B) = "<<result.SrootB<<"  S/sqrt(S+B) = "<<result.SrootSB<<endl;
    }

    /*
    Float_t iSig = ((TH1*)fmc_ggH400->Get("met3_over_qt_6")) -> Integral(7,18); //from 0.6 to 1.8
    Float_t iBkg = (hs_met_over_qt[6]->Sum()) -> Integral(7,18);
    Float_t SB      = 0.1*iSig/iBkg;        //Higgs signal is 10*real in the histograms
    Float_t SrootB  = 0.1*iSig/sqrt(iBkg);
    Float_t SrootSB = 0.1*iSig/sqrt(iBkg + 0.1*iSig);
    cout<<"After cuts of projMET/qT from 0.6 to 1.8 "<<endl;
    cout<<iSig<<"  "<<iBkg<<" SB = "<<SB<<" SrootB = "<<SrootB<<" SrootSB = "<<SrootSB<<endl;
    cout<<endl;
        
    result =  calculateOptimalCuts(3, hs_mt[7]->Sum(), (TH1*)fmc_ggH400->Get("mt2_7"));
    cout<<"from "<<100+10*(result.bin1-1)<<" to "<< 100+10*result.bin2<<"\n S/B = "<<result.SB<<"  S/sqrt(B) = "<<result.SrootB<<"  S/sqrt(S+B) = "<<result.SrootSB<<endl;

    result =  calculateOptimalCuts(3, hs_mt[6]->Sum(), (TH1*)fmc_ggH400->Get("mt2_6"));
    cout<<"from "<<100+10*(result.bin1-1)<<" to "<< 100+10*result.bin2<<"\n S/B = "<<result.SB<<"  S/sqrt(B) = "<<result.SrootB<<"  S/sqrt(S+B) = "<<result.SrootSB<<endl;

    */

  }
  
  //-----------------------------------------------//
  // Zjets  backgroung from Photon samples (Anton) //
  //-----------------------------------------------//
  
  if(doTest){
    TString testpath("~/afs/public_html/test/");  
    
    drawMuliPlot("pfMET", 1, 0.001, 1000000, 0,5, hs_met_et[6], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p01.png");
    drawMuliPlot("MT", 1, 0.001, 1000000, 0,5, hs_mt[6], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p02.png");

    cout<<"dbg"<<endl;

    drawMuliPlot("MT 250", 0, 0.001, 10, 0,5, hs_mt[7], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p03.png");
    drawMuliPlot("MT 300", 0, 0.001, 10, 0,5, hs_mt[8], c2, leg01, list_overlay, intLumi); 
    c2 -> SaveAs(testpath+"p04.png");
       
    drawMuliPlot("N b-jets", 0, 0.0, 30, 0.,5, hs_jet_b_N[9], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(testpath+"p05.png");
    drawMuliPlot("N jets", 0, 0.0, 20, 0.,5, hs_jet_N[9], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(testpath+"p06.png");


    //PrintYields(ph_mt[nn], (TH1*)fmc_ggH400 -> Get(Form("mt2_%i",nn))->Clone(), (TH1*)fData ->Get(Form("mt2_%i",nn))->Clone(), sel, nn, hPath);

    // TH1* sig  =  fmc_ggH400 -> Get("mt2_7") -> Clone();
    //TH1* data =  fData ->Get("mt2_7") -> Clone();
  
    // for(Int_t nn=3; nn<=10; nn++)
    //  PrintYields(hs_mt[nn], (TH1*)fmc_ggH400 -> Get(Form("mt2_%i",nn))->Clone(), (TH1*)fData ->Get(Form("mt2_%i",nn))->Clone(), sel, nn, hPath);

    //for(Int_t nn=4; nn<8; nn++)
    
    //drawMuliPlot("projMET/q_{T}", 1, 0.001, 1000000, 0,3, hs_met_over_qt[5], c2, leg01, list_overlay, intLumi);
    //c2 -> SaveAs(testpath+"p02.png");

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
  

  if(makeZjetsQt){
    
    TFile *fA_Zj = new TFile(Form("./forAnton_Zjets_%i.root",sel), "RECREATE");
    TFile *fA_Da = new TFile(Form("./forAnton_Data_%i.root",sel), "RECREATE");
    TFile *fA_Ds = new TFile(Form("./forAnton_Data_sbtr_%i.root",sel), "RECREATE");
    
    fA_Zj -> cd();
    TH1F *Zjets_qt   = fmc_Zjets->Get("di_qt_4")   -> Clone();
    TH1F *Zjets_mass = fmc_Zjets->Get("di_mass_4") -> Clone();
    TH1F *Zjets_met  = fmc_Zjets->Get("met3_et_5") -> Clone();
    TH1F *Zjets_met2 = fmc_Zjets->Get("met2_et_5") -> Clone();
    TH1F *Zjets_mt5  = fmc_Zjets->Get("mt2_5")   -> Clone();
    TH1F *Zjets_mt6  = fmc_Zjets->Get("mt2_6")   -> Clone();
    TH1F *Zjets_nj5  = fmc_Zjets->Get("jet_N_5") -> Clone();
    TH1F *Zjets_nj6  = fmc_Zjets->Get("jet_N_6") -> Clone();
    TH1F *Zjets_mOq5 = fmc_Zjets->Get("met2_over_qt_5") -> Clone();
    TH1F *Zjets_mOq6 = fmc_Zjets->Get("met2_over_qt_6") -> Clone();
    Zjets_qt   -> Write();
    Zjets_mass -> Write();
    
    fA_Da -> cd();
    TH1F *da_qt   = fData->Get("di_qt_4")   -> Clone();
    TH1F *da_mass = fData->Get("di_mass_4") -> Clone();
    TH1F *da_met  = fData->Get("met3_et_5") -> Clone();
    TH1F *da_met2 = fData->Get("met2_et_5") -> Clone();
    TH1F *da_mt5  = fData->Get("mt2_5")     -> Clone();
    TH1F *da_mt6  = fData->Get("mt2_6")     -> Clone();
    TH1F *da_nj5  = fData->Get("jet_N_5")   -> Clone();
    TH1F *da_nj6  = fData->Get("jet_N_6")   -> Clone();
    TH1F *da_mOq5 = fData->Get("met2_over_qt_5") -> Clone();
    TH1F *da_mOq6 = fData->Get("met2_over_qt_6") -> Clone();
    da_qt   -> Write();
    da_mass -> Write();
    
    
    fA_Ds -> cd();
    TH1F *ds_qt   = (TH1F*)ph_di_qt[4]->Sum()->Clone();
    TH1F *ds_mass = (TH1F*)ph_di_mass[4]->Sum()->Clone();
    da_qt   -> Add(ds_qt,-1);
    da_mass -> Add(ds_mass,-1);

    da_qt   -> Write();
    da_mass -> Write();

    fA_Ds -> Close();
    fA_Da -> Close();
    fA_Zj -> Close();
  }


  // Mll in electons, EE vs EB
  if(doEBEE){
    drawMuliPlot("Leptons in Barrel, |#eta|<1.444,  M(ll)", 1, 0.001, 1000000, 0,4, hs_di_mass_EB[F0], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(testpath+"p01.png");
    drawMuliPlot("Leptons in Endcap, |#eta|>1.566, M(ll)", 1, 0.001, 1000000, 0,4, hs_di_mass_EE[F0], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(testpath+"p02.png");
    drawMuliPlot("Leptons in EB/EE, mixed, M(ll)", 1, 0.001, 1000000, 0,4, hs_di_mass_EX[F0], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(testpath+"p03.png");
  }

  if(doOverview){
    drawMuliPlot("projMET", 1, 0.001, 1000000, 0,5, hs_met_et[F0], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(imgpath+"ov01.png");
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
    
    drawMuliPlot("M(ll)", 1, 0.001, 1000000, 0,2, hs_di_mass[F0], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(imgpath+"ov07.png");
    drawMuliPlot("MT", 1, 0.001, 1000000, 0,2, hs_mt[F0], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(imgpath+"ov08.png");
    
    drawMuliPlot("M(ll)", 1, 0.001, 1000000, 0,5, hs_di_mass[F1], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(imgpath+"ov09.png");
    drawMuliPlot("MT", 1, 0.001, 1000000, 0,5, hs_mt[F1], c2, leg01, list_overlay, intLumi);
    c2 -> SaveAs(imgpath+"ov10.png");
    
  /*    
      drawMuliPlot("M(ll)", 1, 0.001, 10000000, 0,2, hs_di_mass[F2], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov11.png");
      drawMuliPlot("MT", 1, 0.001, 1000000, 0,2, hs_mt[F2], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov12.png");
      
      drawMuliPlot("N jets", 1, 0.001, 1000000, 0,5, hs_jet_N[F0], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov13.png");
      drawMuliPlot("N b-jets", 1, 0.001, 1000000, 0,5, hs_jet_b_N[F0], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov14.png");
      
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
      drawMuliPlot("#Delta#phi(MET, closest jet), p_{T}>20, |#eta|<2.4", 1, 0.001, 10000000, 0,2, hs_met_dPhiClosJet1[F0], c2, leg01, list_overlay, intLumi);
      c2 -> SaveAs(imgpath+"ov22.png");
      
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


  /*

  c1 -> SetLogy(0);
  fmc_ggH200 -> cd();
  met3_et_ovQt_9 -> Draw("colz");
  met3_et_ovQt_9 -> SetTitle(";projMET; projMET/q_{T}");
  c1 -> SaveAs(imgpath+"ov05.png");

  hs_met_et_ovQt[F1] -> Draw("colz");
  hs_met_et_ovQt[F1] -> SetTitle(";projMET; projMET/q_{T}");
    c1 -> SaveAs(imgpath+"ov06.png");
  */

  cout<<"\n\n   -- end of job  --"<<endl;
  
}

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




//void *THStack::PrintYields(Int_t sel, Int_t ver, string option="tex")
void PrintYields(THStack *stack[], TH1 *signal, TH1 *data, Int_t sel, Int_t num, TString ver, string option="tex")
{
  Int_t bins = 0;
  ofstream oo;
  if (sel==1)  oo.open(Form("yields_%s_muons.txt",ver.Data()), ofstream::out);
  if (sel==2)  oo.open(Form("yields_%s_electrons.txt",ver.Data()), ofstream::out);
  oo.precision(2); oo.setf(ios::fixed, ios::floatfield);

  TList * mylist = (TList*)stack[0]->GetHists();
  TIter next(mylist);
  TH1 *hh  = (TH1*) mylist -> First() ->Clone();
  bins = hh -> GetNbinsX();
  if(option=="tex")    oo<<hh->Integral(0,bins+1)<<"\t& ";
  if(option=="twiki")  oo<<hh->Integral(0,bins+1)<<"\t& ";

  TObject *obj; 
  while ((obj = next()))
    {
      // cout<<obj->GetName()<<endl;
      if(obj == mylist->First()) continue;
      hh = (TH1*)obj;     
      if(option=="tex")    oo<<hh->Integral(0,bins+1)<<"\t& ";
      if(option=="twiki")  oo<<hh->Integral(0,bins+1)<<"\t& ";
    }

  TH1* sig = (TH1*)signal->Clone();
  TH1* bkg = (TH1*)stack->Sum();

  Int_t sig_nBins = sig->GetNbinsX();
  Int_t bkg_nBins = bkg->GetNbinsX();
  Int_t nBins=0;
  if (sig_nBins!=bkg_nBins) cout<<" IN calculateOptimalCuts \n WARNING:  different number of bis"<<endl;
  else nBins = sig_nBins;

  Float_t iSig = signal -> Integral(0,nBins+1);  //Count overflows
  Float_t iBkg = bkg -> Integral(0,nBins+1);
  Float_t SB      = 0.1*iSig/iBkg;        //Higgs signal is 10*real in the histograms
  Float_t SrootB  = 0.1*iSig/sqrt(iBkg);
  Float_t SrootSB = 0.1*iSig/sqrt(iBkg + 0.1*iSig);

  if(option=="tex") {
    oo<<"\t"<<data->Integral(0,nBins+1)<<"\t&";

    oo<<"\t"<<iBkg<<"\t& ";
    
    oo<<"\t"<<0.1*iSig<<"\t&";

    oo.precision(3);
    oo<<"\t"<<SB<<"\t&"<<SrootB<<"\t \\\\ \n \\hline"<<endl;
  }

  oo<<"end of printing yields"<<endl;
  oo.close();
}



void PrintYields(TList *bgList, TList *sigList, TFile *dataFile, Int_t sel, TString path, string option="tex")
{
  Int_t bins = 0;
  ofstream oo;
  oo.open(Form("%s/yields_%s.txt",path.Data(), option.c_str()), ofstream::out);
  oo.precision(2); oo.setf(ios::fixed, ios::floatfield);

  TString beginLine("");
  TString endLine("");
  TString separator("");
  TString title("");
  if(option=="tex"){
    title = " sel & WW & WZ   \\\\ \\hline";
    beginLine = "   ";
    endLine   = "\\\\ \\hline";
    separator = "\t &";
  }
  if(option=="twiki"){
    title = "| sel | *top*  | *ttbar*  | *WZ*  | *WW* | *ZZ*  | *Zjets*  | *Data*  | *Total bg*  | *higgs* | *S/B* |";
    beginLine = "|   ";
    endLine   = "\t |";
    separator = "\t |";
  }

  oo<<title<<endl;

  for(Int_t j = 6; j<12; j++)
    {
      Int_t size = bgList->GetSize();
      TFile *ff[10];    TH1 *hh[10];
      if(size>10) {cout<<"To many plots to overlay"<<endl; return 0;}

      oo<<beginLine;
      oo<<j<<separator;

      Float_t total_bg = 0, total_bgError;

      for(Int_t n=0; n<size; n++)
	{
	  ff[n] = (TFile*)bgList->At(n);
	  hh[n] = (TH1*)ff[n]->Get( Form("Andrey/mt2_%i",j) )->Clone();
	  cout<<n<<" file  "<<ff[n] -> GetName()<<"   histoname: "<<hh[n]->GetName()<<endl;
	
	  bins = hh[n] -> GetNbinsX();
	  oo<<hh[n]->Integral(0,bins+1)<<separator;
	  total_bg += hh[n]->Integral(0,bins+1);
	  total_bgError += 0;
	}      


      TH1* data = (TH1*)dataFile->Get(  Form("Andrey/mt2_%i",j) )->Clone();
      TFile* sigFile  = (TFile*)sigList->First();//->Get(  Form("Andrey/mt2_%i",j) )->Clone();
      TH1* sig  = (TH1*)sigFile->Get(Form("Andrey/mt2_%i",j))->Clone();
      
      Int_t sig_nBins = sig->GetNbinsX();
      Int_t bkg_nBins = hh[0]->GetNbinsX();
      Int_t nBins=0;
      if (sig_nBins!=bkg_nBins) cout<<" IN yields \n WARNING:  different number of bis"<<endl;
      else nBins = sig_nBins;
      
      Float_t iSig = sig -> Integral(0,nBins+1);  //Count overflows
      Float_t iBkg = total_bg;
      Float_t SB   = 0.1*iSig/iBkg;        //Higgs signal is 10*real in the histograms
      //Float_t SrootB  = 0.1*iSig/sqrt(iBkg);
      //Float_t SrootSB = 0.1*iSig/sqrt(iBkg + 0.1*iSig);
      
      oo.precision(0);
      oo<<data->Integral(0,nBins+1);
      oo.precision(2);
      oo<<separator<<iBkg<<separator<<0.1*iSig<<separator;
      oo.precision(3);
      oo<<SB<<endLine<<endl;  
      
  }



  oo<<"end of printing yields"<<endl;
  oo.close();
  
}


 
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
