void plot(){
  gROOT->ProcessLine(".L ../data/tdrstyle.C");
  setTDRStyle();
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLegendBorderSize(0);
  gROOT->ForceStyle();
  cout.precision(3); cout.setf(ios::fixed, ios::floatfield);
  TH1::SetDefaultSumw2(kTRUE);

#define nCuts 10
#define F0 6
#define F1 5
#define F2 2

  TString histoPath("v03");
  //TString sample("none");
  Int_t sel = 1; //1-muon, 2 -electron

  //Types of met: met - pfMet, met1 - type1 corrected, met2 - pfMet passed Noise filters, 
  //met3 - projMet, met4 - puProj corrected met (those two are passed Noise filters) 
  TString metType("met4"); 


  if (sel==1)  TString imgpath("~/afs/public_html/higgs/overview/muon/");  
  if (sel==2)  TString imgpath("~/afs/public_html/higgs/overview/electron/");  
  //TString imgpath(Form("~/afs/public_html/higgs/%s/",sample.Data()));  

  TFile* fmc_ZllG      = new TFile(Form("./%s/dir_%i_MC_ZllG_/hhhh.root", histoPath.Data(), sel));
  TFile* fmc_ZeeGamma  = new TFile(Form("./%s/dir_%i_MC_ZeeGamma_/hhhh.root", histoPath.Data(), sel));
  TFile* fmc_Wjets     = new TFile(Form("./%s/dir_%i_MC_Wjets_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_ZZtoAny   = new TFile(Form("./%s/dir_%i_MC_ZZtoAny_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_ggH400    = new TFile(Form("./%s/dir_%i_MC_ggH400_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_tS        = new TFile(Form("./%s/dir_%i_MC_tSchannel_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_tT        = new TFile(Form("./%s/dir_%i_MC_tTchannel_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_tW        = new TFile(Form("./%s/dir_%i_MC_tWchannel_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_tt        = new TFile(Form("./%s/dir_%i_MC_tt2l2nu2b_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_DYmumu    = new TFile(Form("./%s/dir_%i_MC_DYmumu_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_DYee      = new TFile(Form("./%s/dir_%i_MC_DYee_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_WW        = new TFile(Form("./%s/dir_%i_MC_WW_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_WZ        = new TFile(Form("./%s/dir_%i_MC_WZ_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_Zbb0  = new TFile(Form("./%s/dir_%i_MC_Zbb0_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_Zbb1  = new TFile(Form("./%s/dir_%i_MC_Zbb1_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_Zbb2  = new TFile(Form("./%s/dir_%i_MC_Zbb2_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_Zbb3  = new TFile(Form("./%s/dir_%i_MC_Zbb3_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_Zcc0  = new TFile(Form("./%s/dir_%i_MC_Zcc0_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_Zcc1  = new TFile(Form("./%s/dir_%i_MC_Zcc1_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_Zcc2  = new TFile(Form("./%s/dir_%i_MC_Zcc2_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_Zcc3  = new TFile(Form("./%s/dir_%i_MC_Zcc3_/hhhh.root",histoPath.Data(), sel));
  //TFile* fmc_  = new TFile(Form("./%s/dir_%i__/hhhh.root",histoPath.Data(), sel));

  TFile* fda_2011A_DoubleMu  = new TFile(Form("./%s/dir_1_2011A_May10_DoubleMu_/hhhh.root",histoPath.Data()));
  TFile* fda_2011A_DoubleEl  = new TFile(Form("./%s/dir_2_2011A_May10_DoubleElectron_/hhhh.root",histoPath.Data()));
  


  // TFile* fmc_ZeeGamma = new TFile(Form("./%s/histogramsMC_%s.root", histoPath.Data(), tune.Data()));

  TH1F *ZeeGamma_met_phi[nCuts], *ZeeGamma_met_et[nCuts], *ZeeGamma_met_over_qt[nCuts];
  TH1F *ZeeGamma_met1_phi[nCuts], *ZeeGamma_met1_et[nCuts], *ZeeGamma_met1_over_qt[nCuts];

  TH1F *ZllG_met_phi[nCuts], *ZllG_met_et[nCuts], *ZllG_met_over_qt[nCuts];
  TH1F *Wjets_met_phi[nCuts], *Wjets_met_et[nCuts], *Wjets_met_over_qt[nCuts];
  TH1F *ZZtoAny_met_phi[nCuts], *ZZtoAny_met_et[nCuts], *ZZtoAny_met_over_qt[nCuts];
  TH1F *ggH400_met_phi[nCuts], *ggH400_met_et[nCuts], *ggH400_met_over_qt[nCuts];
  TH1F *t_all_met_phi[nCuts], *t_all_met_et[nCuts], *t_all_met_over_qt[nCuts];
  TH1F *tS_met_phi[nCuts], *tS_met_et[nCuts], *tS_met_over_qt[nCuts];
  TH1F *tT_met_phi[nCuts], *tT_met_et[nCuts], *tT_met_over_qt[nCuts];
  TH1F *tW_met_phi[nCuts], *tW_met_et[nCuts], *tW_met_over_qt[nCuts];
  TH1F *tt_met_phi[nCuts], *tt_met_et[nCuts], *tt_met_over_qt[nCuts];
  TH1F *data_met_phi[nCuts], *data_met_et[nCuts], *data_met_over_qt[nCuts];
  TH1F *Zjets_met_phi[nCuts], *Zjets_met_et[nCuts], *Zjets_met_over_qt[nCuts];
  TH1F *DYmumu_met_phi[nCuts], *DYmumu_met_et[nCuts], *DYmumu_met_over_qt[nCuts];
  TH1F *DYee_met_phi[nCuts], *DYee_met_et[nCuts], *DYee_met_over_qt[nCuts];
  TH1F *WW_met_phi[nCuts], *WW_met_et[nCuts], *WW_met_over_qt[nCuts];
  TH1F *WZ_met_phi[nCuts], *WZ_met_et[nCuts], *WZ_met_over_qt[nCuts];

  TH1F *Zjets_met_phi[nCuts], *Zjets_met_et[nCuts], *Zjets_met_over_qt[nCuts];
  TH1F *Zbb0_met_phi[nCuts], *Zbb0_met_et[nCuts], *Zbb0_met_over_qt[nCuts];
  TH1F *Zbb1_met_phi[nCuts], *Zbb1_met_et[nCuts], *Zbb1_met_over_qt[nCuts];
  TH1F *Zbb2_met_phi[nCuts], *Zbb2_met_et[nCuts], *Zbb2_met_over_qt[nCuts];
  TH1F *Zbb3_met_phi[nCuts], *Zbb3_met_et[nCuts], *Zbb3_met_over_qt[nCuts];
  TH1F *Zcc0_met_phi[nCuts], *Zcc0_met_et[nCuts], *Zcc0_met_over_qt[nCuts];
  TH1F *Zcc1_met_phi[nCuts], *Zcc1_met_et[nCuts], *Zcc1_met_over_qt[nCuts];
  TH1F *Zcc2_met_phi[nCuts], *Zcc2_met_et[nCuts], *Zcc2_met_over_qt[nCuts];
  TH1F *Zcc3_met_phi[nCuts], *Zcc3_met_et[nCuts], *Zcc3_met_over_qt[nCuts];
  //TH1F *_met_phi[nCuts], *_met_et[nCuts], *_met_over_qt[nCuts];
  //TH1F *_met_phi[nCuts], *_met_et[nCuts], *_met_over_qt[nCuts];
  //TH1F *_met_phi[nCuts], *_met_et[nCuts], *_met_over_qt[nCuts];

  for(Int_t n=0; n<nCuts; n++)
    {
      ZeeGamma_met_phi[n]     = (TH1F*)fmc_ZeeGamma->Get(Form("%s_phi_%i",metType.Data(),n));
      ZeeGamma_met_et[n]      = (TH1F*)fmc_ZeeGamma->Get(Form("%s_et_%i",metType.Data(),n));
      ZeeGamma_met_over_qt[n] = (TH1F*)fmc_ZeeGamma->Get(Form("%s_over_qt_%i",metType.Data(),n));

      ZllG_met_phi[n]     = (TH1F*)fmc_ZllG->Get(Form("%s_phi_%i",metType.Data(),n));
      ZllG_met_et[n]      = (TH1F*)fmc_ZllG->Get(Form("%s_et_%i",metType.Data(),n));
      ZllG_met_over_qt[n] = (TH1F*)fmc_ZllG->Get(Form("%s_over_qt_%i",metType.Data(),n));

      Wjets_met_phi[n]     = (TH1F*)fmc_Wjets->Get(Form("%s_phi_%i",metType.Data(),n));
      Wjets_met_et[n]      = (TH1F*)fmc_Wjets->Get(Form("%s_et_%i",metType.Data(),n));
      Wjets_met_over_qt[n] = (TH1F*)fmc_Wjets->Get(Form("%s_over_qt_%i",metType.Data(),n));

      ZZtoAny_met_phi[n]     = (TH1F*)fmc_ZZtoAny->Get(Form("%s_phi_%i",metType.Data(),n));
      ZZtoAny_met_et[n]      = (TH1F*)fmc_ZZtoAny->Get(Form("%s_et_%i",metType.Data(),n));
      ZZtoAny_met_over_qt[n] = (TH1F*)fmc_ZZtoAny->Get(Form("%s_over_qt_%i",metType.Data(),n));

      tS_met_phi[n]     = (TH1F*)fmc_tS->Get(Form("%s_phi_%i",metType.Data(),n));
      tS_met_et[n]      = (TH1F*)fmc_tS->Get(Form("%s_et_%i",metType.Data(),n));
      tS_met_over_qt[n] = (TH1F*)fmc_tS->Get(Form("%s_over_qt_%i",metType.Data(),n));

      tT_met_phi[n]     = (TH1F*)fmc_tT->Get(Form("%s_phi_%i",metType.Data(),n));
      tT_met_et[n]      = (TH1F*)fmc_tT->Get(Form("%s_et_%i",metType.Data(),n));
      tT_met_over_qt[n] = (TH1F*)fmc_tT->Get(Form("%s_over_qt_%i",metType.Data(),n));

      tW_met_phi[n]     = (TH1F*)fmc_tW->Get(Form("%s_phi_%i",metType.Data(),n));
      tW_met_et[n]      = (TH1F*)fmc_tW->Get(Form("%s_et_%i",metType.Data(),n));
      tW_met_over_qt[n] = (TH1F*)fmc_tW->Get(Form("%s_over_qt_%i",metType.Data(),n));

      tt_met_phi[n]     = (TH1F*)fmc_tt->Get(Form("%s_phi_%i",metType.Data(),n));
      tt_met_et[n]      = (TH1F*)fmc_tt->Get(Form("%s_et_%i",metType.Data(),n));
      tt_met_over_qt[n] = (TH1F*)fmc_tt->Get(Form("%s_over_qt_%i",metType.Data(),n));

      ggH400_met_phi[n]     = (TH1F*)fmc_ggH400->Get(Form("%s_phi_%i",metType.Data(),n));
      ggH400_met_et[n]      = (TH1F*)fmc_ggH400->Get(Form("%s_et_%i",metType.Data(),n));
      ggH400_met_over_qt[n] = (TH1F*)fmc_ggH400->Get(Form("%s_over_qt_%i",metType.Data(),n));

      if (sel==1)
	{
	  data_met_phi[n]     = (TH1F*)fda_2011A_DoubleMu->Get(Form("%s_phi_%i",metType.Data(),n));
	  data_met_et[n]      = (TH1F*)fda_2011A_DoubleMu->Get(Form("%s_et_%i",metType.Data(),n));
	  data_met_over_qt[n] = (TH1F*)fda_2011A_DoubleMu->Get(Form("%s_over_qt_%i",metType.Data(),n));
	}
      if(sel==2)
	{
	  data_met_phi[n]     = (TH1F*)fda_2011A_DoubleEl->Get(Form("%s_phi_%i",metType.Data(),n));
	  data_met_et[n]      = (TH1F*)fda_2011A_DoubleEl->Get(Form("%s_et_%i",metType.Data(),n));
	  data_met_over_qt[n] = (TH1F*)fda_2011A_DoubleEl->Get(Form("%s_over_qt_%i",metType.Data(),n));
	}

      DYmumu_met_phi[n]     = (TH1F*)fmc_DYmumu->Get(Form("%s_phi_%i",metType.Data(),n));
      DYmumu_met_et[n]      = (TH1F*)fmc_DYmumu->Get(Form("%s_et_%i",metType.Data(),n));
      DYmumu_met_over_qt[n] = (TH1F*)fmc_DYmumu->Get(Form("%s_over_qt_%i",metType.Data(),n));

      DYee_met_phi[n]     = (TH1F*)fmc_DYee->Get(Form("%s_phi_%i",metType.Data(),n));
      DYee_met_et[n]      = (TH1F*)fmc_DYee->Get(Form("%s_et_%i",metType.Data(),n));
      DYee_met_over_qt[n] = (TH1F*)fmc_DYee->Get(Form("%s_over_qt_%i",metType.Data(),n));

      WW_met_phi[n]       = (TH1F*)fmc_WW->Get(Form("%s_phi_%i",metType.Data(),n));
      WW_met_et[n]        = (TH1F*)fmc_WW->Get(Form("%s_et_%i",metType.Data(),n));
      WW_met_over_qt[n]   = (TH1F*)fmc_WW->Get(Form("%s_over_qt_%i",metType.Data(),n));

      WZ_met_phi[n]       = (TH1F*)fmc_WZ->Get(Form("%s_phi_%i",metType.Data(),n));
      WZ_met_et[n]        = (TH1F*)fmc_WZ->Get(Form("%s_et_%i",metType.Data(),n));
      WZ_met_over_qt[n]   = (TH1F*)fmc_WZ->Get(Form("%s_over_qt_%i",metType.Data(),n));

      Zbb0_met_phi[n]     = (TH1F*)fmc_Zbb0->Get(Form("%s_phi_%i",metType.Data(),n));
      Zbb0_met_et[n]      = (TH1F*)fmc_Zbb0->Get(Form("%s_et_%i",metType.Data(),n));
      Zbb0_met_over_qt[n] = (TH1F*)fmc_Zbb0->Get(Form("%s_over_qt_%i",metType.Data(),n));

      Zbb1_met_phi[n]     = (TH1F*)fmc_Zbb1->Get(Form("%s_phi_%i",metType.Data(),n));
      Zbb1_met_et[n]      = (TH1F*)fmc_Zbb1->Get(Form("%s_et_%i",metType.Data(),n));
      Zbb1_met_over_qt[n] = (TH1F*)fmc_Zbb1->Get(Form("%s_over_qt_%i",metType.Data(),n));

      Zbb2_met_phi[n]     = (TH1F*)fmc_Zbb2->Get(Form("%s_phi_%i",metType.Data(),n));
      Zbb2_met_et[n]      = (TH1F*)fmc_Zbb2->Get(Form("%s_et_%i",metType.Data(),n));
      Zbb2_met_over_qt[n] = (TH1F*)fmc_Zbb2->Get(Form("%s_over_qt_%i",metType.Data(),n));

      Zbb3_met_phi[n]     = (TH1F*)fmc_Zbb3->Get(Form("%s_phi_%i",metType.Data(),n));
      Zbb3_met_et[n]      = (TH1F*)fmc_Zbb3->Get(Form("%s_et_%i",metType.Data(),n));
      Zbb3_met_over_qt[n] = (TH1F*)fmc_Zbb3->Get(Form("%s_over_qt_%i",metType.Data(),n));

      Zcc0_met_phi[n]     = (TH1F*)fmc_Zcc0->Get(Form("%s_phi_%i",metType.Data(),n));
      Zcc0_met_et[n]      = (TH1F*)fmc_Zcc0->Get(Form("%s_et_%i",metType.Data(),n));
      Zcc0_met_over_qt[n] = (TH1F*)fmc_Zcc0->Get(Form("%s_over_qt_%i",metType.Data(),n));

      Zcc1_met_phi[n]     = (TH1F*)fmc_Zcc1->Get(Form("%s_phi_%i",metType.Data(),n));
      Zcc1_met_et[n]      = (TH1F*)fmc_Zcc1->Get(Form("%s_et_%i",metType.Data(),n));
      Zcc1_met_over_qt[n] = (TH1F*)fmc_Zcc1->Get(Form("%s_over_qt_%i",metType.Data(),n));

      Zcc2_met_phi[n]     = (TH1F*)fmc_Zcc2->Get(Form("%s_phi_%i",metType.Data(),n));
      Zcc2_met_et[n]      = (TH1F*)fmc_Zcc2->Get(Form("%s_et_%i",metType.Data(),n));
      Zcc2_met_over_qt[n] = (TH1F*)fmc_Zcc2->Get(Form("%s_over_qt_%i",metType.Data(),n));

      Zcc3_met_phi[n]     = (TH1F*)fmc_Zcc3->Get(Form("%s_phi_%i",metType.Data(),n));
      Zcc3_met_et[n]      = (TH1F*)fmc_Zcc3->Get(Form("%s_et_%i",metType.Data(),n));
      Zcc3_met_over_qt[n] = (TH1F*)fmc_Zcc3->Get(Form("%s_over_qt_%i",metType.Data(),n));

      /*

      _met_phi[n]     = (TH1F*)fmc_->Get(Form("%s_phi_%i",metType.Data(),n));
      _met_et[n]      = (TH1F*)fmc_->Get(Form("%s_et_%i",metType.Data(),n));
      _met_over_qt[n] = (TH1F*)fmc_->Get(Form("%s_over_qt_%i",metType.Data(),n));
      */

    }

  Float_t intLumi = 191.0;
  Float_t ev_ZeeGamma  = ZeeGamma_met_et[0] -> GetEntries();
  Float_t ev_ZllG      = 165000.;
  Float_t ev_Wjets     = Wjets_met_et[0]    -> GetEntries();
  Float_t ev_tt        = tt_met_et[0]       -> GetEntries();
  Float_t ev_tS        = tS_met_et[0]       -> GetEntries();
  Float_t ev_tT        = tT_met_et[0]       -> GetEntries();
  Float_t ev_tW        = tW_met_et[0]       -> GetEntries();
  Float_t ev_ZZtoAny   = ZZtoAny_met_et[0]  -> GetEntries();
  Float_t ev_ggH400    = ggH400_met_et[0]   -> GetEntries();
  Float_t ev_DYmumu = 1984154.;
  Float_t ev_DYee   = 1992276.;
  Float_t ev_WW     = 110000.;
  Float_t ev_WZ     = 110000.;
  Float_t ev_Zbb0   = 347393.;
  Float_t ev_Zbb1   = 2025811.;
  Float_t ev_Zbb2   = 10842.;
  Float_t ev_Zbb3   = 10898.;
  Float_t ev_Zcc0   = 438067.;
  Float_t ev_Zcc1   = 184365.;
  Float_t ev_Zcc2   = 10735.;
  Float_t ev_Zcc3   = 10177.;


  cout<<"ZeeG: "<<ev_ZeeGamma<<"  Wjets: "<<ev_Wjets<<"  ZZtoAny: "<<ev_ZZtoAny<<" H400: "<<ev_ggH400<<endl;
  cout<<"tt: "<<ev_tt<<"    tT: "<<ev_tT<<"    tS: "<<ev_tS<<"    tW: "<<ev_tW<<endl;
 
  Float_t sigma_ZeeGamma  = 33.73; 
  Float_t sigma_ZllG      = 14.3;
  Float_t sigma_Wjets     = 31314.0;
  Float_t sigma_tt        = 16.5;
  Float_t sigma_tS        = 1.36;
  Float_t sigma_tT        = 20.9;
  Float_t sigma_tW        = 10.6;
  Float_t sigma_ZZtoAny   = 5.9;
  Float_t sigma_ggH400    = 0.02229;
  Float_t sigma_DYmumu    = 1666.;
  Float_t sigma_DYee      = 1666.;
  Float_t sigma_WW        = 4.51;
  Float_t sigma_WZ        = 0.596;

  Float_t sigma_Zbb0      = 2.92;
  Float_t sigma_Zbb1      = 1.65;
  Float_t sigma_Zbb2      = 0.623;
  Float_t sigma_Zbb3      = 0.274;
  Float_t sigma_Zcc0      = 2.92;
  Float_t sigma_Zcc1      = 1.63;
  Float_t sigma_Zcc2      = 0.627;
  Float_t sigma_Zcc3      = 0.281;
  // Float_t sigma_ = ;

  Float_t sc_ZllG      = intLumi*sigma_ZeeGamma/ev_ZeeGamma;
  Float_t sc_ZeeGamma  = intLumi*sigma_ZeeGamma/ev_ZeeGamma;
  Float_t sc_Wjets     = intLumi*sigma_Wjets/ev_Wjets;
  Float_t sc_tt        = intLumi*sigma_tt/ev_tt;
  Float_t sc_tS        = intLumi*sigma_tS/ev_tS;
  Float_t sc_tT        = intLumi*sigma_tT/ev_tT;
  Float_t sc_tW        = intLumi*sigma_tW/ev_tW;
  Float_t sc_ZZtoAny   = intLumi*sigma_ZZtoAny/ev_ZZtoAny;
  Float_t sc_ggH400    = 10*intLumi*sigma_ggH400/ev_ggH400;
  Float_t sc_DYmumu    = intLumi*sigma_DYmumu/ev_DYmumu;  
  Float_t sc_DYee      = intLumi*sigma_DYee/ev_DYee;  
  Float_t sc_WW        = intLumi*sigma_WW/ev_WW;  
  Float_t sc_WZ        = intLumi*sigma_WZ/ev_WZ;  
  
  Float_t sc_Zbb0  = intLumi*sigma_Zbb0/ev_Zbb0;
  Float_t sc_Zbb1  = intLumi*sigma_Zbb1/ev_Zbb1;
  Float_t sc_Zbb2  = intLumi*sigma_Zbb2/ev_Zbb2;
  Float_t sc_Zbb3  = intLumi*sigma_Zbb3/ev_Zbb3;
  Float_t sc_Zcc0  = intLumi*sigma_Zcc0/ev_Zcc0;
  Float_t sc_Zcc1  = intLumi*sigma_Zcc1/ev_Zcc1;
  Float_t sc_Zcc2  = intLumi*sigma_Zcc2/ev_Zcc2;
  Float_t sc_Zcc3  = intLumi*sigma_Zcc3/ev_Zcc3;
  
  cout<<sc_tt<<endl;
  THStack *hs_met_et[nCuts], *hs_met_over_qt[nCuts];

  for(Int_t n=0; n<nCuts; n++)
    {
      ZeeGamma_met_et[n]      -> Scale(sc_ZeeGamma);
      ZeeGamma_met_over_qt[n] -> Scale(sc_ZeeGamma);

      ZllG_met_et[n]      -> Scale(sc_ZllG);
      ZllG_met_over_qt[n] -> Scale(sc_ZllG);

      Wjets_met_et[n]      -> Scale(sc_Wjets);
      Wjets_met_over_qt[n] -> Scale(sc_Wjets);

      tt_met_et[n]      -> Scale(sc_tt);
      tt_met_over_qt[n] -> Scale(sc_tt);

      tT_met_et[n]      -> Scale(sc_tT);
      tT_met_over_qt[n] -> Scale(sc_tT);
      tS_met_et[n]      -> Scale(sc_tS);
      tS_met_over_qt[n] -> Scale(sc_tS);
      tW_met_et[n]      -> Scale(sc_tW);
      tW_met_over_qt[n] -> Scale(sc_tW);

      t_all_met_over_qt[n] = (TH1F*)tS_met_over_qt[n];
      t_all_met_over_qt[n] -> Add(tT_met_over_qt[n]);  
      t_all_met_over_qt[n] -> Add(tW_met_over_qt[n]);
      t_all_met_over_qt[n] -> SetLineColor(kGreen+2);
      t_all_met_et[n]  = (TH1F*)tS_met_et[n];
      t_all_met_et[n]  -> Add(tT_met_et[n]);
      t_all_met_et[n]  -> Add(tW_met_et[n]);

      //t_all_met_phi[n] = (TH1F*)tS_met_phi[n];

      ZZtoAny_met_et[n]      -> Scale(sc_ZZtoAny);
      ZZtoAny_met_over_qt[n] -> Scale(sc_ZZtoAny);

      DYmumu_met_et[n]      -> Scale(sc_DYmumu);
      DYmumu_met_over_qt[n] -> Scale(sc_DYmumu);
      DYee_met_et[n]      -> Scale(sc_DYee);
      DYee_met_over_qt[n] -> Scale(sc_DYee);
      WW_met_et[n]        -> Scale(sc_WW);
      WW_met_over_qt[n]   -> Scale(sc_WW);

      WZ_met_et[n]        -> Scale(sc_WZ);
      WZ_met_over_qt[n]   -> Scale(sc_WZ);

      Zbb0_met_et[n]      -> Scale(sc_Zbb0);
      Zbb0_met_over_qt[n] -> Scale(sc_Zbb0);
      Zbb1_met_et[n]      -> Scale(sc_Zbb1);
      Zbb1_met_over_qt[n] -> Scale(sc_Zbb1);
      Zbb2_met_et[n]      -> Scale(sc_Zbb2);
      Zbb2_met_over_qt[n] -> Scale(sc_Zbb2);
      Zbb3_met_et[n]      -> Scale(sc_Zbb3);
      Zbb3_met_over_qt[n] -> Scale(sc_Zbb3);
      Zcc0_met_et[n]      -> Scale(sc_Zcc0);
      Zcc0_met_over_qt[n] -> Scale(sc_Zcc0);
      Zcc1_met_et[n]      -> Scale(sc_Zcc1);
      Zcc1_met_over_qt[n] -> Scale(sc_Zcc1);
      Zcc2_met_et[n]      -> Scale(sc_Zcc2);
      Zcc2_met_over_qt[n] -> Scale(sc_Zcc2);
      Zcc3_met_et[n]      -> Scale(sc_Zcc3);
      Zcc3_met_over_qt[n] -> Scale(sc_Zcc3);
     
      /*
	Zjets_met_et[n] = (TH1F*)Zbb0_met_et[n];
      Zjets_met_et[n] -> Add(Zbb1_met_et[n]);  
      Zjets_met_et[n] -> Add(Zbb2_met_et[n]);
      Zjets_met_et[n] -> Add(Zbb3_met_et[n]);
      Zjets_met_et[n] -> Add(Zcc0_met_et[n]);
      Zjets_met_et[n] -> Add(Zcc1_met_et[n]);
      Zjets_met_et[n] -> Add(Zcc2_met_et[n]);
      Zjets_met_et[n] -> Add(Zcc3_met_et[n]);
      Zjets_met_over_qt[n] = (TH1F*)Zbb0_met_over_qt[n];
      Zjets_met_over_qt[n] -> Add(Zbb1_met_over_qt[n]);  
      Zjets_met_over_qt[n] -> Add(Zbb2_met_over_qt[n]);
      Zjets_met_over_qt[n] -> Add(Zbb3_met_over_qt[n]);
      Zjets_met_over_qt[n] -> Add(Zcc0_met_over_qt[n]);
      Zjets_met_over_qt[n] -> Add(Zcc1_met_over_qt[n]);
      Zjets_met_over_qt[n] -> Add(Zcc2_met_over_qt[n]);
      Zjets_met_over_qt[n] -> Add(Zcc3_met_over_qt[n]);
      */

      if(sel==1)  
	{
	  Zjets_met_et[n] = (TH1F*)DYmumu_met_et[n];
	  Zjets_met_over_qt[n] = (TH1F*)DYmumu_met_over_qt[n];
	  //Zjets_met_et[n] -> Add(DYmumu_met_et[n]);
	  //Zjets_met_over_qt[n] -> Add(DYmumu_met_over_qt[n]);
	}
      if(sel==2) 
	{
	  Zjets_met_et[n] = (TH1F*)DYee_met_et[n];
	  Zjets_met_over_qt[n] = (TH1F*)DYee_met_over_qt[n];
	  //Zjets_met_et[n] -> Add(DYee_met_et[n]);
	  //Zjets_met_over_qt[n] -> Add(DYee_met_over_qt[n]);
	}
      
      Zjets_met_et[n]      -> SetFillColor(kRed+1);
      Zjets_met_et[n]      -> SetLineColor(kGreen+2);
      Zjets_met_over_qt[n] -> SetFillColor(kRed+1);
      Zjets_met_over_qt[n] -> SetLineColor(kGreen+2);

      ZeeGamma_met_et[n] -> SetFillColor(kYellow+1);
      //ZeeGamma_met_et[n] -> SetLineColor(kBlue);
      ZeeGamma_met_over_qt[n] -> SetFillColor(kYellow+1);
      //ZeeGamma_met_over_qt[n] -> SetLineColor(kBlue);

      ZllG_met_et[n]      -> SetLineColor(kYellow);
      ZllG_met_over_qt[n] -> SetLineColor(kYellow);
      Wjets_met_et[n]     -> SetFillColor(kCyan+1);
      Wjets_met_et[n]     -> SetLineColor(kGreen-2);
      Wjets_met_over_qt[n] -> SetFillColor(kCyan);
      Wjets_met_over_qt[n] -> SetLineColor(kGreen-2);
      tt_met_over_qt[n]    -> SetFillColor(kGreen+2);
      tt_met_over_qt[n]    -> SetLineColor(kOrange+1);
      tt_met_et[n]         -> SetFillColor(kGreen+2);
      tt_met_et[n]         -> SetLineColor(kOrange+1);
      t_all_met_et[n]      -> SetFillColor(kOrange-3);
      t_all_met_et[n]      -> SetLineColor(kBlue);
      t_all_met_over_qt[n] -> SetFillColor(kOrange-3);
      t_all_met_over_qt[n] -> SetLineColor(kBlue);
      ZZtoAny_met_et[n]      -> SetFillColor(kMagenta);
      ZZtoAny_met_et[n]      -> SetLineColor(kBlue);
      ZZtoAny_met_over_qt[n] -> SetFillColor(kMagenta);
      ZZtoAny_met_over_qt[n] -> SetLineColor(kBlue);

      WW_met_et[n]        -> SetFillColor(kSpring);
      WW_met_et[n]        -> SetLineColor(kBlue-2);
      WW_met_over_qt[n]   -> SetFillColor(kSpring);
      WW_met_over_qt[n]   -> SetLineColor(kBlue-2);
      WZ_met_et[n]        -> SetFillColor(kBlue-2);
      WZ_met_et[n]        -> SetLineColor(kRed+2);
      WZ_met_over_qt[n]   -> SetFillColor(kBlue-2);
      WZ_met_over_qt[n]   -> SetLineColor(kRed+2);


      hs_met_et[n] = new THStack(Form("hs_met_et_%i",n),"Stacked met");
      hs_met_et[n]->Add(t_all_met_et[n]);
      hs_met_et[n]->Add(tt_met_et[n]);
        if(sel==1)
	{
	}
      //hs_met_et[n]->Add(ZllG_met_et[n]);
      hs_met_et[n]->Add(WZ_met_et[n]);
      hs_met_et[n]->Add(WW_met_et[n]);
      hs_met_et[n]->Add(ZZtoAny_met_et[n]);
      //hs_met_et[n]->Add(ZllG_met_et[n]);
      hs_met_et[n]->Add(Zjets_met_et[n]);
      hs_met_et[n]->Add(Wjets_met_et[n]);
      //hs_met_et[n]->Add(_met_et[n]);

      hs_met_over_qt[n] = new THStack(Form("hs_met_over_qt_%i",n),"Stacked met");
      hs_met_over_qt[n]->Add(t_all_met_over_qt[n]);
      hs_met_over_qt[n]->Add(tt_met_over_qt[n]);
      if(sel==1)
	{
	}

      //hs_met_over_qt[n]->Add(ZllG_met_over_qt[n]);
      hs_met_over_qt[n]->Add(WZ_met_over_qt[n]);
      hs_met_over_qt[n]->Add(WW_met_over_qt[n]);
      hs_met_over_qt[n]->Add(ZZtoAny_met_over_qt[n]);
      //hs_met_over_qt[n]->Add(ZllG_met_over_qt[n]);
      hs_met_over_qt[n]->Add(Zjets_met_over_qt[n]);
      hs_met_over_qt[n]->Add(Wjets_met_over_qt[n]);
      //hs_met_over_qt[n]->Add(_met_over_qt[n]);


      ggH400_met_et[n]      -> Scale(sc_ggH400);
      ggH400_met_over_qt[n] -> Scale(sc_ggH400);

      ggH400_met_et[n]      -> SetLineColor(kCyan);
      ggH400_met_over_qt[n] -> SetLineColor(kCyan);


      data_met_over_qt[n]   -> SetMarkerStyle(20);
      data_met_et[n]        -> SetMarkerStyle(20);
      data_met_phi[n]       -> SetMarkerStyle(20);
      /*
      */
      //_met_et[n]      -> Scale(sc_);
      //_met_over_qt[n] -> Scale(sc_);

      // _met_et[n]      -> Scale(sc_);
      //_met_over_qt[n] -> Scale(sc_);

    }


  hs_met_et[F0] -> Draw();
  hs_met_et[F0]   -> SetTitle("; pfMET; Events");
  if(metType=="met3")  hs_met_et[F0]   -> SetTitle("; projMET; Events");
  if(metType=="met4")  hs_met_et[F0]   -> SetTitle("; PU corr projMET; Events");
  hs_met_et[F0] -> SetMinimum(0.001);
  data_met_et[F0] -> Draw("same e1pl");
  //data_met_et[F0] -> SetMarkerStyle(20);
  ggH400_met_et[F0] -> Draw("same");
  ZllG_met_et[F0] -> Draw("same");

  leg01 = new TLegend(0.53,0.7,0.95,0.95);
  leg01 -> SetNColumns(2);
  leg01 -> SetTextSize(0.04);
  leg01->AddEntry(data_met_et[F0],     "Data","epl");
  leg01->AddEntry(Zjets_met_et[F0],    "Z + jets","f");
  leg01->AddEntry(ZZtoAny_met_et[F0],  "ZZ","f");
  leg01->AddEntry(Wjets_met_et[F0],    "W + jets","f");
  leg01->AddEntry(WW_met_et[F0],       "WW","f");
  leg01->AddEntry(t_all_met_et[F0],    "t#rightarrow l#nub","f");
  leg01->AddEntry(WZ_met_et[F0],       "WZ","f");
  leg01->AddEntry(tt_met_et[F0],       "tt#rightarrow 2l2#nu2b","f");
  leg01->AddEntry(ZllG_met_et[F0], "Z#gamma#rightarrow ll#gamma","f");
  //if(sel==1) leg01->AddEntry(ZmumuGamma_met_et[F0], "Z#gamma#rightarrow #mu#mu#gamma","f");

  leg01->AddEntry(ggH400_met_et[F0],   "10xH400","f");
  leg01->SetFillColor(kWhite);
  leg01 -> Draw();

  c1 -> SetLogy(1);
  c1 -> SaveAs(imgpath+"ov01.png");


  hs_met_over_qt[F0] -> Draw();
  hs_met_over_qt[F0] -> SetTitle("; pfMET/q_{T}; Events");
  if(metType=="met3")  hs_met_over_qt[F0]   -> SetTitle("; projMET/q_{T}; Events");
  if(metType=="met4")  hs_met_over_qt[F0]   -> SetTitle("; PU corr projMET/q_{T}; Events");
  hs_met_over_qt[F0] -> SetMinimum(0.001);
  data_met_over_qt[F0] -> Draw("same e1pl");
  //data_met_over_qt[F0] -> SetMarkerStyle(20);
  ggH400_met_over_qt[F0] -> Draw("same");
  ZllG_met_over_qt[F0] -> Draw("same");
  leg01 -> Draw();
  c1 -> SaveAs(imgpath+"ov02.png");

  hs_met_et[F1]     -> Draw();
  hs_met_et[F1]     -> SetTitle("; pfMET; Events");
  if(metType=="met3")  hs_met_et[F1]   -> SetTitle("; projMET; Events");
  if(metType=="met4")  hs_met_et[F1]   -> SetTitle("; PU corr projMET; Events");
  hs_met_et[F1]     -> SetMinimum(0.001);
  hs_met_et[F1]     -> SetMaximum(100000);
  data_met_et[F1]   -> Draw("same e1pl");
  ggH400_met_et[F1] -> Draw("same");
  ZllG_met_et[F1] -> Draw("same");
  leg01 -> Draw(); 
  c1 -> SaveAs(imgpath+"ov03.png");


  hs_met_over_qt[F1] -> Draw();
  hs_met_over_qt[F1] -> SetTitle("; pfMET/q_{T}; Events");
  if(metType=="met3")  hs_met_over_qt[F1]   -> SetTitle("; projMET/q_{T}; Events");
  if(metType=="met4")  hs_met_over_qt[F1]   -> SetTitle("; PU corr projMET/q_{T}; Events");
  hs_met_over_qt[F1] -> SetMinimum(0.001);
  hs_met_over_qt[F1] -> SetMaximum(100000);
  data_met_over_qt[F1] -> Draw("same e1pl");
  //data_met_over_qt[F1] -> SetMarkerStyle(20);
  ggH400_met_over_qt[F1] -> Draw("same");
  ZllG_met_over_qt[F1] -> Draw("same");
  leg01 -> Draw();
  c1 -> SaveAs(imgpath+"ov04.png");

  hs_met_et[F2]     -> Draw();
  hs_met_et[F2]     -> SetTitle("; pfMET; Events");
  if(metType=="met3")  hs_met_et[F2]   -> SetTitle("; projMET; Events");
  if(metType=="met4")  hs_met_et[F2]   -> SetTitle("; PU corr projMET; Events");
  hs_met_et[F2]     -> SetMinimum(0.001);
  hs_met_et[F2]     -> SetMaximum(100000);
  data_met_et[F2]   -> Draw("same e1pl");
  //data_met_et[F2]   -> SetMarkerStyle(20);
  ggH400_met_et[F2] -> Draw("same");
  ZllG_met_et[F2] -> Draw("same");
  leg01 -> Draw(); 
  c1 -> SaveAs(imgpath+"ov05.png");

  hs_met_over_qt[F2] -> Draw();
  hs_met_over_qt[F2] -> SetTitle("; pfMET/q_{T}; Events");
  if(metType=="met3")  hs_met_over_qt[F2]   -> SetTitle("; projMET/q_{T}; Events");
  if(metType=="met4")  hs_met_over_qt[F2]   -> SetTitle("; PU corr projMET/q_{T}; Events");
  hs_met_over_qt[F2] -> SetMinimum(0.001);
  hs_met_over_qt[F2] -> SetMaximum(100000);
  data_met_over_qt[F2] -> Draw("same e1pl");
  //data_met_over_qt[F2] -> SetMarkerStyle(20);
  ggH400_met_over_qt[F2] -> Draw("same");
  ZllG_met_over_qt[F2] -> Draw("same");
  leg01 -> Draw();
  c1 -> SaveAs(imgpath+"ov06.png");

  /*

  TH1F * temp = (TH1F*) data_met_et[F2]->Clone();
  //TH1F * temp = (TH1F*)hs_met_et[F2] -> GetStack() -> Clone();
  //temp -> SetLineColor(kGray);
  temp -> Draw("hist");
  data_met_et[F2]   -> Draw("same e1pl");
  ggH400_met_et[F2] -> Draw("same");
  // leg01 -> Draw();
  c1 -> SaveAs(imgpath+"ov06.png");

  TCanvas *c2 = new TCanvas("c2","ratios", 0,0, 600,1000);
  c2 -> Divide(1,2);
  c2 -> cd(1);
  c2_1 -> SetLogy();
  hs_met_et[F1] -> Draw("A");
  hs_met_et[F1] -> GetHistogram() -> Draw("hist");
  data_met_et[F1]   -> Draw("same e1pl");
     //  ggH400_met_et[F1] -> Draw("same");
  c2 -> cd(2);
  //temp -> Divide(Zjets_met_et[F1]); 
  temp -> Divide((TH1F*)hs_met_et[F1] -> GetHistogram(), data_met_et[F1]); 
  //temp -> Divide(data_met_et[F1]); 
  temp -> Draw("ep");
  c2 -> SaveAs(imgpath+"ov05.png");

  */
  /*
  ZZtoAny_met_phi[F1]  -> Draw();
  ZZtoAny_met_phi[F1]  -> SetMinimum(0);
  ZZtoAny_met_phi[F1]  -> SetTitle("Met phi; #phi(MET); Events");
  Wjets_met_phi[F1]    -> Draw("same");
  ZeeGamma_met_phi[F1] -> Draw("same"); 
  tt_met_phi[F1]       -> Draw("same");
  ggH400_met_phi[F1]   -> Draw("same");
  data_met_phi[F1]     -> Draw("same");
  leg01 -> Draw();
  */

}
