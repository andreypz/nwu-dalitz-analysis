void plot(){
  // gROOT->ProcessLine(".L ./drawOverflow.C");
  gROOT->ProcessLine(".L ../data/tdrstyle.C");
  setTDRStyle();
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLegendBorderSize(0);
  gROOT->ForceStyle();
  cout.precision(7); cout.setf(ios::fixed, ios::floatfield);
  TH1::SetDefaultSumw2(kTRUE);

#define nCuts 10
#define F0 6
#define F1 5
#define F2 2

  TString histoPath("v06");
  //TString sample("none");
  Int_t sel = 1; //1-muon, 2 -electron

  //Types of met: met - pfMet, met1 - type1 corrected, met2 - pfMet passed Noise filters, 
  //met3 - projMet, met4 - puProj corrected met (those two are passed Noise filters) 
  TString metType("met3");   TString mtType("mt0"); 


  if (sel==1)  TString imgpath("~/afs/public_html/higgs/overview/muon/");  
  if (sel==2)  TString imgpath("~/afs/public_html/higgs/overview/electron/");  
  //TString imgpath(Form("~/afs/public_html/higgs/%s/",sample.Data()));  

  TFile* fmc_ZllG      = new TFile(Form("./%s/dir_%i_MC_ZllG_/hhhh.root", histoPath.Data(), sel));
  //TFile* fmc_ZeeGamma  = new TFile(Form("./%s/dir_%i_MC_ZeeGamma_/hhhh.root", histoPath.Data(), sel));
  TFile* fmc_Wjets     = new TFile(Form("./%s/dir_%i_MC_Wjets_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_ZZtoAny   = new TFile(Form("./%s/dir_%i_MC_ZZtoAny_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_tS        = new TFile(Form("./%s/dir_%i_MC_tSchannel_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_tT        = new TFile(Form("./%s/dir_%i_MC_tTchannel_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_tW        = new TFile(Form("./%s/dir_%i_MC_tWchannel_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_tt        = new TFile(Form("./%s/dir_%i_MC_tt2l2nu2b_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_DYmumu    = new TFile(Form("./%s/dir_%i_MC_DYmumu_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_DYee      = new TFile(Form("./%s/dir_%i_MC_DYee_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_DYtautau  = new TFile(Form("./%s/dir_%i_MC_DYtautau_/hhhh.root",histoPath.Data(), sel));
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

  TFile* fmc_ggH200    = new TFile(Form("./%s/dir_%i_MC_ggH400_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_ggH400    = new TFile(Form("./%s/dir_%i_MC_ggH400_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_ggH600    = new TFile(Form("./%s/dir_%i_MC_ggH400_/hhhh.root",histoPath.Data(), sel));
  //TFile* fmc_  = new TFile(Form("./%s/dir_%i__/hhhh.root",histoPath.Data(), sel));

  TFile* fda_2011A_DoubleMu  = new TFile(Form("./%s/dir_1_2011A_May10_DoubleMu_/hhhh.root",histoPath.Data()));
  TFile* fda_2011A_DoubleEl  = new TFile(Form("./%s/dir_2_2011A_May10_DoubleElectron_/hhhh.root",histoPath.Data()));
  
  TFile *fr = new TFile(Form("./forRadek_%i.root",sel), "RECREATE");

  TH1F *data_met_phi[nCuts], *data_met_et[nCuts], *data_met_over_qt[nCuts], *data_di_qt[nCuts], *data_mt[nCuts];

  TH1F *ZllG_met_phi[nCuts], *ZllG_met_et[nCuts], *ZllG_met_over_qt[nCuts], *ZllG_di_qt[nCuts], *ZllG_mt[nCuts];
  TH1F *Wjets_met_phi[nCuts], *Wjets_met_et[nCuts], *Wjets_met_over_qt[nCuts], *Wjets_di_qt[nCuts], *Wjets_mt[nCuts];
  TH1F *ZZtoAny_met_phi[nCuts], *ZZtoAny_met_et[nCuts], *ZZtoAny_met_over_qt[nCuts], *ZZtoAny_di_qt[nCuts], *ZZtoAny_mt[nCuts];
  TH1F *ggH200_met_phi[nCuts], *ggH200_met_et[nCuts], *ggH200_met_over_qt[nCuts], *ggH200_di_qt[nCuts], *ggH200_mt[nCuts];
  TH1F *ggH400_met_phi[nCuts], *ggH400_met_et[nCuts], *ggH400_met_over_qt[nCuts], *ggH400_di_qt[nCuts], *ggH400_mt[nCuts];
  TH1F *ggH600_met_phi[nCuts], *ggH600_met_et[nCuts], *ggH600_met_over_qt[nCuts], *ggH600_di_qt[nCuts], *ggH600_mt[nCuts];
  TH1F *t_all_met_phi[nCuts], *t_all_met_et[nCuts], *t_all_met_over_qt[nCuts], *t_all_di_qt[nCuts], *t_all_mt[nCuts];
  TH1F *tS_met_phi[nCuts], *tS_met_et[nCuts], *tS_met_over_qt[nCuts], *tS_di_qt[nCuts], *tS_mt[nCuts];
  TH1F *tT_met_phi[nCuts], *tT_met_et[nCuts], *tT_met_over_qt[nCuts], *tT_di_qt[nCuts], *tT_mt[nCuts];
  TH1F *tW_met_phi[nCuts], *tW_met_et[nCuts], *tW_met_over_qt[nCuts], *tW_di_qt[nCuts], *tW_mt[nCuts];
  TH1F *tt_met_phi[nCuts], *tt_met_et[nCuts], *tt_met_over_qt[nCuts], *tt_di_qt[nCuts], *tt_mt[nCuts];
  TH1F *WW_met_phi[nCuts], *WW_met_et[nCuts], *WW_met_over_qt[nCuts], *WW_di_qt[nCuts], *WW_mt[nCuts];
  TH1F *WZ_met_phi[nCuts], *WZ_met_et[nCuts], *WZ_met_over_qt[nCuts], *WZ_di_qt[nCuts], *WZ_mt[nCuts];

  TH2F *data_met_et_ovQt[nCuts], *ZllG_met_et_ovQt[nCuts], *Wjets_met_et_ovQt[nCuts],*ZZtoAny_met_et_ovQt[nCuts],*ggH200_met_et_ovQt[nCuts],*ggH400_met_et_ovQt[nCuts], *ggH600_met_et_ovQt[nCuts],*t_all_met_et_ovQt[nCuts],*tS_met_et_ovQt[nCuts],*tT_met_et_ovQt[nCuts],*tW_met_et_ovQt[nCuts],*tt_met_et_ovQt[nCuts],*WW_met_et_ovQt[nCuts],*WZ_met_et_ovQt[nCuts],*DYmumu_met_et_ovQt[nCuts],*DYee_met_et_ovQt[nCuts],*DYtautau_met_et_ovQt[nCuts],*Zjets_met_et_ovQt[nCuts],*Zbb0_met_et_ovQt[nCuts],*Zbb1_met_et_ovQt[nCuts],*Zbb2_met_et_ovQt[nCuts],*Zbb3_met_et_ovQt[nCuts],*Zcc0_met_et_ovQt[nCuts],*Zcc1_met_et_ovQt[nCuts],*Zcc2_met_et_ovQt[nCuts],*Zcc3_met_et_ovQt[nCuts];
  //,*Z_met_et_ovQt[nCuts],*_met_et_ovQt[nCuts],*_met_et_ovQt[nCuts], 


  TH1F *DYmumu_met_phi[nCuts], *DYmumu_met_et[nCuts], *DYmumu_met_over_qt[nCuts], *DYmumu_di_qt[nCuts], *DYmumu_mt[nCuts];
  TH1F *DYee_met_phi[nCuts], *DYee_met_et[nCuts], *DYee_met_over_qt[nCuts], *DYee_di_qt[nCuts], *DYee_mt[nCuts];
  TH1F *DYtautau_met_phi[nCuts], *DYtautau_met_et[nCuts], *DYtautau_met_over_qt[nCuts], *DYtautau_di_qt[nCuts], *DYtautau_mt[nCuts];
  TH1F *Zjets_met_phi[nCuts], *Zjets_met_et[nCuts], *Zjets_met_over_qt[nCuts], *Zjets_di_qt[nCuts], *Zjets_mt[nCuts];

  TH1F *Zbb0_met_phi[nCuts], *Zbb0_met_et[nCuts], *Zbb0_met_over_qt[nCuts], *Zbb0_di_qt[nCuts], *Zbb0_mt[nCuts];
  TH1F *Zbb1_met_phi[nCuts], *Zbb1_met_et[nCuts], *Zbb1_met_over_qt[nCuts], *Zbb1_di_qt[nCuts], *Zbb1_mt[nCuts];
  TH1F *Zbb2_met_phi[nCuts], *Zbb2_met_et[nCuts], *Zbb2_met_over_qt[nCuts], *Zbb2_di_qt[nCuts], *Zbb2_mt[nCuts];
  TH1F *Zbb3_met_phi[nCuts], *Zbb3_met_et[nCuts], *Zbb3_met_over_qt[nCuts], *Zbb3_di_qt[nCuts], *Zbb3_mt[nCuts];
  TH1F *Zcc0_met_phi[nCuts], *Zcc0_met_et[nCuts], *Zcc0_met_over_qt[nCuts], *Zcc0_di_qt[nCuts], *Zcc0_mt[nCuts];
  TH1F *Zcc1_met_phi[nCuts], *Zcc1_met_et[nCuts], *Zcc1_met_over_qt[nCuts], *Zcc1_di_qt[nCuts], *Zcc1_mt[nCuts];
  TH1F *Zcc2_met_phi[nCuts], *Zcc2_met_et[nCuts], *Zcc2_met_over_qt[nCuts], *Zcc2_di_qt[nCuts], *Zcc2_mt[nCuts];
  TH1F *Zcc3_met_phi[nCuts], *Zcc3_met_et[nCuts], *Zcc3_met_over_qt[nCuts], *Zcc3_di_qt[nCuts], *Zcc3_mt[nCuts];
  //TH1F *_met_phi[nCuts], *_met_et[nCuts], *_met_over_qt[nCuts], *_mt[nCuts];
  //TH1F *_met_phi[nCuts], *_met_et[nCuts], *_met_over_qt[nCuts], *_mt[nCuts];
  //TH1F *_met_phi[nCuts], *_met_et[nCuts], *_met_over_qt[nCuts], *_mt[nCuts];
 

  for(Int_t n=0; n<nCuts; n++)
    {
  
      ZllG_mt[n]          = (TH1F*)fmc_ZllG->Get(Form("%s_%i",mtType.Data(),n));
      ZllG_di_qt[n]       = (TH1F*)fmc_ZllG->Get(Form("di_qt_%i",n));
      ZllG_met_phi[n]     = (TH1F*)fmc_ZllG->Get(Form("%s_phi_%i",metType.Data(),n));
      ZllG_met_et[n]      = (TH1F*)fmc_ZllG->Get(Form("%s_et_%i",metType.Data(),n));
      ZllG_met_over_qt[n] = (TH1F*)fmc_ZllG->Get(Form("%s_over_qt_%i",metType.Data(),n));
      ZllG_met_et_ovQt[n] = (TH2F*)fmc_ZllG->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      Wjets_mt[n]          = (TH1F*)fmc_Wjets->Get(Form("%s_%i",mtType.Data(),n));
      Wjets_di_qt[n]       = (TH1F*)fmc_Wjets->Get(Form("di_qt_%i",n));
      Wjets_met_phi[n]     = (TH1F*)fmc_Wjets->Get(Form("%s_phi_%i",metType.Data(),n));
      Wjets_met_et[n]      = (TH1F*)fmc_Wjets->Get(Form("%s_et_%i",metType.Data(),n));
      Wjets_met_over_qt[n] = (TH1F*)fmc_Wjets->Get(Form("%s_over_qt_%i",metType.Data(),n));
      Wjets_met_et_ovQt[n] = (TH2F*)fmc_Wjets->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      ZZtoAny_mt[n]          = (TH1F*)fmc_ZZtoAny->Get(Form("%s_%i",mtType.Data(),n));
      ZZtoAny_di_qt[n]       = (TH1F*)fmc_ZZtoAny->Get(Form("di_qt_%i",n));
      ZZtoAny_met_phi[n]     = (TH1F*)fmc_ZZtoAny->Get(Form("%s_phi_%i",metType.Data(),n));
      ZZtoAny_met_et[n]      = (TH1F*)fmc_ZZtoAny->Get(Form("%s_et_%i",metType.Data(),n));
      ZZtoAny_met_over_qt[n] = (TH1F*)fmc_ZZtoAny->Get(Form("%s_over_qt_%i",metType.Data(),n));
      ZZtoAny_met_et_ovQt[n] = (TH2F*)fmc_ZZtoAny->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      tS_mt[n]          = (TH1F*)fmc_tS->Get(Form("%s_%i",mtType.Data(),n));
      tS_di_qt[n]       = (TH1F*)fmc_tS->Get(Form("di_qt_%i",n));
      tS_met_phi[n]     = (TH1F*)fmc_tS->Get(Form("%s_phi_%i",metType.Data(),n));
      tS_met_et[n]      = (TH1F*)fmc_tS->Get(Form("%s_et_%i",metType.Data(),n));
      tS_met_over_qt[n] = (TH1F*)fmc_tS->Get(Form("%s_over_qt_%i",metType.Data(),n));
      tS_met_et_ovQt[n] = (TH2F*)fmc_tS->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      tT_mt[n]          = (TH1F*)fmc_tT->Get(Form("%s_%i",mtType.Data(),n));
      tT_di_qt[n]       = (TH1F*)fmc_tT->Get(Form("di_qt_%i",n));
      tT_met_phi[n]     = (TH1F*)fmc_tT->Get(Form("%s_phi_%i",metType.Data(),n));
      tT_met_et[n]      = (TH1F*)fmc_tT->Get(Form("%s_et_%i",metType.Data(),n));
      tT_met_over_qt[n] = (TH1F*)fmc_tT->Get(Form("%s_over_qt_%i",metType.Data(),n));
      tT_met_et_ovQt[n] = (TH2F*)fmc_tT->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      tW_mt[n]          = (TH1F*)fmc_tW->Get(Form("%s_%i",mtType.Data(),n));
      tW_di_qt[n]       = (TH1F*)fmc_tW->Get(Form("di_qt_%i",n));
      tW_met_phi[n]     = (TH1F*)fmc_tW->Get(Form("%s_phi_%i",metType.Data(),n));
      tW_met_et[n]      = (TH1F*)fmc_tW->Get(Form("%s_et_%i",metType.Data(),n));
      tW_met_over_qt[n] = (TH1F*)fmc_tW->Get(Form("%s_over_qt_%i",metType.Data(),n));
      tW_met_et_ovQt[n] = (TH2F*)fmc_tW->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      tt_mt[n]          = (TH1F*)fmc_tt->Get(Form("%s_%i",mtType.Data(),n));
      tt_di_qt[n]       = (TH1F*)fmc_tt->Get(Form("di_qt_%i",n));
      tt_met_phi[n]     = (TH1F*)fmc_tt->Get(Form("%s_phi_%i",metType.Data(),n));
      tt_met_et[n]      = (TH1F*)fmc_tt->Get(Form("%s_et_%i",metType.Data(),n));
      tt_met_over_qt[n] = (TH1F*)fmc_tt->Get(Form("%s_over_qt_%i",metType.Data(),n));
      tt_met_et_ovQt[n] = (TH2F*)fmc_tt->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      ggH200_mt[n]          = (TH1F*)fmc_ggH200->Get(Form("%s_%i",mtType.Data(),n));
      ggH200_di_qt[n]       = (TH1F*)fmc_ggH200->Get(Form("di_qt_%i",n));
      ggH200_met_phi[n]     = (TH1F*)fmc_ggH200->Get(Form("%s_phi_%i",metType.Data(),n));
      ggH200_met_et[n]      = (TH1F*)fmc_ggH200->Get(Form("%s_et_%i",metType.Data(),n));
      ggH200_met_over_qt[n] = (TH1F*)fmc_ggH200->Get(Form("%s_over_qt_%i",metType.Data(),n));
      ggH200_met_et_ovQt[n] = (TH2F*)fmc_ggH200->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      ggH400_mt[n]          = (TH1F*)fmc_ggH400->Get(Form("%s_%i",mtType.Data(),n));
      ggH400_di_qt[n]       = (TH1F*)fmc_ggH400->Get(Form("di_qt_%i",n));
      ggH400_met_phi[n]     = (TH1F*)fmc_ggH400->Get(Form("%s_phi_%i",metType.Data(),n));
      ggH400_met_et[n]      = (TH1F*)fmc_ggH400->Get(Form("%s_et_%i",metType.Data(),n));
      ggH400_met_over_qt[n] = (TH1F*)fmc_ggH400->Get(Form("%s_over_qt_%i",metType.Data(),n));
      ggH400_met_et_ovQt[n] = (TH2F*)fmc_ggH400->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      ggH600_mt[n]          = (TH1F*)fmc_ggH600->Get(Form("%s_%i",mtType.Data(),n));
      ggH600_di_qt[n]       = (TH1F*)fmc_ggH600->Get(Form("di_qt_%i",n));
      ggH600_met_phi[n]     = (TH1F*)fmc_ggH600->Get(Form("%s_phi_%i",metType.Data(),n));
      ggH600_met_et[n]      = (TH1F*)fmc_ggH600->Get(Form("%s_et_%i",metType.Data(),n));
      ggH600_met_over_qt[n] = (TH1F*)fmc_ggH600->Get(Form("%s_over_qt_%i",metType.Data(),n));
      ggH600_met_et_ovQt[n] = (TH2F*)fmc_ggH600->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      if (sel==1)
	{
	  data_mt[n]          = (TH1F*)fda_2011A_DoubleMu->Get(Form("%s_%i",mtType.Data(),n));
	  data_di_qt[n]       = (TH1F*)fda_2011A_DoubleMu->Get(Form("di_qt_%i",n));
	  data_met_phi[n]     = (TH1F*)fda_2011A_DoubleMu->Get(Form("%s_phi_%i",metType.Data(),n));
	  data_met_et[n]      = (TH1F*)fda_2011A_DoubleMu->Get(Form("%s_et_%i",metType.Data(),n));
	  data_met_over_qt[n] = (TH1F*)fda_2011A_DoubleMu->Get(Form("%s_over_qt_%i",metType.Data(),n));
	  data_met_et_ovQt[n] = (TH2F*)fda_2011A_DoubleMu->Get(Form("%s_et_ovQt_%i",metType.Data(),n));
	}
      if(sel==2)
	{
	  data_mt[n]          = (TH1F*)fda_2011A_DoubleEl->Get(Form("%s_%i",mtType.Data(),n));
	  data_di_qt[n]       = (TH1F*)fda_2011A_DoubleEl->Get(Form("di_qt_%i",n));
	  data_met_phi[n]     = (TH1F*)fda_2011A_DoubleEl->Get(Form("%s_phi_%i",metType.Data(),n));
	  data_met_et[n]      = (TH1F*)fda_2011A_DoubleEl->Get(Form("%s_et_%i",metType.Data(),n));
	  data_met_over_qt[n] = (TH1F*)fda_2011A_DoubleEl->Get(Form("%s_over_qt_%i",metType.Data(),n));
	  data_met_et_ovQt[n] = (TH2F*)fda_2011A_DoubleEl->Get(Form("%s_et_ovQt_%i",metType.Data(),n));
	}

      DYmumu_mt[n]          = (TH1F*)fmc_DYmumu->Get(Form("%s_%i",mtType.Data(),n));
      DYmumu_di_qt[n]       = (TH1F*)fmc_DYmumu->Get(Form("di_qt_%i",n));
      DYmumu_met_phi[n]     = (TH1F*)fmc_DYmumu->Get(Form("%s_phi_%i",metType.Data(),n));
      DYmumu_met_et[n]      = (TH1F*)fmc_DYmumu->Get(Form("%s_et_%i",metType.Data(),n));
      DYmumu_met_over_qt[n] = (TH1F*)fmc_DYmumu->Get(Form("%s_over_qt_%i",metType.Data(),n));
      DYmumu_met_et_ovQt[n] = (TH2F*)fmc_DYmumu->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      DYee_mt[n]          = (TH1F*)fmc_DYee->Get(Form("%s_%i",mtType.Data(),n));
      DYee_di_qt[n]       = (TH1F*)fmc_DYee->Get(Form("di_qt_%i",n));
      DYee_met_phi[n]     = (TH1F*)fmc_DYee->Get(Form("%s_phi_%i",metType.Data(),n));
      DYee_met_et[n]      = (TH1F*)fmc_DYee->Get(Form("%s_et_%i",metType.Data(),n));
      DYee_met_over_qt[n] = (TH1F*)fmc_DYee->Get(Form("%s_over_qt_%i",metType.Data(),n));
      DYee_met_et_ovQt[n] = (TH2F*)fmc_DYee->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      DYtautau_mt[n]          = (TH1F*)fmc_DYtautau->Get(Form("%s_%i",mtType.Data(),n));
      DYtautau_di_qt[n]       = (TH1F*)fmc_DYtautau->Get(Form("di_qt_%i",n));
      DYtautau_met_phi[n]     = (TH1F*)fmc_DYtautau->Get(Form("%s_phi_%i",metType.Data(),n));
      DYtautau_met_et[n]      = (TH1F*)fmc_DYtautau->Get(Form("%s_et_%i",metType.Data(),n));
      DYtautau_met_over_qt[n] = (TH1F*)fmc_DYtautau->Get(Form("%s_over_qt_%i",metType.Data(),n));
      DYtautau_met_et_ovQt[n] = (TH2F*)fmc_DYtautau->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      WW_mt[n]            = (TH1F*)fmc_WW->Get(Form("%s_%i",mtType.Data(),n));
      WW_di_qt[n]         = (TH1F*)fmc_WW->Get(Form("di_qt_%i",n));
      WW_met_phi[n]       = (TH1F*)fmc_WW->Get(Form("%s_phi_%i",metType.Data(),n));
      WW_met_et[n]        = (TH1F*)fmc_WW->Get(Form("%s_et_%i",metType.Data(),n));
      WW_met_over_qt[n]   = (TH1F*)fmc_WW->Get(Form("%s_over_qt_%i",metType.Data(),n));
      WW_met_et_ovQt[n]   = (TH2F*)fmc_WW->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      WZ_mt[n]            = (TH1F*)fmc_WZ->Get(Form("%s_%i",mtType.Data(),n));
      WZ_di_qt[n]         = (TH1F*)fmc_WZ->Get(Form("di_qt_%i",n));
      WZ_met_phi[n]       = (TH1F*)fmc_WZ->Get(Form("%s_phi_%i",metType.Data(),n));
      WZ_met_et[n]        = (TH1F*)fmc_WZ->Get(Form("%s_et_%i",metType.Data(),n));
      WZ_met_over_qt[n]   = (TH1F*)fmc_WZ->Get(Form("%s_over_qt_%i",metType.Data(),n));
      WZ_met_et_ovQt[n]   = (TH2F*)fmc_WZ->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      Zbb0_mt[n]          = (TH1F*)fmc_Zbb0->Get(Form("%s_%i",mtType.Data(),n));
      Zbb0_di_qt[n]       = (TH1F*)fmc_Zbb0->Get(Form("di_qt_%i",n));
      Zbb0_met_phi[n]     = (TH1F*)fmc_Zbb0->Get(Form("%s_phi_%i",metType.Data(),n));
      Zbb0_met_et[n]      = (TH1F*)fmc_Zbb0->Get(Form("%s_et_%i",metType.Data(),n));
      Zbb0_met_over_qt[n] = (TH1F*)fmc_Zbb0->Get(Form("%s_over_qt_%i",metType.Data(),n));
      Zbb0_met_et_ovQt[n] = (TH2F*)fmc_Zbb0->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      Zbb1_mt[n]          = (TH1F*)fmc_Zbb1->Get(Form("%s_%i",mtType.Data(),n));
      Zbb1_di_qt[n]       = (TH1F*)fmc_Zbb1->Get(Form("di_qt_%i",n));
      Zbb1_met_phi[n]     = (TH1F*)fmc_Zbb1->Get(Form("%s_phi_%i",metType.Data(),n));
      Zbb1_met_et[n]      = (TH1F*)fmc_Zbb1->Get(Form("%s_et_%i",metType.Data(),n));
      Zbb1_met_over_qt[n] = (TH1F*)fmc_Zbb1->Get(Form("%s_over_qt_%i",metType.Data(),n));
      Zbb1_met_et_ovQt[n] = (TH2F*)fmc_Zbb1->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      Zbb2_mt[n]          = (TH1F*)fmc_Zbb2->Get(Form("%s_%i",mtType.Data(),n));
      Zbb2_di_qt[n]       = (TH1F*)fmc_Zbb2->Get(Form("di_qt_%i",n));
      Zbb2_met_phi[n]     = (TH1F*)fmc_Zbb2->Get(Form("%s_phi_%i",metType.Data(),n));
      Zbb2_met_et[n]      = (TH1F*)fmc_Zbb2->Get(Form("%s_et_%i",metType.Data(),n));
      Zbb2_met_over_qt[n] = (TH1F*)fmc_Zbb2->Get(Form("%s_over_qt_%i",metType.Data(),n));
      Zbb2_met_et_ovQt[n] = (TH2F*)fmc_Zbb2->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      Zbb3_mt[n]          = (TH1F*)fmc_Zbb3->Get(Form("%s_%i",mtType.Data(),n));
      Zbb3_di_qt[n]       = (TH1F*)fmc_Zbb3->Get(Form("di_qt_%i",n));
      Zbb3_met_phi[n]     = (TH1F*)fmc_Zbb3->Get(Form("%s_phi_%i",metType.Data(),n));
      Zbb3_met_et[n]      = (TH1F*)fmc_Zbb3->Get(Form("%s_et_%i",metType.Data(),n));
      Zbb3_met_over_qt[n] = (TH1F*)fmc_Zbb3->Get(Form("%s_over_qt_%i",metType.Data(),n));
      Zbb3_met_et_ovQt[n] = (TH2F*)fmc_Zbb3->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      Zcc0_mt[n]          = (TH1F*)fmc_Zcc0->Get(Form("%s_%i",mtType.Data(),n));
      Zcc0_di_qt[n]       = (TH1F*)fmc_Zcc0->Get(Form("di_qt_%i",n));
      Zcc0_met_phi[n]     = (TH1F*)fmc_Zcc0->Get(Form("%s_phi_%i",metType.Data(),n));
      Zcc0_met_et[n]      = (TH1F*)fmc_Zcc0->Get(Form("%s_et_%i",metType.Data(),n));
      Zcc0_met_over_qt[n] = (TH1F*)fmc_Zcc0->Get(Form("%s_over_qt_%i",metType.Data(),n));
      Zcc0_met_et_ovQt[n] = (TH2F*)fmc_Zcc0->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      Zcc1_mt[n]          = (TH1F*)fmc_Zcc1->Get(Form("%s_%i",mtType.Data(),n));
      Zcc1_di_qt[n]       = (TH1F*)fmc_Zcc1->Get(Form("di_qt_%i",n));
      Zcc1_met_phi[n]     = (TH1F*)fmc_Zcc1->Get(Form("%s_phi_%i",metType.Data(),n));
      Zcc1_met_et[n]      = (TH1F*)fmc_Zcc1->Get(Form("%s_et_%i",metType.Data(),n));
      Zcc1_met_over_qt[n] = (TH1F*)fmc_Zcc1->Get(Form("%s_over_qt_%i",metType.Data(),n));
      Zcc1_met_et_ovQt[n] = (TH2F*)fmc_Zcc1->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      Zcc2_mt[n]          = (TH1F*)fmc_Zcc2->Get(Form("%s_%i",mtType.Data(),n));
      Zcc2_di_qt[n]       = (TH1F*)fmc_Zcc2->Get(Form("di_qt_%i",n));
      Zcc2_met_phi[n]     = (TH1F*)fmc_Zcc2->Get(Form("%s_phi_%i",metType.Data(),n));
      Zcc2_met_et[n]      = (TH1F*)fmc_Zcc2->Get(Form("%s_et_%i",metType.Data(),n));
      Zcc2_met_over_qt[n] = (TH1F*)fmc_Zcc2->Get(Form("%s_over_qt_%i",metType.Data(),n));
      Zcc2_met_et_ovQt[n] = (TH2F*)fmc_Zcc2->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      Zcc3_mt[n]          = (TH1F*)fmc_Zbb3->Get(Form("%s_%i",mtType.Data(),n));
      Zcc3_di_qt[n]       = (TH1F*)fmc_Zcc3->Get(Form("di_qt_%i",n));
      Zcc3_met_phi[n]     = (TH1F*)fmc_Zcc3->Get(Form("%s_phi_%i",metType.Data(),n));
      Zcc3_met_et[n]      = (TH1F*)fmc_Zcc3->Get(Form("%s_et_%i",metType.Data(),n));
      Zcc3_met_over_qt[n] = (TH1F*)fmc_Zcc3->Get(Form("%s_over_qt_%i",metType.Data(),n));
      Zcc3_met_et_ovQt[n] = (TH2F*)fmc_Zcc3->Get(Form("%s_et_ovQt_%i",metType.Data(),n));

      /*

      _met_phi[n]     = (TH1F*)fmc_->Get(Form("%s_phi_%i",metType.Data(),n));
      _met_et[n]      = (TH1F*)fmc_->Get(Form("%s_et_%i",metType.Data(),n));
      _met_over_qt[n] = (TH1F*)fmc_->Get(Form("%s_over_qt_%i",metType.Data(),n));
      */

    }


  Float_t intLumi = 191.0;
  Float_t ev_ZllG     = 165000.;
  Float_t ev_Wjets    = Wjets_met_et[0]    -> GetEntries();
  Float_t ev_tt       = tt_met_et[0]       -> GetEntries();
  Float_t ev_tS       = tS_met_et[0]       -> GetEntries();
  Float_t ev_tT       = tT_met_et[0]       -> GetEntries();
  Float_t ev_tW       = tW_met_et[0]       -> GetEntries();
  Float_t ev_ZZtoAny  = ZZtoAny_met_et[0]  -> GetEntries();
  Float_t ev_DYmumu   = 1984154.;
  Float_t ev_DYee     = 1992276.;
  Float_t ev_DYtautau = 1995369.;
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
  Float_t ev_ggH200   = 109274; //ok crab_0_110619_122140/
  Float_t ev_ggH400   = ggH400_met_et[0]   -> GetEntries();
  Float_t ev_ggH600   = 109255; //ok crab_0_110619_122331/ 


  cout<<"ZllG: "<<ev_ZllG<<"  Wjets: "<<ev_Wjets<<"  ZZtoAny: "<<ev_ZZtoAny<<" H400: "<<ev_ggH400<<endl;
  cout<<"tt: "<<ev_tt<<"    tT: "<<ev_tT<<"    tS: "<<ev_tS<<"    tW: "<<ev_tW<<endl;
 
  //Float_t sigma_ZeeGamma  = 33.73; 

  Float_t sigma_ggH200    = 0.06277;
  Float_t sigma_ggH400    = 0.02229;
  Float_t sigma_ggH600    = 0.004232;

  Float_t sigma_ZllG      = 14.3;
  Float_t sigma_Wjets     = 31314.0;
  Float_t sigma_tt        = 16.5;
  Float_t sigma_tS        = 1.36;
  Float_t sigma_tT        = 20.9;
  Float_t sigma_tW        = 10.6;
  Float_t sigma_ZZtoAny   = 5.9;
  Float_t sigma_DYmumu    = 1666.;
  Float_t sigma_DYee      = 1666.;
  Float_t sigma_DYtautau  = 1666.;
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

  Float_t sc_ZllG      = intLumi*sigma_ZllG/ev_ZllG;
  // Float_t sc_ZeeGamma  = intLumi*sigma_ZeeGamma/ev_ZeeGamma;
  Float_t sc_Wjets     = intLumi*sigma_Wjets/ev_Wjets;
  Float_t sc_tt        = intLumi*sigma_tt/ev_tt;
  Float_t sc_tS        = intLumi*sigma_tS/ev_tS;
  Float_t sc_tT        = intLumi*sigma_tT/ev_tT;
  Float_t sc_tW        = intLumi*sigma_tW/ev_tW;
  Float_t sc_ZZtoAny   = intLumi*sigma_ZZtoAny/ev_ZZtoAny;
  Float_t sc_DYmumu    = intLumi*sigma_DYmumu/ev_DYmumu;  
  Float_t sc_DYee      = intLumi*sigma_DYee/ev_DYee;  
  Float_t sc_DYtautau  = intLumi*sigma_DYtautau/ev_DYtautau;  
  Float_t sc_WW        = intLumi*sigma_WW/ev_WW;  
  Float_t sc_WZ        = intLumi*sigma_WZ/ev_WZ;  
  Float_t sc_ggH200    = 10*intLumi*sigma_ggH200/ev_ggH200;
  Float_t sc_ggH400    = 10*intLumi*sigma_ggH400/ev_ggH400;
  Float_t sc_ggH600    = 10*intLumi*sigma_ggH600/ev_ggH600;
  
  Float_t sc_Zbb0  = intLumi*sigma_Zbb0/ev_Zbb0;
  Float_t sc_Zbb1  = intLumi*sigma_Zbb1/ev_Zbb1;
  Float_t sc_Zbb2  = intLumi*sigma_Zbb2/ev_Zbb2;
  Float_t sc_Zbb3  = intLumi*sigma_Zbb3/ev_Zbb3;
  Float_t sc_Zcc0  = intLumi*sigma_Zcc0/ev_Zcc0;
  Float_t sc_Zcc1  = intLumi*sigma_Zcc1/ev_Zcc1;
  Float_t sc_Zcc2  = intLumi*sigma_Zcc2/ev_Zcc2;
  Float_t sc_Zcc3  = intLumi*sigma_Zcc3/ev_Zcc3;
  
  
 //cout<<sc_tt<<endl;
  THStack *hs_met_et[nCuts], *hs_met_over_qt[nCuts], *hs_di_qt[nCuts], *hs_met_et_ovQt[nCuts], *hs_mt[nCuts];
  TH1F *sum_met_et[nCuts], *sum_met_over_qt[nCuts], *sum_di_qt[nCuts], *sum_mt[nCuts];
  TH2F *sum_met_et_ovQt[nCuts];

  for(Int_t n=0; n<nCuts; n++)
    {
      if(histoPath!="v07" && histoPath!="v08")
	{
      ZllG_mt[n]          -> Scale(sc_ZllG);
      ZllG_di_qt[n]       -> Scale(sc_ZllG);
      ZllG_met_et[n]      -> Scale(sc_ZllG);
      ZllG_met_over_qt[n] -> Scale(sc_ZllG);
      ZllG_met_et_ovQt[n] -> Scale(sc_ZllG);

      Wjets_mt[n]          -> Scale(sc_Wjets);
      Wjets_di_qt[n]       -> Scale(sc_Wjets);
      Wjets_met_et[n]      -> Scale(sc_Wjets);
      Wjets_met_over_qt[n] -> Scale(sc_Wjets);
      Wjets_met_et_ovQt[n] -> Scale(sc_Wjets);

      tt_mt[n]          -> Scale(sc_tt);
      tt_di_qt[n]       -> Scale(sc_tt);
      tt_met_et[n]      -> Scale(sc_tt);
      tt_met_over_qt[n] -> Scale(sc_tt);
      tt_met_et_ovQt[n] -> Scale(sc_tt);

      tT_mt[n]          -> Scale(sc_tT);
      tT_di_qt[n]       -> Scale(sc_tT);
      tT_met_et[n]      -> Scale(sc_tT);
      tT_met_over_qt[n] -> Scale(sc_tT);
      tT_met_et_ovQt[n] -> Scale(sc_tT);
      tS_mt[n]          -> Scale(sc_tS);
      tS_di_qt[n]       -> Scale(sc_tS);
      tS_met_et[n]      -> Scale(sc_tS);
      tS_met_over_qt[n] -> Scale(sc_tS);
      tS_met_et_ovQt[n] -> Scale(sc_tS);
      tW_mt[n]          -> Scale(sc_tW);
      tW_di_qt[n]       -> Scale(sc_tW);
      tW_met_et[n]      -> Scale(sc_tW);
      tW_met_over_qt[n] -> Scale(sc_tW);
      tW_met_et_ovQt[n] -> Scale(sc_tW);
	}
      t_all_mt[n]  = (TH1F*)tS_mt[n];
      t_all_mt[n]  -> Add(tT_mt[n]);
      t_all_mt[n]  -> Add(tW_mt[n]);

      t_all_di_qt[n]  = (TH1F*)tS_di_qt[n];
      t_all_di_qt[n]  -> Add(tT_di_qt[n]);
      t_all_di_qt[n]  -> Add(tW_di_qt[n]);
      
      t_all_met_over_qt[n] = (TH1F*)tS_met_over_qt[n];
      t_all_met_over_qt[n] -> Add(tT_met_over_qt[n]);  
      t_all_met_over_qt[n] -> Add(tW_met_over_qt[n]);

      t_all_met_et[n]  = (TH1F*)tS_met_et[n];
      t_all_met_et[n]  -> Add(tT_met_et[n]);
      t_all_met_et[n]  -> Add(tW_met_et[n]);

      t_all_met_et_ovQt[n]  = (TH2F*)tS_met_et_ovQt[n];
      t_all_met_et_ovQt[n]  -> Add(tT_met_et_ovQt[n]);
      t_all_met_et_ovQt[n]  -> Add(tW_met_et_ovQt[n]);

      //t_all_met_phi[n] = (TH1F*)tS_met_phi[n];


      if(histoPath!="v07" && histoPath!="v08")
	{
      ZZtoAny_mt[n]          -> Scale(sc_ZZtoAny);
      ZZtoAny_di_qt[n]       -> Scale(sc_ZZtoAny);
      ZZtoAny_met_et[n]      -> Scale(sc_ZZtoAny);
      ZZtoAny_met_over_qt[n] -> Scale(sc_ZZtoAny);
      ZZtoAny_met_et_ovQt[n] -> Scale(sc_ZZtoAny);


      DYmumu_mt[n]          -> Scale(sc_DYmumu);
      DYmumu_di_qt[n]       -> Scale(sc_DYmumu);
      DYmumu_met_et[n]      -> Scale(sc_DYmumu);
      DYmumu_met_over_qt[n] -> Scale(sc_DYmumu);
      DYmumu_met_et_ovQt[n] -> Scale(sc_DYmumu);
      DYee_mt[n]            -> Scale(sc_DYee);
      DYee_di_qt[n]         -> Scale(sc_DYee);
      DYee_met_et[n]        -> Scale(sc_DYee);
      DYee_met_over_qt[n]   -> Scale(sc_DYee);
      DYee_met_et_ovQt[n]   -> Scale(sc_DYee);
      DYtautau_mt[n]          -> Scale(sc_DYtautau);
      DYtautau_di_qt[n]       -> Scale(sc_DYtautau);
      DYtautau_met_et[n]      -> Scale(sc_DYtautau);
      DYtautau_met_over_qt[n] -> Scale(sc_DYtautau);
      DYtautau_met_et_ovQt[n] -> Scale(sc_DYtautau);

      WW_mt[n]            -> Scale(sc_WW);
      WW_di_qt[n]         -> Scale(sc_WW);
      WW_met_et[n]        -> Scale(sc_WW);
      WW_met_over_qt[n]   -> Scale(sc_WW);
      WW_met_et_ovQt[n]   -> Scale(sc_WW);

      WZ_mt[n]            -> Scale(sc_WZ);
      WZ_di_qt[n]         -> Scale(sc_WZ);
      WZ_met_et[n]        -> Scale(sc_WZ);
      WZ_met_over_qt[n]   -> Scale(sc_WZ);
      WZ_met_et_ovQt[n]   -> Scale(sc_WZ);

      Zbb0_mt[n]          -> Scale(sc_Zbb0);
      Zbb0_di_qt[n]       -> Scale(sc_Zbb0);
      Zbb0_met_et[n]      -> Scale(sc_Zbb0);
      Zbb0_met_over_qt[n] -> Scale(sc_Zbb0);
      Zbb0_met_et_ovQt[n] -> Scale(sc_Zbb0);
      Zbb1_mt[n]          -> Scale(sc_Zbb1);
      Zbb1_di_qt[n]       -> Scale(sc_Zbb1);
      Zbb1_met_et[n]      -> Scale(sc_Zbb1);
      Zbb1_met_over_qt[n] -> Scale(sc_Zbb1);
      Zbb1_met_et_ovQt[n] -> Scale(sc_Zbb1);
      Zbb2_mt[n]          -> Scale(sc_Zbb2);
      Zbb2_di_qt[n]       -> Scale(sc_Zbb2);
      Zbb2_met_et[n]      -> Scale(sc_Zbb2);
      Zbb2_met_over_qt[n] -> Scale(sc_Zbb2);
      Zbb2_met_et_ovQt[n] -> Scale(sc_Zbb2);
      Zbb3_mt[n]          -> Scale(sc_Zbb3);
      Zbb3_di_qt[n]       -> Scale(sc_Zbb3);
      Zbb3_met_et[n]      -> Scale(sc_Zbb3);
      Zbb3_met_over_qt[n] -> Scale(sc_Zbb3);
      Zbb3_met_et_ovQt[n] -> Scale(sc_Zbb3);
      Zcc0_mt[n]          -> Scale(sc_Zcc0);
      Zcc0_di_qt[n]       -> Scale(sc_Zcc0);
      Zcc0_met_et[n]      -> Scale(sc_Zcc0);
      Zcc0_met_over_qt[n] -> Scale(sc_Zcc0);
      Zcc0_met_et_ovQt[n] -> Scale(sc_Zcc0);
      Zcc1_mt[n]          -> Scale(sc_Zcc1);
      Zcc1_di_qt[n]       -> Scale(sc_Zcc1);
      Zcc1_met_et[n]      -> Scale(sc_Zcc1);
      Zcc1_met_over_qt[n] -> Scale(sc_Zcc1);
      Zcc1_met_et_ovQt[n] -> Scale(sc_Zcc1);
      Zcc2_mt[n]          -> Scale(sc_Zcc2);
      Zcc2_di_qt[n]       -> Scale(sc_Zcc2);
      Zcc2_met_et[n]      -> Scale(sc_Zcc2);
      Zcc2_met_over_qt[n] -> Scale(sc_Zcc2);
      Zcc2_met_et_ovQt[n] -> Scale(sc_Zcc2);
      Zcc3_mt[n]          -> Scale(sc_Zcc3);
      Zcc3_di_qt[n]       -> Scale(sc_Zcc3);
      Zcc3_met_et[n]      -> Scale(sc_Zcc3);
      Zcc3_met_over_qt[n] -> Scale(sc_Zcc3);
      Zcc3_met_et_ovQt[n] -> Scale(sc_Zcc3);
     
	}      
	  Zjets_mt[n]     = (TH1F*)DYmumu_mt[n];
	  Zjets_di_qt[n]  = (TH1F*)DYmumu_di_qt[n];
	  Zjets_met_et[n] = (TH1F*)DYmumu_met_et[n];
	  Zjets_met_over_qt[n] = (TH1F*)DYmumu_met_over_qt[n];
	  Zjets_met_et_ovQt[n] = (TH2F*)DYmumu_met_et_ovQt[n];
     
	
	  Zjets_mt[n]    -> Add(DYee_mt[n]);	 
	  Zjets_di_qt[n]  -> Add(DYee_di_qt[n]);
	  Zjets_met_et[n] -> Add(DYee_met_et[n]);
	  Zjets_met_over_qt[n] -> Add(DYee_met_over_qt[n]);
	  Zjets_met_et_ovQt[n] -> Add(DYee_met_et_ovQt[n]);
	  //Zjets_met_et[n] -> Add(DYee_met_et[n]);

	  Zjets_mt[n]    -> Add(DYtautau_mt[n]);	 
      Zjets_di_qt[n]       -> Add(DYtautau_di_qt[n]);
      Zjets_met_et[n]      -> Add(DYtautau_met_et[n]);
      Zjets_met_et_ovQt[n] -> Add(DYtautau_met_et_ovQt[n]);

      Zjets_mt[n]       -> SetFillColor(kRed+1);
      Zjets_mt[n]       -> SetLineColor(kGreen+2);
      Zjets_di_qt[n]       -> SetFillColor(kRed+1);
      Zjets_di_qt[n]       -> SetLineColor(kGreen+2);
      Zjets_met_et[n]      -> SetFillColor(kRed+1);
      Zjets_met_et[n]      -> SetLineColor(kGreen+2);
      Zjets_met_over_qt[n] -> SetFillColor(kRed+1);
      Zjets_met_over_qt[n] -> SetLineColor(kGreen+2);

      ZllG_mt[n]       -> SetLineColor(kYellow);
      ZllG_di_qt[n]       -> SetLineColor(kYellow);
      ZllG_met_et[n]      -> SetLineColor(kYellow);
      ZllG_met_over_qt[n] -> SetLineColor(kYellow);
      Wjets_mt[n]      -> SetFillColor(kCyan+1);
      Wjets_mt[n]     -> SetLineColor(kGreen-2);
      Wjets_di_qt[n]      -> SetFillColor(kCyan+1);
      Wjets_di_qt[n]     -> SetLineColor(kGreen-2);
      Wjets_met_et[n]     -> SetFillColor(kCyan+1);
      Wjets_met_et[n]     -> SetLineColor(kGreen-2);
      Wjets_met_over_qt[n] -> SetFillColor(kCyan);
      Wjets_met_over_qt[n] -> SetLineColor(kGreen-2);
      tt_mt[n]    -> SetFillColor(kGreen+2);
      tt_mt[n]    -> SetFillColor(kGreen+2);
      tt_met_over_qt[n]    -> SetFillColor(kGreen+2);
      tt_met_over_qt[n]    -> SetFillColor(kGreen+2);
      tt_met_over_qt[n]    -> SetLineColor(kOrange+1);
      tt_di_qt[n]         -> SetFillColor(kGreen+2);
      tt_di_qt[n]         -> SetLineColor(kOrange+1);
      tt_met_et[n]         -> SetFillColor(kGreen+2);
      tt_met_et[n]         -> SetLineColor(kOrange+1);
  
      t_all_mt[n]      -> SetFillColor(kOrange-3);
      t_all_mt[n]      -> SetLineColor(kBlue);
      t_all_di_qt[n]      -> SetFillColor(kOrange-3);
      t_all_di_qt[n]      -> SetLineColor(kBlue);
      t_all_met_et[n]      -> SetFillColor(kOrange-3);
      t_all_met_et[n]      -> SetLineColor(kBlue);
      t_all_met_over_qt[n] -> SetFillColor(kOrange-3);
      t_all_met_over_qt[n] -> SetLineColor(kBlue);
      ZZtoAny_mt[n]      -> SetFillColor(kMagenta);
      ZZtoAny_mt[n]      -> SetLineColor(kBlue);
      ZZtoAny_di_qt[n]      -> SetFillColor(kMagenta);
      ZZtoAny_di_qt[n]      -> SetLineColor(kBlue);
      ZZtoAny_met_et[n]      -> SetFillColor(kMagenta);
      ZZtoAny_met_et[n]      -> SetLineColor(kBlue);
      ZZtoAny_met_over_qt[n] -> SetFillColor(kMagenta);
      ZZtoAny_met_over_qt[n] -> SetLineColor(kBlue);

      WW_mt[n]        -> SetFillColor(kSpring);
      WW_mt[n]        -> SetLineColor(kBlue-2);
      WW_di_qt[n]        -> SetFillColor(kSpring);
      WW_di_qt[n]        -> SetLineColor(kBlue-2);
      WW_met_et[n]        -> SetFillColor(kSpring);
      WW_met_et[n]        -> SetLineColor(kBlue-2);
      WW_met_over_qt[n]   -> SetFillColor(kSpring);
      WW_met_over_qt[n]   -> SetLineColor(kBlue-2);
      WZ_mt[n]        -> SetFillColor(kBlue-2);
      WZ_mt[n]        -> SetLineColor(kRed+2);
      WZ_di_qt[n]        -> SetFillColor(kBlue-2);
      WZ_di_qt[n]        -> SetLineColor(kRed+2);
      WZ_met_et[n]        -> SetFillColor(kBlue-2);
      WZ_met_et[n]        -> SetLineColor(kRed+2);
      WZ_met_over_qt[n]   -> SetFillColor(kBlue-2);
      WZ_met_over_qt[n]   -> SetLineColor(kRed+2);


      hs_mt[n] = new THStack(Form("hs_mt_%i",n),"Stacked MT");
      hs_mt[n]->Add(t_all_mt[n]);
      hs_mt[n]->Add(tt_mt[n]);
      hs_mt[n]->Add(WZ_mt[n]);
      hs_mt[n]->Add(WW_mt[n]);
      hs_mt[n]->Add(ZZtoAny_mt[n]);
      //hs_mt[n]->Add(ZllG_mt[n]);
      hs_mt[n]->Add(Zjets_mt[n]);
      hs_mt[n]->Add(Wjets_mt[n]);
      //hs_mt[n]->Add(_mt[n]);

      sum_mt[n] = (TH1F*)t_all_mt[n]->Clone();
      sum_mt[n]->Add(tt_mt[n]);
      sum_mt[n]->Add(WZ_mt[n]);
      sum_mt[n]->Add(WW_mt[n]);
      sum_mt[n]->Add(ZZtoAny_mt[n]);
      sum_mt[n]->Add(Zjets_mt[n]);
      sum_mt[n]->Add(Wjets_mt[n]);
      //sum_mt[n]->Add(_mt[n]);

      hs_di_qt[n] = new THStack(Form("hs_di_qt_%i",n),"Stacked dilepton qT");
      hs_di_qt[n]->Add(t_all_di_qt[n]);
      hs_di_qt[n]->Add(tt_di_qt[n]);
      hs_di_qt[n]->Add(WZ_di_qt[n]);
      hs_di_qt[n]->Add(WW_di_qt[n]);
      hs_di_qt[n]->Add(ZZtoAny_di_qt[n]);
      //hs_di_qt[n]->Add(ZllG_di_qt[n]);
      hs_di_qt[n]->Add(Zjets_di_qt[n]);
      hs_di_qt[n]->Add(Wjets_di_qt[n]);
      //hs_di_qt[n]->Add(_di_qt[n]);

      sum_di_qt[n] = (TH1F*)t_all_di_qt[n]->Clone();
      sum_di_qt[n]->Add(tt_di_qt[n]);
      sum_di_qt[n]->Add(WZ_di_qt[n]);
      sum_di_qt[n]->Add(WW_di_qt[n]);
      sum_di_qt[n]->Add(ZZtoAny_di_qt[n]);
      sum_di_qt[n]->Add(Zjets_di_qt[n]);
      sum_di_qt[n]->Add(Wjets_di_qt[n]);
      //sum_di_qt[n]->Add(_di_qt[n]);

      hs_met_et_ovQt[n] = new THStack(Form("hs_met_et_ovQt_%i",n),"Stacked dilepton qT");
      hs_met_et_ovQt[n]->Add(t_all_met_et_ovQt[n]);
      hs_met_et_ovQt[n]->Add(tt_met_et_ovQt[n]);
      hs_met_et_ovQt[n]->Add(WZ_met_et_ovQt[n]);
      hs_met_et_ovQt[n]->Add(WW_met_et_ovQt[n]);
      hs_met_et_ovQt[n]->Add(ZZtoAny_met_et_ovQt[n]);
      //hs_met_et_ovQt[n]->Add(ZllG_met_et_ovQt[n]);
      hs_met_et_ovQt[n]->Add(Zjets_met_et_ovQt[n]);
      hs_met_et_ovQt[n]->Add(Wjets_met_et_ovQt[n]);
      //hs_met_et_ovQt[n]->Add(_met_et_ovQt[n]);

      sum_met_et_ovQt[n] = (TH2F*)t_all_met_et_ovQt[n]->Clone();
      sum_met_et_ovQt[n]->Add(tt_met_et_ovQt[n]);
      sum_met_et_ovQt[n]->Add(WZ_met_et_ovQt[n]);
      sum_met_et_ovQt[n]->Add(WW_met_et_ovQt[n]);
      sum_met_et_ovQt[n]->Add(ZZtoAny_met_et_ovQt[n]);
      sum_met_et_ovQt[n]->Add(Zjets_met_et_ovQt[n]);
      sum_met_et_ovQt[n]->Add(Wjets_met_et_ovQt[n]);
      //sum_met_et_ovQt[n]->Add(_met_et_ovQt[n]);


      hs_met_et[n] = new THStack(Form("hs_met_et_%i",n),"Stacked met");
      hs_met_et[n]->Add(t_all_met_et[n]);
      hs_met_et[n]->Add(tt_met_et[n]);
      hs_met_et[n]->Add(WZ_met_et[n]);
      hs_met_et[n]->Add(WW_met_et[n]);
      hs_met_et[n]->Add(ZZtoAny_met_et[n]);
      //hs_met_et[n]->Add(ZllG_met_et[n]);
      hs_met_et[n]->Add(Zjets_met_et[n]);
      hs_met_et[n]->Add(Wjets_met_et[n]);
      //hs_met_et[n]->Add(_met_et[n]);

      sum_met_et[n] = (TH1F*)t_all_met_et[n]->Clone();
      sum_met_et[n]->Add(tt_met_et[n]);
      sum_met_et[n]->Add(WZ_met_et[n]);
      sum_met_et[n]->Add(WW_met_et[n]);
      sum_met_et[n]->Add(ZZtoAny_met_et[n]);
      sum_met_et[n]->Add(Zjets_met_et[n]);
      sum_met_et[n]->Add(Wjets_met_et[n]);
      //sum_met_et[n]->Add(_met_et[n]);

      hs_met_over_qt[n] = new THStack(Form("hs_met_over_qt_%i",n),"Stacked met");
      hs_met_over_qt[n]->Add(t_all_met_over_qt[n]);
      hs_met_over_qt[n]->Add(tt_met_over_qt[n]);
      hs_met_over_qt[n]->Add(WZ_met_over_qt[n]);
      hs_met_over_qt[n]->Add(WW_met_over_qt[n]);
      hs_met_over_qt[n]->Add(ZZtoAny_met_over_qt[n]);
      //hs_met_over_qt[n]->Add(ZllG_met_over_qt[n]);
      hs_met_over_qt[n]->Add(Zjets_met_over_qt[n]);
      hs_met_over_qt[n]->Add(Wjets_met_over_qt[n]);
      //hs_met_over_qt[n]->Add(_met_over_qt[n]);

      sum_met_over_qt[n] = (TH1F*)t_all_met_over_qt[n]->Clone();
      sum_met_over_qt[n]->Add(tt_met_over_qt[n]);
      sum_met_over_qt[n]->Add(WZ_met_over_qt[n]);
      sum_met_over_qt[n]->Add(WW_met_over_qt[n]);
      sum_met_over_qt[n]->Add(ZZtoAny_met_over_qt[n]);
      sum_met_over_qt[n]->Add(Zjets_met_over_qt[n]);
      sum_met_over_qt[n]->Add(Wjets_met_over_qt[n]);
      //sum_met_over_qt[n]->Add(_met_over_qt[n]);

      if(histoPath!="v07" && histoPath!="v08")
	{
      ggH200_mt[n]          -> Scale(sc_ggH200);
      ggH200_di_qt[n]       -> Scale(sc_ggH200);
      ggH200_met_et_ovQt[n] -> Scale(sc_ggH200);
      ggH200_met_et[n]      -> Scale(sc_ggH200);
      ggH200_met_over_qt[n] -> Scale(sc_ggH200);

      ggH400_mt[n]          -> Scale(sc_ggH400);
      ggH400_di_qt[n]       -> Scale(sc_ggH400);
      ggH400_met_et_ovQt[n] -> Scale(sc_ggH400);
      ggH400_met_et[n]      -> Scale(sc_ggH400);
      ggH400_met_over_qt[n] -> Scale(sc_ggH400);

      ggH600_mt[n]          -> Scale(sc_ggH600);
      ggH600_di_qt[n]       -> Scale(sc_ggH600);
      ggH600_met_et_ovQt[n] -> Scale(sc_ggH600);
      ggH600_met_et[n]      -> Scale(sc_ggH600);
      ggH600_met_over_qt[n] -> Scale(sc_ggH600);

	}
      ggH200_mt[n]          -> SetLineColor(kYellow);
      ggH200_di_qt[n]       -> SetLineColor(kYellow);
      ggH200_met_et[n]      -> SetLineColor(kYellow);
      ggH200_met_over_qt[n] -> SetLineColor(kYellow);


      ggH400_mt[n]          -> SetLineColor(kCyan);
      ggH400_di_qt[n]       -> SetLineColor(kCyan);
      ggH400_met_et[n]      -> SetLineColor(kCyan);
      ggH400_met_over_qt[n] -> SetLineColor(kCyan);



      ggH600_mt[n]          -> SetLineColor(kCyan);
      ggH600_di_qt[n]       -> SetLineColor(kCyan);
      ggH600_met_et[n]      -> SetLineColor(kCyan);
      ggH600_met_over_qt[n] -> SetLineColor(kCyan);

      data_mt[n]            -> SetMarkerStyle(20);
      data_di_qt[n]         -> SetMarkerStyle(20);
      data_met_et[n]        -> SetMarkerStyle(20);
      data_met_over_qt[n]   -> SetMarkerStyle(20);
      data_met_phi[n]       -> SetMarkerStyle(20);
	  

    }



  hs_met_et[F0] -> Draw("hist");

  TCanvas *c2 = new TCanvas("c2","example",600,700);
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  c2 -> cd();
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  pad1 -> SetLogy();
  hs_met_et[F0] -> Draw("hist");
  hs_met_et[F0]   -> SetTitle("; pfMET; Events");
  hs_met_et[F0] -> SetMaximum(100000);
  hs_met_et[F0] -> SetMinimum(0.001);
  data_met_et[F0]   -> Draw("same e1pl");
  ggH400_met_et[F0] -> Draw("same hist");
  ZllG_met_et[F0]   -> Draw("same hist");

  leg01 = new TLegend(0.53,0.7,0.95,0.95);
  leg01->SetNColumns(2);
  leg01->SetTextSize(0.04);
  leg01->AddEntry(data_met_et[F0],    "Data","epl");
  leg01->AddEntry(Zjets_met_et[F0],   "Z + jets","f");
  leg01->AddEntry(ZZtoAny_met_et[F0], "ZZ","f");
  leg01->AddEntry(Wjets_met_et[F0],   "W + jets","f");
  leg01->AddEntry(WW_met_et[F0],      "WW","f");
  leg01->AddEntry(t_all_met_et[F0],   "t#rightarrow l#nub","f");
  leg01->AddEntry(WZ_met_et[F0],      "WZ","f");
  leg01->AddEntry(tt_met_et[F0],      "tt#rightarrow 2l2#nu2b","f");
  leg01->AddEntry(ZllG_met_et[F0],    "Z#gamma#rightarrow ll#gamma","f");
  leg01->AddEntry(ggH400_met_et[F0],  "10xH400","f");
  leg01->SetFillColor(kWhite);
  leg01 -> Draw();

  c2 -> cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetBottomMargin(0.25);
  pad2->SetTopMargin(0);
  pad2->Draw();
  pad2->cd();
  TH1F * h1 = (TH1F*) data_met_et[F0] -> Clone("temp1");
  TH1F * h2 = (TH1F*) sum_met_et[F0] -> Clone("temp2");
  //h1->SetStats(0);
  h1->Divide(h2);
  h1 -> SetTitle("; pfMET; Data/MC");
  if(metType=="met3") h1 -> SetTitle("; projMET; Data/MC");
  if(metType=="met4") h1 -> SetTitle("; PU corr projMET; Data/MC");
  h1 -> SetMaximum(5);
  h1 -> SetMinimum(0);
  h1->GetYaxis()->SetNdivisions(206);
  h1->GetYaxis()->SetTitleOffset(0.4);
  h1->SetTitleSize(0.1,"XYZ");
  h1 ->SetLabelSize(0.1,"XY");
  h1->Draw("ep");

  c2 -> SaveAs(imgpath+"ov01.png");

  c2 -> cd();
  pad1 -> cd();
  hs_met_over_qt[F0] -> Draw("hist");
  hs_met_over_qt[F0] -> SetTitle("; pfMET/q_{T}; Events");
  hs_met_over_qt[F0] -> SetMinimum(0.001);
  hs_met_over_qt[F0] -> SetMaximum(100000);
  data_met_over_qt[F0]   -> Draw("same e1pl");
  ggH400_met_over_qt[F0] -> Draw("same hist");
  ZllG_met_over_qt[F0]   -> Draw("same hist");
  leg01 -> Draw();

  pad2 -> cd();

  TH1F* h1 = (TH1F*) data_met_over_qt[F0]->Clone("");
  TH1F* h2 = (TH1F*) sum_met_over_qt[F0]-> Clone("");
  h1 -> Divide(h2);
  h1 -> SetTitle("; pfMET/q_{T}; Data/MC");
  if(metType=="met3") h1 -> SetTitle("; projMET/q_{T}; Data/MC");
  if(metType=="met4") h1 -> SetTitle("; PU corr projMET/q_{T}; Data/MC");
  h1 -> SetMaximum(5);
  h1 -> SetMinimum(0);
  h1->GetYaxis()->SetNdivisions(206);
  h1->GetYaxis()->SetTitleOffset(0.4);
  h1->SetTitleSize(0.1,"XYZ");
  h1 ->SetLabelSize(0.1,"XY");
  h1 -> Draw("ep");
  c2 -> SaveAs(imgpath+"ov02.png");


  pad1->cd();
  hs_met_et[F1] -> Draw("hist");
  hs_met_et[F1]   -> SetTitle("; pfMET; Events");
  hs_met_et[F1] -> SetMaximum(100000);
  hs_met_et[F1] -> SetMinimum(0.001);
  data_met_et[F1]   -> Draw("same e1pl");
  ggH400_met_et[F1] -> Draw("same hist");
  ZllG_met_et[F1]   -> Draw("same hist");
  leg01 -> Draw();

  pad2->cd();
  TH1F * h1 = (TH1F*) data_met_et[F1]->Clone("");
  TH1F * h2 = (TH1F*) sum_met_et[F1]-> Clone("");
  h1->Divide(h2);
  h1 -> SetTitle("; pfMET; Data/MC");
  if(metType=="met3") h1 -> SetTitle("; projMET; Data/MC");
  if(metType=="met4") h1 -> SetTitle("; PU corr projMET; Data/MC");
  h1 -> SetMaximum(5);
  h1 -> SetMinimum(0);
  h1->GetYaxis()->SetNdivisions(206);
  h1->GetYaxis()->SetTitleOffset(0.4);
  h1->SetTitleSize(0.1,"XYZ");
  h1 ->SetLabelSize(0.1,"XY");
  h1->Draw("ep");

  c2-> SaveAs(imgpath+"ov03.png");
 
  pad1 -> cd();
  hs_met_over_qt[F1] -> Draw("hist");
  hs_met_over_qt[F1] -> SetTitle("; pfMET/q_{T}; Events");
  hs_met_over_qt[F1] -> SetMinimum(0.001);
  hs_met_over_qt[F1] -> SetMaximum(100000);
  data_met_over_qt[F1]   -> Draw("same e1pl");
  ggH400_met_over_qt[F1] -> Draw("same hist");
  ZllG_met_over_qt[F1]   -> Draw("same hist");
  leg01 -> Draw();

  pad2 -> cd();
  h1 = (TH1F*) data_met_over_qt[F1]->Clone("temp1");
  h2 = (TH1F*) sum_met_over_qt[F1]-> Clone("temp2");
  h1 -> Divide(h2);
  h1 -> SetTitle("; pfMET/q_{T}; Data/MC");
  if(metType=="met3") h1 -> SetTitle("; projMET/q_{T}; Data/MC");
  if(metType=="met4") h1 -> SetTitle("; PU corr projMET/q_{T}; Data/MC");
  h1 -> SetMaximum(5);
  h1 -> SetMinimum(0);
  h1->GetYaxis()->SetNdivisions(206);
  h1->GetYaxis()->SetTitleOffset(0.4);
  h1->SetTitleSize(0.1,"XYZ");
  h1 ->SetLabelSize(0.1,"XY");
  h1 -> Draw("ep");

  c2 -> SaveAs(imgpath+"ov04.png");


  pad1 -> cd();
  hs_di_qt[F0] -> Draw("hist");
  hs_di_qt[F0] -> SetTitle("; ; Events");
  hs_di_qt[F0] -> SetMinimum(0.001);
  hs_di_qt[F0] -> SetMaximum(100000);
  data_di_qt[F0]   -> Draw("same e1pl");
  ggH400_di_qt[F0] -> Draw("same hist");
  ZllG_di_qt[F0]   -> Draw("same hist");
  leg01 -> Draw();

  pad2 -> cd();
  h1 = (TH1F*) data_di_qt[F0]->Clone("temp1");
  h2 = (TH1F*) sum_di_qt[F0]-> Clone("temp2");
  h1 -> Divide(h2);
  h1 -> SetTitle(";q_{T} (di-lepton p_{T}); Data/MC");
  h1 -> SetMaximum(10);
  h1 -> SetMinimum(0);
  h1->GetYaxis()->SetNdivisions(206);
  h1->GetYaxis()->SetTitleOffset(0.4);
  h1->SetTitleSize(0.1,"XYZ");
  h1 ->SetLabelSize(0.1,"XY");
  h1 -> Draw("ep");
  c2 -> SaveAs(imgpath+"ov05.png");


  pad1 -> cd();
  hs_di_qt[F1] -> Draw("hist");
  hs_di_qt[F1] -> SetTitle("; ; Events");
  hs_di_qt[F1] -> SetMinimum(0.001);
  hs_di_qt[F1] -> SetMaximum(100000);
  data_di_qt[F1]   -> Draw("same e1pl");
  ggH400_di_qt[F1] -> Draw("same hist");
  ZllG_di_qt[F1]   -> Draw("same hist");
  leg01 -> Draw();

  pad2 -> cd();
  h1 = (TH1F*) data_di_qt[F1]->Clone("temp1");
  h2 = (TH1F*) sum_di_qt[F1]-> Clone("temp2");
  h1 -> Divide(h2);
  h1 -> SetTitle(";q_{T} (di-lepton p_{T}); Data/MC");
  h1 -> SetMaximum(3);
  h1 -> SetMinimum(0);
  h1->GetYaxis()->SetNdivisions(206);
  h1->GetYaxis()->SetTitleOffset(0.4);
  h1->SetTitleSize(0.1,"XYZ");
  h1 ->SetLabelSize(0.1,"XY");
  h1 -> Draw("ep");
  c2 -> SaveAs(imgpath+"ov06.png");



  pad1 -> cd();
  hs_mt[F0] -> Draw("hist");
  hs_mt[F0] -> SetTitle("; ; Events");
  hs_mt[F0] -> SetMinimum(0.001);
  hs_mt[F0] -> SetMaximum(100000);
  data_mt[F0]   -> Draw("same e1pl");
  ggH400_mt[F0] -> Draw("same hist");
  ZllG_mt[F0]   -> Draw("same hist");
  leg01 -> Draw();

  pad2 -> cd();
  h1 = (TH1F*) data_mt[F0]->Clone("temp1");
  h2 = (TH1F*) sum_mt[F0]-> Clone("temp2");
  h1 -> Divide(h2);
  h1 -> SetTitle(";MT (using pfMET); Data/MC");
  h1 -> SetMaximum(5);
  h1 -> SetMinimum(0);
  h1->GetYaxis()->SetNdivisions(206);
  h1->GetYaxis()->SetTitleOffset(0.4);
  h1->SetTitleSize(0.1,"XYZ");
  h1 ->SetLabelSize(0.1,"XY");
  h1 -> Draw("ep");
  c2 -> SaveAs(imgpath+"ov07.png");


  pad1 -> cd();
  hs_mt[F1] -> Draw("hist");
  hs_mt[F1] -> SetTitle("; ; Events");
  hs_mt[F1] -> SetMinimum(0.001);
  hs_mt[F1] -> SetMaximum(100000);
  data_mt[F1]   -> Draw("same e1pl");
  ggH400_mt[F1] -> Draw("same hist");
  ZllG_mt[F1]   -> Draw("same hist");
  leg01 -> Draw();

  pad2 -> cd();
  h1 = (TH1F*) data_mt[F1]->Clone("temp1");
  h2 = (TH1F*) sum_mt[F1]-> Clone("temp2");
  h1 -> Divide(h2);
  h1 -> SetTitle(";MT (using pfMET); Data/MC");
  h1 -> SetMaximum(5);
  h1 -> SetMinimum(0);
  h1->GetYaxis()->SetNdivisions(206);
  h1->GetYaxis()->SetTitleOffset(0.4);
  h1->SetTitleSize(0.1,"XYZ");
  h1 ->SetLabelSize(0.1,"XY");
  h1 -> Draw("ep");
  c2 -> SaveAs(imgpath+"ov08.png");


  c1 ->cd();

 

  sum_met_et_ovQt[F0] -> Draw("colz");
  sum_met_et_ovQt[F0] -> SetTitle(";MET; MET/q_{T}");
  c1 -> SaveAs(imgpath+"ov09.png");

  ggH200_met_et_ovQt[F0] -> Draw("colz");
  ggH200_met_et_ovQt[F0] -> SetTitle(";MET; MET/q_{T}");
  c1 -> SaveAs(imgpath+"ov10.png");

  //ggH400_met_et_ovQt[F0] -> Draw("colz");
  //ggH400_met_et_ovQt[F0] -> SetTitle(";MET; MET/q_{T}");
  //c1 -> SaveAs(imgpath+"ov10.png");


  sum_met_et_ovQt[F1] -> Draw("colz");
  sum_met_et_ovQt[F1] -> SetTitle(";MET; MET/q_{T}");
  c1 -> SaveAs(imgpath+"ov11.png");

  ggH200_met_et_ovQt[F1] -> Draw("colz");
  ggH200_met_et_ovQt[F1] -> SetTitle(";MET; MET/q_{T}");
  c1 -> SaveAs(imgpath+"ov12.png");

  //ggH400_met_et_ovQt[F1] -> Draw("colz");
  //ggH400_met_et_ovQt[F1] -> SetTitle(";MET; MET/q_{T}");
  //c1 -> SaveAs(imgpath+"ov12.png");



  Double_t StoB = 0;
  Double_t signal = 0.1*ggH400_met_et[F0] ->Integral(8,41);
  Double_t backgr = sum_met_et[F0]->Integral(8,41);
  StoB = signal/backgr;
  cout<<"S/B = "<<StoB<<"  sig: "<<signal<<"  bg: "<<backgr<<endl;
  
  signal = 0.1*ggH400_met_et[F1] ->Integral(8,41);
  backgr = sum_met_et[F1]->Integral(8,41);
  StoB = signal/backgr;
  cout<<"S/B = "<<StoB<<"  sig: "<<signal<<"  bg: "<<backgr<<endl;
 
  signal = 0.1*ggH400_met_et_ovQt[F1]->Integral(8, 41, 0,41);
  backgr = sum_met_et_ovQt[F1]->Integral(8, 41, 0,41);
  StoB = signal/backgr;
  cout<<"S/B = "<<StoB<<"  sig: "<<signal<<"  bg: "<<backgr<<endl;

  signal = 0.1*ggH400_met_et_ovQt[F0]->Integral(8, 41, 0,41);
  backgr = sum_met_et_ovQt[F0]->Integral(8, 41, 0,41);
  StoB = signal/backgr;
  cout<<"S/B = "<<StoB<<"  sig: "<<signal<<"  bg: "<<backgr<<endl;



  //Float_t signal1 = 0.1*ggH200_met_et[F1] ->Integral(8,41);
  //Float_t signal1 = 0.1*ggH400_met_et[F1] ->Integral(8,41);
  Float_t signal1 = 0.1*ggH600_met_et[F1] ->Integral(8,41);
  //Float_t signal2 = 0.1*ggH200_met_et[F1] ->Integral(0,41);
  //Float_t signal2 = 0.1*ggH400_met_et[F1] ->Integral(0,41);
  Float_t signal2 = 0.1*ggH600_met_et[F1] ->Integral(0,41);
  Float_t backgr1 = sum_met_et[F1]->Integral(8,41);
  Float_t backgr2 = sum_met_et[F1]->Integral(0,41);
  Float_t eff_sig = signal1/signal2;
  Float_t eff_bkg = backgr1/backgr2;

  //StoB = signal/backgr;
  cout<<"sig eff:  "<<eff_sig<<"   bkg eff: "<<eff_bkg<<endl; 



  //Zjets_di_qt[F1] -> Draw("hist");
  //Zjets_mt[F1] -> Draw("hist");
  //Zjets_met_et[F1] -> Draw("hist");
  //c1 -> SaveAs(imgpath+"ov07.png");
  Zjets_di_qt[F1] -> Write();
  Zjets_mt[F1] -> Write();
  Zjets_met_et[F1] -> Write();

 


  fr -> Close();


  cout<<n<<" end dbg"<<endl;
  
}
