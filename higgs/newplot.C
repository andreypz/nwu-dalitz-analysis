#include <iostream>

using namespace std;

void newplot(Int_t sel =1, TString hPath ="v00") {
  gROOT->LoadMacro("./makePlot.C");
  gROOT->LoadMacro("./merge.C");
  //gROOT->ProcessLine(".L ../data/tdrstyle.C");
  
  TString ssel("none");

  if (sel==1)  ssel = "muon";
  if (sel==2)  ssel = "electron";

  TString histoPath = Form("%s/%s", hPath.Data(), ssel.Data());
  cout<<histoPath.Data()<<endl;

  TStopwatch timer;	
  timer.Start();

  TFile* fda_2011A_DoubleMu_May10    = new TFile(Form("./%s/muon/hhhh_DoubleMu_May10.root",hPath.Data()));
  TFile* fda_2011A_DoubleMu_PromptV4 = new TFile(Form("./%s/muon/hhhh_DoubleMu_PromptV4.root",hPath.Data()));
  TFile* fda_2011A_DoubleMu_Aug05    = new TFile(Form("./%s/muon/hhhh_DoubleMu_Aug05.root",hPath.Data()));
  TFile* fda_2011A_DoubleMu_PromptV6 = new TFile(Form("./%s/muon/hhhh_DoubleMu_PromptV6.root",hPath.Data()));

  TFile* fda_2011A_DoubleEl_May10      = new TFile(Form("./%s/electron/hhhh_DoubleElectron_May10.root",hPath.Data()));
  TFile* fda_2011A_DoubleEl_PromptV4_1 = new TFile(Form("./%s/electron/hhhh_DoubleElectron_PromptV4_1.root",hPath.Data()));
  TFile* fda_2011A_DoubleEl_PromptV4_2 = new TFile(Form("./%s/electron/hhhh_DoubleElectron_PromptV4_2.root",hPath.Data()));
  TFile* fda_2011A_DoubleEl_PromptV4_3 = new TFile(Form("./%s/electron/hhhh_DoubleElectron_PromptV4_3.root",hPath.Data()));
  TFile* fda_2011A_DoubleEl_PromptV4_4 = new TFile(Form("./%s/electron/hhhh_DoubleElectron_PromptV4_4.root",hPath.Data()));
  TFile* fda_2011A_DoubleEl_Aug05      = new TFile(Form("./%s/electron/hhhh_DoubleElectron_Aug05.root",hPath.Data()));
  TFile* fda_2011A_DoubleEl_PromptV6   = new TFile(Form("./%s/electron/hhhh_DoubleElectron_PromptV6.root",hPath.Data()));


  TFile* fmc_DYmumu_1     = new TFile(Form("./%s/muon/hhhh_DYToMuMu_1.root",hPath.Data() ));
  TFile* fmc_DYmumu_2     = new TFile(Form("./%s/muon/hhhh_DYToMuMu_2.root",hPath.Data() ));
  TFile* fmc_DYmumu_3     = new TFile(Form("./%s/muon/hhhh_DYToMuMu_3.root",hPath.Data() ));
  TFile* fmc_DYmumu_4     = new TFile(Form("./%s/muon/hhhh_DYToMuMu_4.root",hPath.Data() ));
 
  TFile* fmc_DYee_1     = new TFile(Form("./%s/electron/hhhh_DYToEE_1.root",hPath.Data() ));
  TFile* fmc_DYee_2     = new TFile(Form("./%s/electron/hhhh_DYToEE_2.root",hPath.Data() ));
  TFile* fmc_DYee_3     = new TFile(Form("./%s/electron/hhhh_DYToEE_3.root",hPath.Data() ));
  TFile* fmc_DYee_4     = new TFile(Form("./%s/electron/hhhh_DYToEE_4.root",hPath.Data() ));

 // TFile* fmc_Wjets     = new TFile(Form("./%s/hhhh_Wjets.root",histoPath.Data() ));

  TFile* fmc_ZZtoAny   = new TFile(Form("./%s/hhhh_ZZ.root",histoPath.Data() ));
  TFile* fmc_tt_1      = new TFile(Form("./%s/hhhh_ttbar_1.root",histoPath.Data() ));
  TFile* fmc_tt_2      = new TFile(Form("./%s/hhhh_ttbar_2.root",histoPath.Data() ));
  TFile* fmc_tt_3      = new TFile(Form("./%s/hhhh_ttbar_3.root",histoPath.Data() ));
  TFile* fmc_tt_4      = new TFile(Form("./%s/hhhh_ttbar_4.root",histoPath.Data() ));
  TFile* fmc_WW        = new TFile(Form("./%s/hhhh_WW.root",histoPath.Data() ));
  TFile* fmc_WZ        = new TFile(Form("./%s/hhhh_WZ.root",histoPath.Data() ));
  TFile* fmc_tW        = new TFile(Form("./%s/hhhh_tW.root",histoPath.Data() ));


  TFile *m_Zjets = new TFile(Form("./m_Zjets_%i.root",sel), "RECREATE");
  TFile *m_Top   = new TFile(Form("./m_Top_%i.root",sel), "RECREATE");
  TFile *m_ttbar = new TFile(Form("./m_ttbar_%i.root",sel), "RECREATE");
  TFile *m_Data  = new TFile(Form("./m_Data_%i.root",sel), "RECREATE");
  
  list_Data = new TList();
  if(sel==1){
    list_Data -> Add(fda_2011A_DoubleMu_May10);
    list_Data -> Add(fda_2011A_DoubleMu_PromptV4);
    list_Data -> Add(fda_2011A_DoubleMu_Aug05);
    list_Data -> Add(fda_2011A_DoubleMu_PromptV6);
  }
  if(sel==2){
    list_Data -> Add(fda_2011A_DoubleEl_May10);
    list_Data -> Add(fda_2011A_DoubleEl_PromptV4_1);
    list_Data -> Add(fda_2011A_DoubleEl_PromptV4_2);
    list_Data -> Add(fda_2011A_DoubleEl_PromptV4_3);
    list_Data -> Add(fda_2011A_DoubleEl_PromptV4_4);
    list_Data -> Add(fda_2011A_DoubleEl_Aug05);
    list_Data -> Add(fda_2011A_DoubleEl_PromptV6);
  } 
  MergeRootfile(m_Data, list_Data, kBlack, 0);
  
  
  list_Zjets = new TList();
  
  if(sel==1){
    list_Zjets->Add(fmc_DYmumu_1);
    list_Zjets->Add(fmc_DYmumu_2);
    list_Zjets->Add(fmc_DYmumu_3);
    list_Zjets->Add(fmc_DYmumu_4);
  }
  if(sel==2){
    list_Zjets->Add(fmc_DYee_1);
    list_Zjets->Add(fmc_DYee_2);
    list_Zjets->Add(fmc_DYee_3);
    list_Zjets->Add(fmc_DYee_4);
  }
  //list_Zjets->Add(fmc_DYtautau);
  
  MergeRootfile(m_Zjets, list_Zjets, kGreen+2, kRed+1);


  list_ttbar = new TList();
  list_ttbar -> Add(fmc_tt_1);
  list_ttbar -> Add(fmc_tt_2);
  list_ttbar -> Add(fmc_tt_3);
  list_ttbar -> Add(fmc_tt_4);
  MergeRootfile(m_ttbar, list_ttbar, kOrange+1, kGreen+2);


  list_Top = new TList();
  //list_Top -> Add(fmc_tS);
  //list_Top -> Add(fmc_tT);
  list_Top -> Add(fmc_tW);
  MergeRootfile(m_Top, list_Top, kBlue, kOrange-3);


  m_Data  -> Close();
  m_Top -> Close();
  m_Zjets -> Close();
  m_ttbar -> Close();

  gROOT->ProcessLine(Form(".x makePlot.C(%i, \"%s\")", sel, histoPath.Data() ));

  cout << "\n\nDone!" << endl;
  cout << "CPU Time : " << timer.CpuTime() << endl;
  cout << "RealTime : " << timer.RealTime() << endl;
  
}

