#include <iostream>

using namespace std;

void newplot(Int_t sel =1, TString hPath ="v30") {
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

  TFile* fda_2011A_DoubleMu_May10  = new TFile(Form("./%s/hhhh_DoubleMu_May10.root",histoPath.Data()));
  TFile* fda_2011A_DoubleEl_May10  = new TFile(Form("./%s/hhhh_DoubleMu_May10.root",histoPath.Data()));
  TFile* fda_2011A_Prompt_v4_DoubleMu  = new TFile(Form("./%s/hhhh_DoubleMu_PromptV4.root",histoPath.Data()));
  TFile* fda_2011A_Prompt_v4_DoubleEl  = new TFile(Form("./%s/hhhh_DoubleMu_PromptV4.root",histoPath.Data()));


  TFile* fmc_DYmumu_1     = new TFile(Form("./%s/hhhh_DYToMuMu_1.root",histoPath.Data() ));
  TFile* fmc_DYmumu_2     = new TFile(Form("./%s/hhhh_DYToMuMu_2.root",histoPath.Data() ));
  TFile* fmc_DYmumu_3     = new TFile(Form("./%s/hhhh_DYToMuMu_3.root",histoPath.Data() ));
  TFile* fmc_DYmumu_4     = new TFile(Form("./%s/hhhh_DYToMuMu_4.root",histoPath.Data() ));
  // TFile* fmc_Wjets     = new TFile(Form("./%s/hhhh_Wjets.root",histoPath.Data() ));

  TFile* fmc_ZZtoAny   = new TFile(Form("./%s/hhhh_ZZ.root",histoPath.Data() ));
  TFile* fmc_tt_1      = new TFile(Form("./%s/hhhh_ttbar_1.root",histoPath.Data() ));
  TFile* fmc_tt_2      = new TFile(Form("./%s/hhhh_ttbar_2.root",histoPath.Data() ));
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
    list_Data -> Add(fda_2011A_Prompt_v4_DoubleMu);
  }
  if(sel==2){
    list_Data -> Add(fda_2011A_DoubleEl_May10);
    list_Data -> Add(fda_2011A_Prompt_v4_DoubleEl);
  } 
  MergeRootfile(m_Data, list_Data, kBlack, 0);
  
  
  list_Zjets = new TList();
  
  list_Zjets->Add(fmc_DYmumu_1);
  list_Zjets->Add(fmc_DYmumu_2);
  list_Zjets->Add(fmc_DYmumu_3);
  list_Zjets->Add(fmc_DYmumu_4);
  //list_Zjets->Add(fmc_DYee);
  //list_Zjets->Add(fmc_DYtautau);
  
  MergeRootfile(m_Zjets, list_Zjets, kGreen+2, kRed+1);

  //Zjets -> ls();

  list_ttbar = new TList();
  list_ttbar -> Add(fmc_tt_1);
  list_ttbar -> Add(fmc_tt_2);
  MergeRootfile(m_ttbar, list_ttbar, kOrange+1, kGreen+2);


  list_Top = new TList();
  //list_Top -> Add(fmc_tS);
  //list_Top -> Add(fmc_tT);
  list_Top -> Add(fmc_tW);
  MergeRootfile(m_Top, list_Top, kBlue, kOrange-3);

  m_Top -> Close();

  m_Data  -> Close();
  m_Zjets -> Close();
  m_ttbar -> Close();

  gROOT->ProcessLine(Form(".x makePlot.C(%i, \"%s\")", sel, histoPath.Data() ));

  cout << "\n\nDone!" << endl;
  cout << "CPU Time : " << timer.CpuTime() << endl;
  cout << "RealTime : " << timer.RealTime() << endl;
  
}

