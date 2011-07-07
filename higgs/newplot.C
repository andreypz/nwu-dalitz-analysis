#include <iostream>

using namespace std;

void newplot(Int_t sel =1, TString hPath ="v16") {
  gROOT->LoadMacro("./makePlot.C");
  gROOT->LoadMacro("./merge.C");
  //gROOT->ProcessLine(".L ../data/tdrstyle.C");
  
  TString histoPath = hPath.Data();
  cout<<histoPath.Data()<<endl;

  TStopwatch timer;	
  timer.Start();

  TFile* fmc_tS        = new TFile(Form("./%s/dir_%i_MC_tSchannel_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_tT        = new TFile(Form("./%s/dir_%i_MC_tTchannel_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_tW        = new TFile(Form("./%s/dir_%i_MC_tWchannel_/hhhh.root",histoPath.Data(), sel));

  TFile* fmc_DYmumu    = new TFile(Form("./%s/dir_%i_MC_DYmumu_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_DYee      = new TFile(Form("./%s/dir_%i_MC_DYee_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_DYtautau  = new TFile(Form("./%s/dir_%i_MC_DYtautau_/hhhh.root",histoPath.Data(), sel));

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

  TFile* fda_2011A_Prompt_v4_DoubleMu  = new TFile(Form("./%s/dir_1_2011A_Prompt_v4_DoubleMu_/hhhh.root",histoPath.Data()));
  TFile* fda_2011A_Prompt_v4_DoubleEl  = new TFile(Form("./%s/dir_2_2011A_Prompt_v4_DoubleElectron_/hhhh.root",histoPath.Data()));

  TFile *m_Zjets = new TFile(Form("./m_Zjets_%i.root",sel), "RECREATE");
  TFile *m_ZQQ   = new TFile(Form("./m_ZQQ_%i.root",sel), "RECREATE");
  TFile *m_Top   = new TFile(Form("./m_Top_%i.root",sel), "RECREATE");

  TFile *m_Data   = new TFile(Form("./m_Data_%i.root",sel), "RECREATE");
  
  
  list_Data = new TList();
  if(sel==1){
    list_Data ->Add(fda_2011A_DoubleMu);
    list_Data ->Add(fda_2011A_Prompt_v4_DoubleMu);
  }
  if(sel==2){
    list_Data ->Add(fda_2011A_DoubleEl);
    list_Data ->Add(fda_2011A_Prompt_v4_DoubleEl);
  } 
  MergeRootfile(m_Data, list_Data, kBlack, 0);
  

  list_Zjets = new TList();
  
  list_Zjets->Add(fmc_DYmumu);
  list_Zjets->Add(fmc_DYee);
  list_Zjets->Add(fmc_DYtautau);
  
  MergeRootfile(m_Zjets, list_Zjets, kGreen+2, kRed+1);

  //Zjets -> ls();

  list_Top = new TList();
  list_Top -> Add(fmc_tS);
  list_Top -> Add(fmc_tT);
  list_Top -> Add(fmc_tW);
  MergeRootfile(m_Top, list_Top, kBlue, kOrange-3);

  list_ZQQ = new TList();
  list_ZQQ -> Add(fmc_Zbb0);
  list_ZQQ -> Add(fmc_Zbb1);
  list_ZQQ -> Add(fmc_Zbb2);
  list_ZQQ -> Add(fmc_Zbb3);
  list_ZQQ -> Add(fmc_Zcc0);
  list_ZQQ -> Add(fmc_Zcc1);
  list_ZQQ -> Add(fmc_Zcc2);
  list_ZQQ -> Add(fmc_Zcc3);
  MergeRootfile(m_ZQQ, list_ZQQ, kYellow, kBlack);

  m_Zjets -> Close();
  m_ZQQ -> Close();
  m_Top -> Close();

  gROOT->ProcessLine(Form(".x makePlot.C(%i, \"%s\")", sel,histoPath.Data() ));

  cout << "\n\nDone!" << endl;
  cout << "CPU Time : " << timer.CpuTime() << endl;
  cout << "RealTime : " << timer.RealTime() << endl;
  
}

