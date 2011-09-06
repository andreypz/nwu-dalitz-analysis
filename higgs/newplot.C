#include <iostream>
#include "TChain.h"

using namespace std;

void newplot(Int_t sel =1, TString hPath ="v00", Int_t doMerge=0) {
  gROOT->LoadMacro("./makePlot.C");
  gROOT->LoadMacro("./merge.C");
  //gROOT->ProcessLine(".L ../data/tdrstyle.C");
  
  TStopwatch timer;	
  timer.Start();

  TString ssel("none");
  if (sel==1)  ssel = "muon";
  if (sel==2)  ssel = "electron";

  TString histoPath = Form("%s/%s", hPath.Data(), ssel.Data());
  cout<<histoPath.Data()<<endl;

  if(doMerge==0){
    gROOT->ProcessLine(Form(".x makePlot.C(%i, \"%s\")", sel, hPath.Data() ));
    
    cout << "\n\nDone!" << endl;
    cout << "CPU Time : " << timer.CpuTime() << endl;
    cout << "RealTime : " << timer.RealTime() << endl;
    return;
  }
  //Else: do the merging first:

  system (Form("rm ./%s/m_*.root ", hPath.Data()));

  system (Form("hadd ./%s/m_DataPh_1.root ./%s/muGamma/hhhh_Pho*", hPath.Data(), hPath.Data()));
  system (Form("hadd ./%s/m_DataPh_2.root ./%s/eGamma/hhhh_Pho*", hPath.Data(), hPath.Data()));

  system (Form("hadd ./%s/m_Data_1.root ./%s/muon/hhhh_Doub*", hPath.Data(), hPath.Data()));
  system (Form("hadd ./%s/m_Data_2.root ./%s/electron/hhhh_Doub*", hPath.Data(), hPath.Data()));

  system (Form("hadd ./%s/m_Top_1.root ./%s/muon/hhhh_tW*", hPath.Data(), hPath.Data()));
  system (Form("hadd ./%s/m_Top_2.root ./%s/electron/hhhh_tW*", hPath.Data(), hPath.Data()));

  system (Form("hadd ./%s/m_ttbar_1.root ./%s/muon/hhhh_tt*", hPath.Data(), hPath.Data()));
  system (Form("hadd ./%s/m_ttbar_2.root ./%s/electron/hhhh_tt*", hPath.Data(), hPath.Data()));

  system (Form("hadd ./%s/m_Zjets_1.root ./%s/muon/hhhh_DY*", hPath.Data(), hPath.Data()));
  system (Form("hadd ./%s/m_Zjets_2.root ./%s/electron/hhhh_DY*", hPath.Data(), hPath.Data()));


  Float_t photonLumi = 201.2. + 928.2;// + 407.5 +450.6;
  TFile *m_Zjets  = new TFile(Form("./%s/m_Zjets_%i.root", hPath.Data(), sel), "UPDATE");
  TFile *m_Top    = new TFile(Form("./%s/m_Top_%i.root", hPath.Data(), sel), "UPDATE");
  TFile *m_ttbar  = new TFile(Form("./%s/m_ttbar_%i.root", hPath.Data(), sel), "UPDATE");
  TFile *m_Data   = new TFile(Form("./%s/m_Data_%i.root", hPath.Data(), sel), "UPDATE");
  TFile *m_DataPh = new TFile(Form("./%s/m_DataPh_%i.root", hPath.Data(), sel), "UPDATE");
  
  RescaleToLumiAndColors(m_DataPh, photonLumi, 1000, kRed, -1);

  RescaleToLumiAndColors(m_Zjets, 1000,1000, kGreen+2, kRed+1);
  RescaleToLumiAndColors(m_ttbar, 1000,1000, kOrange+1, kGreen+2);
  RescaleToLumiAndColors(m_Top,   1000,1000,  kBlue, kOrange-3);

  /*
  TFile* fda_DoubleMu_May10    = new TFile(Form("./%s/muon/hhhh_DoubleMu_May10.root",hPath.Data()));
  TFile* fda_DoubleMu_PromptV4 = new TFile(Form("./%s/muon/hhhh_DoubleMu_PromptV4.root",hPath.Data()));
  TFile* fda_DoubleMu_Aug05    = new TFile(Form("./%s/muon/hhhh_DoubleMu_Aug05.root",hPath.Data()));
  TFile* fda_DoubleMu_PromptV6 = new TFile(Form("./%s/muon/hhhh_DoubleMu_PromptV6.root",hPath.Data()));

  TFile* fda_DoubleEl_May10      = new TFile(Form("./%s/electron/hhhh_DoubleElectron_May10.root",hPath.Data()));
  TFile* fda_DoubleEl_PromptV4_1 = new TFile(Form("./%s/electron/hhhh_DoubleElectron_PromptV4_1.root",hPath.Data()));
  TFile* fda_DoubleEl_PromptV4_2 = new TFile(Form("./%s/electron/hhhh_DoubleElectron_PromptV4_2.root",hPath.Data()));
  TFile* fda_DoubleEl_PromptV4_3 = new TFile(Form("./%s/electron/hhhh_DoubleElectron_PromptV4_3.root",hPath.Data()));
  TFile* fda_DoubleEl_PromptV4_4 = new TFile(Form("./%s/electron/hhhh_DoubleElectron_PromptV4_4.root",hPath.Data()));
  TFile* fda_DoubleEl_Aug05      = new TFile(Form("./%s/electron/hhhh_DoubleElectron_Aug05.root",hPath.Data()));
  TFile* fda_DoubleEl_PromptV6   = new TFile(Form("./%s/electron/hhhh_DoubleElectron_PromptV6.root",hPath.Data()));

  TFile *fmc_DYmumu[4], *fmc_DYee[4];
  for(Int_t i=0; i<4; i++){
     fmc_DYmumu[i]    = new TFile(Form("./%s/muon/hhhh_DYToMuMu_%i.root",hPath.Data(), i+1 ));
     fmc_DYee[i]      = new TFile(Form("./%s/electron/hhhh_DYToEE_%i.root",hPath.Data(), i+1 ));
  }
  TString gammasel("0");
  if (sel==1)  gammasel = "muGamma";
  if (sel==2)  gammasel = "eGamma";
  
  TString hGammaPath = Form("%s/%s", hPath.Data(), gammasel.Data());
  TFile*  fda_Photon_May10_1    = new TFile(Form("./%s/hhhh_Photon_May10_1.root",hGammaPath.Data()));
  TFile*  fda_Photon_May10_2    = new TFile(Form("./%s/hhhh_Photon_May10_2.root",hGammaPath.Data()));
  TFile*  fda_Photon_May10_3    = new TFile(Form("./%s/hhhh_Photon_May10_3.root",hGammaPath.Data()));
  TFile*  fda_Photon_Aug05_1    = new TFile(Form("./%s/hhhh_Photon_Aug05_1.root",hGammaPath.Data()));
  TFile*  fda_Photon_Aug05_2    = new TFile(Form("./%s/hhhh_Photon_Aug05_2.root",hGammaPath.Data()));
  TFile*  fda_Photon_Aug05_3    = new TFile(Form("./%s/hhhh_Photon_Aug05_3.root",hGammaPath.Data()));
  TFile*  fda_Photon_PromptV4_1 = new TFile(Form("./%s/hhhh_Photon_PromptV4_1.root",hGammaPath.Data()));
  TFile*  fda_Photon_PromptV4_2 = new TFile(Form("./%s/hhhh_Photon_PromptV4_2.root",hGammaPath.Data()));
  TFile*  fda_Photon_PromptV4_3 = new TFile(Form("./%s/hhhh_Photon_PromptV4_3.root",hGammaPath.Data()));
  TFile*  fda_Photon_PromptV4_4 = new TFile(Form("./%s/hhhh_Photon_PromptV4_4.root",hGammaPath.Data()));
  TFile*  fda_Photon_PromptV4_5 = new TFile(Form("./%s/hhhh_Photon_PromptV4_5.root",hGammaPath.Data()));
  TFile*  fda_Photon_PromptV6_1 = new TFile(Form("./%s/hhhh_Photon_PromptV6_1.root",hGammaPath.Data()));
  TFile*  fda_Photon_PromptV6_2 = new TFile(Form("./%s/hhhh_Photon_PromptV6_2.root",hGammaPath.Data()));
  TFile*  fda_Photon_PromptV6_3 = new TFile(Form("./%s/hhhh_Photon_PromptV6_3.root",hGammaPath.Data()));
  TFile*  fda_Photon_PromptV6_4 = new TFile(Form("./%s/hhhh_Photon_PromptV6_4.root",hGammaPath.Data()));
  TFile*  fda_Photon_PromptV6_5 = new TFile(Form("./%s/hhhh_Photon_PromptV6_5.root",hGammaPath.Data()));

  TFile*  fda_Gamma_PromptV6_1 = new TFile(Form("./%s/gamma/hhhh_Photon_PromptV6_1.root",hPath.Data()));
  TFile*  fda_Gamma_PromptV6_2 = new TFile(Form("./%s/gamma/hhhh_Photon_PromptV6_2.root",hPath.Data()));
  TFile*  fda_Gamma_PromptV6_3 = new TFile(Form("./%s/gamma/hhhh_Photon_PromptV6_3.root",hPath.Data()));
  TFile*  fda_Gamma_PromptV6_4 = new TFile(Form("./%s/gamma/hhhh_Photon_PromptV6_4.root",hPath.Data()));
  TFile*  fda_Gamma_PromptV6_5 = new TFile(Form("./%s/gamma/hhhh_Photon_PromptV6_5.root",hPath.Data()));

 // TFile* fmc_Wjets     = new TFile(Form("./%s/hhhh_Wjets.root",histoPath.Data() ));

  TFile* fmc_tt_1      = new TFile(Form("./%s/hhhh_ttbar_1.root",histoPath.Data() ));
  TFile* fmc_tt_2      = new TFile(Form("./%s/hhhh_ttbar_2.root",histoPath.Data() ));
  TFile* fmc_tt_3      = new TFile(Form("./%s/hhhh_ttbar_3.root",histoPath.Data() ));
  TFile* fmc_tt_4      = new TFile(Form("./%s/hhhh_ttbar_4.root",histoPath.Data() ));
  TFile* fmc_ZZtoAny   = new TFile(Form("./%s/hhhh_ZZ.root",histoPath.Data() ));
  TFile* fmc_WW        = new TFile(Form("./%s/hhhh_WW.root",histoPath.Data() ));
  TFile* fmc_WZ        = new TFile(Form("./%s/hhhh_WZ.root",histoPath.Data() ));
  TFile* fmc_tW        = new TFile(Form("./%s/hhhh_tW.root",histoPath.Data() ));


  TFile *m_Zjets  = new TFile(Form("./%s/m_Zjets_%i.root", histoPath.Data(), sel), "RECREATE");
  TFile *m_Top    = new TFile(Form("./%s/m_Top_%i.root", histoPath.Data(), sel), "RECREATE");
  TFile *m_ttbar  = new TFile(Form("./%s/m_ttbar_%i.root", histoPath.Data(), sel), "RECREATE");
  TFile *m_Data   = new TFile(Form("./%s/m_Data_%i.root", histoPath.Data(), sel), "RECREATE");
  TFile *m_DataPh = new TFile(Form("./%s/m_DataPh_%i.root", histoPath.Data(), sel), "RECREATE");
  
  list_Data = new TList();
  if(sel==1){
    list_Data -> Add(fda_DoubleMu_May10);
    list_Data -> Add(fda_DoubleMu_PromptV4);
    //   list_Data -> Add(fda_DoubleMu_Aug05);
    //list_Data -> Add(fda_DoubleMu_PromptV6);
  }
  if(sel==2){
    list_Data -> Add(fda_DoubleEl_May10);
    list_Data -> Add(fda_DoubleEl_PromptV4_1);
    list_Data -> Add(fda_DoubleEl_PromptV4_2);
    list_Data -> Add(fda_DoubleEl_PromptV4_3);
    list_Data -> Add(fda_DoubleEl_PromptV4_4);
    list_Data -> Add(fda_DoubleEl_Aug05);
    list_Data -> Add(fda_DoubleEl_PromptV6);
  } 
  //MergeRootfile(m_Data, list_Data, kBlack, 0);

  list_DataPh = new TList();
  list_DataPh -> Add(fda_Photon_May10_1);
  list_DataPh -> Add(fda_Photon_May10_2);
  list_DataPh -> Add(fda_Photon_May10_3);
  list_DataPh -> Add(fda_Photon_PromptV4_1);
  list_DataPh -> Add(fda_Photon_PromptV4_2);
  list_DataPh -> Add(fda_Photon_PromptV4_3);
  list_DataPh -> Add(fda_Photon_PromptV4_4);
  list_DataPh -> Add(fda_Photon_PromptV4_5);
  list_DataPh -> Add(fda_Photon_Aug05_1);
  list_DataPh -> Add(fda_Photon_Aug05_2);
  list_DataPh -> Add(fda_Photon_Aug05_3);
  list_DataPh -> Add(fda_Photon_PromptV6_1);
  list_DataPh -> Add(fda_Photon_PromptV6_2);
  list_DataPh -> Add(fda_Photon_PromptV6_3);
  list_DataPh -> Add(fda_Photon_PromptV6_4);
  list_DataPh -> Add(fda_Photon_PromptV6_5);
  //    MergeRootfile(m_DataPh, list_DataPh, kRed+1, 0);

  //RescaleToLumi(m_DataPh, photonLumi, 1000);
  
  list_Zjets = new TList();
  if(sel==1){
    for(Int_t i=0; i<4; i++)
      list_Zjets->Add(fmc_DYmumu[i]);
  }
  if(sel==2){
    for(Int_t i=0; i<4; i++)
      list_Zjets->Add(fmc_DYee[i]);
  }
  //list_Zjets->Add(fmc_DYtautau);
  
  //  MergeRootfile(m_Zjets, list_Zjets, kGreen+2, kRed+1);

  list_ttbar = new TList();
  list_ttbar -> Add(fmc_tt_1);
  list_ttbar -> Add(fmc_tt_2);
  list_ttbar -> Add(fmc_tt_3);
  list_ttbar -> Add(fmc_tt_4);
  // MergeRootfile(m_ttbar, list_ttbar, kOrange+1, kGreen+2);


  list_Top = new TList();
  //list_Top -> Add(fmc_tS);
  //list_Top -> Add(fmc_tT);
  list_Top -> Add(fmc_tW);
  // MergeRootfile(m_Top, list_Top, kBlue, kOrange-3);
*/


  m_DataPh -> Close();
  m_Data   -> Close();
  m_Top    -> Close();
  m_Zjets  -> Close();
  m_ttbar  -> Close();

  gROOT->ProcessLine(Form(".x makePlot.C(%i, \"%s\")", sel, hPath.Data() ));

  cout << "\n\nDone!" << endl;
  cout << "CPU Time : " << timer.CpuTime() << endl;
  cout << "RealTime : " << timer.RealTime() << endl;
  
}

