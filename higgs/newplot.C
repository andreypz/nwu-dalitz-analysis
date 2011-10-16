#include <iostream>
#include "TChain.h"

using namespace std;

void newplot(Int_t sel =1, TString hPath ="v00", Int_t doMerge=0) {
  gROOT->LoadMacro("./makePlot.C");
  gROOT->LoadMacro("./utils.C");
  
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

  system (Form("rm ./%s/m_*%i.root ", hPath.Data(), sel));

//  system (Form("hadd ./%s/m_DataPh_%i.root ./%s/muGamma/hhhh_Pho*", hPath.Data(), sel, hPath.Data()));
  
  system (Form("hadd ./%s/m_Data_%i.root  ./%s/hhhh_Doub*.root ",     hPath.Data(), sel, histoPath.Data()));
  system (Form("hadd ./%s/m_Top_%i.root   ./%s/hhhh_t*W.root ", hPath.Data(), sel,  histoPath.Data()));
  system (Form("hadd ./%s/m_ttbar_%i.root ./%s/hhhh_tt*",       hPath.Data(), sel, histoPath.Data()));
  if (sel==1) system (Form("hadd ./%s/m_Zjets_%i.root ./%s/hhhh_DYToMuMu*",       hPath.Data(), sel, histoPath.Data()));
  if (sel==2) system (Form("hadd ./%s/m_Zjets_%i.root ./%s/hhhh_DYToEE*",       hPath.Data(), sel, histoPath.Data()));
  //system (Form("hadd ./%s/m_Zjets_%i.root ./%s/hhhh_DY*",       hPath.Data(), sel, histoPath.Data()));


  //Float_t photonLumi = 201.2. + 928.2;// + 407.5 +450.6;
  TFile *m_Zjets  = new TFile(Form("./%s/m_Zjets_%i.root",  hPath.Data(), sel), "UPDATE");
  TFile *m_Top    = new TFile(Form("./%s/m_Top_%i.root",    hPath.Data(), sel), "UPDATE");
  TFile *m_ttbar  = new TFile(Form("./%s/m_ttbar_%i.root",  hPath.Data(), sel), "UPDATE");
  TFile *m_Data   = new TFile(Form("./%s/m_Data_%i.root",   hPath.Data(), sel), "UPDATE");
//  TFile *m_DataPh = new TFile(Form("./%s/m_DataPh_%i.root", hPath.Data(), sel), "UPDATE");
  
  //RescaleToLumiAndColors(m_DataPh, photonLumi, 1000, kRed, -1);
  RescaleToLumiAndColors(m_ttbar, 1000,1000, kOrange+1, kGreen+2);
  RescaleToLumiAndColors(m_Top,   1000,1000, kBlue, kOrange-3);
  RescaleToLumiAndColors(m_Zjets, 1000,1000, kGreen+2, kRed+1);


  //m_DataPh -> Close();
  m_Data   -> Close();
  m_Top    -> Close();
  m_Zjets  -> Close();
  m_ttbar  -> Close();

  gROOT->ProcessLine(Form(".x makePlot.C(%i, \"%s\")", sel, hPath.Data() ));

  cout << "\n\nDone!" << endl;
  cout << "CPU Time : " << timer.CpuTime() << endl;
  cout << "RealTime : " << timer.RealTime() << endl;
  
}

