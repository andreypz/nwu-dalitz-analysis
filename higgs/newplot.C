#include <iostream>
#include "TChain.h"

using namespace std;

void newplot(Int_t sel = 1, TString hPath ="v00", Int_t doMerge=0) {
  gROOT->LoadMacro("./makePlot.C");
  gROOT->LoadMacro("./utils.C");
  
  TStopwatch timer;	
  timer.Start();

  TString ssel("none"), gsel("none");
  if (sel==1)  {ssel = "muon";     gsel ="muGamma";}
  if (sel==2)  {ssel = "electron"; gsel ="eGamma";}

  TString histoPath = Form("%s/%s", hPath.Data(), ssel.Data());
  TString gammaPath = Form("%s/%s", hPath.Data(), gsel.Data());
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

  system (Form("hadd ./%s/m_DataPh_%i.root ./%s/hhhh_*Photon*", hPath.Data(), sel, gammaPath.Data()));
  
  system (Form("hadd ./%s/m_Data_%i.root  ./%s/hhhh_Doub*.root ", hPath.Data(), sel, histoPath.Data()));
  system (Form("hadd ./%s/m_Top_%i.root   ./%s/hhhh_t*W.root ",   hPath.Data(), sel,  histoPath.Data()));
  system (Form("hadd ./%s/m_ttbar_%i.root ./%s/hhhh_tt*",         hPath.Data(), sel, histoPath.Data()));

  //system (Form("hadd ./%s/m_libZjets_%i.root ./%s/lib5EtaBoost/hhhh_DYjets*",       hPath.Data(), sel, histoPath.Data()));
 // system (Form("hadd ./%s/m_libZjets_%i.root ./%s/lib3Zjets/hhhh_DYjets*",       hPath.Data(), sel, histoPath.Data()));
  system (Form("hadd ./%s/m_Zjets_%i.root ./%s/hhhh_DYjets*",       hPath.Data(), sel, histoPath.Data()));

  Float_t photonLumi = 1;
  if(sel==1) photonLumi = 215.1 + 927.6 + 370.9 + 663.0 + 2511; //double mu     
  if(sel==2) photonLumi = 215.1 + 789.2 + 313.2 + 662.2 + 2511; //double ele    
 

  TFile *m_libZjets = new TFile(Form("./%s/m_libZjets_%i.root",  hPath.Data(), sel), "UPDATE");
  TFile *m_Zjets    = new TFile(Form("./%s/m_Zjets_%i.root",  hPath.Data(), sel), "UPDATE");
  TFile *m_Top      = new TFile(Form("./%s/m_Top_%i.root",    hPath.Data(), sel), "UPDATE");
  TFile *m_ttbar    = new TFile(Form("./%s/m_ttbar_%i.root",  hPath.Data(), sel), "UPDATE");
  TFile *m_Data     = new TFile(Form("./%s/m_Data_%i.root",   hPath.Data(), sel), "UPDATE");
  TFile *m_DataPh   = new TFile(Form("./%s/m_DataPh_%i.root", hPath.Data(), sel), "UPDATE");
  

  m_Data     -> Close();

  RescaleToLumiAndColors(m_DataPh,1, photonLumi, 1000, kRed+1, kRed-4, 3325);
  m_DataPh   -> Close();
  RescaleToLumiAndColors(m_ttbar,1, 1000,1000, kMagenta+1, kBlue-3, 1001);
  m_ttbar    -> Close();
  RescaleToLumiAndColors(m_Top,1,   1000,1000, kOrange+9, kOrange+6,1001);
  m_Top      -> Close();
  RescaleToLumiAndColors(m_Zjets,1, 1000,1000, kRed+2, kRed+1,3004);
  m_Zjets    -> Close();
 // RescaleToLumiAndColors(m_libZjets,1, 1000,1000, kGreen+2, kRed+1,3357);
  m_libZjets -> Close();




  gROOT->ProcessLine(Form(".x makePlot.C(%i, \"%s\")", sel, hPath.Data() ));

  cout << "\n\nDone!" << endl;
  cout << "CPU Time : " << timer.CpuTime() << endl;
  cout << "RealTime : " << timer.RealTime() << endl;
  
}

