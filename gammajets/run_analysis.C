#include <iostream>

using namespace std;

void run_analysis(UInt_t sel = 3, string aa = "none", string bb="none") {
  gROOT->LoadMacro("../src/TCPrimaryVtx.cc+");
  gROOT->LoadMacro("../src/TCMuon.cc+");
  gROOT->LoadMacro("../src/TCElectron.cc+");
  gROOT->LoadMacro("../src/TCJet.cc+");
  gROOT->LoadMacro("../src/TCMET.cc+");
  gROOT->LoadMacro("../src/TCPhoton.cc+");
  gROOT->LoadMacro("../src/TCTrigger.cc+");
  
  TChain* fChain = new TChain("higgs/rTree");


  TString selection("none");
  if (sel==1) selection = "muon";
  if (sel==2) selection = "electron";
  if (sel==3) selection = "photon";


// Radek's ntuples
//fChain->Add("/uscms_data/d2/radek/HZZAna/Photon2011A-May10ReReco-v2-r1/sample/*root");


// ntuples remade by Nate
///fChain->Add("/uscms_data/d2/aa/ntuples/gamma/*root");


// new ntuples production from Radek
// Prompt 
fChain->Add("/uscms_data/d2/radek/HZZAna/Photon2011A-PromptReco-v1-r2/samples/*root");

// re-reco
fChain->Add("/uscms_data/d2/radek/HZZAna/Photon2011A-May10ReReco-v3-r2/samples/*root");


  TStopwatch timer;
  timer.Start();

  ofstream params;
  params.open("./params.txt");
  params<<selection<<endl;
  params.close();

  fChain->Process("analyzer_higgs.C+");

  cout << "\n\nDone!" << endl;
  cout << "CPU Time : " << timer.CpuTime() << endl;
  cout << "RealTime : " << timer.RealTime() << endl;
  
}

