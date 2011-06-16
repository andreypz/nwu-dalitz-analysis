#include <iostream>

using namespace std;

void run_analysis(UInt_t sel = 1, string aa = "none", string bb="none") {
  gROOT->LoadMacro("../src/TCPrimaryVtx.cc+");
  gROOT->LoadMacro("../src/TCMuon.cc+");
  gROOT->LoadMacro("../src/TCElectron.cc+");
  gROOT->LoadMacro("../src/TCJet.cc+");
  gROOT->LoadMacro("../src/TCMET.cc+");
  
  TChain* fChain = new TChain("higgs/rTree");
  
  TString selection("none");
  if (sel==1) selection = "muon";
  if (sel==2) selection = "electron";

  TString pathToFiles1, pathToFiles2;
  //TString sample1("2010A");
  //TString sample1("2010B");
  //TString sample1("2011A_May10_DoubleElectron");
  //TString sample1("2011A_v1_DoubleMu");
  TString sample1("MC_DYmumu");
  //TString sample1("DYee");
  //TString sample1("ggH400");

  TString sample2("2011A_v2_DoubleMu");
  TString sample2("none");
  if (aa!="none")  sample1 = aa;  
  if (bb!="none")  sample2 = bb;

  pathToFiles1 = Form("~/nobackup/higgsfiles4/%s",sample1.Data());
  pathToFiles2 = Form("~/nobackup/higgsfiles4/%s",sample2.Data());
  //pathToFiles = Form("~/nobackup/higgsfiles3/%s",sample.Data());
 
  system (Form("ls %s/out_higgs_* > files.txt", pathToFiles1.Data()));
  if(sample2!="none")  system (Form("ls %s/out_higgs_* >> files.txt", pathToFiles2.Data()));

  ifstream sourceFiles("files.txt");

  string line;	
  while(getline(sourceFiles, line))
    {
      cout<<"Added file: "<<line<<endl;
      fChain->Add(line.c_str());
    }
  
  sourceFiles.close();
  
  //fChain->Add("./out_higgs.root");  

  TStopwatch timer;
  timer.Start();

  cout<<sample1<<"  and "<<sample2<<endl;

  ofstream params;
  params.open("./params.txt");
  params<<selection<<endl;
  params.close();

  fChain->Process("analyzer_higgs.C+");

  cout << "\n\nDone!" << endl;
  cout << "CPU Time : " << timer.CpuTime() << endl;
  cout << "RealTime : " << timer.RealTime() << endl;
  
  /*

  //TFile* fMC   = new TFile("./histogramsMC_D6T.root");

  gROOT -> SetStyle("Plain");
  gStyle -> SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetPaintTextFormat("4.2f"); 

  TString imgpath("~/afs/public_html/dndeta/");
  //TString imgpath("~/afs/public_html/test/");
  

  */
  
}

