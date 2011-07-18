#include <iostream>
#include <string>

using namespace std;
//
//pair<string,string> getColors(string);
//Float_t getCS(string);



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
  TString sample1("MC_ggH400");
  //TString sample1("MC_DYmumu");
  //TString sample1("MC_ZZtoAny");

  TString sample2("2011A_v2_DoubleMu");
  TString sample2("none");
  if (aa!="none")  sample1 = aa;  
  if (bb!="none")  sample2 = bb;

  pathToFiles1 = Form("~/nobackup/higgsfiles4/%s",sample1.Data());
  pathToFiles2 = Form("~/nobackup/higgsfiles4/%s",sample2.Data());
 
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
  
  cout<<sample1<<"  and "<<sample2<<endl;

  gROOT->LoadMacro("./xSecAndColors.C");
  Float_t cs;
  //TString *fillColor, *lineColor; 
  Int_t fillColor, lineColor;
  cs = getCS(sample1.Data());
  fillColor = (TString*)getColors(sample1.Data(), 1);
  lineColor = (TString*)getColors(sample1.Data(), 2);
   //cout<<lineColor->Data()<<endl; 
  
  ofstream params;   //Parameters as an input to main analyzer!
  params.open("./params.txt");
  params<<selection<<endl;  //muon or electron selection
  params<<sample1<<endl;    // sample  (e.g. MC_DYmumu)
  params<<1000000<<endl;  //total number of events in the sample
  params<<cs<<endl;   //cross section for this ssample
  params<<fillColor<<endl;  //Color for fill histograms
  params<<lineColor<<endl;  //Colors for line
  //params<<fillColor->Data()<<endl;  //Color for fill histograms
  //params<<lineColor->Data()<<endl;  //Colors for line
  params.close();


//Ready to run the analyzer!
  TStopwatch timer;
  timer.Start();

  fChain->Process("analyzer_higgs.C+");

  cout << "\n\nDone!" << endl;
  cout << "CPU Time : " << timer.CpuTime() << endl;
  cout << "RealTime : " << timer.RealTime() << endl;
    
}

