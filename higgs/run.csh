#!/bin/csh

echo "Start runninng"

set dir=`echo $1 | cut -d _ -f 1 `

cp template_higgsAnalyzer.C higgsAnalyzer.C

sed -i "s/SUFFIX/$1/g" higgsAnalyzer.C
sed -i "s/TRIGGER/$2/g" higgsAnalyzer.C
sed -i "s/SELECTION/$4/g" higgsAnalyzer.C

cat > run.C << +EOF
    

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;

void run() {

gROOT->LoadMacro("../src/TCJet.cc+");
gROOT->LoadMacro("../src/TCMET.cc+");
gROOT->LoadMacro("../src/TCElectron.cc+");
gROOT->LoadMacro("../src/TCMuon.cc+");
gROOT->LoadMacro("../src/TCTau.cc+");
gROOT->LoadMacro("../src/TCPhoton.cc+");
gROOT->LoadMacro("../src/TCGenJet.cc+");
gROOT->LoadMacro("../src/TCGenParticle.cc+");
gROOT->LoadMacro("../src/TCPrimaryVtx.cc+");
gROOT->LoadMacro("../src/TCTrigger.cc+");
gROOT->LoadMacro("../src/TCTriggerObject.cc+");
gROOT->LoadMacro("../src/WeightUtils.cc+");
     
TChain* fChain = new TChain("ntupleProducer/eventTree");

ifstream sourceFiles("./sourceFiles/$3.txt");
string line;
int  count = 0;
  cout<<"Adding files from $3 to chain..."<<endl;

	  while (sourceFiles >> line) {
	    fChain->Add(line.c_str());      
  ++count;
}
cout<<count<<" files added!"<<endl;
sourceFiles.close();


  TString sample1("$1");
  
cout<<"sample:  "<<sample1<<"  selection:  $4"<<endl;

gROOT->LoadMacro("./xSecAndColors.C");
Float_t cs, Ne;
Int_t fillColor, lineColor;
fillColor = (Int_t)getXsecOrColors(sample1.Data(), 1);
lineColor = (Int_t)getXsecOrColors(sample1.Data(), 2);
cs = getXsecOrColors(sample1.Data(), 3);
Ne = getXsecOrColors(sample1.Data(), 4);

ofstream params;   //Parameters as an input to main analyzer!                                                                                              
params.open("./params.txt");
params<<"$4"<<endl;  //muon or electron selection                                                                                                     
params<<sample1<<endl;    // sample  (e.g. MC_DYmumu)                                                                                                      
params<<Ne<<endl;  //total number of events in the sample
params<<cs<<endl;   //cross section for this ssample                                                                                                       
params<<fillColor<<endl;  //Color for fill histograms                                                                                                      
params<<lineColor<<endl;  //Colors for line
params.close();

  TStopwatch timer;
  timer.Start();

  fChain->Process("higgsAnalyzer.C+");

  cout << "\n\nDone!" << endl;
  cout << "CPU Time : " << timer.CpuTime() << endl;
  cout << "RealTime : " << timer.RealTime() << endl;
cout << "\n";
}
				  
+EOF

root -l -b -q run.C

rm run.C
mv a_higgsHistograms.root hhhh_$3.root

echo "End runninng"
