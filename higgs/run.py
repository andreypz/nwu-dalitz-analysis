#!/usr/bin/env python

import sys,os
#sys.argv.append( '-b' )

from ROOT import *

import shutil

import config as c
#shutil.copy2('template_higgsAnalyzer.C', 'higgsAnalyzer.C')


period    = "2011"
sample    = "ggHZZ400"
selection = "muon"
trigger   = "0"
isbatch   = ""
nargs = len(sys.argv)
print sys.argv[0], nargs
if nargs>=1:
    sourcefilename = sample
if nargs>=2:
    sample    = sys.argv[1]
    sourcefilename = sample
if nargs>=3:
    selection = sys.argv[2]
if nargs>=4:
    sourcefilename =  sys.argv[3]
if nargs>=5:
    trigger   = sys.argv[4] 
if nargs>=6:
    period    = sys.argv[5]      
if nargs>=7:
    isbatch   = sys.argv[6]

print sample, selection, sourcefilename, trigger, period, isbatch 


tempfile = open("template_higgsAnalyzer.C","r")
whole_thing = tempfile.read()
whole_thing = whole_thing.replace("SUFFIX", sample)
whole_thing = whole_thing.replace("TRIGGER", str(trigger))
whole_thing = whole_thing.replace("SELECTION",selection)
whole_thing = whole_thing.replace("PERIOD",period)
tempfile.close()

cfile = open("higgsAnalyzer.C","w")
cfile.write(whole_thing)
cfile.close()

if(isbatch == "b"):
    sourceFiles = "./input.txt"
    print "Do batch, source files: \n", sourceFiles
else:
    sourceFiles = "./sourceFiles/"+sourcefilename+".txt"
    print "Run locally, source files: \n", sourceFiles



gSystem.Load("../plugins/libShapeLine.so");

gROOT.LoadMacro("../src/TCJet.cc+");
gROOT.LoadMacro("../src/TCMET.cc+");
gROOT.LoadMacro("../src/TCElectron.cc+");
gROOT.LoadMacro("../src/TCMuon.cc+");
gROOT.LoadMacro("../src/TCTau.cc+");
gROOT.LoadMacro("../src/TCPhoton.cc+");
gROOT.LoadMacro("../src/TCGenJet.cc+");
gROOT.LoadMacro("../src/TCGenParticle.cc+");
gROOT.LoadMacro("../src/TCPrimaryVtx.cc+");
#gROOT.LoadMacro("../src/TCTrigger.cc+");
gROOT.LoadMacro("../src/TCTriggerObject.cc+");
gROOT.LoadMacro("../plugins/WeightUtils.cc+");
gROOT.LoadMacro("../plugins/TriggerSelector.cc+");
gROOT.LoadMacro("../plugins/ZedEventsLibrary.cc+");

fChain = TChain("ntupleProducer/eventTree");

f = open(sourceFiles,"r")
nfiles=0

for l in f.readlines():
    line = l.split('\n')[0]
    if line.strip():
        print line
        fChain.Add(line)      
        nfiles+=1
print nfiles, " files added!"
#fChain.Add("/eos/uscms/store/user/bpollack/May15/MC/ZZJets/nuTuple_100_1_RGT.root")
#fChain.Add("/eos/uscms/store/user/bpollack/Apr17/MC/WZJets/nuTuple_9_1_IRt.root")

para = c.SampleParams()
xsec_and_colors = para.xsec_and_colors()

lineColor = xsec_and_colors[sample][0]
fillColor = xsec_and_colors[sample][1]
cs = xsec_and_colors[sample][2]
Ne = xsec_and_colors[sample][3]


params = open("./params.txt","w")
params.write(selection+" "+sample+" "+str(Ne)+" "+str(cs)+" "+str(fillColor)+" "+str(lineColor))
params.close();

timer = TStopwatch()
timer.Start()

fChain.Process("higgsAnalyzer.C+")

print "Done!", "CPU Time: ", timer.CpuTime(), "RealTime : ", timer.RealTime() 

"""

run = '''
using namespace std;

void run() {

gSystem->Load("../plugins/libShapeLine.so");

gROOT->LoadMacro("../src/TCJet.cc+");
gROOT->LoadMacro("../src/TCMET.cc+");
gROOT->LoadMacro("../src/TCElectron.cc+");
gROOT->LoadMacro("../src/TCMuon.cc+");
gROOT->LoadMacro("../src/TCTau.cc+");
gROOT->LoadMacro("../src/TCPhoton.cc+");
gROOT->LoadMacro("../src/TCGenJet.cc+");
gROOT->LoadMacro("../src/TCGenParticle.cc+");
gROOT->LoadMacro("../src/TCPrimaryVtx.cc+");
//#gROOT->LoadMacro("../src/TCTrigger.cc+");
gROOT->LoadMacro("../src/TCTriggerObject.cc+");
gROOT->LoadMacro("../plugins/WeightUtils.cc+");
gROOT->LoadMacro("../plugins/TriggerSelector.cc+");
gROOT->LoadMacro("../plugins/ZedEventsLibrary.cc+");

TChain* fChain = new TChain("ntupleProducer/eventTree");

ifstream sourceFiles("SOURCEFILES");
string line;
int  count = 0;
while (sourceFiles >> line) {
  fChain->Add(line.c_str());
  ++count;
}

fChain->Process("higgsAnalyzer.C+");  
}
'''

run = run.replace("SOURCEFILES",sourceFiles)

ff = open("run.C","w")
ff.write(run)
ff.close()

os.system('root -b -q run.C')
"""

os.system('mv a_higgsHistograms.root hhhh_$3.root')

print "End runninng"
