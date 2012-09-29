#!/usr/bin/env python

import sys,os
#sys.argv.append( '-b' )

from ROOT import *

import shutil

import config as c


period    = "2011"
sample    = "ggHZZ250"
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

print "Running on sample", sample, "with selection", selection, "in sourse ", sourcefilename, "and trigger", trigger, period, isbatch 


#tempfile = open("template_vbfZAnalyzer.C","r")
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

gROOT.LoadMacro("../plugins/rochcor.cc+");

gROOT.LoadMacro("../src/TCPhysObject.cc+");
gROOT.LoadMacro("../src/TCJet.cc+");
gROOT.LoadMacro("../src/TCMET.cc+");
gROOT.LoadMacro("../src/TCElectron.cc+");
gROOT.LoadMacro("../src/TCMuon.cc+");
gROOT.LoadMacro("../src/TCTau.cc+");
gROOT.LoadMacro("../src/TCPhoton.cc+");
gROOT.LoadMacro("../src/TCGenJet.cc+");
gROOT.LoadMacro("../src/TCGenParticle.cc+");
gROOT.LoadMacro("../src/TCPrimaryVtx.cc+");
gROOT.LoadMacro("../src/TCTriggerObject.cc+");
gROOT.LoadMacro("../plugins/WeightUtils.cc+");
gROOT.LoadMacro("../plugins/TriggerSelector.cc+");
gROOT.LoadMacro("../plugins/ZedEventsLibrary.cc+");
gROOT.LoadMacro("../plugins/HistManager.cc+");

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


timer = TStopwatch()
timer.Start()

fChain.Process("higgsAnalyzer.C+")

print "Done!", "CPU Time: ", timer.CpuTime(), "RealTime : ", timer.RealTime() 

os.system('mv a_higgsHistograms.root hhhh_'+sourcefilename+'.root')

print "End runninng"
