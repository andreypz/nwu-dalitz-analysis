#!/usr/bin/env python
import sys,os
#sys.argv.append( '-b' )

from optparse import OptionParser
parser = OptionParser(usage="usage: %prog [options -e], -m], -p 2011] samplename sourcefile")
parser.add_option("-e", "--ele", dest="electron", action="store_true", default=False, help="Use electron selection")
parser.add_option("-p", "--period", dest="period", default="2012", help="Set data period (2011/2012)")
parser.add_option("-b", "--batch", dest="batch", action="store_true", default=False, help="Run in batch")
parser.add_option("--ss", dest="same_sign", action="store_true", default=False, help="Run same sign selection")

from ROOT import *
import shutil

(options, args) = parser.parse_args()

if len(args) < 2:
    parser.print_usage()
    exit(1)

sample    = args[0]
sourcefilename = args[1]

period    = options.period
selection = "muon"
do_ss = "kFALSE"
if options.same_sign:
    do_ss = "kTRUE"
    print "Doing same sign analysis"

if options.electron :
    print "options.electorn???", options.electron
    selection="electron"

isbatch   = options.batch
    
if(isbatch):
    sourceFiles = "./input.txt"
    print "Do batch, source files: \n", sourceFiles
else:
    sourceFiles = "./sourceFiles/"+sourcefilename+".txt"
    print "Run locally, source files: \n", sourceFiles



print "Running on sample", sample, "with selection", selection, "in source ", sourcefilename, ", period =", period, " batch =",  isbatch


#tempfile = open("template_vbfZAnalyzer.C","r")
tempfile = open("template_higgsAnalyzer.C","r")
whole_thing = tempfile.read()
whole_thing = whole_thing.replace("SUFFIX",sample)
whole_thing = whole_thing.replace("SELECTION",selection)
whole_thing = whole_thing.replace("PERIOD",period)
whole_thing = whole_thing.replace('"DOSAMESIGN"',do_ss)
tempfile.close()

cfile = open("higgsAnalyzer.C","w")
cfile.write(whole_thing)
cfile.close()

gSystem.Load("../plugins/lib/libShapeLine.so");
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

#os.system('mv a_higgsHistograms.root hhhh_'+sourcefilename+'.root')

print "End runninng"
os.system("rm zgamma.C")
