#!/usr/bin/env python
import sys,os
#sys.argv.append( '-b' )

from optparse import OptionParser
parser = OptionParser(usage="usage: %prog [options -e], -m], -p 2011] sample sourcefile")
parser.add_option("-e", "--ele", dest="electron", action="store_true", default=False,
                  help="Use electron selection (by default it will run muon selection)")
parser.add_option("-t","--trigger", dest="trigger", type="string", default="none",
                  help="Select a trigger to run. options are: double-mu, single-mu, mu-pho, pho, single-el")
parser.add_option("-p", "--period", dest="period", default="2012", help="Set data period (2011/2012)")
parser.add_option("-b", "--batch", dest="batch", action="store_true", default=False,
                  help="Run in batch (uning the scripts in batch_condor).That would imply the use of input.txt")
parser.add_option("-c", "--clean", dest="clean", action="store_true", default=False, help="Clean the libs")

(options, args) = parser.parse_args()

if options.clean:
    os.system("rm ../src/*.so ../src/*.d ../plugins/*.so ../plugins/*.d")
    print "All cleaned"
    exit(0)
    
if len(args) < 2:
    parser.print_usage()
    exit(1)

sample    = args[0]
sourcefilename = args[1] 

period    = options.period
selection = "mu"
if options.electron:
    print "options.electorn", options.electron
    selection="el"
trigger = options.trigger

isbatch   = options.batch    
if(isbatch):
    sourceFiles = "./input.txt"
    print "Do batch, source files: \n", sourceFiles
else:
    sourceFiles = "./"+sourcefilename
    print "Run locally, source files: \n", sourceFiles

    
print "Running on sample", sample, "with selection", selection, "and trigger", trigger,\
      "in source ", sourcefilename, ", period =", period, " batch =",  isbatch


tempfile = open("template_zgamma.C","r")
whole_thing = tempfile.read()
whole_thing = whole_thing.replace("@SELECTION",selection)
whole_thing = whole_thing.replace("@TRIGGER",trigger)
tempfile.close()

cfile = open("zgamma.C","w")
cfile.write(whole_thing)
cfile.close()

from ROOT import *

#gSystem.Load("../plugins/lib/libShapeLine.so");
#gROOT.LoadMacro("../plugins/rochcor.cc+");
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

#gROOT.LoadMacro("../plugins/WeightUtils.cc+");
gROOT.LoadMacro("../plugins/TriggerSelector.cc+");
gROOT.LoadMacro("../plugins/ZGAngles.cc+");
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

fChain.Process("zgamma.C+")

print "Done!", "CPU Time: ", timer.CpuTime(), "RealTime : ", timer.RealTime()

os.system('mv a_'+selection+'_higgsHistograms.root hhhh_'+sourcefilename+'.root')

print "End runninng"
os.system("rm zgamma.C")
