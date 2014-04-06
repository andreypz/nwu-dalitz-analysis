#!/usr/bin/env python
import sys,os
#sys.argv.append( '-b' )

from optparse import OptionParser
parser = OptionParser(usage="usage: %prog [options -e], -m], -p 2011] sample sourcefile")
parser.add_option("-s", '--sel', dest="sel", type="string", default='mugamma',
                  help="Selection to be used. Options are: '4mu','2e2mu', 'zee','mugamma', 'egamma'")
parser.add_option("-t","--trigger", dest="trigger", type="string", default="mugamma",
                  help="Select a trigger to run. options are: mumu, single-mu, mugamma, pho, single-el")
parser.add_option("-p", "--period", dest="period", default="2012", help="Set data period (2011/2012)")
parser.add_option("-b", "--batch", dest="batch", action="store_true", default=False,
                  help="Run in batch (uning the scripts in batch_condor).That would imply the use of input.txt")
parser.add_option("-c", "--clean", dest="clean", action="store_true", default=False, help="Clean the libs")
parser.add_option("-g", "--gen", dest="gen", action="store_true", default=False, help="Make gen selection")

(options, args) = parser.parse_args()

if options.clean:
    os.system("rm ../src/*.so ../src/*.d ../plugins/*.so ../plugins/*.d ./*.so")
    print "All cleaned"
    exit(0)

if len(args) < 2:
    parser.print_usage()
    exit(1)

sample    = args[0]
sourcefilename = args[1]

period    = options.period
trigger   = options.trigger
selection = "mugamma"
ana = "zgamma"

selection=options.sel
if selection=='egamma':
    ana="egamma"
if selection in ['4mu','2e2mu']:
    ana='fourLeptons'

isbatch   = options.batch

if(isbatch):
    sourceFiles = "./input.txt"
    print "Do batch, source files: \n", sourceFiles
else:
    sourceFiles = "./"+sourcefilename
    print "Run locally, source files: \n", sourceFiles


print "Running ", ana, " analyzer, on sample", sample, "with selection=", selection, ", trigger=", trigger,\
      "in source=", sourcefilename, ", period=", period, " batch=",  isbatch, " and gen selection=", options.gen

from ROOT import *

gROOT.SetMacroPath(".:../src/:../plugins/:../interface/")

gROOT.LoadMacro("HistManager.cc+")
gROOT.LoadMacro("ZGAngles.cc+")
gROOT.LoadMacro("WeightUtils.cc+")
gROOT.LoadMacro("TriggerSelector.cc+")

gROOT.LoadMacro("rochcor2012v2.C+")
#gSystem.Load("lib/libShapeLine.so")
gROOT.LoadMacro("TCPhysObject.cc+")
gROOT.LoadMacro("TCJet.cc+")
gROOT.LoadMacro("TCMET.cc+")
gROOT.LoadMacro("TCEGamma.cc+")
gROOT.LoadMacro("TCTrack.cc+")
gROOT.LoadMacro("TCElectron.cc+")
gROOT.LoadMacro("TCPhoton.cc+")
gROOT.LoadMacro("TCMuon.cc+")
gROOT.LoadMacro("TCTau.cc+")
gROOT.LoadMacro("TCGenJet.cc+")
gROOT.LoadMacro("TCGenParticle.cc+")
gROOT.LoadMacro("TCPrimaryVtx.cc+")
gROOT.LoadMacro("TCTriggerObject.cc+")
gROOT.LoadMacro("ObjectID.cc+")
gROOT.LoadMacro("HistMaker.cc+")
#gROOT.LoadMacro("PhosphorCorrectorFunctor.cc+")

fChain = TChain("ntupleProducer/eventTree")

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

fChain.Process(ana+".C+", "%s %s %s %s" % (sample,selection,trigger,str(options.gen).lower()))

os.system('mv a_'+selection+'_higgsHistograms.root hhhh_'+sourcefilename+'.root')

print "Done!", "CPU Time: ", timer.CpuTime(), "RealTime : ", timer.RealTime()

print "End runninng"
