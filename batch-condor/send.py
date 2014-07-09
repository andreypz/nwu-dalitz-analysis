#!/usr/bin/env python

import BatchMaster as b
import sys, os
cfg = b.JobConfig

from optparse import OptionParser
parser = OptionParser(usage="usage: %prog [options -e], --data, --bkg, --sig], -p 2011] version")
parser.add_option("-c", "--clean", dest="clean",  action="store_true", default=False,
                  help="Clean the directory with histogram output.")
parser.add_option("-s", '--sel', dest="sel", type="string", default='mugamma',
                  help="Selection to be used. Options are: '4mu','2e2mu', 'zee','mugamma', 'egamma'")
parser.add_option("-t","--trigger", dest="trigger", type="string", default="mugamma",
                  help="Select a trigger to run. It only matters here for MC samples. For Data the right trigger is hardcoded.\
Options are: mumu, single-mu, mugamma, pho, single-el")
parser.add_option("--data", dest="data", action="store_true", default=False, help="Run over the data sample")
parser.add_option("--sig",  dest="sig",  action="store_true", default=False, help="Run over the signal sample")
parser.add_option("--bkg",  dest="bkg",  action="store_true", default=False, help="Run over the background samples")
parser.add_option("--qcd",  dest="qcd",  action="store_true", default=False, help="Run over the QCD samples")
parser.add_option("--test", dest="test", action="store_true", default=False, help="Run over the test sample")
parser.add_option("--all",  dest="all",  action="store_true", default=False, help="Run over all the samples!")
parser.add_option("--gen",  dest="gen",  action="store_true", default=False, help="Do gen-level analysis.")

(options, args) = parser.parse_args()


if len(args) < 1:
  parser.print_usage()
  exit(1)

''' Specify parameters '''
EOS         = '/eos/uscms/store/user'
outputPath  = EOS+'/andreypz/batch_output/zgamma/8TeV'
executable  = 'batchJob.sh'
selection   = "mugamma"

DIR      = EOS+'/lpchzg'
DIRNATE  = EOS+'/naodell'
DIRBRIAN = EOS+'/bpollack/storage'
whereWeRun = ' '

if '/tthome' in os.getcwd():
  print 'We are running at NWU, yhaa'
  print os.getcwd()
  whereWeRun+='nwu'
  DIR   = '/tthome/share/noobs/'
  DIRME = '/tthome/andrey'
  DIRNATE     = '/tthome/naodell/storage/data'
  DIRBRIAN    = '/tthome/bpollack/storage'
  outputPath  = '/tthome/andrey/batch_output/zgamma/8TeV'

selection = options.sel
if selection not in ['mugamma','egamma','jp-mugamma','zee','mumu','nate']:
  print 'this selection is not supported:',selection
  sys.exit(0)

trig      = ' '+options.trigger+' '

version   = args[0]
period    = '2012'
doTest    = options.test
doData    = options.data
doBG      = options.bkg
doSignal  = options.sig
doQCD     = options.qcd

if options.gen and not options.sig:
  print "We only will do gen level analysis on Signal MC sample for now"
  sys.exit(0)
gen = " 0"
if options.gen:
  gen = " gen"
if options.all:
  doData = 1
  doBG   = 1
  doSignal = 1

'''
Set job configurations.  The order of arguments is:
(Dataset, path to data, number of jobs, arguments to pass to executable, output directory name)
'''

test = []
data = []
bg   = []
signal = []

test.extend([
    #cfg('DYJetsPow20-RD1',  DIR+'/nuTuples_v9.8_8TeV/MC/DYToMuMu_M-20_RD1', 25, 'DYJets '+selection+trig+'  2012 0' + whereWeRun ),
    #cfg('ZGDalitz',     DIR+'/nuTuples_v9.8_8TeV/dalitz2/DYtoMuMuGamma',     1, 'ZG     '+selection+trig+'  2012 0' + whereWeRun ),
    #cfg('dal-mad120', DIR+'/nuTuples_v9.6_8TeV/dalitz/ggHiggsToMuMuGamma_MH120',1, 'dalitz '+selection+trig+period+gen+whereWeRun),
    #cfg('MuEG_Run2012D',  DIRNATE+'/nuTuples_v9.6_8TeV/Data/MuEG_Run2012D', 15, 'DATA '+selection+' mugamma '+' 2012 0' + whereWeRun),
    cfg('ggH-mad150', DIR+'/nuTuples_v9.8_8TeV/dalitz2/ggHiggsToMuMuGamma_MH150',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
    ])


if period =="2012":

  # Data configs

  if selection in ['4mu', '2e2mu','mumu','nate']:
    data.extend([
        cfg('DoubleMu_Run2012A',  DIR+'/nuTuples_v9.8_8TeV/Data/DoubleMu_Run2012A',   5, 'DATA '+selection+' mumu '+' 2012 0' + whereWeRun),
        cfg('DoubleMu_Run2012B',  DIR+'/nuTuples_v9.8_8TeV/Data/DoubleMu_Run2012B',   8, 'DATA '+selection+' mumu '+' 2012 0' + whereWeRun),
        cfg('DoubleMu_Run2012C',  DIR+'/nuTuples_v9.8_8TeV/Data/DoubleMu_Run2012C',  12, 'DATA '+selection+' mumu '+' 2012 0' + whereWeRun),
        cfg('DoubleMu_Run2012D',  DIR+'/nuTuples_v9.8_8TeV/Data/DoubleMu_Run2012D',  16, 'DATA '+selection+' mumu '+' 2012 0' + whereWeRun),
        ])

  if selection in ['zee','2e2mu']:
    data.extend([
        cfg('DoubleElectron_Run2012A',  DIRBRIAN+'/nuTuples_v9.6_8TeV/Data/DoubleElectron_Run2012A',   5, 'DATA '+selection+' ee '+' 2012 0' + whereWeRun),
        cfg('DoubleElectron_Run2012B',  DIRBRIAN+'/nuTuples_v9.6_8TeV/Data/DoubleElectron_Run2012B',   8, 'DATA '+selection+' ee '+' 2012 0' + whereWeRun),
        cfg('DoubleElectron_Run2012C',  DIRBRIAN+'/nuTuples_v9.6_8TeV/Data/DoubleElectron_Run2012C',  12, 'DATA '+selection+' ee '+' 2012 0' + whereWeRun),
        cfg('DoubleElectron_Run2012D',  DIRBRIAN+'/nuTuples_v9.6_8TeV/Data/DoubleElectron_Run2012D',  20, 'DATA '+selection+' ee '+' 2012 0' + whereWeRun),
        ])


  if selection in ['mugamma','jp-mugamma','2e2mu']:
    data.extend([
        cfg('MuEG_Run2012A',  DIR+'/nuTuples_v9.8_8TeV/Data/MuEG_Run2012A',  5, 'DATA '+selection+' mugamma '+' 2012 0' + whereWeRun),
        cfg('MuEG_Run2012B',  DIR+'/nuTuples_v9.8_8TeV/Data/MuEG_Run2012B',  8, 'DATA '+selection+' mugamma '+' 2012 0' + whereWeRun),
        cfg('MuEG_Run2012C',  DIR+'/nuTuples_v9.8_8TeV/Data/MuEG_Run2012C', 12, 'DATA '+selection+' mugamma '+' 2012 0' + whereWeRun),
        cfg('MuEG_Run2012D',  DIR+'/nuTuples_v9.8_8TeV/Data/MuEG_Run2012D', 15, 'DATA '+selection+' mugamma '+' 2012 0' + whereWeRun),
        ])

  if selection in ['elgamma']:
    data.extend([
        cfg('DoublePhoton_Run2012A', dCache+'/andreypz/nuTuples_v9_8TeV/Photon/Run2012A-22Jan2013',         10, 'DATA '+selection+' pho '+' 2012 0' + whereWeRun),
        cfg('DoublePhoton_Run2012B', dCache+'/andreypz/nuTuples_v9_8TeV/DoublePhoton/Run2012B-22Jan2013',   15, 'DATA '+selection+' pho '+' 2012 0' + whereWeRun),
        cfg('DoublePhoton_Run2012C', dCache+'/andreypz/nuTuples_v9_8TeV/DoublePhoton/Run2012C-22Jan2013',   20, 'DATA '+selection+' pho '+' 2012 0' + whereWeRun),
        cfg('DoublePhoton_Run2012D', dCache+'/andreypz/nuTuples_v9_8TeV/DoublePhoton/Run2012D-22Jan2013',   30, 'DATA '+selection+' pho '+' 2012 0' + whereWeRun),
        ])

# Background configs

  bg.extend([
      #cfg('DYJets10',  DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/DYJetsToLL_M-10To50filter', 15, 'DYJets ' +selection+trig+'  2012 0' + whereWeRun ),
      #cfg('DYJets50',  DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/DYJetsToLL_M-50',            25, 'DYJets ' +selection+trig+'  2012 0' + whereWeRun ),
      #cfg('ZGToLLG',   DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/ZGToLLG',                10, 'ZG     ' +selection+trig+'  2012 0' + whereWeRun ),

      cfg('DYJets50-RD1',  DIR+'/nuTuples_v9.8_8TeV/MC/DYJetsToLL_M-50_RD1',  25, 'DYJets ' +selection+trig+'  2012 0' + whereWeRun ),
      cfg('ZGToLLG-RD1',   DIR+'/nuTuples_v9.8_8TeV/MC/ZGToLLG_RD1',          10, 'ZG     ' +selection+trig+'  2012 0' + whereWeRun ),
      cfg('DYJetsPow20-RD1',  DIR+'/nuTuples_v9.8_8TeV/MC/DYToMuMu_M-20_RD1',    25, 'DYJets ' +selection+trig+'  2012 0' + whereWeRun ),

      ])


  if selection in ['mugamma','jp-mugamma','4mu','2e2mu','apz']:
    bg.extend([
        cfg('DYJetsDalitz', DIR+'/nuTuples_v9.8_8TeV/dalitz2/DYtoMuMuJet',    1, 'DYJets '+selection+trig+'  2012 0' + whereWeRun ),
        cfg('ZGDalitz',     DIR+'/nuTuples_v9.8_8TeV/dalitz2/DYtoMuMuGamma',  1, 'ZG     '+selection+trig+'  2012 0' + whereWeRun ),
        #cfg('ZGDalitz-OLD', DIRME+'/nuTuples_v9.6_8TeV/dalitz/DYtoMuMuGamma', 1, 'ZG     '+selection+trig+'  2012 0' + whereWeRun ),
        ])

  if doQCD:
   bg.extend([
        cfg('QCD_Pt_20_30',   DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/QCD_Pt_20_30_EMEnriched',  5, 'QCD '+selection+trig+'  2012 0' + whereWeRun ),
        cfg('QCD_Pt_30_80',   DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/QCD_Pt_30_80_EMEnriched',  5, 'QCD '+selection+trig+'  2012 0' + whereWeRun ),
        cfg('QCD_Pt_80_170',  DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/QCD_Pt_80_170_EMEnriched', 5, 'QCD '+selection+trig+'  2012 0' + whereWeRun ),
        cfg('QCD_Pt_170_250', DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/QCD_Pt_170_250_EMEnriched',5, 'QCD '+selection+trig+'  2012 0' + whereWeRun ),
        cfg('QCD_Pt_250_350', DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/QCD_Pt_250_350_EMEnriched',5, 'QCD '+selection+trig+'  2012 0' + whereWeRun ),
        cfg('QCD_Pt_350',     DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/QCD_Pt_350_EMEnriched',    5, 'QCD '+selection+trig+'  2012 0' + whereWeRun ),
        ])

# Signal MC configs
  signal.extend([
      cfg('ggHZG-125',  DIR+'/nuTuples_v9.8_8TeV/MC/ggHZG_M125_RD1',1, 'HZG '+selection+trig+' '+period+gen + whereWeRun),
      cfg('ggH-mad125', DIR+'/nuTuples_v9.8_8TeV/dalitz2/ggHiggsToMuMuGamma_MH125',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
      ])

  if selection in ['mugamma','jp-mugamma']:
    signal.extend([
        cfg('ZtoJPsiGamma',     DIR+'/nuTuples_v9.8_8TeV/dalitz2/ZtoJPsiGamma-MuMuGamma',1, 'zjp '+selection+trig+period+gen + whereWeRun),
        cfg('HiggsToJPsiGamma', DIR+'/nuTuples_v9.8_8TeV/dalitz2/HiggsToJPsiGamma',      1, 'hjp '+selection+trig+period+gen + whereWeRun),
        ])
  if selection in ['mugamma']:
    signal.extend([
        cfg('ggH-mad120', DIR+'/nuTuples_v9.8_8TeV/dalitz2/ggHiggsToMuMuGamma_MH120',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        cfg('ggH-mad130', DIR+'/nuTuples_v9.8_8TeV/dalitz2/ggHiggsToMuMuGamma_MH130',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        cfg('ggH-mad135', DIR+'/nuTuples_v9.8_8TeV/dalitz2/ggHiggsToMuMuGamma_MH135',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        cfg('ggH-mad140', DIR+'/nuTuples_v9.8_8TeV/dalitz2/ggHiggsToMuMuGamma_MH140',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        cfg('ggH-mad145', DIR+'/nuTuples_v9.8_8TeV/dalitz2/ggHiggsToMuMuGamma_MH145',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        cfg('ggH-mad150', DIR+'/nuTuples_v9.8_8TeV/dalitz2/ggHiggsToMuMuGamma_MH150',1, 'dalitz '+selection+trig+period+gen + whereWeRun),

        cfg('vbf-mad120', DIR+'/nuTuples_v9.8_8TeV/dalitz2/vbfHiggsToMuMuGamma_MH120',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        cfg('vbf-mad125', DIR+'/nuTuples_v9.8_8TeV/dalitz2/vbfHiggsToMuMuGamma_MH125',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        cfg('vbf-mad130', DIR+'/nuTuples_v9.8_8TeV/dalitz2/vbfHiggsToMuMuGamma_MH130',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        cfg('vbf-mad135', DIR+'/nuTuples_v9.8_8TeV/dalitz2/vbfHiggsToMuMuGamma_MH135',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        cfg('vbf-mad140', DIR+'/nuTuples_v9.8_8TeV/dalitz2/vbfHiggsToMuMuGamma_MH140',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        #cfg('vbf-mad145', DIR+'/nuTuples_v9.8_8TeV/dalitz2/vbfHiggsToMuMuGamma_MH145',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        #cfg('vbf-mad150', DIR+'/nuTuples_v9.8_8TeV/dalitz2/vbfHiggsToMuMuGamma_MH150',1, 'dalitz '+selection+trig+period+gen + whereWeRun),

        #cfg('vh-mad120', DIR+'/nuTuples_v9.8_8TeV/dalitz2/VHiggsToMuMuGamma_MH120',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        #cfg('vh-mad125', DIR+'/nuTuples_v9.8_8TeV/dalitz2/VHiggsToMuMuGamma_MH125',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        #cfg('vh-mad130', DIR+'/nuTuples_v9.8_8TeV/dalitz2/VHiggsToMuMuGamma_MH130',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        #cfg('vh-mad135', DIR+'/nuTuples_v9.8_8TeV/dalitz2/VHiggsToMuMuGamma_MH135',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        cfg('vh-mad140', DIR+'/nuTuples_v9.8_8TeV/dalitz2/VHiggsToMuMuGamma_MH140',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        cfg('vh-mad145', DIR+'/nuTuples_v9.8_8TeV/dalitz2/VHiggsToMuMuGamma_MH145',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        cfg('vh-mad150', DIR+'/nuTuples_v9.8_8TeV/dalitz2/VHiggsToMuMuGamma_MH150',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        ])

  elif selection in ["elgamma"]:
    signal.extend([
        cfg('ggH-mad120', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH120', 1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        cfg('ggH-mad125', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH125', 1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        cfg('ggH-mad130', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH130', 1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        cfg('ggH-mad135', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH135', 1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        cfg('ggH-mad140', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH140', 1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        cfg('ggH-mad145', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH145', 1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        cfg('ggH-mad150', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH150', 1, 'dalitz '+selection+trig+period+gen + whereWeRun),
        ])

else:
  print "Only 2012! Other periods are not supported"

inputSamples = []

if doTest:
  inputSamples.extend(test)
if doData:
  inputSamples.extend(data)
if doBG:
  inputSamples.extend(bg)
if doSignal:
  inputSamples.extend(signal)

if len(inputSamples) is not 0:
  batcher = b.BatchMaster(inputSamples, outputPath, shortQueue = False, stageDir = '../StageBatch_'+version,
                          executable = executable, prefix = version + '/' + selection +"_"+ period)
  print "Submitting to batch"
  batcher.submit_to_batch()
