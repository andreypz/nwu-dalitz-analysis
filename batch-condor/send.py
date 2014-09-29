#!/usr/bin/env python

import BatchMaster as b
import sys, os
cfg = b.JobConfig

from optparse import OptionParser
parser = OptionParser(usage="usage: %prog [options -e], --data, --bkg, --sig], -p 2011] version")
parser.add_option("-c", "--clean", dest="clean",  action="store_true", default=False,
                  help="Clean the directory with histogram output.")
parser.add_option("-s", '--sel', dest="sel", type="string", default='mugamma',
                  help="Selection to be used. Options are: '4mu','2e2mu', 'zee','mugamma', 'elgamma'")
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

(opt, args) = parser.parse_args()


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

selection = opt.sel
if selection not in ['mugamma','elgamma','jp-mugamma','zee','mumu','nate']:
  print 'this selection is not supported:',selection
  sys.exit(0)

trig      = ' '+opt.trigger+' '

version   = args[0]
period    = '2012'
doTest    = opt.test
doData    = opt.data
doBG      = opt.bkg
doSignal  = opt.sig
doQCD     = opt.qcd

if opt.gen and not doSignal and not doTest:
  print "We only will do gen level analysis on Signal MC sample for now"
  sys.exit(0)
gen = " 0"
if opt.gen:
  gen = " gen"
if opt.all:
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
    #cfg('ZGDalitz',     DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/DYtoMuMuGamma',     1, 'ZG     '+selection+trig+'  2012 0' + whereWeRun ),
    #cfg('ggH-mad125', DIR+'/nuTuples_v9.8.1_8TeV/dalitz-ee/ggHiggsToEEGamma_MH125', 1, 'dalitz '+selection+trig+period+gen + whereWeRun),
    #cfg('ggH-mad125', DIR+'/nuTuples_v9.8.1_8TeV/dalitz-el/ggHiggsToEEGamma_MH125', 1, 'dalitz '+selection+trig+period+gen + whereWeRun),
    #cfg('ggH-mad125', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/ggHiggsToMuMuGamma_MH125',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
    cfg('ggHZG-125',  DIR+'/nuTuples_v9.8_8TeV/MC/ggHZG_M125_RD1',1, 'HZG '+selection+trig+' '+period+gen + whereWeRun),
    #cfg('ggH-mad150', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/ggHiggsToMuMuGamma_MH150',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
    #cfg('HiggsToJPsiGamma', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/HiggsToJPsiGamma',      1, 'hjp '+selection+trig+period+gen + whereWeRun),
  ])


if period =="2012":

  # Data configs

  if selection in ['4mu', '2e2mu','mumu','nate']:
    data.extend([
        cfg('DoubleMu_Run2012A',  DIR+'/nuTuples_v9.8_8TeV/Data/DoubleMu_Run2012A',   5, 'DATA '+selection+' mumu '+' 2012A 0' + whereWeRun),
        cfg('DoubleMu_Run2012B',  DIR+'/nuTuples_v9.8_8TeV/Data/DoubleMu_Run2012B',   8, 'DATA '+selection+' mumu '+' 2012B 0' + whereWeRun),
        cfg('DoubleMu_Run2012C',  DIR+'/nuTuples_v9.8_8TeV/Data/DoubleMu_Run2012C',  12, 'DATA '+selection+' mumu '+' 2012C 0' + whereWeRun),
        cfg('DoubleMu_Run2012D',  DIR+'/nuTuples_v9.8_8TeV/Data/DoubleMu_Run2012D',  16, 'DATA '+selection+' mumu '+' 2012D 0' + whereWeRun),
        ])

  if selection in ['zee','2e2mu']:
    data.extend([
        cfg('DoubleElectron_Run2012A',  DIRBRIAN+'/nuTuples_v9.6_8TeV/Data/DoubleElectron_Run2012A',   5, 'DATA '+selection+' ee '+' 2012 0' + whereWeRun),
        cfg('DoubleElectron_Run2012B',  DIRBRIAN+'/nuTuples_v9.6_8TeV/Data/DoubleElectron_Run2012B',   8, 'DATA '+selection+' ee '+' 2012 0' + whereWeRun),
        cfg('DoubleElectron_Run2012C',  DIRBRIAN+'/nuTuples_v9.6_8TeV/Data/DoubleElectron_Run2012C',  12, 'DATA '+selection+' ee '+' 2012 0' + whereWeRun),
        cfg('DoubleElectron_Run2012D',  DIRBRIAN+'/nuTuples_v9.6_8TeV/Data/DoubleElectron_Run2012D',  20, 'DATA '+selection+' ee '+' 2012 0' + whereWeRun),
        ])

  if selection in ['elgamma']:
    data.extend([
        cfg('Photon_Run2012A',        DIR+'/nuTuples_v9.8_8TeV/Data/Photon_2012A',        10, 'DATA '+selection+trig+' 2012 0' + whereWeRun),
        cfg('DoublePhoton_Run2012B',  DIR+'/nuTuples_v9.8_8TeV/Data/DoublePhoton_2012B',  15, 'DATA '+selection+trig+' 2012 0' + whereWeRun),
        cfg('DoublePhoton_Run2012C',  DIR+'/nuTuples_v9.8_8TeV/Data/DoublePhoton_2012C',  30, 'DATA '+selection+trig+' 2012 0' + whereWeRun),
        cfg('DoublePhoton_Run2012D',  DIR+'/nuTuples_v9.8_8TeV/Data/DoublePhoton_2012D',  30, 'DATA '+selection+trig+' 2012 0' + whereWeRun)
        #cfg('DoubleElectron_Run2012A',  DIR+'/nuTuples_v9.8_8TeV/Data/DoubleElectron_Run2012A',   5, 'DATA '+selection+trig+' 2012 0' + whereWeRun),
        ])

  if selection in ['mugamma','jp-mugamma','2e2mu']:
    data.extend([
        cfg('MuEG_Run2012A', DIR+'/nuTuples_v9.8_8TeV/Data/MuEG_Run2012A',  5, 'DATA '+selection+' mugamma '+' 2012A 0'+ whereWeRun),
        cfg('MuEG_Run2012B', DIR+'/nuTuples_v9.8_8TeV/Data/MuEG_Run2012B',  8, 'DATA '+selection+' mugamma '+' 2012B 0'+ whereWeRun),
        cfg('MuEG_Run2012C', DIR+'/nuTuples_v9.8_8TeV/Data/MuEG_Run2012C', 12, 'DATA '+selection+' mugamma '+' 2012C 0'+ whereWeRun),
        cfg('MuEG_Run2012D', DIR+'/nuTuples_v9.8_8TeV/Data/MuEG_Run2012D', 15, 'DATA '+selection+' mugamma '+' 2012D 0'+ whereWeRun),
        ])

# Background configs

  params = selection+trig+'  2012 0' + whereWeRun;
  bg.extend([
      #cfg('DYJets10',  DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/DYJetsToLL_M-10To50filter', 15, 'DYJets ' +params),
      #cfg('DYJets50',  DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/DYJetsToLL_M-50',            25, 'DYJets ' +params),
      #cfg('ZGToLLG',   DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/ZGToLLG',                10, 'ZG     ' +params),

      #cfg('DYJets50-RD1',   DIR+'/nuTuples_v9.8_8TeV/MC/DYJetsToLL_M-50_RD1',  25, 'DYJets ' +params),
      #cfg('ZGToLLG-RD1',    DIR+'/nuTuples_v9.8_8TeV/MC/ZGToLLG_RD1',          10, 'ZG     ' +params),
      #cfg('DYJetsPow20-RD1',DIR+'/nuTuples_v9.8_8TeV/MC/DYToMuMu_M-20_RD1',    35, 'DYJets ' +params),

      ])

  if selection in ['mugamma','jp-mugamma','4mu','2e2mu','apz']:
    bg.extend([
        cfg('DYJetsDalitz', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/DYtoMuMuJet',    1, 'DYJets '+params),
        cfg('ZGDalitz',     DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/DYtoMuMuGamma',  1, 'ZG     '+params),
        #cfg('ZGDalitz-OLD', DIRME+'/nuTuples_v9.6_8TeV/dalitz/DYtoMuMuGamma', 1, 'ZG     '+params),
        ])
  elif selection in ['elgamma']:
    bg.extend([
        cfg('DYJetsDalitz', DIR+'/nuTuples_v9.8_8TeV/dalitz-el/DYtoEEJet',    1, 'DYJets '+params),
        cfg('ZGDalitz',     DIR+'/nuTuples_v9.8_8TeV/dalitz-el/DYtoEEGamma',  1, 'ZG     '+params),
        #cfg('ZGDalitz-OLD', DIRME+'/nuTuples_v9.6_8TeV/dalitz/DYtoMuMuGamma', 1, 'ZG     '+params),
        ])

  if doQCD:
    params = 'QCD '+selection+trig+'  2012 0' + whereWeRun;
    bg.extend([
        cfg('QCD_EM_Pt_20to30',   DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/QCD_Pt_20_30_EMEnriched',  15, params),
        cfg('QCD_EM_Pt_30to80',   DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/QCD_Pt_30_80_EMEnriched',  15, params),
        cfg('QCD_EM_Pt_80to170',  DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/QCD_Pt_80_170_EMEnriched', 15, params),
        cfg('QCD_EM_Pt_170to250', DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/QCD_Pt_170_250_EMEnriched',20, params),
        cfg('QCD_EM_Pt_250to350', DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/QCD_Pt_250_350_EMEnriched',20, params),
        cfg('QCD_EM_Pt_350',      DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/QCD_Pt_350_EMEnriched',    20, params),
        ])
    if selection in ['mugamma','jp-mugamma']:
      bg.extend([
          cfg('QCD_Mu_Pt_30to50',  DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/QCD_Pt_30to50_bEnriched_MuEnrichedPt_14', 15, params),
          cfg('QCD_Mu_Pt_50to150', DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/QCD_Pt_50to150_bEnriched_MuEnrichedPt_14',15, params),
          cfg('QCD_Mu_Pt_150',     DIR+'/nuTuples_v9.8_8TeV/MC_skimmed/QCD_Pt_150_bEnriched_MuEnrichedPt_14',    15, params),
          ])


# Signal MC configs

  signal.extend([
      cfg('ggHZG-125',  DIR+'/nuTuples_v9.8_8TeV/MC/ggHZG_M125_RD1',1, 'HZG '+selection+trig+' '+period+gen + whereWeRun),
      ])


  if selection in ['mugamma','jp-mugamma']:
    params = selection+trig+period+gen + whereWeRun;
    signal.extend([
        cfg('ZtoJPsiGamma',     DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/ZtoJPsiGamma-MuMuGamma',1, 'zjp '+params),
        cfg('HiggsToJPsiGamma', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/HiggsToJPsiGamma',      1, 'hjp '+params),
        cfg('JPsiToMuMu-S10',   DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/JPsiToMuMu_S10',        1,  'jp '+params),
        cfg('JPsiToMuMu_Pt20',  DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/JPsiToMuMu_Pt20_RD2',   1,  'jp '+params),

        ])
  if selection in ['mugamma']:
    params = 'dalitz '+selection+trig+period+gen + whereWeRun;
    signal.extend([
        cfg('ggH-mad120', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/ggHiggsToMuMuGamma_MH120',1, params),
        cfg('ggH-mad125', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/ggHiggsToMuMuGamma_MH125',1, params),
        cfg('ggH-mad130', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/ggHiggsToMuMuGamma_MH130',1, params),
        cfg('ggH-mad135', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/ggHiggsToMuMuGamma_MH135',1, params),
        cfg('ggH-mad140', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/ggHiggsToMuMuGamma_MH140',1, params),
        cfg('ggH-mad145', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/ggHiggsToMuMuGamma_MH145',1, params),
        cfg('ggH-mad150', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/ggHiggsToMuMuGamma_MH150',1, params),

        cfg('ggH-mad120-rd1', DIR+'/nuTuples_v9.8_8TeV/dalitz-rd1/ggHiggsToMuMuGamma_MH120',1, params),
        cfg('ggH-mad125-rd1', DIR+'/nuTuples_v9.8_8TeV/dalitz-rd1/ggHiggsToMuMuGamma_MH125',1, params),
        cfg('ggH-mad130-rd1', DIR+'/nuTuples_v9.8_8TeV/dalitz-rd1/ggHiggsToMuMuGamma_MH130',1, params),
        cfg('ggH-mad135-rd1', DIR+'/nuTuples_v9.8_8TeV/dalitz-rd1/ggHiggsToMuMuGamma_MH135',1, params),
        cfg('ggH-mad140-rd1', DIR+'/nuTuples_v9.8_8TeV/dalitz-rd1/ggHiggsToMuMuGamma_MH140',1, params),
        cfg('ggH-mad145-rd1', DIR+'/nuTuples_v9.8_8TeV/dalitz-rd1/ggHiggsToMuMuGamma_MH145',1, params),
        cfg('ggH-mad150-rd1', DIR+'/nuTuples_v9.8_8TeV/dalitz-rd1/ggHiggsToMuMuGamma_MH150',1, params),

        cfg('vbf-mad120', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/vbfHiggsToMuMuGamma_MH120',1, params),
        cfg('vbf-mad125', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/vbfHiggsToMuMuGamma_MH125',1, params),
        cfg('vbf-mad130', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/vbfHiggsToMuMuGamma_MH130',1, params),
        cfg('vbf-mad135', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/vbfHiggsToMuMuGamma_MH135',1, params),
        cfg('vbf-mad140', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/vbfHiggsToMuMuGamma_MH140',1, params),
        cfg('vbf-mad145', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/vbfHiggsToMuMuGamma_MH145',1, params),
        cfg('vbf-mad150', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/vbfHiggsToMuMuGamma_MH150',1, params),

        cfg('vh-mad120', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/VHiggsToMuMuGamma_MH120',1, params),
        cfg('vh-mad125', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/VHiggsToMuMuGamma_MH125',1, params),
        cfg('vh-mad130', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/VHiggsToMuMuGamma_MH130',1, params),
        cfg('vh-mad135', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/VHiggsToMuMuGamma_MH135',1, params),
        cfg('vh-mad140', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/VHiggsToMuMuGamma_MH140',1, params),
        cfg('vh-mad145', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/VHiggsToMuMuGamma_MH145',1, params),
        cfg('vh-mad150', DIR+'/nuTuples_v9.8_8TeV/dalitz-mu/VHiggsToMuMuGamma_MH150',1, params),
        ])

  elif selection in ["elgamma"]:
    params = 'dalitz '+selection+trig+period+gen + whereWeRun;

    signal.extend([
        cfg('ggH-mad120', DIR+'/nuTuples_v9.8_8TeV/dalitz-el/ggHiggsToEEGamma_MH120', 1, params),
        cfg('ggH-mad125', DIR+'/nuTuples_v9.8_8TeV/dalitz-el/ggHiggsToEEGamma_MH125', 1, params),
        cfg('ggH-mad130', DIR+'/nuTuples_v9.8_8TeV/dalitz-el/ggHiggsToEEGamma_MH130', 1, params),
        cfg('ggH-mad135', DIR+'/nuTuples_v9.8_8TeV/dalitz-el/ggHiggsToEEGamma_MH135', 1, params),
        cfg('ggH-mad140', DIR+'/nuTuples_v9.8_8TeV/dalitz-el/ggHiggsToEEGamma_MH140', 1, params),
        cfg('ggH-mad145', DIR+'/nuTuples_v9.8_8TeV/dalitz-el/ggHiggsToEEGamma_MH145', 1, params),
        cfg('ggH-mad150', DIR+'/nuTuples_v9.8_8TeV/dalitz-el/ggHiggsToEEGamma_MH150', 1, params),

        cfg('vbf-mad120', DIR+'/nuTuples_v9.8_8TeV/dalitz-el/vbfHiggsToEEGamma_MH120', 1, params),
        cfg('vbf-mad125', DIR+'/nuTuples_v9.8_8TeV/dalitz-el/vbfHiggsToEEGamma_MH125', 1, params),
        cfg('vbf-mad130', DIR+'/nuTuples_v9.8_8TeV/dalitz-el/vbfHiggsToEEGamma_MH130', 1, params),
        cfg('vbf-mad135', DIR+'/nuTuples_v9.8_8TeV/dalitz-el/vbfHiggsToEEGamma_MH135', 1, params),
        cfg('vbf-mad140', DIR+'/nuTuples_v9.8_8TeV/dalitz-el/vbfHiggsToEEGamma_MH140', 1, params),
        cfg('vbf-mad145', DIR+'/nuTuples_v9.8_8TeV/dalitz-el/vbfHiggsToEEGamma_MH145', 1, params),
        cfg('vbf-mad150', DIR+'/nuTuples_v9.8_8TeV/dalitz-el/vbfHiggsToEEGamma_MH150', 1, params),
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
