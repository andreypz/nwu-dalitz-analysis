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
parser.add_option("--test", dest="test", action="store_true", default=False, help="Run over the test sample")
parser.add_option("--all",  dest="all",  action="store_true", default=False, help="Run over all the samples!")
parser.add_option("--gen",  dest="gen",  action="store_true", default=False, help="Do gen-level analysis.")

(options, args) = parser.parse_args()


if len(args) < 1:
    parser.print_usage()
    exit(1)

''' Specify parameters '''
dCache      = '/pnfs/cms/WAX/11/store/user'
EOS         = '/eos/uscms/store/user'
outputPath  = EOS+'/andreypz/batch_output/zgamma/8TeV'
executable  = 'batchJob.sh'
selection = "mugamma"

DIR = EOS+'/lpchzg'
DIRNATE  = EOS+'/naodell'
DIRBRIAN = EOS+'/bpollack'
whereWeRun = ' '
if '/tthome' in os.getcwd():
    print 'We are running at NWU, yhaa'
    print os.getcwd()
    whereWeRun+='nwu'
    DIR = '/tthome/andrey'
    outputPath  = '/tthome/andrey/batch_output/zgamma/8TeV'
    DIRNATE  = '/tthome/naodell/storage/data'
    DIRBRIAN = '/tthome/bpollack/storage'


selection = options.sel
trig      = ' '+options.trigger+' '

version    = args[0]
period    = '2012'
doTest    = options.test
doData    = options.data
doBG      = options.bkg
doSignal  = options.sig
if options.gen and not options.sig:
    print "We only will do gen level analysis on Signal MC sample for now"
    sys.exit(0)
gen=" 0"
if options.gen:
    gen=" gen"
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
bg = []
signal = []



test.extend([
    cfg('dal-mad120', DIR+'/nuTuples_v9.6_8TeV/dalitz/ggHiggsToMuMuGamma_MH120',1, 'dalitz '+selection+trig+period+gen+whereWeRun),
])


if period =="2012":

    if selection in ['4mu', '2e2mu']:
        data.extend([
                cfg('DoubleMu_Run2012A',  DIRBRIAN+'/nuTuples_v9.6_8TeV/Data/DoubleMu_Run2012A',   5, 'DATA '+selection+' mumu '+' 2012 0' + whereWeRun),
                cfg('DoubleMu_Run2012B',  DIRBRIAN+'/nuTuples_v9.6_8TeV/Data/DoubleMu_Run2012B',  10, 'DATA '+selection+' mumu '+' 2012 0' + whereWeRun),
                cfg('DoubleMu_Run2012C',  DIRBRIAN+'/nuTuples_v9.6_8TeV/Data/DoubleMu_Run2012C',  15, 'DATA '+selection+' mumu '+' 2012 0' + whereWeRun),
                cfg('DoubleMu_Run2012D',  DIRBRIAN+'/nuTuples_v9.6_8TeV/Data/DoubleMu_Run2012D',  20, 'DATA '+selection+' mumu '+' 2012 0' + whereWeRun),
            ])

    if selection == ['zee','2e2mu']:
        data.extend([
                cfg('DoubleElectron_Run2012A',  DIRBRIAN+'/nuTuples_v9.6_8TeV/Data/DoubleElectron_Run2012A',   5, 'DATA '+selection+' ee '+' 2012 0' + whereWeRun),
                cfg('DoubleElectron_Run2012B',  DIRBRIAN+'/nuTuples_v9.6_8TeV/Data/DoubleElectron_Run2012B',   5, 'DATA '+selection+' ee '+' 2012 0' + whereWeRun),
                cfg('DoubleElectron_Run2012C',  DIRBRIAN+'/nuTuples_v9.6_8TeV/Data/DoubleElectron_Run2012C',  10, 'DATA '+selection+' ee '+' 2012 0' + whereWeRun),
                cfg('DoubleElectron_Run2012D',  DIRBRIAN+'/nuTuples_v9.6_8TeV/Data/DoubleElectron_Run2012D',  20, 'DATA '+selection+' ee '+' 2012 0' + whereWeRun),
                ])


    if selection in ['mugamma']:
        data.extend([
                ])

    if selection in ['mugamma','apz','2e2mu']:
        data.extend([
                cfg('MuEG_Run2012A',  DIRNATE+'/nuTuples_v9.6_8TeV/Data/MuEG_Run2012A',  5, 'DATA '+selection+' mugamma '+' 2012 0' + whereWeRun),
                cfg('MuEG_Run2012B',  DIRNATE+'/nuTuples_v9.6_8TeV/Data/MuEG_Run2012B',  8, 'DATA '+selection+' mugamma '+' 2012 0' + whereWeRun),
                cfg('MuEG_Run2012C',  DIRNATE+'/nuTuples_v9.6_8TeV/Data/MuEG_Run2012C',  8, 'DATA '+selection+' mugamma '+' 2012 0' + whereWeRun),
                cfg('MuEG_Run2012D',  DIRNATE+'/nuTuples_v9.6_8TeV/Data/MuEG_Run2012D',  15,'DATA '+selection+' mugamma '+' 2012 0' + whereWeRun),
                ])

    if selection == 'elgamma':
        data.extend([
                cfg('DoublePhoton_Run2012A', dCache+'/andreypz/nuTuples_v9_8TeV/Photon/Run2012A-22Jan2013',         10, 'DATA '+selection+' pho '+' 2012 0' + whereWeRun),
                cfg('DoublePhoton_Run2012B', dCache+'/andreypz/nuTuples_v9_8TeV/DoublePhoton/Run2012B-22Jan2013',   15, 'DATA '+selection+' pho '+' 2012 0' + whereWeRun),
                cfg('DoublePhoton_Run2012C', dCache+'/andreypz/nuTuples_v9_8TeV/DoublePhoton/Run2012C-22Jan2013',   20, 'DATA '+selection+' pho '+' 2012 0' + whereWeRun),
                cfg('DoublePhoton_Run2012D', dCache+'/andreypz/nuTuples_v9_8TeV/DoublePhoton/Run2012D-22Jan2013',   30, 'DATA '+selection+' pho '+' 2012 0' + whereWeRun),
            ])


    bg.extend([
            cfg('DYjets50',  DIRBRIAN+'/nuTuples_v9.6_8TeV/MC/DYJetsToLL_M-50_RD1',  20, 'DYjets50 ' +selection+trig+'  2012 0' + whereWeRun ),
            cfg('ZGToLLG',   DIRBRIAN+'/nuTuples_v9.6_8TeV/MC/ZGToLLG_RD1',          5, 'ZGToLLG  ' +selection+trig+'  2012 0' + whereWeRun ),
        ])

    #if selection == 'elgamma':
    #    bg.extend([
    #            ])

    if selection in ['mugamma','4mu','2e2mu','apz']:
        bg.extend([
                cfg('DYjetDalitz', DIR+'/nuTuples_v9.6_8TeV/dalitz/DYtoMuMuJet', 1, 'DY '+selection+trig+'  2012 0' + whereWeRun ),
                ])

    signal.extend([
        cfg('ggHZG-125', DIRBRIAN+'/nuTuples_v9.6_8TeV/MC/ggHZG_M125_RD1',1, 'HZG '+selection+trig+' '+period+gen + whereWeRun),
        ])

    if selection in ["mugamma"]:
        signal.extend([

            cfg('ggH-mad120', DIR+'/nuTuples_v9.6_8TeV/dalitz/ggHiggsToMuMuGamma_MH120',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
            cfg('ggH-mad125', DIR+'/nuTuples_v9.6_8TeV/dalitz/ggHiggsToMuMuGamma_MH125',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
            cfg('ggH-mad130', DIR+'/nuTuples_v9.6_8TeV/dalitz/ggHiggsToMuMuGamma_MH130',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
            cfg('ggH-mad135', DIR+'/nuTuples_v9.6_8TeV/dalitz/ggHiggsToMuMuGamma_MH135',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
            cfg('ggH-mad140', DIR+'/nuTuples_v9.6_8TeV/dalitz/ggHiggsToMuMuGamma_MH140',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
            cfg('ggH-mad145', DIR+'/nuTuples_v9.6_8TeV/dalitz/ggHiggsToMuMuGamma_MH145',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
            cfg('ggH-mad150', DIR+'/nuTuples_v9.6_8TeV/dalitz/ggHiggsToMuMuGamma_MH150',1, 'dalitz '+selection+trig+period+gen + whereWeRun),

            cfg('vbf-mad120', DIR+'/nuTuples_v9.6_8TeV/dalitz/vbfHiggsToMuMuGamma_MH120',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
            cfg('vbf-mad125', DIR+'/nuTuples_v9.6_8TeV/dalitz/vbfHiggsToMuMuGamma_MH125',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
            cfg('vbf-mad130', DIR+'/nuTuples_v9.6_8TeV/dalitz/vbfHiggsToMuMuGamma_MH130',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
            cfg('vbf-mad135', DIR+'/nuTuples_v9.6_8TeV/dalitz/vbfHiggsToMuMuGamma_MH135',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
            cfg('vbf-mad140', DIR+'/nuTuples_v9.6_8TeV/dalitz/vbfHiggsToMuMuGamma_MH140',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
            cfg('vbf-mad145', DIR+'/nuTuples_v9.6_8TeV/dalitz/vbfHiggsToMuMuGamma_MH145',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
            cfg('vbf-mad150', DIR+'/nuTuples_v9.6_8TeV/dalitz/vbfHiggsToMuMuGamma_MH150',1, 'dalitz '+selection+trig+period+gen + whereWeRun),

            cfg('vh-mad120', DIR+'/nuTuples_v9.6_8TeV/dalitz/VHiggsToMuMuGamma_MH120',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
            cfg('vh-mad125', DIR+'/nuTuples_v9.6_8TeV/dalitz/VHiggsToMuMuGamma_MH125',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
            cfg('vh-mad130', DIR+'/nuTuples_v9.6_8TeV/dalitz/VHiggsToMuMuGamma_MH130',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
            cfg('vh-mad135', DIR+'/nuTuples_v9.6_8TeV/dalitz/VHiggsToMuMuGamma_MH135',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
            cfg('vh-mad140', DIR+'/nuTuples_v9.6_8TeV/dalitz/VHiggsToMuMuGamma_MH140',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
            cfg('vh-mad145', DIR+'/nuTuples_v9.6_8TeV/dalitz/VHiggsToMuMuGamma_MH145',1, 'dalitz '+selection+trig+period+gen + whereWeRun),
            cfg('vh-mad150', DIR+'/nuTuples_v9.6_8TeV/dalitz/VHiggsToMuMuGamma_MH150',1, 'dalitz '+selection+trig+period+gen + whereWeRun),


            ])
    elif selection=="elgamma":
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
