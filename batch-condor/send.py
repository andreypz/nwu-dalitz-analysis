#!/usr/bin/env python

import BatchMaster as b
import sys, os
cfg = b.JobConfig

from optparse import OptionParser
parser = OptionParser(usage="usage: %prog [options -e], --data, --bkg, --sig], -p 2011] version")
parser.add_option("-c", "--clean", dest="clean",  action="store_true", default=False,
                  help="Clean the directory with histogram output.")
parser.add_option("-e", "--elgamma", dest="elgamma", action="store_true", default=False,
                  help="Use electron selection (by default it will run muon selection)")
parser.add_option("--mugamma",   dest="mugamma",  action="store_true", default=False,
                  help="Use Mu+Photon trigger (for running on MuEG path)")
parser.add_option("--mumu",      dest="mumu",     action="store_true", default=False, help="Use Double-Mu trigger")
parser.add_option("--singlemu",  dest="singlemu", action="store_true", default=False, help="Use Iso-mu trigger")
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
    

if options.elgamma:
    selection="elgamma"
if options.mugamma:
    selection="mugamma"
if options.mumu:
    selection="mumu"
if options.singlemu:
    selection="single-mu"
                                            
    
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
    #cfg('SingleMu_Run2012A',        dCache+'/andreypz/nuTuples_v6_8TeV/SingleMu/Run2012A-22Jan2013',         5, 'DATA '+selection+' 2012 0' + whereWeRun),
    #cfg('MuEG_Run2012A',        dCache+'/andreypz/nuTuples_v6_8TeV/MuEG/Run2012A-22Jan2013',         5, 'DATA mugamma 2012 0' + whereWeRun),
    cfg('dal-mad120', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToMuMuGamma_MH120',1, 'dalitz mugamma '+period+gen),

])


if period =="2012":
    
    if selection == 'mumu':
        data.extend([
            cfg('DoubleMu_Run2012A',  DIR+'/bpollack/V09_05_8TeV/DoubleMu/Run2012A',  5, 'DATA '+selection+' 2012 0' + whereWeRun),
            cfg('DoubleMu_Run2012B',  DIR+'/bpollack/V09_05_8TeV/DoubleMu/Run2012B',  8, 'DATA '+selection+' 2012 0' + whereWeRun),
            cfg('DoubleMu_Run2012C',  DIR+'/bpollack/V09_05_8TeV/DoubleMu/Run2012C', 10, 'DATA '+selection+' 2012 0' + whereWeRun),
            cfg('DoubleMu_Run2012D',  DIR+'/bpollack/V09_05_8TeV/DoubleMu/Run2012D', 10, 'DATA '+selection+' 2012 0' + whereWeRun),
            ])

    if selection == 'single-mu':
        data.extend([
            cfg('SingleMu_Run2012A',  dCache+'/andreypz/nuTuples_v6_8TeV/SingleMu/Run2012A-22Jan2013',     5, 'DATA '+selection+' 2012 0' + whereWeRun),
            cfg('SingleMu_Run2012B',  dCache+'/andreypz/nuTuples_v6_8TeV/SingleMu/Run2012B-22Jan2013',     8, 'DATA '+selection+' 2012 0' + whereWeRun),
            cfg('SingleMu_Run2012C',  dCache+'/andreypz/nuTuples_v6_8TeV/SingleMu/Run2012C-22Jan2013-v1',  8, 'DATA '+selection+' 2012 0' + whereWeRun),
            cfg('SingleMu_Run2012D',  dCache+'/andreypz/nuTuples_v6_8TeV/SingleMu/Run2012D-22Jan2013',    15, 'DATA '+selection+' 2012 0' + whereWeRun),
            ])
        
    if selection == 'mugamma':
    
        data.extend([
            
            cfg('MuEG_Run2012A',  DIRNATE+'/nuTuples_v9.6_8TeV/Data/MuEG_Run2012A',  5, 'DATA '+selection+' 2012 0' + whereWeRun),
            cfg('MuEG_Run2012B',  DIRNATE+'/nuTuples_v9.6_8TeV/Data/MuEG_Run2012B',  8, 'DATA '+selection+' 2012 0' + whereWeRun),
            cfg('MuEG_Run2012C',  DIRNATE+'/nuTuples_v9.6_8TeV/Data/MuEG_Run2012C',  8, 'DATA '+selection+' 2012 0' + whereWeRun),
            cfg('MuEG_Run2012D',  DIRNATE+'/nuTuples_v9.6_8TeV/Data/MuEG_Run2012D',  15,'DATA '+selection+' 2012 0' + whereWeRun),

            ])
    
    if selection == 'elgamma':
        data.extend([
            cfg('DoublePhoton_Run2012A', dCache+'/andreypz/nuTuples_v9_8TeV/Photon/Run2012A-22Jan2013',         10, 'DATA '+selection+' 2012 0' + whereWeRun),
            cfg('DoublePhoton_Run2012B', dCache+'/andreypz/nuTuples_v9_8TeV/DoublePhoton/Run2012B-22Jan2013',   15, 'DATA '+selection+' 2012 0' + whereWeRun),
            cfg('DoublePhoton_Run2012C', dCache+'/andreypz/nuTuples_v9_8TeV/DoublePhoton/Run2012C-22Jan2013',   20, 'DATA '+selection+' 2012 0' + whereWeRun),
            cfg('DoublePhoton_Run2012D', dCache+'/andreypz/nuTuples_v9_8TeV/DoublePhoton/Run2012D-22Jan2013',   30, 'DATA '+selection+' 2012 0' + whereWeRun),
            ])
        
        
    bg.extend([
        cfg('DYjets50',  DIR+'/bpollack/V09_05_8TeV/MC/DYJetsToLL_M-50',  1, 'DYjets50 ' +selection+' '+period),
        ])

    if selection == 'elgamma':
        bg.extend([
            #cfg('DYjets10', dCache+'/andreypz/nuTuples_v9_8TeV/DYJets_M10To50', 20, 'DYjets10 '+selection+' '+period),
            #cfg('DYjets50',  dCache+'/andreypz/nuTuples_v9_8TeV/DYJets_M0To50',  1, 'DYjets0 ' +selection+' '+period),
            ])

    if selection == 'mugamma':
        bg.extend([
            #cfg('DYdalitz', dCache+'/andreypz/nutuples_v9.4_8TeV/DYMuMuGammaDalitz', 1, 'DY '+selection+' '+period),
            ])
    """
    cfg('WZJetsTo3LNu',  dCache+'/andreypz/nuTuples_v5_8TeV/WZJetsTo3LNu',  1, 'WZ '+selection+' '+period),
    cfg('WWJetsTo2L2Nu', dCache+'/andreypz/nuTuples_v5_8TeV/WWJetsTo2L2Nu', 1, 'WW '+selection+' '+period),
    cfg('ZZJetsTo2L2Nu', dCache+'/andreypz/nuTuples_v5_8TeV/ZZJetsTo2L2Nu', 1, 'ZZ '+selection+' '+period),
    cfg('ttbar',         dCache+'/andreypz/nuTuples_v5_8TeV/TTJets',       10, 'ttbar '+selection+' '+period),
    cfg('tW',            dCache+'/andreypz/nuTuples_v5_8TeV/tW',            1, 'tW '+selection+' '+period),
    cfg('tbarW',         dCache+'/andreypz/nuTuples_v5_8TeV/tbarW',         1, 'tbarW '+selection+' '+period),
    cfg('DYjets',        dCache+'/andreypz/nuTuples_v5_8TeV/DYjets',       30, 'DYjets '+selection+' '+period),
    cfg('DYjets10',      dCache+'/naodell/nuTuples_v5_8TeV/DYJetsToLL_M-10To50', 20, 'DYjets10 '+selection+' '+period),
    cfg('vbfZ',          dCache+'/andreypz/nuTuples_v5_8TeV/Zvbf',          1, 'vbfZ '+selection+' '+period),
    
    cfg('WJetsToLNu',    dCache+'/naodell/nuTuples_v5_8TeV/WJetsToLNu',     5, 'WJetsToLNu '+selection+' '+period),
    cfg('WZJetsTo2L2Q',  dCache+'/naodell/nuTuples_v5_8TeV/WZJetsTo2L2Q',   1, 'WZJetsTo2L2Q '+selection+' '+period),
        cfg('ZZJetsTo2L2Q',  dCache+'/naodell/nuTuples_v5_8TeV/ZZJetsTo2L2Q',   1, 'ZZJetsTo2L2Q '+selection+' '+period),
        """
        
        
    signal.extend([
        cfg('ggHZG-125', DIRBRIAN+'/nuTuples_v9.6_8TeV/MC/ggHZG_M125_RD1',1, 'HZG '+selection+' '+period+gen + whereWeRun),
        ])
    
    if selection in ["mumu","mugamma","single-mu"]:
        signal.extend([
            
            cfg('ggH-mad120', DIR+'/nuTuples_v9.6_8TeV/dalitz/ggHiggsToMuMuGamma_MH120',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('ggH-mad125', DIR+'/nuTuples_v9.6_8TeV/dalitz/ggHiggsToMuMuGamma_MH125',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('ggH-mad130', DIR+'/nuTuples_v9.6_8TeV/dalitz/ggHiggsToMuMuGamma_MH130',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('ggH-mad135', DIR+'/nuTuples_v9.6_8TeV/dalitz/ggHiggsToMuMuGamma_MH135',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('ggH-mad140', DIR+'/nuTuples_v9.6_8TeV/dalitz/ggHiggsToMuMuGamma_MH140',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('ggH-mad145', DIR+'/nuTuples_v9.6_8TeV/dalitz/ggHiggsToMuMuGamma_MH145',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('ggH-mad150', DIR+'/nuTuples_v9.6_8TeV/dalitz/ggHiggsToMuMuGamma_MH150',1, 'dalitz '+selection+' '+period+gen + whereWeRun),

            cfg('vbf-mad120', DIR+'/nuTuples_v9.6_8TeV/dalitz/vbfHiggsToMuMuGamma_MH120',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('vbf-mad125', DIR+'/nuTuples_v9.6_8TeV/dalitz/vbfHiggsToMuMuGamma_MH125',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('vbf-mad130', DIR+'/nuTuples_v9.6_8TeV/dalitz/vbfHiggsToMuMuGamma_MH130',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('vbf-mad135', DIR+'/nuTuples_v9.6_8TeV/dalitz/vbfHiggsToMuMuGamma_MH135',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('vbf-mad140', DIR+'/nuTuples_v9.6_8TeV/dalitz/vbfHiggsToMuMuGamma_MH140',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('vbf-mad145', DIR+'/nuTuples_v9.6_8TeV/dalitz/vbfHiggsToMuMuGamma_MH145',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('vbf-mad150', DIR+'/nuTuples_v9.6_8TeV/dalitz/vbfHiggsToMuMuGamma_MH150',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            
            cfg('vh-mad120', DIR+'/nuTuples_v9.6_8TeV/dalitz/VHiggsToMuMuGamma_MH120',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('vh-mad125', DIR+'/nuTuples_v9.6_8TeV/dalitz/VHiggsToMuMuGamma_MH125',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('vh-mad130', DIR+'/nuTuples_v9.6_8TeV/dalitz/VHiggsToMuMuGamma_MH130',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('vh-mad135', DIR+'/nuTuples_v9.6_8TeV/dalitz/VHiggsToMuMuGamma_MH135',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('vh-mad140', DIR+'/nuTuples_v9.6_8TeV/dalitz/VHiggsToMuMuGamma_MH140',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('vh-mad145', DIR+'/nuTuples_v9.6_8TeV/dalitz/VHiggsToMuMuGamma_MH145',1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('vh-mad150', DIR+'/nuTuples_v9.6_8TeV/dalitz/VHiggsToMuMuGamma_MH150',1, 'dalitz '+selection+' '+period+gen + whereWeRun),


            ])
    elif selection=="elgamma":
        signal.extend([
            cfg('ggH-mad120', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH120', 1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('ggH-mad125', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH125', 1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('ggH-mad130', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH130', 1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('ggH-mad135', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH135', 1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('ggH-mad140', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH140', 1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('ggH-mad145', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH145', 1, 'dalitz '+selection+' '+period+gen + whereWeRun),
            cfg('ggH-mad150', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH150', 1, 'dalitz '+selection+' '+period+gen + whereWeRun),
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
    print "Submitting to batch?"
    batcher.submit_to_batch()
