#!/usr/bin/env python

import BatchMaster as b
import sys
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
#outputPath  = dCache+'/andreypz/batch_output/zgamma/8TeV'
outputPath  = EOS+'/andreypz/batch_output/zgamma/8TeV'
executable  = 'batchJob.sh'
selection = "mugamma"

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
gen=" "
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

test.extend([
    #cfg('SingleMu_Run2012A',        dCache+'/andreypz/nuTuples_v6_8TeV/SingleMu/Run2012A-22Jan2013',         5, 'DATA '+selection+' 2012'),
    #cfg('MuEG_Run2012A',        dCache+'/andreypz/nuTuples_v6_8TeV/MuEG/Run2012A-22Jan2013',         5, 'DATA mugamma 2012'),
    cfg('dal-mad120', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToMuMuGamma_MH120',1, 'dalitz mugamma '+period+gen),

])


if period =="2012":
    
    if selection == 'mumu':
        data.extend([
            cfg('DoubleMu_Run2012A',  EOS+'/bpollack/V09_04_8TeV/DoubleMu/Run2012A',  5, 'DATA '+selection+' 2012'),
            cfg('DoubleMu_Run2012B',  EOS+'/bpollack/V09_04_8TeV/DoubleMu/Run2012B',  8, 'DATA '+selection+' 2012'),
            cfg('DoubleMu_Run2012C',  EOS+'/bpollack/V09_04_8TeV/DoubleMu/Run2012C', 10, 'DATA '+selection+' 2012'),
            cfg('DoubleMu_Run2012D',  EOS+'/bpollack/V09_04_8TeV/DoubleMu/Run2012D', 10, 'DATA '+selection+' 2012'),
            ])

    if selection == 'single-mu':
        data.extend([
            cfg('SingleMu_Run2012A',  dCache+'/andreypz/nuTuples_v6_8TeV/SingleMu/Run2012A-22Jan2013',     5, 'DATA '+selection+' 2012'),
            cfg('SingleMu_Run2012B',  dCache+'/andreypz/nuTuples_v6_8TeV/SingleMu/Run2012B-22Jan2013',     8, 'DATA '+selection+' 2012'),
            cfg('SingleMu_Run2012C',  dCache+'/andreypz/nuTuples_v6_8TeV/SingleMu/Run2012C-22Jan2013-v1',  8, 'DATA '+selection+' 2012'),
            cfg('SingleMu_Run2012D',  dCache+'/andreypz/nuTuples_v6_8TeV/SingleMu/Run2012D-22Jan2013',    15, 'DATA '+selection+' 2012'),
            ])
        
    if selection == 'mugamma':
    
        data.extend([
            cfg('MuEG_Run2012A',  dCache+'/andreypz/nutuples_v9.3_8TeV/MuEG/Run2012A-22Jan2013',  5, 'DATA '+selection+' 2012'),
            cfg('MuEG_Run2012B',  dCache+'/andreypz/nutuples_v9.3_8TeV/MuEG/Run2012B-22Jan2013',  8, 'DATA '+selection+' 2012'),
            cfg('MuEG_Run2012C',  dCache+'/andreypz/nutuples_v9.3_8TeV/MuEG/Run2012C-22Jan2013',  8, 'DATA '+selection+' 2012'),
            cfg('MuEG_Run2012D',  dCache+'/andreypz/nutuples_v9.3_8TeV/MuEG/Run2012D-22Jan2013',  15,'DATA '+selection+' 2012'),

            ])
    
    if selection == 'elgamma':
        data.extend([
            cfg('DoublePhoton_Run2012A', dCache+'/andreypz/nuTuples_v9_8TeV/Photon/Run2012A-22Jan2013',         10, 'DATA '+selection+' 2012'),
            cfg('DoublePhoton_Run2012B', dCache+'/andreypz/nuTuples_v9_8TeV/DoublePhoton/Run2012B-22Jan2013',   15, 'DATA '+selection+' 2012'),
            cfg('DoublePhoton_Run2012C', dCache+'/andreypz/nuTuples_v9_8TeV/DoublePhoton/Run2012C-22Jan2013',   20, 'DATA '+selection+' 2012'),
            cfg('DoublePhoton_Run2012D', dCache+'/andreypz/nuTuples_v9_8TeV/DoublePhoton/Run2012D-22Jan2013',   30, 'DATA '+selection+' 2012'),
            ])
        
        
    bg = []
    
    bg.extend([
        cfg('DYjets50',  EOS+'/bpollack/V09_04_8TeV/MC/DYJetsToLL_M-50',  1, 'DYjets50 ' +selection+' '+period),
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
        
        

    signal = []

    signal.extend([
        cfg('ggHZG-125', EOS+'/bpollack/V09_04_8TeV/MC/ggHZG_M125_S10',1, 'HZG '+selection+' '+period+gen),
        ])
    
    if selection in ["mumu","mugamma","single-mu"]:
        signal.extend([
            
            #cfg('dal-MCFM125',dCache+'/andreypz/nuTuples_v9.4_8TeV/MCFM_dalitz_v2',    1, 'dalitz '+selection+' '+period),
            cfg('dal-mad120', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToMuMuGamma_MH120',1, 'dalitz '+selection+' '+period+gen),
            cfg('dal-mad125', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToMuMuGamma_MH125',1, 'dalitz '+selection+' '+period+gen),
            cfg('dal-mad130', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToMuMuGamma_MH130',1, 'dalitz '+selection+' '+period+gen),
            cfg('dal-mad135', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToMuMuGamma_MH135',1, 'dalitz '+selection+' '+period+gen),
            cfg('dal-mad140', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToMuMuGamma_MH140',1, 'dalitz '+selection+' '+period+gen),
            cfg('dal-mad145', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToMuMuGamma_MH145',1, 'dalitz '+selection+' '+period+gen),
            cfg('dal-mad150', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToMuMuGamma_MH150',1, 'dalitz '+selection+' '+period+gen),

            cfg('vbf-mad120', dCache+'/andreypz/nuTuples_v9.4_8TeV/vbfHiggsToMuMuGamma_MH120',1, 'dalitz '+selection+' '+period+gen),
            cfg('vbf-mad125', dCache+'/andreypz/nuTuples_v9.4_8TeV/vbfHiggsToMuMuGamma_MH125',1, 'dalitz '+selection+' '+period+gen),
            cfg('vbf-mad130', dCache+'/andreypz/nuTuples_v9.4_8TeV/vbfHiggsToMuMuGamma_MH130',1, 'dalitz '+selection+' '+period+gen),
            cfg('vbf-mad135', dCache+'/andreypz/nuTuples_v9.4_8TeV/vbfHiggsToMuMuGamma_MH135',1, 'dalitz '+selection+' '+period+gen),
            cfg('vbf-mad140', dCache+'/andreypz/nuTuples_v9.4_8TeV/vbfHiggsToMuMuGamma_MH140',1, 'dalitz '+selection+' '+period+gen),
            cfg('vbf-mad145', dCache+'/andreypz/nuTuples_v9.4_8TeV/vbfHiggsToMuMuGamma_MH145',1, 'dalitz '+selection+' '+period+gen),
            cfg('vbf-mad150', dCache+'/andreypz/nuTuples_v9.4_8TeV/vbfHiggsToMuMuGamma_MH150',1, 'dalitz '+selection+' '+period+gen),

            ])
    elif selection=="elgamma":
        signal.extend([
            cfg('dal-mad120', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH120', 1, 'dalitz '+selection+' '+period+gen),
            cfg('dal-mad125', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH125', 1, 'dalitz '+selection+' '+period+gen),
            cfg('dal-mad130', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH130', 1, 'dalitz '+selection+' '+period+gen),
            cfg('dal-mad135', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH135', 1, 'dalitz '+selection+' '+period+gen),
            cfg('dal-mad140', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH140', 1, 'dalitz '+selection+' '+period+gen),
            cfg('dal-mad145', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH145', 1, 'dalitz '+selection+' '+period+gen),
            cfg('dal-mad150', dCache+'/andreypz/nuTuples_v9.4_8TeV/HiggsToEEGamma_MH150', 1, 'dalitz '+selection+' '+period+gen),
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
                            executable = executable, selection = version + '/' + selection +"_"+ period)
    print "Submitting to batch?"
    batcher.submit_to_batch()
            
