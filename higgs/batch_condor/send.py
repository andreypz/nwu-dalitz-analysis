#!/usr/bin/env python

import BatchMaster as b
import sys

nargs = len(sys.argv)
print "Running", sys.argv[0], nargs


''' Specify parameters '''
dCache      = '/pnfs/cms/WAX/11/store/user' 
EOS = '/eos/uscms/store/user'
outputPath  = '/uscms_data/d2/andreypz/cmssw/higgs2/CMSSW_4_4_4/src/NWU/Higgs/higgs/batch_condor'

selection = 'muon'
if nargs>1 and sys.argv[1]=="ele":
    selection="electron"
    
period    = '2011'
doTest    = 1
doData    = 1
doBG      = 1
doSignal  = 0

''' 
    Set job configurations.  The order of arguments is:
    (Dataset, path to data, number of jobs, arguments to pass to executable, output directory name)
'''

test = []
test.extend([
#b.JobConfig('vbfZ', dCache+'/naodell/March18/MC/VBFZ', 10, 'vbfZ 0 '+selection+' '+period, selection),
#b.JobConfig('WW', '/eos/uscms/store/user/bpollack/May15/MC/WWJets/', 1, 'WW 0 '+selection+' '+period, selection),
#b.JobConfig('ggHZZ250', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ250', 1, 'ggHZZ250 0 '+selection+' '+period, selection),
#b.JobConfig('ggHZZ125', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ125', 1, 'ggHZZ125 0 '+selection+' '+period, selection),
#b.JobConfig('DoubleMu_Run2011A', dCache+'/andreypz/nuTuples_v2_7TeV/DoubleMu_HZZ_Run2011A', 25, 'DATA 4,5,8 muon 2011A', selection),
#b.JobConfig('DoubleMu_Run2011B', dCache+'/andreypz/nuTuples_v2_7TeV/DoubleMu_HZZ_Run2011B', 25, 'DATA 4,5,8 muon 2011B', selection),
b.JobConfig('ggHZZ125', dCache+'/andreypz/nuTuples_v2_7TeV/ggH125', 1, 'ggHZZ125 0 '+selection+' '+period, selection),
#b.JobConfig('DYjets', dCache+'/andreypz/nuTuples_v2_7TeV/DYjets', 30, 'DYjets 0 '+selection+' '+period, selection),
])



if selection == 'muon':
    data = []
    data.extend([
        #b.JobConfig('DoubleMu_Run2011A', dCache+'/andreypz/nuTuples_v2_7TeV/DoubleMu_HZZ_Run2011A', 25, 'DATA 4,5,8 muon 2011A', selection),
        b.JobConfig('DoubleMu_Run2011B', EOS+'/andreypz/nuTuples_v2_7TeV/DoubleMu_HZZ_Run2011B', 25, 'DATA 4,5,8 muon 2011B', selection),
        ])

                            


if selection == 'electron':
    data = []
    data.extend([
        b.JobConfig('DoubleEle_Run2011A', dCache+'/andreypz/nuTuples_v1_7TeV/DoubleElectron_HZZ_Run2011A', 30, 'DATA 16,17,18 electron 2011A', selection),
        b.JobConfig('DoubleEle_Run2011B', dCache+'/andreypz/nuTuples_v1_7TeV/DoubleElectron_HZZ_Run2011B', 30, 'DATA 16,17,18 electron 2011B', selection),
        ])



if selection == 'gamma' or selection == 'muGamma' or selection == 'eGamma':
    data = []
    if period in ['2011']:
        data.extend([
        ])

if selection == 'muEG':
    data = []
    if period in ['2011',]:
        data.extend([
         ])
 
bg = []
bg.extend([
#b.JobConfig('DYjets', dCache+'/andreypz/nuTuples_v2_7TeV/DYjets', 30, 'DYjets 0 '+selection+' '+period, selection),
b.JobConfig('ZZ', dCache+'/naodell/nuTuples_v2_7TeV/ZZJetsTo2L2Nu', 1, 'ZZ 0 '+selection+' '+period, selection),
b.JobConfig('WZ', dCache+'/naodell/nuTuples_v2_7TeV/WZJetsTo3LNu', 1, 'WZ 0 '+selection+' '+period, selection),
b.JobConfig('WW', dCache+'/naodell/nuTuples_v2_7TeV/WWJetsTo2L2Nu', 1, 'WW 0 '+selection+' '+period, selection),
b.JobConfig('ttbar', dCache+'/andreypz/nuTuples_v2_7TeV/TTJets', 10, 'ttbar 0 '+selection+' '+period, selection),
b.JobConfig('tW', dCache+'/andreypz/nuTuples_v2_7TeV/tW', 1, 'tW 0 '+selection+' '+period, selection),
b.JobConfig('tbarW', dCache+'/andreypz/nuTuples_v2_7TeV/tbarW', 1, 'tbarW 0 '+selection+' '+period, selection),
#b.JobConfig('ZZpythia', '/eos/uscms/store/user/andreypz/nuTuples_v1_7TeV/ZZ_v3', 1, 'ZZ 0 '+selection+' '+period, selection),
])


signal = []

signal.extend([
b.JobConfig('ggHZZ125', dCache+'/andreypz/nuTuples_v2_7TeV/ggH125', 1, 'ggHZZ125 0 '+selection+' '+period, selection),
#b.JobConfig('ggHZZ200', dCache+'/andreypz/nuTuples_v2_7TeV/ggH200', 1, 'ggHZZ200 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ250', dCache+'/andreypz/nuTuples_v2_7TeV/ggH250', 1, 'ggHZZ250 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ300', dCache+'/andreypz/nuTuples_v2_7TeV/ggH300', 1, 'ggHZZ300 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ350', dCache+'/andreypz/nuTuples_v2_7TeV/ggH350', 1, 'ggHZZ350 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ400', dCache+'/andreypz/nuTuples_v2_7TeV/ggH400', 1, 'ggHZZ400 0 '+selection+' '+period, selection),


b.JobConfig('ggHWW125', dCache+'/andreypz/nuTuples_v2_7TeV/ggHWW130', 1, 'ggHWW125 0 '+selection+' '+period, selection),
#b.JobConfig('ggHWW200', dCache+'/andreypz/nuTuples_v2_7TeV/ggHWW200', 1, 'ggHWW200 0 '+selection+' '+period, selection),
b.JobConfig('ggHWW250', dCache+'/andreypz/nuTuples_v2_7TeV/ggHWW250', 1, 'ggHWW250 0 '+selection+' '+period, selection),
b.JobConfig('ggHWW300', dCache+'/andreypz/nuTuples_v2_7TeV/ggHWW300', 1, 'ggHWW300 0 '+selection+' '+period, selection),
b.JobConfig('ggHWW350', dCache+'/andreypz/nuTuples_v2_7TeV/ggHWW350', 1, 'ggHWW350 0 '+selection+' '+period, selection),
b.JobConfig('ggHWW400', dCache+'/andreypz/nuTuples_v2_7TeV/ggHWW400', 1, 'ggHWW400 0 '+selection+' '+period, selection),

b.JobConfig('VBFHZZ125', dCache+'/andreypz/nuTuples_v2_7TeV/VBFHZZ125', 1, 'VBFHZZ125 0 '+selection+' '+period, selection),
b.JobConfig('VBFHZZ200', dCache+'/andreypz/nuTuples_v2_7TeV/VBFHZZ200', 1, 'VBFHZZ200 0 '+selection+' '+period, selection),
b.JobConfig('VBFHZZ250', dCache+'/andreypz/nuTuples_v2_7TeV/VBFHZZ250', 1, 'VBFHZZ250 0 '+selection+' '+period, selection),
b.JobConfig('VBFHZZ300', dCache+'/andreypz/nuTuples_v2_7TeV/VBFHZZ300', 1, 'VBFHZZ300 0 '+selection+' '+period, selection),
b.JobConfig('VBFHZZ350', dCache+'/andreypz/nuTuples_v2_7TeV/VBFHZZ350', 1, 'VBFHZZ350 0 '+selection+' '+period, selection),
b.JobConfig('VBFHZZ400', dCache+'/andreypz/nuTuples_v2_7TeV/VBFHZZ400', 1, 'VBFHZZ400 0 '+selection+' '+period, selection),

])

btemp = b.BatchMaster([],outputPath)
btemp.CreateWorkingDir("code_dir")

if doData:
    batcher = b.BatchMaster(data, outputPath)
    batcher.SubmitToLPC()
if doBG:
    batcher = b.BatchMaster(bg, outputPath)
    batcher.SubmitToLPC()
if doSignal:
    batcher = b.BatchMaster(signal, outputPath)
    batcher.SubmitToLPC()
if doTest:
    batcher = b.BatchMaster(test, outputPath)
    batcher.SubmitToLPC()
