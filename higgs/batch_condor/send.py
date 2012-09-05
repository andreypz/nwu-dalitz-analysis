#! /usr/bin/env python
import BatchMaster as b
import sys

nargs = len(sys.argv)
print "Running", sys.argv[0], nargs


''' Specify parameters '''
dCache      = '/pnfs/cms/WAX/11/store/user' 
outputPath  = '/uscms_data/d2/andreypz/cmssw/higgs2/CMSSW_4_4_4/src/NWU/Higgs/higgs/batch_condor'

selection = 'muon'
if nargs>1 and sys.argv[1]=="ele":
    selection="electron"
    
period    = '2011'
doTest    = 1
doData    = 0
doBG      = 0
doSignal  = 0

''' 
    Set job configurations.  The order of arguments is:
    (Dataset, path to data, number of jobs, arguments to pass to executable, output directory name)
'''

test = []
test.extend([
b.JobConfig('vbfZ', dCache+'/naodell/March18/MC/VBFZ', 10, 'vbfZ 0 '+selection+' '+period, selection),
#b.JobConfig('WW', '/eos/uscms/store/user/bpollack/May15/MC/WWJets/', 1, 'WW 0 '+selection+' '+period, selection),
#b.JobConfig('ggHZZ300', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ300', 1, 'ggHZZ300 0 '+selection+' '+period, selection),
#b.JobConfig('DoubleEle_Run2011A', dCache+'/andreypz/nuTuples_v1_7TeV/DoubleElectron_HZZ_Run2011A', 30, 'DATA 16,17,18 electron 2011A', selection),
#b.JobConfig('DoubleEle_Run2011B', dCache+'/andreypz/nuTuples_v1_7TeV/DoubleElectron_HZZ_Run2011B', 30, 'DATA 16,17,18 electron 2011B', selection),
#b.JobConfig('DYjets', dCache+'/andreypz/nuTuples_v1_7TeV/DYJetsToLL', 30, 'DYjets 0 '+selection+' '+period, selection),
])



if selection == 'muon':
    data = []
    data.extend([
        b.JobConfig('DoubleMu_Run2011A', dCache+'/andreypz/nuTuples_v1_7TeV/DoubleMu_HZZ_Run2011A', 25, 'DATA 4,5,8 muon 2011A', selection),
        b.JobConfig('DoubleMu_Run2011B', dCache+'/andreypz/nuTuples_v1_7TeV/DoubleMu_HZZ_Run2011B', 25, 'DATA 4,5,8 muon 2011B', selection),
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
        b.JobConfig('Photon_Aug05', dCache+'/radek/Data3/Photon_Run2011A_05Aug2011-v1', 20, 'DATA  14,15,16,17,18,19,20,21,22,23 '+selection+' 2011', selection)
        ])

if selection == 'muEG':
    data = []
    if period in ['2011',]:
        data.extend([
        b.JobConfig('Aug05', dCache+'/andreypz/Oct21/MuEG_Aug05', 10, 'DATA  24,25,27,28 muEG 2011', selection),
         ])
 
bg = []
bg.extend([
b.JobConfig('DYjets', dCache+'/andreypz/nuTuples_v1_7TeV/DYJetsToLL', 30, 'DYjets 0 '+selection+' '+period, selection),
b.JobConfig('ZZ', dCache+'/andreypz/nuTuples_v1_7TeV/ZZJetsTo2L2Nu', 1, 'ZZ 0 '+selection+' '+period, selection),
#b.JobConfig('ZZ', dCache+'/andreypz/nuTuples_v1_7TeV/ZZ_v3', 1, 'ZZ 0 '+selection+' '+period, selection),
b.JobConfig('WW', '/eos/uscms/store/user/bpollack/May15/MC/WWJets/', 1, 'WW 0 '+selection+' '+period, selection),
b.JobConfig('WZ', '/eos/uscms/store/user/bpollack/Apr17/MC/WZJets', 1, 'WZ 0 '+selection+' '+period, selection),
b.JobConfig('ttbar', dCache+'/andreypz/nuTuples_v1_7TeV/TTJets', 10, 'ttbar 0 '+selection+' '+period, selection),
b.JobConfig('tW', dCache+'/radek/MC3b/T_tW-channel-DR-r2', 1, 'tW 0 '+selection+' '+period, selection),
b.JobConfig('tbarW', dCache+'/radek/MC3b/Tbar_tW-channel-DR-r2', 1, 'tbarW 0 '+selection+' '+period, selection),
#b.JobConfig('ggWW', dCache+'/andreypz/nuTuples_v1_7TeV/GluGluToWWTo4L', 1, 'ggWW 0 '+selection+' '+period, selection),
b.JobConfig('ZZpythia', '/eos/uscms/store/user/andreypz/nuTuples_v1_7TeV/ZZ_v3', 1, 'ZZ 0 '+selection+' '+period, selection),
])


signal = []

signal.extend([
b.JobConfig('ggHZZ125', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ125', 1, 'ggHZZ125 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ200', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ200', 1, 'ggHZZ200 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ250', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ250', 1, 'ggHZZ250 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ300', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ300', 1, 'ggHZZ300 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ350', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ350', 1, 'ggHZZ350 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ400', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ400', 1, 'ggHZZ400 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ450', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ450', 1, 'ggHZZ450 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ500', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ500', 1, 'ggHZZ500 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ550', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ550', 1, 'ggHZZ550 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ600', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ600', 1, 'ggHZZ600 0 '+selection+' '+period, selection),

#b.JobConfig('ggHWW125', dCache+'/andreypz/nuTuples_v1_7TeV/ggHWW125', 1, 'ggHWW125 0 '+selection+' '+period, selection),
b.JobConfig('ggHWW200', dCache+'/andreypz/nuTuples_v1_7TeV/ggHWW200', 1, 'ggHWW200 0 '+selection+' '+period, selection),

b.JobConfig('VBFHZZ125', dCache+'/andreypz/nuTuples_v1_7TeV/VBFHZZ125', 1, 'VBFHZZ125 0 '+selection+' '+period, selection),
b.JobConfig('VBFHZZ200', dCache+'/andreypz/nuTuples_v1_7TeV/VBFHZZ200', 1, 'VBFHZZ200 0 '+selection+' '+period, selection),
b.JobConfig('VBFHZZ250', dCache+'/andreypz/nuTuples_v1_7TeV/VBFHZZ250', 1, 'VBFHZZ250 0 '+selection+' '+period, selection),

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
