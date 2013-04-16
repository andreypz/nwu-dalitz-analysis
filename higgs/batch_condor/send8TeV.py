#!/usr/bin/env python

import BatchMaster as b
import sys

nargs = len(sys.argv)
print "Running", sys.argv[0], nargs


''' Specify parameters '''
dCache      = '/pnfs/cms/WAX/11/store/user' 
EOS         = '/eos/uscms/store/user'
outputPath  = EOS+'/andreypz/batch_output/8TeV/'

selection = 'muon'
if nargs>1 and sys.argv[1]=="ele":
    selection="electron"
    
period    = '2012'
doTest    = 0
doData    = 1
doBG      = 1
doSignal  = 1

''' 
    Set job configurations.  The order of arguments is:
    (Dataset, path to data, number of jobs, arguments to pass to executable, output directory name)
'''

test = []
test.extend([
#b.JobConfig('ZZ', dCache+'/andreypz/nuTuples_v5_8TeV/ZZJetsTo2L2Nu', 10, 'ZZ '+selection+' '+period, selection),
#b.JobConfig('DoubleMu_Run2012A', dCache+'/devildog/nuTuples_v1_8TeV/DoubleMu_Run2012A', 1, 'DATA muon 2012', selection),
#b.JobConfig('DoubleMu_Run2012B', dCache+'/devildog/nuTuples_v1_8TeV/DoubleMu_Run2012B', 10, 'DATA muon 2012', selection),
#b.JobConfig('DoubleMu_Run2012A', dCache+'/devildog/nuTuples_v1_8TeV/DoubleMu_Run2012A', 10, 'DATA muon 2012', selection),
#b.JobConfig('DoubleMu_Run2012B', dCache+'/devildog/nuTuples_v1_8TeV/DoubleMu_Run2012B', 10, 'DATA muon 2012', selection),
#b.JobConfig('vbfZ',          dCache+'/andreypz/nuTuples_v5_8TeV/Zvbf',          1, 'vbfZ '+selection+' '+period, selection),
#b.JobConfig('DYjets10',      dCache+'/naodell/nuTuples_v5_8TeV/DYJetsToLL_M-10To50', 20, 'DYjets10 '+selection+' '+period, selection),
#b.JobConfig('WZJetsTo3LNu',  dCache+'/andreypz/nuTuples_v5_8TeV/WZJetsTo3LNu',  1, 'WZ '+selection+' '+period, selection),
#b.JobConfig('DYjets',        dCache+'/andreypz/nuTuples_v5_8TeV/DYjets',       30, 'DYjets '+selection+' '+period, selection),
    b.JobConfig('DoubleMu_Run2012D',        dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012D_v2',            20, 'DATA muon 2012', selection),
])



if selection == 'muon':
    data = []
    data.extend([
        b.JobConfig('DoubleMu_Run2012A',        dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012A_v2',         5, 'DATA muon 2012', selection),
        b.JobConfig('DoubleMu_Run2012B',        dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012B_v2',        10, 'DATA muon 2012', selection),
        b.JobConfig('DoubleMu_Run2012C_prompt', dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012C_Prompt_v2', 10, 'DATA muon 2012', selection),
        b.JobConfig('DoubleMu_Run2012C_recover',dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012C_Recover_v2', 1, 'DATA muon 2012', selection),
        b.JobConfig('DoubleMu_Run2012C_24Aug',  dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012C_24Aug_v2',   1, 'DATA muon 2012', selection),
        b.JobConfig('DoubleMu_Run2012D',        dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012D_v2',        10, 'DATA muon 2012', selection),
        #b.JobConfig('DoubleMu_Run2012A', dCache+'/devildog/nuTuples_v1_8TeV/DoubleMu_Run2012A', 10, 'DATA muon 2012', selection),
        #b.JobConfig('DoubleMu_Run2012B', dCache+'/devildog/nuTuples_v1_8TeV/DoubleMu_Run2012B', 30, 'DATA muon 2012', selection),
        ])

                            


if selection == 'electron':
    data = []
    data.extend([
        b.JobConfig('DoubleEle_Run2012A',        dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012A_v2',         10, 'DATA electron 2012', selection),
        b.JobConfig('DoubleEle_Run2012B',        dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012B_v2',         20, 'DATA electron 2012', selection),
        b.JobConfig('DoubleEle_Run2012C_prompt', dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012C_Prompt_v2',  20, 'DATA electron 2012', selection),
        b.JobConfig('DoubleEle_Run2012C_recover',dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012C_Recover_v2', 20, 'DATA electron 2012', selection),
        b.JobConfig('DoubleEle_Run2012C_24Aug',  dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012C_24Aug_v2',   20, 'DATA electron 2012', selection),
        b.JobConfig('DoubleEle_Run2012D',        dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012D_v2',         20, 'DATA electron 2012', selection),

        #b.JobConfig('DoubleEle_Run2012A', dCache+'/naodell/nuTuples_v5_8TeV/DoubleElectron_Run2012A', 10, 'DATA electron 2012', selection),
        #b.JobConfig('DoubleEle_Run2012B', dCache+'/naodell/nuTuples_v5_8TeV/DoubleElectron_Run2012B', 30, 'DATA electron 2012', selection),
        ])



bg = []
bg.extend([
b.JobConfig('WZJetsTo3LNu',  dCache+'/andreypz/nuTuples_v5_8TeV/WZJetsTo3LNu',  1, 'WZ '+selection+' '+period, selection),
b.JobConfig('WWJetsTo2L2Nu', dCache+'/andreypz/nuTuples_v5_8TeV/WWJetsTo2L2Nu', 1, 'WW '+selection+' '+period, selection),
b.JobConfig('ZZJetsTo2L2Nu', dCache+'/andreypz/nuTuples_v5_8TeV/ZZJetsTo2L2Nu', 1, 'ZZ '+selection+' '+period, selection),
b.JobConfig('ttbar',         dCache+'/andreypz/nuTuples_v5_8TeV/TTJets',       10, 'ttbar '+selection+' '+period, selection),
b.JobConfig('tW',            dCache+'/andreypz/nuTuples_v5_8TeV/tW',            1, 'tW '+selection+' '+period, selection),
b.JobConfig('tbarW',         dCache+'/andreypz/nuTuples_v5_8TeV/tbarW',         1, 'tbarW '+selection+' '+period, selection),
b.JobConfig('DYjets',        dCache+'/andreypz/nuTuples_v5_8TeV/DYjets',       30, 'DYjets '+selection+' '+period, selection),
b.JobConfig('DYjets10',      dCache+'/naodell/nuTuples_v5_8TeV/DYJetsToLL_M-10To50', 20, 'DYjets10 '+selection+' '+period, selection),
b.JobConfig('vbfZ',          dCache+'/andreypz/nuTuples_v5_8TeV/Zvbf',          1, 'vbfZ '+selection+' '+period, selection),

b.JobConfig('WJetsToLNu',    dCache+'/naodell/nuTuples_v5_8TeV/WJetsToLNu',     5, 'WJetsToLNu '+selection+' '+period, selection),
b.JobConfig('WZJetsTo2L2Q',  dCache+'/naodell/nuTuples_v5_8TeV/WZJetsTo2L2Q',   1, 'WZJetsTo2L2Q '+selection+' '+period, selection),
b.JobConfig('ZZJetsTo2L2Q',  dCache+'/naodell/nuTuples_v5_8TeV/ZZJetsTo2L2Q',   1, 'ZZJetsTo2L2Q '+selection+' '+period, selection),

])


signal = []

signal.extend([
b.JobConfig('ggHZZ125', dCache+'/andreypz/nuTuples_v5_8TeV/ggH130', 1, 'ggHZZ125 '+selection+' '+period, selection),
b.JobConfig('ggHZZ200', dCache+'/andreypz/nuTuples_v5_8TeV/ggH200', 1, 'ggHZZ200 '+selection+' '+period, selection),
#b.JobConfig('ggHZZ250', dCache+'/andreypz/nuTuples_v5_8TeV/ggH250', 1, 'ggHZZ250 '+selection+' '+period, selection),
#b.JobConfig('ggHZZ300', dCache+'/andreypz/nuTuples_v5_8TeV/ggH300', 1, 'ggHZZ300 '+selection+' '+period, selection),
#b.JobConfig('ggHZZ350', dCache+'/andreypz/nuTuples_v5_8TeV/ggH350', 1, 'ggHZZ350 '+selection+' '+period, selection),
#b.JobConfig('ggHZZ400', dCache+'/andreypz/nuTuples_v5_8TeV/ggH400', 1, 'ggHZZ400 '+selection+' '+period, selection),

#b.JobConfig('ggHWW125', dCache+'/andreypz/nuTuples_v5_8TeV/ggHWW130', 1, 'ggHWW125 '+selection+' '+period, selection),
#b.JobConfig('ggHWW200', dCache+'/andreypz/nuTuples_v5_8TeV/ggHWW200', 1, 'ggHWW200 '+selection+' '+period, selection),
#b.JobConfig('ggHWW250', dCache+'/andreypz/nuTuples_v5_8TeV/ggHWW250', 1, 'ggHWW250 '+selection+' '+period, selection),
#b.JobConfig('ggHWW300', dCache+'/andreypz/nuTuples_v5_8TeV/ggHWW300', 1, 'ggHWW300 '+selection+' '+period, selection),
#b.JobConfig('ggHWW350', dCache+'/andreypz/nuTuples_v5_8TeV/ggHWW350', 1, 'ggHWW350 '+selection+' '+period, selection),
#b.JobConfig('ggHWW400', dCache+'/andreypz/nuTuples_v5_8TeV/ggHWW400', 1, 'ggHWW400 '+selection+' '+period, selection),

b.JobConfig('VBFHZZ125', dCache+'/andreypz/nuTuples_v5_8TeV/VBFHZZ130', 1, 'VBFHZZ125 '+selection+' '+period, selection),
b.JobConfig('VBFHZZ200', dCache+'/andreypz/nuTuples_v5_8TeV/VBFHZZ200', 1, 'VBFHZZ200 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ250', dCache+'/andreypz/nuTuples_v5_8TeV/VBFHZZ250', 1, 'VBFHZZ250 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ300', dCache+'/andreypz/nuTuples_v5_8TeV/VBFHZZ300', 1, 'VBFHZZ300 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ350', dCache+'/andreypz/nuTuples_v5_8TeV/VBFHZZ350', 1, 'VBFHZZ350 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ400', dCache+'/andreypz/nuTuples_v5_8TeV/VBFHZZ400', 1, 'VBFHZZ400 '+selection+' '+period, selection),

])

# This is just to creare a working(code) directory 
btemp = b.BatchMaster([],outputPath)
btemp.CreateWorkingDir("code_dir")

if doData:
    batcher = b.BatchMaster(data, outputPath, True)
    batcher.SubmitToLPC()
if doBG:
    batcher = b.BatchMaster(bg, outputPath)
    batcher.SubmitToLPC()
if doSignal:
    batcher = b.BatchMaster(signal, outputPath, True)
    batcher.SubmitToLPC()
if doTest:
    batcher = b.BatchMaster(test, outputPath, True)
    batcher.SubmitToLPC()
