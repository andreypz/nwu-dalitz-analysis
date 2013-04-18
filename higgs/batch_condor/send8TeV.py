#!/usr/bin/env python

import BatchMaster as b
import sys
cfg = b.JobConfig

nargs = len(sys.argv)
print "Running", sys.argv[0], nargs


''' Specify parameters '''
dCache      = '/pnfs/cms/WAX/11/store/user' 
EOS         = '/eos/uscms/store/user'
outputPath  = EOS+'/andreypz/batch_output/8TeV'
executable  = 'batchJob.csh'

selection = 'muon'
if nargs>1 and sys.argv[1]=="ele":
    selection="electron"
    
version   = "v80"
period    = '2012'
doTest    = 1
doData    = 0
doBG      = 0
doSignal  = 0

''' 
    Set job configurations.  The order of arguments is:
    (Dataset, path to data, number of jobs, arguments to pass to executable, output directory name)
'''

test = []
data = []
test.extend([
#cfg('ZZ', dCache+'/andreypz/nuTuples_v5_8TeV/ZZJetsTo2L2Nu', 10, 'ZZ '+selection+' '+period),
#cfg('DoubleMu_Run2012A', dCache+'/devildog/nuTuples_v1_8TeV/DoubleMu_Run2012A', 1, 'DATA muon 2012'),
#cfg('DoubleMu_Run2012B', dCache+'/devildog/nuTuples_v1_8TeV/DoubleMu_Run2012B', 10, 'DATA muon 2012'),
#cfg('DoubleMu_Run2012A', dCache+'/devildog/nuTuples_v1_8TeV/DoubleMu_Run2012A', 10, 'DATA muon 2012'),
#cfg('DoubleMu_Run2012B', dCache+'/devildog/nuTuples_v1_8TeV/DoubleMu_Run2012B', 10, 'DATA muon 2012'),
#cfg('vbfZ',          dCache+'/andreypz/nuTuples_v5_8TeV/Zvbf',          1, 'vbfZ '+selection+' '+period),
#cfg('DYjets10',      dCache+'/naodell/nuTuples_v5_8TeV/DYJetsToLL_M-10To50', 20, 'DYjets10 '+selection+' '+period),
#cfg('WZJetsTo3LNu',  dCache+'/andreypz/nuTuples_v5_8TeV/WZJetsTo3LNu',  1, 'WZ '+selection+' '+period),
#cfg('DYjets',        dCache+'/andreypz/nuTuples_v5_8TeV/DYjets',       30, 'DYjets '+selection+' '+period),
    #cfg('DoubleMu_Run2012C_24Aug',  dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012C_24Aug_v2',   5, 'DATA muon 2012'),
    cfg('DoubleEle_Run2012D',        dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012D',         20, 'DATA electron 2012'),
])



if selection == 'muon':
    data.extend([
        cfg('DoubleMu_Run2012A',        dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012A_v2',         5, 'DATA muon 2012'),
        cfg('DoubleMu_Run2012B',        dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012B_v2',        10, 'DATA muon 2012'),
        cfg('DoubleMu_Run2012C_prompt', dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012C_Prompt_v2', 10, 'DATA muon 2012'),
        cfg('DoubleMu_Run2012C_recover',dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012C_Recover_v2', 2, 'DATA muon 2012'),
        cfg('DoubleMu_Run2012C_24Aug',  dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012C_24Aug_v2',   2, 'DATA muon 2012'),
        cfg('DoubleMu_Run2012D',        dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012D_v2',        10, 'DATA muon 2012'),
        #cfg('DoubleMu_Run2012A', dCache+'/devildog/nuTuples_v1_8TeV/DoubleMu_Run2012A', 10, 'DATA muon 2012'),
        #cfg('DoubleMu_Run2012B', dCache+'/devildog/nuTuples_v1_8TeV/DoubleMu_Run2012B', 30, 'DATA muon 2012'),
        ])

if selection == 'electron':
    data.extend([
        cfg('DoubleEle_Run2012A',        dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012A_v2',         10, 'DATA electron 2012'),
        cfg('DoubleEle_Run2012B',        dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012B_v2',         20, 'DATA electron 2012'),
        cfg('DoubleEle_Run2012C_prompt', dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012C_Prompt_v2',  20, 'DATA electron 2012'),
        cfg('DoubleEle_Run2012C_recover',dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012C_Recover_v2', 20, 'DATA electron 2012'),
        cfg('DoubleEle_Run2012C_24Aug',  dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012C_24Aug_v2',   20, 'DATA electron 2012'),
        cfg('DoubleEle_Run2012D',        dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012D',         20, 'DATA electron 2012'),

        #cfg('DoubleEle_Run2012A', dCache+'/naodell/nuTuples_v5_8TeV/DoubleElectron_Run2012A', 10, 'DATA electron 2012'),
        #cfg('DoubleEle_Run2012B', dCache+'/naodell/nuTuples_v5_8TeV/DoubleElectron_Run2012B', 30, 'DATA electron 2012'),
        ])



bg = []
bg.extend([
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
])


signal = []

signal.extend([
cfg('ggHZZ125', dCache+'/andreypz/nuTuples_v5_8TeV/ggH130', 1, 'ggHZZ125 '+selection+' '+period),
cfg('ggHZZ200', dCache+'/andreypz/nuTuples_v5_8TeV/ggH200', 1, 'ggHZZ200 '+selection+' '+period),
#cfg('ggHZZ250', dCache+'/andreypz/nuTuples_v5_8TeV/ggH250', 1, 'ggHZZ250 '+selection+' '+period),
#cfg('ggHZZ300', dCache+'/andreypz/nuTuples_v5_8TeV/ggH300', 1, 'ggHZZ300 '+selection+' '+period),
#cfg('ggHZZ350', dCache+'/andreypz/nuTuples_v5_8TeV/ggH350', 1, 'ggHZZ350 '+selection+' '+period),
#cfg('ggHZZ400', dCache+'/andreypz/nuTuples_v5_8TeV/ggH400', 1, 'ggHZZ400 '+selection+' '+period),

#cfg('ggHWW125', dCache+'/andreypz/nuTuples_v5_8TeV/ggHWW130', 1, 'ggHWW125 '+selection+' '+period),
#cfg('ggHWW200', dCache+'/andreypz/nuTuples_v5_8TeV/ggHWW200', 1, 'ggHWW200 '+selection+' '+period),
#cfg('ggHWW250', dCache+'/andreypz/nuTuples_v5_8TeV/ggHWW250', 1, 'ggHWW250 '+selection+' '+period),
#cfg('ggHWW300', dCache+'/andreypz/nuTuples_v5_8TeV/ggHWW300', 1, 'ggHWW300 '+selection+' '+period),
#cfg('ggHWW350', dCache+'/andreypz/nuTuples_v5_8TeV/ggHWW350', 1, 'ggHWW350 '+selection+' '+period),
#cfg('ggHWW400', dCache+'/andreypz/nuTuples_v5_8TeV/ggHWW400', 1, 'ggHWW400 '+selection+' '+period),

cfg('VBFHZZ125', dCache+'/andreypz/nuTuples_v5_8TeV/VBFHZZ130', 1, 'VBFHZZ125 '+selection+' '+period),
cfg('VBFHZZ200', dCache+'/andreypz/nuTuples_v5_8TeV/VBFHZZ200', 1, 'VBFHZZ200 '+selection+' '+period),
#cfg('VBFHZZ250', dCache+'/andreypz/nuTuples_v5_8TeV/VBFHZZ250', 1, 'VBFHZZ250 '+selection+' '+period),
#cfg('VBFHZZ300', dCache+'/andreypz/nuTuples_v5_8TeV/VBFHZZ300', 1, 'VBFHZZ300 '+selection+' '+period),
#cfg('VBFHZZ350', dCache+'/andreypz/nuTuples_v5_8TeV/VBFHZZ350', 1, 'VBFHZZ350 '+selection+' '+period),
#cfg('VBFHZZ400', dCache+'/andreypz/nuTuples_v5_8TeV/VBFHZZ400', 1, 'VBFHZZ400 '+selection+' '+period),

])

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
    batcher = b.BatchMaster(inputSamples, outputPath, shortQueue = False, stageDir = '../../batchStage', executable = executable, selection = version + '/' + selection +"_"+ period)
    print "Submitting to batch?"
    batcher.submit_to_batch()
            
