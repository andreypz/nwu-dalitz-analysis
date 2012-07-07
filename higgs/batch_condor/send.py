#! /usr/bin/env python
import BatchMaster as b

''' Specify parameters '''
dCache      = '/pnfs/cms/WAX/11/store/user' 
outputPath  = '/uscms_data/d2/andreypz/cmssw/higgs2/CMSSW_4_4_4/src/NWU/Higgs/higgs/batch_condor'

selection = 'muon'
period    = '2011'
doData    = False
doBG      = True
doSignal  = False

''' 
    Set job configurations.  The order of arguments is:
    (Dataset, path to data, number of jobs, arguments to pass to executable, output directory name)
'''

if selection == 'muon':
    data = []
    data.append(b.JobConfig('DoubleMu_Run2011A', dCache+'/andreypz/nuTuples_v1_7TeV/DoubleMu_HZZ_Run2011A', 40, 'DATA 4,5,8 muon 2011', selection))
    data.append(b.JobConfig('DoubleMu_Run2011B', dCache+'/andreypz/nuTuples_v1_7TeV/DoubleMu_HZZ_Run2011B', 40, 'DATA 4,5,8 muon 2011', selection))

                            


if selection == 'electron':
    data = []
    data.append(b.JobConfig('DoubleEle_Run2011A', dCache+'/andreypz/nuTuples_v1_7TeV/DoubleElectron_HZZ_Run2011A', 30, 'DATA 16,17,18 electron 2011', selection))
    data.append(b.JobConfig('DoubleEle_Run2011B', dCache+'/andreypz/nuTuples_v1_7TeV/DoubleElectron_HZZ_Run2011B', 30, 'DATA 16,17,18 electron 2011', selection))



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
b.JobConfig('ZZ', dCache+'/andreypz/nuTuples_v1_7TeV/ZZJetsTo2L2Nu', 1, 'ZZ 0 '+selection+' '+period, selection),
#b.JobConfig('WW', dCache+'/andreypz/nuTuples_v1_7TeV/WWTo2L2Nu', 1, 'WW 0 '+selection+' '+period, selection),
#b.JobConfig('WZ', dCache+'/andreypz/nuTuples_v1_7TeV/WZTo3LNu', 1, 'WZ 0 '+selection+' '+period, selection),
#b.JobConfig('ttbar', dCache+'/andreypz/nuTuples_v1_7TeV/TTJets', 20, 'ttbar 0 '+selection+' '+period, selection),
#b.JobConfig('tW', dCache+'/andreypz/nuTuples_v1_7TeV/T_tW-channel-DR-r2', 1, 'tW 0 '+selection+' '+period, selection),
#b.JobConfig('tbarW', dCache+'/andreypz/nuTuples_v1_7TeV/Tbar_tW-channel-DR-r2', 1, 'tW 0 '+selection+' '+period, selection),
#b.JobConfig('DYjets', dCache+'/andreypz/nuTuples_v1_7TeV/MC/DYJetsToLL', 10, 'DYjets 0 '+selection+' '+period, selection)
#b.JobConfig('ggWW', dCache+'/andreypz/nuTuples_v1_7TeV/GluGluToWWTo4L', 1, 'ggWW 0 '+selection+' '+period, selection),
])


signal = []

signal.extend([
b.JobConfig('ggHZZ200', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ200', 1, 'ggHZZ200 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ250', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ250', 1, 'ggHZZ250 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ300', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ300', 1, 'ggHZZ300 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ350', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ350', 1, 'ggHZZ350 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ400', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ400', 1, 'ggHZZ400 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ450', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ450', 1, 'ggHZZ450 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ500', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ500', 1, 'ggHZZ500 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ550', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ550', 1, 'ggHZZ550 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ600', dCache+'/andreypz/nuTuples_v1_7TeV/ggHZZ600', 1, 'ggHZZ600 0 '+selection+' '+period, selection)
#b.JobConfig('ggHWW250', dCache+'/andreypz/nuTuples_v1_7TeV/ggHWW_M250', 1, 'ggHWW250 0 '+selection+' '+period, selection),
#b.JobConfig('ggHWW300', dCache+'/andreypz/nuTuples_v1_7TeV/ggHWW_M300', 1, 'ggHWW300 0 '+selection+' '+period, selection),
#b.JobConfig('ggHWW350', dCache+'/andreypz/nuTuples_v1_7TeV/ggHWW_M350', 1, 'ggHWW350 0 '+selection+' '+period, selection),
#b.JobConfig('ggHWW400', dCache+'/andreypz/nuTuples_v1_7TeV/ggHWW_M400', 1, 'ggHWW400 0 '+selection+' '+period, selection),
#b.JobConfig('ggHWW450', dCache+'/andreypz/nuTuples_v1_7TeV/ggHWW_M450', 1, 'ggHWW450 0 '+selection+' '+period, selection),
#b.JobConfig('ggHWW500', dCache+'/andreypz/nuTuples_v1_7TeV/ggHWW_M500', 1, 'ggHWW500 0 '+selection+' '+period, selection),
#b.JobConfig('ggHWW550', dCache+'/andreypz/nuTuples_v1_7TeV/ggHWW_M550', 1, 'ggHWW550 0 '+selection+' '+period, selection),
#b.JobConfig('ggHWW600', dCache+'/andreypz/nuTuples_v1_7TeV/ggHWW_M600', 1, 'ggHWW600 0 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ250', dCache+'/andreypz/nuTuples_v1_7TeV/VBFHZZ_M250', 1, 'VBFHZZ250 0 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ300', dCache+'/andreypz/nuTuples_v1_7TeV/VBFHZZ_M300', 1, 'VBFHZZ300 0 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ350', dCache+'/andreypz/nuTuples_v1_7TeV/VBFHZZ_M350', 1, 'VBFHZZ350 0 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ400', dCache+'/andreypz/nuTuples_v1_7TeV/VBFHZZ_M400', 1, 'VBFHZZ400 0 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ450', dCache+'/andreypz/nuTuples_v1_7TeV/VBFHZZ_M450', 1, 'VBFHZZ450 0 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ500', dCache+'/andreypz/nuTuples_v1_7TeV/VBFHZZ_M500', 1, 'VBFHZZ500 0 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ550', dCache+'/andreypz/nuTuples_v1_7TeV/VBFHZZ_M550', 1, 'VBFHZZ550 0 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ600', dCache+'/andreypz/nuTuples_v1_7TeV/VBFHZZ_M600', 1, 'VBFHZZ600 0 '+selection+' '+period, selection)
])

if doData:
    batcher = b.BatchMaster(data, outputPath)
    batcher.SubmitToLPC()
if doBG:
    batcher = b.BatchMaster(bg, outputPath)
    batcher.SubmitToLPC()
if doSignal:
    batcher = b.BatchMaster(signal, outputPath)
    batcher.SubmitToLPC()
