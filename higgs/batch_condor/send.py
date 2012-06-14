#! /usr/bin/env python
import BatchMaster as b

''' Specify parameters '''
dCache      = '/pnfs/cms/WAX/11/store/user' 
outputPath  = '/uscms_data/d2/andreypz/cmssw/higgs2/CMSSW_4_4_4/src/NWU/Higgs/higgs/batch_condor'

selection = 'muon'
period    = '2011'
doData    = False
doBG      = False
doSignal  = True

''' 
    Set job configurations.  The order of arguments is:
    (Dataset, path to data, number of jobs, arguments to pass to executable, output directory name)
'''

if selection == 'muon':
    data = []
    if period in ['2011A', 'Combined']:
        data.append(b.JobConfig('Run2011A', dCache+'/naodell/April15/DoubleMu_Run2011A', 40, 'DATA 4,5,8  muon 2011A', selection))

    if period in ['2011B', 'Combined']:
        data.append(b.JobConfig('Run2011B', dCache+'/naodell/April15/DoubleMu_Run2011B', 40, 'DATA 4,5,8 muon 2011B', selection))

                            


if selection == 'electron':
    data = []
    if period in ['2011A', 'Combined']:
        data.append(b.JobConfig('DoubleEle_Run2011A', dCache+'/naodell/April15/DoubleElectron_Run2011A', 30, 'DATA 16,17,18 electron 2011B', selection))

    if period in ['2011B', 'Combined']:
        data.append(b.JobConfig('DoubleEle_Run2011B', dCache+'/naodell/April15/DoubleElectron_Run2011B', 30, 'DATA 16,17,18 electron 2011B', selection))



if selection == 'gamma' or selection == 'muGamma' or selection == 'eGamma':
    data = []
    if period in ['2011A', 'Combined']:
        data.extend([
        b.JobConfig('Photon_Aug05', dCache+'/radek/Data3/Photon_Run2011A_05Aug2011-v1', 20, 'DATA  14,15,16,17,18,19,20,21,22,23 '+selection+' 2011A', selection),
        b.JobConfig('Photon_May10', dCache+'/radek/Data3/Photon_Run2011A_May10ReReco-v1', 10, 'DATA  14,15,16,17,18,19,20,21,22,23 '+selection+' 2011A', selection),
        b.JobConfig('Photon_PromptV4', dCache+'/radek/Data3/Photon_Run2011A_PromptReco-v4-r2', 20, 'DATA  14,15,16,17,18,19,20,21,22,23 '+selection+' 2011A', selection),
        b.JobConfig('Photon_PromptV6', dCache+'/radek/Data3/Photon_Run2011A_PromptReco-v6-r2', 10, 'DATA  14,15,16,17,18,19,20,21,22,23 '+selection+' 2011A', selection),
        ])
    if period in ['2011B', 'Combined']:
        data.append(b.JobConfig('Photon_Run2011B', dCache+'/naodell/Nov09/Photon_Run2011B', 30, 'DATA 15,16,17,18,19,20,21,22,23,24,25 '+selection+' 2011B', selection))

if selection == 'muEG':
    data = []
    if period in ['2011A', 'Combined']:
        data.extend([
        b.JobConfig('Aug05', dCache+'/andreypz/Oct21/MuEG_Aug05', 10, 'DATA  24,25,27,28 muEG 2011A', selection),
        b.JobConfig('May10', dCache+'/andreypz/Oct21/MuEG_May10', 10, 'DATA  24,25,27,28 muEG 2011A', selection),
        b.JobConfig('PromptV4', dCache+'/andreypz/Oct21/MuEG_PromptV4', 20, 'DATA  24,25,27,28 muEG 2011A', selection),
        b.JobConfig('PromptV6', dCache+'/andreypz/Oct21/MuEG_PromptV6', 20, 'DATA  24,25,27,28 muEG 2011A', selection)
        ])
    if period in ['2011B', 'Combined']:
        data.append(b.JobConfig('PromptV6', dCache+'/naodell/Nov09/MuEG_Run2011B', 20, 'DATA  24,25,27,28 muEG 2011B', selection))

bg = []
bg.extend([
#b.JobConfig('ZZ', dCache+'/bpollack/Oct09/MC/ZZJetsTo2L2Nu', 1, 'ZZ 0 '+selection+' '+period, selection),
#b.JobConfig('ggWW', dCache+'/bpollack/Oct09/MC/GluGluToWWTo4L', 1, 'ggWW 0 '+selection+' '+period, selection),
#b.JobConfig('WW', dCache+'/bpollack/Oct09/MC/WWJetsTo2L2Nu', 1, 'WW 0 '+selection+' '+period, selection),
#b.JobConfig('WZ', dCache+'/bpollack/Oct09/MC/WZJetsTo3LNu', 1, 'WZ 0 '+selection+' '+period, selection),
#b.JobConfig('ttbar', dCache+'/radek/MC3b/TTTo2L2Nu2B_7TeV-powheg-pythia6', 20, 'ttbar 0 '+selection+' '+period, selection),
#b.JobConfig('tW', dCache+'/radek/MC3b/T_tW-channel-DR-r2', 1, 'tW 0 '+selection+' '+period, selection),
#b.JobConfig('tbarW', dCache+'/radek/MC3b/Tbar_tW-channel-DR-r2', 1, 'tW 0 '+selection+' '+period, selection),
b.JobConfig('DYjets', dCache+'/bpollack/Apr17/MC/ZJetsV2/', 10, 'DYjets 0 '+selection+' '+period, selection)
#b.JobConfig('ZZ', dCache+'/bpollack/MC/ZZTo2L2Nu', 5, 'ZZ 0 '+selection+' '+period, selection)
#b.JobConfig('WW', dCache+'/bpollack/MC/WWTo2L2Nu', 4, 'WW 0 '+selection+' '+period, selection),
#b.JobConfig('WZ', dCache+'/bpollack/MC/WZTo3LNu', 4, 'WZ 0 '+selection+' '+period, selection),
])

dCache = "/uscms/home/bpollack/nobackup"
signal = []

signal.extend([
b.JobConfig('ggHZZ250', dCache+'/Oct09/MC/ggHZZ_M250', 1, 'ggHZZ250 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ300', dCache+'/Oct09/MC/ggHZZ_M300', 1, 'ggHZZ300 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ350', dCache+'/Oct09/MC/ggHZZ_M350', 1, 'ggHZZ350 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ400', dCache+'/Oct09/MC/ggHZZ_M400', 1, 'ggHZZ400 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ450', dCache+'/Oct09/MC/ggHZZ_M450', 1, 'ggHZZ450 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ500', dCache+'/Oct09/MC/ggHZZ_M500', 1, 'ggHZZ500 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ550', dCache+'/Oct09/MC/ggHZZ_M550', 1, 'ggHZZ550 0 '+selection+' '+period, selection),
b.JobConfig('ggHZZ600', dCache+'/Oct09/MC/ggHZZ_M600', 1, 'ggHZZ600 0 '+selection+' '+period, selection)
#b.JobConfig('ggHWW250', dCache+'/Oct09/MC/ggHWW_M250', 1, 'ggHWW250 0 '+selection+' '+period, selection),
#b.JobConfig('ggHWW300', dCache+'/Oct09/MC/ggHWW_M300', 1, 'ggHWW300 0 '+selection+' '+period, selection),
#b.JobConfig('ggHWW350', dCache+'/Oct09/MC/ggHWW_M350', 1, 'ggHWW350 0 '+selection+' '+period, selection),
#b.JobConfig('ggHWW400', dCache+'/Oct09/MC/ggHWW_M400', 1, 'ggHWW400 0 '+selection+' '+period, selection),
#b.JobConfig('ggHWW450', dCache+'/Oct09/MC/ggHWW_M450', 1, 'ggHWW450 0 '+selection+' '+period, selection),
#b.JobConfig('ggHWW500', dCache+'/Oct09/MC/ggHWW_M500', 1, 'ggHWW500 0 '+selection+' '+period, selection),
#b.JobConfig('ggHWW550', dCache+'/Oct09/MC/ggHWW_M550', 1, 'ggHWW550 0 '+selection+' '+period, selection),
#b.JobConfig('ggHWW600', dCache+'/Oct09/MC/ggHWW_M600', 1, 'ggHWW600 0 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ250', dCache+'/Oct09/MC/VBFHZZ_M250', 1, 'VBFHZZ250 0 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ300', dCache+'/Oct09/MC/VBFHZZ_M300', 1, 'VBFHZZ300 0 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ350', dCache+'/Oct09/MC/VBFHZZ_M350', 1, 'VBFHZZ350 0 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ400', dCache+'/Oct09/MC/VBFHZZ_M400', 1, 'VBFHZZ400 0 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ450', dCache+'/Oct09/MC/VBFHZZ_M450', 1, 'VBFHZZ450 0 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ500', dCache+'/Oct09/MC/VBFHZZ_M500', 1, 'VBFHZZ500 0 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ550', dCache+'/Oct09/MC/VBFHZZ_M550', 1, 'VBFHZZ550 0 '+selection+' '+period, selection),
#b.JobConfig('VBFHZZ600', dCache+'/Oct09/MC/VBFHZZ_M600', 1, 'VBFHZZ600 0 '+selection+' '+period, selection)
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
