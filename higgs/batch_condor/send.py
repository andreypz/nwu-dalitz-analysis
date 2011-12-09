#! /usr/bin/env python
import BatchMaster as b

''' Specify parameters '''
dCache      = '/pnfs/cms/WAX/11/store/user'
outputPath  = '/uscms_data/d2/andreypz/cmssw/higgs7/CMSSW_4_2_8/src/NWU/Higgs/higgs/batch_condor'

selection = 'muGamma'
period    = 'Combined'
doData    = True
doBG      = False
doSignal  = False   

''' 
    Set job configurations.  The order of arguments is:
    (Dataset, path to data, number of jobs, arguments to pass to executable, output directory name)
'''

if selection == 'muon':
    data = []
    if period in ['2011A', 'Combined']:
        data.extend([
        b.JobConfig('Aug05', dCache+'/andreypz/Data5/DoubleMu_Aug05', 4, 'DATA 4 muon 2011A', selection) 
#        b.JobConfig('May10', dCache+'/andreypz/Data5/DoubleMu_May10', 4, 'DATA 4,8 muon 2011A', selection), 
#        b.JobConfig('PromptV4', dCache+'/andreypz/Data5/DoubleMu_PromptV4', 20, 'DATA 4,8 muon 2011A', selection), 
#        b.JobConfig('PromptV6', dCache+'/andreypz/Data5/DoubleMu_PromptV6', 20, 'DATA 4,5 muon 2011A', selection)
        ])
#    if period in ['2011B', 'Combined']:
#        data.append(b.JobConfig('Run2011B', dCache+'/radek/Data3/DoubleMu_Run2011B_HZZ_PromptSkim_v1_r2', 30, 'DATA 4,5 muon 2011B', selection))

if selection == 'electron':
    data = []
    if period in ['2011A', 'Combined']:
        data.extend([
        b.JobConfig('Aug05', dCache+'/naodell/Oct21/DoubleElectron_Aug05', 4, 'DATA 12,13,14 electron 2011A', selection),
        b.JobConfig('May10', dCache+'/andreypz/Data5/DoubleElectron_May10', 4, 'DATA 12,13 electron 2011A', selection),
        b.JobConfig('PromptV4', dCache+'/naodell/Oct21/DoubleElectron_PromptV4', 20, 'DATA 12,13,14 electron 2011A', selection),
        b.JobConfig('PromptV6', dCache+'/naodell/Oct21/DoubleElectron_PromptV6', 20, 'DATA 12,13,14 electron 2011A', selection)
        ])
    if period in ['2011B', 'Combined']:
        data.append(b.JobConfig('Run2011B', dCache+'/radek/Data3/DoubleElectron_Run2011B_PromptReco_v1', 30, 'DATA 14 electron 2011B', selection))

if selection == 'gamma' or selection == 'muGamma' or selection == 'eGamma':
    data = []
    if period in ['2011A', 'Combined']:
        data.extend([
#        b.JobConfig('Aug05', dCache+'/radek/Data3/Photon_Run2011A_05Aug2011-v1', 20, 'PhotonJets  14,15,16,17,18,19,20,21,22,23 '+selection+' 2011A', selection)
#        b.JobConfig('May10', dCache+'/radek/Data3/Photon_Run2011A_May10ReReco-v1', 10, 'PhotonJets  14,15,16,17,18,19,20,21,22,23 '+selection+' 2011A', selection),
 #       b.JobConfig('PromptV4', dCache+'/radek/Data3/Photon_Run2011A_PromptReco-v4-r2', 20, 'PhotonJets  14,15,16,17,18,19,20,21,22,23 '+selection+' 2011A', selection),
  #      b.JobConfig('PromptV6', dCache+'/radek/Data3/Photon_Run2011A_PromptReco-v6-r2', 10, 'PhotonJets  14,15,16,17,18,19,20,21,22,23 '+selection+' 2011A', selection),
        ])
    if period in ['2011B', 'Combined']:
        data.append(b.JobConfig('Run2011B', dCache+'/naodell/Nov09/Photon_Run2011B', 30, 'PhotonJets 15,16,17,18,19,20,21,22,23,24,25 '+selection+' 2011B', selection))

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
#b.JobConfig('ZZ', dCache+'/bpollack/MC/ZZTo2L2Nu', 5, 'ZZ 0 '+selection+' '+period, selection)
#b.JobConfig('WW', dCache+'/bpollack/MC/WWTo2L2Nu', 4, 'WW 0 '+selection+' '+period, selection),
#b.JobConfig('WZ', dCache+'/bpollack/MC/WZTo3LNu', 4, 'WZ 0 '+selection+' '+period, selection),
b.JobConfig('ZZ', dCache+'/bpollack/Oct09/MC/ZZJetsTo2L2Nu', 5, 'ZZ 0 '+selection+' '+period, selection)
#b.JobConfig('GluGluWW', dCache+'/bpollack/Oct09/MC/GluGluToWWTo4L', 4, 'GluGluWW 0 '+selection+' '+period, selection),
##b.JobConfig('WWJets', dCache+'/bpollack/Oct09/MC/WWJetsTo2L2Nu', 4, 'WWJets 0 '+selection+' '+period, selection),
#b.JobConfig('WZJets', dCache+'/bpollack/Oct09/MC/WZJetsTo3LNu', 4, 'WZJets 0 '+selection+' '+period, selection),
#b.JobConfig('ttbar', dCache+'/radek/MC3b/TTTo2L2Nu2B_7TeV-powheg-pythia6', 20, 'ttbar 0 '+selection+' '+period, selection),
#b.JobConfig('tW', dCache+'/radek/MC3b/T_tW-channel-DR-r2', 3, 'tW 0 '+selection+' '+period, selection),
#b.JobConfig('tbarW', dCache+'/radek/MC3b/Tbar_tW-channel-DR-r2', 3, 'tW 0 '+selection+' '+period, selection),
#b.JobConfig('ZJets', dCache+'/naodell/Oct09/MC/DYJetsToLL', 30, 'ZJets 0 '+selection+' '+period, selection)
])

signal = []
signal.extend([
b.JobConfig('HZZ250', dCache+'/bpollack/Oct09/MC/ggHZZ_M250', 2, 'HZZ250 0 '+selection+' '+period, selection),
b.JobConfig('HZZ300', dCache+'/bpollack/Oct09/MC/ggHZZ_M300', 2, 'HZZ300 0 '+selection+' '+period, selection),
b.JobConfig('HZZ350', dCache+'/bpollack/Oct09/MC/ggHZZ_M350', 2, 'HZZ350 0 '+selection+' '+period, selection),
b.JobConfig('HZZ400', dCache+'/bpollack/Oct09/MC/ggHZZ_M400', 2, 'HZZ400 0 '+selection+' '+period, selection),
b.JobConfig('HZZ450', dCache+'/bpollack/Oct09/MC/ggHZZ_M450', 2, 'HZZ450 0 '+selection+' '+period, selection),
b.JobConfig('HZZ500', dCache+'/bpollack/Oct09/MC/ggHZZ_M500', 2, 'HZZ500 0 '+selection+' '+period, selection),
b.JobConfig('HZZ550', dCache+'/bpollack/Oct09/MC/ggHZZ_M550', 2, 'HZZ550 0 '+selection+' '+period, selection),
b.JobConfig('HZZ600', dCache+'/bpollack/Oct09/MC/ggHZZ_M600', 2, 'HZZ600 0 '+selection+' '+period, selection),
b.JobConfig('HWW250', dCache+'/bpollack/Oct09/MC/ggHWW_M250', 2, 'HWW250 0 '+selection+' '+period, selection),
b.JobConfig('HWW300', dCache+'/bpollack/Oct09/MC/ggHWW_M300', 2, 'HWW300 0 '+selection+' '+period, selection),
b.JobConfig('HWW350', dCache+'/bpollack/Oct09/MC/ggHWW_M350', 2, 'HWW350 0 '+selection+' '+period, selection),
b.JobConfig('HWW400', dCache+'/bpollack/Oct09/MC/ggHWW_M400', 2, 'HWW400 0 '+selection+' '+period, selection),
b.JobConfig('HWW450', dCache+'/bpollack/Oct09/MC/ggHWW_M450', 2, 'HWW450 0 '+selection+' '+period, selection),
b.JobConfig('HWW500', dCache+'/bpollack/Oct09/MC/ggHWW_M500', 2, 'HWW500 0 '+selection+' '+period, selection),
b.JobConfig('HWW550', dCache+'/bpollack/Oct09/MC/ggHWW_M550', 2, 'HWW550 0 '+selection+' '+period, selection),
b.JobConfig('HWW600', dCache+'/bpollack/Oct09/MC/ggHWW_M600', 2, 'HWW600 0 '+selection+' '+period, selection),
b.JobConfig('VBF250', dCache+'/bpollack/Oct09/MC/VBFHZZ_M250', 2, 'VBF250 0 '+selection+' '+period, selection),
b.JobConfig('VBF300', dCache+'/bpollack/Oct09/MC/VBFHZZ_M300', 2, 'VBF300 0 '+selection+' '+period, selection),
b.JobConfig('VBF350', dCache+'/bpollack/Oct09/MC/VBFHZZ_M350', 2, 'VBF350 0 '+selection+' '+period, selection),
b.JobConfig('VBF400', dCache+'/bpollack/Oct09/MC/VBFHZZ_M400', 2, 'VBF400 0 '+selection+' '+period, selection),
b.JobConfig('VBF450', dCache+'/bpollack/Oct09/MC/VBFHZZ_M450', 2, 'VBF450 0 '+selection+' '+period, selection),
b.JobConfig('VBF500', dCache+'/bpollack/Oct09/MC/VBFHZZ_M500', 2, 'VBF500 0 '+selection+' '+period, selection),
b.JobConfig('VBF550', dCache+'/bpollack/Oct09/MC/VBFHZZ_M550', 2, 'VBF550 0 '+selection+' '+period, selection),
b.JobConfig('VBF600', dCache+'/bpollack/Oct09/MC/VBFHZZ_M600', 2, 'VBF600 0 '+selection+' '+period, selection)
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
