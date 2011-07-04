import FWCore.ParameterSet.Config as cms

process = cms.Process("HiggsZZllnunu")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1600) )
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(500)

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag  = 'START41_V0::All'
process.GlobalTag.globaltag  = 'GR_R_41_V0::All'
#process.GlobalTag.globaltag  = 'FT_R_42_V13A::All'

# Hcal Noise filter
#process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minRatio = cms.double(-999)
process.HBHENoiseFilterResultProducer.maxRatio = cms.double(999)
process.HBHENoiseFilterResultProducer.minHPDHits = cms.int32(17)
process.HBHENoiseFilterResultProducer.minRBXHits = cms.int32(999)
process.HBHENoiseFilterResultProducer.minHPDNoOtherHits = cms.int32(10)
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(10)
process.HBHENoiseFilterResultProducer.minHighEHitTime = cms.double(-9999.0)
process.HBHENoiseFilterResultProducer.maxHighEHitTime = cms.double(9999.0)
process.HBHENoiseFilterResultProducer.maxRBXEMF = cms.double(-999.0)
process.HBHENoiseFilterResultProducer.minNumIsolatedNoiseChannels = cms.int32(9999)
process.HBHENoiseFilterResultProducer.minIsolatedNoiseSumE = cms.double(9999)
process.HBHENoiseFilterResultProducer.minIsolatedNoiseSumEt = cms.double(9999)
process.HBHENoiseFilterResultProducer.useTS4TS5 = cms.bool(True)



#Ecal noise filter
###############################################################################
process.load('PhysicsTools/EcalAnomalousEventFilter/ecalanomalouseventfilter_cfi')
process.EcalAnomalousEventFilter.FilterAlgo= cms.untracked.string("TuningMode")
#process.EcalAnomalousEventFilter.FilterAlgo= cms.untracked.string("FilterMode")
process.EcalAnomalousEventFilter.cutBoundEnergyDeadCellsEB=cms.untracked.double(10)
process.EcalAnomalousEventFilter.cutBoundEnergyDeadCellsEE=cms.untracked.double(10)
process.EcalAnomalousEventFilter.cutBoundEnergyGapEB=cms.untracked.double(100)
process.EcalAnomalousEventFilter.cutBoundEnergyGapEE=cms.untracked.double(100)
#process.EcalAnomalousEventFilter.limitFilterToEB=cms.untracked.bool(True)
#process.EcalAnomalousEventFilter.limitFilterToEE=cms.untracked.bool(True)
process.EcalAnomalousEventFilter.enableGap=cms.untracked.bool(False)

process.BE1214 = process.EcalAnomalousEventFilter.clone()
process.BE1214.limitDeadCellToChannelStatusEB = cms.vint32(12,14)
process.BE1214.limitDeadCellToChannelStatusEE = cms.vint32(12,14)
###############################################################################
process.load('JetMETAnalysis.ecalDeadCellTools.RA2TPfilter_cff')

ecalDead = cms.Sequence(process.BE1214)


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#'file:/uscms/home/andreypz/nobackup/mc_minbias.root'
#    '/store/mc/Fall10/QCD_Pt_80to120_TuneZ2_7TeV_pythia6/GEN-SIM-RECO/START38_V12-v1/0000/FE6C9B24-41CB-DF11-9325-00E08178C0CD.root'
#'file:/uscms/home/andreypz/nobackup/mc_minbias_summer11.root'

#'/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/819/FC9C0351-8753-E011-885F-003048F024DE.root',
#'/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/819/CA3D9364-5C53-E011-AAA3-0030487CD7E0.root',
#'/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/819/6AF10E63-5C53-E011-9512-003048F024DE.root',
#'/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/815/C0B6309B-EB53-E011-A13B-000423D94E70.root',

#'/store/data/Run2011A/DoubleMu/AOD/PromptReco-v1/000/160/819/48968165-8753-E011-A317-003048F118C2.root'
#'/store/relval/CMSSW_4_1_4/RelValTTbar_Tauola/GEN-SIM-RECO/START311_V2_PU_E7TeV_AVE_2_BX156-v1/0025/EE9DC9BB-3162-E011-8154-0018F3D09648.root',
#'/store/relval/CMSSW_4_1_4/RelValTTbar_Tauola/GEN-SIM-RECO/START311_V2_PU_E7TeV_AVE_2_BX156-v1/0025/947911E1-3362-E011-A8D4-00304867915A.root',
#'/store/relval/CMSSW_4_1_4/RelValTTbar_Tauola/GEN-SIM-RECO/START311_V2_PU_E7TeV_AVE_2_BX156-v1/0025/5A48C829-3162-E011-8467-002618943865.root',

#'/store/data/Run2011A/DoubleElectron/RECO/May10ReReco-v1/0005/FA14BEF5-4C7C-E011-98D7-001A928116B8.root',
#'/store/data/Run2011A/DoubleElectron/RECO/May10ReReco-v1/0005/F8CF75C9-8B7B-E011-89B7-00304867902E.root',
#'/store/data/Run2011A/DoubleElectron/RECO/May10ReReco-v1/0005/F8BAA7A5-887B-E011-8E51-002618943957.root',
#'/store/data/Run2011A/DoubleElectron/RECO/May10ReReco-v1/0005/F859FB7E-8A7B-E011-ACC9-002618943836.root',

'/store/data/Run2011A/DoubleMu/RECO/PromptReco-v4/000/165/121/1A873C93-A381-E011-902F-0030487CD710.root',

#'/store/data/Run2011A/DoubleElectron/AOD/May10ReReco-v1/0005/F0B77A28-897B-E011-86FC-0026189437FC.root',
#'/store/data/Run2011A/DoubleElectron/AOD/May10ReReco-v1/0005/EE63D495-8B7B-E011-9E52-001A92811700.root',
#'/store/data/Run2011A/DoubleElectron/AOD/May10ReReco-v1/0005/DEAFCEF8-887B-E011-9BF2-00261894397F.root',

#'/store/mc/Spring11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0005/08A906D0-1D55-E011-BEDC-003048679214.root',
#'/store/mc/Spring11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0005/089269D6-1B55-E011-9F14-00261894385D.root',
#'/store/mc/Spring11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0005/08017A2D-1855-E011-B844-00304867D838.root',
#'/store/mc/Spring11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0005/044EA4DB-1E55-E011-9DF7-00261894385D.root',

)
)

# Select primary vertices
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi")
process.offlinePrimaryVerticesDAWithBS = process.offlinePrimaryVerticesDA.clone()
process.offlinePrimaryVerticesDAWithBS.useBeamConstraint = cms.bool(True)
process.offlinePrimaryVerticesDAWithBS.TkClusParameters.TkDAClusParameters.Tmin= cms.double(4.)
process.offlinePrimaryVerticesDAWithBS.TkClusParameters.TkDAClusParameters.vertexSize= cms.double(0.01)

#For eleID maps:
process.load("ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi")
process.simpleEleId80relIso = process.simpleCutBasedElectronID.clone()
process.simpleEleId80relIso.electronQuality = "80relIso"

#from RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_cfi import *
#from RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi import *



process.demo0 = cms.EDAnalyzer('Dummy',)
process.demo1 = process.demo0.clone()
process.demo2 = process.demo0.clone()
process.demo3 = process.demo0.clone()

process.load("RecoJets.Configuration.RecoPFJets_cff")
process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)
process.kt6PFJets.rParam = cms.double(0.6)

process.ak5PFJetsHZZ2l2nu = process.ak5PFJets.clone()
process.ak5PFJetsHZZ2l2nu.doAreaFastjet = True
process.ak5PFJetsHZZ2l2nu.Rho_EtaMax = cms.double(2.5)


process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.correctedJets = process.ak5PFJetsL1FastL2L3.clone()
process.correctedJets.src = cms.InputTag('ak5PFJetsHZZ2l2nu')

### To get b-tags from ak5PFJets
#process.load("RecoBTag.Configuration.RecoBTag_cff")

process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")
process.ak5JetTracksAssociatorAtVertex.jets = cms.InputTag("ak5PFJetsHZZ2l2nu")
process.ak5JetTracksAssociatorAtCaloFace.jets = cms.InputTag("ak5PFJetsHZZ2l2nu")
process.ak5JetExtender.jets = cms.InputTag("ak5PFJetsHZZ2l2nu")
#process.ak5JetTracksAssociatorAtVertex.jets = cms.InputTag("correctedJets")
#process.ak5JetTracksAssociatorAtCaloFace.jets = cms.InputTag("correctedJets")
#process.ak5JetExtender.jets = cms.InputTag("correctedJets")

from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5PFJet
process.metJESCorPFAK5 = metJESCorAK5PFJet.clone()
process.metJESCorPFAK5.inputUncorJetsLabel = "ak5PFJets"
process.metJESCorPFAK5.metType = "PFMET"
process.metJESCorPFAK5.inputUncorMetLabel = "pfMet"
process.metJESCorPFAK5.useTypeII = False
process.metJESCorPFAK5.jetPTthreshold = cms.double(10.0)
process.metJESCorPFAK5.corrector = cms.string('ak5PFL2L3')

process.higgs = cms.EDAnalyzer('HiggsZZllnunuNtuplyzer',
                               RecoJetTag        =    cms.untracked.InputTag("ak5PFJetsHZZ2l2nu"),
#                               RecoJetTag        =    cms.untracked.InputTag("ak5PFJets"),
                               RecoMETTag        =    cms.untracked.InputTag("pfMet"),
                               ElectronTag       =    cms.untracked.InputTag("gsfElectrons"),
                               MuonTag           =    cms.untracked.InputTag("muons"),
                               electronIDMap     =    cms.untracked.InputTag("simpleEleId80relIso"),                               
#                               electronIDMap     =    cms.untracked.InputTag("eidTight"),                               
	                       rhoCorrTag	 =    cms.untracked.InputTag("kt6PFJets", "rho", "HiggsZZllnunu"),
                               hcalFilterTag     = cms.untracked.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"),  
                               ecalAnomalousFilterTag = cms.untracked.InputTag("BE1214","anomalousECALVariables")	
)

process.load("NWU/Higgs/hzzSkim_cff")
#process.load("Configuration/Skimming/PDWG_HZZSkim_cff")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.TimerService = cms.Service("TimerService", useCPUtime = cms.untracked.bool(True))
process.pts = cms.EDFilter("PathTimerInserter")
process.PathTimerService = cms.Service("PathTimerService")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('out_higgs.root')
)

process.hltHighLevel = cms.EDFilter("HLTHighLevel",
                                    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                    HLTPaths = cms.vstring(
    #'HLT_I*','HLT_Mu15*','HLT_Mu2*','HLT_Mu30_*','HLT_Mu*',  # single muons
    #'HLT_DoubleMu4*','HLT_DoubleMu6*','HLT_DoubleMu7*','HLT_L2DoubleMu3*','HLT_L2DoubleMu23*','HLT_T*', # double muons
    #'HLT_D*','HLT_Ele17_*_Ele*','HLT_Ele*_SC*','HLT_T*','HLT_P*'  # double electron
    '*' #Any, will sort this out later (in analyzer) 
    ),
                                    
                                    eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
                                    andOr = cms.bool(True),             # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
                                    throw = cms.bool(False)    # throw exception on unknown path names
                                    )


preSequence = cms.Sequence(process.offlinePrimaryVerticesDAWithBS 
                           * process.simpleEleId80relIso 
                           #  * process.correctedJets
                           * process.kt6PFJets 
                           * process.ak5PFJetsHZZ2l2nu 
                           * process.metJESCorPFAK5 
                           * process.ak5JetTracksAssociatorAtVertex 
                           * process.impactParameterTagInfos 
                           * ( process.trackCountingHighEffBJetTags +
                               process.trackCountingHighPurBJetTags +
                               process.jetProbabilityBJetTags       +
                               process.jetBProbabilityBJetTags) 
                           * process.HBHENoiseFilterResultProducer
                           * ecalDead
                           )

process.dummy0 = cms.Path(process.demo0)
process.diSequence = cms.Path((process.goodHzzMuons + process.goodHzzElectrons)
                      *(process.diHzzMuons + process.diHzzElectrons + process.crossHzzLeptons ) )
                
process.diMuonFilter     = cms.Path(process.diHzzMuonsFilter      *process.demo1*preSequence *process.higgs)
process.diElectronFilter = cms.Path(process.diHzzElectronsFilter  *process.demo2*preSequence *process.higgs)
process.EleMuFilter      = cms.Path(process.crossHzzLeptonsFilter *process.demo3*preSequence *process.higgs)


#process.p = cms.Path(process.hltHighLevel * preSequence * process.higgs) # No Ecal (Run on AOD)
#(process.zzdiMuonSequence + process.zzdiElectronSequence + process.zzeleMuSequence)

###############################################################################
#process.out = cms.OutputModule("PoolOutputModule",
#process.out = cms.OutputModule("AsciiOutputModule",
#                               SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('diMuonFilter','diElectronFilter')
#                                                                 #fileName = cms.untracked.string('aatest.root'),
#                                                                 #                               outputCommands = cms.untracked.vstring('drop *', )
#                                                                 )
#                               )
#process.e = cms.EndPath(process.out * process.higgs)
#############################################
