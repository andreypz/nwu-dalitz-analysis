// -*- C++ -*-
//
// Package:    HiggsZZllnunuNtuplyzer
// Class:      HiggsZZllnunuNtuplyzer
// 
/**\class HiggsZZllnunuNtuplyzer HiggsZZllnunuNtuplyzer.cc NWU/Higgs/src/HiggsZZllnunuNtuplyzer.cc

 Description: Higgs to ZZ to 2l2nu analysis

 Implementation: Reco to Root ntuplyzer.

*/
//
// Original Author:  Andrey Pozdnyakov
//         Created:  Fri Mar 25 08:07:03 CDT 2011
// $Id: HiggsZZllnunuNtuplyzer.cc,v 1.3 2011/08/14 18:52:42 andrey Exp $
//

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// JEC
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "JetMETCorrections/Type1MET/interface/Type1MET.h"
#include "JetMETCorrections/Type1MET/interface/Type1METAlgo.h"
#include "DataFormats/AnomalousEcalDataFormats/interface/AnomalousECALVariables.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"

// Need for HLT trigger info:
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"


//ROOT stuff
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TVector3.h"
#include "TClonesArray.h"

// ntuple storage classes
#include "TCPrimaryVtx.h"
#include "TCJet.h"
#include "TCMET.h"
#include "TCMuon.h"
#include "TCElectron.h"

//#include "CMGTools/HtoZZ2l2nu/interface/ReducedMETComputer.h"

//typedef struct {Float_t x,y,z;} POINT;

using namespace edm;
using namespace std;
using namespace reco;

class HiggsZZllnunuNtuplyzer : public edm::EDAnalyzer {
   public:
      explicit HiggsZZllnunuNtuplyzer(const edm::ParameterSet&);
      ~HiggsZZllnunuNtuplyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  // ----------member data ---------------------------
  //-- Functions ---
  float sumPtSquared(const Vertex& v);
  bool isFilteredOutScraping(const edm::Event& iEvent, const edm::EventSetup& iSetup, int numtrack=10, double thresh=0.25);
  bool triggerDecision(edm::Handle<edm::TriggerResults> &hltR, int iTrigger);

  //-- Inputs --  
  edm::InputTag recoJetTag_, recoMETTag_;
  edm::InputTag muonTag_, electronTag_, electronIDMap_;
  edm::InputTag rhoCorrTag_, hcalFilterTag_, ecalAnomalousFilterTag_;
  edm::InputTag hltTriggerResults_, hltName_;
  bool useEcalFilter_;

//--- Variables ---
  Bool_t isRealData;
  UInt_t eventNumber, runNumber, lumiSection, bunchCross;
  Float_t rhoFactor;

//--- Objects ---  
  TriggerNames triggerNames;

  //ReducedMETComputer rmet_;
  
  struct {Float_t x,y,z;} beamSpot;

  //TClonesArray *beamSpot;
  // static POINT beamSpot;
  //TClonesArray *primaryVtx;
  TClonesArray *primaryVtxWithBS, *primaryVtxDAWithBS, *primaryVtxDA;
  TClonesArray *recoJets, *recoMET, *corrMET, *recoMuons, *recoElectrons;

  Bool_t isNoiseHcal, isDeadEcalCluster;
  Bool_t isCSCTightHalo, isCSCLooseHalo, isHcalTightHalo, isHcalLooseHalo, isEcalTightHalo, isEcalLooseHalo;
  Bool_t isGlobalTightHalo, isGlobalLooseHalo;
  Bool_t isScraping;

  UInt_t triggerStatusMuon, triggerStatusElectron;
 
  TTree * rTree;
  edm::Service<TFileService> fs;
};


HiggsZZllnunuNtuplyzer::HiggsZZllnunuNtuplyzer(const edm::ParameterSet& iConfig)
{

  //rmet_(1.0,1.0,0.0,0.0,1.0);

  recoJetTag_       = iConfig.getUntrackedParameter<edm::InputTag>("RecoJetTag");
  recoMETTag_       = iConfig.getUntrackedParameter<edm::InputTag>("RecoMETTag");
  muonTag_          = iConfig.getUntrackedParameter<edm::InputTag>("MuonTag");
  electronTag_      = iConfig.getUntrackedParameter<edm::InputTag>("ElectronTag");
  electronIDMap_    = iConfig.getUntrackedParameter<edm::InputTag>("electronIDMap");
  rhoCorrTag_       = iConfig.getUntrackedParameter<edm::InputTag>("rhoCorrTag");
  hcalFilterTag_    = iConfig.getUntrackedParameter<edm::InputTag>("hcalFilterTag");
  useEcalFilter_ = iConfig.getUntrackedParameter<bool>("UseEcalFilter","false");
  ecalAnomalousFilterTag_ = iConfig.getUntrackedParameter<edm::InputTag>("ecalAnomalousFilterTag");
  
  //hltTriggerResults_ =  iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults","TriggerResults");
  //hltName_          = iConfig.getUntrackedParameter<edm::InputTag>("hltName","HLT");
  //triggerResultsTag_
}


HiggsZZllnunuNtuplyzer::~HiggsZZllnunuNtuplyzer()
{
}

void HiggsZZllnunuNtuplyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std; 
   using namespace reco;

   eventNumber  = iEvent.id().event(); 
   runNumber    = iEvent.id().run();
   lumiSection  = (unsigned int)iEvent.luminosityBlock();
   bunchCross   = (unsigned int)iEvent.bunchCrossing();
   isRealData   = iEvent.isRealData();
   if (!isRealData) bunchCross = 111; //just because in MC it's not filled 

   edm::Handle<reco::BeamSpot> beamSpotHandle;
   iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
   reco::BeamSpot vertexBeamSpot = *beamSpotHandle;

   //TVector3* bs = new ((*beamSpot)[0]) TVector3;
   //bs -> SetXYZ(vertexBeamSpot.position().x(), vertexBeamSpot.position().y(), vertexBeamSpot.position().z());
   beamSpot.x = vertexBeamSpot.position().x();
   beamSpot.y = vertexBeamSpot.position().y();
   beamSpot.z = vertexBeamSpot.position().z();

  //----------------------------------------------------
  //---------------- Filling Primary Verticies -------------
  //--------------------------------------------------------

  //Handle<reco::VertexCollection> primaryVtcs;
  //iEvent.getByLabel("offlinePrimaryVertices", primaryVtcs);
  Handle<reco::VertexCollection> primaryVtcsWithBS;
  iEvent.getByLabel("offlinePrimaryVerticesWithBS", primaryVtcsWithBS);
  Handle<reco::VertexCollection> primaryVtcsDAWithBS;
  iEvent.getByLabel("offlinePrimaryVerticesDAWithBS", primaryVtcsDAWithBS);
  Handle<reco::VertexCollection> primaryVtcsDA;
  iEvent.getByLabel("VerticesDA", primaryVtcsDA);

  
  Int_t  nPrimVtx = 0;
  //Double_t chi2OverNdof = 1000000;
  // Not using this collection for now  
  for(VertexCollection::const_iterator vtx_iter = primaryVtcsDA->begin(); vtx_iter!= primaryVtcsDA->end(); ++vtx_iter)
    {
      reco::Vertex myVtx = reco::Vertex(*vtx_iter);
      TCPrimaryVtx* vtxCon = new ((*primaryVtxDA)[nPrimVtx]) TCPrimaryVtx;
      vtxCon->SetPosition(myVtx.x(), myVtx.y(), myVtx.z());
      vtxCon->SetNDof(myVtx.ndof());
      vtxCon->SetChi2(myVtx.chi2());
      vtxCon->SetIsFake(myVtx.isFake());
      vtxCon->SetNtracks(myVtx.nTracks(0)); //0 is the minWeight, default is 0.5
      if (myVtx.nTracks(0)!=myVtx.tracksSize()) LogWarning(" Vertex nTracks: ")<<"is something wrong here?";
      vtxCon->SetSumPt2Trks(sumPtSquared(myVtx));
      ++nPrimVtx;
    }// end of Primary Vertecies loop
  
  
//PV with Beam Spot constrain

  nPrimVtx=0;
  for(VertexCollection::const_iterator vtx_iter = primaryVtcsWithBS->begin(); vtx_iter!= primaryVtcsWithBS->end(); ++vtx_iter)
    {

      reco::Vertex myVtx = reco::Vertex(*vtx_iter);
      TCPrimaryVtx* vtxCon = new ((*primaryVtxWithBS)[nPrimVtx]) TCPrimaryVtx;
    /*
      vtxCon->SetPosition(myVtx.x(), myVtx.y(), myVtx.z());
      vtxCon->SetNDof(myVtx.ndof());
      vtxCon->SetChi2(myVtx.chi2());
      vtxCon->SetIsFake(myVtx.isFake());
      vtxCon->SetNtracks(myVtx.nTracks(0)); //0 is the minWeight, default is 0.5
      if (myVtx.nTracks(0)!=myVtx.tracksSize()) LogWarning(" Vertex nTracks: ")<<"is something wrong here?";
      vtxCon->SetSumPt2Trks(sumPtSquared(myVtx));
*/

      ++nPrimVtx;
    }// end of Primary Vertecies loop


  //DA PV with Beam Spot constrain
  nPrimVtx=0;
  for(VertexCollection::const_iterator vtx_iter = primaryVtcsDAWithBS->begin(); vtx_iter!= primaryVtcsDAWithBS->end(); ++vtx_iter)
    {
      reco::Vertex myVtx = reco::Vertex(*vtx_iter);
      TCPrimaryVtx* vtxCon = new ((*primaryVtxDAWithBS)[nPrimVtx]) TCPrimaryVtx;
      /*
      vtxCon->SetPosition(myVtx.x(), myVtx.y(), myVtx.z());
      vtxCon->SetNDof(myVtx.ndof());
      vtxCon->SetChi2(myVtx.chi2());
      vtxCon->SetIsFake(myVtx.isFake());
      vtxCon->SetNtracks(myVtx.nTracks(0)); //0 is the minWeight, default is 0.5
      if (myVtx.nTracks(0)!=myVtx.tracksSize()) LogWarning(" Vertex nTracks: ")<<"is something wrong here?";
      vtxCon->SetSumPt2Trks(sumPtSquared(myVtx));
      */

      ++nPrimVtx;
    }// end of Primary Vertecies loop


  //------ End of PV filling ----------------

  ///---------------------------
  //---------  Jets -----------
  //---------------------------

  const JetCorrector* correctorL1L2L3  = JetCorrector::getJetCorrector("ak5PFL1FastL2L3",iSetup);
  //const JetCorrector* correctorL1  = JetCorrector::getJetCorrector("ak5PFL1Offset",iSetup);
  //const JetCorrector* correctorL2  = JetCorrector::getJetCorrector("ak5PFL2Relative",iSetup);
  //const JetCorrector* correctorL3  = JetCorrector::getJetCorrector("ak5PFL3Absolute",iSetup);
  const JetCorrector* correctorRes = JetCorrector::getJetCorrector("ak5PFResidual", iSetup);

  //Getting rho factor
  Handle<double> rhoCorr;
  iEvent.getByLabel(rhoCorrTag_, rhoCorr);
  rhoFactor = (Float_t)(*rhoCorr);

  Handle<reco::PFJetCollection> PFJets;
  iEvent.getByLabel(recoJetTag_, PFJets);

  //These Btag collections are produced from recoJetTag_ of PFJets, therefore they should be exactly the same (pt, eta, etc.)
  edm::Handle<reco::JetTagCollection> bTagHandle1;
  iEvent.getByLabel("trackCountingHighEffBJetTags", bTagHandle1);
  reco::JetTagCollection::const_iterator jet_b_1 = bTagHandle1->begin();

  edm::Handle<reco::JetTagCollection> bTagHandle2;
  iEvent.getByLabel("trackCountingHighPurBJetTags", bTagHandle2);
  reco::JetTagCollection::const_iterator jet_b_2 = bTagHandle2->begin();


  //std::vector<LorentzVector> jetmomenta;

  Int_t cc= 0;
  Int_t pfCount = 0;
  for (PFJetCollection::const_iterator jet_iter = PFJets->begin(); jet_iter!= PFJets->end(); ++jet_iter) 
    {
      reco::PFJet myJet  = reco::PFJet(*jet_iter);
      reco::PFJet corJet = reco::PFJet(*jet_iter);

      int index = jet_iter - PFJets->begin();
      edm::RefToBase<reco::Jet> jetRef(edm::Ref<PFJetCollection>(PFJets,index));
      double scale123 = correctorL1L2L3->correction(*jet_iter,jetRef,iEvent,iSetup);

      corJet.scaleEnergy(scale123);
      
      Int_t nLowPtJets = 0;
      if (corJet.pt() > 20) 
	{
	  // Don't need very low-pt jets
	  //cout<<"dbg  "<<pfCount<<endl;
	  TCJet* jetCon = new ((*recoJets)[pfCount]) TCJet;
  
	  Float_t dR1 = deltaR(myJet, *(jet_b_1->first));
	  Float_t dR2 = deltaR(myJet, *(jet_b_2->first));
       

	  if (dR1>0.001 || dR2>0.001)
	    {
	      //jetCon->SetBDiscrTrkCountHiEff(-100);
              //jetCon->SetBDiscrTrkCountHiPure(-100);

	      LogWarning("bTags")<<"something is wrong with b-tag jets!"<<
		"\n cor jets: "<<cc<<"   pt: "<< corJet.pt()<<"  eta: "<<corJet.eta()<<"  phi: "<<corJet.phi()<<
		"\n uncorrec: "<<cc<<"   pt: "<< myJet.pt()<<"  eta: "<<myJet.eta()<<"  phi: "<<myJet.phi()<<
		"\n b  b1 "<<cc<<"    pt: "<<jet_b_1->first->pt()<<"   eta "<<jet_b_1->first->eta()<<"   phi "<<jet_b_1->first->phi()<<"   b: "<<jet_b_1->second<<
		"\n b  b2 "<<cc<<"    pt: "<<jet_b_2->first->pt()<<"   eta "<<jet_b_2->first->eta()<<"   phi "<<jet_b_2->first->phi()<<"   b: "<<jet_b_2->second<<
		"\n"<<dR1<<"    "<<dR2<<
		"\n size:  pf: "<<PFJets->size()<<"    b1: "<<bTagHandle1->size()<<"  b2: "<<bTagHandle2->size();

	    }
	  cc++;

	  /*
	  jetCon->SetBDiscrTrkCountHiEff(jet_b_1->second);
	  jetCon->SetBDiscrTrkCountHiPure(jet_b_2->second);
	  
	  
	  jetCon->SetJetCorr(1, scale123); 	  
	  if (isRealData) 
	    {
	      float scaleRes = correctorRes->correction(corJet);
	      jetCon->SetJetCorr(4, scaleRes);
	    }

	  jetCon->SetP4(myJet.px(), myJet.py(), myJet.pz(), myJet.energy());
	  jetCon->SetVtx(-999.0, -999.0, -999.0);
	  jetCon->SetP4(myJet.px(), myJet.py(), myJet.pz(), myJet.energy());
	  jetCon->SetVtx(-999.0, -999.0, -999.0);
	  jetCon->SetChHadFrac(myJet.chargedHadronEnergyFraction());
	  jetCon->SetNeuHadFrac(myJet.neutralHadronEnergyFraction());
	  jetCon->SetChEmFrac(myJet.chargedEmEnergyFraction());
	  jetCon->SetNeuEmFrac(myJet.neutralEmEnergyFraction());
	  jetCon->SetNumConstit(myJet.chargedMultiplicity() + myJet.neutralMultiplicity());
	  jetCon->SetNumChPart(myJet.chargedMultiplicity());	  
	  */
	  //jetmomenta.push_back(myJet.p4());

	  pfCount++;  
	}
      else nLowPtJets++;
      
      jet_b_1++;
      jet_b_2++;

      //LogInfo("Higgs: ")<<" number of jets:  pt>5 "<< pfCount<<"  low pt<5 "<<nLowPtJets;
      
    }//End of pfJets
  
  //cout<<"dbg: out of pfJet loop"<<endl;
  //------- End of Jets----------------

  //----------------------------------
  //----------MET------------------
  //--------------------------------
  Handle<PFMETCollection> pfMEThandle;
  iEvent.getByLabel("pfMet", pfMEThandle);
  reco::PFMET iMET = pfMEThandle->front();
 
  Int_t metCount=0;
  //for (PFMETCollection::const_iterator iMET = MET->begin(); iMET != MET->end(); ++iMET) 
  //{}
  TCMET* metCon = new ((*recoMET)[metCount]) TCMET;
  /*
  metCon->SetSumEt(iMET.sumEt());
  metCon->SetMet(iMET.et());
  metCon->SetPhi(iMET.phi());
  metCon->SetPhotonEtFraction(iMET.photonEtFraction());
  metCon->SetElectronEtFraction(iMET.electronEtFraction());
  metCon->SetMuonEtFraction(iMET.muonEtFraction());
  metCon->SetNeutralHadronEtFraction(iMET.neutralHadronEtFraction());
  metCon->SetChargedHadronEtFraction(iMET.chargedHadronEtFraction());
  metCon->SetHFHadronEtFraction(iMET.HFHadronEtFraction());
  metCon->SetHFEMEtFraction(iMET.HFEMEtFraction());
  //metCount++;
  */
  
  //Type 1 corrected MET
  Handle<METCollection> pfMEThandle1;
  iEvent.getByLabel("metJESCorPFAK5", pfMEThandle1);
  reco::MET iMET1 = pfMEThandle1->front();
  //reco::PFMET iMET1;
  //iEvent.getByLabel("metJESCorPFAK5", iMET1);
  //reco::PFMET iMET1 = pfMEThandle1->front();
  
  TCMET* metCon1 = new ((*corrMET)[0]) TCMET;
  /*
  metCon1->SetSumEt(iMET1.sumEt());
  metCon1->SetMet(iMET1.et());
  metCon1->SetPhi(iMET1.phi());
  */
  //metCon1->SetPhotonEtFraction(iMET1.photonEtFraction());
  //metCon1->SetElectronEtFraction(iMET1.electronEtFraction());
  //metCon1->SetMuonEtFraction(iMET1.muonEtFraction());
  //metCon1->SetNeutralHadronEtFraction(iMET1.neutralHadronEtFraction());
  //metCon1->SetChargedHadronEtFraction(iMET1.chargedHadronEtFraction());
  //metCon1->SetHFHadronEtFraction(iMET1.HFHadronEtFraction());
  //metCon1->SetHFEMEtFraction(iMET1.HFEMEtFraction());
  
  //if(metCount>1)  LogError("MET: ")<<"Two many met objects! "<<metCount<<endl;

  //cout<<"dbg: out of pfJet loop"<<endl;
  //---------------End of MET -------------


  //----------------------------------------
  //------------ Muons --------------------
  //----------------------------------------
  Handle<MuonCollection> muons;
  iEvent.getByLabel(muonTag_, muons);

  int muCount = 0;
  for (MuonCollection::const_iterator mu = muons->begin(); mu != muons->end(); ++mu) 
    {
    if (mu->isGlobalMuon() && mu->pt() > 9)
      {
	TCMuon* muCon = new ((*recoMuons)[muCount]) TCMuon;
	/*
	muCon->Setp4(mu->px(), mu->py(), mu->pz(), mu->p());
	muCon->SetVtx(mu->globalTrack()->vx(),mu->globalTrack()->vy(),mu->globalTrack()->vz());
	muCon->SetCharge(mu->charge());
	muCon->SetisGLB(mu->isGlobalMuon());
	muCon->SetisTRK(mu->isTrackerMuon());
	muCon->SetnPXLHits(mu->globalTrack()->hitPattern().numberOfValidPixelHits());
	muCon->SetnTRKHits(mu->globalTrack()->hitPattern().numberOfValidTrackerHits());
	muCon->SetnValidMuHits(mu->globalTrack()->hitPattern().numberOfValidMuonHits());
	muCon->SetnMatchSeg(mu->numberOfMatches());
	muCon->SetNormChi2(mu->globalTrack()->normalizedChi2());
	muCon->SetCaloComp(mu->caloCompatibility());
	muCon->SetSegComp(muon::segmentCompatibility(*mu));
	muCon->SetEMIso(mu->isolationR03().emEt);
	muCon->SetHADIso(mu->isolationR03().hadEt);
	muCon->SetTRKIso(mu->isolationR03().sumPt);
	*/
	muCount++;
      }
  }
  //-----------End of muons ----------------------

  //Reduced Met
  /*
  //Double_t lepton1pterr = 1;
  //Double_t lepton2pterr = 1;

  //TLorentzVector themet(0,0,0,0);

  //  if (muCount>=2)  rmet_.compute(recoMuons[0]->P4(),lepton1pterr, recoMuons[1]->P4(),lepton2pterr, jetmomenta, themet );
  float reducedMET=rmet_.reducedMET();
  LorentzVector rmetP(rmet_.reducedMETComponents().first,rmet_.reducedMETComponents().second,0,reducedMET);
  std::pair<double, double> rmet_dil = rmet_.dileptonProjComponents();
  std::pair<double, double> rmet_delta = rmet_.dileptonPtCorrComponents();
  std::pair<double, double> rmet_rjet = rmet_.sumJetProjComponents();
  std::pair<double, double> rmet_met = rmet_.metProjComponents();
  std::pair<double, double> rmet_r = rmet_.recoilProjComponents();
  */

  
  //----------------------------------------
  //------------Electrons ------------------
  //----------------------------------------
  Handle<edm::ValueMap<float> > eIDValueMap;
  iEvent.getByLabel( electronIDMap_ , eIDValueMap );
  const edm::ValueMap<float> & eIDmap = * eIDValueMap ;

  Handle<GsfElectronCollection> electrons;
  iEvent.getByLabel(electronTag_, electrons);

  int elCount = 0;
  for (unsigned int i = 0; i < electrons->size(); i++)
    {
      edm::Ref<reco::GsfElectronCollection> electronRef(electrons,i);
      Int_t id = eIDmap[electronRef];
      if ( (id==5 || id==7)  && electronRef->pt() > 9)
	{
	  TCElectron* elCon = new ((*recoElectrons)[elCount]) TCElectron;
	  /*
	    elCon->Setp4(electronRef->px(),electronRef->py(),electronRef->pz(),electronRef->p());
	  elCon->SetVtx(electronRef->gsfTrack()->vx(),electronRef->gsfTrack()->vy(),electronRef->gsfTrack()->vz());
	  elCon->SetCharge(electronRef->charge());
	  elCon->SetNormChi2(electronRef->gsfTrack()->normalizedChi2());
	  elCon->SetEMIso(electronRef->dr03EcalRecHitSumEt());
	  elCon->SetHADIso(electronRef->dr03HcalTowerSumEt());
	  elCon->SetTRKIso(electronRef->dr03TkSumPt());

	  //double superClusterEta = electronRef->superCluster()->eta();

	  elCon->SetHadOverEm(electronRef->hadronicOverEm());
	  elCon->SetDPhiSuperCluster(electronRef->deltaPhiSuperClusterTrackAtVtx());
	  elCon->SetDEtaSuperCluster(electronRef->deltaEtaSuperClusterTrackAtVtx());
	  elCon->SetSigmaIetaIeta(electronRef->sigmaIetaIeta());
	  */
	  elCount++;
	}
  }

  Handle<bool> hcalNoiseFilterHandle;
  iEvent.getByLabel(hcalFilterTag_, hcalNoiseFilterHandle);
  if (hcalNoiseFilterHandle.isValid())  isNoiseHcal = !(Bool_t)(*hcalNoiseFilterHandle);

  //cout<<"noisy hcal? "<<isNoiseHcal<<endl;

  isDeadEcalCluster = kFALSE;
  //edm::InputTag ecalAnomalousFilterTag_("BE1214","anomalousECALVariables");
  Handle<AnomalousECALVariables> anomalousECALvarsHandle;
  iEvent.getByLabel(ecalAnomalousFilterTag_, anomalousECALvarsHandle);
  AnomalousECALVariables anomalousECALvars;
  if (anomalousECALvarsHandle.isValid()) {
    anomalousECALvars = *anomalousECALvarsHandle;
    isDeadEcalCluster = anomalousECALvars.isDeadEcalCluster();
  } else {
    // LogWarning("ECAL dead cluster filter: ") <<" Anomalous ECAL Vars not valid/found: " ;
  } 
  
  //if (isDeadEcalCluster) cout<<"dead ecal?  "<<isDeadEcalCluster<<endl;


  edm::Handle<BeamHaloSummary> TheBeamHaloSummary;
  iEvent.getByLabel("BeamHaloSummary",TheBeamHaloSummary);

  const BeamHaloSummary TheSummary = (*TheBeamHaloSummary.product() );
  //if( TheSummary.CSCTightHaloId())
  // cout << "This event has been identified as a halo event with tight CSC-based Halo Id" << endl;
  // if( TheSummary.CSCLooseHaloId())
  //  cout << "This event has been identified as a halo event with loose CSC-based Halo Id" << endl;

  isCSCTightHalo = TheSummary.CSCTightHaloId();
  isCSCLooseHalo = TheSummary.CSCLooseHaloId();

  /*// These filters are not validated so far
  isHcalTightHalo = TheSummary.HcalTightHaloId();
  isHcalLooseHalo = TheSummary.HcalLooseHaloId();

  isEcalTightHalo = TheSummary.EcalTightHaloId();
  isEcalLooseHalo = TheSummary.EcalLooseHaloId();

  isGlobalTightHalo = TheSummary.EcalTightHaloId();
  isGlobalLooseHalo = TheSummary.EcalLooseHaloId();
  */

  isScraping = isFilteredOutScraping(iEvent, iSetup, 10, 0.25); 

  /*-----------------*/
  /*---- Triggers ---*/
  /*-----------------*/

  triggerStatusMuon = 0x0;
  triggerStatusElectron = 0x0;

  if(isRealData)
    {
      edm::Handle<TriggerResults> hltR;
      edm::InputTag triggerResultsTag = InputTag("TriggerResults","","HLT");
      iEvent.getByLabel(triggerResultsTag, hltR);
      
      
      vector<string>  hlNames;
      
      const TriggerNames & triggerNames = iEvent.triggerNames(*hltR);
      hlNames = triggerNames.triggerNames();
      
      bool triggerPassed = false;
      
      for (uint i=0; i<hlNames.size(); ++i) 
	{
	  triggerPassed = triggerDecision(hltR, i);
	  if(triggerPassed)
	    {
	      if (hlNames[i].find("HLT_DoubleMu3")!=string::npos)
		triggerStatusMuon |= 0x01 << 0;
	      if (hlNames[i].find("HLT_DoubleMu4")!=string::npos)
		triggerStatusMuon |= 0x01 << 1;
	      if (hlNames[i].find("HLT_DoubleMu5")!=string::npos)
		triggerStatusMuon |= 0x01 << 2;
	      if (hlNames[i].find("HLT_DoubleMu6")!=string::npos)
		triggerStatusMuon |= 0x01 << 3;
	      

	      if (hlNames[i].find("HLT_Ele15_LW")!=string::npos)
		triggerStatusElectron |= 0x01 << 0;
	      if (hlNames[i].find("HLT_Ele15_SW")!=string::npos)
		triggerStatusElectron |= 0x01 << 1;
	      if (hlNames[i].find("HLT_Ele17_SW")!=string::npos)
		triggerStatusElectron |= 0x01 << 2;
	    }
	} 
    }







  //- Fill the Ntuple and clear the memory ---//
  rTree-> Fill();
    


  //primaryVtx         -> Clear("C");
  primaryVtxWithBS   -> Clear("C");
  primaryVtxDAWithBS -> Clear("C");
  primaryVtxDA       -> Clear("C");
  
  recoJets      -> Clear("C");
  recoMET       -> Clear("C");
  corrMET       -> Clear("C");
  recoElectrons -> Clear("C");
  recoMuons     -> Clear("C");
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
HiggsZZllnunuNtuplyzer::beginJob()
{
  //primaryVtx       = new TClonesArray("TCPrimaryVtx");
  primaryVtxWithBS   = new TClonesArray("TCPrimaryVtx");
  primaryVtxDAWithBS = new TClonesArray("TCPrimaryVtx");
  primaryVtxDA       = new TClonesArray("TCPrimaryVtx");

  recoJets         = new TClonesArray("TCJet");
  recoMET          = new TClonesArray("TCMET");
  corrMET          = new TClonesArray("TCMET");
  recoElectrons    = new TClonesArray("TCElectron");
  recoMuons        = new TClonesArray("TCMuon");

  //beamSpot         = new TClonesArray("TVector3");
 
  rTree = fs->make<TTree>("rTree","rTree");
  rTree->Branch("isRealData",&isRealData,"isRealData/O");
  rTree->Branch("eventNumber",&eventNumber, "eventNumber/i");
  rTree->Branch("runNumber",&runNumber, "runNumber/i");
  rTree->Branch("lumiSection",&lumiSection, "lumiSection/i");
  rTree->Branch("bunchCross",&bunchCross, "bunchCross/i");

  rTree->Branch("rhoFactor",&rhoFactor, "rhoFactor/F");
  rTree->Branch("isNoiseHcal",&isNoiseHcal,"isNoiseHcal/O");
  rTree->Branch("isDeadEcalCluster",&isDeadEcalCluster,"isDeadEcalCluster/O");

  rTree->Branch("beamSpot", &beamSpot, "x:y:z");
  //rTree->Branch("beamSpot", &beamSpot, 6400, 1);
  //rTree->Branch("primaryVtx",&primaryVtx, 6400, 0);
  rTree->Branch("primaryVtxWithBS",&primaryVtxWithBS, 6400, 0);
  rTree->Branch("primaryVtxDAWithBS",&primaryVtxDAWithBS, 6400, 0);
  rTree->Branch("primaryVtxDA",&primaryVtxDA, 6400, 0);

  rTree->Branch("recoJets",&recoJets, 6400, 0);
  rTree->Branch("recoElectrons",&recoElectrons, 6400, 0);
  rTree->Branch("recoMuons",&recoMuons, 6400, 0);
  rTree->Branch("recoMET",&recoMET, 6400, 0);
  rTree->Branch("corrMET",&corrMET, 6400, 0);


  rTree->Branch("isCSCTightHalo",&isCSCTightHalo,"isCSCTightHalo/O");
  rTree->Branch("isCSCLooseHalo",&isCSCLooseHalo,"isCSCLooseHalo/O");

  /*
  rTree->Branch("isHcalTightHalo",&isHcalTightHalo,"isHcalTightHalo/O");
  rTree->Branch("isHcalLooseHalo",&isHcalLooseHalo,"isHcalLooseHalo/O");

  rTree->Branch("isEcalTightHalo",&isEcalTightHalo,"isEcalTightHalo/O");
  rTree->Branch("isEcalLooseHalo",&isEcalLooseHalo,"isEcalLooseHalo/O");

  rTree->Branch("isGlobalTightHalo",&isGlobalTightHalo,"isGlobalTightHalo/O");
  rTree->Branch("isGlobalLooseHalo",&isGlobalLooseHalo,"isGlobalLooseHalo/O");
  */

  rTree->Branch("isScraping",&isScraping,"isScraping/O");

  rTree->Branch("triggerStatusMuon",&triggerStatusMuon, "triggerStatusMuon/i");
  rTree->Branch("triggerStatusElectron",&triggerStatusElectron, "triggerStatusElectron/i");
 
  //rTree->Branch("hltPrescale",hltPrescale, "hltPrescale[32]/i");

  //rTree->Branch("lumiDeadCount",&lumiDeadCount, "lumiDeadCount/f");
  //rTree->Branch("lumiLiveFrac",&lumiLiveFrac, "lumiLiveFrac/f");
  //rTree->Branch("intDelLumi",&intDelLumi, "intDelLumi/f");
  //rTree->Branch("ptHat",&ptHat, "ptHat/f");
  //rTree->Branch("qScale", &qScale, "qScale/f");
  //rTree->Branch("evtWeight", &evtWeight, "evtWeight/f");
  //rTree->Branch("crossSection", &crossSection, "crossSection/f");


}

// ------------ method called once each job just after ending the event loop  ------------
void 
HiggsZZllnunuNtuplyzer::endJob() {
}

float HiggsZZllnunuNtuplyzer::sumPtSquared(const Vertex & v)
{
  float sum = 0.;
  float pT;

  for (Vertex::trackRef_iterator it = v.tracks_begin(); it != v.tracks_end(); it++) {
    pT = (**it).pt();
    float epT=(**it).ptError(); pT=pT>epT ? pT-epT : 0;

    sum += pT*pT;
  }

  return sum;
}



bool HiggsZZllnunuNtuplyzer::isFilteredOutScraping( const edm::Event& iEvent, const edm::EventSetup& iSetup, int numtrack, double thresh)
{

  bool accepted = false;
  float fraction = 0;  
  // get GeneralTracks collection

  edm::Handle<reco::TrackCollection> tkRef;
  iEvent.getByLabel("generalTracks",tkRef);    
  const reco::TrackCollection* tkColl = tkRef.product();
 
  int numhighpurity=0;
  reco::TrackBase::TrackQuality _trackQuality = reco::TrackBase::qualityByName("highPurity");

  if(tkColl->size()>(UInt_t)numtrack){ 
    reco::TrackCollection::const_iterator itk = tkColl->begin();
    reco::TrackCollection::const_iterator itk_e = tkColl->end();
    for(;itk!=itk_e;++itk){
      if(itk->quality(_trackQuality)) numhighpurity++;
    }
    fraction = (float)numhighpurity/(float)tkColl->size();
    if(fraction>thresh) accepted=true;
  }else{
    //if less than 10 Tracks accept the event anyway    
    accepted= true;
  }
    
  /*
  if (debugOn) {
    int ievt = iEvent.id().event();
    int irun = iEvent.id().run();
    int ils = iEvent.luminosityBlock();
    int bx = iEvent.bunchCrossing();
    
    std::cout << "FilterOutScraping_debug: Run " << irun << " Event " << ievt << " Lumi Block " << ils << " Bunch Crossing " << bx << " Fraction " << fraction << " NTracks " << tkColl->size() << " Accepted " << accepted << std::endl;
  }
  */

  return !accepted;  //iif filtered out it's not accepted.

}

bool HiggsZZllnunuNtuplyzer::triggerDecision(edm::Handle<edm::TriggerResults> &hltR, int iTrigger)
{
  bool triggerPassed = false;
  if(hltR->wasrun(iTrigger) &&
     hltR->accept(iTrigger) &&
     !hltR->error(iTrigger) ){
    triggerPassed = true;
  }
  return triggerPassed;
}


//define this as a plug-in
DEFINE_FWK_MODULE(HiggsZZllnunuNtuplyzer);
