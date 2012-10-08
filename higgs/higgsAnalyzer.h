//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 13 14:46:29 2012 by ROOT version 5.27/06b
// from TTree eventTree/eventTree
// found on file: /uscms/home/andreypz/nobackup/nuTuple_new.root
//////////////////////////////////////////////////////////

#ifndef higgsAnalyzer_h
#define higgsAnalyzer_h

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TString.h>
#include <TKey.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TVector2.h>
#include <TProfile.h>
#include <TRandom3.h>

#include "../src/TCPhysObject.h"
#include "../src/TCJet.h"
#include "../src/TCMET.h"
#include "../src/TCElectron.h"
#include "../src/TCMuon.h"
#include "../src/TCTau.h"
#include "../src/TCPhoton.h"
#include "../src/TCGenParticle.h"
#include "../src/TCGenJet.h"
#include "../src/TCPrimaryVtx.h"
#include "../src/TCTriggerObject.h"

#include "../plugins/MetDefinitions.h"
#include "../plugins/WeightUtils.h"
#include "../plugins/TriggerSelector.h"
#include "../plugins/ZedEventsLibrary.h"
#include "../plugins/rochcor.h"
#include "../plugins/HistManager.h"


#define nC 35  //nCuts in the analysis. make plots after each cut

struct muIdAndIsoCuts{
  Bool_t IsPF;
  Bool_t IsGLB;
  Float_t ptErrorOverPt;
  Float_t NumberOfValidMuonHits;
  Float_t NumberOfValidTrackerHits;
  Float_t NumberOfValidPixelHits;
  Float_t NumberOfMatches;
  Float_t NormalizedChi2;
  Float_t dxy;
  Float_t dz;
  Float_t chIso04;
  Float_t nhIso04;
  Float_t phIso04;
  Float_t pfIso04;
};

struct elIdAndIsoCuts{
  //broken into [0] barrel and [1] endcap
  float ptErrorOverPt[2];
  float dEtaIn[2];
  float dPhiIn[2];
  float sigmaIetaIeta[2];
  float HadOverEm[2];
  float dxy[2];
  float dz[2];
  float fabsEPDiff[2];
  float ConversionMissHits[2];
  float PassedConversionProb[2];
  float pfIso04[2];
};

class higgsAnalyzer : public TSelector {

 private:

  TFile* histoFile;
  TTree* thisTree;

  muIdAndIsoCuts muIdAndIsoCutsTight, muIdAndIsoCutsLoose, muIdAndIsoCutsSoft;
  elIdAndIsoCuts elIdAndIsoCutsTight, elIdAndIsoCutsLoose;
  //Variables to Fill into histos. Have to be global
  
  Float_t qT, diEta, diPhi, diAngle, deltaPhiLep, Mll, Mll_EB, Mll_EE, Mll_EX;  
  Float_t MET, pfMET, pfMET1, puCorrMET, projMET, ZprojMET, redMET1, redMET2, compMET;
  Float_t pfMET_lg, pfMET_recoil;
  Float_t dPhiMetDiLep, metProjOnQt, metPerpQt, metOverQt;

  Float_t MET_phi, MET1_phi;
  Float_t MT, MT1, MTZ, pTll;
  Int_t nVtx, nVtxTotal;
  Float_t nDofVtx1, nDofVtx2;
  Int_t nJets, nJets15, nJetsEta24, nJetsB , nJetsBssv, nJetsB25, nJetsB30;
  Int_t nGamma;
  Float_t dPhiClos1, dPhiClos2;
  Float_t lep1_eta, lep1_phi, lep1_pt;
  Float_t lep2_eta, lep2_phi, lep2_pt;
  Float_t lep_dPhi, lep_dEta, lep_dR, lep_angle, lep_angleLog, lep_ptRatio;
  Float_t dRjetlep1, dRjetlep2, ptLeadJet, etaLeadJet, phiLeadJet, ptLeadBJet;
  Float_t ptTrailJet, etaTrailJet, phiTrailJet;
  Float_t deltaEtaDiJet, massDiJet, zeppDiJetDiLep, zeppDiJetDiLep_y;
  Float_t hig_M, hig_Pt;

  ofstream nout[nC], fout[nC], ffout, ncout;
  TTree * _mvaTree;

  HistManager *hists;
  rochcor *roch;
  ZedEventsLibrary *zLib;
  WeightUtils *weighter;
  TriggerSelector *triggerSelector;

  float nEventsWeighted[nC];
  
  UInt_t nEvents[nC];
  UInt_t nEventsPassNoiseFilter[6][nC]; //1-Hcal, 2-Ecal, 3- Scraping, 4-CSCTight, 5-CSCLoose,   0-passed all
  UInt_t totEvents;  

public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types                                                                                                               
   TClonesArray    *recoJets;
   //TClonesArray    *recoJPT;
   TClonesArray    *recoElectrons;
   TClonesArray    *recoMuons;
   TClonesArray    *recoTaus;
   TClonesArray    *recoPhotons;
   TCMET           *recoMET;
   TClonesArray    *triggerObjects;
   TClonesArray    *genJets;
   TClonesArray    *genParticles;
   TClonesArray    *primaryVtx;
   TVector3        *beamSpot;
   Int_t           nPUVertices;
   Float_t         nPUVerticesTrue;
   Bool_t          isRealData;
   UInt_t          runNumber;
   ULong64_t       eventNumber;
   UInt_t          lumiSection;
   UInt_t          bunchCross;
   Float_t         ptHat;
   Float_t         qScale;
   Float_t         evtWeight;
   Float_t         rhoFactor;
   Float_t         rho25Factor;
   Float_t         rhoMuFactor;
   ULong64_t       triggerStatus;
   UInt_t          hltPrescale[64];
   Bool_t          NoiseFilters_isScraping;
   Bool_t          NoiseFilters_isNoiseHcalHBHE;
   Bool_t          NoiseFilters_isNoiseHcalLaser;
   Bool_t          NoiseFilters_isNoiseEcalTP;
   Bool_t          NoiseFilters_isNoiseEcalBE;
   Bool_t          NoiseFilters_isCSCTightHalo;
   Bool_t          NoiseFilters_isCSCLooseHalo;

   // List of branches
   TBranch        *b_recoJets;   //!
   TBranch        *b_recoElectrons;   //!
   TBranch        *b_recoMuons;   //!
   TBranch        *b_recoTaus;   //!
   TBranch        *b_recoPhotons;   //!
   TBranch        *b_recoMET;   //!
   TBranch        *b_triggerObjects;   //!
   TBranch        *b_genJets;   //!
   TBranch        *b_genParticles;   //!
   TBranch        *b_primaryVtx;   //!
   TBranch        *b_beamSpot;   //!
   TBranch        *b_nPUVertices;   //!
   TBranch        *b_nPUVerticesTrue;   //! 
   TBranch        *b_isRealData;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!  
   TBranch        *b_lumiSection;   //! 
   TBranch        *b_bunchCross;   //!
   TBranch        *b_ptHat;   //!   
   TBranch        *b_qScale;   //!
   TBranch        *b_evtWeight;   //! 
   TBranch        *b_rhoFactor;   //! 
   TBranch        *b_rho25Factor;   //!
   TBranch        *b_rhoMuFactor;   //!
   TBranch        *b_triggerStatus;   //!
   TBranch        *b_hltPrescale;   //!
   TBranch        *b_NoiseFilters;   //!

   higgsAnalyzer(TTree * /*tree*/ =0) { }
   virtual ~higgsAnalyzer() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   //virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   //virtual void    SlaveTerminate();
   virtual void    Terminate();



   //virtual bool    CosmicMuonFilter(TCMuon , TCMuon );
   virtual float   CalculateTransMass(TLorentzVector p1, TLorentzVector p2);
   virtual float   CalculateTransMassAlt(TLorentzVector p1, TLorentzVector p2);
   virtual float   DeltaPhiJetMET(TLorentzVector , std::vector<TLorentzVector> ); 

   virtual bool PassMuonIdAndIso(TCMuon *l, muIdAndIsoCuts c, TVector3 *pv);
   virtual bool PassElectronIdAndIso(TCElectron *l, elIdAndIsoCuts c, TVector3 *pv);
  
   virtual float EffAreaMuon(TCMuon *lep, TString, Bool_t, Int_t);
   virtual float EffAreaElectron(TCElectron *lep, TString, Bool_t, Int_t);
   void DumpMuonInfo(TCMuon *lep, TVector3 *v, Float_t, Float_t, Float_t);
   void DumpElectronInfo(TCElectron *lep, TVector3 *v, Float_t, Float_t, Float_t);
   void FillHistosBasic(Int_t, Double_t);
   void FillHistosFull(Int_t, Double_t);
   void CountEvents(Int_t);
   void PrintOut(Int_t, Bool_t isShort=kFALSE);
   void PrintOutNoisy(Int_t);
   // pair<int, int> GetQtBin(Float_t, Float_t);

   ClassDef(higgsAnalyzer,0);
};

#endif

#ifdef higgsAnalyzer_cxx
void higgsAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   recoJets = 0;
   //recoJPT = 0;
   recoElectrons = 0;
   recoMuons = 0;
   recoTaus = 0;
   recoPhotons = 0;
   recoMET = 0;
   //recoMETNoPU = 0;
   triggerObjects = 0;
   genJets = 0;
   genParticles = 0;
   primaryVtx = 0;
   beamSpot = 0;
   // Set branch addresses and branch pointers

   totEvents = 0;
   if (!tree) return;
   thisTree    = tree; 
   fChain = tree;
   fChain->SetMakeClass(1);


   fChain->SetBranchAddress("recoJets", &recoJets, &b_recoJets);
   //fChain->SetBranchAddress("recoJPT", &recoJPT, &b_recoJPT);
   fChain->SetBranchAddress("recoElectrons", &recoElectrons, &b_recoElectrons);
   fChain->SetBranchAddress("recoMuons", &recoMuons, &b_recoMuons);
   fChain->SetBranchAddress("recoTaus", &recoTaus, &b_recoTaus);
   fChain->SetBranchAddress("recoPhotons", &recoPhotons, &b_recoPhotons);
   fChain->SetBranchAddress("recoMET", &recoMET, &b_recoMET);
   //fChain->SetBranchAddress("recoMETNoPU", &recoMETNoPU, &b_recoMETNoPU);
   fChain->SetBranchAddress("triggerObjects", &triggerObjects, &b_triggerObjects);
   fChain->SetBranchAddress("genJets", &genJets, &b_genJets);
   fChain->SetBranchAddress("genParticles", &genParticles, &b_genParticles);
   fChain->SetBranchAddress("primaryVtx", &primaryVtx, &b_primaryVtx);
   fChain->SetBranchAddress("beamSpot", &beamSpot, &b_beamSpot);
   fChain->SetBranchAddress("nPUVertices", &nPUVertices, &b_nPUVertices);
   fChain->SetBranchAddress("nPUVerticesTrue", &nPUVerticesTrue, &b_nPUVerticesTrue);
   fChain->SetBranchAddress("isRealData", &isRealData, &b_isRealData);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
   fChain->SetBranchAddress("bunchCross", &bunchCross, &b_bunchCross);
   fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat);
   fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
   fChain->SetBranchAddress("evtWeight", &evtWeight, &b_evtWeight);
   fChain->SetBranchAddress("rhoFactor", &rhoFactor, &b_rhoFactor);
   fChain->SetBranchAddress("rho25Factor", &rho25Factor, &b_rho25Factor);
   fChain->SetBranchAddress("rhoMuFactor", &rhoMuFactor, &b_rhoMuFactor);
   fChain->SetBranchAddress("triggerStatus", &triggerStatus, &b_triggerStatus);
   fChain->SetBranchAddress("hltPrescale", hltPrescale, &b_hltPrescale);
   fChain->SetBranchAddress("NoiseFilters", &NoiseFilters_isScraping, &b_NoiseFilters);

   //Cuts

   //electrons are two types: Barrel/Endcap 
   elIdAndIsoCutsTight.ptErrorOverPt[0] = 0.8;
   elIdAndIsoCutsTight.dxy[0]           = 0.02;
   elIdAndIsoCutsTight.dz[0]            = 0.1;
   elIdAndIsoCutsTight.dPhiIn[0]        = 0.06;
   elIdAndIsoCutsTight.dEtaIn[0]        = 0.004;
   elIdAndIsoCutsTight.sigmaIetaIeta[0] = 0.01;
   elIdAndIsoCutsTight.HadOverEm[0]     = 0.12;
   elIdAndIsoCutsTight.dxy[0]           = 0.02;
   elIdAndIsoCutsTight.dz[0]            = 0.1;
   elIdAndIsoCutsTight.fabsEPDiff[0]    = 0.05;
   elIdAndIsoCutsTight.pfIso04[0]       = 0.15;

   elIdAndIsoCutsTight.ptErrorOverPt[1] = 0.8;
   elIdAndIsoCutsTight.dxy[1]           = 0.02;
   elIdAndIsoCutsTight.dz[1]            = 0.1;
   elIdAndIsoCutsTight.dPhiIn[1]        = 0.03;
   elIdAndIsoCutsTight.dEtaIn[1]        = 0.007;
   elIdAndIsoCutsTight.sigmaIetaIeta[1] = 0.03;
   elIdAndIsoCutsTight.HadOverEm[1]     = 0.10;
   elIdAndIsoCutsTight.dxy[1]           = 0.02;
   elIdAndIsoCutsTight.dz[1]            = 0.1;
   elIdAndIsoCutsTight.fabsEPDiff[1]    = 0.05;
   elIdAndIsoCutsTight.pfIso04[1]       = 0.15;


   elIdAndIsoCutsLoose.ptErrorOverPt[0] = 99999;
   elIdAndIsoCutsLoose.dxy[0]           = 0.02;
   elIdAndIsoCutsLoose.dz[0]            = 0.1;
   elIdAndIsoCutsLoose.dPhiIn[0]        = 0.06;
   elIdAndIsoCutsLoose.dEtaIn[0]        = 0.004;
   elIdAndIsoCutsLoose.sigmaIetaIeta[0] = 0.01;
   elIdAndIsoCutsLoose.HadOverEm[0]     = 0.12;
   elIdAndIsoCutsLoose.dxy[0]           = 0.02;
   elIdAndIsoCutsLoose.dz[0]            = 0.1;
   elIdAndIsoCutsLoose.fabsEPDiff[0]    = 99999;
   elIdAndIsoCutsLoose.pfIso04[0]       = 0.15;

   elIdAndIsoCutsLoose.ptErrorOverPt[1] = 99999;
   elIdAndIsoCutsLoose.dxy[1]           = 0.02;
   elIdAndIsoCutsLoose.dz[1]            = 0.1;
   elIdAndIsoCutsLoose.dPhiIn[1]        = 0.03;
   elIdAndIsoCutsLoose.dEtaIn[1]        = 0.007;
   elIdAndIsoCutsLoose.sigmaIetaIeta[1] = 0.03;
   elIdAndIsoCutsLoose.HadOverEm[1]     = 0.12;
   elIdAndIsoCutsLoose.dxy[1]           = 0.02;
   elIdAndIsoCutsLoose.dz[1]            = 0.1;
   elIdAndIsoCutsLoose.fabsEPDiff[1]    = 99999;
   elIdAndIsoCutsLoose.pfIso04[1]       = 0.15;


   muIdAndIsoCutsTight.ptErrorOverPt            = 0.1;
   muIdAndIsoCutsTight.NumberOfValidMuonHits    = 0;
   muIdAndIsoCutsTight.NumberOfValidTrackerHits = 10;
   muIdAndIsoCutsTight.NumberOfValidPixelHits   = 0;
   muIdAndIsoCutsTight.NumberOfMatches          = 1;
   muIdAndIsoCutsTight.NormalizedChi2           = 10;
   muIdAndIsoCutsTight.dxy             = 0.02;
   muIdAndIsoCutsTight.dz              = 0.05;
   muIdAndIsoCutsTight.pfIso04         = 0.2;

   muIdAndIsoCutsLoose.ptErrorOverPt            = 9999;
   muIdAndIsoCutsLoose.NumberOfValidMuonHits    = -1;
   muIdAndIsoCutsLoose.NumberOfValidTrackerHits = 2;
   muIdAndIsoCutsLoose.NumberOfValidPixelHits   = -1;
   muIdAndIsoCutsLoose.NumberOfMatches          = -1;
   muIdAndIsoCutsLoose.NormalizedChi2           = 9999;
   muIdAndIsoCutsLoose.dxy            = 0.2;
   muIdAndIsoCutsLoose.dz             = 0.2;
   muIdAndIsoCutsLoose.pfIso04        = 9999;

   muIdAndIsoCutsSoft.ptErrorOverPt            = 9999;
   muIdAndIsoCutsSoft.NumberOfValidMuonHits    = -1;
   muIdAndIsoCutsSoft.NumberOfValidTrackerHits = 10;
   muIdAndIsoCutsSoft.NumberOfValidPixelHits   = -1;
   muIdAndIsoCutsSoft.NumberOfMatches          = -1;
   muIdAndIsoCutsSoft.NormalizedChi2           = 9999;
   muIdAndIsoCutsSoft.dxy                 = 0.3;
   muIdAndIsoCutsSoft.dz                  = 3.;
   muIdAndIsoCutsSoft.pfIso04             = 0.2;
}

Bool_t higgsAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

  UInt_t initEvents;
  TFile   *inFile     = thisTree->GetCurrentFile();
  TTree   *jobTree    = (TTree*)inFile->Get("ntupleProducer/jobTree");
  TBranch *eventsInFile    =  jobTree->GetBranch("nEvents");

  eventsInFile->SetAddress(&initEvents);
  eventsInFile->GetEntry(0);

  totEvents += initEvents;
  cout<<" Events in file: "<<initEvents<<"  total processed: "<<totEvents<<endl;

  MET = 0;
  FillHistosBasic(0, initEvents);

  return kTRUE;
}

#endif // #ifdef higgsAnalyzer_cxx
