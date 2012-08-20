// $Id: higgsAnalyzer.h,v 1.24 2012/08/04 23:36:49 andrey Exp $

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

#include "TClonesArray.h"
#include "../src/TCJet.h"
#include "../src/TCMET.h"
#include "../src/TCElectron.h"
#include "../src/TCMuon.h"
#include "../src/TCTau.h"
#include "../src/TCPhoton.h"
#include "../src/TCGenParticle.h"
#include "../src/TCGenJet.h"
#include "../src/TCPrimaryVtx.h"
//#include "../src/TCTrigger.h"
#include "../src/TCTriggerObject.h"

#include "../plugins/WeightUtils.h"
#include "../plugins/TriggerSelector.h"
#include "../plugins/ZedEventsLibrary.h"
#include "../plugins/rochcor.h"

#define nC 30  //nCuts in the analysis. make plots after each cut

class higgsAnalyzer : public TSelector {

 private:

  TFile* histoFile;
  //Variables to Fill into histos. Have to be global

  TString samp;  
  Float_t qT, diEta, diPhi, Mll, Mll_EB, Mll_EE, Mll_EX;
  
  Float_t MET, pfMET, pfMET1, puCorrMET, projMET, ZprojMET, redMET1, redMET2, compMET;
  Float_t pfMET_lg, pfMET_recoil;
  Float_t MET_phi, MET1_phi;
  Float_t MT, MT1, MTZ, pTll;
  Float_t METqt, MET1qt, projMETqt;
  Int_t nVtx, nVtxTotal;
  Float_t nDofVtx1, nDofVtx2;
  Int_t nJets, nJetsEta24, nJetsB , nJetsBssv, nJetsB25, nJetsB30;
  Int_t nGamma;
  Float_t dPhiClos1, dPhiClos2;
  Float_t lep1_eta, lep1_phi, lep1_pt;
  Float_t lep2_eta, lep2_phi, lep2_pt;
  Float_t dRjetlep1, dRjetlep2, ptLeadJet, etaLeadJet, phiLeadJet, ptLeadBJet;
  Float_t hig_M, hig_Pt;

  ofstream nout[nC], fout[nC], ffout, ncout;
  TTree * _kinTree;

  TH1F *evt_byCut, *evt_byCut_raw;
  TH1F *mva_discr[3];
  TH2F *evt_libQt;
  TH1F *mtZ[nC], *mt0[nC], *mt1[nC], *mt2[nC], *mt3[nC], *mt4[nC];
  TH1F *met0_phi[nC], *met0_et[nC], *met0_over_qt[nC];
  TH1F *met1_phi[nC], *met1_et[nC], *met1_over_qt[nC];
  TH1F *met2_phi[nC], *met2_et[nC], *met2_over_qt[nC];
  TH1F *met3_phi[nC], *met3_et[nC], *met3_over_qt[nC];
  TH1F *met4_phi[nC], *met4_et[nC], *met4_over_qt[nC], *met4_puSig[nC];
  TH1F *met5_et[nC],*met6_et[nC],*met7_et[nC],*met8_et[nC], *met9_et[nC], *met10_et[nC];
  TH1F *l1_phi[nC], *l1_eta[nC], *l1_pt[nC];
  TH1F *l2_phi[nC], *l2_eta[nC], *l2_pt[nC];
  TH1F *btag_hp[nC];
  TH1F *di_qt[nC], *di_eta[nC], *di_phi[nC], *di_mass[nC], *di_mass_EB[nC], *di_mass_EE[nC], *di_mass_EX[nC];
  TH1F *di_dPhiMet[nC];
  TH1F *jet_N[nC], *jet_N24[nC], *jet_dRlep1[nC], *jet_dRlep2[nC], *jet_pt[nC], *jet_eta[nC], *jet_phi[nC];
  TH1F *jet_b_N[nC], *jet_b_Nssv[nC], *jet_b_N25[nC],*jet_b_N30[nC], *jet_b_pt[nC];
  TH1F *ph_nGamma[nC];
  TH1F *l0_dPhi[nC], *l0_dEta[nC], *l0_dR[nC], *l0_ptRatio[nC];
  TH1F *met1_dPhiLeadJet1[nC], *met1_dPhiLeadJet2[nC], *met1_dPhiClosJet1[nC], *met1_dPhiClosJet2[nC];
  TH1F *met1_lg[nC], *met1_recoil_lg[nC];
//  TH1F *met2_dPhiLeadJet1[nC], *met2_dPhiLeadJet2[nC], *met2_dPhiClosJet1[nC], *met2_dPhiClosJet2[nC];

  TH2F *met0_et_ovQt[nC], *met1_et_ovQt[nC], *met2_et_ovQt[nC], *met3_et_ovQt[nC], *met4_et_ovQt[nC];
  TH2F *mtZ_met2[nC], *mt2_met2[nC];
  TH2F *mtZ_met3[nC], *mt2_met3[nC];

  TH1F *vtx_nPV_tot[nC], *vtx_nPV_raw[nC], *vtx_nPV_weight[nC], *vtx_ndof_1[nC], *vtx_ndof_2[nC];

  TH1F *run_events[nC];
  TH1F *evt_weight[nC];

  TH1F *higgs_pt[nC], *higgs_mass[nC], *higgs_w_pt[nC], *higgs_w_mass[nC];

  rochcor *roch;
  ZedEventsLibrary *zLib;
  WeightUtils *weighter;
  TriggerSelector *triggerSelector;


 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  
  // Declaration of leaf types                                                                                                                                 
  TClonesArray    *recoJets;
  TClonesArray    *recoJPT;
  TClonesArray    *recoElectrons;
  TClonesArray    *recoMuons;
  TClonesArray    *pfMuons;
  TClonesArray    *recoTaus;
  TClonesArray    *recoPhotons;
  TClonesArray    *pfPhotons;
  TCMET           *recoMET;
  TCMET           *recoMETNoPU;
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
  Bool_t          isScraping;
  Bool_t          isNoiseHcal;
  Bool_t          isDeadEcalCluster;
  Bool_t          isCSCTightHalo;
  Bool_t          isCSCLooseHalo;
  Float_t         ptHat;
  Float_t         qScale;
  Float_t         evtWeight;
  Float_t         rhoFactor;
  ULong64_t       triggerStatus;
  int          hltPrescale[64];
  
  // List of branches                                                                                                                                          
  TBranch        *b_recoJets;   //! 
  TBranch        *b_recoJPT;   //!
  TBranch        *b_recoElectrons;   //!
  TBranch        *b_recoMuons;   //!                                                                                                                           
  TBranch        *b_pfMuons;   //!                                                                                                                             
  TBranch        *b_recoTaus;   //!                                                                                                                            
  TBranch        *b_recoPhotons;   //!                                                                                                                         
  TBranch        *b_pfPhotons;   //!                                                                                                                           
  TBranch        *b_recoMET;   //!                                                                                                                             
  TBranch        *b_recoMETNoPU;   //!                                                                                                                         
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
  TBranch        *b_isScraping;   //!                                                                                                                          
  TBranch        *b_isNoiseHcal;   //!                                                                                                                         
  TBranch        *b_isDeadEcalCluster;   //!                                                                                                                   
  TBranch        *b_isCSCTightHalo;   //!                                                                                                                      
  TBranch        *b_isCSCLooseHalo;   //!                                                                                                                      
  TBranch        *b_ptHat;   //!                                                                                                                               
  TBranch        *b_qScale;   //!                                                                                                                              
  TBranch        *b_evtWeight;   //!                                                                                                                           
  TBranch        *b_rhoFactor;   //!                                                                                                                           
  TBranch        *b_triggerStatus;   //!                                                                                                                       
  TBranch        *b_hltPrescale;   //!                                                                                                                         

  
  //For counting events
  //        int          nEvents[16];
  //For counting weighted events
  float nEventsWeighted[nC];
  
  UInt_t nEvents[nC];
  UInt_t nEventsPassNoiseFilter[6][nC]; //1-Hcal, 2-Ecal, 3- Scraping, 4-CSCTight, 5-CSCLoose,   0-passed all
  
  higgsAnalyzer(TTree * /*tree*/ =0) { }
  virtual ~higgsAnalyzer() { }
  virtual int     Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  //virtual void    SlaveBegin(TTree *tree) { TString option = GetOption();};
  virtual void    Init(TTree *tree);
  virtual bool    Notify();
  virtual bool    Process(Long64_t entry);
  virtual int     GetEntry(Long64_t entry, int getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate() {};
  virtual void    Terminate();
  
  virtual bool    CosmicMuonFilter(TCMuon , TCMuon );
  virtual float   CalculateTransMass(TLorentzVector p1, TLorentzVector p2);
  virtual float   CalculateTransMassAlt(TLorentzVector p1, TLorentzVector p2);
  virtual float   DeltaPhiJetMET(TLorentzVector , std::vector<TLorentzVector> ); 

  
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
  recoJPT = 0;
  recoElectrons = 0;
  recoMuons = 0;
  pfMuons = 0;
  recoTaus = 0;
  recoPhotons = 0;
  pfPhotons = 0;
  recoMET = 0;
  recoMETNoPU = 0;
  triggerObjects = 0;
  genJets = 0;
  genParticles = 0;
  primaryVtx = 0;
  beamSpot = 0;


  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress("recoJets", &recoJets, &b_recoJets);
  fChain->SetBranchAddress("recoJPT", &recoJPT, &b_recoJPT);
  fChain->SetBranchAddress("recoElectrons", &recoElectrons, &b_recoElectrons);
  fChain->SetBranchAddress("recoMuons", &recoMuons, &b_recoMuons);
  fChain->SetBranchAddress("pfMuons", &pfMuons, &b_pfMuons);
  fChain->SetBranchAddress("recoTaus", &recoTaus, &b_recoTaus);
  fChain->SetBranchAddress("recoPhotons", &recoPhotons, &b_recoPhotons);
  fChain->SetBranchAddress("pfPhotons", &pfPhotons, &b_pfPhotons);
  fChain->SetBranchAddress("recoMET", &recoMET, &b_recoMET);
  fChain->SetBranchAddress("recoMETNoPU", &recoMETNoPU, &b_recoMETNoPU);
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
  fChain->SetBranchAddress("isScraping", &isScraping, &b_isScraping);
  fChain->SetBranchAddress("isNoiseHcal", &isNoiseHcal, &b_isNoiseHcal);
  fChain->SetBranchAddress("isDeadEcalCluster", &isDeadEcalCluster, &b_isDeadEcalCluster);
  fChain->SetBranchAddress("isCSCTightHalo", &isCSCTightHalo, &b_isCSCTightHalo);
  fChain->SetBranchAddress("isCSCLooseHalo", &isCSCLooseHalo, &b_isCSCLooseHalo);
  fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat);
  fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
  fChain->SetBranchAddress("evtWeight", &evtWeight, &b_evtWeight);
  fChain->SetBranchAddress("rhoFactor", &rhoFactor, &b_rhoFactor);
  fChain->SetBranchAddress("triggerStatus", &triggerStatus, &b_triggerStatus);
  fChain->SetBranchAddress("hltPrescale", hltPrescale, &b_hltPrescale);

}

bool higgsAnalyzer::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

#endif 
