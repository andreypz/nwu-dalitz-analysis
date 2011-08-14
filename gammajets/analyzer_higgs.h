//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 12 10:25:05 2011 by ROOT version 5.27/06b
// from TTree rTree/rTree
// found on file: out_higgs.root
//////////////////////////////////////////////////////////

#ifndef analyzer_higgs_h
#define analyzer_higgs_h

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include "../src/TCPrimaryVtx.h"
#include "../src/TCMuon.h"
#include "../src/TCElectron.h"
#include "../src/TCJet.h"
#include "../src/TCMET.h"
#include "../src/TCPhoton.h"
#include "../src/TCTrigger.h"

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>

#define nCuts 10


class analyzer_higgs : public TSelector {
  ofstream nout[nCuts], fout[nCuts], ffout, ncout;
  
  double dz(TCPrimaryVtx *primVtx, TVector3* trackVtx, TVector3* trackMom);
  double dxy(TCPrimaryVtx *primVtx, TVector3* trackVtx, TVector3* trackMom);
  double higgsMT(Float_t pTll, Float_t Mll, Float_t MET, TVector2 ZptXY, TVector2 METXY);
  double projectedMET(Float_t MET, Float_t MET_phi, TLorentzVector Lepton1, TLorentzVector Lepton2);
  double pileUpCorrectedMET(Float_t projMET, Int_t nVtx);

  void CountEvents(UInt_t nEvents[], UInt_t nEventsPassNoiseFilter[]);
  
public :
  TFile* histoFile;
  
  TH1F *met_phi[nCuts], *met_et[nCuts], *met_over_qt[nCuts];
  TH1F *met1_phi[nCuts], *met1_et[nCuts], *met1_over_qt[nCuts];
  TH1F *met2_phi[nCuts], *met2_et[nCuts], *met2_over_qt[nCuts];
  TH1F *met3_phi[nCuts], *met3_et[nCuts], *met3_over_qt[nCuts];
  TH1F *met4_phi[nCuts], *met4_et[nCuts], *met4_over_qt[nCuts], *met4_sig[nCuts];
  TH1F *mu1_phi[nCuts], *mu1_eta[nCuts], *mu1_pt[nCuts];
  TH1F *mu2_phi[nCuts], *mu2_eta[nCuts], *mu2_pt[nCuts];
  TH1F *btag_hp[nCuts];

  TH1F *h_hlt_fired, *h_hlt_fired2;
  TH2F *h_hlt_fired_prescale, *h_hlt_fired_prescale2;
  TH1F *h_pho_pt, *h_pho_mass;
  TH1F *h_phow_pt, *h_phow_mass;
  
  TH1F *h_pho_eta, *h_pho_phi;
  TH1F *h_phow_eta, *h_phow_phi; 

  TH1F *h_pho1_pt, *h_pho1_mass;
  TH1F *h_phow1_pt, *h_phow1_mass;

  TH1F *h_pho1_eta, *h_pho1_phi;
  TH1F *h_phow1_eta, *h_phow1_phi;

  TH1F *h_pho2_pt, *h_pho2_mass;
  TH1F *h_phow2_pt, *h_phow2_mass;

  TH1F *h_pho2_eta, *h_pho2_phi;
  TH1F *h_phow2_eta, *h_phow2_phi;

  // histograms for event-level variable
  TH1F *h_phow1_mt[12], *h_phow1_met[12], *h_phow1_projMet[12], *h_phow1_metOverQt[12], *h_phow1_projMetOverQt[12], *h_phow1_nJets[12];
  TH1F *h_phow2_mt[12], *h_phow2_met[12], *h_phow2_projMet[12], *h_phow2_metOverQt[12], *h_phow2_projMetOverQt[12], *h_phow2_nJets[12];



  // histograms that hold the weights for photons to be applied for
  // bg estimation in ee and mumu. The combined histogram "ll" is 
  // obtained after combining the samples and needs additional weight 
  // factor for ee and mumu. First Check the ee, mumu... ll if thre is 
  // strong reason for it

  TH1F *h1_photonWeightEe, *h1_photonWeightMumu, *h1_photonWeightLl;

  Float_t photonWeightEe(Float_t pt);     // weight for the photons in ee bg estimate 
  Float_t photonWeightMumu(Float_t pt);   //  same for mumu
 

  // histograms to hold the ee, mumu mass spectra observed in data
  // when "photon mass" is generated directly from the histograms
 TH1F *h1_photonMassFromEe, *h1_photonMassFromMumu; 





   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Bool_t          isRealData;
   Bool_t          isDeadEcalCluster;
   Bool_t          isNoiseHcal;

   Bool_t          isScraping;
   Bool_t          isCSCTightHalo;
   Bool_t          isCSCLooseHalo;
   UInt_t triggerStatusMuon;
   UInt_t triggerStatusElectron;

   UInt_t          eventNumber;
   UInt_t          runNumber;
   UInt_t          lumiSection;
   UInt_t          bunchCross;
   Float_t         rhoFactor;
   Float_t         beamSpot_x;
   Float_t         beamSpot_y;
   Float_t         beamSpot_z;
   //TClonesArray    *primaryVtx;
   TClonesArray    *primaryVtxWithBS;
   TClonesArray    *primaryVtxDAWithBS;
   TClonesArray    *recoJets;
   TClonesArray    *recoElectrons;
   TClonesArray    *recoMuons;
   TClonesArray    *recoMET;
   TClonesArray    *corrMET;
   TClonesArray    *recoPhotons; 
   TClonesArray    *hltTriggerResults;; 

   // List of branches
   TBranch        *b_isRealData;   //!
   TBranch        *b_isDeadEcalCluster;   //!
   TBranch        *b_isNoiseHcal;   //!
   TBranch          *b_isScraping;
   TBranch          *b_isCSCTightHalo;
   TBranch         *b_isCSCLooseHalo;
   TBranch  *b_triggerStatusMuon;
   TBranch  *b_triggerStatusElectron;

   TBranch        *b_eventNumber;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_lumiSection;   //!
   TBranch        *b_bunchCross;   //!
   TBranch        *b_rhoFactor;   //!
   TBranch        *b_beamSpot;   //!
   //TBranch        *b_primaryVtx;   //!
   TBranch        *b_primaryVtxWithBS;   //!
   TBranch        *b_primaryVtxDAWithBS;   //!
   TBranch        *b_recoJets;   //!
   TBranch        *b_recoElectrons;   //!
   TBranch        *b_recoMuons;   //!
   TBranch        *b_recoMET;   //!
   TBranch        *b_corrMET;   //!
   TBranch        *b_recoPhotons;   //!  
   TBranch        *b_hltTriggerResults;   //!  

   analyzer_higgs(TTree * /*tree*/ =0) { }
   virtual ~analyzer_higgs() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(analyzer_higgs,0);
};

#endif

#ifdef analyzer_higgs_cxx
void analyzer_higgs::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   //primaryVtx = 0;
   primaryVtxWithBS = 0;
   primaryVtxDAWithBS = 0;
   recoJets = 0;
   recoElectrons = 0;
   recoMuons = 0;
   recoMET = 0;
   corrMET = 0;
   recoPhotons = 0;
   hltTriggerResults = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("isRealData", &isRealData, &b_isRealData);
   fChain->SetBranchAddress("isDeadEcalCluster", &isDeadEcalCluster, &b_isDeadEcalCluster);
   fChain->SetBranchAddress("isNoiseHcal", &isNoiseHcal, &b_isNoiseHcal);
   fChain->SetBranchAddress("isCSCTightHalo", &isCSCTightHalo, &b_isCSCTightHalo);
   fChain->SetBranchAddress("isCSCLooseHalo", &isCSCLooseHalo, &b_isCSCLooseHalo);
   fChain->SetBranchAddress("isScraping", &isScraping, &b_isScraping);
   fChain->SetBranchAddress("triggerStatusMuon", &triggerStatusMuon, &b_triggerStatusMuon);
   fChain->SetBranchAddress("triggerStatusElectron", &triggerStatusElectron, &b_triggerStatusElectron);

   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
   fChain->SetBranchAddress("bunchCross", &bunchCross, &b_bunchCross);
   fChain->SetBranchAddress("rhoFactor", &rhoFactor, &b_rhoFactor);
   fChain->SetBranchAddress("beamSpot", &beamSpot_x, &b_beamSpot);
   //fChain->SetBranchAddress("primaryVtx", &primaryVtx, &b_primaryVtx);
   fChain->SetBranchAddress("primaryVtxWithBS", &primaryVtxWithBS, &b_primaryVtxWithBS);
   fChain->SetBranchAddress("primaryVtxDAWithBS", &primaryVtxDAWithBS, &b_primaryVtxDAWithBS);
   fChain->SetBranchAddress("recoJets", &recoJets, &b_recoJets);
   fChain->SetBranchAddress("recoElectrons", &recoElectrons, &b_recoElectrons);
   fChain->SetBranchAddress("recoMuons", &recoMuons, &b_recoMuons);
   fChain->SetBranchAddress("recoMET", &recoMET, &b_recoMET);
   fChain->SetBranchAddress("corrMET", &corrMET, &b_corrMET);
   fChain->SetBranchAddress("recoPhotons", &recoPhotons, &b_recoPhotons);
   fChain->SetBranchAddress("hltTriggerResults", &hltTriggerResults, &b_hltTriggerResults);
}

Bool_t analyzer_higgs::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef analyzer_higgs_cxx
