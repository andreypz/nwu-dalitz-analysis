//////////////////////////////////////////////////////////// This class has been automatically generated on
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

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>

#include <TKey.h>
#include <TList.h>
//#include <TFileIter.h>                                                                                                                                     
#include <TH2.h>
#include <TStyle.h>
#include <TString.h>
#include <TMath.h>
//#include <TClass.h>

#define nC 10  //nCuts in the analysis. make plots after each cut

class analyzer_higgs : public TSelector {
  ofstream nout[nC], fout[nC], ffout, ncout;
  
  double dz(TCPrimaryVtx *primVtx, TVector3* trackVtx, TVector3* trackMom);
  double dxy(TCPrimaryVtx *primVtx, TVector3* trackVtx, TVector3* trackMom);
  double higgsMT(Float_t pTll, Float_t Mll, Float_t MET, TVector2 ZptXY, TVector2 METXY);
  double projectedMET(Float_t MET, Float_t MET_phi, TLorentzVector Lepton1, TLorentzVector Lepton2);
  double pileUpCorrectedMET(Float_t projMET, Int_t nVtx);

  float getNevents(string, TH1F* );
  void scaleAndColor(TString , Float_t , Float_t , Float_t , Int_t , Int_t );
//  void CountEvents(UInt_t nEvents[], UInt_t nEventsPassNoiseFilter[]);
  
public :
  TFile* histoFile;
  
  TH1F *mtZ[nC], *mt0[nC], *mt1[nC], *mt2[nC], *mt3[nC], *mt4[nC];
  TH1F *met0_phi[nC], *met0_et[nC], *met0_over_qt[nC];
  TH1F *met1_phi[nC], *met1_et[nC], *met1_over_qt[nC];
  TH1F *met2_phi[nC], *met2_et[nC], *met2_over_qt[nC];
  TH1F *met3_phi[nC], *met3_et[nC], *met3_over_qt[nC];
  TH1F *met4_phi[nC], *met4_et[nC], *met4_over_qt[nC], *met4_puSig[nC];
  TH1F *mu1_phi[nC], *mu1_eta[nC], *mu1_pt[nC];
  TH1F *mu2_phi[nC], *mu2_eta[nC], *mu2_pt[nC];
  TH1F *btag_hp[nC];
  TH1F *di_qt[nC], *di_mass[nC], *di_mass_EB[nC], *di_mass_EE[nC], *di_mass_EX[nC];
  TH1F *jet_N[nC], *jet_pt[nC];
  TH1F *jet_b_N[nC], *jet_b_pt[nC];

  TH1F *met2_dPhiLeadJet1[nC], *met2_dPhiLeadJet2[nC], *met2_dPhiClosJet1[nC], *met2_dPhiClosJet2[nC];

  TH2F *met0_et_ovQt[nC], *met1_et_ovQt[nC], *met2_et_ovQt[nC], *met3_et_ovQt[nC], *met4_et_ovQt[nC]; 
  TH2F *mtZ_met2[nC], *mt2_met2[nC];
  TH2F *mtZ_met3[nC], *mt2_met3[nC];

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
