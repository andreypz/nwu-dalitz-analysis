// $Id: higgsAnalyzer.h,v 1.3 2011/08/25 10:18:29 andrey

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
#include <TH2.h>
#include <TStyle.h>
#include <TString.h>
#include <TKey.h>
#include <TList.h>


#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile.h"

#include "../src/TCJet.h"
#include "../src/TCMET.h"
#include "../src/TCElectron.h"
#include "../src/TCMuon.h"
#include "../src/TCTau.h"
#include "../src/TCPhoton.h"
#include "../src/TCGenJet.h"
#include "../src/TCPrimaryVtx.h"
#include "../src/TCTrigger.h"
#include "../src/TCTriggerObject.h"

#define nC 20  //nCuts in the analysis. make plots after each cut

class higgsAnalyzer : public TSelector {

 private:

  TFile* histoFile;
  //Variables to Fill into histos. Have to be global
  
  Float_t qT, diEta, Mll, Mll_EB, Mll_EE, Mll_EX;
  Float_t MET, MET_phi, MET1, MET1_phi;
  Float_t puCorrMET;
  Float_t MT, MT1, MTZ, pTll;
  Float_t METqt, MET1qt;
  Int_t nVtx;
  Float_t nDofVtx1, nDofVtx2;
  Int_t nJets, nJetsB, nJetsB25, nJetsB30,;
  /* jetCount, muCountLoose, eleCountLoose;
  Float_t dPhiClos1,dPhiClos2, dPhiLead1, dPhiLead2;
  */
  
  ofstream nout[nC], fout[nC], ffout, ncout;
  TTree * cutTree;

  TH1F *mtZ[nC], *mt0[nC], *mt1[nC], *mt2[nC], *mt3[nC], *mt4[nC];
  TH1F *met0_phi[nC], *met0_et[nC], *met0_over_qt[nC];
  TH1F *met1_phi[nC], *met1_et[nC], *met1_over_qt[nC];
  TH1F *met2_phi[nC], *met2_et[nC], *met2_over_qt[nC];
  TH1F *met3_phi[nC], *met3_et[nC], *met3_over_qt[nC];
  TH1F *met4_phi[nC], *met4_et[nC], *met4_over_qt[nC], *met4_puSig[nC];
  TH1F *mu1_phi[nC], *mu1_eta[nC], *mu1_pt[nC];
  TH1F *mu2_phi[nC], *mu2_eta[nC], *mu2_pt[nC];
  TH1F *btag_hp[nC];
  TH1F *di_qt[nC], *di_eta[nC], *di_mass[nC], *di_mass_EB[nC], *di_mass_EE[nC], *di_mass_EX[nC];
  TH1F *jet_N[nC], *jet_dRlep1[nC], *jet_dRlep2[nC], *jet_pt[nC];
  TH1F *jet_b_N[nC], *jet_b_N25[nC],*jet_b_N30[nC], *jet_b_pt[nC];

  TH1F *met2_dPhiLeadJet1[nC], *met2_dPhiLeadJet2[nC], *met2_dPhiClosJet1[nC], *met2_dPhiClosJet2[nC];

  TH2F *met0_et_ovQt[nC], *met1_et_ovQt[nC], *met2_et_ovQt[nC], *met3_et_ovQt[nC], *met4_et_ovQt[nC];
  TH2F *mtZ_met2[nC], *mt2_met2[nC];
  TH2F *mtZ_met3[nC], *mt2_met3[nC];

  TH1F *vtx_nPV_raw[nC], *vtx_nPV_weight[nC], *vtx_ndof_1[nC], *vtx_ndof_2[nC];

		//miscellaneous
		TH1D*     h1_ptHat;
		TH1D*     h1_triggerStatus;
		TH1D*     h1_goodRuns;
		TH1D*     h1_acceptanceByCut;
		TH1D*     h1_acceptanceByCutRaw;
		TH1D*     h1_pvMult;
		TH1D*     h1_simVertexMult;
		TH1D*     h1_eventWeight;
		TH2D*     h2_nEventsByHMass;
		TProfile* p1_nVtcs;

		TH1D*     h1_bJetMultPreVeto;
		TH1D*     h1_bJetMultPostVeto;

		TH1D*     h1_Met;

	public :
		TTree          *fChain;   //!pointer to the analyzed TTree or TChain

		// Declaration of leaf types
		TClonesArray  *recoJets;
		TClonesArray  *recoMET;
		TClonesArray  *genJets;
		TClonesArray  *recoMuons;
		TClonesArray  *recoElectrons;
		TClonesArray  *recoPhotons;
		TClonesArray  *primaryVtx;
		int           eventNumber;
		int           runNumber;
		int           lumiSection;
		int           nPUVertices;
		double        ptHat;
		float         rhoFactor;
		unsigned int  triggerStatus;
		int           hltPrescale[64];
		bool          isRealData, isNoiseHcal, isDeadEcalCluster, isScraping, isCSCTightHalo;

		// List of branches
		TBranch        *b_recoJets;   //!
		TBranch        *b_recoMET;   //!
		TBranch        *b_genJets;   //!
		TBranch        *b_recoMuons;
		TBranch        *b_recoElectrons;
		TBranch        *b_recoPhotons;
		TBranch        *b_primaryVtx;   //!
		TBranch        *b_eventNumber;   //!
		TBranch        *b_runNumber;   //!
		TBranch        *b_lumiSection;   //!
		TBranch        *b_rhoFactor;   //!
		TBranch        *b_ptHat;   //!
		TBranch        *b_triggerStatus;   //!
		TBranch        *b_hltPrescale;   //!
		TBranch        *b_nPUVertices;   //!
		TBranch        *b_isRealData;   //!
		TBranch        *b_isNoiseHcal;   //!
		TBranch        *b_isDeadEcalCluster;   //!
		TBranch        *b_isScraping;   //!
		TBranch        *b_isCSCTightHalo;   //!

		//For counting events
		//        int          nEvents[16];
		//For counting weighted events
		float        nEventsWeighted[nC];
		
		UInt_t nEvents[nC];
		UInt_t nEventsPassNoiseFilter[6][nC]; //1-Hcal, 2-Ecal, 3- Scraping, 4-CSCTight, 5-CSCLoose,   0-passed all
		
/*
		//MET corrections
		struct metCorrs {
			float Sigma0;
			float SigmaPU;
			float Shift;
		} metCorr;

		std::map<string, metCorrs> metCorrMap;
*/

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

		virtual float   Dz(TVector3 objVtx, TLorentzVector objP4, TVector3 vtx);
		virtual float   Dxy(TVector3 objVtx, TLorentzVector objP4, TVector3 vtx);
		virtual void    PostSelectionYieldCounter(float nEventsPS, TLorentzVector metP4, TLorentzVector zP4, float dPhiJetMet, float evtWeight);
		//virtual void    MetPlusZPlots(TLorentzVector metP4, TLorentzVector ZP4, float evtWeight);
		//virtual void    MetPlusLeptonPlots(TLorentzVector metP4, TLorentzVector p1, TLorentzVector p2, float evtWeight);
		//virtual void    LeptonBasicPlots(TLorentzVector p1, TLorentzVector p2, float evtWeight);
		//virtual void    DileptonBasicPlots(TLorentzVector ZP4, float evtWeight);
		virtual bool    CosmicMuonFilter(TCMuon , TCMuon );
		virtual float   CalculateTransMass(TLorentzVector p1, TLorentzVector p2);
		virtual float   CalculateTransMassAlt(TLorentzVector p1, TLorentzVector p2);
		virtual float   DeltaPhiJetMET(TLorentzVector , std::vector<TLorentzVector> ); 
		virtual float   GetEventWeight(int, TLorentzVector, TLorentzVector);
		virtual float   PUCorrectedMET(float, int, std::string);
		virtual float   GetPhotonMass();

		virtual TLorentzVector GetReducedMET(TLorentzVector sumJet, TLorentzVector lep1, TLorentzVector lep2, TLorentzVector metP4, int version);



		float getNevents(string, TH1F* );
		void scaleAndColor(TString , Float_t , Float_t , Float_t , Int_t , Int_t );

		//float weightZZ(Float_t );
		//float puReweighting(Int_t );

		void FillHistos(Int_t, Double_t);
		void FillHistosNoise(Int_t, Double_t);
		void CountEvents(Int_t);
		void PrintOut(Int_t);
		void PrintOutNoisy(Int_t);


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
	recoMET = 0;
	genJets = 0;
	primaryVtx = 0;
	recoMuons = 0;
	recoElectrons = 0;
	recoPhotons = 0;
	// Set branch addresses and branch pointers
	if (!tree) return;
	fChain = tree;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("recoJets", &recoJets, &b_recoJets);
	fChain->SetBranchAddress("recoElectrons", &recoElectrons, &b_recoElectrons);
	fChain->SetBranchAddress("recoMuons", &recoMuons, &b_recoMuons);
	fChain->SetBranchAddress("recoPhotons", &recoPhotons, &b_recoPhotons);
	fChain->SetBranchAddress("recoMET", &recoMET, &b_recoMET);
	fChain->SetBranchAddress("genJets", &genJets, &b_genJets);
	fChain->SetBranchAddress("primaryVtx", &primaryVtx, &b_primaryVtx);

	fChain->SetBranchAddress("isRealData", &isRealData, &b_isRealData);
	fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
	fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
	fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
	fChain->SetBranchAddress("rhoFactor", &rhoFactor, &b_rhoFactor);
	fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat);
	fChain->SetBranchAddress("nPUVertices", &nPUVertices, &b_nPUVertices);
	fChain->SetBranchAddress("triggerStatus", &triggerStatus, &b_triggerStatus);
	fChain->SetBranchAddress("hltPrescale",&hltPrescale, &b_hltPrescale);

	fChain->SetBranchAddress("isNoiseHcal", &isNoiseHcal, &b_isNoiseHcal);
	fChain->SetBranchAddress("isDeadEcalCluster", &isDeadEcalCluster, &b_isDeadEcalCluster);
	fChain->SetBranchAddress("isScraping", &isScraping, &b_isScraping);
	fChain->SetBranchAddress("isCSCTightHalo", &isCSCTightHalo, &b_isCSCTightHalo);

	// ADD THESE!
	//eventTree->Branch("recoTaus",&recoTaus, 6400, 0);
	//eventTree->Branch("beamSpot", &beamSpot, 6400, 0);
	//eventTree->Branch("bunchCross",&bunchCross, "bunchCross/i");
	//eventTree->Branch("qScale", &qScale, "qScale/F");
	//eventTree->Branch("evtWeight", &evtWeight, "evtWeight/F");

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
