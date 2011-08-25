#define higgsAnalyzer_cxx

#include "higgsAnalyzer.h"

using namespace std;

/////////////////////////////
//Specify parameters here. //
/////////////////////////////

string  selection      = "muon";
int     JC_LVL         = 3;
int     trigger[]      = {4,8};
string  suffix         = "DATA";


UInt_t verboseLvl = 1;

/////////////////
//Analysis cuts//
/////////////////

float   jetPtCut[]     = {20.};
float   muPtCut[]      = {20., 10.};
float   elePtCut[]     = {20., 10.};
float   zMassCut[]     = {76.2, 106.2};
float   metCut[]       = {83., 9999.};
float   mtCut[]        = {242., 320.};
float   metPlusQtCut[] = {0., 9999.};
float   metByQtCut[]   = {0., 9999.};
float   bJetPtCut      = 20.;
float   qtCut          = 25.;
float   dPhiJetMETCut  = 0.28;
int     nJetsCut[]     = {0, 99};
Float_t cut_vz = 24, cut_vd0 = 2, cut_vndof = 4;  //PV filter cuts


// Cuts for mass points in the PAS.

float dPhiMinCut[] = {0.62, 0.28, 0.14, TMath::Pi(), TMath::Pi(), TMath::Pi(), TMath::Pi(), TMath::Pi()};
float metMinCut[]  = {69., 83., 97., 112., 126., 141., 155., 170.};
float mtMinCut[]   = {216., 242., 267., 292., 315., 336., 357., 377.};
float mtMaxCut[]   = {272., 320., 386., 471., 540., 600., 660., 720.};

///////////////////////////
//Resources for weighting//
///////////////////////////

TFile reweightFile = TFile("../data/Reweight.root", "OPEN");
TH1D  *h_eGammaReweight  = (TH1D*)reweightFile.Get("h1_eGammaWeight");
TH1D  *h_muGammaReweight = (TH1D*)reweightFile.Get("h1_muGammaWeight");
TH1D  *h_puReweight      = (TH1D*)reweightFile.Get("h1_puReweightHist");

//metCorrMap["pfMet"] = (metCorr(1., 1., 0.));

// Do something about these: should just have one sort condition function
bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());} 
bool MuonSortCondition(const TCMuon& m1, const TCMuon& m2) {return (m1.Pt() > m2.Pt());}
bool ElectronSortCondition(const TCElectron& e1, const TCElectron& e2) {return (e1.Pt() > e2.Pt());}

void higgsAnalyzer::Begin(TTree * /*tree*/) 
{
    TString option = GetOption();
    TH1::SetDefaultSumw2(kTRUE);

    histoFile = new TFile("higgsHistograms.root", "RECREATE");
    histoFile->mkdir("Misc", "Misc");
    histoFile->mkdir("Lepton", "Lepton");
    histoFile->mkdir("Jet", "Jet");
    histoFile->mkdir("MET", "MET");
    histoFile->mkdir("MET+Lepton", "MET+Lepton");
    histoFile->mkdir("2D", "2D");
    histoFile->mkdir("TESTS", "TESTS");

    histoFile->mkdir("Andrey", "Andrey");

    histoFile->cd("Andrey");
    for(Int_t n=0; n<nC; n++)
      met0_et[n]       = new TH1F(Form("met0_et_%i",n), "met0_et", 40, 0,400);

    for(Int_t n=5; n<nC; n++)
      {
	//met0_et_ovQt[n]  = new TH2F(Form("met0_et_ovQt_%i",n), "pfMET vs Met/qt", 40, 0,400, 40, 0,4);
	//met1_et_ovQt[n]  = new TH2F(Form("met1_et_ovQt_%i",n), "MET1 vs Met/qt", 40, 0,400, 40, 0,4);
	met2_et_ovQt[n]  = new TH2F(Form("met2_et_ovQt_%i",n), "pfMet noise vs Met/qt", 40, 0,400, 40, 0,4);
	met3_et_ovQt[n]  = new TH2F(Form("met3_et_ovQt_%i",n), "projMET vs Met/qt", 40, 0,400, 40, 0,4);
	met4_et_ovQt[n]  = new TH2F(Form("met4_et_ovQt_%i",n), "puCorrMET vs Met/qt", 40, 0,400, 40, 0,4);

	met1_phi[n]     = new TH1F(Form("met1_phi_%i",n), "met1_phi", 40, -TMath::Pi(), TMath::Pi());
	met1_et[n]      = new TH1F(Form("met1_et_%i",n), "met1_et", 40, 0,400);
	met1_over_qt[n] = new TH1F(Form("met1_over_qt_%i",n), "met1_over_qt", 40, 0,4);

	met2_phi[n]     = new TH1F(Form("met2_phi_%i",n), "met2_phi", 40, -TMath::Pi(), TMath::Pi());
	met2_et[n]      = new TH1F(Form("met2_et_%i",n), "met2_et", 40, 0,400);
	met2_over_qt[n] = new TH1F(Form("met2_over_qt_%i",n), "met2_over_qt", 40, 0,4);

	met3_phi[n]     = new TH1F(Form("met3_phi_%i",n), "met3_phi", 40, -TMath::Pi(), TMath::Pi());
	met3_et[n]      = new TH1F(Form("met3_et_%i",n), "met3_et", 40, 0,400);
	met3_over_qt[n] = new TH1F(Form("met3_over_qt_%i",n), "met2_over_qt", 40, 0,4);

	met4_phi[n]     = new TH1F(Form("met4_phi_%i",n), "met4_phi", 40, -TMath::Pi(), TMath::Pi());
	met4_et[n]      = new TH1F(Form("met4_et_%i",n), "met4_et", 42, -20,400);
	//met4_over_qt[n] = new TH1F(Form("met4_over_qt_%i",n), "met4_over_qt", 42, -0.2,4);
	//	met4_puSig[n]   = new TH1F(Form("met4_sig_%i",n), "met4_sig", 42, -0.4,8);

	//met2_dPhiLeadJet1[n] =  new TH1F(Form("met2_dPhiLeadJet1_%i",n), "delta phi to lead jets", 40, 0, TMath::Pi());
	//met2_dPhiLeadJet2[n] =  new TH1F(Form("met2_dPhiLeadJet2_%i",n), "delta phi to lead jets", 40, 0, TMath::Pi());
	//met2_dPhiClosJet1[n] =  new TH1F(Form("met2_dPhiClosJet1_%i",n), "delta phi to closest jets", 40, 0, TMath::Pi());
	//met2_dPhiClosJet2[n] =  new TH1F(Form("met2_dPhiClosJet2_%i",n), "delta phi to closest jets", 40, 0, TMath::Pi());


	//mu1_phi[n] = new TH1F(Form("mu1_phi_%i",n), "mu1_phi", 40, -TMath::Pi(), TMath::Pi());
	//mu1_eta[n] = new TH1F(Form("mu1_eta_%i",n), "mu1_eta", 40, -2.6, 2.6);
	//mu1_pt[n]  = new TH1F(Form("mu1_pt_%i",n), "mu1_pt", 40, 0, 200);
	//mu2_phi[n] = new TH1F(Form("mu2_phi_%i",n), "mu2_phi", 40, -TMath::Pi(), TMath::Pi());
	//mu2_eta[n] = new TH1F(Form("mu2_eta_%i",n), "mu2_eta", 40, -2.6, 2.6);
	//mu2_pt[n]  = new TH1F(Form("mu2_pt_%i",n), "mu2_pt", 40, 0, 200);
	btag_hp[n] = new TH1F(Form("btag_hp_%i",n), "btag_hp", 40, -10, 10);

	//mtZ[n]     = new TH1F(Form("mtZ_%i",n), "MTZ pf", 50, 100, 600);
	mt0[n]     = new TH1F(Form("mt0_%i",n), "MT pf", 50, 100, 600);
	mt1[n]     = new TH1F(Form("mt1_%i",n), "MT type1", 50, 100, 600);
	mt2[n]     = new TH1F(Form("mt2_%i",n), "MT noise", 50, 100, 600);
	mt3[n]     = new TH1F(Form("mt3_%i",n), "MT proj", 50, 100, 600);
	mt4[n]     = new TH1F(Form("mt4_%i",n), "MT pu corr", 50, 100, 600);

	//mtZ_met2[n]   = new TH2F(Form("mtZ_met2_%i",n), "MTZ pf vs pfMet noise", 50, 100, 600, 40, 0,400);
	mt2_met2[n]   = new TH2F(Form("mt2_met2_%i",n), "MT pf vs pfMet noise", 50, 100, 600, 40, 0,400);

	//mtZ_met3[n]   = new TH2F(Form("mtZ_met3_%i",n), "MTZ pf vs projMet noise", 50, 100, 600, 40, 0,400);
	//mt2_met3[n]   = new TH2F(Form("mt2_met3_%i",n), "MT pf vs projfMet noise", 50, 100, 600, 40, 0,400);

	di_qt[n]      = new TH1F(Form("di_qt_%i",n), "di-lepton qt", 40, 0,400);
	di_eta[n]     = new TH1F(Form("di_eta_%i",n), "di-lepton Eta", 40, -2.6,2.6);
	di_mass[n]    = new TH1F(Form("di_mass_%i",n), "di-lepton Mass", 50, 70,120);
	di_mass_EB[n] = new TH1F(Form("di_mass_EB_%i",n), "di-lepton Mass in Ecal Barrel", 50, 70,120);
	di_mass_EE[n] = new TH1F(Form("di_mass_EE_%i",n), "di-lepton Mass in Ecal Endcap", 50, 70,120);
	di_mass_EX[n] = new TH1F(Form("di_mass_EX_%i",n), "di-lepton Mass in Ecal E/B mix", 50, 70,120);

	jet_N[n]      = new TH1F(Form("jet_N_%i",n), "Number of jets", 20,0,20);
	jet_dRlep1[n] = new TH1F(Form("jet_dRlep1_%i",n), "dR(jet, lepton1)", 50, 0,1);
	jet_dRlep2[n] = new TH1F(Form("jet_dRlep2_%i",n), "dR(jet, lepton2)", 50, 0,1);
	jet_b_N[n]    = new TH1F(Form("jet_b_N_%i",n), "Number of b-jets", 10,0,10);

	jet_pt[n]   = new TH1F(Form("jet_pt_%i",n), "Pt of jets, eta<2.4", 40,0,200);
	jet_b_pt[n] = new TH1F(Form("jet_b_pt_%i",n), "Pt of b-jets, eta<2.4", 40,0,200);


	vtx_nPV_raw[n]    = new  TH1F(Form("vtx_nPV_raw_%i",n), "Number of PVs, raw", 0,0,20);
	vtx_nPV_weight[n] = new  TH1F(Form("vtx_nPV_weight_%i",n), "Number of PVs, reweighted", 0,0,20);
	vtx_ndof_1[n]     = new  TH1F(Form("vtx_ndof_1_%i",n), "Ndof of first PV", 0,0,100);
	vtx_ndof_2[n]     = new  TH1F(Form("vtx_ndof_2_%i",n), "Ndof of second PV", 0,0,100);

      }
    ffout.open("./events_printout_DATA_final.txt",ofstream::out);
    ncout.open("./counts_for_tex_DATA.txt",ofstream::out);

    ffout.precision(4); ffout.setf(ios::fixed, ios::floatfield);
    for(Int_t i=6; i<nC; i++)
      {
	fout[i].open(Form("./events_printout_DATA_%i.txt",i),ofstream::out);
	fout[i].precision(3); fout[i].setf(ios::fixed, ios::floatfield);
	fout[i]<<
	  "******\t*********\t*****\t********\t************\t*********\t********\t*****\t*\n"<<
	  "* Run \t* Event  \t* LS \t*  m(ll)\t*   MT(llvv)\t* PFMET  \t* PT(ll)\t* rho\t*\n"<<
	  "******\t*********\t*****\t********\t************\t*********\t********\t*****\t*\n";
	nout[i].open(Form("./events_filtered_DATA_%i.txt",i),ofstream::out);
	nout[i].precision(3); fout[i].setf(ios::fixed, ios::floatfield);
	nout[i]<<
	  "******\t*********\t*********\t******************************************\n"<<
	  "* Run \t* Event  \t*  LS   \t*  PFMET * Hcal * Ecal * scrap * CSCTight * \n"<<
	  "******\t*********\t*********\t******************************************\n";
      }

    histoFile->cd("Misc");
    h1_ptHat                        = new TH1D("h1_ptHat_DATA", "ptHat", 37, 15.0, 200.0);
    h1_triggerStatus                = new TH1D("h1_triggerStatus_DATA", "Triggers", 32, 0.5, 32.5);
    h1_goodRuns                     = new TH1D("h1_goodRuns_DATA", "Number of events passing selection cuts by run;Run number; nEvents", 9200, 160000., 175000.);
    h1_acceptanceByCut              = new TH1D("h1_acceptanceByCut_DATA", "Weighted number of events passing cuts by cut; cut; N_{evts}", 16, -0.5, 15.5);
    h1_acceptanceByCutRaw           = new TH1D("h1_acceptanceByCutRaw_DATA", "Raw number of events passing cuts; cut; N_{evts}", 16, -0.5, 15.5);
    h1_pvMult                       = new TH1D("h1_pvMult_DATA", "Multiplicity of PVs", 25, 0.5, 25.5);
    h1_simVertexMult                = new TH1D("h1_simVertexMult_DATA", "Multiplicity of simulated vertices", 25, -0.5, 24.5);
    h1_eventWeight                  = new TH1D("h1_eventWeight_DATA", "event weight", 100, 0., 2.);
    h2_nEventsByHMass               = new TH2D("h2_nEventsByHMass_DATA", "", 10, 0., 10., 10, 0., 10.);
    p1_nVtcs                        = new TProfile("p1_nVtcs", "Average number of vertices per run; Run Number; nVertices", 8700.0, 135000.0, 144200.0, 0.0, 6.0);

    histoFile->cd("Lepton");
    h1_leadLeptonPt                 = new TH1D("h1_leadLeptonPt_DATA", "p_{T} leading lepton;p_{T};N_{evts}", 48, 10., 250.);
    h1_leadLeptonEta                = new TH1D("h1_leadLeptonEta_DATA", "#eta leading lepton;#eta;N_{evts}", 25, -2.5, 2.5);
    h1_leadLeptonPhi                = new TH1D("h1_leadLeptonPhi_DATA", "#phi leading lepton;#phi;N_{evts}", 18, -TMath::Pi(), TMath::Pi());
    h1_trailingLeptonPt             = new TH1D("h1_trailingLeptonPt_DATA", "p_{T} trailing lepton;p_{T};N_{evts}", 38, 10., 200.);
    h1_trailingLeptonEta            = new TH1D("h1_trailingLeptonEta_DATA", "#eta trailing lepton;#eta;N_{evts}", 25, -2.5, 2.5);
    h1_trailingLeptonPhi            = new TH1D("h1_trailingLeptonPhi_DATA", "#phi trailing lepton;#phi;N_{evts}", 18, -TMath::Pi(), TMath::Pi());

    h1_diLeptonTransMass            = new TH1D("h1_diLeptonTransMass_DATA", "M_{T,ll};M_{T,ll};N_{evts}", 100, 55., 255.);
    h1_diLeptonMass                 = new TH1D("h1_diLeptonMass_DATA", "M_{ll}; M_{ll};N_{evts}", 40, 70., 110.);
    h1_diLeptonQt                   = new TH1D("h1_diLeptonQt_DATA", "q_{T};Q_{T};N_{evts}", 50, 0., 500.);
    h1_diLeptonPtRatio              = new TH1D("h1_diLeptonPtRatio_DATA", "dilepton p_{T} ratio;p_{T,2}/p_{T,1};N_{evts}", 25, 0., 1.);
    h1_diLeptonDeltaEta             = new TH1D("h1_diLeptonDeltaEta_DATA", "dilepton #Delta#eta;#Delta#eta;N_{evts}", 25, 0., 2.5);
    h1_diLeptonDeltaPhi             = new TH1D("h1_diLeptonDeltaPhi_DATA", "dilepton #Delta#phi;#Delta#phi;N_{evts}", 18, 0., TMath::Pi());
    h1_diLeptonDeltaR               = new TH1D("h1_diLeptonDeltaR_DATA", "dilepton #Delta R;#Delta R;N_{evts}", 18, 0., 4.5);  	

    histoFile->cd("Jet");
    h1_leadJetPt                  = new TH1D("h1_leadJetPt_DATA", "p_{T} of lead jet;p_{T};N_{evts}", 26, 10., 250.);
    h1_leadJetEta                 = new TH1D("h1_leadJetEta_DATA", "#eta of lead jet;#eta;N_{evts}", 25, -2.5, 2.5);
    h1_leadJetPhi                 = new TH1D("h1_leadJetPhi_DATA", "#phi of lead jet;#phi;N_{evts}", 18, -TMath::Pi(), TMath::Pi());
    h1_tailJetPt                  = new TH1D("h1_tailJetPt_DATA", "p_{T} of tail jet;p_{T};N_{evts}", 26, 10., 250.);
    h1_tailJetEta                 = new TH1D("h1_tailJetEta_DATA", "#eta of tail jet;#eta;N_{evts}", 25, -2.5, 2.5);
    h1_tailJetPhi                 = new TH1D("h1_tailJetPhi_DATA", "#phi of tail jet;#phi;N_{evts}", 18, -TMath::Pi(), TMath::Pi());
    h1_nearestJetEta              = new TH1D("h1_nearestJetEta_DATA", "#eta of jet_{nearest};#eta;N_{evts}", 25, -2.5, 2.5);
    h1_jetMult                    = new TH1D("h1_jetMult_DATA", "Multiplicity of jets;N_{jets};N_{evts}", 11, -0.5, 10.5);
    h1_leadBJetPt                 = new TH1D("h1_leadBJetPt_DATA", "p_{T} of lead b-jet;p_{T};N_{evts}", 36, 10., 190.);
    h1_leadBJetEta                = new TH1D("h1_leadBJetEta_DATA", "#eta of lead b-jet;#eta;N_{evts}", 25, -2.5, 2.5);
    h1_leadBJetPhi                = new TH1D("h1_leadBJetPhi_DATA", "#phi of lead b-jet;#phi;N_{evts}", 18, -TMath::Pi(), TMath::Pi());
    h1_bJetMultPostVeto           = new TH1D("h1_bJetMultPostVeto_DATA", "Multiplicity of b-jets (post-veto);N_{jets};N_{evts}", 11, -0.5, 10.5);
    h1_bJetMultPreVeto            = new TH1D("h1_bJetMultPreVeto_DATA", "Multiplicity of b-jets (pre-veto);N_{jets};N_{evts}", 11, -0.5, 10.5);

    histoFile->cd("MET");
    h1_Met                        = new TH1D("h1_Met_DATA", "MET;MET;N_{evts}", 50, 0., 250.);
    h1_MetPhi                     = new TH1D("h1_MetPhi_DATA", "#phi MET;#phi;N_{evts}", 18, -TMath::Pi(), TMath::Pi());
    h1_MetSumEt                   = new TH1D("h1_MetSumEt_DATA", "#Sigma E_{T} of MET;#Sigma E_{T};N_{evts}", 75, 0., 1500.);
    h1_ReducedMET                 = new TH1D("h1_ReducedMET_DATA", " Reduced MET;ReducedMET;N_{evts}", 50, 0., 250.);
    h1_ReducedMETTransverse       = new TH1D("h1_ReducedMETTransverse_DATA", " Reduced MET (Transverse);TransRedMET;N_{evts}", 40, -50., 50.);
    h1_ReducedMETLongitudinal     = new TH1D("h1_ReducedMETLongitudinal_DATA", " Reduced MET (longitudinal);TransRedMET;N_{evts}", 40, -50., 50.);
    h1_MetNearestJetDeltaPhi      = new TH1D("h1_MetNearestJetDeltaPhi_DATA", "#Delta#phi(j_{nearest}, MET);#Delta#phi;N_{evts}", 18, 0., TMath::Pi());
    h1_ProjectedMet               = new TH1D("h1_ProjectedMet_DATA", "ProjMET; ProjMET; N_{evts}", 50, 0., 250.);

    histoFile->cd("MET+Lepton");
    h1_MetOverQt                  = new TH1D("h1_MetOverQt_DATA", "MET/Q_{T,ll};MET/q_{T};N_{evts}", 45, 0., 9.);
    h1_MetPlusQtMagnitude         = new TH1D("h1_MetPlusQtMagnitude_DATA", "|MET + q_{T}|;|MET + Q_{T,ll}|;N_{evts}", 60, 0., 300.);
    h1_MetQtDeltaPhi              = new TH1D("h1_MetQtDeltaPhi_DATA", "#Delta#phi(q_{T}, MET);#Delta#phi;N_{evts}", 18, 0., TMath::Pi());
    h1_MetLeadLeptonDeltaPhi      = new TH1D("h1_MetLeadLeptonDeltaPhi_DATA", "#Delta#phi(l1, MET);#Delta#phi;N_{evts}", 18, 0., TMath::Pi());
    h1_MetTrailingLeptonDeltaPhi  = new TH1D("h1_MetTrailingLeptonDeltaPhi_DATA", "#Delta#phi(l2, MET);#Delta#phi;N_{evts}", 18, 0., TMath::Pi());
    h1_ProjMetByQt                = new TH1D("h1_ProjMetByQt_DATA", "Longitudinal MET by q_{T};MET_{longitudinal};N_{evts}", 52, -180., 80.);
    h1_OrthoMetByQt               = new TH1D("h1_OrthoMetByQt_DATA", "Transverse MET by q_{T};MET_{Transverse};N_{evts}", 64, -100., 80.);
    h1_ZPlusJetMHT                = new TH1D("h1_ZPlusJetMHT_DATA", "MET of jets + dilepton;MET;N_{evts}", 80, 0., 400.);
    h1_MetDileptonMT              = new TH1D("h1_MetDileptonMT_DATA", "M_{T};M_{T};N_{evts}", 90, 50., 500.);
    h1_MetLeadLeptonMT            = new TH1D("h1_MetLeadLeptonMT_DATA", "M_{T,l1};M_{T,l1};N_{evts}", 50, 0., 250.);
    h1_MetTrailingLeptonMT        = new TH1D("h1_MetTrailingLeptonMT_DATA", "M_{T,l2};M_{T,l2};N_{evts}", 50, 0., 250.);

    histoFile->cd("2D");
    h2_MetByMetOverQt             = new TH2D("h2_MetByMetOverQt_DATA", "MET/q_{T} vs. MET; MET; MET/q_{T}", 40, 0., 200., 40, 0., 4.);
    h2_MetByMetPlusQtMag          = new TH2D("h2_MetByMetPlusQtMag_DATA", "|MET + q_{T}| vs. MET; MET; MET/q_{T}", 40, 0., 200., 60, 0., 300.);
    h2_JetMultVsPVMult            = new TH2D("h2_JetMultVsPVMult_DATA", "Jet multiplicity vs. N_{PV};N_{PV};N_{jets}", 20, 0.5, 20.5, 11, -0.5, 10.5);
    h2_MetLeptonMT                = new TH2D("h2_MetLeptonMT_DATA", "MT(MET,lepton);MT(MET, l1);MT(MET, l2)", 80, 0., 400., 80, 0., 400.);

    histoFile->cd("TESTS");
    h1_TESTS[0]                   = new TH1D("h1_TEST1", "Jet multiplicity;N_{evts};N_{jets}", 10, -0.5, 9.5);
    h1_TESTS[1]                   = new TH1D("h1_TEST2", "no PU, associated", 10, -0.5, 9.5);
    h1_TESTS[2]                   = new TH1D("h1_TEST3", "PU, not associated", 10, -0.5, 9.5);
    h1_TESTS[3]                   = new TH1D("h1_TEST4", "no PU, not associated", 10, -0.5, 9.5);
    h1_TESTS[4]                   = new TH1D("h1_TEST5", "MET (corrected);MET;N_{evts}", 50, 0., 250.);
    h1_TESTS[6]                   = new TH1D("h1_TEST7", "sum pt fraction;#Sigma p_{T};N_{jets}", 30, 0., 1.2);
    h1_TESTS[7]                   = new TH1D("h1_TEST8", "jet-vertex index", 10, -0.5, 9.5);
    h1_TESTS[8]                   = new TH1D("h1_TEST9", "jet d_{xy}", 200, -1., 1.);
    h1_TESTS[9]                   = new TH1D("h1_TEST10", "jet d_{z}", 200, -1., 1.);

    h2_TESTS[0]                   = new TH2D("h2_TEST1", "jet d_{xy} vs. sumPt cut", 10, 0., 1., 200, -1., 1.);
    h2_TESTS[1]                   = new TH2D("h2_TEST2", "jet d_{z} vs. sumPt cut", 10, 0., 1., 200, -1., 1.);

}

bool higgsAnalyzer::Process(Long64_t entry)
{  
    GetEntry(entry);
    CountEvents(0);
    ++nEventsWeighted[0];
    MET = 0;
    FillHistosNoise(0, 1);

    if (nEvents[0] % (int)1e5 == 0) cout<<nEvents[13]<<" events passed of "<<nEvents[0]<<" checked!"<<endl;

    if (selection == "gamma"   && (eventNumber % 3) != 0) return kTRUE;
    if (selection == "eGamma"  && (eventNumber % 3) != 1) return kTRUE;
    if (selection == "muGamma" && (eventNumber % 3) != 2) return kTRUE;

    //////////////////
    //Trigger status//
    //////////////////

    for(int i = 0; i < 32; ++i) {
        unsigned int iHLT = 0x0; 
        iHLT = 0x01 << i;  
        if ((triggerStatus & iHLT) == iHLT) h1_triggerStatus->Fill(i+1);  
    } 

    bool triggerPass = false; 
    if (trigger[0] != 0) { 
        //unsigned int iHLT = 0x0; 
        //for (int i = 0; i < (int)(sizeof(trigger)/sizeof(int)); ++i) iHLT |= 0x01 << (trigger[i] - 1);  
        //if ((triggerStatus & iHLT) == iHLT) triggerPass = true; 

        for (int i = 0; i < (int)(sizeof(trigger)/sizeof(int)); ++i) {
            unsigned int iHLT = 0x01 << (trigger[i] - 1);  
            if ((selection == "muon" || selection == "electron")  && hltPrescale[i] != 1) continue;
            if ((triggerStatus & iHLT) == iHLT) triggerPass = true; 
        }
    } else {triggerPass = true;}
    

    if (!triggerPass) return kTRUE;

    CountEvents(1);
    ++nEventsWeighted[1];
    FillHistosNoise(1, 1);

    ////////////////////////////
    //Check the event vertices//
    ////////////////////////////


    if (!isRealData) h1_simVertexMult->Fill(nPUVertices);
    // if (primaryVtx->GetSize() == 0) return kTRUE;

    //---------------------------------------------------
    //--- Primary vertex filter -------------------------
    //-------------------------------------------------
  
    Int_t PVind[2] = {-1, -1}; //indecies for two best PVs
    Bool_t vertexFilter = kFALSE;
    nVtx = 0;
    for (Int_t vv=0; vv<primaryVtx->GetSize(); vv++)
      {
	TCPrimaryVtx* pVtx = (TCPrimaryVtx*) primaryVtx->At(vv);
	//if(!pVtx->isValid()) cout<<"PV is not valid"<<endl;
	//if(pVtx->isFake()) cout<<"PV is fake"<<endl;

	if( 
	   !pVtx->IsFake() && 
	    pVtx->NDof() > cut_vndof                      //4
	    && fabs(pVtx->Position().z()) <= cut_vz      //15
	    && fabs(pVtx->Position().Perp()) <= cut_vd0   //Not Mag()!; 2
	    )
	  {
	    vertexFilter = kTRUE;

	    if (PVind[0]==-1) PVind[0] = vv; //Set first main PV
	    if (PVind[1]==-1 && PVind[0]!=vv) PVind[1] = vv; //Second PV
	    nVtx++;
	  }
      }

    TCPrimaryVtx *mainPrimaryVertex = 0, *secondPrimaryVertex = 0;
    if (PVind[0]!=-1) mainPrimaryVertex   = (TCPrimaryVtx*)(primaryVtx->At(PVind[0]));
    if (PVind[1]!=-1) secondPrimaryVertex = (TCPrimaryVtx*)(primaryVtx->At(PVind[1]));  

    // Apply the PV filters here! -------------
    if (!vertexFilter) return kTRUE;  

    TVector3* pvPosition = new TVector3();
    *pvPosition = mainPrimaryVertex->Position();


    //////////////////
    // Data quality //
    //////////////////

    //if (isNoiseHcal || isDeadEcalCluster || isScraping || isCSCTightHalo) return kTRUE;// <-- Not yet in ntuples

    ///////////////
    // electrons //
    ///////////////

    vector<TLorentzVector> extraLeptons;
    vector<TCElectron> electrons;
    //int eleCount = 0;

    for (int i = 0; i <  recoElectrons->GetSize(); ++i) {
        TCElectron* thisElec = (TCElectron*) recoElectrons->At(i);    

        if (fabs(thisElec->Eta()) > 2.5 || !thisElec->PassID(95)) continue;

        float eleISOendcap = (thisElec->TrkIso() + thisElec->EmIso() + thisElec->HadIso() - rhoFactor*TMath::Pi()*0.09)/thisElec->Pt(); 
        float eleISObarrel = (thisElec->TrkIso() + TMath::Max(0.0, thisElec->EmIso()-1.0)+ thisElec->HadIso() - rhoFactor*TMath::Pi()*0.09)/thisElec->Pt(); 
        bool  elecPass = false;

        if (
                (fabs(thisElec->Eta()) < 1.442     
                 && eleISObarrel                        < 0.1
                 && thisElec->SigmaIetaIeta()           < 0.01  
                 && fabs(thisElec->DphiSuperCluster())  < 0.06  
                 && fabs(thisElec->DetaSuperCluster())  < 0.004 
                 && thisElec->HadOverEm()               < 0.04      
                ) ||
                (fabs(thisElec->Eta()) >  1.556  
                 && eleISOendcap                        < 0.1 
                 && thisElec->SigmaIetaIeta()           < 0.03  
                 && fabs(thisElec->DphiSuperCluster())  < 0.03  
                 && fabs(thisElec->DetaSuperCluster())  < 0.007 
                 && thisElec->HadOverEm()               < 0.15
                )) elecPass = true; 

        if (elecPass) { 
            if (
                    thisElec->Pt() > elePtCut[0]  
                    && thisElec->PassConversion(80)
                    && fabs(thisElec->Dxy(pvPosition)) < 0.02
                    && fabs(thisElec->Dz(pvPosition))  < 0.1
               ) electrons.push_back(*thisElec);			
        } else {
            if (
                    thisElec->Pt() > elePtCut[1]  
                    && thisElec->PassConversion(95)
                    && fabs(thisElec->Dxy(pvPosition)) < 0.02
                    && fabs(thisElec->Dz(pvPosition)) < 0.1
               ) extraLeptons.push_back(thisElec->P4());			
        }
    } 

    sort(electrons.begin(), electrons.end(), ElectronSortCondition);

    ///////////
    // muons //
    ///////////

    vector<TCMuon> muons;

    for (int i = 0; i < recoMuons->GetSize(); ++ i) {
        TCMuon* thisMuon = (TCMuon*) recoMuons->At(i);    

        if (!(fabs(thisMuon->Eta()) < 2.4	
                    && thisMuon->IsGLB() 
                    && thisMuon->IsTRK()
                    && thisMuon->PtError()/thisMuon->Pt() < 0.1
             )) continue;

        if (
                thisMuon->Pt() > muPtCut[0]
                && thisMuon->NumberOfValidMuonHits()    > 0
                && thisMuon->NumberOfValidTrackerHits() > 10
                && thisMuon->NumberOfValidPixelHits()   > 0
                && thisMuon->NumberOfMatches() > 1
                && thisMuon->NormalizedChi2()  < 10
                && fabs(thisMuon->Dxy(pvPosition)) < 0.02
                && fabs(thisMuon->Dz(pvPosition))  < 0.1 
                && (thisMuon->TrkIso() + thisMuon->HadIso())/thisMuon->Pt() < 0.15
                //&& (thisMuon->TrkIso() + thisMuon->HadIso() + thisMuon->EmIso() - rhoFactor*TMath::Pi()*0.09)/thisMuon->Pt() < 0.15
           ) {
            muons.push_back(*thisMuon);
        } else if (
                thisMuon->Pt() > muPtCut[1] 
                && thisMuon->Pt() <= muPtCut[0] 
                && thisMuon->NumberOfValidMuonHits()    > 0
                && thisMuon->NumberOfValidTrackerHits() > 2 
                && thisMuon->NumberOfValidPixelHits()   > 0
                && thisMuon->NumberOfMatches() > 0
                && thisMuon->NormalizedChi2()  < 100
                && fabs(thisMuon->Dxy(pvPosition)) < 0.02
                && fabs(thisMuon->Dz(pvPosition))  < 0.1
                //&& (thisMuon->TrkIso() + thisMuon->HadIso() + thisMuon->EmIso() - rhoFactor*TMath::Pi()*0.09)/thisMuon->Pt() < 0.15
                ) {
            extraLeptons.push_back(thisMuon->P4());
        }
    } 

    sort(muons.begin(), muons.end(), MuonSortCondition);
    sort(extraLeptons.begin(), extraLeptons.begin(), P4SortCondition);

    /////////////
    // photons //
    /////////////

    vector<TLorentzVector> photons; 
    if (selection == "eGamma" || selection == "muGamma" || selection == "gamma") {
        for (Int_t i = 0; i < recoPhotons->GetSize(); ++i)
        {
            TCPhoton* thisPhoton = (TCPhoton*) recoPhotons->At(i);

            if ( thisPhoton->EmIso()                < (4.2 + 0.006 *thisPhoton->Pt())
                    &&	thisPhoton->HadIso()        < (2.2 + 0.0025*thisPhoton->Pt())
                    &&	thisPhoton->TrkIso()        < (2.0 + 0.001 *thisPhoton->Pt())
                    &&	thisPhoton->HadOverEm()     < 0.05 
                    &&	thisPhoton->SigmaIEtaIEta() > 0.001
                    &&	thisPhoton->SigmaIEtaIEta() < 0.013
                    &&	thisPhoton->TrackVeto() == 0
                    &&	(thisPhoton->EtaSupercluster() < 1.442 || thisPhoton->EtaSupercluster() > 1.556) 
               ) photons.push_back(thisPhoton->P4());
        }
        sort(photons.begin(), photons.end(), P4SortCondition);
    }


    //////////
    // jets //
    //////////

    vector<TLorentzVector> jetP4, bJetP4, testJetP4;
    TLorentzVector sumJetP4;
    int jetCount = 0;
    float fwdJetSumPt = 0.;
    int fwdJetCount = 0;

    nJets=0; nJetsB=0;
    for (int i = 0; i < recoJets->GetSize(); ++i) {
      TCJet* thisJet = (TCJet*) recoJets->At(i);     
      if (fabs(thisJet->P4().Eta()) < 2.4) {
	
	// Prevent lepton overlap //
	bool leptonOverlap = false;
	for (int j = 0; j < (int)muons.size(); ++j) if (thisJet->P4().DeltaR(muons[j].P4()) < 0.4) leptonOverlap = true;
	for (int j = 0; j < (int)electrons.size(); ++j) if (thisJet->P4().DeltaR(electrons[j].P4()) < 0.4) leptonOverlap = true;
	
	if (selection == "eGamma" || selection == "muGamma") {
	  for (int j = 0; j < (int)photons.size(); ++j) if (thisJet->P4().DeltaR(photons[j]) < 0.4) leptonOverlap = true;
	}
	if (leptonOverlap) continue;
	if (
	    thisJet->NumConstit() > 1
	    //&& thisJet->NumChPart() > 0
	    && thisJet->ChHadFrac() > 0
	    && thisJet->ChEmFrac() < 1
	    && (thisJet->NeuHadFrac() + thisJet->NeuEmFrac()) < 0.9
	    && thisJet->P4(JC_LVL).Pt() > jetPtCut[0]
	    ) {
	  if (thisJet->BDiscrTCHE() > 2. && thisJet->P4(JC_LVL).Pt() > bJetPtCut){
	    bJetP4.push_back(thisJet->P4(JC_LVL));
	  }
	  
	  testJetP4.push_back(thisJet->P4(JC_LVL));
	  //h1_TESTS[7]->Fill(thisJet->VtxIndex());
	  
	  if (
	      thisJet->VtxNTracks() > 0
	      && thisJet->VtxSumPtFrac() > 0. 
	      && thisJet->VtxSumPtIndex() == 1
	      ) {
	    jetP4.push_back(thisJet->P4(JC_LVL));
	    sumJetP4 += thisJet->P4();
	  }
	  ++jetCount;
	}
      } else {
	if (
	    thisJet->P4(JC_LVL).Eta() < 5.0
	    && thisJet->P4(JC_LVL).Pt() > 15. //jetPtCut[0]
	    ) {
	  //jetP4.push_back(thisJet->P4(JC_LVL)); 
	  //sumJetP4 += thisJet->P4();
	  fwdJetSumPt += thisJet->Pt(JC_LVL);
	  ++fwdJetCount;
	  ++jetCount;
	}
      }
    }

    sort(jetP4.begin(), jetP4.end(), P4SortCondition);
    sort(bJetP4.begin(), bJetP4.end(), P4SortCondition);

    nJets  = jetP4.size();
    nJetsB = bJetP4.size();

    if (!isRealData) {
        if (nPUVertices == 0) {
            h1_TESTS[1]->Fill(jetP4.size());
            h1_TESTS[3]->Fill(testJetP4.size());
        } else {
            h1_TESTS[0]->Fill(jetP4.size());
            h1_TESTS[2]->Fill(testJetP4.size());
        }
    }


    /////////
    // MET //
    /////////

    TCMET* met = (TCMET*) recoMET;
    TLorentzVector metP4, reducedMetP4;
    //metP4.SetPtEtaPhiE(met->CorrectedMet(), 0, met->CorrectedPhi(), met->CorrectedMet());
    //metP4.SetPtEtaPhiE(puCorrectedMet, 0, met->Phi(), puCorrectedMet);
    metP4.SetPtEtaPhiE(met->Met(), 0, met->Phi(), met->Met());

    ////////////////////////
    // Analysis selection //
    ////////////////////////

    float evtWeight = 1;

    TLorentzVector ZP4;
    MET        = metP4.Pt();
    puCorrMET  = PUCorrectedMET(met->Met(), nVtx, "pfMet");
    ///////////////////////
    // Cosmics rejection //
    ///////////////////////

    if (muons.size() > 1) if (CosmicMuonFilter(muons[0], muons[1])) return kTRUE;
    CountEvents(2);
    nEventsWeighted[2] += evtWeight;
    FillHistosNoise(2, evtWeight);


    float deltaPhiJetMET = 2*TMath::Pi();
    if (jetP4.size() > 0) deltaPhiJetMET= DeltaPhiJetMET(metP4, jetP4);


    TLorentzVector Lepton1(0.,0.,0.,0.), Lepton2(0.,0.,0.,0.);
   if (selection == "electron") {
        /////////////////////
        // 2 good elctrons //
        /////////////////////

        if (electrons.size() < 2) return kTRUE;
        if (electrons[0].Charge() == electrons[1].Charge()) return kTRUE;

        ZP4          = electrons[0].P4() + electrons[1].P4();
        evtWeight    = GetEventWeight(nPUVertices, electrons[0].P4(), electrons[1].P4());
        reducedMetP4 = GetReducedMET(sumJetP4,  electrons[0].P4(), electrons[1].P4(), metP4, 2);

	Lepton1 = electrons[0].P4();
	Lepton2 = electrons[1].P4();

    } else if (selection == "muon") {

        ////////////////////////////////
        // 2 oppositely charged muons //
        ////////////////////////////////

        if (muons.size() < 2) return kTRUE;
        if (muons[0].Charge() == muons[1].Charge()) return kTRUE;

        ZP4          = muons[0].P4() + muons[1].P4();
        evtWeight    = GetEventWeight(nPUVertices, muons[0].P4(), muons[1].P4());
        reducedMetP4 = GetReducedMET(sumJetP4,  muons[0].P4(), muons[1].P4(), metP4, 2);

	Lepton1 = muons[0].P4();
	Lepton2 = muons[1].P4();

    } else if (selection == "eGamma" || selection == "muGamma" || selection == "gamma") {

        /////////////////////////////////////
        // Photons in place of Z candidate //
        /////////////////////////////////////

        if (photons.size() < 1) return kTRUE;
        ZP4.SetPtEtaPhiM(photons[0].Pt(), photons[0].Eta(), photons[0].Phi(), GetPhotonMass());
        evtWeight  = GetEventWeight(primaryVtx->GetSize(), ZP4, TLorentzVector(0, 0, 0, 0));
        reducedMetP4 = GetReducedMET(sumJetP4,  photons[0], TLorentzVector(0, 0, 0, 0), metP4, 2);
    } else {
        return kTRUE;
    }


    CountEvents(3);
    nEventsWeighted[3] += evtWeight;
    FillHistosNoise(3, evtWeight);

    ////////////
    // Z mass //
    ////////////

    if (ZP4.M() < zMassCut[0] || ZP4.M() > zMassCut[1]) return kTRUE;  
    CountEvents(4);
    nEventsWeighted[4] += evtWeight;
    FillHistosNoise(4, evtWeight);

    /////////////////////
    // 3rd lepton veto //
    /////////////////////

    if (extraLeptons.size() > 0) return kTRUE;
    if ((electrons.size() > 0 || muons.size() > 2) && selection == "muon") return kTRUE;
    if ((muons.size() > 0 || electrons.size() > 2) && selection == "electron") return kTRUE;
    if ((muons.size() > 0 || electrons.size() > 0) && (selection == "eGamma" || selection == "muGamma" || selection == "gamma")) return kTRUE;
    CountEvents(5);
    nEventsWeighted[5] += evtWeight;
    FillHistosNoise(5, evtWeight);

    //////////
    // Z Qt //
    //////////

    if (ZP4.Pt() < qtCut) return kTRUE;  

    Bool_t passBasic = kTRUE;

    /////////////////////////////////////////////////
    // Variables for FillHistos() function (Andrey)//
    /////////////////////////////////////////////////
    MT      = CalculateTransMass(metP4, ZP4);
    METqt   = metP4.Pt()/ZP4.Pt();
    MET_phi = metP4.Phi();
    Mll     = ZP4.M(), Mll_EE=0, Mll_EB=0, Mll_EX=0;
    qT      = ZP4.Pt();
    diEta   = ZP4.Eta();

    if(fabs(Lepton1.Eta()) < 1.444 && fabs(Lepton2.Eta())<1.444) Mll_EB = Mll;
    else  if(fabs(Lepton1.Eta()) > 1.566 && fabs(Lepton2.Eta())>1.566)  Mll_EE = Mll;
    else  Mll_EX = Mll;
    nDofVtx1 = mainPrimaryVertex->NDof();
    if(PVind[1]!=-1) nDofVtx2 = secondPrimaryVertex->NDof();


    ////////////////
    // b-jet veto //
    ////////////////
    if (bJetP4.size() > 0) {
        h1_leadBJetPt->Fill(bJetP4[0].Pt(), evtWeight);
        h1_leadBJetEta->Fill(bJetP4[0].Eta(), evtWeight);
        h1_leadBJetPhi->Fill(bJetP4[0].Phi(), evtWeight);
    }

    Bool_t passBveto = kFALSE;
    if (nJetsB == 0)
      {
	passBveto = kTRUE;
	CountEvents(6);
	nEventsWeighted[6] += evtWeight;
	FillHistos(6, evtWeight);
	FillHistosNoise(6, evtWeight);
    
	//Yields for all Higgs masses//
	PostSelectionYieldCounter(nEventsWeighted[6], metP4, ZP4, deltaPhiJetMET, evtWeight); 
    }

  
    ////////////////////////
    // DeltaPhi(MET, jet) //
    ////////////////////////

    if (passBasic && deltaPhiJetMET > dPhiMinCut[1] &&  MET>metMinCut[1]  && (MT > mtMinCut[1] && MT < mtMaxCut[1]))
      {
	//NO b-veto, only basics cuts. For N-b jets distribution. For Higgs  300
	CountEvents(7);
	nEventsWeighted[7] += evtWeight;
       	FillHistos(7, evtWeight);
	FillHistosNoise(7, evtWeight);
	h1_bJetMultPreVeto->Fill(bJetP4.size(), evtWeight);
      }


    if (passBveto && deltaPhiJetMET > dPhiMinCut[0] &&  MET>metMinCut[0]  && (MT > mtMinCut[0] && MT < mtMaxCut[0]))
      {
	//For Higgs  250
	CountEvents(8);
	nEventsWeighted[8] += evtWeight;
	FillHistos(8, evtWeight);
	FillHistosNoise(8, evtWeight);
      }


    if (passBveto && deltaPhiJetMET > dPhiMinCut[1] &&  MET>metMinCut[1]  && (MT > mtMinCut[1] && MT < mtMaxCut[1]))
      {
	//For Higgs  300
	CountEvents(9);
	nEventsWeighted[9] += evtWeight;
       	FillHistos(9, evtWeight);
	FillHistosNoise(9, evtWeight);
      }

    if (passBveto && deltaPhiJetMET > dPhiMinCut[2] &&  MET>metMinCut[2]  && (MT > mtMinCut[2] && MT < mtMaxCut[2]))
      {
	//For Higgs  350
	CountEvents(10);
	nEventsWeighted[10] += evtWeight;
       	FillHistos(10, evtWeight);
	FillHistosNoise(10, evtWeight);

      }
    if (passBveto && deltaPhiJetMET > dPhiMinCut[3] &&  MET>metMinCut[3]  && (MT > mtMinCut[3] && MT < mtMaxCut[3]))
      {
	//For Higgs  400
	CountEvents(11);
	nEventsWeighted[11] += evtWeight;
       	FillHistos(11, evtWeight);
	FillHistosNoise(11, evtWeight);
      }
    if (passBveto && deltaPhiJetMET > dPhiMinCut[4] &&  MET>metMinCut[4]  && (MT > mtMinCut[4] && MT < mtMaxCut[4]))
      {
	//For Higgs  450
	CountEvents(12);
	nEventsWeighted[12] += evtWeight;
       	FillHistos(12, evtWeight);
    	FillHistosNoise(12, evtWeight);
      }


    /*
    ////////////
    // MET/QT //
    ////////////

    if (met->Met()/ZP4.Pt() < metByQtCut[0] || met->Met()/ZP4.Pt() > metByQtCut[1]) return kTRUE;
    CountEvents(11);
    nEventsWeighted[11] += evtWeight;

    //////////////
    // MET + QT //
    //////////////

    if ((metP4 + ZP4).Pt() < metPlusQtCut[0] || (metP4 + ZP4).Pt() > metPlusQtCut[1]) return kTRUE;
    CountEvents(12);
    nEventsWeighted[12] += evtWeight;

    //////////////////////
    // Jet Multiplicity //
    //////////////////////

    if ((int)jetP4.size() < nJetsCut[0] || (int)jetP4.size() > nJetsCut[1]) return kTRUE;
    CountEvents(13);
    nEventsWeighted[13] += evtWeight;
    */

    ////////////////////////////
    // Fill lepton histograms //
    ////////////////////////////

    if (selection == "muon") {

        LeptonBasicPlots(muons[0].P4(), muons[1].P4(), evtWeight);
        MetPlusLeptonPlots(metP4, muons[0].P4(), muons[1].P4(), evtWeight);

    } else if (selection == "electron") {

        LeptonBasicPlots(electrons[0].P4(), electrons[1].P4(), evtWeight);
        MetPlusLeptonPlots(metP4, electrons[0].P4(), electrons[1].P4(), evtWeight);
    } 

    DileptonBasicPlots(ZP4, evtWeight);

    ///////////////////////
    // MET distributions //
    ///////////////////////

    MetPlusZPlots(metP4, ZP4, evtWeight);
    h1_Met->Fill(metP4.Pt(), evtWeight);
    h1_MetPhi->Fill(met->Phi(), evtWeight);
    h1_MetSumEt->Fill(met->SumEt(), evtWeight);
    h1_ZPlusJetMHT->Fill((sumJetP4 + ZP4).Pt(), evtWeight);

    h1_ReducedMET->Fill(reducedMetP4.Pt());
    h1_ReducedMETTransverse->Fill(reducedMetP4.Px()); 
    h1_ReducedMETLongitudinal->Fill(reducedMetP4.Py());

    h1_MetNearestJetDeltaPhi->Fill(deltaPhiJetMET, evtWeight);

    /////////////////////////
    // Fill jet histograms //
    /////////////////////////

    h2_JetMultVsPVMult->Fill(primaryVtx->GetSize(), jetP4.size(), evtWeight);
    h1_jetMult->Fill(jetP4.size(), evtWeight);
    h1_bJetMultPostVeto->Fill(bJetP4.size(), evtWeight);

    if (jetP4.size() > 0) {
        h1_leadJetPt->Fill(jetP4[0].Pt(), evtWeight);
        h1_leadJetEta->Fill(jetP4[0].Eta(), evtWeight);
        h1_leadJetPhi->Fill(jetP4[0].Phi(), evtWeight);
    } 
    if (jetP4.size() > 1) {	
        h1_tailJetPt->Fill(jetP4[jetP4.size()-1].Pt(), evtWeight);
        h1_tailJetEta->Fill(jetP4[jetP4.size()-1].Eta(), evtWeight);
        h1_tailJetPhi->Fill(jetP4[jetP4.size()-1].Phi(), evtWeight);
    }

    //////////
    // misc //
    //////////

    h1_eventWeight->Fill(evtWeight);
    h1_pvMult->Fill(primaryVtx->GetSize(), evtWeight);
    if (!isRealData) h1_ptHat->Fill(ptHat, evtWeight);

    if (isRealData) {
        h1_goodRuns->Fill(runNumber);
        p1_nVtcs->Fill(runNumber, primaryVtx->GetSize());
    }

    return kTRUE;
}

void higgsAnalyzer::Terminate()
{

    cout<<"\nRunning over "<<suffix<<" dataset with "<<selection<<" selection."<<"\n"<<endl;
    cout<<"| CUT DESCRIPTION                    |\t"<< "\t|"<<endl;
    cout<<"| Initial number of events:          |\t"<< nEvents[0]  <<"\t|"<<nEventsWeighted[0]  <<"\t|"<<endl;
    cout<<"| Pass HLT selection:                |\t"<< nEvents[1]  <<"\t|"<<nEventsWeighted[1]  <<"\t|"<<endl;
    cout<<"| Cosmics veto:                      |\t"<< nEvents[2]  <<"\t|"<<nEventsWeighted[2]  <<"\t|"<<endl;
    cout<<"| Two oppositely charged leptons:    |\t"<< nEvents[3]  <<"\t|"<<nEventsWeighted[3]  <<"\t|"<<endl;
    cout<<"| Z mass window:                     |\t"<< nEvents[4]  <<"\t|"<<nEventsWeighted[4]  <<"\t|"<<endl;
    cout<<"| Third lepton veto:                 |\t"<< nEvents[5]  <<"\t|"<<nEventsWeighted[5]  <<"\t|"<<endl;
    cout<<"| b-jet veto:                        |\t"<< nEvents[6]  <<"\t|"<<nEventsWeighted[6]  <<"\t|"<<endl;
    cout<<"| Z Qt :                             |\t"<< nEvents[7]  <<"\t|"<<nEventsWeighted[7]  <<"\t|"<<endl;
    cout<<"| deltaPhi(MET, jet):                |\t"<< nEvents[8]  <<"\t|"<<nEventsWeighted[8]  <<"\t|"<<endl;
    cout<<"| MET:                               |\t"<< nEvents[9] <<"\t|"<<nEventsWeighted[9] <<"\t|"<<endl;
    cout<<"| MT:                                |\t"<< nEvents[10]  <<"\t|"<<nEventsWeighted[10]  <<"\t|"<<endl;
    cout<<"| MET/QT:                            |\t"<< nEvents[11] <<"\t|"<<nEventsWeighted[11] <<"\t|"<<endl;
    cout<<"| MET + QT:                          |\t"<< nEvents[12] <<"\t|"<<nEventsWeighted[12] <<"\t|"<<endl;
    cout<<"| Jet Multiplicity:                  |\t"<< nEvents[13] <<"\t|"<<nEventsWeighted[13] <<"\t|"<<endl;

    for (int i = 0; i < 16; ++i) {
        h1_acceptanceByCut->SetBinContent(i+1, nEventsWeighted[i]);
        h1_acceptanceByCut->SetBinError(i+1, sqrt(nEventsWeighted[i]));
        h1_acceptanceByCutRaw->SetBinContent(i+1, nEvents[i]);
        h1_acceptanceByCutRaw->SetBinError(i+1, sqrt(nEvents[i]));
    }


    //Normalization to 1fb and coloring the histograms
    string sel, sample;
    Float_t CS, nEv;  Int_t lColor, fColor;

    ifstream params;
    params.open("./params.txt");  //This is the file where parameters are stored
    params>>sel;
    params>>sample;  params>>nEv;  params>>CS;  params>>lColor;  params>>fColor;
    cout<<"sample: "<<sample<<"  cs: "<<CS<<"  events: "<<nEv<<" colors: "<<lColor<<"  "<<fColor<<endl;
    //ncout<<endl;

    //Getting total number of events fro ma histogram filled at the beginning
    Float_t nEv1 = getNevents(sample.c_str(), met0_et[0]);

    histoFile->Write();
    histoFile->cd("Andrey");
    scaleAndColor(sample.c_str(), CS, nEv1, 1000.0, lColor, fColor); //1fb

    histoFile->Close();  
}

void higgsAnalyzer::MetPlusZPlots(TLorentzVector metP4, TLorentzVector ZP4, float evtWeight)
{

    h1_MetOverQt->Fill(metP4.Pt()/ZP4.Pt(), evtWeight);
    h1_MetPlusQtMagnitude->Fill((metP4 + ZP4).Pt(), evtWeight);
    h1_MetQtDeltaPhi->Fill(fabs(metP4.DeltaPhi(ZP4)), evtWeight);
    h2_MetByMetOverQt->Fill(metP4.Pt(), metP4.Pt()/ZP4.Pt(), evtWeight);
    h2_MetByMetPlusQtMag->Fill(metP4.Pt(), (metP4 + ZP4).Pt(), evtWeight);

    //MT
    h1_MetDileptonMT->Fill(CalculateTransMass(metP4, ZP4), evtWeight);

    //ProjectedMET
    h1_ProjMetByQt->Fill(metP4.Pt()*cos(metP4.DeltaPhi(ZP4)), evtWeight);
    h1_OrthoMetByQt->Fill(metP4.Pt()*sin(metP4.DeltaPhi(ZP4)), evtWeight);
    if (fabs(metP4.DeltaPhi(ZP4)) < TMath::Pi()/2) {
        h1_ProjectedMet->Fill(metP4.Pt() * sin(fabs(metP4.DeltaPhi(ZP4))), evtWeight);
    } else {
        h1_ProjectedMet->Fill(metP4.Pt(), evtWeight);
    }
}

void higgsAnalyzer::MetPlusLeptonPlots(TLorentzVector metP4, TLorentzVector p1, TLorentzVector p2, float evtWeight)
{
    h1_MetLeadLeptonMT->Fill(CalculateTransMassAlt(metP4, p1), evtWeight);
    h1_MetTrailingLeptonMT->Fill(CalculateTransMassAlt(metP4, p1), evtWeight);
    h2_MetLeptonMT->Fill(CalculateTransMassAlt(metP4, p1), CalculateTransMassAlt(metP4, p2), evtWeight);
    h1_MetLeadLeptonDeltaPhi->Fill(fabs(metP4.DeltaPhi(p1)), evtWeight);
    h1_MetTrailingLeptonDeltaPhi->Fill(fabs(metP4.DeltaPhi(p2)), evtWeight);
}

void higgsAnalyzer::LeptonBasicPlots(TLorentzVector p1, TLorentzVector p2, float evtWeight)
{
    h1_leadLeptonPt->Fill(p1.Pt(), evtWeight);     
    h1_leadLeptonEta->Fill(p1.Eta(), evtWeight);    
    h1_leadLeptonPhi->Fill(p1.Phi(), evtWeight);    
    h1_trailingLeptonPt->Fill(p2.Pt(), evtWeight); 
    h1_trailingLeptonEta->Fill(p2.Eta(), evtWeight);
    h1_trailingLeptonPhi->Fill(p2.Phi(), evtWeight);
    h1_diLeptonPtRatio->Fill(p2.Pt()/p1.Pt(), evtWeight); 
    h1_diLeptonDeltaEta->Fill(fabs(p2.Eta() - p1.Eta()), evtWeight);
    h1_diLeptonDeltaPhi->Fill(fabs(p2.DeltaPhi(p1)), evtWeight);
    h1_diLeptonDeltaR->Fill(p2.DeltaR(p1), evtWeight);  
}

void higgsAnalyzer::DileptonBasicPlots(TLorentzVector ZP4, float evtWeight)
{
    h1_diLeptonTransMass->Fill(ZP4.Mt(), evtWeight);
    h1_diLeptonMass->Fill(ZP4.M(), evtWeight);     
    h1_diLeptonQt->Fill(ZP4.Pt(), evtWeight);       
}

float higgsAnalyzer::CalculateTransMass(TLorentzVector p1, TLorentzVector p2)
{
    //float transE    = sqrt(p1.Pt()*p1.Pt() + pow(91.2,2)) + sqrt(p2.Pt()*p2.Pt() + pow(91.2,2)); //p2.M()*p2.M());
    float transE    = sqrt(p1.Pt()*p1.Pt() + p2.M()*p2.M()) + sqrt(p2.Pt()*p2.Pt() +p2.M()*p2.M());
    float transPt   = (p1 + p2).Pt();
    float transMass = sqrt(transE*transE - transPt*transPt);

    return transMass;
}

float higgsAnalyzer::CalculateTransMassAlt(TLorentzVector p1, TLorentzVector p2)
{
    float transMass = sqrt(2*p2.Pt()*p1.Pt() * (1 - cos(fabs(p2.DeltaPhi(p1)))));
    return transMass;
}

bool higgsAnalyzer::CosmicMuonFilter(TCMuon muon1, TCMuon muon2)
{
    float dimuonAngle = muon1.P4().Angle(muon2.P4().Vect());
    if (TMath::Pi() - dimuonAngle  < 0.05)
        return true;
    else
        return false;
}

float higgsAnalyzer::DeltaPhiJetMET(TLorentzVector metP4, vector<TLorentzVector> jets)
{
    float minDeltaPhi = TMath::Pi();
    TLorentzVector nearestJet;
    for (int i = 0; i < (int)jets.size(); ++i) {
        if (fabs(metP4.DeltaPhi(jets[i])) < minDeltaPhi) {
            minDeltaPhi = fabs(metP4.DeltaPhi(jets[i]));
            nearestJet  = jets[i];
        }
    }
    //h1_nearestJetEta->Fill(nearestJet.Eta(), evtWeight);
    return fabs(nearestJet.DeltaPhi(metP4));
}

float higgsAnalyzer::GetEventWeight(int nPV, TLorentzVector l1, TLorentzVector l2)
{
    float weight = 1.; float puWeight = 1.; float triggerWeight = 1.; float zzWeight = 1.;

    if (!isRealData) {
        //PU reweighting parameterized by number of PVs
        puWeight = h_puReweight->GetBinContent(nPV+1); 

        //weighting for the muon trigger;include parameterization for muons
        if (selection == "muon") triggerWeight = 0.9;
        if (selection == "electron") triggerWeight = 0.99;

        //include dynamic scaling of ZZ samples
        if (suffix == "ZZ") {
            float ZP4 = (l1 + l2).Pt();
            zzWeight = 0.12 + (1.108 + 0.002429*ZP4 + (1.655e-6)*pow(ZP4, 2));  
        }
        weight *= puWeight;
        weight *= triggerWeight;
        weight *= zzWeight;

    } else if (selection == "eGamma") {
        int binNumber = floor(l1.Pt()/10. + 1);
        weight = h_eGammaReweight->GetBinContent(binNumber);
    } else if (selection == "muGamma") {
        int binNumber = floor(l1.Pt()/10. + 1);
        weight = h_muGammaReweight->GetBinContent(binNumber);
    }

    return weight;
}

float higgsAnalyzer::PUCorrectedMET(float met, int nPV, string metType)
{
    Float_t sigma0=0, shift=0, sigmaT=0, shift0=0, puCorMet=0;
    if (metType == "pfMet") {
        if (isRealData) {
            sigma0 = 5.868;
            shift0 = 8.231;
            shift  = 0.5761*(nPV - 1);
            sigmaT = sqrt(pow(sigma0, 2) + (nPV - 1)*pow(2.112, 2));
        } else {
            sigma0 = 4.569;
            shift0 = 6.282;
            shift  = 0.8064*(nPV - 1);
            sigmaT = sqrt(pow(sigma0, 2) + (nPV - 1)*pow(2.387, 2));
        }
    }
    puCorMet = (met - (shift0 + shift))*sigma0/sigmaT + shift0;
    return puCorMet;
}

TLorentzVector higgsAnalyzer::GetReducedMET(TLorentzVector sumJet, TLorentzVector lep1, TLorentzVector lep2, TLorentzVector metP4, int version) 
{
    TLorentzVector Q  = lep1 + lep2;
    float bisectorPhi = min(lep1.Phi(), lep2.Phi()) + lep1.DeltaPhi(lep2)/2;

    TVector2 projDilepton(Q.Px(), Q.Py());
    TVector2 projSumJet(sumJet.Px(), sumJet.Py());
    TVector2 projMET(metP4.Pt()*cos(metP4.Phi()), metP4.Pt()*sin(metP4.Phi()));

    //TVector2 delta;
    //TLorentzVector Thrust = lep1 - lep2; 
    //if (fabs(lep1.DeltaPhi(lep2)) > TMath::Pi()/2) {
    //	delta.Set(0., projDilepton.Rotate(-Thrust.Phi()).Py() + projDilepton.Py());
    //} else {
    //	delta.Set(0, projDilepton.Rotate(-Q.Phi()).Py() + projDilepton.Py());
    //}

    if (version == 1) {
        projDilepton = projDilepton.Rotate(-bisectorPhi);
        projSumJet   = projSumJet.Rotate(-bisectorPhi);
        projMET      = projMET.Rotate(-bisectorPhi);
    } else if (version == 2) {
        projDilepton = projDilepton.Rotate(-Q.Phi());
        projSumJet   = projSumJet.Rotate(-Q.Phi());
        projMET      = projMET.Rotate(-Q.Phi());
    }

    TVector2 unclustered = -1*projMET;                        //projDilepton - 1.*(projMET + projDilepton) + delta;
    TVector2 clustered   = projDilepton + 1.*projSumJet;   // + delta;

    TVector2 reducedMET = TVector2((fabs(unclustered.Px()) < fabs(clustered.Px()) ? unclustered.Px() : clustered.Px()),
            (fabs(unclustered.Py()) < fabs(clustered.Py()) ? unclustered.Py() : clustered.Py()));

    return TLorentzVector(reducedMET.Px(), reducedMET.Py(), 0, reducedMET.Mod());
}

float higgsAnalyzer::Dz(TVector3 objVtx, TLorentzVector objP4, TVector3 vtx)
{
    float vx  = objVtx.x(), vy = objVtx.y(), vz = objVtx.z();
    float px  = objP4.Px(), py = objP4.Py();
    float pz  = objP4.Pz(), pt = objP4.Pt();
    float pvx = vtx.x(), pvy = vtx.y(), pvz = vtx.z();

    float dZ =  (vz-pvz)-((vx-pvx)*px +(vy-pvy)*py)/pt*(pz/pt);
    return dZ;
}

float higgsAnalyzer::Dxy(TVector3 objVtx, TLorentzVector objP4, TVector3 vtx)
{
    float vx  = objVtx.x(), vy = objVtx.y();
    float px  = objP4.Px(), py = objP4.Py(), pt = objP4.Pt();
    float pvx = vtx.x(), pvy = vtx.y();

    double dXY =  (-(vx-pvx)*py + (vy-pvy)*px)/pt;
    return dXY;
}

void higgsAnalyzer::PostSelectionYieldCounter(float nEventsPS, TLorentzVector metP4, TLorentzVector ZP4, float dPhiJetMet, float evtWeight)
{
    for (int i = 0; i < 8; ++i) {

        h2_nEventsByHMass->Fill(i+0.5, 1., nEventsPS);

        if (dPhiJetMet < dPhiMinCut[i]) continue;
        h2_nEventsByHMass->Fill(i+0.5, 2., evtWeight);

        if (CalculateTransMass(metP4, ZP4) < mtMinCut[i] || CalculateTransMass(metP4, ZP4) > mtMaxCut[i]) continue;
        h2_nEventsByHMass->Fill(i+0.5, 3., evtWeight);

        if (metP4.Pt() < metMinCut[i]) continue;
        h2_nEventsByHMass->Fill(i+0.5, 4., evtWeight);
    }
}

float higgsAnalyzer::GetPhotonMass()
{
    float photonMass = 91.2;
    if (selection == "eGamma") {
        TH1D  *h_mass = (TH1D*)reweightFile.Get("h1_diElectronMass");
        photonMass = h_mass->GetRandom();
    }
    if (selection == "muGamma") {
        TH1D  *h_mass = (TH1D*)reweightFile.Get("h1_diMuonMass");
        photonMass = h_mass->GetRandom();
    }
    return photonMass;
}

void higgsAnalyzer::PrintOutNoisy(Int_t num){
    nout[num]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<MET<<"\t* "
	     <<endl;//       <<isNoiseHcal<<"\t* "<<isDeadEcalCluster<<"\t* "<<isScraping<<"\t* "<<isCSCTightHalo<<"\t* "<<endl;
}

void higgsAnalyzer::PrintOut(Int_t num){
  fout[num]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<Mll<<"\t* "<<MT<<"\t* "<<MET<<"\t* "
	   <<endl;//<<pTll<<"\t* "<<rhoFactor<<"\t* "<<muCountLoose<<"\t* "<<eleCountLoose<<endl;
}

void higgsAnalyzer::scaleAndColor(TString sample1, Float_t cs, Float_t nTotEv, Float_t lumi, Int_t line, Int_t fill)
{
  Float_t scaleFactor      = lumi*cs/nTotEv;
  Bool_t isData = kFALSE;
  if(sample1.Contains("DATA")) {scaleFactor = 1.0; isData = kTRUE;}
  if(sample1.Contains("HZZ")) scaleFactor*=10;

  //  cout<<sample1<<"   Rescaling all the histograms by: "<<scaleFactor<<endl;                                                                              
  //  file -> cd();                                                                                                                                          
  TIter nextkey( gDirectory->GetListOfKeys() );
  //TIter nextkey( file->GetListOfKeys() );                                                                                                                  
  TKey *key, *oldkey=0;

  while ( (key = (TKey*)nextkey()))
    {
      if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;
      TObject *obj = key->ReadObj();
      if ( ! (obj->IsA()->InheritsFrom( TH1::Class() )) ) continue;
      TH1 *hist = (TH1*)key->ReadObj();
      hist -> Scale(scaleFactor);
      if(line!=-1)      hist -> SetLineColor(line);
      if(fill!=-1)      hist -> SetFillColor(fill);
      if(isData)        hist -> SetMarkerStyle(20);
      hist->SetLineWidth(2);
      hist->Write("",kWriteDelete); //overwrite the object in the file                                                                                       
      delete obj;
      //cout<<sample1<<"  nEv: "<<nTotEv<<"   cs: "<<cs<<"   Rescaling all the histogram "<<hist->GetName()<<" by: "<<scaleFactor<<endl;                     
    }

}

float higgsAnalyzer::getNevents(string dataset, TH1F* h)
{

  Int_t nEv2 = 0;
  //Get total number of Events  from the histogram (needs to be filled before any cuts applied)
  nEv2 = h       -> GetEntries();
  //nEv2 = 10000;

  if(dataset=="DATA")          nEv2 = -1.;
  else if(dataset=="DYmumu")   nEv2 = 29.5e6;
  else if(dataset=="DYee")     nEv2 = 29.5e6;
  else if(dataset=="DYtautau") nEv2 = 29.5e6;
  else if(dataset=="ttbar")    nEv2 = 10.34e6;
  else if(dataset=="tW")       nEv2 = 795e3;
  else if(dataset=="WW")       nEv2 = 220e3;
  else if(dataset=="ZZ")       nEv2 = 220e3;
  else if(dataset=="WZ")       nEv2 = 205e3;
  else if(dataset=="Wjets")    nEv2 = 31314.0;
  else if(dataset=="ggHZZ140") nEv2 = 99.99e3;
  else if(dataset=="ggHZZ200") nEv2 = 99.62e3;
  else if(dataset=="ggHZZ250") nEv2 = 99.99e3;
  else if(dataset=="ggHZZ300") nEv2 = 96.99e3;
  else if(dataset=="ggHZZ350") nEv2 = 96.99e3;
  else if(dataset=="ggHZZ400") nEv2 = 93.90e3;
  else if(dataset=="ggHZZ600") nEv2 = 93.90e3;
  else if(dataset=="VBFHZZ200") nEv2  = 93.90e3;

  cout<<"Getting events for sample: "<<dataset<<"  :  "<<nEv2<<endl;
  return nEv2;
}

void higgsAnalyzer::FillHistosNoise(Int_t num, Double_t weight){

  // met0_over_qt[num] -> Fill(METqt, weight);
  met0_et[num]      -> Fill(MET, weight);
  //met0_et_ovQt[num] -> Fill(MET, METqt, weight);
  //met0_phi[num]     -> Fill(MET_phi, weight);

  //met1_over_qt[num] -> Fill(MET1qt, weight);
  //met1_et[num]      -> Fill(MET1, weight);
  //met1_et_ovQt[num] -> Fill(MET1, MET1qt, weight);
  //met1_phi[num]     -> Fill(MET1_phi, weight);

  //mt0[num]          -> Fill(MT, weight);

}

void higgsAnalyzer::FillHistos(Int_t num, Double_t weight){

  met2_over_qt[num] -> Fill(METqt, weight);
  met2_et[num]      -> Fill(MET, weight);
  met2_et_ovQt[num] -> Fill(MET, METqt, weight);
  met2_phi[num]     -> Fill(MET_phi, weight);

  met3_et[num]      -> Fill(puCorrMET, weight);

  //met4_over_qt[num] -> Fill(puCorrMETqt, weight);
  //met4_et[num]      -> Fill(puCorrMET, weight);
  //met4_et_ovQt[num] -> Fill(puCorrMET, puCorrMETqt, weight);
  //met4_puSig[num]   -> Fill(puSigMET, weight);

  mt2[num]      -> Fill(MT, weight);
  mt2_met2[num] -> Fill(MT, MET, weight);

  //mt2_met3[num] -> Fill(MT, projMET, weight);
  //mtZ_met3[num] -> Fill(MTZ, projMET, weight);


  di_qt[num]   -> Fill(qT, weight);
  di_eta[num]  -> Fill(diEta);
  di_mass[num] -> Fill(Mll, weight);
  if(Mll_EB!=0) di_mass_EB[num]  -> Fill(Mll_EB, weight);
  if(Mll_EE!=0) di_mass_EE[num]  -> Fill(Mll_EE, weight);
  if(Mll_EX!=0) di_mass_EX[num]  -> Fill(Mll_EX, weight);

  jet_N[num]    -> Fill(nJets, weight);
  jet_b_N[num]  -> Fill(nJetsB, weight);

  vtx_nPV_raw[num]    -> Fill(nVtx);
  vtx_nPV_weight[num] -> Fill(nVtx, weight);
  vtx_ndof_1[num]     -> Fill(nDofVtx1, weight);
  vtx_ndof_2[num]     -> Fill(nDofVtx2, weight);


  /*
  if(jetCount>0)
    {
      met2_dPhiClosJet1[num] -> Fill(dPhiClos1, weight);
      met2_dPhiClosJet2[num] -> Fill(dPhiClos2, weight);
      met2_dPhiLeadJet1[num] -> Fill(dPhiLead1, weight);
      met2_dPhiLeadJet2[num] -> Fill(dPhiLead2, weight);
    }
  */

}

void higgsAnalyzer::CountEvents(Int_t num)
{
  nEvents[num]++;
  /*
  if(!isNoiseHcal)        nEventsPassNoiseFilter[1][num]++;
  if(!isDeadEcalCluster)  nEventsPassNoiseFilter[2][num]++;
  if(!isScraping)         nEventsPassNoiseFilter[3][num]++;
  if(!isCSCTightHalo)     nEventsPassNoiseFilter[4][num]++;
  if(!isCSCLooseHalo)     nEventsPassNoiseFilter[5][num]++;
  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo)
    nEventsPassNoiseFilter[0][num]++;
  */
}
