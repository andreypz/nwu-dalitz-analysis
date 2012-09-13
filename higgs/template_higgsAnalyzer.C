// $Id: template_higgsAnalyzer.C,v 1.51 2012/09/13 16:14:46 andrey Exp $

#define higgsAnalyzer_cxx

#include "higgsAnalyzer.h"
#include <string>
#include "../plugins/MetDefinitions.h"
using namespace std;

// This needs to be un-commented in order to run MVA cuts
//#define USE_MVA  

/////////////////////////////
//Specify parameters here. //
/////////////////////////////

const string  selection      = "muon";
const string  period         = "2011";
const int     JC_LVL         = 0;  //No JEC for Pat jets (they are already applied)
const int     trigger[]      = {0};
const TString suffix("TEST");

vector<int> triggers (trigger, trigger + sizeof(trigger)/sizeof(int));

Bool_t makeMvaTree = 0;
const UInt_t verboseLvl  = 0;
const Bool_t doZlibrary  = 0, isFromData=0;

/////////////////
//Analysis cuts//
/////////////////

const Float_t jetPtCut[]     = {30., 15.};
const Float_t bJetPtCut      = 30.;
const Float_t muPtCut[]      = {20., 10.};
const Float_t elePtCut[]     = {20., 10.};
const Float_t photonPtCut[]  = {25., 1e9};
const Float_t zMassCut1[]     = {76, 106};
const Float_t zMassCut2[]     = {70, 106};
const Float_t qtCut          = 55.;
const Int_t   nJetsCut[]     = {0, 99};
const Float_t cut_vz = 24, cut_vd0 = 2, cut_vndof = 4;  //PV filter cuts

// Cuts for mass points in the PAS. 200, 250,300,350 etc
const Float_t dPhiMinCut = 0.5;
const Float_t metMinCut[] = {75,   75,   85,  85,  90, 110,  110, 120, 120}; 
const Float_t mtMinCut[]  = {175, 225,  250,  300,  325,  350,  400, 450, 450};
const Float_t mtMaxCut[]  = {275, 325,  350,  400,  425,  500,  650, 700, 9999};

Float_t mva_weight, mva_nJets15;

// Do something about these: should just have one sort condition function
bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());} 

TRandom3 *myRandom;


// ************* MVA ********* //

// stuff for BDT selection ----------------------------------------------------------------------------
// inclusions made with the 'define'
// to isolate it from the rest of the code

#ifdef USE_MVA

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

// TMVA weights directory
const TString weightsDir = "../data/mvaWeights";

#define N_DISCR_METHODS 3
// here we will use only BDTG... but will keep the structure 
enum DISCR_METHOD {
  MLPBNN, BDTG, BDT
};
const TString discrMethodName[3] = {
  "MLPBNN", "BDTG", "BDT"
};

const TString discrSampleName = "allBg";
// multiplicities of jets for the discriminants                                                                                                                 
#define  N_DISCR_JET_MULTI 3

enum DISCR_JET_MULTI {
  D_0J, D_1J, D_GT1J
};
const TString discrJetMultiName [N_DISCR_JET_MULTI] = {"0j", "1j", "gt1j"};


#define N_HIGGS_MASSES 3
const Int_t mvaHiggsMassPoint[N_HIGGS_MASSES] = {125,200,250};

const Float_t bdtCut[N_HIGGS_MASSES][3] = 
  {
    {-0.15, -0.07, -0.05},  //m125
    {0.0, 0.06, 0.05},      //m200
    {0.072, 0.092, 0.082},     //m250
  };

//[jet multi bin]
TMVA::Reader* tmvaReader[N_HIGGS_MASSES][3];

//Variables for the MVA tree
#endif
///****** End of MVA related stuff *************//



//For Synchronization
//The events that Jucub cuts, but I keep:
//UInt_t myEVTS[] = {};
//Visa versa
const UInt_t hisEVTS[] = {3912507,3537222, 2260583, 4048415, 1972146, 1294493, 2879208, 65882, 696163, 4124071, 3263050, 640272};


void higgsAnalyzer::Begin(TTree * /*tree*/) 
{
    TString option = GetOption();
    TH1::SetDefaultSumw2(kTRUE);

    cout<<"\n***      Begin the Analyzer      ****"<<endl;
    myRandom = new TRandom3();

    //Rochester Corrections for muons:

    if (selection=="muon")
      roch = new rochcor();    

    // Initialize utilities and selectors here //
    zLib            = new ZedEventsLibrary(selection, isFromData);
    weighter        = new WeightUtils(suffix.Data(), period, selection, isRealData);
    triggerSelector = new TriggerSelector(selection, period, triggers);

    histoFile = new TFile("a_higgsHistograms.root", "RECREATE");
    histoFile->mkdir("Andrey", "Andrey");
    histoFile->mkdir("mvaTree", "mvaTree");
    histoFile->cd("Andrey");

    //In the Title of this histogram I encode the sample!  

    for(Int_t m=0;m<3;m++) //For Higgs masses
      for(Int_t j=0;j<3;j++) //For three bins in jet multiplicity
	mva_discr[m][j] = new TH1F(Form("mva_discr_mh%i_%i",m,j), "MVA discriminator output", 50, -0.6,0.3);

    evt_byCut = new TH1F("evt_byCut", "TEST", nC, 0,nC);
    evt_byCut_raw = new TH1F("evt_byCut_raw", "TEST", nC, 0,nC);
    evt_libQt = new TH2F("evt_libQt", "Events in a library", 4, 0,4,  10, 0,10);

    if (doZlibrary)  zLib->CreateLibrary();
    
    for (Int_t n=0; n<nC; n++)
      {
      //cout<<"zeroing the array: "<<n<<endl;
      nEvents[n]=0;
      nEventsWeighted[n]=0;
    }

    for(Int_t n=0; n<nC; n++)
      {
	met0_et[n]       = new TH1F(Form("met0_et_%i",n), "met0_et", 40, 0,400);
	evt_weight[n]      = new  TH1F(Form("evt_weight_%i",n), "Weights", 50,0,10);
      }

    for(Int_t n=4; n<nC; n++)
      {
	//cout<<n<<"  dbg hists"<<endl;

	met1_et[n]      = new TH1F(Form("met1_et_%i",n), "met1_et", 80, 0,400);
	met2_et[n]      = new TH1F(Form("met2_et_%i",n), "met2_et", 80, 0,400);
	met3_et[n]      = new TH1F(Form("met3_et_%i",n), "met3_et", 84, -20,400);
	met4_et[n]      = new TH1F(Form("met4_et_%i",n), "met4_et", 80, 0,400);
	met5_et[n]      = new TH1F(Form("met5_et_%i",n), "met5_et", 80, 0,400);
	met6_et[n]      = new TH1F(Form("met6_et_%i",n), "met6_et", 80, 0,400);
	met7_et[n]      = new TH1F(Form("met7_et_%i",n), "met7_et", 80, 0,400);
	met8_et[n]      = new TH1F(Form("met8_et_%i",n), "met8_et", 80, 0,400);
	met9_et[n]      = new TH1F(Form("met9_et_%i",n), "met9_et", 80, 0,400);
	met10_et[n]     = new TH1F(Form("met10_et_%i",n), "met10_et", 80, 0,400);

	met1_overQt[n] = new TH1F(Form("met1_overQt_%i",n), "met1_overQt", 40, 0,4);
	met1_projOnQt[n]    = new TH1F(Form("met1_projOnQt_%i",n), "Met projected on QT", 50, -200,100);
	met1_perpQt[n]      = new TH1F(Form("met1_perpQt_%i",n), "Met perpendicular to QT", 50, 0,400);
	met1_phi[n]     = new TH1F(Form("met1_phi_%i",n), "met1_phi", 40, -TMath::Pi(), TMath::Pi());


	//met0_et_ovQt[n]  = new TH2F(Form("met0_et_ovQt_%i",n), "pfMET vs Met/qt", 40, 0,400, 40, 0,4);
	//met1_et_ovQt[n]  = new TH2F(Form("met1_et_ovQt_%i",n), "MET1 vs Met/qt", 40, 0,400, 40, 0,4);
	//met2_et_ovQt[n]  = new TH2F(Form("met2_et_ovQt_%i",n), "pfMet noise vs Met/qt", 40, 0,400, 40, 0,4);
	//met3_et_ovQt[n]  = new TH2F(Form("met3_et_ovQt_%i",n), "projMET vs Met/qt", 40, 0,400, 40, 0,4);
	//met4_et_ovQt[n]  = new TH2F(Form("met4_et_ovQt_%i",n), "puCorrMET vs Met/qt", 40, 0,400, 40, 0,4);

	//met1_lg[n]         = new TH1F(Form("met1_lg_%i",n), "met1_lg", 50,-200,100);
	met1_recoil_lg[n]  = new TH1F(Form("met1_recoil_lg_%i",n), "met1_recoil_lg", 50,-400,100);

	//met2_phi[n]     = new TH1F(Form("met2_phi_%i",n), "met2_phi", 40, -TMath::Pi(), TMath::Pi());
	//met2_over_qt[n] = new TH1F(Form("met2_over_qt_%i",n), "met2_over_qt", 40, 0,4);

	//met3_phi[n]     = new TH1F(Form("met3_phi_%i",n), "met3_phi", 40, -TMath::Pi(), TMath::Pi());
	//met3_over_qt[n] = new TH1F(Form("met3_over_qt_%i",n), "met2_over_qt", 40, 0,4);

	//met4_phi[n]     = new TH1F(Form("met4_phi_%i",n), "met4_phi", 40, -TMath::Pi(), TMath::Pi());
	//met4_over_qt[n] = new TH1F(Form("met4_over_qt_%i",n), "met4_over_qt", 42, -0.2,4);
	//	met4_puSig[n]   = new TH1F(Form("met4_sig_%i",n), "met4_sig", 42, -0.4,8);

	//met1_dPhiLeadJet1[n] =  new TH1F(Form("met1_dPhiLeadJet1_%i",n), "delta phi to lead jets", 50, 0, TMath::Pi());
	//met1_dPhiLeadJet2[n] =  new TH1F(Form("met1_dPhiLeadJet2_%i",n), "delta phi to lead jets", 50, 0, TMath::Pi());
	met1_dPhiClosJet1[n] =  new TH1F(Form("met1_dPhiClosJet1_%i",n), "delta phi to closest jets", 50, 0, TMath::Pi());
	//met1_dPhiClosJet2[n] =  new TH1F(Form("met1_dPhiClosJet2_%i",n), "delta phi to closest jets", 50, 0, TMath::Pi());


	//met2_dPhiLeadJet1[n] =  new TH1F(Form("met2_dPhiLeadJet1_%i",n), "delta phi to lead jets", 50, 0, TMath::Pi());
	//met2_dPhiLeadJet2[n] =  new TH1F(Form("met2_dPhiLeadJet2_%i",n), "delta phi to lead jets", 50, 0, TMath::Pi());
	//met2_dPhiClosJet1[n] =  new TH1F(Form("met2_dPhiClosJet1_%i",n), "delta phi to closest jets", 50, 0, TMath::Pi());
	//met2_dPhiClosJet2[n] =  new TH1F(Form("met2_dPhiClosJet2_%i",n), "delta phi to closest jets", 50, 0, TMath::Pi());


	l1_phi[n] = new TH1F(Form("l1_phi_%i",n), "l1_phi", 50, -TMath::Pi(), TMath::Pi());
	l1_eta[n] = new TH1F(Form("l1_eta_%i",n), "l1_eta", 52, -2.6, 2.6);
	l1_pt[n]  = new TH1F(Form("l1_pt_%i",n), "l1_pt", 50, 0, 200);
	l2_phi[n] = new TH1F(Form("l2_phi_%i",n), "l2_phi", 50, -TMath::Pi(), TMath::Pi());
	l2_eta[n] = new TH1F(Form("l2_eta_%i",n), "l2_eta", 52, -2.6, 2.6);
	l2_pt[n]  = new TH1F(Form("l2_pt_%i",n), "l2_pt", 50, 0, 200);

	l0_angle[n]    = new TH1F(Form("l0_angle_%i",n), "l0_angle", 50, 0, TMath::Pi());
	l0_angleLog[n] = new TH1F(Form("l0_angleLog_%i",n), "l0_angleLog", 50, -5, 1);
	l0_dPhi[n]    = new TH1F(Form("l0_dPhi_%i",n), "l0_dPhi", 50, 0, TMath::Pi());
	l0_dEta[n]    = new TH1F(Form("l0_dEta_%i",n), "l0_dEta", 50, 0, 4);
	l0_dR[n]      = new TH1F(Form("l0_dR_%i",n), "l0_dR", 50, 0, 4);
	l0_ptRatio[n] = new TH1F(Form("l0_ptRatio_%i",n), "l0_ptRatio", 50, 0, 1);

	btag_hp[n] = new TH1F(Form("btag_hp_%i",n), "btag_hp", 40, -10, 10);

	//mtZ[n]     = new TH1F(Form("mtZ_%i",n), "MTZ pf", 50, 100, 600);
	mt0[n]     = new TH1F(Form("mt0_%i",n), "MT pf", 50, 100, 600);
	mt1[n]     = new TH1F(Form("mt1_%i",n), "MT type1", 50, 100, 600);
	mt2[n]     = new TH1F(Form("mt2_%i",n), "MT noise", 50, 100, 600);
	mt3[n]     = new TH1F(Form("mt3_%i",n), "MT proj", 50, 100, 600);
	mt4[n]     = new TH1F(Form("mt4_%i",n), "MT pu corr", 50, 100, 600);

	//mtZ_met2[n]   = new TH2F(Form("mtZ_met2_%i",n), "MTZ pf vs pfMet noise", 50, 100, 600, 40, 0,400);
	//mt2_met2[n]   = new TH2F(Form("mt2_met2_%i",n), "MT pf vs pfMet noise", 50, 100, 600, 40, 0,400);

	//mtZ_met3[n]   = new TH2F(Form("mtZ_met3_%i",n), "MTZ pf vs projMet noise", 50, 100, 600, 40, 0,400);
	//mt2_met3[n]   = new TH2F(Form("mt2_met3_%i",n), "MT pf vs projfMet noise", 50, 100, 600, 40, 0,400);

	di_qt[n]      = new TH1F(Form("di_qt_%i",n), "di-lepton qt", 80, 0,400);
	di_eta[n]     = new TH1F(Form("di_eta_%i",n), "di-lepton Eta", 100, -5,5);
	di_phi[n]     = new TH1F(Form("di_phi_%i",n), "di-lepton Phi", 50,  -TMath::Pi(), TMath::Pi());
	di_mass[n]    = new TH1F(Form("di_mass_%i",n), "di-lepton Mass", 50, 70,120);
	di_mass_EB[n] = new TH1F(Form("di_mass_EB_%i",n), "di-lepton Mass in Ecal Barrel", 50, 70,120);
	di_mass_EE[n] = new TH1F(Form("di_mass_EE_%i",n), "di-lepton Mass in Ecal Endcap", 50, 70,120);
	di_mass_EX[n] = new TH1F(Form("di_mass_EX_%i",n), "di-lepton Mass in Ecal E/B mix", 50, 70,120);
	di_dPhiMet[n] = new TH1F(Form("di_dPhiMet_%i",n), "di_dPhiMet", 50, 0, TMath::Pi());

	jet_N[n]      = new TH1F(Form("jet_N_%i",n), "Number of jets", 20,0,20);
	jet_N15[n]    = new TH1F(Form("jet_N15_%i",n), "Number of jets pt>15", 20,0,20);
	jet_N24[n]    = new TH1F(Form("jet_N24_%i",n), "Number of jets in |eta|<2.4", 20,0,20);
	jet_dRlep1[n] = new TH1F(Form("jet_dRlep1_%i",n), "dR(jet, lepton1)", 50, 0,5);
	jet_dRlep2[n] = new TH1F(Form("jet_dRlep2_%i",n), "dR(jet, lepton2)", 50, 0,5);
	jet_pt1[n]     = new TH1F(Form("jet_pt1_%i",n),  "Pt of leading jet", 50,0,400);
	jet_eta1[n]    = new TH1F(Form("jet_eta1_%i",n), "Eta of leading jet", 50,-5,5);
	jet_phi1[n]    = new TH1F(Form("jet_phi1_%i",n), "Phi of leading jet", 50, -TMath::Pi(), TMath::Pi());
	jet_pt2[n]     = new TH1F(Form("jet_pt2_%i",n),  "Pt of second jet", 50,0,400);
	jet_eta2[n]    = new TH1F(Form("jet_eta2_%i",n), "Eta of second jet", 50,-5,5);
	jet_phi2[n]    = new TH1F(Form("jet_phi2_%i",n), "Phi of second jet", 50, -TMath::Pi(), TMath::Pi());

	jet_diM[n]      = new TH1F(Form("jet_diM_%i",n), "Mass of di-jet system", 50, 0,2000);
	jet_deltaEta[n] = new TH1F(Form("jet_deltaEta_%i",n), "Delta eta of two leading jets", 50, 0,10);
	jet_zeppZ[n]    = new TH1F(Form("jet_zeppZ_%i",n), "Zeppenfeld variable wrt Z", 50,-5,5);
        jet_zeppZy[n]   = new TH1F(Form("jet_zeppZy_%i",n), "Zeppenfeld-y variable wrt Z", 50,-5,5);

	jet_b_N[n]    = new TH1F(Form("jet_b_N_%i",n), "Number of b-jets pt>default", 10,0,10);
	jet_b_Nssv[n] = new TH1F(Form("jet_b_Nssv_%i",n), "Number of b-jets pt>default", 10,0,10);
	jet_b_N25[n]  = new TH1F(Form("jet_b_N25_%i",n), "Number of b-jets pt>25", 10,0,10);
	jet_b_N30[n]  = new TH1F(Form("jet_b_N30_%i",n), "Number of b-jets pt>30", 10,0,10);
	jet_b_pt[n]   = new TH1F(Form("jet_b_pt_%i",n), "Pt of b-jets, eta<2.4", 50,0,400);


	vtx_nPV_tot[n]    = new  TH1F(Form("vtx_nPV_tot_%i",n), "Total Number of PVs, raw", 20,0,20);
	vtx_nPV_raw[n]    = new  TH1F(Form("vtx_nPV_raw_%i",n), "Number of PVs, raw", 20,0,20);
	vtx_nPV_weight[n] = new  TH1F(Form("vtx_nPV_weight_%i",n), "Number of PVs, reweighted", 20,0,20);
	vtx_ndof_1[n]     = new  TH1F(Form("vtx_ndof_1_%i",n), "Ndof of first PV", 50,0,150);
	vtx_ndof_2[n]     = new  TH1F(Form("vtx_ndof_2_%i",n), "Ndof of second PV", 50,0,150);

	ph_nGamma[n]      = new  TH1F(Form("ph_nGamma_%i",n), "photon multiplicity (pt>25)", 10,0,10);

	run_events[n]     =  new  TH1F(Form("run_events_%i",n), "Events per run", 18000, 160000., 178000.);


	if (suffix.Contains("ggHZZ") || suffix.Contains("ggHWW")) {
	  higgs_pt[n]     = new  TH1F(Form("higgs_pt_%i",n), "Higgs pt", 50,0,500);
	  higgs_w_pt[n]   = new  TH1F(Form("higgs_w_pt_%i",n), "Higgs pt, wighted", 50,0,500);
	  
	  higgs_mass[n]   = new  TH1F(Form("higgs_mass_%i",n), "Higgs mass", 150, 0,1500);
	  higgs_w_mass[n] = new  TH1F(Form("higgs_w_mass_%i",n), "Higgs mass, weighted", 150, 0,1500);
	  //higgs_mass[n] = new  TH1F(Form("higgs_mass_%i",n), "Higgs mass", 70,100,800);
	}
       

       }
    if (verboseLvl>0){
      ffout.open("./events_printout_TEST_final.txt",ofstream::out);
      ncout.open("./counts_for_tex_TEST.txt",ofstream::out);
      
      ffout.precision(4); ffout.setf(ios::fixed, ios::floatfield);
      for(Int_t i=0; i<nC; i++)
	{
	  if (i!=2 && i!=29) continue;
	  fout[i].open(Form("./events_printout_TEST_%i.txt",i),ofstream::out);
	  fout[i].precision(3); fout[i].setf(ios::fixed, ios::floatfield);
	  fout[i]<<
	    "*********\t*********\t*****\t***\t****\t****\t************\t*********\t********\t*****\t*\n"<<
	    "* Run    \t* Event  \t* LS \t* nVtx\t* nJets\t* m(ll)\t*   MT(llvv)\t* PFMET  \t* PT(ll)\t* rho\t*\n"<<
	    "*********\t*********\t*****\t****\t****\t***\t************\t*********\t********\t*****\t*\n";
	  if (verboseLvl>2){
	    nout[i].open(Form("./events_filtered_TEST_%i.txt",i),ofstream::out);
	    nout[i].precision(3); fout[i].setf(ios::fixed, ios::floatfield);
	    nout[i]<<
	      "*********\t*********\t*********\t******************************************\n"<<
	      "* Run    \t* Event  \t*  LS   \t*  PFMET * Hcal * Ecal * scrap * CSCTight * \n"<<
	      "*********\t*********\t*********\t******************************************\n";
	  }
      }
    }


    if (suffix.Contains("DATA"))
      makeMvaTree = kFALSE; //Dont produce mvaTrees for Data!

    histoFile->cd("mvaTree");    
    if (makeMvaTree){
      TString mvaTreeName = "mvaTree";
      _mvaTree = new TTree(mvaTreeName, "Tree with kinematic info");
      
      _mvaTree->Branch("lep1Pt", &lep1_pt, "lep1Pt/F");
      _mvaTree->Branch("lep2Pt", &lep2_pt,  "lep2Pt/F");
      _mvaTree->Branch("diLepPt", &qT, "diLepPt/F");
      _mvaTree->Branch("diLepM", &Mll, "diLepM/F");
      _mvaTree->Branch("diLepEta", &diEta, "diLepEta/F");
      _mvaTree->Branch("lepDeltaPhi", &lep_dPhi, "lepDeltaPhi/F");
      _mvaTree->Branch("lepDeltaEta", &lep_dEta, "lepDeltaEta/F");
      _mvaTree->Branch("lepDeltaR", &lep_dR, "lepDeltaR/F");
      _mvaTree->Branch("lepAngle", &lep_angle, "lepAngle/F");
      _mvaTree->Branch("lepPtRatio", &lep_ptRatio, "lepPtRatio/F");
      _mvaTree->Branch("met", &MET, "met/F");
      _mvaTree->Branch("metOverQt", &metOverQt, "metOverQt/F");
      _mvaTree->Branch("metProjOnQt", &metProjOnQt, "metProjOnQt/F");
      _mvaTree->Branch("metPerpQt", &metPerpQt, "metPerpQt/F");
      _mvaTree->Branch("dPhiMetDiLep", &dPhiMetDiLep, "dPhiMetDiLep/F");
      _mvaTree->Branch("dPhiJetMet", &dPhiClos1, "dPhiJetMet/F");
      _mvaTree->Branch("mt", &MT, "mt/F");
      _mvaTree->Branch("nJets", &nJets,  "nJets/I");
      _mvaTree->Branch("nJets15", &nJets15,  "nJets15/I");
      _mvaTree->Branch("deltaEtaDiJet,", &deltaEtaDiJet, "deltaEtaDiJet/F");
      _mvaTree->Branch("massDiJet,", &massDiJet, "massDiJet/F");
      _mvaTree->Branch("zeppDiJetDiLep,", &zeppDiJetDiLep, "zeppDiJetDiLep/F");
      _mvaTree->Branch("weight", &mva_weight, "weight/F");
      _mvaTree->Branch("evNumber", &eventNumber, "evNumber/l");
    }


    //cout<<"dbg   End of Begin job"<<endl;


    // -----------------  MVA stuff ---------------------

#ifdef USE_MVA

    // set up the TMVA readers (the input variables have been set already as global)
    for (Int_t jm = 0; jm <3; ++jm)
      {
	for (int mh = 0; mh < N_HIGGS_MASSES; ++mh) {
	  tmvaReader[mh][jm] = new TMVA::Reader("!Color:!Silent");
	  
	  //add  variables... some are exclusive to particular sample/jet multi discriminators
	  tmvaReader[mh][jm]->AddVariable("lep1Pt", &lep1_pt);
	  tmvaReader[mh][jm]->AddVariable("lep2Pt", &lep2_pt);
	  tmvaReader[mh][jm]->AddVariable("diLepPt", &qT);
	  tmvaReader[mh][jm]->AddVariable("diLepM", &Mll);
	  tmvaReader[mh][jm]->AddVariable("diLepEta", &diEta);
	  tmvaReader[mh][jm]->AddVariable("lepDeltaPhi", &lep_dPhi);
	  tmvaReader[mh][jm]->AddVariable("lepDeltaEta", &lep_dEta);
	  tmvaReader[mh][jm]->AddVariable("lepDeltaR", &lep_dR);
	  tmvaReader[mh][jm]->AddVariable("lepAngle", &lep_angle);
	  tmvaReader[mh][jm]->AddVariable("lepPtRatio", &lep_ptRatio);
	  
	  //tmvaReader[mh][jm]->AddVariable("deltaPhiLep", &deltaPhiLep);
	  tmvaReader[mh][jm]->AddVariable("met", &MET);
	  tmvaReader[mh][jm]->AddVariable("metOverQt", &metOverQt);
	  tmvaReader[mh][jm]->AddVariable("metProjOnQt", &metProjOnQt);
	  tmvaReader[mh][jm]->AddVariable("metPerpQt", &metPerpQt);
	  tmvaReader[mh][jm]->AddVariable("dPhiMetDiLep", &dPhiMetDiLep);
	  tmvaReader[mh][jm]->AddVariable("mt", &MT);
	  tmvaReader[mh][jm]->AddVariable("nJets15", &mva_nJets15);
	  //tmvaReader[mh][jm]->AddVariable("dPhiJetMet", &dPhiJetMet,);
	  if(jm==2){
	    tmvaReader[mh][jm]->AddVariable("deltaEtaDiJet", &deltaEtaDiJet);
	    tmvaReader[mh][jm]->AddVariable("massDiJet", &massDiJet);
	    tmvaReader[mh][jm]->AddVariable("zeppDiJetDiLep", &zeppDiJetDiLep);
	  }
	}
      }
              
    // Book the methods
    for (Int_t jm = 0; jm <3; ++jm)
      {
	int discr = BDT;

	for (int mh = 0; mh < N_HIGGS_MASSES; ++mh) {

	  TString label = TString::Format("%s_%s_MVA_hzz%i_%s", discrMethodName[discr].Data(), discrSampleName.Data(),
					  mvaHiggsMassPoint[mh],discrJetMultiName[jm].Data());

	  //discr_allBgPhoton_hzz250_0j___BDT.weights.xml
	  TString weightFile = TString::Format("%s/discr_%s_hzz%i_%s___%s.weights.xml",
					   weightsDir.Data(), discrSampleName.Data(), 
					   mvaHiggsMassPoint[mh],discrJetMultiName[jm].Data(),discrMethodName[discr].Data());

	  tmvaReader[mh][jm]->BookMVA(label.Data(), weightFile.Data());
	  
	} //hm loop

      }// loop over jet multiplicity bin


#endif
    // ------------------ End of MVA stuff ------------

}

bool higgsAnalyzer::Process(Long64_t entry)
{  
    GetEntry(entry);
    CountEvents(0);
    ++nEventsWeighted[0];
    MET = 0;
    FillHistosBasic(0, 1);

    //Int_t reject[5] = {0,0,0,0,0};

    //cout<<"event #"<<nEvents[0]<<endl;
    if (nEvents[0] == 1) weighter->SetDataBit(isRealData);

    //if(nEvents[0]>100) Abort("\t\t  ** 200 EVENTS PASSED, FINISH   ** ");

    /*
    for(Int_t ev=0; ev<10;ev++)
      {
	if (eventNumber==hisEVTS[ev]){
	  cout<<eventNumber<<"  Found an event in  beginnings!"<<endl;
	  break;
	}
      }
    */

    //cout<<"dbg"<<endl;    
    if (nEvents[0] % (int)5e4 == 0) cout<<nEvents[3]<<" events passed of "<<nEvents[0]<<" checked! (at Z-peak cut)"<<endl;

    //if (selection == "gamma"   && (eventNumber % 3) != 0) return kTRUE;
    //if (selection == "eGamma"  && (eventNumber % 3) != 1) return kTRUE;
    //if (selection == "muGamma" && (eventNumber % 3) != 2) return kTRUE;

    //////////////////
    //Trigger status//
    //////////////////

    /*
      for(int i = 0; i < 32; ++i) {
      unsigned int iHLT = 0x0; 
      iHLT = 0x01 << i;  
      if ((triggerStatus & iHLT) == iHLT) h1_triggerStatus->Fill(i+1);  
    } 
    */
    
    Bool_t triggerPass   = triggerSelector->SelectTriggers(triggerStatus, hltPrescale);
    if (!triggerPass) return kTRUE;
    // Double electron workaround.  Gets rid of hopelessly prescaled events of July 20-26, 2011
    if (selection == "electron" && (runNumber > 171200 && runNumber < 171600)) return kTRUE;

    Int_t  eventPrescale = triggerSelector->GetEventPrescale();

    CountEvents(1);
    ++nEventsWeighted[1];
    FillHistosBasic(1, 1);

    //cout<<"dbg"<<endl;    

    ////////////////////////////
    //Check the event vertices//
    ////////////////////////////

    //  if (!isRealData) h1_simVertexMult->Fill(nPUVertices);
    //---------------------------------------------------
    //--- Primary vertex filter -------------------------
    //-------------------------------------------------
    // if (primaryVtx->GetSize() == 0) return kTRUE;
    nVtxTotal = primaryVtx->GetSize();  

    Int_t PVind[2] = {-1, -1}; //indecies for two best PVs
    Bool_t vertexFilter = kFALSE;
    nVtx = 0;
    for (Int_t vv=0; vv<primaryVtx->GetSize(); vv++)
      {
	TCPrimaryVtx* pVtx = (TCPrimaryVtx*) primaryVtx->At(vv);
	//if(!pVtx->isValid()) cout<<"PV is not valid"<<endl;
	//if(pVtx->isFake()) cout<<"PV is fake"<<endl;

	if( 
	   //!pVtx->IsFake() && 
	    pVtx->NDof() > cut_vndof                      //4
	    && fabs(pVtx->z()) <= cut_vz      //15
	    && fabs(pVtx->Perp()) <= cut_vd0   //Not Mag()!; 2
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
    /*
    for(Int_t ev=0; ev<10;ev++)
      {
	if (eventNumber==hisEVTS[ev]){
	  cout<<eventNumber<<"  After vertex filter"<<endl;
	  break;
	}
      }
    */
    TVector3* pvPosition = new TVector3();
    pvPosition = mainPrimaryVertex;
    

    ///////////////
    // electrons //
    ///////////////
    vector<TCPhysObject> looseLeptons;
    vector<TCPhysObject> electrons;
    //int eleCount = 0;

    for (int i = 0; i <  recoElectrons->GetSize(); ++i) {
        TCElectron* thisElec = (TCElectron*) recoElectrons->At(i);    

        if (fabs(thisElec->Eta()) > 2.5) continue;

        float eleISOendcap = (thisElec->IsoMap("SumPt_R03") + thisElec->IsoMap("EmIso_R03") + thisElec->IsoMap("HadIso_R03") - rhoFactor*TMath::Pi()*0.09)/thisElec->Pt(); 
        float eleISObarrel = (thisElec->IsoMap("SumPt_R03") +  TMath::Max(0.0, thisElec->IsoMap("EmIso_R03")-1.0) + thisElec->IsoMap("EmIso_R03") + thisElec->IsoMap("HadIso_R03") - rhoFactor*TMath::Pi()*0.09)/thisElec->Pt(); 
        bool  elecPass = false;

        if (
                (fabs(thisElec->Eta()) < 1.442     
                 && eleISObarrel                        < 0.1
                 && thisElec->SigmaIetaIeta()           < 0.01  
                 && fabs(thisElec->DphiSuperCluster())  < 0.06  
                 && fabs(thisElec->DetaSuperCluster())  < 0.004 
                 && thisElec->HadOverEm()               < 0.12      
                ) ||
                (fabs(thisElec->Eta()) >  1.556  
                 && eleISOendcap                        < 0.1 
                 && thisElec->SigmaIetaIeta()           < 0.03  
                 && fabs(thisElec->DphiSuperCluster())  < 0.03  
                 && fabs(thisElec->DetaSuperCluster())  < 0.007 
                 && thisElec->HadOverEm()               < 0.1
                )) elecPass = true; 

        if (elecPass) { 
            if (
                    thisElec->Pt() > 10.0
                    && thisElec->PassConversion(80)
                    && fabs(thisElec->Dxy(pvPosition)) < 0.02
                    && fabs(thisElec->Dz(pvPosition))  < 0.1
               ) electrons.push_back(*thisElec);			
        } else if (
		   (
		    (fabs(thisElec->Eta()) < 1.499
		     && thisElec->SigmaIetaIeta()           < 0.01
		     && fabs(thisElec->DphiSuperCluster())  < 0.8
		     && fabs(thisElec->DetaSuperCluster())  < 0.007
		     && thisElec->HadOverEm()               < 0.15
		     ) ||
		    (fabs(thisElec->Eta()) >  1.499
		     && thisElec->SigmaIetaIeta()           < 0.03
		     && fabs(thisElec->DphiSuperCluster())  < 0.7
		     && fabs(thisElec->DetaSuperCluster())  < 0.01
		     && thisElec->HadOverEm()               < 0.99 //noc cut!
		     )
		    )
		   && thisElec->Pt() > 10
		   && thisElec->PassConversion(95)
		   && fabs(thisElec->Dxy(pvPosition)) < 0.04
		   && fabs(thisElec->Dz(pvPosition)) < 0.2
		   ) looseLeptons.push_back(*thisElec);
	
    } 

    sort(electrons.begin(), electrons.end(), P4SortCondition);

    ///////////
    // muons //
    ///////////

    vector<TCPhysObject> muons;
    Int_t softMuons = 0;
    for (int i = 0; i < recoMuons->GetSize(); ++ i) {
        TCMuon* thisMuon = (TCMuon*) recoMuons->At(i);    

	/*
	for(Int_t ev=0; ev<10;ev++)
	  {
	    if (eventNumber==hisEVTS[ev]){
	      cout<<i<<"  A muon with pt="<<thisMuon->Pt()<<"  eta="<<thisMuon->Eta()<<"  charge="<<thisMuon->Charge()<<endl;
	      cout<<"\t\t  thisMuon->NumberOfValidTrackerHits() = "<<thisMuon->NumberOfValidTrackerHits()<<endl;
	      cout<<"\t\t  fabs(thisMuon->Dxy(pvPosition))   "<<fabs(thisMuon->Dxy(pvPosition))<<endl;
	      cout<<"\t\t  fabs(thisMuon->Dz(pvPosition))   "<<fabs(thisMuon->Dz(pvPosition))<<endl;
	      cout<<"\t\t  rhoFactor   "<<rhoFactor<<endl;
	      cout<<"\t\t  (thisMuon->TrkIso() + thisMuon->HadIso() + thisMuon->EmIso() - rhoFactor*TMath::Pi()*0.09)/thisMuon->Pt()   "<<(thisMuon->TrkIso() + thisMuon->HadIso() + thisMuon->EmIso() - rhoFactor*TMath::Pi()*0.09)/thisMuon->Pt()<<endl;
	      cout<<"\t\t  thisMuon->PtError()/thisMuon->Pt()   "<<thisMuon->PtError()/thisMuon->Pt()<<endl;
	      cout<<"\t\t  thisMuon->NumberOfValidTrackerHits()   "<<thisMuon->NumberOfValidTrackerHits()<<endl;
	      cout<<"\t\t  thisMuon->NumberOfValidPixelHits()   "<<thisMuon->NumberOfValidPixelHits()<<endl;
	      cout<<"\t\t  thisMuon->NumberOfMatches()   "<<thisMuon->NumberOfMatches()<<endl;
	      cout<<"\t\t  thisMuon->NormalizedChi2()   "<<thisMuon->NormalizedChi2()<<endl;
	      //cout<<"\t\t     "<<<<endl;
	      //cout<<"\t\t     "<<<<endl;

	      break;
	    }
      }
	    */

	
        if (!(fabs(thisMuon->Eta()) < 2.4	
	      && thisMuon->IsGLB() 
	      && thisMuon->IsTRK()
	      )) continue;
	
        if (
	    thisMuon->Pt() > 3.0
	    && thisMuon->NumberOfValidTrackerHits() > 10
	    && fabs(thisMuon->Dxy(pvPosition)) < 0.3
	    && fabs(thisMuon->Dz(pvPosition))  < 0.30 
	    && (thisMuon->IsoMap("SumPt_R03") + thisMuon->IsoMap("EmIso_R03") + thisMuon->IsoMap("HadIso_R03")  - rhoFactor*TMath::Pi()*0.09)/thisMuon->Pt() < 0.15
	    
	    ) {
	  softMuons++;
	  
	  if(thisMuon->Pt()>10
	     && thisMuon->PtError()/thisMuon->Pt() < 0.1
	     && thisMuon->NumberOfValidMuonHits()    > 0
	     && thisMuon->NumberOfValidTrackerHits() > 10
	     && thisMuon->NumberOfValidPixelHits()   > 0
	     && thisMuon->NumberOfMatches() > 1
	     && thisMuon->NormalizedChi2()  < 10.0
	     && fabs(thisMuon->Dxy(pvPosition)) < 0.02
	     && fabs(thisMuon->Dz(pvPosition))  < 0.1 
	     && (thisMuon->IsoMap("SumPt_R03") + thisMuon->IsoMap("EmIso_R03") + thisMuon->IsoMap("HadIso_R03") - rhoFactor*TMath::Pi()*0.09)/thisMuon->Pt() < 0.15
	     )	  muons.push_back(*thisMuon);
     
	} else if (
		   thisMuon->Pt() > 10.
		   //&& thisMuon->Pt() <= muPtCut[0] 
		   && thisMuon->NumberOfValidMuonHits()    > 0
		   && thisMuon->NumberOfValidTrackerHits() > 2 
		   && thisMuon->NumberOfValidPixelHits()   > 0
		   && thisMuon->NumberOfMatches() > 0
		   && thisMuon->NormalizedChi2()  < 9999.
		   && fabs(thisMuon->Dxy(pvPosition)) < 0.2
		   && fabs(thisMuon->Dz(pvPosition))  < 0.2
		   //&& (thisMuon->TrkIso() + thisMuon->HadIso() + thisMuon->EmIso() - rhoFactor*TMath::Pi()*0.09)/thisMuon->Pt() < 0.15
		   ) {
	  looseLeptons.push_back(*thisMuon);
        }
    } 

    sort(muons.begin(), muons.end(), P4SortCondition);
    sort(looseLeptons.begin(), looseLeptons.end(), P4SortCondition);

    /////////////
    // photons //
    /////////////
    nGamma = 0;
    vector<TCPhysObject> photons; 
    if (selection == "eGamma" || selection == "muGamma" || selection == "gamma") {
        for (Int_t i = 0; i < recoPhotons->GetSize(); ++i)
        {
	  TCPhoton* thisPhoton = (TCPhoton*) recoPhotons->At(i);
	  //cout<<i+1<<" dbg  pt: "<<thisPhoton->Pt()<<"  eta: "<< thisPhoton->Eta()<<endl;
	  if (thisPhoton->Pt() < photonPtCut[0] || thisPhoton->Pt() > photonPtCut[1]) continue;	  
	  if ( 
	      thisPhoton->Pt() >10.0
	      //&& thisPhoton->EmIso()         < (4.2 + 0.006 *thisPhoton->Pt())
	      //&& thisPhoton->HadIso()        < (2.2 + 0.0025*thisPhoton->Pt())
	      //&& thisPhoton->TrkIso()        < (2.0 + 0.001 *thisPhoton->Pt())
	      && thisPhoton->HadOverEm()     < 0.05 
	      && thisPhoton->SigmaIEtaIEta() > 0.001
	      && thisPhoton->SigmaIEtaIEta() < 0.013
	      && thisPhoton->TrackVeto() == 0
	      && fabs(thisPhoton->Eta()) < 1.442
	       )  
	    {
	      photons.push_back(*thisPhoton);

	      //for gamma multiplicity plot:
	      if(thisPhoton->Pt() > 25) nGamma++;	
	    }

	}
	sort(photons.begin(), photons.end(), P4SortCondition);
	  
        if (photons.size() > 0 && eventPrescale > 1) {
	  if (!triggerSelector->PhotonTriggerBins(photons[0].Pt(), true)) return kTRUE;
        }
    }


    //////////
    // Jets //
    //////////

    vector<TLorentzVector> jetP4, bJetP4, softJetP4, bSSVJetP4, testJetP4;
    TLorentzVector sumJetP4(0,0,0,0);
    int jetCount = 0;
    float fwdJetSumPt = 0.;
    int fwdJetCount = 0;

    nJets=0, nJetsEta24=0, nJetsB=0, nJetsBssv=0;
    for (int i = 0; i < recoJets->GetSize(); ++i) {
      TCJet* thisJet = (TCJet*) recoJets->At(i);     

      // Prevent lepton overlap //
      bool leptonOverlap = false;
      for (int j = 0; j < (int)muons.size(); ++j)     if (thisJet->DeltaR(muons[j]) < 0.4) leptonOverlap = true;
      for (int j = 0; j < (int)electrons.size(); ++j) if (thisJet->DeltaR(electrons[j]) < 0.4) leptonOverlap = true;
      
      if (selection == "eGamma" || selection == "muGamma") {
	for (int j = 0; j < (int)photons.size(); ++j) if (thisJet->DeltaR(photons[j]) < 0.4) leptonOverlap = true;
      }
      if (leptonOverlap) continue;
      
      if (fabs(thisJet->Eta()) < 2.4) {
	if (
	    // Loose ID jets. Those are used fo deltaPhi(Jet, Met) cut
	    thisJet->NumConstit()    > 1
	    && thisJet->NeuHadFrac() < 0.99
	    && thisJet->NeuEmFrac()  < 0.99
	    && thisJet->ChEmFrac()   < 0.99
	    && thisJet->ChHadFrac()  > 0.0
	    && thisJet->NumChPart()  > 0.0

	    //&& thisJet->P4(JC_LVL).Pt() > jetPtCut[0]
	    ) {
	  
	  if (thisJet->Pt() > jetPtCut[0]) {
	    if (thisJet->BDiscriminatorMap("JPB") > 0.275 && thisJet->Pt() > bJetPtCut) 
	      //if (thisJet->BDiscrTCHE() > 2. && thisJet->P4(JC_LVL).Pt() > bJetPtCut) 
	      bJetP4.push_back(*thisJet);
	    
	    jetP4.push_back(*thisJet);
	    sumJetP4 += TLorentzVector(*thisJet);
	    ++jetCount;
	    ++nJetsEta24;
	  }
	  else if (thisJet->Pt() > jetPtCut[1]) softJetP4.push_back(*thisJet);
	}
    	
      } else if (fabs(thisJet->Eta()) < 4.9) {
    	if (thisJet->Pt() > jetPtCut[0]
    	    && thisJet->NumConstit()    > 1
    	    && thisJet->NeuHadFrac() < 0.99
    	    && thisJet->NeuEmFrac()  < 0.99
    	    ) {
     	  jetP4.push_back(*thisJet); 
    	  sumJetP4 += TLorentzVector(*thisJet);
    	  fwdJetSumPt += thisJet->Pt();
    	  ++fwdJetCount;
    	  ++jetCount;
    	} else if (thisJet->Pt() > jetPtCut[1]) softJetP4.push_back(*thisJet);
      }
    }
    
    sort(jetP4.begin(), jetP4.end(), P4SortCondition);
    sort(bJetP4.begin(), bJetP4.end(), P4SortCondition);
    sort(bSSVJetP4.begin(), bSSVJetP4.end(), P4SortCondition);
    
    nJets15   = softJetP4.size();
    nJets     = jetP4.size();
    nJetsB    = bJetP4.size();
    nJetsBssv = bSSVJetP4.size();
    nJetsB25=0; nJetsB30=0;
    for (Int_t i = 0; i < (Int_t)bJetP4.size(); ++i) {
      if(bJetP4[i].Pt() >25)  nJetsB25++;
      if(bJetP4[i].Pt() >30)  nJetsB30++;
    }

    /////////
    // MET //
    /////////

    TCMET* met = (TCMET*) recoMET;
    TLorentzVector metP4, met1P4, reducedMet1P4, reducedMet2P4;
    metP4.SetPtEtaPhiE(met->Mod(), 0, met->Phi(), met->Mod());
    //met1P4.SetPtEtaPhiE(met->CorrectedMet(), 0, met->CorrectedPhi(), met->CorrectedMet());
    ////////////////////////
    // Analysis selection //
    ////////////////////////

    TLorentzVector ZP4;
    MET       = met->Mod();
    pfMET     = met->Mod();
    //pfMET1    = met->CorrectedMet();
    puCorrMET = PUCorrectedMET(pfMET, nVtx, "pfMet", isRealData);

    Int_t ch1=0, ch2=0;
    TLorentzVector Lepton1(0.,0.,0.,0.), Lepton2(0.,0.,0.,0.);
    if (selection == "electron") {
        /////////////////////
        // 2 good elctrons //
        /////////////////////

        if (electrons.size() < 2) return kTRUE;
    	//opposite charge requirement
    	//if (electrons[0].Charge() == electrons[1].Charge()) return kTRUE;

        ZP4           = electrons[0] + electrons[1];
        reducedMet1P4 = GetReducedMET(sumJetP4,  electrons[0], electrons[1], metP4, 1);
        reducedMet2P4 = GetReducedMET(sumJetP4,  electrons[0], electrons[1], metP4, 2);

	ch1=electrons[0].Charge();
	ch2=electrons[1].Charge();

    	Lepton1 = electrons[0];
    	Lepton2 = electrons[1];

    } else if (selection == "muon") {

        //////////////////////////////////////////
        // 2 (oppositely charged) muons  w/ pt>20 //
        //////////////////////////////////////////

      if (muons.size() < 2) return kTRUE;
      /*
      for(Int_t ev=0; ev<10;ev++)
	{
	  if (eventNumber==hisEVTS[ev]){
	    cout<<eventNumber<<"  n Mu =2 cut"<<endl;
	    break;
	  }
	}
      */
      //opposite charge requirement
      //if (muons[0].Charge() == muons[1].Charge()) reject[3]=1;//return kTRUE;

      /*
      for(Int_t ev=0; ev<10;ev++)
	{
	  if (eventNumber==hisEVTS[ev]){
	    cout<<eventNumber<<" Same charge leptons"<<endl;
	    break;
	  }
	}
      */
	ch1 = muons[0].Charge();
	ch2 = muons[1].Charge();

    	Lepton1 = muons[0];
    	Lepton2 = muons[1];

	//Rochester Muon corrections to make em peak at Z mass correctly
	//cout<<"before: "<<Lepton1.Pt()<<endl;
	// commenting this out for synchronization
	Int_t runopt = 0;
	
	if(isRealData){
	  
	  if(period=="2011A") runopt=0;
	  else if(period=="2011B") runopt=1;
	  else Abort("   * The Run period is not found! (2011A or 2011B)");
	  
	  if (ch1==-1)	
	    roch->momcor_data(Lepton1,Lepton2, 1,0,runopt);
	  else
	    roch->momcor_data(Lepton2,Lepton1, 1,0,runopt);
	}
	else {
	  if (ch1==-1)	
	    roch->momcor_mc(Lepton1,Lepton2, 1,0,0);
	  else
	    roch->momcor_mc(Lepton2,Lepton1, 1,0,0);
	}
	//cout<<"after:  "<<Lepton1.Pt()<<endl;
	

	if (Lepton2.Pt()>Lepton1.Pt()){
	  // it could get messed up after Roch corrections
	  TLorentzVector temp = Lepton1;
	  Lepton1 = Lepton2;
	  Lepton2 = temp;
	}
        ZP4           = Lepton1 + Lepton2;
        reducedMet1P4 = GetReducedMET(sumJetP4, Lepton1, Lepton2, metP4, 1);
        reducedMet2P4 = GetReducedMET(sumJetP4, Lepton1, Lepton2, metP4, 2);

    } else if (selection == "eGamma" || selection == "muGamma" || selection == "gamma") {

        /////////////////////////////////////
        // Photons in place of Z candidate //
        /////////////////////////////////////

        if (photons.size() < 1) return kTRUE;
    	float photonMass = weighter->GetPhotonMass(); 
    	ZP4.SetPtEtaPhiM(photons[0].Pt(), photons[0].Eta(), photons[0].Phi(), photonMass);
        reducedMet1P4 = GetReducedMET(sumJetP4,  photons[0], TLorentzVector(0, 0, 0, 0), metP4, 1);
    	reducedMet2P4 = GetReducedMET(sumJetP4,  photons[0], TLorentzVector(0, 0, 0, 0), metP4, 2);

    	//needed for weights
    	Lepton1 = ZP4;  
    	Lepton2.SetXYZT(0,0,0,0);

    } else {
        return kTRUE;
    }
    //Event Weight
    float eventWeight = 1.0;\
    // Doesn't work for now. Need to check with Brian/Nate
    eventWeight   = weighter->GetTotalWeight(nPUVertices, jetP4.size(), Lepton1, Lepton2);

    if ((selection == "gamma" || selection == "muGamma" || selection == "eGamma") && isRealData) {
      //cout<<"Weight, before/after prescale   "<<eventWeight;
      eventWeight  *= eventPrescale;
      //cout<<"\t\t"<<eventWeight<<endl;
      //cout<<"nPVs = "<<primaryVtx->GetSize()<<"\t gamma pt= "<<ZP4.Pt()<<endl;
    }

    //DeltaPhi
    float deltaPhiJetMET = 2*TMath::Pi();
    if (jetP4.size() > 0)            deltaPhiJetMET = DeltaPhiJetMET(metP4, jetP4);
    else if (softJetP4.size() > 0)   deltaPhiJetMET = DeltaPhiJetMET(metP4, softJetP4);




    /////////////////// 
    // Gen particles //
    ///////////////////
    
    if (!isRealData && (suffix.Contains("HZZ") || suffix.Contains("HWW") )) {
      
      // NLO k-faktors weights
      for (int i = 0; i < genParticles->GetSize(); ++i) {
    	TCGenParticle* thisParticle = (TCGenParticle*) genParticles->At(i);    
	
    	if (thisParticle->GetPDGId() == 25 && thisParticle->Mother()==25)
	  {
	    //cout<<"   * Found a higgs in evt# = "<<eventNumber<<"   weight before = "<<eventWeight<<endl;
	    //cout<<"Mother:  "<<thisParticle->Mother()<<endl;
	    Int_t i_mass = 0;
	    if (suffix.Contains("ggHZZ") || suffix.Contains("ggHWW")) {
	      i_mass = TString(suffix(5,3)).Atoi();
	      
	      hig_Pt = thisParticle->Pt();
	      hig_M  = thisParticle->M();
	      eventWeight *= weighter->GluGluHiggsWeight(hig_Pt, i_mass);
	      
	      //cout<<i_mass<<" Higgs Pt and Mass: "<<hig_Pt<<"  "<<hig_M<<"   weight after = "<<eventWeight<<endl;
	    }
 
	    else if (suffix.Contains("VBFHZZ")){
	      i_mass = TString(suffix(6,3)).Atoi(); //for vbf name
	    }
	    
	    // Higgs mass lineshapes weights
	    //eventWeight *= weighter->HiggsMassLineShapeWeight(hig_M, i_mass);
	    //cout<<"imass: "<<i_mass<<"   genMass from gen Z's: "<<hig_M<<"   eventWeight= "<<eventWeight<<endl; 
	    
    	 }
      }      
    }
   
    //cout<<"dbg"<<endl;    


   if(selection  == "muon" ||  selection  == "electron")
     if (Lepton1.Pt() < 20.0 || Lepton2.Pt() < 20.0) //reject[3]=1;
       return kTRUE;

   /*
   for(Int_t ev=0; ev<10; ev++)
      {
	if (eventNumber==hisEVTS[ev]){
	  cout<<eventNumber<<"  Veto on pt>20 cut?   "<<reject[3]
	      <<"\t\t lep1: pt"<<Lepton1.Pt()<<"  eta: "<<Lepton1.Eta()<<"  phi: "<<Lepton1.Phi()
	      <<"\t\t lep2: pt"<<Lepton2.Pt()<<"  eta: "<<Lepton2.Eta()<<"  phi: "<<Lepton2.Phi()
	      <<endl;
	  break;
	}

      }
    return kTRUE;
   */
   
    CountEvents(2);
    nEventsWeighted[2] += eventWeight;
    FillHistosBasic(2, eventWeight);

    if(verboseLvl>0)
      PrintOut(2,kTRUE);

    /////////////////////
    // 3rd lepton veto //
    /////////////////////

    if ((electrons.size() > 0 || muons.size() > 2) && selection  == "muon") return kTRUE;
    if ((muons.size() > 0 || electrons.size() > 2) && selection  == "electron") return kTRUE;
    if ((muons.size() > 0 || electrons.size() > 0) && (selection == "eGamma" || selection == "muGamma" || selection == "gamma")) return kTRUE;


    //Loose third lepton Veto:
    if( (selection == "muon" || selection == "electron") && looseLeptons.size()>0) return kTRUE;
    //Soft third muon veto:
    if (softMuons > 2) return kTRUE;


    CountEvents(3);
    nEventsWeighted[3] += eventWeight;
    FillHistosBasic(3, eventWeight);

    ////////////
    // Z mass //
    ////////////

    if (ZP4.M() > zMassCut1[0] && ZP4.M() < zMassCut1[1])
      {
	CountEvents(4);
	nEventsWeighted[4] += eventWeight;
	FillHistosBasic(4, eventWeight);
	//FillHistosFull(4, eventWeight);
      }


    /////////////////////////////////////////////////
    // Variables for FillHistos() function (Andrey)//
    /////////////////////////////////////////////////
    qT      = ZP4.Pt();
    Mll     = ZP4.M(); Mll_EE=0; Mll_EB=0; Mll_EX=0;
    diEta   = ZP4.Eta();
    diPhi   = ZP4.Phi();
    MT      = CalculateTransMass(metP4, ZP4);
    MT1     = CalculateTransMass(metP4, ZP4);
    //METqt   = metP4.Pt()/qT;
    MET_phi = metP4.Phi();
    projMET = projectedMET(MET, MET_phi, Lepton1, Lepton2);
    ZprojMET= ZprojectedMET(MET, MET_phi, Lepton1, Lepton2);
    redMET1 = reducedMet1P4.Pt();
    redMET2 = reducedMet2P4.Pt();
    compMET = (sumJetP4 + ZP4).Pt();

    dPhiMetDiLep = fabs(ZP4.DeltaPhi(metP4));
    metProjOnQt  = metP4.Pt()*cos(metP4.DeltaPhi(ZP4)); 
    metPerpQt    = fabs(metP4.Pt()*sin(metP4.DeltaPhi(ZP4))); 

    metOverQt = MET/qT;

    lep_dPhi = fabs(TVector2::Phi_mpi_pi(lep1_phi - lep2_phi));
    lep_dEta = fabs(lep1_eta - lep2_eta);
    lep_dR   = sqrt(pow(lep_dPhi,2) + pow(lep_dEta,2));
    lep_ptRatio  = lep2_pt/lep1_pt;
    lep_angle    = fabs(Lepton1.Angle(Lepton2.Vect()));
    lep_angleLog = log10(TMath::Pi()-lep_angle);
    dPhiClos1 = deltaPhiJetMET;

    if(PVind[0]!=-1) nDofVtx1 = mainPrimaryVertex->NDof();
    if(PVind[1]!=-1) nDofVtx2 = secondPrimaryVertex->NDof();
  

    //////////////////////////////////
    //  di-Leptons library //////////
    ////////////////////////////////

    if (doZlibrary){
      pair<TLorentzVector, TLorentzVector> myLeptons = zLib->GetLeptonsFromLibrary(ZP4);
      Lepton1 = myLeptons.first;
      Lepton2 = myLeptons.second;

      TLorentzVector diLepton = Lepton1+Lepton2;
      //cout<<"lep1 pt= "<<Lepton1.Pt()<<"     lep2 pt= "<<Lepton2.Pt()<<endl;
      
      //replace the Gamma/Z variables with the di-lepton from the library
      qT    = diLepton.Pt(); 
      Mll   = diLepton.M();
      diEta = diLepton.Eta();
      diPhi = diLepton.Phi();

    }
    //----end of library -------//
    //////////////////////////////

    //-------------------------------------------//
    //--- Variables for FillHistos (continue) ---//
    //-------------------------------------------//
    lep1_eta = Lepton1.Eta();
    lep1_phi = Lepton1.Phi();
    lep1_pt  = Lepton1.Pt();

    lep2_eta = Lepton2.Eta();
    lep2_phi = Lepton2.Phi();
    lep2_pt  = Lepton2.Pt();

    //pfMET_lg = metP4.Pt()*cos(metP4.DeltaPhi(ZP4)); 
    pfMET_recoil = -1*(ZP4+metP4).Pt()*cos((ZP4+metP4).DeltaPhi(ZP4));


    if(selection=="muon" || selection=="electron" || ((selection=="muGamma" || selection=="eGamma") && doZlibrary)){
      if(fabs(Lepton1.Eta()) < 1.444 && fabs(Lepton2.Eta())<1.444) Mll_EB = Mll;
      else  if(fabs(Lepton1.Eta()) > 1.566 && fabs(Lepton2.Eta())>1.566)  Mll_EE = Mll;
      else  Mll_EX = Mll;
    }

    if (jetP4.size()!=0){
      dRjetlep1 =  jetP4[0].DeltaR(Lepton1); 
      dRjetlep2 =  jetP4[0].DeltaR(Lepton2);
      ptLeadJet  = jetP4[0].Pt();
      etaLeadJet = jetP4[0].Eta();
      phiLeadJet = jetP4[0].Phi();

      if (jetP4.size()>1){
	ptTrailJet  = jetP4[1].Pt();
	etaTrailJet = jetP4[1].Eta();
	phiTrailJet = jetP4[1].Phi();
        deltaEtaDiJet = fabs(jetP4[0].Eta() - jetP4[1].Eta());
	massDiJet = (jetP4[0] + jetP4[1]).M();
	zeppDiJetDiLep = ZP4.Eta() - 0.5*(jetP4[0].Eta() + jetP4[1].Eta());

        zeppDiJetDiLep_y = ZP4.Rapidity() - 0.5*(jetP4[0].Rapidity() + jetP4[1].Rapidity());
      }
    }
    if (bJetP4.size()!=0)
      ptLeadBJet  = bJetP4[0].Pt();

    ///-------------------/////


    ////////////////
    // b-jet veto //
    ////////////////
    
    Bool_t passBveto = kFALSE;
    if (nJetsB == 0)
      passBveto = kTRUE;


    if (ZP4.M() < zMassCut2[0] || ZP4.M() > zMassCut2[1]) return kTRUE;

    if (ZP4.Pt() > 30 && MET > 30 && passBveto){
      
      mva_nJets15 = (Float_t)nJets15; // need it to be float for mva input
      mva_weight = eventWeight;

      if(makeMvaTree)
	_mvaTree -> Fill();
      

      if (nJets==0){
	CountEvents(18);
	nEventsWeighted[18] += eventWeight;
	FillHistosFull(18, eventWeight);
	FillHistosBasic(18, eventWeight);
      }
      else if(nJets==1){
	CountEvents(19);
	nEventsWeighted[19] += eventWeight;
	FillHistosFull(19, eventWeight);
	FillHistosBasic(19, eventWeight);
      }
      else{
	CountEvents(20);
	nEventsWeighted[20] += eventWeight;
	FillHistosFull(20, eventWeight);
	FillHistosBasic(20, eventWeight);
      }
      // -------------------------- MVA stuff -------------------------------------------
#ifdef USE_MVA

      // ** Only do this for Odd event numbers because even numbers were used for training!!

      Bool_t skipEvent = kFALSE;     
      if(selection=="muon" && !suffix.Contains("DATA")){
	eventWeight*=2.0;
	if (eventNumber%2==0) skipEvent = kTRUE;   
      }

      if(!skipEvent){
	//Apply MVA cuts 
	
	// *additional* cuts for the BDT pre-selection on top of the ones already applied
	// assign values to the variables used in the BDT      
	// get the MVA discriminators for the considered methods
      
	//             [higgs mass points]
	Bool_t passBdtCut[N_HIGGS_MASSES] = {kFALSE, kFALSE, kFALSE};
	Bool_t passAllBdtCuts[N_HIGGS_MASSES] = {kTRUE, kTRUE, kTRUE};
	
	
	//                    [mva method][higgs mass point]
	Float_t tmvaValue[N_DISCR_METHODS][N_HIGGS_MASSES] = {{-999.0, -999.0, -999.00}};
	
	Int_t discrJetMultiInd = D_0J;
	if (nJets == 1) discrJetMultiInd = D_1J;
	else if (nJets > 1) discrJetMultiInd = D_GT1J;
	
	int discr = BDT; // use only this one for now
	//int discr = MLPBNN; // use only this one for now
	
	for (int mh = 0; mh<N_HIGGS_MASSES; ++mh) {
	  TString label = TString::Format("%s_%s_MVA_hzz%i_%s", discrMethodName[discr].Data(), discrSampleName.Data(),
					  mvaHiggsMassPoint[mh], discrJetMultiName[discrJetMultiInd].Data());
	  tmvaValue[discr][mh] = tmvaReader[mh][discrJetMultiInd]->EvaluateMVA(label.Data());
	  
	  mva_discr[mh][discrJetMultiInd] -> Fill(tmvaValue[discr][mh], eventWeight);
	  //cout<<"   ** TMVA value = "<<tmvaValue[discr][mh]<<endl;
	  if (tmvaValue[discr][mh] > bdtCut[mh][discrJetMultiInd]) passBdtCut[mh] = kTRUE;
	  passAllBdtCuts[mh] = (passAllBdtCuts[mh] && passBdtCut[mh]);
	}
	
	// here we can count events, fill histograms etc
	if (passAllBdtCuts[0]) //h125
	  //if (passAllBdtCuts[0] && Mll>76 && Mll<106 && fabs(lep1_eta)<1.4 && fabs(lep2_eta)<1.4)
	  {    
	    CountEvents(21);
	    nEventsWeighted[21] += eventWeight;
	    FillHistosFull(21, eventWeight);
	    FillHistosBasic(21, eventWeight);
	    
	    if(nJets==0){
	      CountEvents(22);
	      nEventsWeighted[22] += eventWeight;
	      FillHistosFull(22, eventWeight);
	      FillHistosBasic(22, eventWeight);
	    }
	    else if(nJets==1){
	      CountEvents(23);
	      nEventsWeighted[23] += eventWeight;
	      FillHistosFull(23, eventWeight);
	      FillHistosBasic(23, eventWeight);
	    }
	    else{
	      CountEvents(24);
	      nEventsWeighted[24] += eventWeight;
	      FillHistosFull(24, eventWeight);
	      FillHistosBasic(24, eventWeight);
	    }
	  }
	
	if (passAllBdtCuts[1]) //h200
	  {    
	    CountEvents(25);
	    nEventsWeighted[25] += eventWeight;
	    FillHistosFull(25, eventWeight);
	    FillHistosBasic(25, eventWeight);
	    
	    if(nJets==0){
	      CountEvents(26);
	      nEventsWeighted[26] += eventWeight;
	      FillHistosFull(26, eventWeight);
	      FillHistosBasic(26, eventWeight);
	    }
	    else if(nJets==1){
	      CountEvents(27);
	      nEventsWeighted[27] += eventWeight;
	      FillHistosFull(27, eventWeight);
	      FillHistosBasic(27, eventWeight);
	    }
	    else{
	      CountEvents(28);
	      nEventsWeighted[28] += eventWeight;
	      FillHistosFull(28, eventWeight);
	      FillHistosBasic(28, eventWeight);
	    }
	  }
	
	if (passAllBdtCuts[2])  //h250
	  {    
	    CountEvents(29);
	    nEventsWeighted[29] += eventWeight;
	    FillHistosFull(29, eventWeight);
	    FillHistosBasic(29, eventWeight);
	    if(verboseLvl>0)
	      PrintOut(29,kTRUE);
	    
	    if(nJets==0){
	      CountEvents(30);
	      nEventsWeighted[30] += eventWeight;
	      FillHistosFull(30, eventWeight);
	      FillHistosBasic(30, eventWeight);
	    }
	    else if(nJets==1){
	      CountEvents(31);
	      nEventsWeighted[31] += eventWeight;
	      FillHistosFull(31, eventWeight);
	      FillHistosBasic(31, eventWeight);
	    }
	    else{
	      CountEvents(32);
	      nEventsWeighted[32] += eventWeight;
	      FillHistosFull(32, eventWeight);
	      FillHistosBasic(32, eventWeight);
	    }
	  }
	
      }

      if(selection=="muon" && !suffix.Contains("DATA"))
	eventWeight/=2.0;
      
#endif
      // --------------------- End of MVA stuff -------------------

    }


    if (ZP4.M() < zMassCut1[0] || ZP4.M() > zMassCut1[1]) return kTRUE;

    /*
    ///////////////////////
    // Cosmics rejection //
    No needed woth the cuts on Mll!
    Bool_t isCosmic = kFALSE;
    if (selection=="muon" && muons.size() > 1 && CosmicMuonFilter(muons[0], muons[1])) isCosmic = kTRUE;
    if(verboseLvl>0)
      if (isCosmic) 
	PrintOut(4,kTRUE);
    */


    if (ZP4.Pt() < qtCut) return kTRUE;
    
    CountEvents(5);
    nEventsWeighted[5] += eventWeight;
    FillHistosFull(5, eventWeight);
    FillHistosBasic(5, eventWeight);
   
    if (deltaPhiJetMET < dPhiMinCut) return kTRUE;
    CountEvents(6);
    nEventsWeighted[6] += eventWeight;
    FillHistosFull(6, eventWeight);
    FillHistosBasic(6, eventWeight);
    

    /*
    for(Int_t ev=0; ev<10;ev++)
      {
	if (eventNumber==hisEVTS[ev]){
	  cout<<eventNumber<<"  After cosmic filter"<<endl;
	  break;
	}
      }
    */


    //Data quality
    //if (isRealData && (isNoiseHcal || isScraping || isCSCTightHalo)) return kTRUE;


    if (MET < 70) return kTRUE;
    CountEvents(7);
    nEventsWeighted[7] += eventWeight;
    FillHistosFull(7, eventWeight);
    FillHistosBasic(7, eventWeight);
    
    if(passBveto){
      CountEvents(8);
      nEventsWeighted[8] += eventWeight;
      FillHistosFull(8, eventWeight);
      FillHistosBasic(8, eventWeight);
      
      PrintOut(8);

    }    
        

    //    Bool_t passBasic = kTRUE;

    /////////////////////////////////////////
    // DeltaPhi(MET, jet); MET and MT cuts //
    /////////////////////////////////////////

    if (passBveto &&  MET>metMinCut[0]  && (MT > mtMinCut[0] && MT < mtMaxCut[0]))
      {
	//For Higgs  200 . The cuts need to be defined
	CountEvents(9);
	nEventsWeighted[9] += eventWeight;
	FillHistosFull(9, eventWeight);
	FillHistosBasic(9, eventWeight);
	PrintOut(9);
      }

    if (passBveto &&  MET>metMinCut[1]  && (MT > mtMinCut[1] && MT < mtMaxCut[1]))
      {
	//For Higgs  250
	CountEvents(10);
	nEventsWeighted[10] += eventWeight;
	FillHistosFull(10, eventWeight);
	FillHistosBasic(10, eventWeight);
      }

    if (passBveto &&  MET>metMinCut[2]  && (MT > mtMinCut[2] && MT < mtMaxCut[2]))
      {
	//For Higgs  300
	CountEvents(11);
	nEventsWeighted[11] += eventWeight;
       	//FillHistosFull(11, eventWeight);
	FillHistosBasic(11, eventWeight);
      }

    if (passBveto &&  MET>metMinCut[3]  && (MT > mtMinCut[3] && MT < mtMaxCut[3]))
      {
	//For Higgs  350
	CountEvents(12);
	nEventsWeighted[12] += eventWeight;
       	//FillHistosFull(12, eventWeight);
	FillHistosBasic(12, eventWeight);
      }
    if (passBveto &&  MET>metMinCut[4]  && (MT > mtMinCut[4] && MT < mtMaxCut[4]))
      {
	//For Higgs  400
	CountEvents(13);
	nEventsWeighted[13] += eventWeight;
       	//FillHistosFull(13, eventWeight);
	FillHistosBasic(13, eventWeight);

      }
    if (passBveto &&  MET>metMinCut[5]  && (MT > mtMinCut[5] && MT < mtMaxCut[5]))
      {
	//For Higgs  450
	CountEvents(14);
	nEventsWeighted[14] += eventWeight;
       	//FillHistosFull(14, eventWeight);
    	FillHistosBasic(14, eventWeight);
      }

    if (passBveto &&  MET>metMinCut[6]  && (MT > mtMinCut[6] && MT < mtMaxCut[6]))
      {
	//For Higgs  500
	CountEvents(15);
	nEventsWeighted[15] += eventWeight;
       	//FillHistosFull(15, eventWeight);
    	FillHistosBasic(15, eventWeight);
      }

    if (passBveto &&  MET>metMinCut[7]  && (MT > mtMinCut[7] && MT < mtMaxCut[7]))
      {
	//For Higgs  550
	CountEvents(16);
	nEventsWeighted[16] += eventWeight;
       	//FillHistosFull(16, eventWeight);
    	FillHistosBasic(16, eventWeight);
      }

    if (passBveto &&  MET>metMinCut[8]  && (MT > mtMinCut[8] && MT < mtMaxCut[8]))
      {
	//For Higgs  600
	CountEvents(17);
	nEventsWeighted[17] += eventWeight;
       	//FillHistosFull(17, eventWeight);
    	FillHistosBasic(17, eventWeight);
      }
    //cout<<"dbg END"<<endl;    

    return kTRUE;
}

void higgsAnalyzer::Terminate()
{
  for (int i = 2; i < nC; ++i) fout[i].close();

  cout<<"\nRunning over "<<suffix<<" dataset with "<<selection<<" selection."<<"\n"<<endl;
  cout<<"| CUT DESCRIPTION                  |\t"<< "\t|"<<endl;
  cout<<"| Initial number of events:        |\t"<< nEvents[0]  <<"\t|"<<nEventsWeighted[0]  <<"\t|"<<endl;
  cout<<"| Pass HLT selection:              |\t"<< nEvents[1]  <<"\t|"<<nEventsWeighted[1]  <<"\t|"<<endl;
  cout<<"| Two leptons (id and iso):        |\t"<< nEvents[2]  <<"\t|"<<nEventsWeighted[2]  <<"\t|"<<endl;
  cout<<"| Third lepton veto:Z mass window  |\t"<< nEvents[3]  <<"\t|"<<nEventsWeighted[3]  <<"\t|"<<endl;
  cout<<"| Z mass window                    |\t"<< nEvents[4]  <<"\t|"<<nEventsWeighted[4]  <<"\t|"<<endl;
  cout<<"| Z Qt > 55 :                      |\t"<< nEvents[5]  <<"\t|"<<nEventsWeighted[5]  <<"\t|"<<endl;
  cout<<"| dPi(Jet, Met) > 0.5              |\t"<< nEvents[6]  <<"\t|"<<nEventsWeighted[6]  <<"\t|"<<endl;
  cout<<"| Met > 70                         |\t"<< nEvents[7]  <<"\t|"<<nEventsWeighted[7]  <<"\t|"<<endl;
  cout<<"| b-jet veto                       |\t"<< nEvents[8]  <<"\t|"<<nEventsWeighted[8]  <<"\t|"<<endl;
  cout<<"| H200                |\t"<< nEvents[9]  <<"\t|"<<nEventsWeighted[9] <<"\t|"<<endl;
  cout<<"| H250                |\t"<< nEvents[10] <<"\t|"<<nEventsWeighted[10]  <<"\t|"<<endl;
  cout<<"| H300                |\t"<< nEvents[11] <<"\t|"<<nEventsWeighted[11] <<"\t|"<<endl;
  cout<<"|                     |\t"<< nEvents[12] <<"\t|"<<nEventsWeighted[12] <<"\t|"<<endl;
  cout<<"| H400                |\t"<< nEvents[13] <<"\t|"<<nEventsWeighted[13] <<"\t|"<<endl;
  cout<<"|                     |\t"<< nEvents[14] <<"\t|"<<nEventsWeighted[14] <<"\t|"<<endl;
  cout<<"| H500                |\t"<< nEvents[15] <<"\t|"<<nEventsWeighted[15] <<"\t|"<<endl;
  cout<<"|                     |\t"<< nEvents[16] <<"\t|"<<nEventsWeighted[16] <<"\t|"<<endl;
  cout<<"| H600                |\t"<< nEvents[17] <<"\t|"<<nEventsWeighted[17] <<"\t|"<<endl;
  cout<<"| presel 0j           |\t"<< nEvents[18] <<"\t|"<<nEventsWeighted[18] <<"\t|"<<endl;
  cout<<"| presel 1j           |\t"<< nEvents[19] <<"\t|"<<nEventsWeighted[19] <<"\t|"<<endl;
  cout<<"| presel >1j          |\t"<< nEvents[20] <<"\t|"<<nEventsWeighted[20] <<"\t|"<<endl;
  cout<<"| MVA for 125         |\t"<< nEvents[21] <<"\t|"<<nEventsWeighted[21] <<"\t|"<<endl;
  cout<<"| MVA for 125 0j      |\t"<< nEvents[22] <<"\t|"<<nEventsWeighted[22] <<"\t|"<<endl;
  cout<<"| MVA for 125 1j      |\t"<< nEvents[23] <<"\t|"<<nEventsWeighted[23] <<"\t|"<<endl;
  cout<<"| MVA for 125 1+j     |\t"<< nEvents[24] <<"\t|"<<nEventsWeighted[24] <<"\t|"<<endl;

  cout<<"| MVA for 200         |\t"<< nEvents[25] <<"\t|"<<nEventsWeighted[25] <<"\t|"<<endl;
  cout<<"| MVA for 200 0j      |\t"<< nEvents[26] <<"\t|"<<nEventsWeighted[26] <<"\t|"<<endl;
  cout<<"| MVA for 200 1j      |\t"<< nEvents[27] <<"\t|"<<nEventsWeighted[27] <<"\t|"<<endl;
  cout<<"| MVA for 200 1+j     |\t"<< nEvents[28] <<"\t|"<<nEventsWeighted[28] <<"\t|"<<endl;

  cout<<"| MVA for 250         |\t"<< nEvents[29] <<"\t|"<<nEventsWeighted[29] <<"\t|"<<endl;
  cout<<"| MVA for 250 0j      |\t"<< nEvents[30] <<"\t|"<<nEventsWeighted[30] <<"\t|"<<endl;
  cout<<"| MVA for 250 1j      |\t"<< nEvents[31] <<"\t|"<<nEventsWeighted[31] <<"\t|"<<endl;
  cout<<"| MVA for 250 1+j     |\t"<< nEvents[32] <<"\t|"<<nEventsWeighted[32] <<"\t|"<<endl;

  
  histoFile->Write();
  histoFile->Close();  
  
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
    float dimuonAngle = muon1.Angle(muon2.Vect());
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
    //return fabs(nearestJet.DeltaPhi(metP4));
    return minDeltaPhi;
}

void higgsAnalyzer::PrintOutNoisy(Int_t num){
  nout[num]<<"* "<<runNumber<<"\t\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<MET<<"\t* "
	   <<endl;//       <<isNoiseHcal<<"\t* "<<isDeadEcalCluster<<"\t* "<<isScraping<<"\t* "<<isCSCTightHalo<<"\t* "<<endl;
}

void higgsAnalyzer::PrintOut(Int_t num, Bool_t isShort){
  if(isShort)
    fout[num]<<" "<<runNumber<<" "<<eventNumber<<endl;
  
  else
    fout[num]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<nVtx<<"\t* "<<nJets<<"\t* "
	     <<Mll<<"\t* "<<MT<<"\t* "<<MET<<"\t* "
	     <<qT<<"\t* "<<rhoFactor<<"\t* "<<endl;//<<muCountLoose<<"\t* "<<eleCountLoose<<endl;

}

void higgsAnalyzer::FillHistosBasic(Int_t num, Double_t weight){
  met0_et[num] -> Fill(MET, weight);
  if(num!=0)
    { 
      evt_byCut -> Fill(num, weight);  
      evt_byCut_raw -> Fill(num, 1);  
    }
  evt_weight[num] -> Fill(weight);
}


void higgsAnalyzer::FillHistosFull(Int_t num, Double_t weight){
  
  met1_et[num]      -> Fill(MET, weight);
  met2_et[num]      -> Fill(pfMET1, weight);
  met3_et[num]      -> Fill(puCorrMET, weight);
  met4_et[num]      -> Fill(projMET, weight);
  met5_et[num]      -> Fill(ZprojMET, weight);
  met6_et[num]      -> Fill(redMET1, weight);
  met7_et[num]      -> Fill(redMET2, weight);
  met8_et[num]      -> Fill(compMET, weight);

  met1_overQt[num] -> Fill(metOverQt, weight);
  //met1_et_ovQt[num] -> Fill(MET, METqt, weight);
  met1_phi[num]     -> Fill(MET_phi, weight);

  l1_eta[num]  -> Fill(lep1_eta, weight);
  l1_phi[num]  -> Fill(lep1_phi, weight);
  l1_pt[num]   -> Fill(lep1_pt, weight);

  l2_eta[num]  -> Fill(lep2_eta, weight);
  l2_phi[num]  -> Fill(lep2_phi, weight);
  l2_pt[num]   -> Fill(lep2_pt, weight);

  l0_angle[num]       -> Fill(lep_angle,weight);
  l0_angleLog[num]    -> Fill(lep_angleLog,weight);
  l0_dPhi[num]        -> Fill(lep_dPhi,weight);
  l0_dEta[num]        -> Fill(lep_dEta,weight);
  l0_dR[num]          -> Fill(lep_dR,weight);
  l0_ptRatio[num]     -> Fill(lep_ptRatio,weight);

  Float_t dPhiMetZ =  fabs(TVector2::Phi_mpi_pi(MET_phi-diPhi));
  di_dPhiMet[num]     -> Fill(dPhiMetZ,weight);
  met1_projOnQt[num]  -> Fill(metProjOnQt,weight);
  met1_perpQt[num]    -> Fill(metPerpQt,weight);
  met1_recoil_lg[num] -> Fill(pfMET_recoil,weight);

  //met4_over_qt[num] -> Fill(puCorrMETqt, weight);
  //met4_et_ovQt[num] -> Fill(puCorrMET, puCorrMETqt, weight);
  //met4_puSig[num]   -> Fill(puSigMET, weight);

  mt2[num]      -> Fill(MT, weight);
  //mt2_met2[num] -> Fill(MT, MET, weight);

  //mt2_met3[num] -> Fill(MT, projMET, weight);
  //mtZ_met3[num] -> Fill(MTZ, projMET, weight);


  di_qt[num]   -> Fill(qT, weight);
  di_eta[num]  -> Fill(diEta, weight);
  di_phi[num]  -> Fill(diPhi, weight);
  di_mass[num] -> Fill(Mll, weight);
  if(Mll_EB!=0) di_mass_EB[num]  -> Fill(Mll_EB, weight);
  if(Mll_EE!=0) di_mass_EE[num]  -> Fill(Mll_EE, weight);
  if(Mll_EX!=0) di_mass_EX[num]  -> Fill(Mll_EX, weight);

  jet_N[num]     -> Fill(nJets, weight);
  jet_N15[num]   -> Fill(nJets15, weight);
  jet_N24[num]   -> Fill(nJetsEta24, weight);
  jet_b_N[num]   -> Fill(nJetsB, weight);
  jet_b_Nssv[num]-> Fill(nJetsBssv, weight);
  jet_b_N25[num] -> Fill(nJetsB25, weight);
  jet_b_N30[num] -> Fill(nJetsB30, weight);


  vtx_nPV_tot[num]    -> Fill(nVtxTotal, weight);
  vtx_nPV_raw[num]    -> Fill(nVtx);
  vtx_nPV_weight[num] -> Fill(nVtx, weight);
  vtx_ndof_1[num]     -> Fill(nDofVtx1, weight);
  vtx_ndof_2[num]     -> Fill(nDofVtx2, weight);

  if(isRealData)  run_events[num] -> Fill(runNumber);


  if(nJets>0)
    {
      jet_pt1[num]     -> Fill(ptLeadJet, weight);
      jet_eta1[num]    -> Fill(etaLeadJet, weight);
      jet_phi1[num]    -> Fill(phiLeadJet, weight);

      met1_dPhiClosJet1[num] -> Fill(dPhiClos1, weight);
      //met1_dPhiClosJet2[num] -> Fill(dPhiClos2, weight);
      //met2_dPhiLeadJet1[num] -> Fill(dPhiLead1, weight);
      //met2_dPhiLeadJet2[num] -> Fill(dPhiLead2, weight);
    }

  if(nJets>1)
    {
      jet_pt2[num]    -> Fill(ptTrailJet, weight);
      jet_eta2[num]   -> Fill(etaTrailJet, weight);
      jet_phi2[num]   -> Fill(phiTrailJet, weight);
      jet_b_pt[num]   -> Fill(ptLeadBJet, weight);
      jet_dRlep1[num] -> Fill(dRjetlep1, weight); 
      jet_dRlep2[num] -> Fill(dRjetlep2, weight); 
      
      jet_diM[num]      -> Fill(massDiJet,weight);
      jet_deltaEta[num] -> Fill(deltaEtaDiJet, weight);
      jet_zeppZ[num]    -> Fill(zeppDiJetDiLep, weight);
      jet_zeppZy[num]   -> Fill(zeppDiJetDiLep_y, weight);

    }
      

  // if(selection =="muGamma" || selection =="eGamma" || selection =="gamma"){
  
  ph_nGamma[num] -> Fill(nGamma, weight);


  if (suffix.Contains("ggHZZ") || suffix.Contains("ggHWW")) {
    higgs_pt[num]   -> Fill(hig_Pt);
    higgs_mass[num] -> Fill(hig_M);
    higgs_w_pt[num]   -> Fill(hig_Pt, weight);
    higgs_w_mass[num] -> Fill(hig_M, weight);
 
  }
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
