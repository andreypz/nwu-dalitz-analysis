// $Id: template_higgsAnalyzer.C,v 1.67 2013/04/18 00:13:48 andrey Exp $

#define higgsAnalyzer_cxx
#include "higgsAnalyzer.h"
#include <string>
#include <vector>

using namespace std;

// This needs to be un-commented in order to run MVA cuts
//#define USE_MVA

/////////////////////////////
//Specify parameters here. //
/////////////////////////////

const string  selection      = "SELECTION";
const string  period         = "PERIOD";
const Bool_t  do_same_sign   = "DOSAMESIGN";
const int     JC_LVL         = 0;  //No JEC for Pat jets (they are already applied)
const TString suffix("SUFFIX");

Bool_t makeMvaTree = 0;
const UInt_t verboseLvl  = 1;
const Bool_t doZlibrary  = 0, isFromData=0;

/////////////////
//Analysis cuts//
/////////////////

const Float_t jetPtCut[]     = {30., 15.};
const Float_t bJetPtCut      = 30.;
const Float_t muPtCut[]     = {20.,10.};
const Float_t elePtCut[]    = {20.,10.};
const Float_t photonPtCut[] = {25.,1e9};
const Float_t zMassCut1[]   = {76,106};
const Float_t zMassCut2[]   = {70,106};
const Float_t	qtCut       = 55.;
const Int_t   nJetsCut[]    = {0,99};
Float_t	cut_vz  = 24., cut_vd0 = 2. ,cut_vndof = 4.;	//PV filter cuts

// Cuts for mass points in the PAS. 200,	250,	300,	350 etc
const Float_t	dPhiMinCut = 0.5;
const Float_t metMinCut[]  = {75,       75,	85,	85,	90,	110,	110,	120,	120};
const Float_t mtMinCut[]   = {175,	225,	250,	300,	325,	350,	400,	450,	450};
const Float_t mtMaxCut[]   = {275,	325,	350,	400,	425,	500,	650,	700,	9999};

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
//const UInt_t hisEVTS[] = {888043, 888050,888055, //mu
//			  242, 305966, 58878818, 58878887, 58878917,  //ele
//			  43345};
//Int_t evSize = sizeof(hisEVTS)/sizeof(int);


//const UInt_t hisEVTS[] = {3912507,3537222, 2260583, 4048415, 1972146, 1294493, 2879208, 65882, 696163, 4124071, 3263050, 640272};


void higgsAnalyzer::Begin(TTree * tree)
{

  vector<string>* triggerNames = 0;
  TFile   *inFile         = tree->GetCurrentFile();
  TTree   *jobTree        = (TTree*)inFile->Get("ntupleProducer/jobTree");
  jobTree->SetBranchAddress("triggerNames", &triggerNames);
  jobTree->GetEntry();

  //for (unsigned i = 0; i < triggerNames->size(); ++i) cout << triggerNames->at(i) << endl;

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
  triggerSelector = new TriggerSelector(selection, period, *triggerNames);

  histoFile = new TFile("a_higgsHistograms.root", "RECREATE");

  histoFile->mkdir("gen_b", "gen_b");
  histoFile->mkdir("Histos", "Histos");
  histoFile->mkdir("DataInfo", "DataInfo");
  histoFile->mkdir("Muons", "Muons");
  histoFile->mkdir("Electrons", "Electrons");
  histoFile->mkdir("mvaTree", "mvaTree");
  histoFile->cd("Histos");

  hists = new HistManager(histoFile);

  //In the Title of this histogram I encode the sample!

  if (doZlibrary)  zLib->CreateLibrary();

  for (Int_t n=0; n<nC; n++)
    {
      //cout<<"zeroing the array: "<<n<<endl;
      nEvents[n]=0;
      nEventsWeighted[n]=0;
    }

  if (verboseLvl>0){
    //ffout.open("./events_printout_SUFFIX_final.txt",ofstream::out);
    //ncout.open("./counts_for_tex_SUFFIX.txt",ofstream::out);

    //ffout.precision(4); ffout.setf(ios::fixed, ios::floatfield);
    for(Int_t i=0; i<nC; i++)
      {
        if (i!=2 && i!=29) continue;
        fout[i].open(Form("./events_printout_SUFFIX_%i.txt",i),ofstream::out);
        fout[i].precision(3); fout[i].setf(ios::fixed, ios::floatfield);
        fout[i]<<
          "*********\t*********\t*****\t***\t****\t****\t************\t*********\t********\t*****\t*\n"<<
          "* Run    \t* Event  \t* LS \t* nVtx\t* nJets\t* m(ll)\t*   MT(llvv)\t* PFMET  \t* PT(ll)\t* rho\t*\n"<<
          "*********\t*********\t*****\t****\t****\t***\t************\t*********\t********\t*****\t*\n";
        if (verboseLvl>2){
          nout[i].open(Form("./events_filtered_SUFFIX_%i.txt",i),ofstream::out);
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

  histoFile->cd();
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

      } //mh loop

    }// loop over jet multiplicity bin


#endif
  // ------------------ End of MVA stuff ------------

}

bool higgsAnalyzer::Process(Long64_t entry)
{
  GetEntry(entry);
  CountEvents(0);
  ++nEventsWeighted[0];

  //MET = 0;
  //FillHistosBasic(0, 1); This is moved to Notify()

  //Int_t reject[5] = {0,0,0,0,0};

  //cout<<"event #  "<<eventNumber<<"  **N events** at cut 0: "<<nEvents[0]<<endl;
  if (nEvents[0] == 1) weighter->SetDataBit(isRealData);

  //if(!isRealData)  return kTRUE;

  //if(nEvents[0]>100) Abort("\t\t  ** 200 EVENTS PASSED, FINISH   ** ");

  /*
    for(Int_t ev=0; ev<evSize;ev++)
    {
    if (eventNumber==hisEVTS[ev]){
    cout<<eventNumber<<"  Found an event in  beginnings!"<<endl;
    break;
    }
    }
  */

  if (nEvents[0] % (int)5e4 == 0) cout<<nEvents[3]<<" events passed of "<<nEvents[0]<<" checked! (at Z-peak cut)"<<endl;

  //return kTRUE;

  //if (selection == "gamma"   && (eventNumber % 3) != 0) return kTRUE;
  //if (selection == "eGamma"  && (eventNumber % 3) != 1) return kTRUE;
  //if (selection == "muGamma" && (eventNumber % 3) != 2) return kTRUE;

  //////////////////
  //Trigger status//
  //////////////////
  Bool_t triggerPass  = 1;
  Bool_t isFound  = 0;
  Int_t prescale = 99;
  if (isRealData)
    {
      if (selection=="muon")
        {
          if(runNumber>1 && runNumber <196046)
            triggerSelector->SelectTrigger("HLT_Mu17_Mu8_v", triggerStatus, hltPrescale, isFound, triggerPass, prescale);
          
          else if(runNumber>=196046 && runNumber <9999999)
            triggerSelector->SelectTrigger("HLT_Mu17_Mu8_v", triggerStatus, hltPrescale, isFound, triggerPass, prescale);
          
        }
      if (selection=="electron")
        {
          if(runNumber>1 && runNumber <999999)
            triggerSelector->SelectTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v", triggerStatus, hltPrescale, isFound, triggerPass, prescale);
          
        }

      if(!isFound)
        Abort("TRG **** Warning ***\n The trigger name you specified is not in the list of trigger names");
      else if  (triggerPass && prescale!=1)
        {
          //cout<<"  Caution **: the trigger  passed but it is Prescaled!  with a prescale =  "<<prescale<<endl;
          //triggerPass = kFALSE;
          //if (prescale!=90)
          //cout<<"Run #  "<<runNumber<<"  **N events** at this cut: "<<nEvents[1]<<"  prescale "<<prescale<<endl;
        }
      
      hists->fillProfile(runNumber, prescale, "run_prescale", "run_prescale", 40000, 200000., 240000, 0,100, 1, "DataInfo");

    }


  if (!triggerPass) return kTRUE;
  // Double electron workaround.  Gets rid of hopelessly prescaled events of July 20-26, 2011
  //if (selection == "electron" && (runNumber > 171200 && runNumber < 171600)) return kTRUE;


  Int_t  eventPrescale = 1;//triggerSelector->GetEventPrescale();

  //cout<<"event #  "<<eventNumber<<"  **N events** at cut 1: "<<nEvents[1]<<endl;

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
  //mainPrimaryVertex   = (TCPrimaryVtx*)(primaryVtx->At(0));
  //if (nVtx>1) secondPrimaryVertex = (TCPrimaryVtx*)(primaryVtx->At(1));

  // Apply the PV filters here! -------------
  if (!vertexFilter) return kTRUE;
  //cout<<"event #  "<<eventNumber<<"  **N events** after vertex: "<<nEvents[1]<<endl;

  CountEvents(1);
  ++nEventsWeighted[1];
  FillHistosBasic(1, 1);


  /*
    for(Int_t ev=0; ev<evSize;ev++)
    {
    if (eventNumber==hisEVTS[ev]){
    cout<<eventNumber<<"  After vertex filter"<<endl;
    break;
    }
    }
  */
  TVector3* pvPosition;// = new TVector3();
  pvPosition = mainPrimaryVertex;


  ///////////////
  // electrons //
  ///////////////
  //vector<TCPhysObject> leptons;
  vector<TCPhysObject> looseLeptons;
  vector<TCPhysObject> electrons;
  //int eleCount = 0;

  //if (recoMuons->GetSize() < 2) return kTRUE;
  //CountEvents(4);

  //if (recoElectrons->GetSize()>1)
  //cout<<"Event with at least 2 electrons event #  "<<eventNumber<<endl;
  for (int i = 0; i <  recoElectrons->GetSize(); ++i) {
    TCElectron* thisElec = (TCElectron*) recoElectrons->At(i);

    /*
      for(Int_t ev=0; ev<evSize;ev++)
      {
      if (eventNumber==hisEVTS[ev]){
      cout<<eventNumber<<"  n Ele =2 cut"<<endl;
      //DumpElectronInfo(thisElec, pvPosition, rhoFactor, rho25Factor,rhoMuFactor);
      break;

      }
      }
    */
    if (fabs(thisElec->Eta()) > 2.5) continue;

    if (thisElec->Pt() > 10.0)
      MakeElectronPlots(thisElec, pvPosition);

    if (thisElec->Pt() > 10.0
        && PassElectronIdAndIso(thisElec, elIdAndIsoCutsTight, pvPosition)
        ) electrons.push_back(*thisElec);
    else if (thisElec->Pt() > 10.0
             && PassElectronIdAndIso(thisElec, elIdAndIsoCutsLoose, pvPosition)
             ) looseLeptons.push_back(*thisElec);
  }

  sort(electrons.begin(), electrons.end(), P4SortCondition);


  //!! low mass resonance rejection !!//
  bool lowMassOS = false;

  for (unsigned i = 1; i < electrons.size(); ++i) {
    for (unsigned j = 0; j < i; ++j) {
      if (
          electrons[i].Type() == electrons[j].Type()
          && electrons[i].Charge() != electrons[j].Charge()
          && (electrons[i] + electrons[j]).M() < 12
          ) lowMassOS = true;
    }
  }


  if (electrons.size()==3 && electrons[0].Pt()>20 ){
    if (fabs(electrons[0].Charge() + electrons[1].Charge() + electrons[2].Charge())==1)
      {
        Float_t triM = (electrons[0]+electrons[1]+electrons[2]).M();
        hists->fill1DHist(triM, "ele_triM","M(eee)", 60, 0,300, 1, "Electrons");
        if (!lowMassOS)
          hists->fill1DHist(triM, "ele_triM_noLowMass","M(eee)", 60, 0,300, 1, "Electrons");
      }
  }

  ///////////
  // muons //
  ///////////

  vector<TCPhysObject> muons;
  Int_t softMuons = 0;
  for (int i = 0; i < recoMuons->GetSize(); ++ i) {
    TCMuon* thisMuon = (TCMuon*) recoMuons->At(i);

	if (!(fabs(thisMuon->Eta()) < 2.4
	      && thisMuon->IsGLB()
	      && thisMuon->IsTRK()
	      )) continue;

    if (thisMuon->Pt() > 10.0)
      MakeMuonPlots(thisMuon, pvPosition);

	if (thisMuon->Pt() > 3.0
	    && PassMuonIdAndIso(thisMuon, muIdAndIsoCutsSoft,pvPosition)
	    ) softMuons++;

	if(thisMuon->Pt()>10
	   && PassMuonIdAndIso(thisMuon, muIdAndIsoCutsTight,pvPosition)
	   )	  {
	  muons.push_back(*thisMuon);
	}
	else if (thisMuon->Pt() > 10.
             && PassMuonIdAndIso(thisMuon, muIdAndIsoCutsLoose,pvPosition)
             ) looseLeptons.push_back(*thisMuon);

  }

  sort(muons.begin(), muons.end(), P4SortCondition);
  sort(looseLeptons.begin(), looseLeptons.end(), P4SortCondition);

  //!! low mass resonance rejection !!//
  lowMassOS = false;
  for (unsigned i = 1; i < muons.size(); ++i) {
    for (unsigned j = 0; j < i; ++j) {
      if (
          muons[i].Type() == muons[j].Type()
          && muons[i].Charge() != muons[j].Charge()
          && (muons[i] + muons[j]).M() < 12
          ) lowMassOS = true;
    }
  }


  if (muons.size()==3 && muons[0].Pt()>20 ){
    if (fabs(muons[0].Charge() + muons[1].Charge() + muons[2].Charge())==1)
      {
        Float_t triM = (muons[0]+muons[1]+muons[2]).M();
        hists->fill1DHist(triM, "mu_triM","M(mumumu)", 60, 0,300, 1, "Muons");
        if (!lowMassOS)
          hists->fill1DHist(triM, "mu_triM_noLowMass","M(mumumu)", 60, 0,300, 1, "Muons");
      }
  }

  /////////////
  // photons //
  /////////////
  nGamma = 0;
  vector<TCPhysObject> photons;
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

  //if (selection == "eGamma" || selection == "muGamma" || selection == "gamma") {
    //if (photons.size() > 0 && eventPrescale > 1) {
    //if (!triggerSelector->PhotonTriggerBins(photons[0].Pt(), true)) return kTRUE;
    //}
  //}


  //////////
  // Jets //
  //////////

  vector<TLorentzVector> jetP4, bJetP4, softJetP4, bSSVJetP4, testJetP4;
  TLorentzVector sumJetP4(0,0,0,0);

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
          && thisJet->NeuHadFrac() < 0.90
          && thisJet->NeuEmFrac()  < 0.90
          && thisJet->ChEmFrac()   < 0.99
          && thisJet->ChHadFrac()  > 0.0
          && thisJet->NumChPart()  > 0.0

          //&& thisJet->P4(JC_LVL).Pt() > jetPtCut[0]
          ) {

        //cout<<"DBG jets  pt="<<thisJet->Pt()<<endl;
        if (thisJet->Pt() > jetPtCut[0]) {
          if (thisJet->BDiscriminatorMap("JPB") > 0.275 && thisJet->Pt() > bJetPtCut)
            //if (thisJet->BDiscrTCHE() > 2. && thisJet->P4(JC_LVL).Pt() > bJetPtCut)
            bJetP4.push_back(*thisJet);

          jetP4.push_back(*thisJet);
          sumJetP4 += TLorentzVector(*thisJet);
          ++nJetsEta24;
        }
        else if (thisJet->Pt() > jetPtCut[1]) softJetP4.push_back(*thisJet);
      }

    } else if (fabs(thisJet->Eta()) < 4.9) {
      if (thisJet->NumConstit()    > 1
          && thisJet->NeuHadFrac() < 0.90
          && thisJet->NeuEmFrac()  < 0.90
          ) {
        if (thisJet->Pt() > jetPtCut[0])
          {
            jetP4.push_back(*thisJet);
            sumJetP4 += TLorentzVector(*thisJet);
            fwdJetSumPt += thisJet->Pt();
            ++fwdJetCount;
          }
        else if (thisJet->Pt() > jetPtCut[1]) softJetP4.push_back(*thisJet);
      }
    }
  }

  sort(jetP4.begin(), jetP4.end(), P4SortCondition);
  sort(bJetP4.begin(), bJetP4.end(), P4SortCondition);
  //sort(bSSVJetP4.begin(), bSSVJetP4.end(), P4SortCondition);


  nJets15   = softJetP4.size();
  nJets     = jetP4.size();
  nJetsB    = bJetP4.size();
  //nJetsBssv = bSSVJetP4.size();
  nJetsB25=0; nJetsB30=0;
  for (Int_t i = 0; i < (Int_t)bJetP4.size(); ++i) {
    if(bJetP4[i].Pt() > 25)  nJetsB25++;
    if(bJetP4[i].Pt() > 30)  nJetsB30++;
  }

  /////////
  // MET //
  /////////

  TCMET* met = (TCMET*) recoMET;
  TLorentzVector metP4, met1P4;//, reducedMet1P4, reducedMet2P4;
  metP4.SetPtEtaPhiE(met->Mod(), 0, met->Phi(), met->Mod());
  //met1P4.SetPtEtaPhiE(met->CorrectedMet(), 0, met->CorrectedPhi(), met->CorrectedMet());

  TLorentzVector ZP4;
  MET       = met->Mod();
  MET_sumEt = met->SumEt();
  pfMET     = met->Mod();
  //pfMET1    = met->CorrectedMet();
  //puCorrMET = PUCorrectedMET(pfMET, nVtx, "pfMet", isRealData);


  Int_t ch1=0, ch2=0;
  TLorentzVector Lepton1(0.,0.,0.,0.), Lepton2(0.,0.,0.,0.);
  if (selection == "electron") {

    //Find out if electrons reconstracted while it's actually a muon!
    //UInt_t initsize =     electrons.size();
    //if(initsize>0)
    //cout<<"DBG  electrons   size = "<<electrons.size()<<endl;
    vector<TCPhysObject>::iterator it = electrons.begin();
    for (; it!= electrons.end(); )
      {
        vector<TCPhysObject>::iterator jt = muons.begin();

        Bool_t isErased = kFALSE;
        for (; jt!= muons.end(); )
          {
            Float_t dR = (*it).DeltaR(*jt);
            jt++;
            hists->fill1DHist(dR, "ele_dRmu","dR(ele, mu)", 40, 0,4, 1, "Electrons");
            //cout<<"DBG  erasing electrons?  dR = "<<dR<<endl;

            if (dR<0.1)
              {
                //cout<<" yes"<<endl;
                // it = electrons.erase(it); //This is actually a muon faking an electron! Remove it.
                // test isErased = kTRUE;
                break;
              }
          }
        if (!isErased)
          ++it; //if it's not erased we need to increment the iterator, otherwise it's infinite loop...
      }

    //if(initsize>electrons.size())
    //cout<<"DBG  electrons   size = "<<electrons.size()<<"   init size = "<<initsize<<endl;


    /////////////////////
    // 2 good elctrons //
    /////////////////////


    if (electrons.size() < 2) return kTRUE;
    //CountEvents(5);

	/*
      for(Int_t ev=0; ev<evSize;ev++)
	  {
      if (eventNumber==hisEVTS[ev]){
      cout<<eventNumber<<" Same charge leptons"<<endl;
      break;
      }
	  }
	*/
    //opposite charge requirement

    Int_t sign = electrons[0].Charge()*electrons[1].Charge();
    if (do_same_sign && sign==-1) return kTRUE;
    if (!do_same_sign && sign==1) return kTRUE;


    ZP4           = electrons[0] + electrons[1];
    //reducedMet1P4 = GetReducedMET(sumJetP4,  electrons[0], electrons[1], metP4, 1);
    //reducedMet2P4 = GetReducedMET(sumJetP4,  electrons[0], electrons[1], metP4, 2);

	ch1=electrons[0].Charge();
	ch2=electrons[1].Charge();

    Lepton1 = electrons[0];
    Lepton2 = electrons[1];

  } else if (selection == "muon") {

    //////////////////////////////////////////
    // 2 (oppositely charged) muons  w/ pt>20 //
    //////////////////////////////////////////

    if (muons.size() < 2) return kTRUE;
    //opposite charge requirement
    //if (muons[0].Charge() == muons[1].Charge()) reject[3]=1;//return kTRUE;

    Int_t sign = muons[0].Charge()*muons[1].Charge();
    if (do_same_sign && sign==-1) return kTRUE;
    if (!do_same_sign && sign==1) return kTRUE;

	ch1 = muons[0].Charge();
	ch2 = muons[1].Charge();

    Lepton1 = muons[0];
    Lepton2 = muons[1];

	//Rochester Muon corrections to make em peak at Z mass correctly
	//cout<<"before: "<<Lepton1.Pt()<<endl;
	// commenting this out for synchronization
	/*
      Int_t runopt = 0;

      if(isRealData){

	  if(period=="2011A") runopt=0;
	  else if(period=="2011B") runopt=1;
	  else runopt=0;

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
      */

	if (Lepton2.Pt()>Lepton1.Pt()){
	  // it could get messed up after Roch corrections
	  TLorentzVector temp = Lepton1;
	  Lepton1 = Lepton2;
	  Lepton2 = temp;
	}
    ZP4           = Lepton1 + Lepton2;
    //reducedMet1P4 = GetReducedMET(sumJetP4, Lepton1, Lepton2, metP4, 1);
    //reducedMet2P4 = GetReducedMET(sumJetP4, Lepton1, Lepton2, metP4, 2);

  } else if (selection == "eGamma" || selection == "muGamma" || selection == "gamma") {

    /////////////////////////////////////
    // Photons in place of Z candidate //
    /////////////////////////////////////

    if (photons.size() < 1) return kTRUE;
    float photonMass = weighter->GetPhotonMass();
    ZP4.SetPtEtaPhiM(photons[0].Pt(), photons[0].Eta(), photons[0].Phi(), photonMass);
    //reducedMet1P4 = GetReducedMET(sumJetP4,  photons[0], TLorentzVector(0, 0, 0, 0), metP4, 1);
    //reducedMet2P4 = GetReducedMET(sumJetP4,  photons[0], TLorentzVector(0, 0, 0, 0), metP4, 2);

    //needed for weights
    Lepton1 = ZP4;
    Lepton2.SetXYZT(0,0,0,0);

  } else {
    return kTRUE;
  }


  //Event Weight
  float eventWeight = 1.0;
  // Doesn't work for now. Need to check with Brian/Nate
  eventWeight   = weighter->GetTotalWeight(nPUVerticesTrue, jetP4.size(), Lepton1, Lepton2);

  //cout<<"eventWeight "<<eventWeight<<"  nPU = "<<nPUVerticesTrue<<endl;

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

  /* a study of Higgs mass line shape
     if (!isRealData) {


     vector<TCPhysObject> genleps, genneutrinos;
     for (int i = 0; i < genParticles->GetSize(); ++i) {
     TCGenParticle* thisParticle = (TCGenParticle*) genParticles->At(i);
     if (thisParticle->GetPDGId() == 23 && thisParticle->Mother()==23)
     {
     hists->fill1DHist(thisParticle->M(), "gen_Z_M", "Mass of gen Z", 50,0,200, eventWeight, "Histos");

     //cout<<eventNumber<<"  Mass = "<<thisParticle->M()<<"    pt= "<<thisParticle->Pt()<<endl;
     }

     if ((abs(thisParticle->GetPDGId()) == 11 || abs(thisParticle->GetPDGId())==13) && thisParticle->Mother()==23)
     {
     //cout<<eventNumber<<"    Lepton from Z, ID= "<<thisParticle->GetPDGId()<<"   "<<thisParticle->M()<<endl;
     genleps.push_back(*thisParticle);
     }
     if ((abs(thisParticle->GetPDGId()) == 12 || abs(thisParticle->GetPDGId())==14 || abs(thisParticle->GetPDGId())==16) && thisParticle->Mother()==23)
     {
     //cout<<eventNumber<<"    Lepton from Z, ID= "<<thisParticle->GetPDGId()<<"   "<<thisParticle->M()<<endl;
     genneutrinos.push_back(*thisParticle);
     }
     }
     if (genneutrinos.size()==2)
     {
     hists->fill1DHist(genneutrinos[0].Eta(), "gen_nu1_eta", "eta of gen lep1", 40,-5,5, eventWeight, "Histos");
     hists->fill1DHist(genneutrinos[1].Eta(), "gen_nu2_eta", "eta of gen lep1", 40,-5,5, eventWeight, "Histos");
     }
     if (genleps.size()==2)
     {
     Float_t genll = (genleps[0]+genleps[1]).M();
     hists->fill1DHist(genll, "gen_ll_M", "Gen level m(ll)", 50,0,200, eventWeight, "Histos");
     hists->fill1DHist(genleps[0].Eta(), "gen_l1_eta", "eta of gen lep1", 40,-5,5, eventWeight, "Histos");
     hists->fill1DHist(genleps[1].Eta(), "gen_l2_eta", "eta of gen lep2", 40,-5,5, eventWeight, "Histos");
     }
     //else cout<<" wrong size of gen leptons collection  "<<genleps.size()<<endl;

     }
  */

  if (!isRealData && (suffix.Contains("HZZ") || suffix.Contains("HWW") ) && period=="2011") {

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
    if (Lepton1.Pt() < 20.0 || Lepton2.Pt() < 20.0)
      return kTRUE;




  //////////////////////////////////
  //  di-Leptons library //////////
  ////////////////////////////////
  /*
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
    */

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



  /////////////////////////////////////////////////
  // Variables for FillHistos() function (Histos)//
  /////////////////////////////////////////////////
  qT      = ZP4.Pt();
  Mll     = ZP4.M(); Mll_EE=0; Mll_EB=0; Mll_EX=0;
  diEta   = ZP4.Eta();
  diPhi   = ZP4.Phi();
  MT      = CalculateTransMass(metP4, ZP4);
  MT1     = CalculateTransMass(metP4, ZP4);
  MET_phi = metP4.Phi();
  projMET = projectedMET(MET, MET_phi, Lepton1, Lepton2);
  //redMET1 = reducedMet1P4.Pt();
  //redMET2 = reducedMet2P4.Pt();

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

    dRjetDiLept = jetP4[0].DeltaR(ZP4);

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

  if(Mll<10) return kTRUE;

  CountEvents(2);
  nEventsWeighted[2] += eventWeight;
  FillHistosBasic(2, eventWeight);
  FillHistosFull(2, eventWeight);

  //    cout<<"event #  "<<eventNumber<<"  **N events** at cut 2: "<<nEvents[2]<<endl;
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
  //FillHistosFull(3, eventWeight);
  //cout<<"event #  "<<eventNumber<<"  **N events** at cut 3: "<<nEvents[3]<<endl;


  return kTRUE;

  ////////////
  // Z mass //
  ////////////

  if (ZP4.M() > zMassCut1[0] && ZP4.M() < zMassCut1[1])
    {
      CountEvents(4);
      nEventsWeighted[4] += eventWeight;
      FillHistosBasic(4, eventWeight);
      FillHistosFull(4, eventWeight);
      if(NoiseFilters_isScraping
         || NoiseFilters_isNoiseHcalHBHE
         || NoiseFilters_isNoiseHcalLaser
         || NoiseFilters_isNoiseEcalTP
         || NoiseFilters_isNoiseEcalBE
         ){
        hists->fill1DHist(4, "evt_byCut_filters", "Filtered events", nC, 0,nC, 1, "Histos");
        hists->fill1DHist(MET, "met_filters_4", "Met", 80, 0,400, eventWeight, "Histos");
      }
    }


  ////////////////
  // b-jet veto //
  ////////////////

  Bool_t passBveto = kFALSE;
  if (nJetsB == 0)
    passBveto = kTRUE;


  if (ZP4.M() < zMassCut2[0] || ZP4.M() > zMassCut2[1]) return kTRUE;

  if (ZP4.Pt() > 15 && MET > 30 && passBveto){

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

        hists->fill1DHist(tmvaValue[discr][mh], Form("mva_discr_mh%i_%i",mh,discrJetMultiInd),"MVA discriminator output", 50, -0.6,0.3, eventWeight, "Histos");
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
          //if(verboseLvl>0)
          // PrintOut(29,kTRUE);

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


  if (MET < 70) return kTRUE;
  CountEvents(5);
  nEventsWeighted[5] += eventWeight;
  FillHistosFull(5, eventWeight);
  FillHistosBasic(5, eventWeight);

  if(NoiseFilters_isScraping
     || NoiseFilters_isNoiseHcalHBHE
     || NoiseFilters_isNoiseHcalLaser
     || NoiseFilters_isNoiseEcalTP
     || NoiseFilters_isNoiseEcalBE
     ) {
    hists->fill1DHist(5, "evt_byCut_filters", "SUFFIX", nC, 0,nC, 1, "Histos");
    hists->fill1DHist(MET, "met_filters_5", "Met", 80, 0,400, eventWeight, "Histos");
  }
  if (deltaPhiJetMET < dPhiMinCut) return kTRUE;
  CountEvents(6);
  nEventsWeighted[6] += eventWeight;
  FillHistosFull(6, eventWeight);
  FillHistosBasic(6, eventWeight);

  if (ZP4.Pt() < qtCut) return kTRUE;
  CountEvents(7);
  nEventsWeighted[7] += eventWeight;
  FillHistosFull(7, eventWeight);
  FillHistosBasic(7, eventWeight);
  //cout<<"event #  "<<eventNumber<<"  **N events** at cut 7: "<<nEvents[7]<<endl;


  //Data quality
  //if (isRealData && (isNoiseHcal || isScraping || isCSCTightHalo)) return kTRUE;

  // B jets and gen particles study
  if (!isRealData) {

    vector<TCGenParticle> gen_b_quarks;
    //vector<TCPhysObject> gen_b_quarks;
    Int_t Nb_fromW =0;
    for (int i = 0; i < genParticles->GetSize(); ++i) {
      TCGenParticle* thisParticle = (TCGenParticle*) genParticles->At(i);

      //if (abs(thisParticle->GetPDGId()) == 24)
      //   cout<<i<<" b gen id "<<""<<thisParticle->GetPDGId()<<"   mother= "<<thisParticle->Mother()<<"  pt="<<thisParticle->Pt()<<endl;


      if (abs(thisParticle->GetPDGId()) == 5 && abs(thisParticle->Mother())==24)
        {
          cout<<i<<" b gen id "<<""<<thisParticle->GetPDGId()<<"   mother= "<<thisParticle->Mother()<<"  pt="<<thisParticle->Pt()<<endl;
          gen_b_quarks.push_back(*thisParticle);
          Nb_fromW++;
        }

      if (abs(thisParticle->GetPDGId()) == 5 &&  thisParticle->GetStatus()==2)
        {
          //cout<<i<<" b gen id "<<""<<thisParticle->GetPDGId()<<"   status "<<thisParticle->GetStatus()<<"   mother= "<<thisParticle->Mother()<<"  pt="<<thisParticle->Pt()<<endl;
          if (gen_b_quarks.size() == 0)
            gen_b_quarks.push_back(*thisParticle);

          Bool_t matched= kFALSE;
          for (UInt_t j=0; j < gen_b_quarks.size(); j++)
            {
              if ( fabs(gen_b_quarks[j].Pt()-thisParticle->Pt())/thisParticle->Pt() < 0.01 && fabs(gen_b_quarks[j].Eta()-thisParticle->Eta()) < 0.01)
                matched = kTRUE;
            }
          if(!matched)//don't want dulicate particles
            gen_b_quarks.push_back(*thisParticle);

        }
    }

    if(gen_b_quarks.size()>2)
      {
        cout<<gen_b_quarks.size() <<"   ev number = "<<eventNumber<<endl;
        for(UInt_t k=0; k<gen_b_quarks.size(); k++)
          cout<<"ID = "<<gen_b_quarks[k].GetPDGId()
              <<"\t mother= "<<gen_b_quarks[k].Mother()
              <<"\t pt    = "<<gen_b_quarks[k].Pt()
              <<"\t eta   = "<<gen_b_quarks[k].Eta()
              <<"\t phi   = "<<gen_b_quarks[k].Phi()
              <<endl;

      }

    Int_t Nb25=0, Nb30=0;
    for (UInt_t j=0; j < gen_b_quarks.size(); j++)
      {
        if (gen_b_quarks[j].Pt()>25 && fabs(gen_b_quarks[j].Eta())<2.5)
          Nb25++;
        if (gen_b_quarks[j].Pt()>30 && fabs(gen_b_quarks[j].Eta())<2.5)
          Nb30++;

        Float_t pt  = gen_b_quarks[j].Pt();
        Float_t eta = gen_b_quarks[j].Eta();
        hists->fill1DHist(pt,  "b_pt_all", "Pt of all b quarks", 50,0,100, 1,"gen_b");
        hists->fill1DHist(eta, "b_eta_all", "Eta of all b qaurks", 50,-5,5, 1,"gen_b");
      }
    hists->fill1DHist(Nb25, "b_Nb25", "Number of b-quarks with pt>25, eta<2.5", 10,0,10, 1,"gen_b");
    hists->fill1DHist(Nb30, "b_Nb30", "Number of b-quarks with pt>30, eta<2.5", 10,0,10, 1,"gen_b");
    hists->fill1DHist(Nb_fromW, "b_Nb_fromW", "Number of b-quarks in the event from a W", 10,0,10, 1,"gen_b");

    /* need more thought into it
    if(nJets==0)
      {
        hists->fill1DHist(Nb25, "b_Nb25_nj0", "Number of b-quarks with pt>25, eta<2.5", 10,0,10, 1,"gen_b");
        hists->fill1DHist(Nb30, "b_Nb30_nj0", "Number of b-quarks with pt>30, eta<2.5", 10,0,10, 1,"gen_b");
      }
    else
      {
        hists->fill1DHist(Nb25, "b_Nb25_nj1", "Number of b-quarks with pt>25, eta<2.5", 10,0,10, 1,"gen_b");
        hists->fill1DHist(Nb30, "b_Nb30_nj1", "Number of b-quarks with pt>30, eta<2.5", 10,0,10, 1,"gen_b");
      }
    */


  }



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

  passBveto = kTRUE;

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
  cout<<"| Ntuple evts; HLT; vertex         |\t"<< nEvents[1]  <<"\t|"<<nEventsWeighted[1]  <<"\t|"<<endl;
  cout<<"| Two leptons (id and iso):        |\t"<< nEvents[2]  <<"\t|"<<nEventsWeighted[2]  <<"\t|"<<endl;
  cout<<"| Third lepton veto:               |\t"<< nEvents[3]  <<"\t|"<<nEventsWeighted[3]  <<"\t|"<<endl;
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


  histoFile->cd();
  histoFile->Write();
  histoFile->Close();

  //triggerSelector -> Delete();
  //weighter        -> Delete();
  //hists           -> Delete();
  //zLib            -> Delete();
  cout<<" ** End of Analyzer **"<<endl;
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

/*
  bool higgsAnalyzer::CosmicMuonFilter(TCMuon muon1, TCMuon muon2)
  {
  float dimuonAngle = muon1.Angle(muon2.Vect());
  if (TMath::Pi() - dimuonAngle  < 0.05)
  return true;
  else
  return false;
  }
*/

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

void higgsAnalyzer::FillHistosBasic(Int_t num, Double_t weight)
{
  hists->fill1DHist(MET,Form("met0_et_%i",num),"pfMet", 40,0,400, weight, "Histos");
  if(num!=0)
    {
      hists->fill1DHist(num, "evt_byCut","SUFFIX", nC,0,nC, weight, "Histos");
      hists->fill1DHist(num, "evt_byCut_raw", "SUFFIX", nC, 0,nC, 1, "Histos");
    }
  hists->fill1DHist(weight, Form("evt_weight_%i",num), "Weights", 50,0,10, 1, "Histos");
  hists->fill1DHist(nPUVerticesTrue, Form("evt_puTrue_raw_%i",num), "PU true raw", 240, 0,60, 1, "Histos");
  hists->fill1DHist(nPUVerticesTrue, Form("evt_puTrue_%i",num), "PU true", 240, 0,60, weight, "Histos");
  hists->fill1DHist(nPUVertices, Form("evt_pu_%i",num), "PU observed", 240, 0,60, weight, "Histos");
}


void higgsAnalyzer::FillHistosFull(Int_t num, Double_t weight){
  hists->fill1DHist(MET, Form("met1_et_%i",num), "met1_et", 80,0,400, weight, "Histos");
  hists->fill1DHist(MET_sumEt, Form("met1_SumEt_%i",num), "met1_SumEt", 50,0,2000, weight, "Histos");
  hists->fill1DHist(pfMET1, Form("met2_et_%i",num), "met2_et", 80,0,400, weight, "Histos");
  hists->fill1DHist(MET_phi, Form("met1_phi_%i",num), "met1_phi", 40, -TMath::Pi(), TMath::Pi(), weight, "Histos");


  hists->fill1DHist(metOverQt, Form("met1_overQt_%i",num), "met over qT", 40,0,4, weight, "Histos");
  hists->fill1DHist(metProjOnQt, Form("met1_projOnQt_%i",num), "met1_projOnQt", 50, -200,100, weight, "Histos");
  hists->fill1DHist(metPerpQt, Form("met1_perpQt_%i",num), "met1_perpQt", 50, 0,400, weight, "Histos");

  hists->fill1DHist(pfMET_recoil, Form("met1_recoil_lg_%i",num), "met1_recoil_lg", 50, -400,100, weight, "Histos");

  hists->fill1DHist(lep1_eta, Form("l1_eta_%i",num), "l1_eta", 52, -2.6, 2.6,  weight, "Histos");
  hists->fill1DHist(lep1_phi, Form("l1_phi_%i",num), "l1_pi", 50, -TMath::Pi(), TMath::Pi(),  weight, "Histos");
  hists->fill1DHist(lep1_pt,  Form("l1_pt_%i",num),  "l1_pt", 50, 0, 200,  weight, "Histos");
  hists->fill1DHist(lep2_eta, Form("l2_eta_%i",num), "l2_eta", 52, -2.6, 2.6,  weight, "Histos");
  hists->fill1DHist(lep2_phi, Form("l2_phi_%i",num), "l2_pi", 50, -TMath::Pi(), TMath::Pi(),  weight, "Histos");
  hists->fill1DHist(lep2_pt,  Form("l2_pt_%i",num),  "l2_pt", 50, 0, 200,  weight, "Histos");

  hists->fill1DHist(lep_angle,  Form("l0_angle_%i",num), "l0_angle", 50, 0, TMath::Pi(),  weight, "Histos");
  hists->fill1DHist(lep_angleLog,  Form("l0_angleLog_%i",num), "l0_angleLog", 50, -5,1,  weight, "Histos");

  hists->fill1DHist(lep_dPhi, Form("l0_dPhi_%i",num), "l0_dPhi", 50, 0, TMath::Pi(),  weight, "Histos");
  hists->fill1DHist(lep_dEta, Form("l0_dEta_%i",num), "l0_dEta", 50, 0, 4,  weight, "Histos");
  hists->fill1DHist(lep_dR,   Form("l0_dR_%i",num),   "l0_dR",   50, 0, 4,  weight, "Histos");
  hists->fill1DHist(lep_ptRatio, Form("l0_ptRatio_%i",num), "l0_ptRatio", 50, 0, 1,  weight, "Histos");

  hists->fill1DHist(MT, Form("mt_%i",num), "mt", 50, 100,600, weight, "Histos");

  hists->fill1DHist(qT, Form("di_qt_%i",num), "di_qt", 80, 0,400, weight, "Histos");
  hists->fill1DHist(diEta, Form("di_eta_%i",num), "di_eta", 100, -5,5, weight, "Histos");
  hists->fill1DHist(diPhi, Form("di_phi_%i",num), "di_phi", 50, -TMath::Pi(), TMath::Pi(), weight, "Histos");
  hists->fill1DHist(Mll, Form("di_mass_%i",num), "di_mass", 50, 70,120, weight, "Histos");

  if(num==2 || num==3)
    hists->fill1DHist(Mll, Form("di_mass_ext_%i",num), "di_mass_ext", 50, 0,200, weight, "Histos");

  if(Mll_EB!=0)
    hists->fill1DHist(Mll_EB, Form("di_mass_EB_%i",num), "di_mass_EB", 50, 70,120, weight, "Histos");
  if(Mll_EE!=0)
    hists->fill1DHist(Mll_EE, Form("di_mass_EE_%i",num), "di_mass_EE", 50, 70,120, weight, "Histos");
  if(Mll_EX!=0)
    hists->fill1DHist(Mll_EX, Form("di_mass_EX_%i",num), "di_mass_EX", 50, 70,120, weight, "Histos");

  Float_t dPhiMetZ =  fabs(TVector2::Phi_mpi_pi(MET_phi-diPhi));
  hists->fill1DHist(dPhiMetZ, Form("di_dPhiMet_%i",num), "dPhiMet", 50, 0, TMath::Pi(), weight, "Histos");


  hists->fill1DHist(nJets, Form("jet_N_%i",num), "jet_N", 20, 0,20, weight, "Histos");
  hists->fill1DHist(nJets15, Form("jet_N15_%i",num), "jet_N15", 20, 0,20, weight, "Histos");
  hists->fill1DHist(nJetsEta24, Form("jet_N24_%i",num), "jet_N24", 20, 0,20, weight, "Histos");

  hists->fill1DHist(nJetsB, Form("jet_b_N_%i",num), "jet_b_N", 10, 0,10, weight, "Histos");
  hists->fill1DHist(nJetsB25, Form("jet_b_N25_%i",num), "jet_b_N25", 10, 0,10, weight, "Histos");

  if(nJets>0)
    {
      hists->fill1DHist(ptLeadJet, Form("jet_pt1_%i",num), "jet_pt1", 50, 0,400, weight, "Histos");
      hists->fill1DHist(etaLeadJet, Form("jet_eta1_%i",num), "jet_eta1", 50, -5,5, weight, "Histos");
      hists->fill1DHist(phiLeadJet, Form("jet_phi1_%i",num), "jet_phi1", 50,  -TMath::Pi(), TMath::Pi(), weight, "Histos");
      hists->fill1DHist(dRjetlep1, Form("jet_dRlep1_%i",num), "jet_dRlep1", 50, 0,5, weight, "Histos");
      hists->fill1DHist(dRjetlep2, Form("jet_dRlep2_%i",num), "jet_dRlep2", 50, 0,5, weight, "Histos");

      hists->fill1DHist(dPhiClos1, Form("met1_dPhiClosJet1_%i",num), "met1_dPhiClosJet1", 50, 0,TMath::Pi(), weight, "Histos");

      hists->fill1DHist(ptLeadBJet, Form("jet_b_pt_%i",num), "jet_b_pt", 50, 0,400, weight, "Histos");

      hists->fill1DHist(dRjetDiLept, Form("jet_dRDiLept_%i",num), "dR (di-elpt,jet)", 50, 0,5, weight, "Histos");


    }

  if(nJets>1)
    {
      hists->fill1DHist(ptTrailJet, Form("jet_pt2_%i",num), "jet_pt2", 50, 0,400, weight, "Histos");
      hists->fill1DHist(etaTrailJet, Form("jet_eta2_%i",num), "jet_eta2", 50, -5,5, weight, "Histos");
      hists->fill1DHist(phiTrailJet, Form("jet_phi2_%i",num), "jet_phi2", 50,  -TMath::Pi(), TMath::Pi(), weight, "Histos");

      hists->fill1DHist(massDiJet, Form("jet_diM_%i",num), "jet_diM", 50, 0,2000, weight, "Histos");
      hists->fill1DHist(deltaEtaDiJet, Form("jet_deltaEta_%i",num), "jet_deltaEta", 50, 0,10, weight, "Histos");
      hists->fill1DHist(zeppDiJetDiLep, Form("jet_zeppZ_%i",num), "jet_zeppZ", 50, -5,5, weight, "Histos");
      hists->fill1DHist(zeppDiJetDiLep_y, Form("jet_zeppZy_%i",num), "jet_zeppZy", 50, -5,5, weight, "Histos");
    }


  // if(selection =="muGamma" || selection =="eGamma" || selection =="gamma"){

  hists->fill1DHist(nGamma, Form("ph_nGamma_%i",num), "ph_nGamma", 20, 0,20, weight, "Histos");

  if(isRealData)
    hists->fill1DHist(runNumber, Form("run_events_%i",num), "run_events", 40000, 200000., 240000, 1, "DataInfo");


  hists->fill1DHist(nVtxTotal, Form("vtx_nPV_tot_%i",num), "vtx_nPV_tot", 40, 0,40, weight, "Histos");
  hists->fill1DHist(nVtx, Form("vtx_nPV_raw_%i",num), "vtx_nPV_raw", 40, 0,40, 1, "Histos");
  hists->fill1DHist(nVtx, Form("vtx_nPV_weight_%i",num), "vtx_nPV_weight", 40, 0,40, weight, "Histos");

  hists->fill1DHist(nDofVtx1, Form("vtx_ndof1_%i",num), "vtx_ndof1", 50, 0,200, weight, "Histos");
  hists->fill1DHist(nDofVtx1, Form("vtx_ndof2_%i",num), "vtx_ndof2", 50, 0,200, weight, "Histos");
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

void higgsAnalyzer::DumpElectronInfo(TCElectron *lep, TVector3 *vert, Float_t rho44, Float_t rho25, Float_t rhoMu)
{
  cout<<"  A Electron with pt="<<lep->Pt()<<"  eta="<<lep->Eta()<<"  charge="<<lep->Charge()<<endl;
  //cout<<"\t\t  lep->NumberOfValidTrackerHits() = "<<lep->NumberOfValidTrackerHits()<<endl;
  cout<<"\t\t  fabs(lep->Dxy(vert))  = "<<fabs(lep->Dxy(vert))<<endl;
  cout<<"\t\t  fabs(lep->Dz(vert))   = "<<fabs(lep->Dz(vert))<<endl;
  cout<<"\t\t  rho44     = "<<rho44<<"\t  rho25    = "<<rho25<<"  "<<rhoMu<<endl;
  //float eleISOendcap = (lep->IsoMap("pfPhotonEt_R03") + lep->IsoMap("pfChargedHadron_R03") + lep->IsoMap("pfNeutralHadron_R03") - rho25*TMath::Pi()*0.09)/lep->Pt();
  //float eleISObarrel = (TMath::Max(0.0, lep->IsoMap("pfPhotonEt_R03")-1.0) + lep->IsoMap("pfChargedHadron_R03")  + lep->IsoMap("pfNeutralHadron_R03") - rho25*TMath::Pi()*0.09)/lep->Pt();

  //cout<<"\t\t  rel isolation (-rho) in barrel = "<<eleISObarrel<<"    and endcup = "<< eleISOendcap<<endl;

  cout<<"\t\t  pfPhotonEt_R03        = "<<lep->IsoMap("pfPhotonEt_R03")<<endl;
  cout<<"\t\t  pfChargedHadron_R03   = "<<lep->IsoMap("pfChargedHadron_R03")<<endl;
  cout<<"\t\t  pfNeutralHadron_R03   = "<<lep->IsoMap("pfNeutralHadron_R03")<<endl;

  cout<<"\t\t  pfChIso_R04         = "<<lep->IsoMap("pfChIso_R04")<<endl;
  cout<<"\t\t  pfNeuIso_R04        = "<<lep->IsoMap("pfNeuIso_R04")<<endl;
  //cout<<"\t\t  pfPhoIso_R04        = "<<lep->IsoMap("pfPhoIso_R04")<<endl;

  cout<<"\t\t  EffArea_R03        = "<<lep->IsoMap("EffArea_R03")<<endl;
  cout<<"\t\t  EffArea_R04        = "<<lep->IsoMap("EffArea_R04")<<endl;
  cout<<"\t\t eff area from the table  = "<<  EffAreaElectron(lep, "2012", 0, 0)<<endl;

  cout<<"\t\t HadOvEm              = "<<lep->HadOverEm()<<endl;
  cout<<"\t\t DphiSuperCluster      = "<<lep->DphiSuperCluster()<<endl;
  cout<<"\t\t DetaSuperCluster      = "<<lep->DetaSuperCluster()<<endl;
  cout<<"\t\t SigmaIetaIeta         = "<<lep->SigmaIetaIeta()<<endl;
  //cout<<"\t\t ConversionFlag        = "<<lep->ConversionFlag()<<endl;


  cout<<"\t\t  PtError()/Pt()        = "<<lep->PtError()/lep->Pt()<<endl;
  cout<<"\t\t  NormalizedChi2()      = "<<lep->NormalizedChi2()<<endl;

}



void higgsAnalyzer::DumpMuonInfo(TCMuon *lep, TVector3 *vert, Float_t rho44, Float_t rho25, Float_t rhoMu)
{
  cout<<"  A Muon with pt="<<lep->Pt()<<"  eta="<<lep->Eta()<<"  charge="<<lep->Charge()<<endl;
  //cout<<"\t\t  lep->NumberOfValidTrackerHits() = "<<lep->NumberOfValidTrackerHits()<<endl;
  cout<<"\t\t  fabs(lep->Dxy(vert))  = "<<fabs(lep->Dxy(vert))<<endl;
  cout<<"\t\t  fabs(lep->Dz(vert))   = "<<fabs(lep->Dz(vert))<<endl;
  cout<<"\t\t  rho44     = "<<rho44<<"\t  rho25    = "<<rho25<<"  "<<rhoMu<<endl;

  //cout<<"\t\t  rel isolation (-rho) in barrel = "<<eleISObarrel<<"    and endcup = "<< eleISOendcap<<endl;

  cout<<"\t\t  pfChargedPt_R04              = "<<lep->IsoMap("pfChargedPt_R04")<<endl;
  cout<<"\t\t  pfChargedHadronPt_R04        = "<<lep->IsoMap("pfChargedHadronPt_R04")<<endl;
  //cout<<"\t\t  pfPhoIso_R04                 = "<<lep->IsoMap("pfPhotonEt_R04")<<endl;
  cout<<"\t\t  pfNeutralHadronEt_R04        = "<<lep->IsoMap("pfNeutralHadronEt_R04")<<endl;

  //cout<<"\t\t tot iso - effrho              = "<<(lep->IsoMap("pfChargedHadronPt_R04") + lep->IsoMap("pfPhoIso_R04") + lep->IsoMap("pfNeutralHadronEt_R04")  - rho25*EffAreaMuon(lep, "2012", 0, 0))/lep->Pt()<<endl;

  cout<<"\t\t  EffArea_R03        = "<<lep->IdMap("EffArea_R03")<<endl;
  cout<<"\t\t  EffArea_R04        = "<<lep->IdMap("EffArea_R04")<<endl;

  cout<<"\t\t eff area from the table  = "<<  EffAreaMuon(lep, "2012", 0, 0)<<endl;

  cout<<"\t\t  PtError()/Pt()        = "<<lep->PtError()/lep->Pt()<<endl;
  cout<<"\t\t  NormalizedChi2()      = "<<lep->NormalizedChi2()<<endl;

}

float higgsAnalyzer::EffAreaElectron(TCElectron *lep, TString dataPeriod, Bool_t isData, Int_t a)
{
  Float_t area=0, area_neu = 0, area_phot=0, area_both=0;

  if (isData && dataPeriod.Contains("2011"))
    {
      if (fabs(lep->Eta())<1.0)         {area_neu = 0.044; area_phot = 0.140; area_both = 0.18;}
      else if  (fabs(lep->Eta())<1.479) {area_neu = 0.065; area_phot = 0.130; area_both = 0.20;}
      else if  (fabs(lep->Eta())<2.0)   {area_neu = 0.068; area_phot = 0.079; area_both = 0.15;}
      else if  (fabs(lep->Eta())<2.2)   {area_neu = 0.057; area_phot = 0.130; area_both = 0.19;}
      else if  (fabs(lep->Eta())<2.3)   {area_neu = 0.058; area_phot = 0.150; area_both = 0.21;}
      else if  (fabs(lep->Eta())<2.4)   {area_neu = 0.061; area_phot = 0.160; area_both = 0.22;}
      else                             {area_neu = 0.110; area_phot = 0.180; area_both = 0.29;}
    }

  if (dataPeriod.Contains("2012"))
    {
      if (fabs(lep->Eta())<1.0)         {area_both = 0.19;}
      else if  (fabs(lep->Eta())<1.479) {area_both = 0.25;}
      else if  (fabs(lep->Eta())<2.0)   {area_both = 0.12;}
      else if  (fabs(lep->Eta())<2.2)   {area_both = 0.21;}
      else if  (fabs(lep->Eta())<2.3)   {area_both = 0.27;}
      else if  (fabs(lep->Eta())<2.4)   {area_both = 0.44;}
      else                             {area_both = 0.52;}
    }

  if (a==0) area = area_both;
  else if  (a==1) area = area_neu;
  else if  (a==2) area = area_phot;

  return area;
}
float higgsAnalyzer::EffAreaMuon(TCMuon *lep, TString dataPeriod, Bool_t isData, Int_t a)
{
  Float_t area=0, area_neu = 0, area_phot=0, area_both=0;

  if (dataPeriod.Contains("2012"))
    {
      if (fabs(lep->Eta())<1.0)         {area_neu = 0.166; area_phot = 0.504; area_both = 0.674;}
      else if  (fabs(lep->Eta())<1.5)   {area_neu = 0.259; area_phot = 0.306; area_both = 0.565;}
      else if  (fabs(lep->Eta())<2.0)   {area_neu = 0.247; area_phot = 0.198; area_both = 0.442;}
      else if  (fabs(lep->Eta())<2.2)   {area_neu = 0.220; area_phot = 0.287; area_both = 0.515;}
      else if  (fabs(lep->Eta())<2.3)   {area_neu = 0.340; area_phot = 0.525; area_both = 0.821;}
      else if  (fabs(lep->Eta())<2.4)   {area_neu = 0.216; area_phot = 0.488; area_both = 0.660;}
      else                              {area_neu = 0.000; area_phot = 0.000; area_both = 0.00;}
    }

  if (dataPeriod.Contains("2011"))
    {
      if (fabs(lep->Eta())<1.0)         {area_both = 0.132;}
      else if  (fabs(lep->Eta())<1.5)   {area_both = 0.120;}
      else if  (fabs(lep->Eta())<2.0)   {area_both = 0.114;}
      else if  (fabs(lep->Eta())<2.2)   {area_both = 0.139;}
      else if  (fabs(lep->Eta())<2.3)   {area_both = 0.168;}
      else if  (fabs(lep->Eta())<2.4)   {area_both = 0.189;}
      else                              {area_both = 0.00;}
    }

  if (a==0) area = area_both;
  else if  (a==1) area = area_neu;
  else if  (a==2) area = area_phot;

  return area;
}



void higgsAnalyzer::MakeElectronPlots(TCElectron *ele, TVector3 *pv)
{
  Float_t eleIso = CalculateElectronIso(ele);

  hists->fill1DHist(ele->Dxy(pv), "ele_dxy", "dxy", 50, 0,0.3, 1, "Electrons");
  hists->fill1DHist(ele->Dz(pv), "ele_dz", "dz", 50, 0,0.3, 1, "Electrons");
  hists->fill1DHist(eleIso, "ele_iso", "Isolation", 50, 0,0.7, 1, "Electrons");
  hists->fill1DHist(ele->PtError()/ele->Pt(), "ele_ptErrorOverPt", "ptErrorOverPt", 50, 0,0.6, 1, "Electrons");
  hists->fill1DHist(ele->M(), "ele_mass", "Mass of the electron", 50, -4,4, 1, "Electrons");

}


float higgsAnalyzer::CalculateElectronIso(TCElectron *lep)
{
  float eleISO = 0;
  if(period=="2012")
    eleISO = (lep->IsoMap("pfChIso_R04") + TMath::Max(0.0, (Double_t)(lep->IsoMap("pfPhoIso_R04") + lep->IsoMap("pfNeuIso_R04") - rho25Factor*lep->IsoMap("EffArea_R04"))))/lep->Pt();
  else if(period=="2011")
    eleISO = (lep->IsoMap("SumPt_R04") + TMath::Max(0.0, (Double_t)(lep->IsoMap("EmIso_R04") + lep->IsoMap("HadIso_R04") - rhoFactor*TMath::Pi()*0.09)))/lep->Pt();
  else cout<<"Ele iso - Period is not defined!"<<period<<endl;

  return eleISO;
}


bool higgsAnalyzer::PassElectronIdAndIso(TCElectron *lep, elIdAndIsoCuts cuts, TVector3 *pv)
{
  bool pass = false;

  Float_t eleISO = CalculateElectronIso(lep);

  if(((fabs(lep->Eta()) < 1.442
       && lep->PtError()/lep->Pt()       < cuts.ptErrorOverPt[0]
       && lep->SigmaIetaIeta()           < cuts.sigmaIetaIeta[0]
       && fabs(lep->DphiSuperCluster())  < cuts.dPhiIn[0]
       && fabs(lep->DetaSuperCluster())  < cuts.dEtaIn[0]
       && lep->HadOverEm()               < cuts.HadOverEm[0]
       && (period=="2011" || lep->IsoMap("fabsEPDiff")      < cuts.fabsEPDiff[0])
       && fabs(lep->Dxy(pv))             < cuts.dxy[0]
       && fabs(lep->Dz(pv))              < cuts.dz[0]
       && eleISO                         < cuts.pfIso04[0]
       )||
      (fabs(lep->Eta()) >  1.556
       && lep->PtError()/lep->Pt()       < cuts.ptErrorOverPt[1]
       && lep->SigmaIetaIeta()           < cuts.sigmaIetaIeta[1]
       && fabs(lep->DphiSuperCluster())  < cuts.dPhiIn[1]
       && fabs(lep->DetaSuperCluster())  < cuts.dEtaIn[1]
       && lep->HadOverEm()               < cuts.HadOverEm[1]
       && (period=="2011" || lep->IsoMap("fabsEPDiff")      < cuts.fabsEPDiff[1])
       && fabs(lep->Dxy(pv))             < cuts.dxy[1]
       && fabs(lep->Dz(pv))              < cuts.dz[1]
       && eleISO                         < cuts.pfIso04[1]
       ))
     && lep->ConversionVeto()

     ) pass = true;

  //pass = true;

  //cout<<"ele pass? "<<pass<<endl;
  return pass;
}

void higgsAnalyzer::MakeMuonPlots(TCMuon *mu, TVector3 *pv)
{
  hists->fill1DHist(mu->NumberOfValidMuonHits(), "mu_NumberOfValidMuonHits", "NumberOfValidMuonHits", 60, 0,60, 1, "Muons");
  hists->fill1DHist(mu->NumberOfValidTrackerHits(), "mu_NumberOfValidTrackerHits", "NumberOfValidTrackerHits", 40, 0,40, 1, "Muons");
  hists->fill1DHist(mu->NumberOfValidPixelHits(), "mu_NumberOfValidPixelHits", "NumberOfValidPixelHits", 20, 0,20, 1, "Muons");
  hists->fill1DHist(mu->NormalizedChi2(), "mu_NormalizedChi2", "NormalizedChi2", 50, 0,12, 1, "Muons");
  hists->fill1DHist(mu->Dxy(pv), "mu_dxy", "dxy", 50, 0,0.3, 1, "Muons");
  hists->fill1DHist(mu->Dz(pv), "mu_dxy", "dz", 50, 0,0.3, 1, "Muons");

  Float_t muIso = CalculateMuonIso(mu);
  hists->fill1DHist(muIso, "mu_iso", "Isoalation", 50, 0,0.7, 1, "Muons");

  hists->fill1DHist(mu->PtError()/mu->Pt(), "mu_ptErrorOverPt", "ptErrorOverPt", 50, 0,0.6, 1, "Muons");
}


float higgsAnalyzer::CalculateMuonIso(TCMuon *lep)
{
  Double_t area = EffAreaMuon(lep, "2012",0,0);
  float muISO = 0;
  if (period=="2012")
    muISO = (lep->IsoMap("pfChargedHadronPt_R04") + TMath::Max(0.0, lep->IsoMap("pfPhotonEt_R04") + lep->IsoMap("pfNeutralHadronEt_R04")  - TMath::Max(0.0, (Double_t)rho25Factor*area)))/lep->Pt();
  else if(period=="2011")
    muISO = (lep->IsoMap("pfChargedHadronPt_R04") + TMath::Max(0.0, lep->IsoMap("pfPhotonEt_R04") + lep->IsoMap("pfNeutralEt_R04")  - TMath::Max(0.0, (Double_t)rhoFactor*TMath::Pi()*0.09)))/lep->Pt();
  else cout<<"Muon iso - Period is not defined!"<<period<<endl;

  return muISO;
}


bool higgsAnalyzer::PassMuonIdAndIso(TCMuon *lep, muIdAndIsoCuts cuts, TVector3 *pv)
{
  bool pass = false;

  Float_t muISO =  CalculateMuonIso(lep);

  if( lep->PtError()/lep->Pt() < cuts.ptErrorOverPt
      && lep->NumberOfValidMuonHits()    > cuts.NumberOfValidMuonHits
      && lep->NumberOfValidTrackerHits() > cuts.NumberOfValidTrackerHits
      && lep->NumberOfValidPixelHits()   > cuts.NumberOfValidPixelHits
      && lep->NumberOfMatches() > cuts.NumberOfMatches
      && lep->NormalizedChi2()  < cuts.NormalizedChi2
      && fabs(lep->Dxy(pv))     < cuts.dxy
      && fabs(lep->Dz(pv))      < cuts.dz
      && muISO                  < cuts.pfIso04
      )
    pass = true;
  //cout<<"muon pass? "<<pass<<endl;

  //pass = false;
  //if (muISO < cuts.pfIso04)
  //  pass = true;

  return pass;
}

