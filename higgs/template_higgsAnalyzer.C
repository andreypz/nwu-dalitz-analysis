// $Id: template_higgsAnalyzer.C,v 1.31 2012/06/22 22:11:44 andrey Exp $

#define higgsAnalyzer_cxx

#include "higgsAnalyzer.h"
#include <string>
#include "../plugins/MetDefinitions.h"
using namespace std;

/////////////////////////////
//Specify parameters here. //
/////////////////////////////

string  selection      = "SELECTION";
string  period         = "PERIOD";
int     JC_LVL         = 4;
int     trigger[]      = {TRIGGER};
string  suffix         = "SUFFIX";

vector<int> triggers (trigger, trigger + sizeof(trigger)/sizeof(int));

UInt_t verboseLvl  = 0;
Bool_t doZlibrary  = 0, isFromData=1;
Bool_t makeKinTree = 0;

/////////////////
//Analysis cuts//
/////////////////

Float_t jetPtCut[]     = {30., 15.};
Float_t bJetPtCut      = 30.;
Float_t muPtCut[]      = {20., 10.};
Float_t elePtCut[]     = {20., 10.};
float   photonPtCut[]  = {25., 1e9};
Float_t zMassCut[]     = {76.2, 106.2};
Float_t qtCut          = 55.;
Int_t   nJetsCut[]     = {0, 99};
Float_t cut_vz = 24, cut_vd0 = 2, cut_vndof = 4;  //PV filter cuts

// Cuts for mass points in the PAS. 250,300,350 etc
Float_t dPhiMinCut[] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
//Float_t dPhiMinCut[] = {0.47, 0.33, 0.21, 0.12, 0.06, 0.01, 0.00, 0.00};
Float_t metMinCut[]  = {70,     79,   95,  115,  134,  150, 159, 160};
Float_t mtMinCut[]   = {222,   264,  298,  327,  354,  382, 413, 452};
Float_t mtMaxCut[]   = {272,   331,  393,  460,  531,  605, 684, 767};

//For Anton's tree:
TLorentzVector  kt_lep1, kt_lep2, kt_met, kt_corrMet, kt_allJets, kt_V;
UInt_t kt_nPV, kt_nSimPV, kt_nJets;
Bool_t kt_isMC, kt_hasBJet;
Float_t kt_dPhiJetMet, kt_weight;
ULong64_t kt_evNumber;

// Do something about these: should just have one sort condition function
bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());} 
bool MuonSortCondition(const TCMuon& m1, const TCMuon& m2) {return (m1.Pt() > m2.Pt());}
bool ElectronSortCondition(const TCElectron& e1, const TCElectron& e2) {return (e1.Pt() > e2.Pt());}
bool PhotonSortCondition(const TCPhoton& g1, const TCPhoton& g2) {return (g1.Pt() > g2.Pt());}

TRandom3 *myRandom;

void higgsAnalyzer::Begin(TTree * /*tree*/) 
{
    TString option = GetOption();
    TH1::SetDefaultSumw2(kTRUE);
    
    cout<<"\n***      Begin the Analyzer      ****"<<endl;
    myRandom = new TRandom3();


    // Initialize utilities and selectors here //
    zLib            = new ZedEventsLibrary(selection, isFromData);
    weighter        = new WeightUtils(suffix, period, selection, isRealData);
    triggerSelector = new TriggerSelector(selection, period, triggers);
    
    histoFile = new TFile("a_higgsHistograms.root", "RECREATE");
    histoFile->mkdir("Andrey", "Andrey");
    histoFile->mkdir("Anton", "Anton");
    histoFile->cd("Andrey");

    //In the Title of this histogram I encode the sample!  
    evt_byCut = new TH1F("evt_byCut", "SUFFIX", 20, 0,20);
    evt_libQt = new TH2F("evt_libQt", "Events in a library", 4, 0,4,  10, 0,10);

    if (doZlibrary)  zLib->CreateLibrary();
    
    for (Int_t n=0; n<nC; n++)
      {
      //cout<<"zeroing the array: "<<n<<endl;
      nEvents[n]=0;
      nEventsWeighted[n]=0;
    }

    for(Int_t n=0; n<nC; n++)
      met0_et[n]       = new TH1F(Form("met0_et_%i",n), "met0_et", 40, 0,400);

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


	//met0_et_ovQt[n]  = new TH2F(Form("met0_et_ovQt_%i",n), "pfMET vs Met/qt", 40, 0,400, 40, 0,4);
	met1_et_ovQt[n]  = new TH2F(Form("met1_et_ovQt_%i",n), "MET1 vs Met/qt", 40, 0,400, 40, 0,4);
	//met2_et_ovQt[n]  = new TH2F(Form("met2_et_ovQt_%i",n), "pfMet noise vs Met/qt", 40, 0,400, 40, 0,4);
	//met3_et_ovQt[n]  = new TH2F(Form("met3_et_ovQt_%i",n), "projMET vs Met/qt", 40, 0,400, 40, 0,4);
	//met4_et_ovQt[n]  = new TH2F(Form("met4_et_ovQt_%i",n), "puCorrMET vs Met/qt", 40, 0,400, 40, 0,4);

	met1_phi[n]     = new TH1F(Form("met1_phi_%i",n), "met1_phi", 40, -TMath::Pi(), TMath::Pi());
	met1_over_qt[n] = new TH1F(Form("met1_over_qt_%i",n), "met1_over_qt", 40, 0,4);

	met1_lg[n]         = new TH1F(Form("met1_lg_%i",n), "met1_lg", 50,-200,100);
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
	mt2_met2[n]   = new TH2F(Form("mt2_met2_%i",n), "MT pf vs pfMet noise", 50, 100, 600, 40, 0,400);

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
	jet_dRlep1[n] = new TH1F(Form("jet_dRlep1_%i",n), "dR(jet, lepton1)", 50, 0,5);
	jet_dRlep2[n] = new TH1F(Form("jet_dRlep2_%i",n), "dR(jet, lepton2)", 50, 0,5);
	jet_pt[n]     = new TH1F(Form("jet_pt_%i",n),  "Pt of leading jet", 50,0,400);
	jet_eta[n]    = new TH1F(Form("jet_eta_%i",n), "Eta of leading jet", 50,-5,5);
	jet_phi[n]    = new TH1F(Form("jet_phi_%i",n), "Phi of leading jet", 50, -TMath::Pi(), TMath::Pi());

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

	evt_weight[n]      = new  TH1F(Form("evt_weight_%i",n), "Weights", 50,0,10);

	run_events[n]     =  new  TH1F(Form("run_events_%i",n), "Events per run", 18000, 160000., 178000.);


	if (suffix.compare(0, 5, "ggHZZ") == 0 || suffix.compare(0, 5, "ggHWW") == 0) {
	  higgs_pt[n]     = new  TH1F(Form("higgs_pt_%i",n), "Higgs pt", 50,0,500);
	  higgs_w_pt[n]   = new  TH1F(Form("higgs_w_pt_%i",n), "Higgs pt, wighted", 50,0,500);

	  higgs_mass[n]   = new  TH1F(Form("higgs_mass_%i",n), "Higgs mass", 150, 0,1500);
	  higgs_w_mass[n] = new  TH1F(Form("higgs_w_mass_%i",n), "Higgs mass, weighted", 150, 0,1500);
	  //higgs_mass[n] = new  TH1F(Form("higgs_mass_%i",n), "Higgs mass", 70,100,800);
	}


       }
    if (verboseLvl>0){
      ffout.open("./events_printout_SUFFIX_final.txt",ofstream::out);
      ncout.open("./counts_for_tex_SUFFIX.txt",ofstream::out);
      
      ffout.precision(4); ffout.setf(ios::fixed, ios::floatfield);
      for(Int_t i=6; i<nC; i++)
	{
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

    histoFile->cd("Anton");    

    if (makeKinTree){
      TString kinTreeName = "kinTree";
      _kinTree = new TTree(kinTreeName, "Tree with kinematic info");
      _kinTree->Branch("V", "TLorentzVector",  &kt_V);
      _kinTree->Branch("lep1", "TLorentzVector",  &kt_lep1);
      _kinTree->Branch("lep2", "TLorentzVector",  &kt_lep2);
      _kinTree->Branch("met", "TLorentzVector",  &kt_met);
      _kinTree->Branch("corrMet", "TLorentzVector",  &kt_corrMet);
      _kinTree->Branch("allJets", "TLorentzVector",  &kt_allJets);
      _kinTree->Branch("nPV", &kt_nPV,  "nPV/i");
      _kinTree->Branch("nSimPV", &kt_nSimPV,  "nSimPV/i");
      _kinTree->Branch("nJets", &kt_nJets,  "nJets/i");
      _kinTree->Branch("isMc", &kt_isMC,  "isMC/O");
      _kinTree->Branch("hasBJet", &kt_hasBJet,  "hasBJet/O");
      _kinTree->Branch("dPhiJetMet", &kt_dPhiJetMet, "dPhiJetMet/F");
      _kinTree->Branch("evNumber", &kt_evNumber, "evNumber/l");
      _kinTree->Branch("weight", &kt_weight, "weight/F");
    }
}
bool higgsAnalyzer::Process(Long64_t entry)
{  
    GetEntry(entry);
    CountEvents(0);
    ++nEventsWeighted[0];
    MET = 0;
    FillHistosNoise(0, 1);

    //cout<<nEvents[0]<<endl;
    if (nEvents[0] == 1) weighter->SetDataBit(isRealData);
    if(nEvents[0]>200) return kTRUE;
    
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
    
    bool triggerPass   = triggerSelector->SelectTriggers(triggerStatus, hltPrescale);
    if (!triggerPass) return kTRUE;
    // Double electron workaround.  Gets rid of hopelessly prescaled events of July 20-26, 2011
    if (selection == "electron" && (runNumber > 171200 && runNumber < 171600)) return kTRUE;

    int  eventPrescale = triggerSelector->GetEventPrescale();

    CountEvents(1);
    ++nEventsWeighted[1];
    FillHistosNoise(1, 1);

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
    

    ///////////////
    // electrons //
    ///////////////
    vector<TLorentzVector> looseLeptons;
    vector<TCElectron> electrons;
    //int eleCount = 0;

    for (int i = 0; i <  recoElectrons->GetSize(); ++i) {
        TCElectron* thisElec = (TCElectron*) recoElectrons->At(i);    

        if (fabs(thisElec->Eta()) > 2.5) continue;

        float eleISOendcap = (thisElec->TrkIso() + thisElec->EmIso() + thisElec->HadIso() - rhoFactor*TMath::Pi()*0.09)/thisElec->Pt(); 
        float eleISObarrel = (thisElec->TrkIso() + TMath::Max(0.0, thisElec->EmIso()-1.0)+ thisElec->HadIso() - rhoFactor*TMath::Pi()*0.09)/thisElec->Pt(); 
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
                    thisElec->Pt() > 10  
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
		   ) looseLeptons.push_back(thisElec->P4());
	
    } 

    sort(electrons.begin(), electrons.end(), ElectronSortCondition);

    ///////////
    // muons //
    ///////////

    vector<TCMuon> muons;
    Int_t softMuons = 0;
    for (int i = 0; i < recoMuons->GetSize(); ++ i) {
        TCMuon* thisMuon = (TCMuon*) recoMuons->At(i);    
	
        if (!(fabs(thisMuon->Eta()) < 2.4	
	      && thisMuon->IsGLB() 
	      && thisMuon->IsTRK()
	      && thisMuon->PtError()/thisMuon->Pt() < 0.1
	      )) continue;
	
        if (
	    thisMuon->Pt() > 3
	    && thisMuon->NumberOfValidTrackerHits() > 10
	    && fabs(thisMuon->Dxy(pvPosition)) < 0.3
	    && fabs(thisMuon->Dz(pvPosition))  < 0.30 
	    && (thisMuon->TrkIso() + thisMuon->HadIso() + thisMuon->EmIso() - rhoFactor*TMath::Pi()*0.09)/thisMuon->Pt() < 0.15
	    ) {
	  softMuons++;

	  if(thisMuon->Pt()>10
	     && thisMuon->NumberOfValidMuonHits()    > 0
	     && thisMuon->NumberOfValidTrackerHits() > 10
	     && thisMuon->NumberOfValidPixelHits()   > 0
	     && thisMuon->NumberOfMatches() > 1
	     && thisMuon->NormalizedChi2()  < 10
	     && fabs(thisMuon->Dxy(pvPosition)) < 0.2
	     && fabs(thisMuon->Dz(pvPosition))  < 0.5 
	     && (thisMuon->TrkIso() + thisMuon->HadIso() + thisMuon->EmIso() - rhoFactor*TMath::Pi()*0.09)/thisMuon->Pt() < 0.15
	     )	  muons.push_back(*thisMuon);
     
	} else if (
		   thisMuon->Pt() > 10 
		   //&& thisMuon->Pt() <= muPtCut[0] 
		   && thisMuon->NumberOfValidMuonHits()    > 0
		   && thisMuon->NumberOfValidTrackerHits() > 2 
		   && thisMuon->NumberOfValidPixelHits()   > 0
		   && thisMuon->NumberOfMatches() > 0
		   && thisMuon->NormalizedChi2()  < 100
		   && fabs(thisMuon->Dxy(pvPosition)) < 0.02
		   && fabs(thisMuon->Dz(pvPosition))  < 0.1
		   //&& (thisMuon->TrkIso() + thisMuon->HadIso() + thisMuon->EmIso() - rhoFactor*TMath::Pi()*0.09)/thisMuon->Pt() < 0.15
		   ) {
	  looseLeptons.push_back(thisMuon->P4());
        }
    } 

    sort(muons.begin(), muons.end(), MuonSortCondition);
    sort(looseLeptons.begin(), looseLeptons.end(), P4SortCondition);

    /////////////
    // photons //
    /////////////
    nGamma = 0;
    vector<TCPhoton> photons; 
    if (selection == "eGamma" || selection == "muGamma" || selection == "gamma") {
        for (Int_t i = 0; i < recoPhotons->GetSize(); ++i)
        {
	  TCPhoton* thisPhoton = (TCPhoton*) recoPhotons->At(i);
	  //cout<<i+1<<" dbg  pt: "<<thisPhoton->Pt()<<"  eta: "<< thisPhoton->Eta()<<endl;
	  if (thisPhoton->Pt() < photonPtCut[0] || thisPhoton->Pt() > photonPtCut[1]) continue;	  
	  if ( 
	      thisPhoton->Pt() >10
	      //&& thisPhoton->EmIso()         < (4.2 + 0.006 *thisPhoton->Pt())
	      //&& thisPhoton->HadIso()        < (2.2 + 0.0025*thisPhoton->Pt())
	      //&& thisPhoton->TrkIso()        < (2.0 + 0.001 *thisPhoton->Pt())
	      && thisPhoton->HadOverEm()     < 0.05 
	      && thisPhoton->SigmaIEtaIEta() > 0.001
	      && thisPhoton->SigmaIEtaIEta() < 0.013
	      && thisPhoton->TrackVeto() == 0
	      && fabs(thisPhoton->EtaSupercluster()) < 1.442
	       )  
	    {
	      photons.push_back(*thisPhoton);

	      //for gamma multiplicity plot:
	      if(thisPhoton->Pt() > 25) nGamma++;	
	    }

	}
	sort(photons.begin(), photons.end(), PhotonSortCondition);
	  
        if (photons.size() > 0 && eventPrescale > 1) {
	  if (!triggerSelector->PhotonTriggerBins(photons[0].Pt(), true)) return kTRUE;
        }
    }


    //////////
    // jets //
    //////////

    vector<TLorentzVector> jetP4, bJetP4, softJetP4, bSSVJetP4, testJetP4;
    TLorentzVector sumJetP4(0,0,0,0);
    int jetCount = 0;
    float fwdJetSumPt = 0.;
    int fwdJetCount = 0;

    nJets=0; nJetsB=0, nJetsBssv=0;
    for (int i = 0; i < recoJets->GetSize(); ++i) {
      TCJet* thisJet = (TCJet*) recoJets->At(i);     

      // Prevent lepton overlap //
      bool leptonOverlap = false;
      for (int j = 0; j < (int)muons.size(); ++j)     if (thisJet->P4().DeltaR(muons[j].P4()) < 0.4) leptonOverlap = true;
      for (int j = 0; j < (int)electrons.size(); ++j) if (thisJet->P4().DeltaR(electrons[j].P4()) < 0.4) leptonOverlap = true;
      
      if (selection == "eGamma" || selection == "muGamma") {
	for (int j = 0; j < (int)photons.size(); ++j) if (thisJet->P4().DeltaR(photons[j].P4()) < 0.4) leptonOverlap = true;
      }
      if (leptonOverlap) continue;
      
      if (fabs(thisJet->P4().Eta()) < 2.4) {
	if (
	    thisJet->NumConstit()   > 1
	    && thisJet->ChHadFrac() > 0
	    && thisJet->ChEmFrac()  < 1
	    && (thisJet->NeuHadFrac() + thisJet->NeuEmFrac()) < 0.9
	    //&& thisJet->P4(JC_LVL).Pt() > jetPtCut[0]
	    ) {
	  if (
	      kTRUE  //Removing jet-vertex assosiation
	      //thisJet->VtxNTracks() > 0
	      //&& thisJet->VtxSumPtFrac() > 0. 
	      //&& thisJet->VtxSumPtIndex() == (UInt_t)(PVind[0]+1) //index in pv collection (passed vtx selection)
	      ) {

	    if (thisJet->P4(JC_LVL).Pt() > jetPtCut[0]) {
	      if (thisJet->BDiscrTCHE() > 2. && thisJet->P4(JC_LVL).Pt() > bJetPtCut) 
		bJetP4.push_back(thisJet->P4(JC_LVL));
	      if (( thisJet->BDiscrTCHE() > 2 || thisJet->BDiscrSSVHE()> 1.74 )&& thisJet->P4(JC_LVL).Pt() > bJetPtCut) 
		bSSVJetP4.push_back(thisJet->P4(JC_LVL));
	      
	      jetP4.push_back(thisJet->P4(JC_LVL));
	      sumJetP4 += thisJet->P4();
	      ++jetCount;
	    }
	    else if (thisJet->P4(JC_LVL).Pt() > jetPtCut[1]) softJetP4.push_back(thisJet->P4(JC_LVL));
	  }
	}
      } else if (fabs(thisJet->P4().Eta()) < 4.8) {
	if (thisJet->P4(JC_LVL).Pt() > jetPtCut[0]) {
     	  jetP4.push_back(thisJet->P4(JC_LVL)); 
	  sumJetP4 += thisJet->P4();
	  fwdJetSumPt += thisJet->Pt(JC_LVL);
	  ++fwdJetCount;
	  ++jetCount;
	} else if (thisJet->P4(JC_LVL).Pt() > jetPtCut[1]) softJetP4.push_back(thisJet->P4(JC_LVL));
      }
    }

    sort(jetP4.begin(), jetP4.end(), P4SortCondition);
    sort(bJetP4.begin(), bJetP4.end(), P4SortCondition);
    sort(bSSVJetP4.begin(), bSSVJetP4.end(), P4SortCondition);

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
    metP4.SetPtEtaPhiE(met->Met(), 0, met->Phi(), met->Met());
    met1P4.SetPtEtaPhiE(met->CorrectedMet(), 0, met->CorrectedPhi(), met->CorrectedMet());
    ////////////////////////
    // Analysis selection //
    ////////////////////////

    TLorentzVector ZP4;
    MET       = met->Met();
    pfMET     = met->Met();
    pfMET1    = met->CorrectedMet();
    puCorrMET = PUCorrectedMET(pfMET, nVtx, "pfMet", isRealData);

    ///////////////////////
    // Cosmics rejection //
    ///////////////////////

    if (selection=="muon" && muons.size() > 1) if (CosmicMuonFilter(muons[0], muons[1])) return kTRUE;
    //Data quality
    if (isRealData && (isNoiseHcal || isScraping || isCSCTightHalo)) return kTRUE;

    TLorentzVector Lepton1(0.,0.,0.,0.), Lepton2(0.,0.,0.,0.);
    if (selection == "electron") {
        /////////////////////
        // 2 good elctrons //
        /////////////////////

        if (electrons.size() < 2) return kTRUE;
	//no opposite charge requirement
	//if (electrons[0].Charge() == electrons[1].Charge()) return kTRUE;

        ZP4           = electrons[0].P4() + electrons[1].P4();
        reducedMet1P4 = GetReducedMET(sumJetP4,  electrons[0].P4(), electrons[1].P4(), metP4, 1);
        reducedMet2P4 = GetReducedMET(sumJetP4,  electrons[0].P4(), electrons[1].P4(), metP4, 2);

	Lepton1 = electrons[0].P4();
	Lepton2 = electrons[1].P4();

    } else if (selection == "muon") {

        //////////////////////////////////////////
        // 2 (oppositely charged) muons  w/ pt>20 //
        //////////////////////////////////////////

        if (muons.size() < 2) return kTRUE;
	//no opposite charge requirement
        //if (muons[0].Charge() == muons[1].Charge()) return kTRUE;

        ZP4           = muons[0].P4() + muons[1].P4();
        reducedMet1P4 = GetReducedMET(sumJetP4,  muons[0].P4(), muons[1].P4(), metP4, 1);
        reducedMet2P4 = GetReducedMET(sumJetP4,  muons[0].P4(), muons[1].P4(), metP4, 2);

	Lepton1 = muons[0].P4();
	Lepton2 = muons[1].P4();

    } else if (selection == "eGamma" || selection == "muGamma" || selection == "gamma") {

        /////////////////////////////////////
        // Photons in place of Z candidate //
        /////////////////////////////////////

        if (photons.size() < 1) return kTRUE;
	float photonMass = weighter->GetPhotonMass(); 
	ZP4.SetPtEtaPhiM(photons[0].Pt(), photons[0].Eta(), photons[0].Phi(), photonMass);
        reducedMet1P4 = GetReducedMET(sumJetP4,  photons[0].P4(), TLorentzVector(0, 0, 0, 0), metP4, 1);
	reducedMet2P4 = GetReducedMET(sumJetP4,  photons[0].P4(), TLorentzVector(0, 0, 0, 0), metP4, 2);

	//needed for weights
	Lepton1 = ZP4;  
	Lepton2.SetXYZT(0,0,0,0);

    } else {
        return kTRUE;
    }
    //Event Weight
    float evtWeight = 1.0;\
    // Doesn't work for now. Need to check with Brian/Nate
    evtWeight   = weighter->GetTotalWeight(nPUVertices, jetP4.size(), Lepton1, Lepton2);

    if ((selection == "gamma" || selection == "muGamma" || selection == "eGamma") && isRealData) {
      //cout<<"Weight, before/after prescale   "<<evtWeight;
      evtWeight  *= eventPrescale;
      //cout<<"\t\t"<<evtWeight<<endl;
      //cout<<"nPVs = "<<primaryVtx->GetSize()<<"\t gamma pt= "<<ZP4.Pt()<<endl;
    }

    //DeltaPhi
    float deltaPhiJetMET = 2*TMath::Pi();
    if (jetP4.size() > 0)            deltaPhiJetMET = DeltaPhiJetMET(metP4, jetP4);
    else if (softJetP4.size() > 0)   deltaPhiJetMET = DeltaPhiJetMET(metP4, softJetP4);



    /////////////////// 
    // Gen particles //
    ///////////////////
    
    if (!isRealData && (suffix.compare(2, 3, "HZZ") == 0 || suffix.compare(2, 3, "HWW") == 0 || suffix.compare(3, 3, "HZZ") == 0)) {
      vector<TLorentzVector> genLeptons;
      vector<TLorentzVector> genNeutrinos;
    
      // NLO k-faktors weights
      for (int i = 0; i < genParticles->GetSize(); ++i) {
	TCGenParticle* thisParticle = (TCGenParticle*) genParticles->At(i);    
	
	if (
	    (fabs(thisParticle->GetPDGId()) == 11 
	     || fabs(thisParticle->GetPDGId()) == 13 
	     || fabs(thisParticle->GetPDGId()) == 15) 
	    && thisParticle->Mother() ==23
	    ) genLeptons.push_back(thisParticle->P4());
	
	if (
	    (fabs(thisParticle->GetPDGId()) == 12 
	     || fabs(thisParticle->GetPDGId()) == 14 
	     || fabs(thisParticle->GetPDGId()) == 16) 
	    && thisParticle->Mother() ==23
	    ) genNeutrinos.push_back(thisParticle->P4());


	//if (thisParticle->GetPDGId() == 25)
	// {
	//
	// }
      }

      if (genLeptons.size() > 1 && genNeutrinos.size() > 1) {
	hig_Pt   = (genLeptons[0]+genLeptons[1]+genNeutrinos[0]+genNeutrinos[1]).Pt();
	hig_M    = (genLeptons[0]+genLeptons[1]+genNeutrinos[0]+genNeutrinos[1]).M();

	Int_t i_mass = 0;
	if (suffix.compare(0, 5, "ggHZZ") == 0 || suffix.compare(0, 5, "ggHWW") == 0) {
	  i_mass = atoi(suffix.substr(5,3).c_str());
	  evtWeight *= weighter->GluGluHiggsWeight(hig_Pt, i_mass);
	} 
	else if (suffix.compare(0, 6, "VBFHZZ") == 0){
	  i_mass = atoi(suffix.substr(6,3).c_str()); //for vbf name
	}
	
	//cout<<"nPU vertic = "<<nPUVertices<<"  evtWeight before reweighting  "<<evtWeight<<endl;
	
	// Higgs mass lineshapes weights
	evtWeight *= weighter->HiggsMassLineShapeWeight(hig_M, i_mass);
	//cout<<"imass: "<<i_mass<<"   genMass from gen Z's: "<<hig_M<<"   evtWeight= "<<evtWeight<<endl;	  

      }
      
      //cout<<"Higgs Pt and Mass form ZZ: "<<hig_Pt<<"  "<<hig_M<<endl;
      
    }
   
    //cout<<"dbg"<<endl;    



    /////////////////////
    // 3rd lepton veto //
    /////////////////////

   if(selection  == "muon" ||  selection  == "electron")
     if (Lepton1.Pt() < 20 || Lepton2.Pt() < 20) return kTRUE;

    if ((electrons.size() > 0 || muons.size() > 2) && selection  == "muon") return kTRUE;
    if ((muons.size() > 0 || electrons.size() > 2) && selection  == "electron") return kTRUE;
    if ((muons.size() > 0 || electrons.size() > 0) && (selection == "eGamma" || selection == "muGamma" || selection == "gamma")) return kTRUE;

    //Loose third lepton Veto:
    //if( (selection == "muon" || selection == "electron") && looseLeptons.size()>0) return kTRUE;

    CountEvents(2);
    nEventsWeighted[2] += evtWeight;
    FillHistosNoise(2, evtWeight);

    ////////////
    // Z mass //
    ////////////

    if (ZP4.M() < zMassCut[0] || ZP4.M() > zMassCut[1]) return kTRUE;  
    CountEvents(3);
    nEventsWeighted[3] += evtWeight;
    FillHistosNoise(3, evtWeight);


    //Soft third muon veto:
    if (softMuons > 2) return kTRUE;

    CountEvents(4);
    nEventsWeighted[4] += evtWeight;

    
    /////////////////////////////////////////////////
    // Variables for FillHistos() function (Andrey)//
    /////////////////////////////////////////////////
    qT      = ZP4.Pt();
    Mll     = ZP4.M(); Mll_EE=0; Mll_EB=0; Mll_EX=0;
    diEta   = ZP4.Eta();
    diPhi   = ZP4.Phi();
    MT      = CalculateTransMass(metP4, ZP4);
    MT1     = CalculateTransMass(metP4, ZP4);
    METqt   = metP4.Pt()/qT;
    MET_phi = metP4.Phi();
    projMET = projectedMET(MET, MET_phi, Lepton1, Lepton2);
    ZprojMET= ZprojectedMET(MET, MET_phi, Lepton1, Lepton2);
    redMET1 = reducedMet1P4.Pt();
    redMET2 = reducedMet2P4.Pt();
    compMET = (sumJetP4 + ZP4).Pt();

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

    pfMET_lg = metP4.Pt()*cos(metP4.DeltaPhi(ZP4)); 
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
    }
    if (bJetP4.size()!=0)
      ptLeadBJet  = bJetP4[0].Pt();

    ///-------------------/////


    FillHistosNoise(4, evtWeight);
    FillHistos(4, evtWeight);
    
    ////////////////
    // b-jet veto //
    ////////////////
    
    Bool_t passBveto = kFALSE;
    if (nJetsB == 0)
      passBveto = kTRUE;

    if (makeKinTree){
      kt_V          = ZP4;     //Z or Gamma
      kt_lep1       = Lepton1;
      kt_lep2       = Lepton2;
      kt_met        = metP4;
      kt_corrMet    = met1P4;
      kt_allJets    = sumJetP4;
      kt_nPV        = nVtx;
      kt_nSimPV     = nPUVertices;
      kt_nJets      = nJets;
      kt_isMC       = (!isRealData);
      kt_hasBJet    = (!passBveto);
      kt_dPhiJetMet = deltaPhiJetMET;
      kt_evNumber   = eventNumber;
      kt_weight     = evtWeight;
      _kinTree -> Fill();
    }

    if (ZP4.Pt() < qtCut) return kTRUE;  
    CountEvents(5);
    nEventsWeighted[5] += evtWeight;
    FillHistosNoise(5, evtWeight);
    

    if(passBveto){
      CountEvents(6);
      nEventsWeighted[6] += evtWeight;
      FillHistos(6, evtWeight);
      FillHistosNoise(6, evtWeight);
    }    
        
    Bool_t passBasic = kTRUE;

    /////////////////////////////////////////
    // DeltaPhi(MET, jet); MET and MT cuts //
    /////////////////////////////////////////

    if (passBveto && deltaPhiJetMET > dPhiMinCut[0] &&  MET>metMinCut[0]  && (MT > mtMinCut[0] && MT < mtMaxCut[0]))
      {
	//For Higgs  200 . The cuts need to be defined
	CountEvents(7);
	nEventsWeighted[7] += evtWeight;
	//FillHistos(7, evtWeight);
	FillHistosNoise(7, evtWeight);
      }

    if (passBveto && deltaPhiJetMET > dPhiMinCut[0] &&  MET>metMinCut[0]  && (MT > mtMinCut[0] && MT < mtMaxCut[0]))
      {
	//For Higgs  250
	CountEvents(8);
	nEventsWeighted[8] += evtWeight;
	FillHistos(8, evtWeight);
	FillHistosNoise(8, evtWeight);
	PrintOut(8);
      }

    if (passBveto && deltaPhiJetMET > dPhiMinCut[1] &&  MET>metMinCut[1]  && (MT > mtMinCut[1] && MT < mtMaxCut[1]))
      {
	//For Higgs  300
	CountEvents(9);
	nEventsWeighted[9] += evtWeight;
       	//FillHistos(9, evtWeight);
	FillHistosNoise(9, evtWeight);
	PrintOut(9);
      }

    if (passBveto && deltaPhiJetMET > dPhiMinCut[2] &&  MET>metMinCut[2]  && (MT > mtMinCut[2] && MT < mtMaxCut[2]))
      {
	//For Higgs  350
	CountEvents(10);
	nEventsWeighted[10] += evtWeight;
       	//FillHistos(10, evtWeight);
	FillHistosNoise(10, evtWeight);
	PrintOut(10);
      }
    if (passBveto && deltaPhiJetMET > dPhiMinCut[3] &&  MET>metMinCut[3]  && (MT > mtMinCut[3] && MT < mtMaxCut[3]))
      {
	//For Higgs  400
	CountEvents(11);
	nEventsWeighted[11] += evtWeight;
       	FillHistos(11, evtWeight);
	FillHistosNoise(11, evtWeight);
	PrintOut(11);

      }
    if (passBveto && deltaPhiJetMET > dPhiMinCut[4] &&  MET>metMinCut[4]  && (MT > mtMinCut[4] && MT < mtMaxCut[4]))
      {
	//For Higgs  450
	CountEvents(12);
	nEventsWeighted[12] += evtWeight;
       	//FillHistos(12, evtWeight);
    	FillHistosNoise(12, evtWeight);
      }

    if (passBveto && deltaPhiJetMET > dPhiMinCut[5] &&  MET>metMinCut[5]  && (MT > mtMinCut[5] && MT < mtMaxCut[5]))
      {
	//For Higgs  500
	CountEvents(13);
	nEventsWeighted[13] += evtWeight;
       	//FillHistos(13, evtWeight);
    	FillHistosNoise(13, evtWeight);
      }

    if (passBveto && deltaPhiJetMET > dPhiMinCut[6] &&  MET>metMinCut[6]  && (MT > mtMinCut[6] && MT < mtMaxCut[6]))
      {
	//For Higgs  550
	CountEvents(14);
	nEventsWeighted[14] += evtWeight;
       	//FillHistos(14, evtWeight);
    	FillHistosNoise(14, evtWeight);
      }

    if (passBveto && deltaPhiJetMET > dPhiMinCut[7] &&  MET>metMinCut[7]  && (MT > mtMinCut[7] && MT < mtMaxCut[7]))
      {
	//For Higgs  600
	CountEvents(15);
	nEventsWeighted[15] += evtWeight;
       	//FillHistos(15, evtWeight);
    	FillHistosNoise(15, evtWeight);
      }
    
    if (passBasic && MET> 70)
      {
	//NO b-veto, only basics cuts. For N-b jets distribution.
	CountEvents(16);
	nEventsWeighted[16] += evtWeight;
       	FillHistos(16, evtWeight);
	FillHistosNoise(16, evtWeight);
      }

    if (passBasic && passBveto && MET> 40)
      {
	CountEvents(17);
	nEventsWeighted[17] += evtWeight;
	FillHistos(17, evtWeight);
	FillHistosNoise(17, evtWeight);
      }
   

    //////////////////////
    // Jet Multiplicity //
    //////////////////////
    /*
    if ((int)jetP4.size() < nJetsCut[0] || (int)jetP4.size() > nJetsCut[1]) return kTRUE;
    CountEvents(13);
    nEventsWeighted[13] += evtWeight;
    */

    //cout<<"dbg END"<<endl;    

    return kTRUE;
}

void higgsAnalyzer::Terminate()
{
  for (int i = 6; i < nC; ++i) fout[i].close();

  cout<<"\nRunning over "<<suffix<<" dataset with "<<selection<<" selection."<<"\n"<<endl;
  cout<<"| CUT DESCRIPTION                    |\t"<< "\t|"<<endl;
  cout<<"| Initial number of events:          |\t"<< nEvents[0]  <<"\t|"<<nEventsWeighted[0]  <<"\t|"<<endl;
  cout<<"| Pass HLT selection:                |\t"<< nEvents[1]  <<"\t|"<<nEventsWeighted[1]  <<"\t|"<<endl;
  cout<<"| Two oppositely charged leptons:    |\t"<< nEvents[2]  <<"\t|"<<nEventsWeighted[2]  <<"\t|"<<endl;
  cout<<"| Z mass window:                     |\t"<< nEvents[3]  <<"\t|"<<nEventsWeighted[3]  <<"\t|"<<endl;
  cout<<"| Third lepton veto:                 |\t"<< nEvents[4]  <<"\t|"<<nEventsWeighted[4]  <<"\t|"<<endl;
  cout<<"| Z Qt :                             |\t"<< nEvents[5]  <<"\t|"<<nEventsWeighted[5]  <<"\t|"<<endl;
  cout<<"| b-jet veto:                        |\t"<< nEvents[6]  <<"\t|"<<nEventsWeighted[6]  <<"\t|"<<endl;
  cout<<"|                            |\t"<< nEvents[7]  <<"\t|"<<nEventsWeighted[7]  <<"\t|"<<endl;
  cout<<"| H250                         |\t"<< nEvents[8]  <<"\t|"<<nEventsWeighted[8]  <<"\t|"<<endl;
  cout<<"| H300                       |\t"<< nEvents[9] <<"\t|"<<nEventsWeighted[9] <<"\t|"<<endl;
  cout<<"|                        |\t"<< nEvents[10]  <<"\t|"<<nEventsWeighted[10]  <<"\t|"<<endl;
  cout<<"| H400                        |\t"<< nEvents[11] <<"\t|"<<nEventsWeighted[11] <<"\t|"<<endl;
  cout<<"|                           |\t"<< nEvents[12] <<"\t|"<<nEventsWeighted[12] <<"\t|"<<endl;
  cout<<"| H500                       |\t"<< nEvents[13] <<"\t|"<<nEventsWeighted[13] <<"\t|"<<endl;
  cout<<"| H550                       |\t"<< nEvents[14] <<"\t|"<<nEventsWeighted[14] <<"\t|"<<endl;
  cout<<"| H600                       |\t"<< nEvents[15] <<"\t|"<<nEventsWeighted[15] <<"\t|"<<endl;
  
      
    //Normalization to 1fb and coloring the histograms
    string sel, sample;
    Float_t CS, nEv;  Int_t lColor, fColor;

    ifstream params;
    params.open("./params.txt");  //This is the file where parameters are stored
    params>>sel;
    params>>sample;  params>>nEv;  params>>CS;  params>>lColor;  params>>fColor;
    cout<<"sample: "<<sample<<"  cs: "<<CS<<"  events: "<<nEv<<" colors: "<<lColor<<"  "<<fColor<<endl;
    
    //histoFile->cd("");
    //samp = "none";
    //samp = sample;
    //cout<<samp<<endl;
    //TObject *obbb = new TObject(samp);
    //histoFile->WriteObjectAny(samp, TClass::TString);

    histoFile->Write();
    histoFile->cd("Andrey");
    scaleAndColor(sample.c_str(), CS, nEv, 1000.0, lColor, fColor); //normalize to 1fb

    histoFile->Close();  

    //reweightFile->Close();  
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
    //return fabs(nearestJet.DeltaPhi(metP4));
    return minDeltaPhi;
}

void higgsAnalyzer::PrintOutNoisy(Int_t num){
  nout[num]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<MET<<"\t* "
	   <<endl;//       <<isNoiseHcal<<"\t* "<<isDeadEcalCluster<<"\t* "<<isScraping<<"\t* "<<isCSCTightHalo<<"\t* "<<endl;
}

void higgsAnalyzer::PrintOut(Int_t num){
  fout[num]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<nVtx<<"\t* "<<nJets<<"\t* "
	   <<Mll<<"\t* "<<MT<<"\t* "<<MET<<"\t* "
	   <<qT<<"\t* "<<rhoFactor<<"\t* "<<endl;//<<muCountLoose<<"\t* "<<eleCountLoose<<endl;
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

      //Adding overflow and underflow bins inside the histogram:
      Int_t nBins    = hist->GetNbinsX();
      Double_t lastBin    = hist->GetBinContent(nBins);
      Double_t ovflBin    = hist->GetBinContent(nBins+1);
      Double_t lastBinErr = hist->GetBinError(nBins);
      Double_t ovflBinErr = hist->GetBinError(nBins+1);
     
      Double_t firstBin    = hist->GetBinContent(1);
      Double_t undflBin    = hist->GetBinContent(0);
      Double_t firstBinErr = hist->GetBinError(1);
      Double_t undflBinErr = hist->GetBinError(0);
  
      hist -> SetBinContent(nBins, lastBin+ovflBin);
      hist -> SetBinError(nBins, sqrt(pow(lastBinErr,2) + pow(ovflBinErr,2)) );
      hist -> SetBinContent(1, firstBin+undflBin);
      hist -> SetBinError(1, sqrt(pow(firstBinErr,2) + pow(undflBinErr,2)) );
      hist -> SetBinContent(0,0);
      hist -> SetBinContent(nBins+1,0);

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

void higgsAnalyzer::FillHistosNoise(Int_t num, Double_t weight){
  met0_et[num]      -> Fill(MET, weight);
}

void higgsAnalyzer::FillHistos(Int_t num, Double_t weight){
  
  met1_et[num]      -> Fill(pfMET, weight);
  met2_et[num]      -> Fill(pfMET1, weight);
  met3_et[num]      -> Fill(puCorrMET, weight);
  met4_et[num]      -> Fill(projMET, weight);
  met5_et[num]      -> Fill(ZprojMET, weight);
  met6_et[num]      -> Fill(redMET1, weight);
  met7_et[num]      -> Fill(redMET2, weight);
  met8_et[num]      -> Fill(compMET, weight);

  met1_over_qt[num] -> Fill(METqt, weight);
  met1_et_ovQt[num] -> Fill(MET, METqt, weight);
  met1_phi[num]     -> Fill(MET_phi, weight);

  l1_eta[num]  -> Fill(lep1_eta, weight);
  l1_phi[num]  -> Fill(lep1_phi, weight);
  l1_pt[num]   -> Fill(lep1_pt, weight);

  l2_eta[num]  -> Fill(lep2_eta, weight);
  l2_phi[num]  -> Fill(lep2_phi, weight);
  l2_pt[num]   -> Fill(lep2_pt, weight);

  Float_t lep_dPhi =  fabs(TVector2::Phi_mpi_pi(lep1_phi - lep2_phi));
  Float_t lep_dEta = fabs(lep1_eta - lep2_eta);
  Float_t lep_dR   = sqrt(pow(lep_dPhi,2) + pow(lep_dEta,2));
  l0_dPhi[num]        -> Fill(lep_dPhi,weight);
  l0_dEta[num]        -> Fill(lep_dEta,weight);
  l0_dR[num]          -> Fill(lep_dR,weight);
  l0_ptRatio[num]     -> Fill(lep2_pt/lep1_pt,weight);

  Float_t dPhiMetZ =  fabs(TVector2::Phi_mpi_pi(MET_phi-diPhi));
  di_dPhiMet[num]     -> Fill(dPhiMetZ,weight);
  met1_lg[num]        -> Fill(pfMET_lg,weight);
  met1_recoil_lg[num] -> Fill(pfMET_recoil,weight);

  //met4_over_qt[num] -> Fill(puCorrMETqt, weight);
  //met4_et_ovQt[num] -> Fill(puCorrMET, puCorrMETqt, weight);
  //met4_puSig[num]   -> Fill(puSigMET, weight);

  mt2[num]      -> Fill(MT, weight);
  mt2_met2[num] -> Fill(MT, MET, weight);

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
  jet_b_N[num]   -> Fill(nJetsB, weight);
  jet_b_Nssv[num]-> Fill(nJetsBssv, weight);
  jet_b_N25[num] -> Fill(nJetsB25, weight);
  jet_b_N30[num] -> Fill(nJetsB30, weight);

  jet_pt[num]     -> Fill(ptLeadJet, weight);
  jet_eta[num]    -> Fill(etaLeadJet, weight);
  jet_phi[num]    -> Fill(phiLeadJet, weight);
  jet_b_pt[num]   -> Fill(ptLeadBJet, weight);
  jet_dRlep1[num] -> Fill(dRjetlep1, weight); 
  jet_dRlep2[num] -> Fill(dRjetlep2, weight); 

  vtx_nPV_tot[num]    -> Fill(nVtxTotal, weight);
  vtx_nPV_raw[num]    -> Fill(nVtx);
  vtx_nPV_weight[num] -> Fill(nVtx, weight);
  vtx_ndof_1[num]     -> Fill(nDofVtx1, weight);
  vtx_ndof_2[num]     -> Fill(nDofVtx2, weight);

  if(isRealData)  run_events[num] -> Fill(runNumber);

  evt_byCut -> Fill(num, weight);  

  evt_weight[num] -> Fill(weight);

  if(nJets>0)
    {
      met1_dPhiClosJet1[num] -> Fill(dPhiClos1, weight);
      //  met1_dPhiClosJet2[num] -> Fill(dPhiClos2, weight);
      //met2_dPhiLeadJet1[num] -> Fill(dPhiLead1, weight);
      //met2_dPhiLeadJet2[num] -> Fill(dPhiLead2, weight);
    }

  // if(selection =="muGamma" || selection =="eGamma" || selection =="gamma"){
  
  ph_nGamma[num] -> Fill(nGamma);


  if (suffix.compare(0, 5, "ggHZZ") == 0 || suffix.compare(0, 5, "ggHWW") == 0) {
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
