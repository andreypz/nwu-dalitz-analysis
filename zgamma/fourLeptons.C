#define fourLeptons_cxx
#include "fourLeptons.h"
#include "../plugins/histo.h"

const UInt_t ntrig = 8;
string myTriggers[ntrig] = {
  "HLT_IsoMu24_v",
  "HLT_IsoMu24_eta2p1_v",
  "HLT_Mu13_Mu8_v",
  "HLT_Mu17_Mu8_v",
  "HLT_Mu17_TkMu8_v",
  "HLT_Mu22_TkMu8_v",
  "HLT_Mu22_Photon22_CaloIdL_v",
  "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v"
};

UInt_t nEventsTrig[nC][ntrig];
float EAPho[7][3] = {
  {0.012,  0.030,   0.148}, //         eta < 1.0
  {0.010,  0.057,   0.130}, // 1.0   < eta < 1.479
  {0.014,  0.039,   0.112}, // 1.479 < eta < 2.0
  {0.012,  0.015,   0.216}, // 2.0   < eta < 2.2
  {0.016,  0.024,   0.262}, // 2.2   < eta < 2.3
  {0.020,  0.039,   0.260}, // 2.3   < eta < 2.4
  {0.012,  0.072,   0.266}  // 2.4   < eta
};
string  myTrigger = "";
Float_t cut_l1pt  = 23;
Float_t cut_l2pt  = 4;
//const Float_t cut_iso_mu1 = 0.4, cut_iso_mu2=0;
Float_t cut_gammapt = 25;
Float_t mllMax = 200;

//Float_t global_Mll = 0;
Bool_t doRochCorr  = 1;
Bool_t makeApzTree = 1;
const UInt_t hisEVTS[] = {9331};
Int_t evSize = sizeof(hisEVTS)/sizeof(int);

bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());}

void fourLeptons::Begin(TTree * tree)
{
  cout<<"\n***      Begin the Analyzer      ****"<<endl;
  TString option = GetOption();

  cout<<"option = "<<option.Data()<<endl;
  TObjArray *args = (TObjArray*)option.Tokenize(" ");
  sample    = (string)((TObjString*)args->At(0))->GetString();
  selection = (string)((TObjString*)args->At(1))->GetString();
  trigger   = (string)((TObjString*)args->At(2))->GetString();
  gen       = (string)((TObjString*)args->At(3))->GetString();
  period    = "2012";


  cout<<sample<<"  "<<selection<<"  "<<trigger<<"  "<<gen<<endl;

  TH1::SetDefaultSumw2(kTRUE);


  vector<string>* triggerNames = 0;
  TFile   *inFile         = tree->GetCurrentFile();
  TTree   *jobTree        = (TTree*)inFile->Get("ntupleProducer/jobTree");
  jobTree->SetBranchAddress("triggerNames", &triggerNames);
  jobTree->GetEntry();

  //cout<<"Printing all the trigger names:"<<endl;
  //for (unsigned i = 0; i < triggerNames->size(); ++i) cout << triggerNames->at(i) << endl;

  myRandom = new TRandom3();

  triggerSelector = new TriggerSelector("", period, *triggerNames);
  weighter = new WeightUtils(sample, period, selection, 0);
  roch = new rochcor2012();

  histoFile = new TFile(Form("a_%s_higgsHistograms.root", selection.c_str()), "RECREATE");

  histoFile->mkdir("apzTree", "apzTree");

  hists = new HistManager(histoFile);

  for (Int_t n1=0; n1<nC; n1++){
    nEvents[n1]=0;
    for (UInt_t n2=0; n2<ntrig; n2++)
      nEventsTrig[n1][n2]=0;
  }

  //Cuts
  //electrons are two types: Barrel/Endcap
  //Loose
  elIdAndIsoCutsLoose.ptErrorOverPt[0] = 99999;
  elIdAndIsoCutsLoose.dEtaIn[0]        = 0.007;
  elIdAndIsoCutsLoose.dPhiIn[0]        = 0.15;
  elIdAndIsoCutsLoose.sigmaIetaIeta[0] = 0.01;
  elIdAndIsoCutsLoose.HadOverEm[0]     = 0.12;
  elIdAndIsoCutsLoose.dxy[0]           = 0.02;
  elIdAndIsoCutsLoose.dz[0]            = 0.2;
  elIdAndIsoCutsLoose.fabsEPDiff[0]    = 0.05;
  elIdAndIsoCutsLoose.pfIso04[0]       = 0.3;

  elIdAndIsoCutsLoose.ptErrorOverPt[1] = 99999;
  elIdAndIsoCutsLoose.dEtaIn[1]        = 0.009;
  elIdAndIsoCutsLoose.dPhiIn[1]        = 0.1;
  elIdAndIsoCutsLoose.sigmaIetaIeta[1] = 0.03;
  elIdAndIsoCutsLoose.HadOverEm[1]     = 0.1;
  elIdAndIsoCutsLoose.dxy[1]           = 0.02;
  elIdAndIsoCutsLoose.dz[1]            = 0.2;
  elIdAndIsoCutsLoose.fabsEPDiff[1]    = 0.05;
  elIdAndIsoCutsLoose.pfIso04[1]       = 0.3;

  //Tight
  elIdAndIsoCutsTight.ptErrorOverPt[0] = 9999.;
  elIdAndIsoCutsTight.dEtaIn[0] = 0.004;
  elIdAndIsoCutsTight.dPhiIn[0] = 0.06;
  elIdAndIsoCutsTight.sigmaIetaIeta[0] = 0.01;
  elIdAndIsoCutsTight.HadOverEm[0] = 0.12;
  elIdAndIsoCutsTight.fabsEPDiff[0] = 0.05;
  elIdAndIsoCutsTight.dxy[0] = 0.02;
  elIdAndIsoCutsTight.dz[0] = 0.1;
  elIdAndIsoCutsTight.pfIso04[0] = 0.15;

  elIdAndIsoCutsTight.ptErrorOverPt[1] = 9999.;
  elIdAndIsoCutsTight.dEtaIn[1] = 0.007;
  elIdAndIsoCutsTight.dPhiIn[1] = 0.03;
  elIdAndIsoCutsTight.sigmaIetaIeta[1] = 0.03;
  elIdAndIsoCutsTight.HadOverEm[1] = 0.10;
  elIdAndIsoCutsTight.fabsEPDiff[1] = 0.05;
  elIdAndIsoCutsTight.dxy[1] = 0.02;
  elIdAndIsoCutsTight.dz[1] = 0.1;
  elIdAndIsoCutsTight.pfIso04[1] = 0.15;

  //Soft muon id (specific for dalitz)
  muIdAndIsoCutsSoft.TrackLayersWithMeasurement = 5;
  muIdAndIsoCutsSoft.PixelLayersWithMeasurement = 1;
  muIdAndIsoCutsSoft.NormalizedChi2_tracker = 3;
  muIdAndIsoCutsSoft.dxy             = 0.2; //3
  muIdAndIsoCutsSoft.dz              = 0.5; //30
  muIdAndIsoCutsSoft.pfIso04         = 0.4;

  //Photon Id cuts
  phIdAndIsoCutsTight.PassedEleSafeVeto[0] = 1;
  phIdAndIsoCutsTight.HadOverEm[0]         = 0.05;
  phIdAndIsoCutsTight.sigmaIetaIeta[0]     = 0.011;
  phIdAndIsoCutsTight.chIso03[0] = 0.7;
  phIdAndIsoCutsTight.nhIso03[0] = 0.4;
  phIdAndIsoCutsTight.phIso03[0] = 0.5;

  phIdAndIsoCutsTight.PassedEleSafeVeto[1] = 1;
  phIdAndIsoCutsTight.HadOverEm[1]         = 0.05;
  phIdAndIsoCutsTight.sigmaIetaIeta[1]     = 0.031;
  phIdAndIsoCutsTight.chIso03[1] = 0.5;
  phIdAndIsoCutsTight.nhIso03[1] = 1.5;
  phIdAndIsoCutsTight.phIso03[1] = 1.0;


  phIdAndIsoCutsHZG.PassedEleSafeVeto[0] = 1;
  phIdAndIsoCutsHZG.HadOverEm[0]         = 0.05;
  phIdAndIsoCutsHZG.sigmaIetaIeta[0]     = 0.011;
  phIdAndIsoCutsHZG.chIso03[0] = 1.5;
  phIdAndIsoCutsHZG.nhIso03[0] = 1.0;
  phIdAndIsoCutsHZG.phIso03[0] = 0.7;

  phIdAndIsoCutsHZG.PassedEleSafeVeto[1] = 1;
  phIdAndIsoCutsHZG.HadOverEm[1]         = 0.05;
  phIdAndIsoCutsHZG.sigmaIetaIeta[1]     = 0.033;
  phIdAndIsoCutsHZG.chIso03[1] = 1.2;
  phIdAndIsoCutsHZG.nhIso03[1] = 1.5;
  phIdAndIsoCutsHZG.phIso03[1] = 1.0;


  //Angles class:
  ang = new ZGAngles();

  fout.open("./out_synch_.txt",ofstream::out);
  fout.precision(3); fout.setf(ios::fixed, ios::floatfield);

  histoFile->cd("apzTree");
  if (makeApzTree){
    _apzTree = new TTree("apzTree", "A tree for studying new particles");
    _apzTree->Branch("dr1234",&apz_dr1234,"dr1234/D");
    _apzTree->Branch("dr12",  &apz_dr12,  "dr12/D");
    _apzTree->Branch("dr34",  &apz_dr34,  "dr34/D");
    _apzTree->Branch("pt12",  &apz_pt12,  "pt12/D");
    _apzTree->Branch("pt34",  &apz_pt34,  "pt34/D");
    _apzTree->Branch("m12",   &apz_m12,   "m12/D");
    _apzTree->Branch("m34",   &apz_m34,   "m34/D");
    _apzTree->Branch("m4l",   &apz_m4l,   "m4l/D");
    _apzTree->Branch("pt1",   &apz_pt1,   "pt1/D");
    _apzTree->Branch("pt2",   &apz_pt2,   "pt2/D");
    _apzTree->Branch("pt3",   &apz_pt3,   "pt3/D");
    _apzTree->Branch("pt4",   &apz_pt4,   "pt4/D");

    _apzTree->Branch("run",  &runNumber,  "run/i");
    _apzTree->Branch("event",&eventNumber,"event/l");
  }

  histoFile->cd();

}

void fourLeptons::SlaveBegin(TTree * /*tree*/)
{   TString option = GetOption(); }

Bool_t fourLeptons::Process(Long64_t entry)
{
  GetEntry(entry);
  CountEvents(0);
  FillHistoCounts(0, 1);
  if (nEvents[0] % (int)5e4 == 0) cout<<nEvents[4]<<" events passed of "<<nEvents[0]<<" checked!"<<endl;

  if (nEvents[0] == 1) weighter->SetDataBit(isRealData);

  //if (nEvents[0]>20) Abort("enough is enosugh!");

  //////////////////
  //Trigger status//
  //////////////////
  Bool_t triggerPass  = 1;
  Bool_t isFound  = 0;
  Int_t  prescale = 99;
  Float_t eventWeight = 1.0;
  for (UInt_t i =0; i<ntrig; i++)
    {
      triggerSelector->SelectTrigger(myTriggers[i], triggerStatus, hltPrescale, isFound, triggerPass, prescale);
      if(triggerPass) nEventsTrig[0][i]++;
    }

  eventWeight = weighter->PUWeight(nPUVerticesTrue);
  //cout<<" weight="<<eventWeight<<endl;

  //----------------------//
  // Selecting triggers   //
  //----------------------//

  Bool_t checkTrigger = kTRUE;
  if (trigger=="mumu")
    {
      myTrigger = "HLT_Mu13_Mu8_v";
      //myTrigger = "HLT_Mu17_Mu8_v";
      cut_l1pt = 15;
      //cut_l2pt = 9;
      //cut_gammapt = 23;
    }
  else if (trigger=="mugamma")
    {
      myTrigger = "HLT_Mu22_Photon22_CaloIdL_v";
      //cut_l1pt = 23;
      //cut_l2pt = 7;
      //cut_gammapt = 25;
    }
  else if (trigger=="single-mu")
    {
      myTrigger = "HLT_IsoMu24_eta2p1_v";
      cut_l1pt = 25;
      cut_l2pt = 4;
      cut_gammapt = 23;
    }
  else if (trigger=="ee")
    {
      myTrigger = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v";
      cut_l1pt = 18;
      cut_l2pt = 4;
      cut_gammapt = 23;
    }
  else
    checkTrigger = kFALSE;
  //cout<<"Other options are not supported!"<<endl;

  //cout<<"My trigger is "<<myTrigger<<endl;


  nVtx = primaryVtx->GetSize();
  // --- Primary vertex ---//

  if (primaryVtx->GetSize()<1) return kTRUE;
  TCPrimaryVtx *mainPrimaryVertex = 0, *secondPrimaryVertex = 0;
  mainPrimaryVertex   = (TCPrimaryVtx*)(primaryVtx->At(0));
  TVector3* pvPosition;// = new TVector3();

  nDofVtx1 = mainPrimaryVertex->NDof();
  if(primaryVtx->GetSize() > 1){
    secondPrimaryVertex = (TCPrimaryVtx*)(primaryVtx->At(1));
    nDofVtx2 = secondPrimaryVertex->NDof();
  }
  pvPosition = mainPrimaryVertex;

  vector<TCElectron> electrons0, electrons, electrons_dalitz, fake_electrons0;
  vector<TCMuon> muons0, muons;
  vector<TCPhoton> photons0, photonsTight, photonsHZG, fake_photons;
  TCPhysObject l1,l2, lPt1, lPt2, lPt3, lPt4;
  TCPhoton gamma0, ufoton, ufelectron;
  //TCPhoton gamma;
  TCPhysObject gamma;

  FillHistoCounts(2, eventWeight);
  CountEvents(2);


  for (Int_t i = 0; i < recoPhotons->GetSize(); ++i) {
    //cout<<"new photon!!!!!!!"<<endl;
    TCPhoton* thisPhoton = (TCPhoton*) recoPhotons->At(i);
    //cout<<"in photon loop event = "<<eventNumber
    //  <<"\n pt="<<thisPhoton->Pt()<<" eta="<<thisPhoton->Eta()<<" phi="<<thisPhoton->Phi()<<" px="<<thisPhoton->Px()<<endl;

    //if (!(fabs(thisPhoton->Eta()) < 2.5)) continue;
    if (!(fabs(thisPhoton->SCEta()) < 2.5)) continue;

    //Double_t corrPhoEnSC   = correctedPhotonEnergy(thisPhoton->SCEnergy(), thisPhoton->SCEta(), thisPhoton->R9(), runNumber,
    //0, "Moriond2014", !isRealData, myRandom);
    Double_t corrPhoEnReco = correctedPhotonEnergy(thisPhoton->E(), thisPhoton->SCEta(), thisPhoton->R9(), runNumber,
						   0, "Moriond2014", !isRealData, myRandom);

    //cout<<runNumber<<"  uncor en="<< thisPhoton->E()<<"  cor en = "<<corrPhoEn<<endl;
    //Double_t scale = 1.0;
    Double_t scale = corrPhoEnReco/thisPhoton->E();
    //Double_t scale = corrPhoEnSC/thisPhoton->E();

    thisPhoton->SetXYZM(scale*thisPhoton->Px(), scale*thisPhoton->Py(),scale*thisPhoton->Pz(),0);

    //cout<<runNumber<<"  uncor en="<< thisPhoton->E()<<"  cor en = "<<corrPhoEn<<endl;
    //thisPhoton->SetPtEtaPhiM(corrPhoEnSC/cosh(SCEta()));

    if(thisPhoton->Pt() > 15){
      if(PassPhotonIdAndIso(thisPhoton, phIdAndIsoCutsHZG, pvPosition)
	 )
	{
	  photonsHZG.push_back(*thisPhoton);
	}

      if(PassPhotonIdAndIso(thisPhoton, phIdAndIsoCutsTight, pvPosition)
	 )
	photonsTight.push_back(*thisPhoton);
    }
  }


  sort(photons0.begin(),     photons0.end(),     P4SortCondition);
  sort(photonsTight.begin(), photonsTight.end(), P4SortCondition);
  sort(photonsHZG.begin(),   photonsHZG.end(),   P4SortCondition);


  //if (!isRealData)
  //eventWeight *= weighter->PhotonSF(gamma);

  for (Int_t i = 0; i <  recoElectrons->GetSize(); ++i) {
    TCElectron* thisElec = (TCElectron*) recoElectrons->At(i);
    if (!(fabs(thisElec->Eta()) < 2.5)) continue;
    if (thisElec->Pt() > 5.0){
      if (
	  PassElectronIdAndIso(thisElec, elIdAndIsoCutsLoose, pvPosition)
	  //PassElectronIdAndIso(thisElec, elIdAndIsoCutsTight, pvPosition)
	  )
	{
	  //cout<<"event = "<<eventNumber
	  //  <<"\n pt="<<thisElec->Pt()<<" eta="<<thisElec->Eta()<<" phi="<<thisElec->Phi()<<" px="<<thisElec->Px()<<endl;

	  electrons0.push_back(*thisElec);
	}
      //else
      //fourLeptons::ElectronDump(*thisElec, elIdAndIsoCutsLoose, pvPosition);

    }
  }
  sort(electrons0.begin(), electrons0.end(), P4SortCondition);
  sort(electrons.begin(),  electrons.end(),  P4SortCondition);


  for (Int_t i = 0; i < recoMuons->GetSize(); ++ i) {
    TCMuon* thisMuon = (TCMuon*) recoMuons->At(i);
    //cout<<"event = "<<eventNumber
    //  <<"\n pt="<<thisMuon->Pt()<<" eta="<<thisMuon->Eta()<<" phi="<<thisMuon->Phi()<<" px="<<thisMuon->Px()<<endl;

    if (!(fabs(thisMuon->Eta()) < 2.4)) continue;

    if(thisMuon->Pt() > 4)
      for(Int_t ev=0; ev<evSize;ev++){
	if (eventNumber==hisEVTS[ev]){cout<<"   mu>4 cut ---> "<<eventNumber<<" <--- Found an event after cut "<<endl;
	  MuonDump(*thisMuon, pvPosition);
	  break;}}

    if(thisMuon->Pt() > 4
       && PassMuonIdAndIso(thisMuon, muIdAndIsoCutsSoft, pvPosition)
       ) {
      if (doRochCorr){
	Float_t qter = 1;
	//cout<<"DBG before rochcor"<<endl;
	if (isRealData)
	  roch->momcor_data(*thisMuon, thisMuon->Charge(), 0, qter);
	else
	  roch->momcor_mc(*thisMuon,  thisMuon->Charge(), 0, qter );

	//cout<<"DBG after rochcor"<<endl;
      }

      muons.push_back(*thisMuon);
    }
  }

  sort(muons.begin(), muons.end(), P4SortCondition);

  hists->fill1DHist(muons.size(),     Form("size_mu_cut%i",  1),";Number of muons",    5,0,5, 1, "Muons");
  //hists->fill1DHist(electrons.size(), Form("size_el_cut%i",  1),";Number of electrons",5,0,5, 1, "Electrons");
  //hists->fill1DHist(electrons0.size(),Form("size_el0_cut%i", 1),";Number of electrons",5,0,5, 1, "Electrons");
  //hists->fill1DHist(photons.size(),   Form("size_ph_cut%i",  1),";Number of photons",  5,0,5, 1, "Photon");
  //hists->fill1DHist(photons0.size(),  Form("size_ph0_cut%i", 1),";Number of photons",  5,0,5, 1, "Photon");


  if (checkTrigger){
    triggerSelector->SelectTrigger(myTrigger, triggerStatus, hltPrescale, isFound, triggerPass, prescale);
    if (!triggerPass) return kTRUE;
  }

  FillHistoCounts(3, eventWeight);
  CountEvents(3);


  if (selection=="2e2mu")
    {

      if (muons.size()<2) return kTRUE;

      Float_t tmpMll = 99999;
      Float_t MAPZ   = 0;

      if (muons.size()==2){
	lPt1 = muons[0];
	lPt2 = muons[1];
      }
      else
	for (UInt_t m1=0; m1 < muons.size(); m1++){
	  for (UInt_t m2=m1; m2 < muons.size(); m2++){
	    if (muons[m1].Charge()*muons[m2].Charge()==1) continue; //same charge - don't care

	    if( fabs((muons[m1]+muons[m2]).M() - MAPZ) < tmpMll)
	      {
		lPt1 = muons[m1];
		lPt2 = muons[m2];
		tmpMll = (muons[m1]+muons[m2]).M();
	      }
	  }
	}

      if (lPt1.Charge()*lPt2.Charge() == 1)
	return kTRUE;//Abort("reco * They are the same charge!");

      //if (fourLeptons::CalculateMuonIso(&muons[0]) > 0.4)
      //return kTRUE;
      //if (lPt2.Pt() > 20 && fourLeptons::CalculateMuonIso(&muons[1]) > 0.4)
      //return kTRUE;

      FillHistoCounts(4, eventWeight);
      CountEvents(4);

      if(electrons0.size()<2) return kTRUE;

      FillHistoCounts(5, eventWeight);
      CountEvents(5);

      lPt3 = electrons0[0];
      lPt4 = electrons0[1];
    }
  else if (selection == "4mu")
    {

      if (muons.size()<4) return kTRUE;

      Float_t tmpMll = 99999;
      Float_t MAPZ   = 0;

      for (UInt_t m1=0; m1 < muons.size(); m1++){
	for (UInt_t m2=m1; m2 < muons.size(); m2++){
	  if (muons[m1].Charge()*muons[m2].Charge()==1) continue; //same charge - don't care

	  if( fabs((muons[m1]+muons[m2]).M() - MAPZ) < tmpMll)
	    {
	      lPt1 = muons[m1];
	      lPt2 = muons[m2];
	      tmpMll = (muons[m1]+muons[m2]).M();
	    }
	}

      }


      Abort("Yet to be done");
    }
  else
    Abort("Selection is not supported");

  if (lPt1.Pt() < 25 || lPt2.Pt() < 4)   return kTRUE;
  if (lPt3.Pt() < 25 || lPt4.Pt() < 5)   return kTRUE;

  FillHistoCounts(6, eventWeight);
  CountEvents(6);


  if (!isfinite(eventWeight))
    cout<<"nan/inf "<<eventWeight<<endl;
  //    if (fabs(eventWeight)>50 || fabs(eventWeight)<0.00001)

  if (!isRealData){
    eventWeight *= weighter->MuonSF(lPt1);
    eventWeight *= weighter->MuonSF(lPt2);
  }


  //for(Int_t ev=0; ev<evSize;ev++){
  //if (eventNumber==hisEVTS[ev]){cout<<eventNumber<<" Found an event after lept selection "<<endl;
  //  cout<<"Mll = "<<Mll<<endl;
  //  break;}
  //}

  Float_t Mllg = 0;//(l1+l2+gamma).M();
  Float_t M4l  = (l1+l2+lPt3+lPt4).M();

  fout<<runNumber<<" "<<eventNumber<<endl;

  Double_t Mll = (lPt1+lPt2).M();
  hists->fill1DHist(Mll,  Form("diLep_mass_low_cut%i",  3), ";M(ll)", 50, 0,20,  1, "");
  hists->fill1DHist(Mll,  Form("diLep_mass_high_cut%i", 3), ";M(ll)", 50, 0,120, 1, "");


  if (Mll > mllMax) return kTRUE;

  FillHistoCounts(7, eventWeight);
  CountEvents(7);

  if (makeApzTree){

    apz_pt1 = lPt1.Pt();
    apz_pt2 = lPt2.Pt();
    apz_pt3 = lPt3.Pt();
    apz_pt4 = lPt4.Pt();

    apz_dr1234 = (lPt1+lPt2).DeltaR(lPt3+lPt4);

    apz_dr12 = lPt1.DeltaR(lPt2);
    apz_dr34 = lPt3.DeltaR(lPt4);

    apz_pt12 = (lPt1+lPt2).Pt();
    apz_pt34 = (lPt3+lPt4).Pt();

    apz_m12 = (lPt1+lPt2).M();
    apz_m34 = (lPt3+lPt4).M();

    apz_m4l = (lPt1+lPt2+lPt3+lPt4).M();

    _apzTree->Fill();

  }
  return kTRUE;


}

void fourLeptons::SlaveTerminate() {}

void fourLeptons::Terminate()
{
  cout<<" ** FOR PAS **"<<endl;
  cout<<"| CUT DESCRIPTION            |\t"<< "\t|"<<endl;
  cout<<"| 0: initial                 |\t"<< nEvents[0]  <<"\t|"<<float(nEvents[0])/nEvents[0]<<"\t|"<<endl;
  cout<<"| 1: n/a                     |\t"<< nEvents[1]  <<"\t|"<<float(nEvents[1])/nEvents[0]<<"\t|"<<endl;
  cout<<"| 2: vtx                     |\t"<< nEvents[2]  <<"\t|"<<float(nEvents[2])/nEvents[1]<<"\t|"<<endl;
  cout<<"| 3: trigger                 |\t"<< nEvents[3]  <<"\t|"<<float(nEvents[3])/nEvents[2]<<"\t|"<<endl;
  cout<<"| 4: 2 muons                 |\t"<< nEvents[4]  <<"\t|"<<float(nEvents[4])/nEvents[3]<<"\t|"<<endl;
  cout<<"| 5: 2 electrons             |\t"<< nEvents[5]  <<"\t|"<<float(nEvents[5])/nEvents[4]<<"\t|"<<endl;
  cout<<"| 6: 2+2 and pt Cuts         |\t"<< nEvents[6]  <<"\t|"<<float(nEvents[6])/nEvents[5]<<"\t|"<<endl;
  cout<<"| 7:   Mll < 25              |\t"<< nEvents[7]  <<"\t|"<<float(nEvents[7])/nEvents[6]<<"\t|"<<endl;
  cout<<"| 8:  n/a                    |\t"<< nEvents[8]  <<"\t|"<<float(nEvents[8])/nEvents[7]<<"\t|"<<endl;

  cout<<"\n\n | Trigger efficiency 3              |\t"<<endl;
  for(UInt_t n=0; n<ntrig; n++){
    UInt_t N=7;
    cout<<n<<"  "<<float(nEventsTrig[N][n])/nEvents[N]<<"   "<<myTriggers[n]<<endl;;
  }

  hists->fill1DHist(-1, "evt_byCut",";cut #;weighted events", nC+1,-1,nC, totEvents, "Counts");
  hists->fill1DHist(-1, "evt_byCut_raw", ";cut #;events",     nC+1,-1,nC, totEvents, "Counts");

  //gDirectory->Print();
  //gDirectory->FindObject("evt_byCut_raw")->Print();
  histoFile->cd();
  //gDirectory->Print();
  //cout<<"Raw events:       "<<((TH1F*)gDirectory->FindObject("evt_byCut_raw"))<<endl;
  //cout<<"Weighted events:  "<<((TH1F*)histoFile->Get("evt_byCut"))->GetBinContent(1)<<endl;
  //cout<<"weight integral = "<<((TH1F*)histoFile->Get("weight__cut5"))->Integral()<<endl;

  histoFile->Write();
  histoFile->Close();
  cout<<" ** End of Analyzer **"<<endl;
}


void fourLeptons::CountEvents(Int_t num)
{
  nEvents[num]++;
}



float fourLeptons::CalculateElectronIso(TCElectron *lep)
{
  float eleISO = 0;
  eleISO = (lep->PfIsoCharged() +
	    TMath::Max(0.0, (Double_t)(lep->PfIsoPhoton() +
				       lep->PfIsoNeutral() -
				       rho25Factor*lep->EffArea())))/lep->Pt();
  return eleISO;
}



bool fourLeptons::PassElectronIdAndIso(TCElectron *lep, elIdAndIsoCuts cuts, TVector3 *pv)
{
  bool pass = false;

  Float_t eleISO = CalculateElectronIso(lep);

  if(((fabs(lep->Eta()) < 1.442
       && lep->PtError()/lep->Pt() < cuts.ptErrorOverPt[0]
       && lep->SigmaIEtaIEta()     < cuts.sigmaIetaIeta[0]
       && fabs(lep->SCDeltaEta())  < cuts.dEtaIn[0]
       && fabs(lep->SCDeltaPhi())  < cuts.dPhiIn[0]
       && lep->HadOverEm()         < cuts.HadOverEm[0]
       && fabs(lep->Dxy(pv))       < cuts.dxy[0]
       && fabs(lep->Dz(pv))        < cuts.dz[0]
       //&& fabs(lep->InverseEnergyMomentumDiff()) < cuts.fabsEPDiff[0]
       //&& eleISO < cuts.pfIso04[0]
       )||
      (fabs(lep->Eta()) > 1.556
       && lep->PtError()/lep->Pt() < cuts.ptErrorOverPt[1]
       && lep->SigmaIEtaIEta()     < cuts.sigmaIetaIeta[1]
       && fabs(lep->SCDeltaEta())  < cuts.dEtaIn[1]
       && fabs(lep->SCDeltaPhi())  < cuts.dPhiIn[1]
       && lep->HadOverEm()         < cuts.HadOverEm[1]
       && fabs(lep->Dxy(pv))       < cuts.dxy[1]
       && fabs(lep->Dz(pv))        < cuts.dz[1]
       // && fabs(lep->InverseEnergyMomentumDiff()) < cuts.fabsEPDiff[1]
       //&& eleISO < cuts.pfIso04[1]
       ))
     //&& lep->PassConversionVeto()
     ) pass = true;

  //pass = true;

  //cout<<"ele pass? "<<pass<<endl;
  return pass;
}


bool fourLeptons::PassElectronIdAndIsoMVA(TCElectron *lep)
{
  bool pass = false;

  Float_t eleISO = CalculateElectronIso(lep);

  if(
     eleISO < 0.15
     && lep->ConversionMissHits()==0
     && lep->PassConversionVeto() == true
     &&
     (
      (lep->Pt()>10 && lep->Pt()<20 &&
       (
	(fabs(lep->Eta()) < 0.8
	 && lep->MvaID() > 0.00)
	||
	(fabs(lep->Eta()) >  0.8 && fabs(lep->Eta()) < 1.479
	 && lep->MvaID() > 0.10)
	||
	(fabs(lep->Eta()) >  1.479 && fabs(lep->Eta()) < 2.5
	 && lep->MvaID() > 0.62)
	)
       )
      ||
      (lep->Pt()>20 &&
       (
	(fabs(lep->Eta()) < 0.8
	 && lep->MvaID() > 0.94)
	||
	(fabs(lep->Eta()) >  0.8 && fabs(lep->Eta()) < 1.479
	 && lep->MvaID() > 0.85)
	||
	(fabs(lep->Eta()) >  1.479 && fabs(lep->Eta()) < 2.5
	 && lep->MvaID() > 0.92)
	)
       )
      )
     )
    pass = true;

  return pass;
}

float fourLeptons::CalculateMuonIso(TCMuon *lep)
{
  float muISO = (lep->PfIsoCharged() + TMath::Max(0.0, lep->PfIsoPhoton() + lep->PfIsoNeutral() - 0.5*lep->PfIsoPU()))/lep->Pt();

  //  Float_t muISO = (lep->IsoMap("pfChargedHadronPt_R04") +
  //TMath::Max(0.0, lep->IsoMap("pfNeutralHadronEt_R04") + lep->IsoMap("pfPhotonEt_R04") - 0.5*lep->IsoMap("pfPUPt_R04")))/lep->Pt();

  return muISO;
}


bool fourLeptons::PassMuonIdAndIso(TCMuon *lep, muIdAndIsoCuts cuts, TVector3 *pv)
{
  bool pass = false;

  //Float_t muISO =  CalculateMuonIso(lep);

  if(1
     && lep->IsPF()
     && ((lep->IsTRK() && lep->NumberOfMatches() > 0) || lep->IsGLB())
     && fabs(lep->Dxy(pv))     < cuts.dxy
     && fabs(lep->Dz(pv))      < cuts.dz
     && (lep->Pt() < 20 || fourLeptons::CalculateMuonIso(lep) < 0.4)
     )
    pass = true;
  /*
  if(1
     && lep->IsPF() && lep->IsTRK()
     && lep->NormalizedChi2_tracker()  < cuts.NormalizedChi2_tracker
     && fabs(lep->Dxy(pv))     < cuts.dxy
     && fabs(lep->Dz(pv))      < cuts.dz
     )
    pass = true;
  */

  return pass;
}

void fourLeptons::CalculatePhotonIso(TCPhoton *ph, float& chIsoCor, float& nhIsoCor, float& phIsoCor){
  float chEA,nhEA,phEA,tmpEta;
  //float chEA,nhEA,phEA,chIsoCor,nhIsoCor,phIsoCor,tmpEta;
   tmpEta = ph->SCEta();

  if (fabs(tmpEta) < 1.0){
    chEA = EAPho[0][0];
    nhEA = EAPho[0][1];
    phEA = EAPho[0][2];
  }else if (fabs(tmpEta) < 1.479){
    chEA = EAPho[1][0];
    nhEA = EAPho[1][1];
    phEA = EAPho[1][2];
  }else if (fabs(tmpEta) < 2.0){
    chEA = EAPho[2][0];
    nhEA = EAPho[2][1];
    phEA = EAPho[2][2];
  }else if (fabs(tmpEta) < 2.2){
    chEA = EAPho[3][0];
    nhEA = EAPho[3][1];
    phEA = EAPho[3][2];
  }else if (fabs(tmpEta) < 2.3){
    chEA = EAPho[4][0];
    nhEA = EAPho[4][1];
    phEA = EAPho[4][2];
  }else if (fabs(tmpEta) < 2.4){
    chEA = EAPho[5][0];
    nhEA = EAPho[5][1];
    phEA = EAPho[5][2];
  }else{
    chEA = EAPho[6][0];
    nhEA = EAPho[6][1];
    phEA = EAPho[6][2];
  }
  chIsoCor = ph->PfIsoCharged()-rhoFactor*chEA;
  nhIsoCor = ph->PfIsoNeutral()-rhoFactor*nhEA;
  phIsoCor = ph->PfIsoPhoton()-rhoFactor*phEA;
  //cout<<chIsoCor<<"  "<<nhIsoCor<<"  "<<phIsoCor<<endl;

}

bool fourLeptons::PassPhotonIdAndIso(TCPhoton *ph, phIdAndIsoCuts cuts, TVector3 *pv)
{
  bool pass = false;
  //Float_t phoISO = 0;
  float tmpEta = ph->SCEta();

  float chIsoCor=0, nhIsoCor=0, phIsoCor=0;
  CalculatePhotonIso(ph, chIsoCor,nhIsoCor,phIsoCor);
  if(
     (fabs(tmpEta)  < 1.442
      && ph->ConversionVeto()       == cuts.PassedEleSafeVeto[0]
      && ph->HadOverEm()             < cuts.HadOverEm[0]
      && ph->SigmaIEtaIEta()         < cuts.sigmaIetaIeta[0]
      && max((double)chIsoCor,0.)    < cuts.chIso03[0]
      && max((double)nhIsoCor,0.)    < cuts.nhIso03[0] + 0.040*ph->Pt()
      && max((double)phIsoCor,0.)    < cuts.phIso03[0] + 0.005*ph->Pt()
      ) ||
     (fabs(tmpEta)  > 1.566
      && ph->ConversionVeto()       == cuts.PassedEleSafeVeto[1]
      && ph->HadOverEm()             < cuts.HadOverEm[1]
      && ph->SigmaIEtaIEta()         < cuts.sigmaIetaIeta[1]

      && max((double)chIsoCor,0.)    < cuts.chIso03[1]
      && max((double)nhIsoCor,0.)    < cuts.nhIso03[1] + 0.040*ph->Pt()
      && max((double)phIsoCor,0.)    < cuts.phIso03[1] + 0.005*ph->Pt()
      )
     ) pass = true;

  return pass;
}


void fourLeptons::FillHistoCounts(Int_t num, Double_t weight)
{
  hists->fill1DHist(num, "evt_byCut",     ";cut;weighted events", nC+1, -1,nC, weight, "Counts");
  hists->fill1DHist(num, "evt_byCut_raw", ";cut;events",          nC+1, -1,nC,      1, "Counts");
}

TCGenParticle * fourLeptons::GetPrimaryAncestor(TCGenParticle *p)
{
  TCGenParticle *a = p;
  while (a->Mother())
    a = a->Mother();
  return a;
}


void fourLeptons::DiscoverGeneology(TCGenParticle *p, ULong64_t ev)
{
  cout<<"---->> event = "<<ev<<"   "<<p->GetPDGId()<<"  st = "<<p->GetStatus()
      <<"\n pt="<<p->Pt()<<" eta="<<p->Eta()<<" phi="<<p->Phi()<<" M="<<p->M()<<endl;
  if(p->Mother()){
    TCGenParticle *m = p->Mother();
    cout<<"    mother = "<<m->GetPDGId()<<"  st="<<m->GetStatus()<<"\n pt="<<m->Pt()<<endl;
    if(m->Mother()){
      TCGenParticle *g = m->Mother();
      cout<<"  grandmother = "<<g->GetPDGId()<<"  st="<<g->GetStatus()<<"\n pt="<<g->Pt()<<endl;

      if(g->Mother()){
	TCGenParticle *gg = g->Mother();
	cout<<"  grand-grandmother = "<<gg->GetPDGId()<<"  st="<<gg->GetStatus()<<"\n pt="<<gg->Pt()<<endl;
	if(gg->Mother()){
	  TCGenParticle *ggg = gg->Mother();
	  cout<<"  grand-grand-grandmother = "<<ggg->GetPDGId()<<"  st="<<ggg->GetStatus()<<"\n pt="<<ggg->Pt()<<endl;

	  if(ggg->Mother()){
	    TCGenParticle *gggg = ggg->Mother();
	    cout<<"  grand-grand-grand-grandmother = "<<gggg->GetPDGId()<<"  st="<<gggg->GetStatus()<<"\n pt="<<gggg->Pt()<<endl;
	  }
	}
      }
    }
  }

}



void fourLeptons::MuonDump(TCMuon mu, TVector3 *pv)
{
  Float_t muISO = fourLeptons::CalculateMuonIso(&mu);
  //Float_t muISO = (mu.IsoMap("pfChargedHadronPt_R04") +
  //               TMath::Max(0.0, mu.IsoMap("pfNeutralHadronEt_R04") + mu.IsoMap("pfPhotonEt_R04") - 0.5*mu.IsoMap("pfPUPt_R04")))/mu.Pt();


  cout  << runNumber << " " << eventNumber << "  pt =" << mu.Pt()
	<< " eta=" << mu.Eta() << " GLB=" << mu.IsGLB() << "  isPF=" << mu.IsPF() << " isTRK="<<mu.IsTRK()<<endl;
  cout  << "chi2_tr="<< mu.NormalizedChi2_tracker() << "  iso=" << muISO <<  "  dxy=" << mu.Dxy(pv) << " dz=" << mu.Dz(pv)	<< "\n  good =" << mu.IsGood()
    //<< "\n  trk =" << mu.TrackLayersWithMeasurement() << "  pix=" <<mu.PixelLayersWithMeasurement()
    //   << "\n " << mu.NormalizedChi2() << " " << mu.NumberOfValidMuonHits() << " " << mu.NumberOfMatchedStations()
    //  << " " << mu.NumberOfValidPixelHits() << " " << mu.TrackLayersWithMeasurement()
	<< endl;
}





void fourLeptons::PhotonDump(TCPhoton pho, phIdAndIsoCuts cuts)
{
  //Float_t phoISO = 0;
  float chIsoCor,nhIsoCor,phIsoCor;
  float tmpEta = pho.SCEta();

  CalculatePhotonIso(&pho, chIsoCor,nhIsoCor,phIsoCor);

  cout  << runNumber << " " << eventNumber << " " << pho.Pt() << " " << pho.Eta() << " " << pho.SCEta() <<endl;
  cout<<"conv:"<<pho.ConversionVeto()<<"  hoem:"<< pho.HadOverEm() <<"  sieie:"<<pho.SigmaIEtaIEta()<<endl;
  cout  <<"\n iso: ch:" << chIsoCor<< "  nh:" << nhIsoCor <<  " ph:"<<phIsoCor<<endl;
  if (fabs(tmpEta)  < 1.442){
    cout<<"\n cuts: "<<cuts.chIso03[0] <<"  "<<cuts.nhIso03[0] + 0.040*pho.Pt()<<"  "<<cuts.phIso03[0] + 0.005*pho.Pt()<<endl;
    cout<<"conv:"<<cuts.PassedEleSafeVeto[0]<<"  hoem:"<< cuts.HadOverEm[0]<<"  sieie:"<<cuts.sigmaIetaIeta[0]<<endl;
  }
  else if (fabs(tmpEta)  > 1.566){
    cout<<"\n cuts: "<<cuts.chIso03[1] <<"  "<<cuts.nhIso03[1] + 0.040*pho.Pt()<<"  "<<cuts.phIso03[1] + 0.005*pho.Pt()<<endl;
    cout<<"conv:"<<cuts.PassedEleSafeVeto[1]<<"  hoem:"<< cuts.HadOverEm[1]<<"  sieie:"<<cuts.sigmaIetaIeta[1]<<endl;
  }

}


void fourLeptons::ElectronDump(TCElectron ele, elIdAndIsoCuts cuts, TVector3 *pv)
{
  cout  << runNumber << " " << eventNumber << " " << ele.Pt() << " " << ele.Eta() << " " << ele.SCEta() <<endl;
  Float_t eleISO = 0;//CalculateElectronIso(*ele);

  cout<<"conv veto                  "<<ele.PassConversionVeto()<<"  "<<""<<endl;
  cout<<"ptErrorOverPt              "<<ele.PtError()/ele.Pt()<<"  "<<cuts.ptErrorOverPt[ele.IsEB()]<<endl;
  cout<<"ele.SigmaIEtaIEta()        "<<ele.SigmaIEtaIEta() <<"  "<<cuts.sigmaIetaIeta[ele.IsEB()] <<endl;
  cout<<"fabs(ele.SCDeltaEta())     "<<fabs(ele.SCDeltaEta()) <<"  "<<cuts.dEtaIn[ele.IsEB()] <<endl;
  cout<<"fabs(ele.SCDeltaPhi())     "<<fabs(ele.SCDeltaPhi()) <<"  "<<cuts.dPhiIn[ele.IsEB()] <<endl;
  cout<<"ele.HadOverEm()            "<<ele.HadOverEm() <<"  "<< cuts.HadOverEm[ele.IsEB()] <<endl;
  cout<<"fabs(ele.Dxy(pv))          "<<fabs(ele.Dxy(pv)) <<"  "<<cuts.dxy[ele.IsEB()] <<endl;
  cout<<"fabs(ele.Dz(pv))           "<<fabs(ele.Dz(pv)) <<"  "<<cuts.dz[ele.IsEB()] <<endl;
  cout<<"InverseEnergyMomentumDiff  "<<fabs(ele.InverseEnergyMomentumDiff()) <<"  "<<cuts.fabsEPDiff[ele.IsEB()] <<endl;
  cout<<"eleISO  "<<eleISO <<"  "<<cuts.pfIso04[ele.IsEB()] <<endl;
  cout<<"***** \n"<<endl;
}
