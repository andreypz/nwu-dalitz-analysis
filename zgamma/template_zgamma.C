#define zgamma_cxx
#include "zgamma.h"
#include <TH2.h>
#include <TStyle.h>


const string period   = "2012";
const string selection = "el";

const UInt_t ntrig = 38;
string myTriggers[ntrig] = {
  "HLT_Mu22_Photon22_CaloIdL_v",
  "HLT_IsoMu24_v",
  "HLT_IsoMu24_eta2p1_v",
  "HLT_Mu13_Mu8_v",
  "HLT_Mu17_Mu8_v",
  "HLT_Mu17_TkMu8_v",
  "HLT_Mu22_TkMu8_v",
  "HLT_Photon90_CaloIdVL_v",
  "HLT_Photon135_v",
  "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_v",
  "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_v",
  "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_v",
  "HLT_Ele8_CaloIdL_CaloIsoVL_v",
  "HLT_Ele17_CaloIdL_CaloIsoVL_v",
  "HLT_Ele27_WP80_v",
  "HLT_Ele22_CaloIdL_CaloIsoVL_v",
  "HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50",
  "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50",
  "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v",
  "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v",
  "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",
  "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v",
  "HLT_Ele15_Ele8_Ele5_CaloIdL_TrkIdVL",
  "HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_CaloIdT_TrkIdVL ",
  "HLT_TripleEle10_CaloIdL_TrkIdVL",
  "HLT_Ele30_CaloIdVT_TrkIdT_PFJet100_PFJet25",
  "HLT_DoubleEle8_CaloIdT_TrkIdVL_Mass8_PFHT175",
  "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60",
  "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70",
  "HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60",
  "HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60",
  "HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60",
  "HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60",
  "HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50",
  "HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85",
  "HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50",
  "HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50",
  "HLT_Photon36_R9Id85_Photon22_R9Id85"};

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

bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());}

void zgamma::Begin(TTree * tree)
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
  //myRandom = new TRandom3();

  triggerSelector = new TriggerSelector("", period, *triggerNames);

  histoFile = new TFile("a_higgsHistograms.root", "RECREATE");
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
  elIdAndIsoCutsLoose.dPhiIn[0]        = 0.1;
  elIdAndIsoCutsLoose.dEtaIn[0]        = 0.01;
  elIdAndIsoCutsLoose.sigmaIetaIeta[0] = 0.02;
  elIdAndIsoCutsLoose.HadOverEm[0]     = 0.3;
  elIdAndIsoCutsLoose.fabsEPDiff[0]    = 99999;
  elIdAndIsoCutsLoose.dxy[0]           = 0.1;
  elIdAndIsoCutsLoose.dz[0]            = 0.2;
  elIdAndIsoCutsLoose.pfIso04[0]       = 0.3;

  elIdAndIsoCutsLoose.ptErrorOverPt[1] = 99999;
  elIdAndIsoCutsLoose.dPhiIn[1]        = 0.1;
  elIdAndIsoCutsLoose.dEtaIn[1]        = 0.03;
  elIdAndIsoCutsLoose.sigmaIetaIeta[1] = 0.1;
  elIdAndIsoCutsLoose.HadOverEm[1]     = 0.3;
  elIdAndIsoCutsLoose.fabsEPDiff[1]    = 99999;
  elIdAndIsoCutsLoose.dxy[1]           = 0.1;
  elIdAndIsoCutsLoose.dz[1]            = 0.3;
  elIdAndIsoCutsLoose.pfIso04[1]       = 0.3;

  //Tight
  elIdAndIsoCutsTight.ptErrorOverPt[0] = 9999.;
  elIdAndIsoCutsTight.dPhiIn[0]        = 0.06;
  elIdAndIsoCutsTight.dEtaIn[0]        = 0.004;
  elIdAndIsoCutsTight.sigmaIetaIeta[0] = 0.01;
  elIdAndIsoCutsTight.HadOverEm[0]     = 0.12;
  elIdAndIsoCutsTight.fabsEPDiff[0]    = 0.05;
  elIdAndIsoCutsTight.dxy[0]           = 0.02;
  elIdAndIsoCutsTight.dz[0]            = 0.1;
  elIdAndIsoCutsTight.pfIso04[0]       = 0.15;

  elIdAndIsoCutsTight.ptErrorOverPt[1] = 9999.;
  elIdAndIsoCutsTight.dPhiIn[1]        = 0.03;
  elIdAndIsoCutsTight.dEtaIn[1]        = 0.007;
  elIdAndIsoCutsTight.sigmaIetaIeta[1] = 0.03;
  elIdAndIsoCutsTight.HadOverEm[1]     = 0.10;
  elIdAndIsoCutsTight.fabsEPDiff[1]    = 0.05;
  elIdAndIsoCutsTight.dxy[1]           = 0.02;
  elIdAndIsoCutsTight.dz[1]            = 0.1;
  elIdAndIsoCutsTight.pfIso04[1]       = 0.15;

  //muon id cuts
  muIdAndIsoCutsTight.ptErrorOverPt            = 9999.;
  muIdAndIsoCutsTight.NumberOfValidMuonHits    = 0;
  muIdAndIsoCutsTight.NumberOfValidTrackerHits = 10;
  muIdAndIsoCutsTight.NumberOfValidPixelHits   = 0;
  muIdAndIsoCutsTight.NumberOfMatches          = 1;
  muIdAndIsoCutsTight.NormalizedChi2           = 10;
  muIdAndIsoCutsTight.dxy             = 0.02;
  muIdAndIsoCutsTight.dz              = 0.05;
  muIdAndIsoCutsTight.pfIso04         = 0.2;


  phIdAndIsoCutsTight.PassedEleSafeVeto[0] = 1;
  phIdAndIsoCutsTight.HadOverEm[0]         = 0.05;
  phIdAndIsoCutsTight.sigmaIetaIeta[0]     = 0.011;
  phIdAndIsoCutsTight.chIso03[0] = 1.5;
  phIdAndIsoCutsTight.nhIso03[0] = 1.0;
  phIdAndIsoCutsTight.phIso03[0] = 0.7;

  phIdAndIsoCutsTight.PassedEleSafeVeto[1] = 1;
  phIdAndIsoCutsTight.HadOverEm[1]         = 0.05;
  phIdAndIsoCutsTight.sigmaIetaIeta[1]     = 0.033;
  phIdAndIsoCutsTight.chIso03[1] = 1.2;
  phIdAndIsoCutsTight.nhIso03[1] = 1.5;
  phIdAndIsoCutsTight.phIso03[1] = 1.0;


}

void zgamma::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();

}

Bool_t zgamma::Process(Long64_t entry)
{

  GetEntry(entry);
  CountEvents(0);
  if (nEvents[0] % (int)5e4 == 0) cout<<nEvents[3]<<" events passed of "<<nEvents[0]<<" checked! (at Z-peak cut)"<<endl;
  
  //////////////////
  //Trigger status//
  //////////////////
  Bool_t triggerPass  = 1;
  Bool_t isFound  = 0;
  Int_t  prescale = 99;


  TCPrimaryVtx *mainPrimaryVertex = 0, *secondPrimaryVertex = 0;
  mainPrimaryVertex   = (TCPrimaryVtx*)(primaryVtx->At(0));
  //secondPrimaryVertex = (TCPrimaryVtx*)(primaryVtx->At(1));  
  TVector3* pvPosition;// = new TVector3();
  pvPosition = mainPrimaryVertex;


  vector<TCPhysObject> electrons, muons, photons;


  for (Int_t i = 0; i < recoPhotons->GetSize(); ++i) {
    //cout<<"new photon!!!!!!!"<<endl;
    TCPhoton* thisPhoton = (TCPhoton*) recoPhotons->At(i);
    if(thisPhoton->Pt() > 15
       && PassPhotonIdAndIso(thisPhoton, phIdAndIsoCutsTight, pvPosition)
       ) photons.push_back(*thisPhoton);
  }
  
  
  TCPhysObject l1,l2, ph;
  if (photons.size()<1) return kTRUE;
  ph = photons[0];

  for (Int_t i = 0; i <  recoElectrons->GetSize(); ++i) {
    TCElectron* thisElec = (TCElectron*) recoElectrons->At(i);

    if (!(fabs(thisElec->Eta()) < 2.5)) continue;

    if (thisElec->Pt() > 8.0
        && PassElectronIdAndIso(thisElec, elIdAndIsoCutsTight, pvPosition)
        //&& PassElectronIdAndIso(thisElec, elIdAndIsoCutsLoose, pvPosition)
        //&& ph.DeltaR(*thisElec) > 0.4 // this should not fake a photon!
        ) electrons.push_back(*thisElec);
  }
  sort(electrons.begin(), electrons.end(), P4SortCondition);
  
  for (Int_t i = 0; i < recoMuons->GetSize(); ++ i) {
    TCMuon* thisMuon = (TCMuon*) recoMuons->At(i);

    if (!(fabs(thisMuon->Eta()) < 2.4)) continue;
    
    if(thisMuon->Pt() > 7
       && PassMuonIdAndIso(thisMuon, muIdAndIsoCutsTight, pvPosition)
       ) muons.push_back(*thisMuon);    
  }

  sort(muons.begin(), muons.end(), P4SortCondition);
  
  
  
  vector<TCPhysObject> gen_lep;
  TCPhysObject gen_ph;

  for (int i = 0; i < genParticles->GetSize(); ++i) {
    TCGenParticle* thisParticle = (TCGenParticle*) genParticles->At(i);
    if (abs(thisParticle->GetPDGId()) == 11 && thisParticle->Mother()==3000001) {
      //cout<<eventNumber<<"  Found an electron from 3000001"<<endl; 
      gen_lep.push_back(*thisParticle);
    }
  }
  if (gen_lep.size()>2) Abort("Can't have more than two leptons from a decay!!!");
  if (gen_lep.size()>=2)
    hists->fill1DHist(gen_lep[0].DeltaR(gen_lep[1]),"gen_ll_deltaR","gen_ll_deltaR",100,0,5, 1,"");



  if(selection=="el"){ //eegamma
    if (electrons.size()<2) return kTRUE;
    l1 = electrons[0];
    l2 = electrons[1];}

  
  if(selection=="mu"){//mumugamma
    if (muons.size()<2) return kTRUE;
    l1 = muons[0];
    l2 = muons[1];
  }


  CountEvents(1);
  for (UInt_t i =0; i<ntrig; i++)
    {
      triggerSelector->SelectTrigger(myTriggers[i], triggerStatus, hltPrescale, isFound, triggerPass, prescale);
      //if(nEvents[0]==0)
      //cout<<i<<"   "<<myTriggers[i]<<"   is Found = "<<isFound<<"   ispassed = "<<triggerPass<<"  prescale = "<<prescale<<endl;
      if(triggerPass) nEventsTrig[1][i]++;
    }

 

  hists->fill1DHist(l1.Pt(),"l1_pt","Leading lepton_pt", 50,0,200, 1,"");
  hists->fill1DHist(l2.Pt(),"l2_pt","Second Lepton pt", 50,0,200, 1,"");
  hists->fill1DHist(ph.Pt(),"ph_pt","Photon pt", 50,0,200, 1,"");

  hists->fill1DHist((l1+l2).M(),"ll_mass","ll_mass",100,0,100, 1,"");
  hists->fill1DHist(l1.DeltaR(l2),"ll_deltaR","ll_deltaR",100,0,5, 1,"");
  hists->fill1DHist(l1.DeltaR(ph),"l1_ph_deltaR","l1_ph_deltaR",100,0,5, 1,"");
  hists->fill1DHist(l2.DeltaR(ph),"l2_ph_deltaR","l2_ph_deltaR",100,0,5, 1,"");

  if (l1.Pt() < 23 || l2.Pt() < 7) return kTRUE;
  if (ph.Pt() < 23) return kTRUE;
  


  CountEvents(2);
  for (UInt_t i =0; i<ntrig; i++)
    {
      triggerSelector->SelectTrigger(myTriggers[i], triggerStatus, hltPrescale, isFound, triggerPass, prescale);
      if(triggerPass) nEventsTrig[2][i]++;
    }


  if (fabs(ph.Eta()) > 1.44) return kTRUE;
  CountEvents(3);

  for (UInt_t i =0; i<ntrig; i++){
    triggerSelector->SelectTrigger(myTriggers[i], triggerStatus, hltPrescale, isFound, triggerPass, prescale);
    if(triggerPass) nEventsTrig[3][i]++;
  }





  return kTRUE;
  
}

void zgamma::SlaveTerminate() {}

void zgamma::Terminate()
{

  cout<<"| CUT DESCRIPTION                  |\t"<< "\t|"<<endl;
  cout<<"| Initial number of events:        |\t"<< nEvents[0]  <<"\t|"<<"\t|"<<endl;
  cout<<"| number of events:        |\t"<< nEvents[1]  <<"\t|"<<"\t|"<<endl;
  cout<<"| number of events:        |\t"<< nEvents[2]  <<"\t|"<<"\t|"<<endl;
  cout<<"| number of events:        |\t"<< nEvents[3]  <<"\t|"<<"\t|"<<endl;


  cout<<"\n\n | Trigger efficiency 2              |\t"<<endl;
  for(UInt_t n=0; n<ntrig; n++){ 
    UInt_t N=2;
    cout<<n<<"  "<<float(nEventsTrig[N][n])/nEvents[N]<<"   "<<myTriggers[n]<<endl;;
  }

  //cout<<"\n\n | Trigger efficiency 3              |\t"<<endl;
  //for(UInt_t n=0; n<ntrig; n++){ 
  //UInt_t N=3;
  //cout<<n<<"  "<<float(nEventsTrig[N][n])/nEvents[N]<<"   "<<myTriggers[n]<<endl;;
  //}

  histoFile->cd();
  histoFile->Write();
  histoFile->Close();
  cout<<" ** End of Analyzer **"<<endl;
}


void zgamma::CountEvents(Int_t num)
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



float zgamma::CalculateElectronIso(TCElectron *lep)
{
  float eleISO = 0;
  if(period=="2012")
    eleISO = (lep->IsoMap("pfChIso_R04") + TMath::Max(0.0, (Double_t)(lep->IsoMap("pfPhoIso_R04") + lep->IsoMap("pfNeuIso_R04") - rho25Factor*lep->IsoMap("EffArea_R04"))))/lep->Pt();
  else if(period=="2011")
    eleISO = (lep->IsoMap("SumPt_R04") + TMath::Max(0.0, (Double_t)(lep->IsoMap("EmIso_R04") + lep->IsoMap("HadIso_R04") - rhoFactor*TMath::Pi()*0.09)))/lep->Pt();
  else cout<<"Ele iso - Period is not defined!"<<period<<endl;

  return eleISO;
}

 
bool zgamma::PassElectronIdAndIso(TCElectron *lep, elIdAndIsoCuts cuts, TVector3 *pv)
{
  bool pass = false;

  Float_t eleISO = CalculateElectronIso(lep);

  if(((fabs(lep->Eta()) < 1.442
       && lep->PtError()/lep->Pt()       < cuts.ptErrorOverPt[0]
       && lep->SigmaIEtaIEta()           < cuts.sigmaIetaIeta[0]
       && fabs(lep->DphiSuperCluster())  < cuts.dPhiIn[0]
       && fabs(lep->DetaSuperCluster())  < cuts.dEtaIn[0]
       && lep->HadOverEm()               < cuts.HadOverEm[0]
       //&& (period=="2011" || lep->IsoMap("fabsEPDiff")      < cuts.fabsEPDiff[0])
       //&& fabs(lep->Dxy(pv))             < cuts.dxy[0]
       //&& fabs(lep->Dz(pv))              < cuts.dz[0]
       //&& eleISO                         < cuts.pfIso04[0]
       )||
      (fabs(lep->Eta()) >  1.556
       && lep->PtError()/lep->Pt()       < cuts.ptErrorOverPt[1]
       && lep->SigmaIEtaIEta()           < cuts.sigmaIetaIeta[1]
       && fabs(lep->DphiSuperCluster())  < cuts.dPhiIn[1]
       && fabs(lep->DetaSuperCluster())  < cuts.dEtaIn[1]
       && lep->HadOverEm()               < cuts.HadOverEm[1]
       //&& (period=="2011" || lep->IsoMap("fabsEPDiff")      < cuts.fabsEPDiff[1])
       //&& fabs(lep->Dxy(pv))             < cuts.dxy[1]
       //&& fabs(lep->Dz(pv))              < cuts.dz[1]
       //&& eleISO                         < cuts.pfIso04[1]
       ))
     //&& lep->ConversionVeto()

     ) pass = true;

  //pass = true;

  //cout<<"ele pass? "<<pass<<endl;
  return pass;
}


float zgamma::CalculateMuonIso(TCMuon *lep)
{
  Double_t area = 0.09;
  //Double_t area = EffAreaMuon(lep, "2012",0,0);
  float muISO = 0;
  if (period=="2012")
    muISO = (lep->IsoMap("pfChargedHadronPt_R04") + TMath::Max(0.0, lep->IsoMap("pfPhotonEt_R04") + lep->IsoMap("pfNeutralHadronEt_R04")  - TMath::Max(0.0, (Double_t)rho25Factor*area)))/lep->Pt();
  else if(period=="2011")
    muISO = (lep->IsoMap("pfChargedHadronPt_R04") + TMath::Max(0.0, lep->IsoMap("pfPhotonEt_R04") + lep->IsoMap("pfNeutralEt_R04")  - TMath::Max(0.0, (Double_t)rhoFactor*TMath::Pi()*0.09)))/lep->Pt();
  else cout<<"Muon iso - Period is not defined!"<<period<<endl;

  return muISO;
}


bool zgamma::PassMuonIdAndIso(TCMuon *lep, muIdAndIsoCuts cuts, TVector3 *pv)
{
  bool pass = false;
  
  Float_t muISO =  CalculateMuonIso(lep);
  
  if( lep->PtError()/lep->Pt() < cuts.ptErrorOverPt
      && lep->NumberOfValidMuonHits()    > cuts.NumberOfValidMuonHits
      && lep->NumberOfValidTrackerHits() > cuts.NumberOfValidTrackerHits
      && lep->NumberOfValidPixelHits()   > cuts.NumberOfValidPixelHits
      && lep->NumberOfMatches() > cuts.NumberOfMatches
      && lep->NormalizedChi2()  < cuts.NormalizedChi2
      //&& fabs(lep->Dxy(pv))     < cuts.dxy
      //&& fabs(lep->Dz(pv))      < cuts.dz
      //&& muISO                  < cuts.pfIso04
      )
    pass = true;
  //cout<<"muon pass? "<<pass<<endl;
  
  //pass = false;
  //if (muISO < cuts.pfIso04)
  //  pass = true;
  
  return pass;
}

void zgamma::CalculatePhotonIso(TCPhoton *ph, float& chIsoCor, float& nhIsoCor, float& phIsoCor){
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

  chIsoCor = ph->IsoMap("chIso03")-rhoFactor*chEA;
  nhIsoCor = ph->IsoMap("nhIso03")-rhoFactor*nhEA;
  phIsoCor = ph->IsoMap("phIso03")-rhoFactor*phEA;

}

bool zgamma::PassPhotonIdAndIso(TCPhoton *ph, phIdAndIsoCuts cuts, TVector3 *pv)
{
  bool pass = false;
  //Float_t phoISO = 0;
  float tmpEta = ph->SCEta();

  float chIsoCor,nhIsoCor,phIsoCor;
  CalculatePhotonIso(ph, chIsoCor,nhIsoCor,phIsoCor);
  if(
     (fabs(tmpEta)  < 1.442
      && ph->ConversionVeto()       == cuts.PassedEleSafeVeto[0]
      && ph->HadOverEm()             < cuts.HadOverEm[0]
      && ph->SigmaIEtaIEta()         < cuts.sigmaIetaIeta[0]
      && max((double)chIsoCor,0.)    < cuts.chIso03[0]
      && max((double)nhIsoCor,0.)    < cuts.nhIso03[0] + 0.04*ph->Pt()
      && max((double)phIsoCor,0.)    < cuts.phIso03[0] + 0.005*ph->Pt() 
      ) ||
     (fabs(tmpEta)  > 1.566
      && ph->ConversionVeto()       == cuts.PassedEleSafeVeto[1]
      && ph->HadOverEm()             < cuts.HadOverEm[1]
      && ph->SigmaIEtaIEta()         < cuts.sigmaIetaIeta[1]

      && max((double)chIsoCor,0.)    < cuts.chIso03[1]
      && max((double)nhIsoCor,0.)    < cuts.nhIso03[1] + 0.04*ph->Pt() 
      && max((double)phIsoCor,0.)    < cuts.phIso03[1] + 0.005*ph->Pt() 
      )
     ) pass = true;
  
  return pass;
}
