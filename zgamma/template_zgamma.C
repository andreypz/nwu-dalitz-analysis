#define zgamma_cxx
#include "zgamma.h"
#include <TH2.h>
#include <TStyle.h>


const string period   = "2012";
const string selection = "@SELECTION";
const string trigger   = "@TRIGGER";

const UInt_t ntrig = 38;
string myTriggers[ntrig] = {
  "HLT_IsoMu24_v",
  "HLT_IsoMu24_eta2p1_v",
  "HLT_Mu13_Mu8_v",
  "HLT_Mu17_Mu8_v",
  "HLT_Mu17_TkMu8_v",
  "HLT_Mu22_TkMu8_v",
  "HLT_Mu22_Photon22_CaloIdL_v",
  "HLT_Photon90_CaloIdVL_v",
  "HLT_Photon135_v",
  "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_v",
  "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_v",
  "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_v",
  "HLT_Ele8_CaloIdL_CaloIsoVL_v",
  "HLT_Ele17_CaloIdL_CaloIsoVL_v",
  "HLT_Ele27_WP80_v",
  "HLT_Ele22_CaloIdL_CaloIsoVL_v",
  "HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v",
  "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v",
  "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v",
  "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v",
  "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",
  "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v",
  "HLT_Ele15_Ele8_Ele5_CaloIdL_TrkIdVL_v",
  "HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_CaloIdT_TrkIdVL_v",
  "HLT_TripleEle10_CaloIdL_TrkIdVL_v",
  "HLT_Ele30_CaloIdVT_TrkIdT_PFJet100_PFJet25_v",
  "HLT_DoubleEle8_CaloIdT_TrkIdVL_Mass8_PFHT175_v",
  "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v",
  "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v",
  "HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v",
  "HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60_v",
  "HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60_v",
  "HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60_v",
  "HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v",
  "HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v",
  "HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v",
  "HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v",
  "HLT_Photon36_R9Id85_Photon22_R9Id85_v"};

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

const bool makeGen = false;
string myTrigger = "";
Float_t cut_l1pt  = 23;
Float_t cut_l2pt  = 9;
Float_t cut_gammapt = 25;

bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());}

void zgamma::Begin(TTree * tree)
{
  cout<<"\n***      Begin the Analyzer      ****"<<endl;
  TString option = GetOption();
  TH1::SetDefaultSumw2(kTRUE);

  vector<string>* triggerNames = 0;
  TFile   *inFile         = tree->GetCurrentFile();
  TTree   *jobTree        = (TTree*)inFile->Get("ntupleProducer/jobTree");
  jobTree->SetBranchAddress("triggerNames", &triggerNames);
  jobTree->GetEntry();

  //cout<<"Printing all the trigger names:"<<endl;
  //for (unsigned i = 0; i < triggerNames->size(); ++i) cout << triggerNames->at(i) << endl;

  //myRandom = new TRandom3();

  triggerSelector = new TriggerSelector("", period, *triggerNames);

  histoFile = new TFile(Form("a_%s_higgsHistograms.root", selection.c_str()), "RECREATE");
  histoFile->mkdir("Counts", "Counts");
  histoFile->mkdir("Muons", "Muons");
  histoFile->mkdir("jpsi", "jpsi");
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

  //Photon Id cuts
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
{   TString option = GetOption(); }

Bool_t zgamma::Process(Long64_t entry)
{
  GetEntry(entry);
  CountEvents(0);
  FillHistoCounts(0, 1);
  if (nEvents[0] % (int)5e4 == 0) cout<<nEvents[3]<<" events passed of "<<nEvents[0]<<" checked!"<<endl;
  
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

  // ----------------------//
  // Gen level particles --//
  //-----------------------//
  vector<TCPhysObject> gen_mu, gen_el;
  TCPhysObject gen_gamma, gen_l1, gen_l2, gen_lPt1, gen_lPt2;
  TCPhysObject gen_higgs, gen_gamma_st1;

  if(!isRealData && makeGen){    
    //Int_t ZID = 3000001; //made-up particle that decays to l+l-
    Int_t ZID = 23;
    for (int i = 0; i < genParticles->GetSize(); ++i) {
      TCGenParticle* thisParticle = (TCGenParticle*) genParticles->At(i);
      //if (thisParticle->Mother()==ZID) {
      if (abs(thisParticle->GetPDGId()) == 11 && thisParticle->GetStatus()==1) 
        gen_el.push_back(*thisParticle);
      if (abs(thisParticle->GetPDGId()) == 13 && (thisParticle->GetStatus()==1 || thisParticle->GetStatus()==2) && (fabs(thisParticle->Grandmother())<7 || thisParticle->Grandmother()==23)) 
        //if (abs(thisParticle->GetPDGId()) == 13 && thisParticle->GetStatus()==1) 
        gen_mu.push_back(*thisParticle);
      //}
      //if (thisParticle->GetPDGId()==22 && thisParticle->GetStatus()==1 && (fabs(thisParticle->Grandmother())<7 || thisParticle->Grandmother()==23) ) // a photon from higgs decay
      if (thisParticle->GetPDGId()==22 && thisParticle->Mother()==25) // a photon from higgs decay
        gen_gamma = *thisParticle;
      if (thisParticle->GetPDGId()==22 && thisParticle->GetStatus()==1) // a photon from higgs decay
        gen_gamma_st1 = *thisParticle;
      if (thisParticle->GetPDGId()==25)
        gen_higgs = *thisParticle;

    }
    sort(gen_el.begin(), gen_el.end(), P4SortCondition);
    sort(gen_mu.begin(), gen_mu.end(), P4SortCondition);
    
    if(selection=="el"){ //eegamma
      if (gen_el.size()<2) return kTRUE;
      //Abort("No gen electrons? What's up with that?");
      gen_lPt1 = gen_el[0];
      gen_lPt2 = gen_el[1];
      
      if (gen_el[0].Charge()==1 && gen_el[1].Charge()==-1){
        gen_l1 = gen_el[0];
        gen_l2 = gen_el[1];
      }
      else if (gen_el[0].Charge()==-1 && gen_el[1].Charge()==1){
        gen_l1 = gen_el[1];
        gen_l2 = gen_el[0];
      }
      else
        Abort("They are the same charge!");
      
      if (gen_el.size()>2) return kTRUE;
      //Abort("Can't have more than two electrons from a decay!!!");
    }
    
    
    for (int i = 0; i < genParticles->GetSize(); ++i) {
      TCGenParticle* thisParticle = (TCGenParticle*) genParticles->At(i);
      //if (abs(thisParticle->GetPDGId()) == 11 && thisParticle->Mother()==3000001) {
      //if (abs(thisParticle->GetPDGId()) == 25)
      //if (abs(thisParticle->GetPDGId()) == 13){
      if (abs(thisParticle->GetPDGId()) == 22 || abs(thisParticle->GetPDGId()) == 13){
        //  (abs(thisParticle->GetPDGId()) == 13 && (thisParticle->GetStatus()==1 || thisParticle->GetStatus()==2))){
        cout<<"event = "<<eventNumber<<"   "<<thisParticle->GetPDGId()
            <<"  st = "<<thisParticle->GetStatus()<<"    mother = "<<thisParticle->Mother()<<"    grandMother = "<<thisParticle->Grandmother()
            <<"\n pt="<<thisParticle->Pt()<<" eta="<<thisParticle->Eta()<<" phi="<<thisParticle->Phi()<<" px="<<thisParticle->Px()<<endl;
        //if ( thisParticle->Grandmother() == thisParticle->GetPDGId())
        //cout<<"\t grandma ?: "<<thisParticle->GetPDGId()<<endl;
        
      }

    }

    if (gen_gamma.E()==0) return kTRUE;
      //Abort(Form("%i: No gen gamma in an event? that sounds bad",eventNumber));//return kTRUE;

    if(selection=="mu"){//mumugamma
      if (gen_mu.size()<2)  return kTRUE;
      
      gen_lPt1 = gen_mu[0];
      gen_lPt2 = gen_mu[1];
      // If there are more muons - it's confusing, so let's just ignore them for now.
      // (they may come from ZH -> ZZgamma process)
      
      
      if (gen_mu[0].Charge()==1 && gen_mu[1].Charge()==-1){
        gen_l1 = gen_mu[0];
        gen_l2 = gen_mu[1];
      }
      else if (gen_mu[0].Charge()==-1 && gen_mu[1].Charge()==1){
        gen_l1 = gen_mu[1];
        gen_l2 = gen_mu[0];
      }
      else
        Abort("They are the same charge!");
      
      if (gen_mu.size()>2) 
        return kTRUE;
      
      //cout<<"\t\t event = "<<eventNumber<<"  "<<endl;
      //cout<<" charge = "<<gen_l1.Charge()<<" pt = "<<gen_l1.Pt()<<" eta="<<gen_l1.Eta()<<" phi="<<gen_l1.Phi()<<endl;
      //cout<<" charge = "<<gen_l2.Charge()<<" pt = "<<gen_l2.Pt()<<" eta="<<gen_l2.Eta()<<" phi="<<gen_l2.Phi()<<endl;
      
      //{
      //}
      //}
      
      //Abort("Can't have more than two muons from a decay!!!");
    }
    
    
    //if(!( gen_lPt1.Pt()>23 && gen_lPt2.Pt()>8 && gen_gamma.Pt()>23)) return kTRUE;
    //if(!( fabs(gen_lPt1.Eta())<2.5 && fabs(gen_lPt2.Eta())<2.5 && fabs(gen_gamma.Eta())<2.5)) return kTRUE;
    
    hists->fill1DHist(gen_l1.DeltaR(gen_l2),   "gen_ll_deltaR","gen_ll_deltaR",100,0,1, 1,"");
    hists->fill1DHist(gen_l1.DeltaR(gen_gamma),"gen_l1gamma_deltaR","gen_l1gamma_deltaR",100,0,5, 1,"");
    hists->fill1DHist(gen_l2.DeltaR(gen_gamma),"gen_l2gamma_deltaR","gen_l2gamma_deltaR",100,0,5, 1,"");

    hists->fill1DHist(gen_higgs.Pt(),     "gen_higgs_pt","Pt of higgs",   100, 0,150,  1, "");
    hists->fill1DHist(gen_higgs.M(),      "gen_higgs_mass","Mass of higgs",   100, 50,180,  1, "");

    hists->fill1DHist(gen_gamma.Pt(),     "gen_gamma_pt","Pt of gamma",   50, 0,150,  1, "");
    hists->fill1DHist(gen_gamma.Eta(),    "gen_gamma_eta","Eta of gamma", 50, -3,3,  1, "");
    hists->fill1DHist(gen_gamma_st1.Pt(), "gen_gamma_st1_pt","Pt of gamma",   50, 0,150,  1, "");
    hists->fill1DHist(gen_gamma_st1.Eta(),"gen_gamma_st1_eta","Eta of gamma", 50, -3,3,  1, "");


  }

  FillHistoCounts(1, eventWeight);
  CountEvents(1);
    
  for (UInt_t i =0; i<ntrig; i++){
    triggerSelector->SelectTrigger(myTriggers[i], triggerStatus, hltPrescale, isFound, triggerPass, prescale);
    if(triggerPass) nEventsTrig[1][i]++;

    if (!isFound)
      cout<<"TRG **** Warning ***\n The trigger name "<<myTriggers[i]<<" is not in the list of trigger names"<<endl;
  }



  //----------------------//
  // Selecting triggers   //
  //----------------------//

  Bool_t checkTrigger = kTRUE;
  if (trigger=="double-mu")
    {
      myTrigger = "HLT_Mu13_Mu8_v";
      //myTrigger = "HLT_Mu17_Mu8_v";
      //cut_l1pt = 15;
      //cut_l2pt = 9;
      //cut_gammapt = 23;
    }
  else if (trigger=="mu-pho")
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
      cut_l2pt = 7;
      cut_gammapt = 23;
    }
  else if (trigger=="pho")
    {
      myTrigger = "HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v";
      cut_l1pt = 38;
      cut_l2pt = 7;
      cut_gammapt = 25;
    }
  else
    checkTrigger = kFALSE;
  //cout<<"Other options are not supported!"<<endl;


  FillHistoCounts(2, eventWeight);
  CountEvents(2);
    

  // --- Primary vertex ---//

  if (primaryVtx->GetSize()<1) return kTRUE;
  TCPrimaryVtx *mainPrimaryVertex = 0;//, *secondPrimaryVertex = 0;
  mainPrimaryVertex   = (TCPrimaryVtx*)(primaryVtx->At(0));

  //secondPrimaryVertex = (TCPrimaryVtx*)(primaryVtx->At(1));  
  TVector3* pvPosition;// = new TVector3();
  pvPosition = mainPrimaryVertex;


  vector<TCPhysObject> electrons, muons, photons;

  for (Int_t i = 0; i < recoPhotons->GetSize(); ++i) {
    //cout<<"new photon!!!!!!!"<<endl;
    TCPhoton* thisPhoton = (TCPhoton*) recoPhotons->At(i);
    //cout<<"event = "<<eventNumber
    //  <<"\n pt="<<thisPhoton->Pt()<<" eta="<<thisPhoton->Eta()<<" phi="<<thisPhoton->Phi()<<" px="<<thisPhoton->Px()<<endl;
    if (!(fabs(thisPhoton->Eta()) < 2.5)) continue;

    if(thisPhoton->Pt() > 15
       && PassPhotonIdAndIso(thisPhoton, phIdAndIsoCutsTight, pvPosition)
       ) photons.push_back(*thisPhoton);
  }
  
  sort(photons.begin(), photons.end(), P4SortCondition);


  for (Int_t i = 0; i <  recoElectrons->GetSize(); ++i) {
    TCElectron* thisElec = (TCElectron*) recoElectrons->At(i);

    if (!(fabs(thisElec->Eta()) < 2.5)) continue;

    if (thisElec->Pt() > 7.0
        //&& PassElectronIdAndIso(thisElec, elIdAndIsoCutsTight, pvPosition)
        && PassElectronIdAndIso(thisElec, elIdAndIsoCutsLoose, pvPosition)
        //&& gamma.DeltaR(*thisElec) > 0.4 // this should not fake a photon!
        ) electrons.push_back(*thisElec);
  }
  sort(electrons.begin(), electrons.end(), P4SortCondition);
  
  for (Int_t i = 0; i < recoMuons->GetSize(); ++ i) {
    TCMuon* thisMuon = (TCMuon*) recoMuons->At(i);
    //cout<<"event = "<<eventNumber
    //  <<"\n pt="<<thisMuon->Pt()<<" eta="<<thisMuon->Eta()<<" phi="<<thisMuon->Phi()<<" px="<<thisMuon->Px()<<endl;

    if (!(fabs(thisMuon->Eta()) < 2.4)) continue;
    
    if(thisMuon->Pt() > 7
       && PassMuonIdAndIso(thisMuon, muIdAndIsoCutsTight, pvPosition)
       ) {
      muons.push_back(*thisMuon);
      MakeMuonPlots(thisMuon, pvPosition);
    }
  }

  sort(muons.begin(), muons.end(), P4SortCondition);
    
  
  TCPhysObject l1,l2, lPt1, lPt2, gamma;
  if(selection=="el"){ //eegamma
    //Only one electron!
    if (electrons.size()<1) return kTRUE;
    lPt1 = electrons[0];
    lPt2 = lPt1;
    
    l1 = electrons[0];
    l2 = l1;
  }

  if(selection=="mu"){//mumugamma
    if (muons.size()<2) return kTRUE;
    lPt1 = muons[0];
    lPt2 = muons[1];
    if (muons[0].Charge()==1 && muons[1].Charge()==-1){
      l1 = muons[0];
      l2 = muons[1];
    }
    else if (muons[0].Charge()==-1 && muons[1].Charge()==1){
      l1 = muons[1];
      l2 = muons[0];
    }
    else{
      cout<<"ch1 = "<<muons[0].Charge()<<"  ch2 = "<<muons[1].Charge()<<endl;
      return kTRUE;//Abort("reco * They are the same charge!");
    }
  }
  Double_t Mll = (l1+l2).M();



  hists->fill1DHist(Mll,  Form("diLep_mass_low_cut%i", 3), ";M(ll)", 50, 0,20,  1, "");
  hists->fill1DHist(Mll,  Form("diLep_mass_high_cut%i", 3),";M(ll)", 50, 0,120, 1, "");


  if (photons.size()<1) return kTRUE;
  gamma = photons[0];

  if(!isRealData && makeGen){
    hists->fill1DHist(gen_l1.DeltaR(l1),"reco_gen_l1_deltaR","reco_gen_l1_deltaR",100,0,5, 1,"");
    hists->fill1DHist(gen_l2.DeltaR(l2),"reco_gen_l2_deltaR","reco_gen_l2_deltaR",100,0,5, 1,"");
    
    hists->fill2DHist(gen_l1.DeltaR(l1), gen_l2.DeltaR(l2),      "reco_gen_2D_ll_deltaR","reco_gen_2D_ll_deltaR",          100,0,5, 100,0,5, 1,"");
    hists->fill2DHist(gen_l1.DeltaR(l1), gen_gamma.DeltaR(gamma),"reco_gen_2D_l1gamma_deltaR","reco_gen_2D_l1gamma_deltaR",100,0,5, 100,0,5, 1,"");
    
    hists->fill1DHist(gen_gamma.DeltaR(gamma),"reco_gen_gamma_deltaR","reco_gen_gamma_deltaR",100,0,5, 1,"");
  }
  


  FillHistosFull(2, eventWeight, l1, l2, lPt1, lPt2, gamma);

  if (lPt1.Pt() < cut_l1pt || lPt2.Pt() < cut_l2pt) return kTRUE;

  FillHistosFull(3, eventWeight, l1, l2, lPt1, lPt2, gamma);
  FillHistoCounts(3, eventWeight);
  CountEvents(3);

  if (Mll > 20) return kTRUE;
  FillHistosFull(4, eventWeight, l1, l2, lPt1, lPt2, gamma);
  FillHistoCounts(4, eventWeight);
  CountEvents(4);

  if (gamma.Pt() < cut_gammapt) return kTRUE;
  FillHistosFull(5, eventWeight, l1, l2, lPt1, lPt2, gamma);
  FillHistoCounts(5, eventWeight);
  CountEvents(5);

  if (lPt1.DeltaR(gamma)<1.0 || lPt2.DeltaR(gamma)<1.0) return kTRUE;

  FillHistosFull(6, eventWeight, l1, l2, lPt1, lPt2, gamma);
  FillHistoCounts(6, eventWeight);
  CountEvents(6);

  if (checkTrigger){
    triggerSelector->SelectTrigger(myTrigger, triggerStatus, hltPrescale, isFound, triggerPass, prescale);
    if (!triggerPass) return kTRUE;
    //  else Abort("Event trigger should mu or el");
  }


  FillHistosFull(7, eventWeight, l1, l2, lPt1, lPt2, gamma);
  FillHistoCounts(7, eventWeight);
  CountEvents(7);



  for (UInt_t i =0; i<ntrig; i++){
    triggerSelector->SelectTrigger(myTriggers[i], triggerStatus, hltPrescale, isFound, triggerPass, prescale);
    if(triggerPass) nEventsTrig[3][i]++;
  }


  if (Mll < 20){
    for (UInt_t i=0; i<ntrig; i++){
      triggerSelector->SelectTrigger(myTriggers[i], triggerStatus, hltPrescale, isFound, triggerPass, prescale);
      if(triggerPass) nEventsTrig[4][i]++;
    }

  }
  

  for (UInt_t i =0; i<ntrig; i++){
    triggerSelector->SelectTrigger(myTriggers[i], triggerStatus, hltPrescale, isFound, triggerPass, prescale);
    if(triggerPass) nEventsTrig[5][i]++;

    Bool_t pa  = 1;
    Bool_t fo  = 0;
    Int_t  pr  = 99;
    triggerSelector->SelectTrigger("HLT_Mu22_Photon22_CaloIdL_v", triggerStatus, hltPrescale, fo, pa, pr);

    if(triggerPass || pa) nEventsTrig[6][i]++;
  }



  if (Mll>2 && Mll<5){ //jpsi window


  }


  /* single ele trigger study
  for (UInt_t i =0; i<ntrig; i++)
    {
      triggerSelector->SelectTrigger(myTriggers[i], triggerStatus, hltPrescale, isFound, triggerPass, prescale);

      Bool_t pa  = 1;
      Bool_t fo  = 0;
      Int_t  pr  = 99;
      triggerSelector->SelectTrigger("HLT_Ele27_WP80_v", triggerStatus, hltPrescale, fo, pa, pr);


      //if(nEvents[0]==0)
      //cout<<i<<"   "<<myTriggers[i]<<"   is Found = "<<isFound<<"   ispassed = "<<triggerPass<<"  prescale = "<<prescale<<endl;
      if(triggerPass) nEventsTrig[2][i]++;
      if(triggerPass || pa) nEventsTrig[5][i]++;
    }
  */


  /*
  CountEvents(3);
  for (UInt_t i =0; i<ntrig; i++) {
      triggerSelector->SelectTrigger(myTriggers[i], triggerStatus, hltPrescale, isFound, triggerPass, prescale);
      if(triggerPass) nEventsTrig[3][i]++;
    }


  if (fabs(gamma.Eta()) > 2.5) return kTRUE;
  CountEvents(4);

  */

  return kTRUE;
  
}

void zgamma::SlaveTerminate() {}

void zgamma::Terminate()
{

  cout<<"| CUT DESCRIPTION             |\t"<< "\t|"<<endl;
  cout<<"| Initial number of events:   |\t"<< nEvents[0]  <<"\t|"<<float(nEvents[0])/nEvents[0]<<"\t|"<<endl;
  cout<<"| Gen particles cuts:         |\t"<< nEvents[1]  <<"\t|"<<float(nEvents[1])/nEvents[0]<<"\t|"<<endl;
  cout<<"| Trigger selection:          |\t"<< nEvents[2]  <<"\t|"<<float(nEvents[2])/nEvents[1]<<"\t|"<<endl;
  cout<<"| pt1,pt2, ptgamma cuts:      |\t"<< nEvents[3]  <<"\t|"<<float(nEvents[3])/nEvents[2]<<"\t|"<<endl;
  cout<<"| M(l,l) <20                  |\t"<< nEvents[4]  <<"\t|"<<float(nEvents[4])/nEvents[3]<<"\t|"<<endl;
  cout<<"| deltaR(lep, gamma) cut      |\t"<< nEvents[5]  <<"\t|"<<float(nEvents[5])/nEvents[4]<<"\t|"<<endl;
  cout<<"| j/psi selection             |\t"<< nEvents[6]  <<"\t|"<<float(nEvents[6])/nEvents[5]<<"\t|"<<endl;


  cout<<"\n\n | Trigger efficiency 2              |\t"<<endl;
  for(UInt_t n=0; n<ntrig; n++){ 
    //UInt_t N=3;
    //cout<<n<<"  "<<float(nEventsTrig[N][n])/nEvents[N]<<"   "<<myTriggers[n]<<endl;;
    cout<<n<<"  "<<float(nEventsTrig[0][n])/nEvents[0]<<"   after sel: "
      //<<nEventsTrig[N][n]<<"/"<<nEvents[N]<<"="
        <<float(nEventsTrig[1][n])/nEvents[1]<<"   "
        <<float(nEventsTrig[2][n])/nEvents[2]<<"   "
        <<float(nEventsTrig[3][n])/nEvents[3]<<"   "
        <<float(nEventsTrig[4][n])/nEvents[4]<<"   "
        <<float(nEventsTrig[5][n])/nEvents[5]<<"   "
        <<float(nEventsTrig[6][n])/nEvents[5]<<"   "
      
        <<myTriggers[n]<<endl;;
  }

  //cout<<"\n\n | Trigger efficiency 3              |\t"<<endl;
  //for(UInt_t n=0; n<ntrig; n++){ 
  //UInt_t N=3;
  //cout<<n<<"  "<<float(nEventsTrig[N][n])/nEvents[N]<<"   "<<myTriggers[n]<<endl;;
  //}

  hists->fill1DHist(-1, "evt_byCut",";cut #;weighted events", nC+1,-1,nC, totEvents, "Counts");
  hists->fill1DHist(-1, "evt_byCut_raw", ";cut #;events",     nC+1,-1,nC, totEvents, "Counts");


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
       //&& lep->IsoMap("fabsEPDiff")      < cuts.fabsEPDiff[0]
       && fabs(lep->Dxy(pv))             < cuts.dxy[0]
       && fabs(lep->Dz(pv))              < cuts.dz[0]
       //&& eleISO                         < cuts.pfIso04[0]
       )||
      (fabs(lep->Eta()) >  1.556
       && lep->PtError()/lep->Pt()       < cuts.ptErrorOverPt[1]
       && lep->SigmaIEtaIEta()           < cuts.sigmaIetaIeta[1]
       && fabs(lep->DphiSuperCluster())  < cuts.dPhiIn[1]
       && fabs(lep->DetaSuperCluster())  < cuts.dEtaIn[1]
       && lep->HadOverEm()               < cuts.HadOverEm[1]
       //&& lep->IsoMap("fabsEPDiff")      < cuts.fabsEPDiff[1]
       && fabs(lep->Dxy(pv))             < cuts.dxy[1]
       && fabs(lep->Dz(pv))              < cuts.dz[1]
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


void zgamma::FillHistoCounts(Int_t num, Double_t weight)
{
  hists->fill1DHist(num, "evt_byCut",";cut #;weighted events", nC+1, -1,nC, weight, "Counts");
  hists->fill1DHist(num, "evt_byCut_raw", ";cut #;events",     nC+1, -1,nC, 1, "Counts");
}

void zgamma::FillHistosFull(Int_t num, Double_t weight,   
                            TCPhysObject l1, TCPhysObject l2, TCPhysObject lPt1, TCPhysObject lPt2, TCPhysObject gamma)
{
  TLorentzVector diLep = l1+l2;
  if(selection=="el") diLep = l1;
  TLorentzVector tri = diLep + gamma;
  
  TLorentzVector gammaCM = gamma;
  TLorentzVector diLepCM = diLep;
  TVector3    b1  = -tri.BoostVector();
  gammaCM.Boost(b1);
  diLepCM.Boost(b1);
  
  
  hists->fill1DHist(tri.M(),     Form("tri_mass_cut%i", num),";M(ll#gamma)",  100, 0,200,  weight, "");
  
  hists->fill1DHist(diLep.M(),   Form("diLep_mass_low_cut%i",  num),";M(ll)", 100, 0,20,  weight, "");
  hists->fill1DHist(diLep.M(),   Form("diLep_mass_high_cut%i", num),";M(ll)", 100, 0,120, weight, "");
  hists->fill1DHist(diLep.M(),   Form("diLep_mass_jpsi_cut%i", num),";M(ll)", 100, 2,5,   weight, "jpsi");
  hists->fill1DHist(diLep.Pt(),  Form("diLep_pt_cut%i",   num),";di-Lepton p_{T}", 50, 0,120, weight, "");
  hists->fill1DHist(diLep.Eta(), Form("diLep_eta_cut%i",  num),";di-Lepton eta",   50, -3.5,3.5, weight, "");
  hists->fill1DHist(diLep.Phi(), Form("diLep_phi_cut%i",  num),";di-Lepton phi",   50, -TMath::Pi(),TMath::Pi(), weight, "");
  hists->fill1DHist(diLepCM.E(), Form("diLep_Ecom_cut%i", num),";E(ll) in CoM",    50, 0,200,  weight, "");
  
  hists->fill1DHist(gammaCM.E(),Form("gamma_Ecom_cut%i", num),";E_{#gamma} in CoM",  50, 0,200,  weight, "");
  hists->fill1DHist(gamma.E(),  Form("gamma_E_cut%i",    num),";Photon Energy", 50, 0,120, weight, "");
  hists->fill1DHist(gamma.Pt(), Form("gamma_pt_cut%i",   num),";Photon p_{T}",  50, 0,100, weight, "");
  hists->fill1DHist(gamma.Eta(),Form("gamma_eta_cut%i",  num),";Photon eta", 50, -3.5,3.5, weight, "");
  hists->fill1DHist(gamma.Phi(),Form("gamma_phi_cut%i",  num),";Photon phi", 50, -TMath::Pi(),TMath::Pi(), weight, "");
  
  hists->fill1DHist(l1.Pt(),    Form("l1_pt_cut%i",  num), ";l+ pt",    50, 0,100, weight, "");
  hists->fill1DHist(l1.Eta(),   Form("l1_eta_cut%i", num), ";l+ eta",   50, -3.5,3.5, weight, "");
  hists->fill1DHist(l1.Phi(),   Form("l1_phi_cut%i", num), ";l+ phi",   50, -TMath::Pi(),TMath::Pi(), weight, "");
  hists->fill1DHist(l2.Pt(),    Form("l2_pt_cut%i",  num), ";l- pt",    50, 0,100, weight, "");
  hists->fill1DHist(l2.Eta(),   Form("l2_eta_cut%i", num), ";l- eta",   50, -3.5,3.5, weight, "");
  hists->fill1DHist(l2.Phi(),   Form("l2_phi_cut%i", num), ";l- phi",   50, -TMath::Pi(),TMath::Pi(), weight, "");
  
  hists->fill1DHist(lPt1.Pt(),  Form("lPt1_pt_cut%i",  num), ";Leading lepton p_{T}", 50, 0,100, weight, "");
  hists->fill1DHist(lPt1.Eta(), Form("lPt1_eta_cut%i", num), ";Leading lepton eta",   50, -3.5,3.5, weight, "");
  hists->fill1DHist(lPt1.Phi(), Form("lPt1_phi_cut%i", num), ";Leading lepton phi",   50, -TMath::Pi(),TMath::Pi(), weight, "");
  hists->fill1DHist(lPt2.Pt(),  Form("lPt2_pt_cut%i",  num), ";Trailing lepton p_{T}",50, 0,100, weight, "");
  hists->fill1DHist(lPt2.Eta(), Form("lPt2_eta_cut%i", num), ";Trailing lepton  eta", 50, -3.5,3.5, weight, "");
  hists->fill1DHist(lPt2.Phi(), Form("lPt2_phi_cut%i", num), ";Trailing lepton  phi", 50, -TMath::Pi(),TMath::Pi(), weight, "");
  
  hists->fill2DHist(lPt1.Pt(), lPt2.Pt(), Form("h2D_Pt2_vs_Pt1_cut%i",     num), ";Leading lepton p_{T};Trailing lepton p_{T}",    50, 0,100, 50,0,100, weight, "");
  hists->fill2DHist(l1.Pt(),   l2.Pt(),   Form("h2D_l1_vs_l2_cut%i",       num), ";l+ pt; l- pt",    50, 0,100, 50,0,100, weight, "");
  hists->fill2DHist(diLep.Pt(),gamma.Pt(),Form("h2D_diLep_vs_gamma_cut%i", num), ";p_{T} of ll system; p_{T} of gamma",    50, 0,100, 50,0,100, weight, "");
  
  hists->fill2DHist(gammaCM.E(), gamma.Pt(),Form("h2D_gamma_Ecom_vs_Pt_cut%i",   num), ";E_{#gamma} in CoM; p_{T}(#gamma)", 50, 0,100, 50,0,100, weight, "");
  hists->fill2DHist(gammaCM.E(), tri.M(),   Form("h2D_gamma_Ecom_vs_triM_cut%i", num), ";E_{#gamma} in CoM; M(ll#gamma)",   50, 0,100, 50,0,200, weight, "");
  
  hists->fill2DHist(diLep.Eta()-gamma.Eta(), gamma.Pt(), Form("h2D_deltaEta_vs_gammaPt_deltaEta_cut%i", num), 
                    ";#Delta#eta(ll, #gamma);p_{T} of {#gamma}",    100, -5,5, 100,0,130, weight, ""); 
  
  hists->fill1DHist(lPt1.DeltaR(lPt2), Form("ll_deltaR_cut%i",       num),";#Delta R(l_{1}, l_{2})",100,0,5, weight,"");
  hists->fill1DHist(lPt1.DeltaR(gamma),Form("l1_gamma_deltaR_cut%i", num),";#Delta R(#gamma,l_{1})",100,0,5, weight,"");
  hists->fill1DHist(lPt2.DeltaR(gamma),Form("l2_gamma_deltaR_cut%i", num),";#Delta R(#gamma,l_{2})",100,0,5, weight,"");
}

void zgamma::MakeMuonPlots(TCMuon *mu, TVector3 *pv)
{
  hists->fill1DHist(mu->NumberOfValidMuonHits(),    "mu_NumberOfValidMuonHits", "NumberOfValidMuonHits", 60, 0,60, 1, "Muons");
  hists->fill1DHist(mu->NumberOfValidTrackerHits(), "mu_NumberOfValidTrackerHits", "NumberOfValidTrackerHits", 40, 0,40, 1, "Muons");
  hists->fill1DHist(mu->NumberOfValidPixelHits(),   "mu_NumberOfValidPixelHits", "NumberOfValidPixelHits", 20, 0,20, 1, "Muons");
  hists->fill1DHist(mu->NormalizedChi2(), "mu_NormalizedChi2", "NormalizedChi2", 50, 0,12, 1, "Muons");
  hists->fill1DHist(mu->Dxy(pv), "mu_dxy", "dxy", 50, 0,0.3, 1, "Muons");
  hists->fill1DHist(mu->Dz(pv),  "mu_dxy", "dz",  50, 0,0.3, 1, "Muons");

  Float_t muIso = CalculateMuonIso(mu);
  hists->fill1DHist(muIso, "mu_iso", "Isolation", 50, 0,0.7, 1, "Muons");

  hists->fill1DHist(mu->PtError()/mu->Pt(), "mu_ptErrorOverPt", "ptErrorOverPt", 50, 0,0.6, 1, "Muons");
}
