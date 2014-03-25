#define egamma_cxx
#include "egamma.h"

const UInt_t ntrig = 15;
string myTriggers[ntrig] = {
  "HLT_Ele27_WP80_v",

  "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v",
  "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v",
  "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",

  "HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v",
  "HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60_v",
  "HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60_v",
  "HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60_v",
  "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v",
  "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v",

  "HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v",
  "HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v",
  "HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v",
  "HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v",
  "HLT_Photon36_R9Id85_Photon22_R9Id85_v"};

const string pho26OR[6] = {
  "HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v",
  "HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60_v",
  "HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60_v",
  "HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60_v",
  "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v",
  "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v"};

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
Float_t cut_l2pt  = 7, cut_l2pt_low = 7;
Float_t cut_gammapt = 25;

Float_t global_Mll = 0;

Int_t  dal=0, nodal=0;
Bool_t makeMvaTree = 0;
Bool_t makeFitTree = 1;
const UInt_t hisEVTS[] = {154950, 155010};
Int_t evSize = sizeof(hisEVTS)/sizeof(int);

bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());}

void egamma::Begin(TTree * tree)
{
  cout<<"\n***      Begin the Analyzer      ****"<<endl;
  TString option = GetOption();

  cout<<"option = "<<option.Data()<<endl;
  TObjArray *args = (TObjArray*)option.Tokenize(" ");
  //TObjArray *args = (TObjArray*)fOption.Tokenize(" ");
  sample    = (string)((TObjString*)args->At(0))->GetString();
  selection = (string)((TObjString*)args->At(1))->GetString();
  trigger   = (string)((TObjString*)args->At(2))->GetString();
  gen       = (string)((TObjString*)args->At(3))->GetString();
  period    = "2012";

  cout<<sample<<selection<<trigger<<gen<<endl;

  if (gen=="true")
    makeGen=true;
  else
    makeGen=false;

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

  histoFile->mkdir("mvaTree", "mvaTree");
  histoFile->mkdir("fitTree", "fitTree");

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
  muIdAndIsoCutsTight.TrackLayersWithMeasurement = 5;
  muIdAndIsoCutsTight.NumberOfValidMuonHits    = 0;
  muIdAndIsoCutsTight.NumberOfValidTrackerHits = 10;
  muIdAndIsoCutsTight.NumberOfValidPixelHits   = 0;
  muIdAndIsoCutsTight.NumberOfMatches          = 1;
  muIdAndIsoCutsTight.NumberOfMatchedStations  = 1;
  muIdAndIsoCutsTight.NormalizedChi2           = 10;
  muIdAndIsoCutsTight.dxy             = 0.02;
  muIdAndIsoCutsTight.dz              = 0.05;
  muIdAndIsoCutsTight.pfIso04         = 0.4;

  //Soft muon id (specific for dalitz)
  muIdAndIsoCutsSoft.TrackLayersWithMeasurement = 5;
  muIdAndIsoCutsSoft.PixelLayersWithMeasurement = 1;
  muIdAndIsoCutsSoft.NormalizedChi2_tracker = 1.8;
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
  //fout<<"runNumber lumiSection eventNumber"<<endl;

  /*
  histoFile->cd("mvaTree");
  if (makeMvaTree){
    _mvaTree = new TTree("mvaTree", "Supercluster variables");
    _mvaTree->Branch("SCPhiWidth",    &mva_SCPhiWidth,    "SCPhiWidth/F");
    _mvaTree->Branch("SCPhiWidth",    &mva_SCPhiWidth,    "SCPhiWidth/F");
    _mvaTree->Branch("SCEtaWidth",    &mva_SCEtaWidth,    "SCEtaWidth/F");
    _mvaTree->Branch("SigmaIEtaIEta", &mva_SigmaIEtaIEta, "SigmaIEtaIEta/F");
    _mvaTree->Branch("SigmaIPhiIPhi", &mva_SigmaIPhiIPhi, "SigmaIPhiIPhi/F");
    _mvaTree->Branch("HadOverEm",     &mva_HadOverEm,     "HadOverEm/F");
    //_mvaTree->Branch("fabsEPDiff",    &mva_fabsEPDiff,    "fabsEPDiff/F");
    //_mvaTree->Branch("ome1x5oe5x5",   &mva_ome1x5oe5x5,   "ome1x5oe5x5/F");
    _mvaTree->Branch("R9",     &mva_R9,     "R9/F");
    _mvaTree->Branch("fbrem",  &mva_fbrem,  "fbrem/F");
    _mvaTree->Branch("SCdPhi", &mva_SCdPhi, "SCdPhi/F");
    _mvaTree->Branch("SCdEta", &mva_SCdEta, "SCdEta/F");
    _mvaTree->Branch("EoP",    &mva_EoP,    "EoP/F");
    //_mvaTree->Branch("", &mva_, "/F");
    //_mvaTree->Branch("", &mva_, "/F");

  }
  */

  histoFile->cd("fitTree");
  if (makeFitTree){
    _fitTree = new TTree("fitTree", "Mllg tree for fitting");
    _fitTree->Branch("m_llg",  &fit_m_llg,  "m_llg/D");
    _fitTree->Branch("m_ll",   &fit_m_ll,   "m_ll/D");
    _fitTree->Branch("ph_eta", &fit_phEta,  "m_phEta/D");
    _fitTree->Branch("isLowPt",&fit_isLowPt,"isLowPt/O");
    _fitTree->Branch("weight", &fit_weight, "weight/D");
  }
  histoFile->cd();

}

void egamma::SlaveBegin(TTree * /*tree*/)
{   TString option = GetOption(); }

Bool_t egamma::Process(Long64_t entry)
{
  GetEntry(entry);
  CountEvents(0);
  FillHistoCounts(0, 1);
  if (nEvents[0] % (int)5e4 == 0) cout<<nEvents[8]<<" events passed of "<<nEvents[0]<<" checked!"<<endl;

  //return kTRUE;
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


  //----------------------//
  // Selecting triggers   //
  //----------------------//

  Bool_t checkTrigger = kTRUE;
  if (trigger=="pho")
    {
      myTrigger = "HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v";
      cut_l1pt = 20;
      cut_l2pt = 1;
      cut_gammapt = 38;
    }
  else if (trigger=="pho26")
    {
      myTrigger = "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v";
      cut_l1pt = 20;
      cut_l2pt = 1;
      //cut_gammapt = 30;
    }
  else
    checkTrigger = kFALSE;
  //cout<<"Other options are not supported!"<<endl;


  // ----------------------//
  // Gen level particles --//
  //-----------------------//
  vector<TCPhysObject> gen_mu, gen_el;
  TCPhysObject gen_gamma, gen_l1, gen_l2, gen_lPt1, gen_lPt2;
  TCPhysObject gen_higgs, gen_gamma_st1;

  Float_t genMll =0;
  Float_t gendR  =0;


  if(!isRealData && makeGen){
    //Int_t ZID = 3000001; //made-up particle that decays to l+l-
    Int_t A = 25;

    if (sample=="DY")
      A=23;

    Int_t fsr_mu_count = 0, fsr_el_count=0;
    Int_t ph=0, h=0;

    for (int i = 0; i < genParticles->GetSize(); ++i) {
      TCGenParticle* thisParticle = (TCGenParticle*) genParticles->At(i);


      //Higgs himself:
      if (thisParticle->GetPDGId()==25 && thisParticle->GetStatus()==3){
	gen_higgs = *thisParticle;
	hists->fill1DHist(thisParticle->M(), "gen_h_mass",";gen mu mass",  200, 124,126, 1,"GEN");
	h++;
      }


      //GEN ELECTRONS
      if (abs(thisParticle->GetPDGId()) == 11 && thisParticle->GetStatus()==1) {
	if (egamma::GetPrimaryAncestor(thisParticle)->GetPDGId()==A)
	  gen_el.push_back(*thisParticle);
	else if (h==0 && thisParticle->Mother()
		 && egamma::GetPrimaryAncestor(thisParticle)->GetPDGId() == thisParticle->GetPDGId())
	  gen_el.push_back(*thisParticle);
      }


      //PHOTON from the Higgs
      if (sample=="dalitz" && thisParticle->GetPDGId()==22 && thisParticle->GetStatus()==1
	  //&& egamma::GetPrimaryAncestor(thisParticle)->GetPDGId()==A
	  && thisParticle->Mother() &&  abs(thisParticle->Mother()->GetPDGId())!=13 && abs(thisParticle->Mother()->GetPDGId())!=11)
	{
	  gen_gamma = *thisParticle;
	  ph++;
	}
      else if (sample=="DY" &&
	       thisParticle->GetPDGId()==22 && thisParticle->GetStatus()==1 && egamma::GetPrimaryAncestor(thisParticle)->GetPDGId()==23)
	gen_gamma = *thisParticle;



      // **** Find FSR PHOTONS *** //

      if (thisParticle->GetPDGId()==22 //&& thisParticle->GetStatus()==1
	  && thisParticle->Mother()
	  //&& egamma::GetPrimaryAncestor(thisParticle)->GetPDGId()==A
	  && abs(thisParticle->Mother()->GetPDGId())==11
	  )
	{
	  hists->fill1DHist(thisParticle->Pt(), "gen_fsr_pho_pt",";gen fsr photon pt",  200, 0,50, 1,"GEN");
	  // if (egamma::GetPrimaryAncestor(thisParticle)->GetPDGId()!=A)
	  //egamma::DiscoverGeneology(thisParticle, eventNumber);

	  Float_t preFSR_pt = thisParticle->Mother()->Pt();
	  for (int j = 0; j < genParticles->GetSize(); ++j) {
	    TCGenParticle* fsr = (TCGenParticle*) genParticles->At(j);
	    if(fsr->GetPDGId()  == thisParticle->Mother()->GetPDGId()
	       && fsr->Mother() == thisParticle->Mother()){

	      //cout<<"   0--- this one --->>>"<<endl;
	      //egamma::DiscoverGeneology(fsr, eventNumber);
	      //cout<<"   <<---0 this one"<<endl;

	      hists->fill1DHist((preFSR_pt - fsr->Pt())/preFSR_pt, "gen_fsr_el_dPt",";gen el (Pt_preFSR - Pt_afterFSR)/pt",  200, -1,1, 1,"GEN");

	      if(thisParticle->Pt() > 1)
		hists->fill1DHist((preFSR_pt - fsr->Pt())/preFSR_pt, "gen_fsr_el_dPt_photonPt1",";gen el (Pt_preFSR - Pt_afterFSR)/pt",  200, -1,1, 1,"GEN");
	      fsr_mu_count++;
	    }
	  }
	}//end of fsr photons loop
    }// end of genparticles loop


    /*if (ph!=1 && sample=="dalitz")
      Abort(Form("ev #%i NONONO There has to be exactly one photon from the Higgs! \n \t\t but there is %i",
		 (int)eventNumber, (int)ph));
    */
    //if (h!=1 && sample=="dalitz")
    //return kTRUE;
    //Abort(Form(" NONONO There has to be exactly one Higgs! \n \t\t but there is %i", h));

    sort(gen_el.begin(), gen_el.end(), P4SortCondition);

    hists->fill1DHist(fsr_el_count, "gen_el_fsrcount",";gen el FSR",  4, 0,4, 1,"GEN");



    if (gen_el.size()!=2)
      return kTRUE;

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

      //cout<<"\t\t event = "<<eventNumber<<"  "<<endl;
      //cout<<" charge = "<<gen_l1.Charge()<<" pt = "<<gen_l1.Pt()<<" eta="<<gen_l1.Eta()<<" phi="<<gen_l1.Phi()<<endl;
      //cout<<" charge = "<<gen_l2.Charge()<<" pt = "<<gen_l2.Pt()<<" eta="<<gen_l2.Eta()<<" phi="<<gen_l2.Phi()<<endl;


    if (gen_el.size()!=2 && gen_mu.size()!=2)
      Abort(Form("  - - WARNING - -\n ev #%i What's up with that?\n There should be exact two leptons from Higgs. Insteat there %i",
		 (int)eventNumber, (int)(gen_el.size() + gen_mu.size())));

    if (gen_gamma.E()==0 && sample=="dalitz")// return kTRUE;
      Abort(Form("%i: No gen gamma in an event? that sounds bad",(int)eventNumber));//return kTRUE;


    //ZGAnlgles:
    double co1,co2,phi,co3;

    if (sample=="dalitz")
      ang->GetAngles(gen_l1, gen_l2, gen_gamma, co1,co2,phi,co3);
    //cout<<eventNumber<<" gen Angles: c1= "<<co1<<"  c2="<<co2<<"   phi="<<phi<<"   coTh="<<co3<<endl;

    hists->fill1DHist(co1, "gen_co1",";gen cos_lp",  100,-1,1, 1,"");
    hists->fill1DHist(co2, "gen_co2",";gen cos_lm",  100,-1,1, 1,"");
    hists->fill1DHist(co3, "gen_co3",";gen cosTheta",100,-1,1, 1,"");
    hists->fill1DHist(phi, "gen_phi",";gen phi lp",  100, -TMath::Pi(), TMath::Pi(), 1,"");

    FillHistoCounts(1, eventWeight);
    CountEvents(1);

    gendR  = gen_l1.DeltaR(gen_l2);
    genMll = (gen_l1+gen_l2).M();

    //if (genMll>2) return kTRUE;


    hists->fill1DHist(genMll,"gen_Mll_0",";gen_Mll",50,0,15, 1,"eff");
    hists->fill1DHist(gendR, "gen_dR_0", ";gen_dR", 50,0,0.3,1,"eff");
    hists->fillProfile(gendR, genMll, "gen_Mll_vs_dR", ";dR(l1,l2);M(l1,l2)", 100, 0, 0.5, 0, 20, 1,"eff");

    //Acceptance study (do not erase!)
    if (sample =="dalitz" && gen_gamma.Pt()>cut_gammapt && fabs(gen_gamma.Eta())<2.5){
      hists->fill1DHist(genMll, "gen_Mll_acc_gamma",";gen_Mll",50,0,15, 1,"eff");
      hists->fill1DHist(gendR,  "gen_dR_acc_gamma", ";gen_dR", 50,0,0.3,1,"eff");

      if(sample=="dalitz" &&
	 gen_lPt1.Pt()>cut_l1pt && fabs(gen_lPt1.Eta())<2.4 &&
	 gen_lPt2.Pt()>cut_l2pt && fabs(gen_lPt2.Eta())<2.4
	 ){
	hists->fill1DHist(genMll,  "gen_Mll_acc_lept",  ";gen_Mll",50,0,15, 1,"eff");
	hists->fill1DHist(gendR,   "gen_dR_acc_lept",   ";gen_dR", 50,0,0.3,1,"eff");

      }
      else if (sample=="dalitz")
	return kTRUE;
    }
    else if (sample=="dalitz")
      return kTRUE;


    if (sample == "dalitz"){
      hists->fill1DHist(gen_l1.DeltaR(gen_gamma),"gen_l1gamma_deltaR","gen_l1gamma_deltaR",100,0,5, 1,"");
      hists->fill1DHist(gen_l2.DeltaR(gen_gamma),"gen_l2gamma_deltaR","gen_l2gamma_deltaR",100,0,5, 1,"");

      hists->fill1DHist(gen_higgs.Pt(),     "gen_higgs_pt","Pt of higgs",     100, 0,150,  1, "");
      hists->fill1DHist(gen_higgs.M(),      "gen_higgs_mass","Mass of higgs", 100, 50,180,  1, "");

      hists->fill1DHist(gen_gamma.Pt(),     "gen_gamma_pt","Pt of gamma",   50, 0,150,  1, "");
      hists->fill1DHist(gen_gamma.Eta(),    "gen_gamma_eta","Eta of gamma", 50, -3,3,  1, "");
      //hists->fill1DHist(gen_gamma_st1.Pt(), "gen_gamma_st1_pt","Pt of gamma",   50, 0,150,  1, "");
      //hists->fill1DHist(gen_gamma_st1.Eta(),"gen_gamma_st1_eta","Eta of gamma", 50, -3,3,  1, "");
    }

  }

  if (checkTrigger){
    Bool_t pass = false;
    for (Int_t i=0; i<6; i++){
      triggerSelector->SelectTrigger(pho26OR[i], triggerStatus, hltPrescale, isFound, triggerPass, prescale);
      if (triggerPass) {pass = true; break;}
    }
    if (!pass) return kTRUE;
    //  else Abort("Event trigger should mu or el");
  }

  FillHistoCounts(2, eventWeight);
  CountEvents(2);


  /*for(Int_t ev=0; ev<evSize;ev++){
    if (eventNumber==hisEVTS[ev]){cout<<eventNumber<<" Found an event after cut 2"<<endl;
    break;}}*/


  // --- Primary vertex ---//

  if (primaryVtx->GetSize()<1) return kTRUE;
  TCPrimaryVtx *mainPrimaryVertex = 0;//, *secondPrimaryVertex = 0;
  mainPrimaryVertex   = (TCPrimaryVtx*)(primaryVtx->At(0));

  //secondPrimaryVertex = (TCPrimaryVtx*)(primaryVtx->At(1));
  TVector3* pvPosition;// = new TVector3();

  pvPosition = mainPrimaryVertex;

  vector<TCElectron> electrons0, electrons, DAlectrons, DAlectronsTr, fake_electrons0;
  vector<TCMuon> muons0, muons;
  vector<TCPhoton> photons0, photons, photonsHZG, fake_photons;
  TCPhysObject l1,l2, lPt1, lPt2;
  TCPhoton gamma0, gamma, ufoton, ufelectron;


  for (Int_t i = 0; i < recoPhotons->GetSize(); ++i) {
    //cout<<"new photon!!!!!!!"<<endl;
    TCPhoton* thisPhoton = (TCPhoton*) recoPhotons->At(i);
    //cout<<"event = "<<eventNumber
    //  <<"\n pt="<<thisPhoton->Pt()<<" eta="<<thisPhoton->Eta()<<" phi="<<thisPhoton->Phi()<<" px="<<thisPhoton->Px()<<endl;


    if(thisPhoton->Pt() > 20)
      for(Int_t ev=0; ev<evSize;ev++){
	if (eventNumber==hisEVTS[ev]){cout<<"   ---> "<<eventNumber<<" <--- Found an event after cut "<<endl;
	  PhotonDump(*thisPhoton, phIdAndIsoCutsHZG);
	  break;}}


    if (!(fabs(thisPhoton->Eta()) < 2.5)) continue;
    //if (!(fabs(thisPhoton->SCEta()) < 1.4442)) continue;


    if (PassPhotonIdAndIso(thisPhoton, phIdAndIsoCutsHZG, pvPosition)
	&& makeGen && gen_lPt1.DeltaR(*thisPhoton) < 0.1 && thisPhoton->Pt()>20
	&& PassPhotonIdAndIso(thisPhoton, phIdAndIsoCutsHZG, pvPosition)
	)
      fake_photons.push_back(*thisPhoton);

    if(thisPhoton->Pt() > 15){
      if (isRealData || !makeGen || (sample=="dalitz" && makeGen && gen_gamma.DeltaR(*thisPhoton) < 0.1) )
	photons0.push_back(*thisPhoton);

      if(PassPhotonIdAndIso(thisPhoton, phIdAndIsoCutsHZG, pvPosition)
	 //&& (isRealData || !makeGen || (sample=="dalitz" && makeGen && gen_gamma.DeltaR(*thisPhoton) < 0.2) )
	 )
	photonsHZG.push_back(*thisPhoton);

      if(PassPhotonIdAndIso(thisPhoton, phIdAndIsoCutsTight, pvPosition)
	 //&& (isRealData || !makeGen || (sample=="dalitz" && makeGen && gen_gamma.DeltaR(*thisPhoton) < 0.2) )
	 )
	photons.push_back(*thisPhoton);
    }
  }

  sort(photons0.begin(), photons0.end(), P4SortCondition);
  sort(photons.begin(),  photons.end(),  P4SortCondition);
  sort(photonsHZG.begin(),   photonsHZG.end(),   P4SortCondition);
  sort(fake_photons.begin(), fake_photons.end(), P4SortCondition);

  //if (fake_photons.size()>1)
  // Abort("Nah, this is too much.");
  if (fake_photons.size()>=1)
    ufoton = fake_photons[0];


  //if (photons.size()>0)
  //gamma = photons[0];

  if (photonsHZG.size()>0)
    gamma = photonsHZG[0];

  if (makeGen){
    if (photonsHZG.size()>0)
      hists->fill1DHist(photonsHZG[0].DeltaR(gen_gamma),"reco_gen_gamma0_deltaR","reco_gen_gamma0_deltaR",100,0,5, 1,"GEN-RECO");
    if (photonsHZG.size()>1)
      hists->fill1DHist(photonsHZG[1].DeltaR(gen_gamma),"reco_gen_gamma1_deltaR","reco_gen_gamma1_deltaR",100,0,5, 1,"GEN-RECO");
    if (photonsHZG.size()>2)
      hists->fill1DHist(photonsHZG[2].DeltaR(gen_gamma),"reco_gen_gamma2_deltaR","reco_gen_gamma2_deltaR",100,0,5, 1,"GEN-RECO");
  }

  //else if (selection=="el")
  // return kTRUE;
  //for electrons, it does not make sense to continue,
  //because an electron selection  would fail without a photon



  for (Int_t i = 0; i <  recoElectrons->GetSize(); ++i) {
    TCElectron* thisElec = (TCElectron*) recoElectrons->At(i);
    if (!(fabs(thisElec->Eta()) < 2.5)) continue;
    if (thisElec->Pt() > 10.0){

      //if(isRealData || !makeGen || (makeGen && (gen_l1.DeltaR(*thisElec) < 0.1 || gen_l2.DeltaR(*thisElec)<0.1)  && thisElec->Pt()>10))
      //{
      //  if (thisElec->MvaID() > -0.20)
      //}

      if (thisElec->Pt() > 30.0 && fabs(thisElec->Eta()) < 1.4442){
	//electrons0.push_back(*thisElec);
	DAlectrons.push_back(*thisElec);
      }
      if(isRealData || !makeGen || (sample=="dalitz" && makeGen && gen_gamma.DeltaR(*thisElec) < 0.1  && thisElec->Pt()>20))
	fake_electrons0.push_back(*thisElec);

      if (
	  PassElectronIdAndIsoMVA(thisElec)
	  //PassElectronIdAndIso(thisElec, elIdAndIsoCutsTight, pvPosition)
	  //&&(isRealData || !makeGen || ( selection=="el" && makeGen && (gen_l1.DeltaR(*thisElec) < 0.2 ||  gen_l2.DeltaR(*thisElec) < 0.2) )) // this should not fake a photon!
	  //&&(selection!="el" || ( gamma.Pt()>0 && gamma.DeltaR(*thisElec) > 0.4) ) // this should not fake a photon!
	  )
	electrons.push_back(*thisElec);

      // if ( PassElectronIdAndIsoDalitz(thisElec, pvPosition, false)
      //   && gamma.Pt()>0 && gamma.DeltaR(*thisElec) > 0.4  // this should not fake a photon!
	//   )
      //DAlectrons.push_back(*thisElec);

      //if ( PassElectronIdAndIsoDalitz(thisElec, true)
      //   )
      //DAlectronsTr.push_back(*thisElec);

    }
  }

  sort(electrons0.begin(), electrons0.end(), P4SortCondition);
  sort(electrons.begin(),  electrons.end(),  P4SortCondition);
  sort(DAlectrons.begin(), DAlectrons.end(), P4SortCondition);
  //sort(DAlectronsTr.begin(),    DAlectronsTr.end(),    P4SortCondition);
  sort(fake_electrons0.begin(), fake_electrons0.end(), P4SortCondition);

  if (makeGen){
    if (electrons.size()>0)
      hists->fill1DHist(electrons[0].DeltaR(gen_lPt1),"reco_gen_el0_deltaR","reco_gen_el0_deltaR",100,0,5, 1,"GEN-RECO");
    if (electrons.size()>1)
      hists->fill1DHist(electrons[1].DeltaR(gen_lPt1),"reco_gen_el1_deltaR","reco_gen_el1_deltaR",100,0,5, 1,"GEN-RECO");
    if (electrons.size()>2)
      hists->fill1DHist(electrons[2].DeltaR(gen_lPt1),"reco_gen_el2_deltaR","reco_gen_el2_deltaR",100,0,5, 1,"GEN-RECO");

    if (DAlectrons.size()>0)
      hists->fill1DHist(DAlectrons[0].DeltaR(gen_lPt1),"reco_gen_dal0_deltaR","reco_gen_dal0_deltaR",100,0,5, 1,"GEN-RECO");
    if (DAlectrons.size()>1)
      hists->fill1DHist(DAlectrons[1].DeltaR(gen_lPt1),"reco_gen_dal1_deltaR","reco_gen_dal1_deltaR",100,0,5, 1,"GEN-RECO");
    if (DAlectrons.size()>2)
      hists->fill1DHist(DAlectrons[2].DeltaR(gen_lPt1),"reco_gen_dal2_deltaR","reco_gen_dal2_deltaR",100,0,5, 1,"GEN-RECO");

    if (DAlectronsTr.size()>0)
      hists->fill1DHist(DAlectronsTr[0].DeltaR(gen_lPt1),"reco_gen_dalTr0_deltaR","reco_gen_dalTr0_deltaR",100,0,5, 1,"GEN-RECO");
    if (DAlectronsTr.size()>1)
      hists->fill1DHist(DAlectronsTr[1].DeltaR(gen_lPt1),"reco_gen_dalTr1_deltaR","reco_gen_dalTr1_deltaR",100,0,5, 1,"GEN-RECO");
    if (DAlectronsTr.size()>2)
      hists->fill1DHist(DAlectronsTr[2].DeltaR(gen_lPt1),"reco_gen_dalTr2_deltaR","reco_gen_dalTr2_deltaR",100,0,5, 1,"GEN-RECO");


    if (DAlectronsTr.size()>0){
      Float_t dalMll = (DAlectronsTr[0].GetTracks()[0] + DAlectronsTr[0].GetTracks()[1]).M();
      hists->fill1DHist(dalMll,"dal_Mll","dal_Mll",100,0,10, 1,"GEN-RECO");
    }
  }

  if (fake_photons.size()==0 && electrons0.size()==0)
    hists->fill1DHist(0, "egamma_reco",";object; reco fraction",   6, 0,6,  1, "Electrons");
  if (electrons0.size()>=1)
    hists->fill1DHist(1, "egamma_reco",";object; reco fraction",   6, 0,6,  1, "Electrons");
  if (fake_photons.size()>=1)
    hists->fill1DHist(2, "egamma_reco",";object; reco fraction",   6, 0,6,  1, "Electrons");
  if (fake_photons.size()>=1 && electrons0.size()>=1)
    hists->fill1DHist(3, "egamma_reco",";object; reco fraction",   6, 0,6,  1, "Electrons");
  if (fake_photons.size()>=1 || electrons0.size()>=1)
    hists->fill1DHist(4, "egamma_reco",";object; reco fraction",   6, 0,6,  1, "Electrons");
  //else


  for (Int_t i = 0; i < recoMuons->GetSize(); ++ i) {
    TCMuon* thisMuon = (TCMuon*) recoMuons->At(i);
    //cout<<"event = "<<eventNumber
    //  <<"\n pt="<<thisMuon->Pt()<<" eta="<<thisMuon->Eta()<<" phi="<<thisMuon->Phi()<<" px="<<thisMuon->Px()<<endl;

    if (!(fabs(thisMuon->Eta()) < 2.4)) continue;

    if(thisMuon->Pt() > 4
       //&& PassMuonIdAndIso(thisMuon, muIdAndIsoCutsTight, pvPosition)
       && PassMuonIdAndIso(thisMuon, muIdAndIsoCutsSoft, pvPosition)
       ) {

      //cout<<"no mAp "<<thisMuon->PfIsoCharged()<<"  "<< thisMuon->PfIsoNeutral()<<"  "<<thisMuon->PfIsoPhoton()<<endl;
      //cout<<"Mapp   "<<thisMuon->IsoMap("pfChargedHadronPt_R04")<<"  "<<thisMuon->IsoMap("pfNeutralHadronEt_R04")<<"  "<<thisMuon->IsoMap("pfPhotonEt_R04")
      //  <<"  "<<thisMuon->IsoMap("pfPUPt_R04")<<endl;

      muons.push_back(*thisMuon);
    }
  }

  sort(muons.begin(), muons.end(), P4SortCondition);


  hists->fill1DHist(muons.size(),     Form("size_mu_cut%i",  1),";Number of muons",    5,0,5, 1, "Muons");
  hists->fill1DHist(electrons.size(), Form("size_el_cut%i",  1),";Number of electrons",5,0,5, 1, "Electrons");
  hists->fill1DHist(electrons0.size(),Form("size_el0_cut%i", 1),";Number of electrons",5,0,5, 1, "Electrons");
  hists->fill1DHist(photons.size(),   Form("size_ph_cut%i",  1),";Number of photons",  5,0,5, 1, "Photon");
  hists->fill1DHist(photons0.size(),  Form("size_ph0_cut%i", 1),";Number of photons",  5,0,5, 1, "Photon");

  /*
    if (electrons.size()>=1)
      MakeElectronPlots(electrons[0], "NewEle-1");
    if (electrons.size()>=2)
      MakeElectronPlots(electrons[1], "NewEle-2");

  */

  //if (photons0.size()<1) return kTRUE;
  //gamma0 = photons0[0];
  //if (gamma0.Pt() < cut_gammapt) return kTRUE;

  if(!isRealData && makeGen){
    hists->fill1DHist(genMll,  "gen_Mll_reco_gamma",  ";gen_Mll",50,0,15, 1,"eff");
    hists->fill1DHist(gendR,   "gen_dR_reco_gamma",   ";gen_dR", 50,0,0.3,1,"eff");
  }



  for(Int_t ev=0; ev<evSize;ev++){
    if (eventNumber==hisEVTS[ev]){cout<<"   ---> "<<eventNumber<<" <--- Found an event after cut "<<endl;
      for (Int_t s = 0; s<photonsHZG.size();s++)
	{
	  cout<<s<<"  photon"<<endl;
	  PhotonDump(photonsHZG[s], phIdAndIsoCutsHZG);
	}
      break;}}

  //if (photonsHZG.size()<1) return kTRUE;
  //gamma = photonsHZG[0];

  if (gamma.Pt() < cut_gammapt)
    //if (gamma.Pt() < cut_gammapt || fabs(gamma.SCEta())>1.4442)
    return kTRUE;



  /*
    if (photonsHZG.size()==1){
    if (gamma.Pt() < cut_gammapt || fabs(gamma.SCEta())>1.4442)
    return kTRUE;
    }
    else if (photonsHZG.size()==2){
    if ((photonsHZG[0].Pt() < cut_gammapt || fabs(photonsHZG[0].SCEta()) > 1.442) &&
    (photonsHZG[1].Pt() < cut_gammapt || fabs(photonsHZG[1].SCEta()) > 1.442))
    return kTRUE;
    }
  */

  if(!isRealData && makeGen){
    hists->fill1DHist(genMll,  "gen_Mll_reco_gamma_iso",  ";gen_Mll",50,0,15, 1,"eff");
    hists->fill1DHist(gendR,   "gen_dR_reco_gamma_iso",   ";gen_dR", 50,0,0.3,1,"eff");

  }


  FillHistoCounts(3, eventWeight);
  CountEvents(3);


  bool isDalitzEle = false;



  if (electrons0.size()>=1){
    if (sample=="dalitz" && makeMvaTree && gendR < 0.05)
      egamma::FillElecronMVATree(electrons0[0]);

    else if (sample=="DY" && makeMvaTree && gendR > 0.1)
      egamma::FillElecronMVATree(electrons0[0]);
  }


  if(!isRealData && makeGen)
    {
      if (electrons0.size()>=1)
	if (electrons0[0].Pt() > 30){
	  hists->fill1DHist(genMll,  "gen_Mll_one_ele_reco", ";gen_Mll",50,0,15, 1,"eff");
	  hists->fill1DHist(gendR,   "gen_dR_one_ele_reco",  ";gen_dR", 50,0,0.3,1,"eff");
	}
      if (electrons.size()>=1)
	if (electrons[0].Pt() > 30){
	  hists->fill1DHist(genMll,  "gen_Mll_one_ele_reco_ID", ";gen_Mll",50,0,15, 1,"eff");
	  hists->fill1DHist(gendR,   "gen_dR_one_ele_reco_ID",  ";gen_dR", 50,0,0.3,1,"eff");
	}
      if (electrons0.size()>=2)
	if (electrons0[0].Pt() > cut_l1pt && electrons0[1].Pt() > cut_l2pt){
	  hists->fill1DHist(genMll,  "gen_Mll_two_ele_reco", ";gen_Mll",50,0,15, 1,"eff");
	  hists->fill1DHist(gendR,   "gen_dR_two_ele_reco",  ";gen_dR", 50,0,0.3,1,"eff");
	}

      if (electrons.size()>=2)
	if (electrons[0].Pt() > cut_l1pt && electrons[1].Pt() > cut_l2pt){
	  hists->fill1DHist(genMll,  "gen_Mll_two_ele_reco_ID", ";gen_Mll",50,0,15, 1,"eff");
	  hists->fill1DHist(gendR,   "gen_dR_two_ele_reco_ID",  ";gen_dR", 50,0,0.3,1,"eff");
	}
    }

  //if (electrons.size()<2) return kTRUE;

  if (DAlectrons.size()>=1)
    {
      //tracks = DAlectrons[0].GetTracks();
      lPt1 = DAlectrons[0];
      lPt2 = DAlectrons[0];
      l1 = lPt1;
      l2 = lPt2;
      //if (tracks.size()<2) return kTRUE;

      /*
	if (tracks[0].Pt()>tracks[1].Pt()){
	lPt1 = tracks[0];
	lPt2 = tracks[1];
	}
	else {
	lPt1 = tracks[1];
	lPt2 = tracks[0];
	}

	if (tracks[0].Charge()==1 && tracks[1].Charge()==-1){
	l1 = tracks[0];
	l2 = tracks[1];
	}
	else if (tracks[0].Charge()==-1 && tracks[1].Charge()==1){
	l1 = tracks[1];
	l2 = tracks[0];
	}
	else {
	cout<<"Track Charges are the same!  "<<tracks[0].Charge() <<"  "<<tracks[1].Charge()<<endl;
	//return kTRUE;
	}
      */
      isDalitzEle = true;
    }
  else
    return kTRUE;

  if (isDalitzEle)
    dal++;
  else
    nodal++;


  Double_t Mllg = (DAlectrons[0] + gamma).M();

  if (DAlectrons.size()>=1)
    {
      MakeElectronPlots(DAlectrons[0], "DalitzEle");
      if (fabs(DAlectrons[0].Eta()) < 1.44)
	MakeElectronPlots(DAlectrons[0], "DalitzEleEB");
      else
	MakeElectronPlots(DAlectrons[0], "DalitzEleEE");
      //cout<<" EEEEE "<<DAlectrons[0].Pt()<<endl;
    }

  MakePhotonPlots(gamma);

  if (DAlectrons.size()<1) return kTRUE;
  MakeElectronPlots(DAlectrons[0], "DalitzEle-Before");

  if(!isRealData && makeGen){
    hists->fill1DHist(genMll,  "gen_Mll_two_ele_reco_crosscheck", ";gen_Mll",50,0,15, 1,"eff");
    hists->fill1DHist(gendR,   "gen_dR_two_ele_reco_crosscheck",  ";gen_dR", 50,0,0.3,1,"eff");
  }

  FillHistosFull(4, eventWeight, l1, l2, lPt1, lPt2, gamma, "", isDalitzEle);
  FillHistoCounts(4, eventWeight);
  CountEvents(4);

  if (DAlectrons[0].GetTracks().size()<2) return kTRUE;

  l1 = DAlectrons[0].GetTracks()[0];
  l2 = DAlectrons[0].GetTracks()[1];

  FillHistosFull(5, eventWeight, l1, l2, lPt1, lPt2, gamma, "", isDalitzEle);
  FillHistoCounts(5, eventWeight);
  CountEvents(5);


  if (!DAlectrons[0].PassConversionVeto())
    return kTRUE;
  if (DAlectrons[0].ConversionMissHits()!=0)
    return kTRUE;
  //else if (lPt1.DeltaR(gamma)<1.0 || lPt2.DeltaR(gamma)<1.0) return kTRUE;

  //TCTrack::ConversionInfo inf = DAlectrons[0].GetTracks()[0].GetConversionInfo();
 //if (inf.lxyPV > 4 || inf.vtxProb < 1e-6 || inf.nHitsMax > 0 ) return kTRUE;


  MakeElectronPlots(DAlectrons[0], "DalitzEle-AfterAll");
  FillHistosFull(6, eventWeight, l1, l2, lPt1, lPt2, gamma, "", isDalitzEle);
  FillHistoCounts(6, eventWeight);
  CountEvents(6);

  if (fabs(gamma.Eta())<1.444)
    FillHistosFull(6, eventWeight, l1, l2, lPt1, lPt2, gamma, "EB", isDalitzEle);
  else
    FillHistosFull(6, eventWeight, l1, l2, lPt1, lPt2, gamma, "EE", isDalitzEle);

  /*
  for (UInt_t i =0; i<ntrig; i++){
    triggerSelector->SelectTrigger(myTriggers[i], triggerStatus, hltPrescale, isFound, triggerPass, prescale);
    if(triggerPass) nEventsTrig[7][i]++;

    if (!isFound)
      cout<<"TRG **** Warning ***\n The trigger name "<<myTriggers[i]<<" is not in the list of trigger names"<<endl;
  }
  */


  Float_t dxy1 = DAlectrons[0].GetTracks()[0].Dxy(pvPosition);
  Float_t dz1  = DAlectrons[0].GetTracks()[0].Dz(pvPosition);
  if (fabs(dxy1) > 0.02 || fabs(dz1) > 0.1 ) return kTRUE;

  FillHistosFull(7, eventWeight, l1, l2, lPt1, lPt2, gamma, "", isDalitzEle);
  FillHistoCounts(7, eventWeight);
  CountEvents(7);


  if(!isRealData && makeGen && sample=="dalitz"){

    hists->fill1DHist(gen_l2.DeltaR(l2),"reco_gen_l2_deltaR","reco_gen_l2_deltaR",100,0,5, 1,"");

    hists->fill2DHist(gen_l1.DeltaR(l1), gen_l2.DeltaR(l2),      "reco_gen_2D_ll_deltaR",     "reco_gen_2D_ll_deltaR",     100,0,5, 100,0,5, 1,"");
    hists->fill2DHist(gen_l1.DeltaR(l1), gen_gamma.DeltaR(gamma),"reco_gen_2D_l1gamma_deltaR","reco_gen_2D_l1gamma_deltaR",100,0,5, 100,0,5, 1,"");

    hists->fill1DHist(gen_gamma.DeltaR(gamma),"reco_gen_gamma_deltaR","reco_gen_gamma_deltaR",100,0,5, 1,"");
  }


  Double_t Mll  = (l1+l2).M();
  hists->fill1DHist(Mll, Form("diLep_mass_low_cut%i", 7), ";M(ll)", 50, 0,20,  1, "");
  hists->fill1DHist(Mll, Form("diLep_mass_high_cut%i",7), ";M(ll)", 50, 0,120, 1, "");

  if (makeGen){
    hists->fill1DHist(Mll/genMll,            "Mll_responce",";Mll_{reco}/Mll_{gen}",100,0,5, 1,"GEN-RECO");
    hists->fill1DHist(lPt1.Pt()/gen_lPt1.Pt(),"pt1_responce",";pt1_{reco}/p1_{gen}",100,0,5, 1,"GEN-RECO");
    hists->fill1DHist(lPt2.Pt()/gen_lPt2.Pt(),"pt2_responce",";pt2_{reco}/p2_{gen}",100,0,5, 1,"GEN-RECO");
  }


  if (!PassElectronIdAndIsoDalitz(&DAlectrons[0], pvPosition, false))
    return kTRUE;

  FillHistosFull(8, eventWeight, l1, l2, lPt1, lPt2, gamma, "", isDalitzEle);
  FillHistoCounts(8, eventWeight);
  CountEvents(8);

  fout<<eventNumber<<endl;


  if (Mll > 20) return kTRUE;

  FillHistosFull(9, eventWeight, l1, l2, lPt1, lPt2, gamma, "", isDalitzEle);
  FillHistoCounts(9, eventWeight);
  CountEvents(9);

  if (Mllg < 70) return kTRUE;
  FillHistosFull(10, eventWeight, l1, l2, lPt1, lPt2, gamma, "", isDalitzEle);
  FillHistoCounts(10, eventWeight);
  CountEvents(10);

  if (Mllg>120 && Mllg<130){
    FillHistosFull(11, eventWeight, l1, l2, lPt1, lPt2, gamma, "", isDalitzEle);
    FillHistoCounts(11, eventWeight);
    CountEvents(11);

    /*
      if (DAlectrons.size()>=1){
      if (fabs(DAlectrons[0].Eta()) < 1.44)
      MakeElectronPlots(DAlectrons[0], "DalitzEleEB-mH");
      else
      MakeElectronPlots(DAlectrons[0], "DalitzEleEE-mH");
      }
    */
  }

  global_Mll = Mll;

  fit_m_llg   = Mllg;
  fit_m_ll    = Mll;
  fit_weight  = eventWeight;
  fit_phEta   = gamma.SCEta();
  fit_isLowPt = false;

  _fitTree -> Fill();



  if(!isRealData && makeGen){
    hists->fill1DHist(genMll,  "gen_Mll_3",";gen_Mll",50,0,15, 1,"eff");
    hists->fill1DHist(gendR,   "gen_dR_3", ";gen_dR", 50,0,0.3,1,"eff");
  }

  //if ((l1+l2).Pt()+gamma.Pt() < 180  && (l1+l2).Pt()>45 && gamma.Pt()>45){
  if (Mll>2.9 && Mll<3.3){ //jpsi window
    FillHistosFull(12, eventWeight, l1, l2, lPt1, lPt2, gamma, "jpsi", isDalitzEle);
    FillHistoCounts(12, eventWeight);
    CountEvents(12);
    //return kTRUE;
  }

  /*
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
  */
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
  for (UInt_t i =0; i<ntrig; i++) {
      triggerSelector->SelectTrigger(myTriggers[i], triggerStatus, hltPrescale, isFound, triggerPass, prescale);
      if(triggerPass) nEventsTrig[3][i]++;
    }

  */

  return kTRUE;

}

void egamma::SlaveTerminate() {}

void egamma::Terminate()
{

  cout<<"| CUT DESCRIPTION             |\t"<< "\t|"<<endl;
  cout<<"| 0: initial          |\t"<< nEvents[0]  <<"\t|"<<float(nEvents[0])/nEvents[0]<<"\t|"<<endl;
  cout<<"| 1: gen matching     |\t"<< nEvents[1]  <<"\t|"<<float(nEvents[1])/nEvents[0]<<"\t|"<<endl;
  cout<<"| 2: trigger          |\t"<< nEvents[2]  <<"\t|"<<float(nEvents[2])/nEvents[1]<<"\t|"<<endl;
  cout<<"| 3: reco gamma iso   |\t"<< nEvents[3]  <<"\t|"<<float(nEvents[3])/nEvents[2]<<"\t|"<<endl;
  cout<<"| 4: ele >30 ineta<1.44|\t"<< nEvents[4]  <<"\t|"<<float(nEvents[4])/nEvents[3]<<"\t|"<<endl;
  cout<<"| 5: nGSF tracks >=2  |\t"<< nEvents[5]  <<"\t|"<<float(nEvents[5])/nEvents[4]<<"\t|"<<endl;
  cout<<"| 6: Conv and missHits|\t"<< nEvents[6]  <<"\t|"<<float(nEvents[6])/nEvents[5]<<"\t|"<<endl;
  cout<<"| 7: d0 dz cuts  |\t"<< nEvents[7]  <<"\t|"<<float(nEvents[7])/nEvents[6]<<"\t|"<<endl;
  cout<<"| 8: DAL ID WP85 |\t"<< nEvents[8]  <<"\t|"<<float(nEvents[8])/nEvents[7]<<"\t|"<<endl;
  cout<<"| 9:  Mll < 20   |\t"<< nEvents[9]  <<"\t|"<<float(nEvents[9])/nEvents[7]<<"\t|"<<endl;
  cout<<"| 10: m(e,g)>70  |\t"<< nEvents[10] <<"\t|"<<float(nEvents[10])/nEvents[7]<<"\t|"<<endl;
  cout<<"| 11: 120<m(e,g)<130  |\t"<< nEvents[11] <<"\t|"<<float(nEvents[11])/nEvents[7]<<"\t|"<<endl;
  cout<<"| 12: jpsi            |\t"<< nEvents[12] <<"\t|"<<float(nEvents[12])/nEvents[7]<<"\t|"<<endl;

  cout<<"dal = "<<dal<<"   nodal = "<<nodal<<"   tot="<<dal+nodal<<endl;

  /*
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
  */

  cout<<"\n\n | Trigger efficiency 3              |\t"<<endl;
  for(UInt_t n=0; n<ntrig; n++){
    UInt_t N=7;
    cout<<n<<"  "<<float(nEventsTrig[N][n])/nEvents[N]<<"   "<<myTriggers[n]<<endl;;
  }

  hists->fill1DHist(-1, "evt_byCut",";cut #;weighted events", nC+1,-1,nC, totEvents, "Counts");
  hists->fill1DHist(-1, "evt_byCut_raw", ";cut #;events",     nC+1,-1,nC, totEvents, "Counts");


  histoFile->cd();
  histoFile->Write();
  histoFile->Close();
  cout<<" ** End of Analyzer **"<<endl;
}


void egamma::CountEvents(Int_t num)
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



float egamma::CalculateElectronIso(TCElectron *lep)
{
  float eleISO = 0;
  /*
if(period=="2012")
    eleISO = (lep->IsoMap("pfChIso_R04") + TMath::Max(0.0, (Double_t)(lep->IsoMap("pfPhoIso_R04") + lep->IsoMap("pfNeuIso_R04") - rho25Factor*lep->IsoMap("EffArea_R04"))))/lep->Pt();
  else if(period=="2011")
    eleISO = (lep->IsoMap("SumPt_R04") + TMath::Max(0.0, (Double_t)(lep->IsoMap("EmIso_R04") + lep->IsoMap("HadIso_R04") - rhoFactor*TMath::Pi()*0.09)))/lep->Pt();
  else cout<<"Ele iso - Period is not defined!"<<period<<endl;
  */

  return eleISO;
}



bool egamma::PassElectronIdAndIso(TCElectron *lep, elIdAndIsoCuts cuts, TVector3 *pv)
{
  bool pass = false;

  //Float_t eleISO = CalculateElectronIso(lep);

  if(((fabs(lep->Eta()) < 1.4442
       && lep->PtError()/lep->Pt() < cuts.ptErrorOverPt[0]
       && lep->SigmaIEtaIEta()     < cuts.sigmaIetaIeta[0]
       && fabs(lep->SCDeltaPhi())  < cuts.dPhiIn[0]
       && fabs(lep->SCDeltaEta())  < cuts.dEtaIn[0]
       && lep->HadOverEm()         < cuts.HadOverEm[0]
       && fabs(lep->Dxy(pv))       < cuts.dxy[0]
       && fabs(lep->Dz(pv))        < cuts.dz[0]
       //&& eleISO                 < cuts.pfIso04[0]
       )||
      (fabs(lep->Eta()) >  1.556
       && lep->PtError()/lep->Pt() < cuts.ptErrorOverPt[1]
       && lep->SigmaIEtaIEta()     < cuts.sigmaIetaIeta[1]
       && fabs(lep->SCDeltaPhi())  < cuts.dPhiIn[1]
       && fabs(lep->SCDeltaEta())  < cuts.dEtaIn[1]
       && lep->HadOverEm()         < cuts.HadOverEm[1]
       && fabs(lep->Dxy(pv))       < cuts.dxy[1]
       && fabs(lep->Dz(pv))        < cuts.dz[1]
       //&& eleISO                 < cuts.pfIso04[1]
       ))
     //&& lep->ConversionVeto()

     ) pass = true;

  //pass = true;

  //cout<<"ele pass? "<<pass<<endl;
  return pass;
}

bool egamma::PassElectronIdAndIsoDalitz(TCElectron *lep, TVector3 *pv, bool useTracks)
{
  bool pass = false;

  //Float_t eleISO = CalculateElectronIso(lep);

  if(((fabs(lep->Eta()) < 1.442

       //Selection from Ming:
       && lep->SigmaIEtaIEta()    < 0.011
       && lep->E2x5()/lep->E5x5() > 0.914
       && fabs(lep->SCEtaWidth()) < 0.0135
       && fabs(lep->SCDeltaEta()) < 0.0071
       && fabs(lep->SCDeltaPhi()) < 0.0375

       && lep->InverseEnergyMomentumDiff() < 0.04

       //&& fabs(lep->Dxy(pv))       < 0.02
       //&& fabs(lep->Dz(pv))        < 0.1

       //Mine:

       //&& fabs(lep->SCDeltaPhi())  > 0.015
       //&& lep->SCPhiWidth() > 0.02
       //&& lep->R9()  <  0.85
       //&& lep->EoP() > 0.9
       //&& lep->MvaID() > -0.20
       //&& eleISO                         < 0.0
       )
      )
     //(fabs(lep->Eta()) >  1.556
       //&& fabs(lep->SCDeltaPhi())  > 0.015
       //&& lep->SCPhiWidth() > 0.02
       //&& lep->EoP() > 1.0
       //&& lep->R9()  <  0.9 && lep->R9()  >  0.2
       //&& lep->MvaID() > -0.20
       //&& lep->IdMap("kfNLayers") >0
       //&& eleISO                         < cuts.pfIso04[1]


     //&&
     //(!useTracks ||
     //lep->GetTracks().size() >= 2
     //)

     //&& lep->PassConversionVeto()

     ) pass = true;

  //pass = true;

  //cout<<"ele pass? "<<pass<<endl;
  return pass;
}



bool egamma::PassElectronIdAndIsoMVA(TCElectron *lep)
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

float egamma::CalculateMuonIso(TCMuon *lep)
{
  //Double_t area = EffAreaMuon(lep, "2012",0,0);
  float muISO =  (lep->PfIsoCharged() + TMath::Max(0.0, lep->PfIsoPhoton() + lep->PfIsoNeutral()  - 0.5*lep->PfIsoPU()))/lep->Pt();
  return muISO;
}


bool egamma::PassMuonIdAndIso(TCMuon *lep, muIdAndIsoCuts cuts, TVector3 *pv)
{
  //Float_t muISO =  zgamma::CalculateMuonIso(lep);
  Float_t muISO =  CalculateMuonIso(lep);
  bool pass = false;
  if(1
     && lep->IsPF() && lep->IsTRK()
     //&& lep->PtError()/lep->Pt() < cuts.ptErrorOverPt
     //&& lep->NumberOfValidMuonHits()    > cuts.NumberOfValidMuonHits
     //&& lep->NumberOfValidTrackerHits() > cuts.NumberOfValidTrackerHits
     //&& lep->NumberOfValidPixelHits()   > cuts.NumberOfValidPixelHits
     //&& lep->NumberOfMatchedStations() > cuts.NumberOfMatchedStations
     //&& lep->NumberOfMatches() > cuts.NumberOfMatches
     //&& lep->NormalizedChi2()  < cuts.NormalizedChi2
     //&& lep->IsGood()
     //&& lep->TrackLayersWithMeasurement() > cuts.TrackLayersWithMeasurement
     //&& lep->PixelLayersWithMeasurement() > cuts.PixelLayersWithMeasurement
     && lep->NormalizedChi2_tracker()  < cuts.NormalizedChi2_tracker
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

void egamma::CalculatePhotonIso(TCPhoton *ph, float& chIsoCor, float& nhIsoCor, float& phIsoCor){
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
  phIsoCor = ph->PfIsoPhoton() -rhoFactor*phEA;
  //cout<<chIsoCor<<"  "<<nhIsoCor<<"  "<<phIsoCor<<endl;

}

bool egamma::PassPhotonIdAndIso(TCPhoton *ph, phIdAndIsoCuts cuts, TVector3 *pv)
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


void egamma::FillHistoCounts(Int_t num, Double_t weight)
{
  hists->fill1DHist(num, "evt_byCut",";cut;weighted events", nC+1, -1,nC, weight, "Counts");
  hists->fill1DHist(num, "evt_byCut_raw", ";cut;events",     nC+1, -1,nC, 1, "Counts");
}

void egamma::FillHistosFull(Int_t num, Double_t weight,
			    TCPhysObject l1, TCPhysObject l2, TCPhysObject lPt1, TCPhysObject lPt2, TCPhysObject gamma, string dir,
			    bool isMergedEle)
{

  char * d = new char [dir.length()+1];
  std::strcpy (d, dir.c_str());

  TLorentzVector diLep = l1+l2;
  //if(isMergedEle) diLep = l1;
  TLorentzVector tri = diLep + gamma;

  TLorentzVector lPt1Gamma = lPt1+gamma;
  TLorentzVector lPt2Gamma = lPt2+gamma;

  TLorentzVector gammaCM = gamma;
  TLorentzVector diLepCM = diLep;
  TVector3    b1  = -tri.BoostVector();
  gammaCM.Boost(b1);
  diLepCM.Boost(b1);


  hists->fill1DHist(tri.M(),     Form("tri_mass_%s_cut%i",         d, num),";M(ll#gamma)",  100, 0,200,  weight, dir);
  hists->fill1DHist(tri.M(),     Form("tri_mass_longTail_%s_cut%i",d, num),";M(ll#gamma)",  100, 0,400,  weight, dir);

  hists->fill1DHist(lPt1Gamma.M(), Form("lPt1Gamma_%s_cut%i",d, num),";M(l1#gamma)",  100, 0,400,  weight, dir);
  hists->fill1DHist(lPt2Gamma.M(), Form("lPt2Gamma_%s_cut%i",d, num),";M(l2#gamma)",  100, 0,400,  weight, dir);

  hists->fill2DHist(lPt1Gamma.M2(), lPt2Gamma.M2(),
		    Form("h2D_dalitzPlot_%s_cut%i",  d, num),";M(l1#gamma)^{2};M(l2#gamma)^{2}", 100,0,20000, 100, 0,10000,  weight, dir);

  double rotation = atan(-1.);
  double sin_rotation = sin(rotation);
  double cos_rotation = cos(rotation);
  double x_prime = (lPt1Gamma.M2()*cos_rotation - lPt2Gamma.M2()*sin_rotation) + (1.-cos_rotation)*pow(125.,2);
  double y_prime = lPt1Gamma.M2()*sin_rotation + lPt2Gamma.M2()*cos_rotation;

  hists->fill2DHist(x_prime, y_prime,
		    Form("h2D_dalitzPlot_rotation_%s_cut%i",  d, num),";M(l1#gamma)^{2};M(l2#gamma)^{2}", 100,0,22000, 100, -14000,5000,  weight, dir);

  hists->fill2DHist(tri.M(), diLep.M(),
		    Form("h2D_tri_vs_diLep_mass0_%s_cut%i",  d, num),";M(llgamma);M(ll)", 100,0,200, 100, 0,20,  weight, dir);
  hists->fill2DHist(tri.M(), diLep.M(),
		    Form("h2D_tri_vs_diLep_mass1_%s_cut%i",  d, num),";M(llgamma);M(ll)", 100,0,200, 100, 0,1.5, weight, dir);
  hists->fill2DHist(tri.M(), diLep.M(),
		    Form("h2D_tri_vs_diLep_mass2_%s_cut%i",  d, num),";M(llgamma);M(ll)", 100,0,200, 100, 1.5,8, weight, dir);
  hists->fill2DHist(tri.M(), diLep.M(),
		    Form("h2D_tri_vs_diLep_mass3_%s_cut%i",  d, num),";M(llgamma);M(ll)", 100,0,200, 100, 8,20,  weight, dir);

  hists->fill1DHist(diLep.M(),   Form("diLep_mass_low1_%s_cut%i",  d, num),";M(ll)", 100, 0,1.5,  weight, dir);
  hists->fill1DHist(diLep.M(),   Form("diLep_mass_low2_%s_cut%i",  d, num),";M(ll)", 100, 1.5,8,  weight, dir);
  hists->fill1DHist(diLep.M(),   Form("diLep_mass_low3_%s_cut%i",  d, num),";M(ll)", 100, 8,20,   weight, dir);
  hists->fill1DHist(diLep.M(),   Form("diLep_mass_jpsi_%s_cut%i",  d, num),";M(ll)", 100, 2.7,3.5,weight, dir);

  hists->fill1DHist(diLep.M(),   Form("diLep_mass_low_%s_cut%i",  d, num),";M(ll)", 100, 0,20,  weight, dir);
  hists->fill1DHist(diLep.M(),   Form("diLep_mass_high_%s_cut%i", d, num),";M(ll)", 100, 0,120, weight, dir);
  //hists->fill1DHist(diLep.M(),   Form("diLep_mass_jpsi_%s_cut%i",d, dir, num),";M(ll)", 100, 2,5,   weight, "jpsi");
  hists->fill1DHist(diLep.Pt(),  Form("diLep_pt_%s_cut%i",   d, num),";di-Lepton p_{T}", 50, 0,120, weight, dir);
  hists->fill1DHist(diLep.Eta(), Form("diLep_eta_%s_cut%i",  d, num),";di-Lepton eta",   50, -3.5,3.5, weight, dir);
  hists->fill1DHist(diLep.Phi(), Form("diLep_phi_%s_cut%i",  d, num),";di-Lepton phi",   50, -TMath::Pi(),TMath::Pi(), weight, dir);
  hists->fill1DHist(diLepCM.E(), Form("diLep_Ecom_%s_cut%i", d, num),";E(ll) in CoM",    50, 0,200,  weight, dir);

  hists->fill1DHist(gammaCM.E(),Form("gamma_Ecom_%s_cut%i", d, num),";E_{#gamma} in CoM",  50, 0,120,  weight, dir);
  hists->fill1DHist(gamma.E(),  Form("gamma_E_%s_cut%i",    d, num),";Photon Energy", 50, 0,200, weight, dir);
  hists->fill1DHist(gamma.Pt(), Form("gamma_pt_%s_cut%i",   d, num),";Photon p_{T}",  50, 0,100, weight, dir);
  hists->fill1DHist(gamma.Eta(),Form("gamma_eta_%s_cut%i",  d, num),";Photon eta", 50, -3.5,3.5, weight, dir);
  hists->fill1DHist(gamma.Phi(),Form("gamma_phi_%s_cut%i",  d, num),";Photon phi", 50, -TMath::Pi(),TMath::Pi(), weight, dir);
  /*
  hists->fill1DHist(l1.Pt(),    Form("l1_pt_%s_cut%i",  d, num), ";l+ pt",    50, 0,100, weight, dir);
  hists->fill1DHist(l1.Eta(),   Form("l1_eta_%s_cut%i", d, num), ";l+ eta",   50, -3.5,3.5, weight, dir);
  hists->fill1DHist(l1.Phi(),   Form("l1_phi_%s_cut%i", d, num), ";l+ phi",   50, -TMath::Pi(),TMath::Pi(), weight, dir);
  hists->fill1DHist(l2.Pt(),    Form("l2_pt_%s_cut%i",  d, num), ";l- pt",    50, 0,100, weight, dir);
  hists->fill1DHist(l2.Eta(),   Form("l2_eta_%s_cut%i", d, num), ";l- eta",   50, -3.5,3.5, weight, dir);
  hists->fill1DHist(l2.Phi(),   Form("l2_phi_%s_cut%i", d, num), ";l- phi",   50, -TMath::Pi(),TMath::Pi(), weight, dir);
  */
  hists->fill1DHist(lPt1.Pt(),  Form("lPt1_pt_%s_cut%i", d, num), ";Leading lepton p_{T}", 50, 0,100,    weight, dir);
  hists->fill1DHist(lPt1.Eta(), Form("lPt1_eta_%s_cut%i",d, num), ";Leading lepton eta",   50, -3.5,3.5, weight, dir);
  hists->fill1DHist(lPt1.Phi(), Form("lPt1_phi_%s_cut%i",d, num), ";Leading lepton phi",   50, -TMath::Pi(),TMath::Pi(), weight, dir);
  hists->fill1DHist(lPt2.Pt(),  Form("lPt2_pt_%s_cut%i", d, num), ";Trailing lepton p_{T}",50, 0,100,    weight, dir);
  hists->fill1DHist(lPt2.Eta(), Form("lPt2_eta_%s_cut%i",d, num), ";Trailing lepton  eta", 50, -3.5,3.5, weight, dir);
  hists->fill1DHist(lPt2.Phi(), Form("lPt2_phi_%s_cut%i",d, num), ";Trailing lepton  phi", 50, -TMath::Pi(),TMath::Pi(), weight, dir);

  //hists->fill2DHist(l1.Pt(),   l2.Pt(),   Form("h2D_l1_vs_l2_%s_cut%i",  d, num),
  //";l+ pt; l- pt",    50, 0,100, 50,0,100, weight, dir);

  hists->fill2DHist(lPt1.Pt(), lPt2.Pt(), Form("h2D_Pt2_vs_Pt1_%s_cut%i", d, num),
		    ";Leading lepton p_{T};Trailing lepton p_{T}",  50, 0,100, 50,0,100, weight, dir);
  hists->fill2DHist(diLep.Pt(),gamma.Pt(),Form("h2D_diLep_vs_gamma_%s_cut%i",d, num),
		    ";p_{T} of ll system; p_{T} of gamma",          50, 0,100, 50,0,100, weight, dir);

  hists->fill2DHist(gammaCM.E(), gamma.Pt(),Form("h2D_gamma_Ecom_vs_Pt_%s_cut%i",  d, num),
		    ";E_{#gamma} in CoM; p_{T}(#gamma)", 50, 0,100, 50,0,100, weight, dir);
  hists->fill2DHist(gammaCM.E(), tri.M(),   Form("h2D_gamma_Ecom_vs_triM_%s_cut%i",d, num),
		    ";E_{#gamma} in CoM; M(ll#gamma)",   50, 0,100, 50,0,200, weight, dir);

  hists->fill2DHist(diLep.Eta()-gamma.Eta(), gamma.Pt(), Form("h2D_deltaEta_vs_gammaPt_deltaEta_%s_cut%i",d, num),
		    ";#Delta#eta(ll, #gamma);p_{T} of {#gamma}",    100, -5,5, 100,0,130, weight, dir);

  hists->fill1DHist(lPt1.DeltaR(lPt2), Form("ll_deltaR_%s_cut%i",        d, num),";#Delta R(l_{1}, l_{2})",100,0,5, weight,dir);
  hists->fill1DHist(lPt1.DeltaR(gamma),Form("lPt1_gamma_deltaR_%s_cut%i",d, num),";#Delta R(#gamma,l_{1})",100,0,5, weight,dir);
  hists->fill1DHist(lPt2.DeltaR(gamma),Form("lPt2_gamma_deltaR_%s_cut%i",d, num),";#Delta R(#gamma,l_{2})",100,0,5, weight,dir);

  //ZGanlgles:
  if(selection!="el"){

    double co1,co2,phi,co3;
    ang->GetAngles(l1,l2,gamma,co1,co2,phi,co3);
    //cout<<eventNumber<<"RECO  Angles: c1= "<<co1<<"  c2="<<co2<<"   phi="<<phi<<"   coco="<<co3<<endl;

    hists->fill1DHist(co1, Form("co1_cut%i", num), "; cos_lp",  100,-1,1, 1,"Angles");
    hists->fill1DHist(co2, Form("co2_cut%i", num), "; cos_lm,", 100,-1,1, 1,"Angles");
    hists->fill1DHist(co3, Form("co3_cut%i", num), "; cosTheta",100,-1,1, 1,"Angles");
    hists->fill1DHist(phi, Form("phi_cut%i", num), "; phi lp",  100, -TMath::Pi(), TMath::Pi(), 1,"Angles");
  }


}

void egamma::MakeMuonPlots(TCMuon mu, TVector3 *pv)
{
  hists->fill1DHist(mu.PixelLayersWithMeasurement(),"mu_PixelLayersWithMeasurement",";PixelLayersWithMeasurement",   15, 0,15, 1, "Muons");
  hists->fill1DHist(mu.TrackLayersWithMeasurement(),"mu_TrackLayersWithMeasurement",";TrackerLayersWithMeasurement", 30, 0,30, 1, "Muons");
  hists->fill1DHist(mu.NumberOfMatchedStations(),   "mu_NumberOfMatchedStations",   ";NumberOfMatchedStations",  10, 0,10, 1, "Muons");
  hists->fill1DHist(mu.NumberOfValidMuonHits(),     "mu_NumberOfValidMuonHits",     ";NumberOfValidMuonHits",    60, 0,60, 1, "Muons");
  hists->fill1DHist(mu.NumberOfValidTrackerHits(),  "mu_NumberOfValidTrackerHits",  ";NumberOfValidTrackerHits", 40, 0,40, 1, "Muons");
  hists->fill1DHist(mu.NumberOfValidPixelHits(),    "mu_NumberOfValidPixelHits",    ";NumberOfValidPixelHits",   15, 0,15, 1, "Muons");
  hists->fill1DHist(mu.NormalizedChi2(),            "mu_NormalizedChi2",            ";NormalizedChi2",         100, 0,12,  1, "Muons");
  hists->fill1DHist(mu.NormalizedChi2_tracker(),    "mu_NormalizedChi2_tracker",    ";NormalizedChi2_tracker", 100, 0,6,   1, "Muons");
  hists->fill1DHist(mu.Dxy(pv), "mu_dxy", ";dxy", 50, 0,0.1, 1, "Muons");
  hists->fill1DHist(mu.Dz(pv),  "mu_dxz", ";dz",  50, 0,0.3, 1, "Muons");

  Float_t muIso = CalculateMuonIso(&mu);
  hists->fill1DHist(muIso, "mu_iso", ";rel isolation, rho corr", 50, 0,0.7, 1, "Muons");
  Float_t iso2 = (mu.PfIsoCharged() + mu.PfIsoNeutral() + mu.PfIsoPhoton())/mu.Pt();
  hists->fill1DHist(iso2, "mu_iso_official", ";pfIsolationR04", 50, 0,0.7, 1, "Muons");

  hists->fill2DHist(iso2, global_Mll,     "mu_iso_vs_Mll",    ";pfIsolationR04;M(l_{1},l_{2})", 50, 0,0.7, 50, 0,20, 1, "Muons");
  hists->fill1DHist(mu.PtError()/mu.Pt(), "mu_ptErrorOverPt", ";ptErrorOverPt", 50, 0,0.2, 1, "Muons");
}

void egamma::MakePhotonPlots(TCPhoton ph)
{
  hists->fill1DHist(ph.R9(),        "ph_R9",       ";R9",        100,    0, 1,    1,"Photon");
  hists->fill1DHist(ph.E2OverE5(),  "ph_E2OverE5", ";E2OverE5",  100,    0, 1.5,  1,"Photon");
  hists->fill1DHist(ph.TrackVeto(), "ph_trackVeto",";track veto",  3,    0, 3,    1,"Photon");
  hists->fill1DHist(ph.HadOverEm(), "ph_HadOverEm",";HadOverEm", 100,    0, 0.05, 1,"Photon");
  //hists->fill1DHist(ph.NormChi2(),  "ph_NormChi2", ";NormChi2",  100,    0, 0.01, 1,"Photon");
  hists->fill1DHist(ph.SCDeltaPhi(),    "ph_SCDeltaPhi",   ";SCDeltaPhi",    100, -0.2, 0.2,  1,"Photon");
  hists->fill1DHist(ph.SCDeltaEta(),    "ph_SCDeltaEta",   ";SCDeltaEta",    100,-0.02, 0.02, 1,"Photon");
  hists->fill1DHist(ph.GetNCrystals(),  "ph_nCrystals",    ";nCrystals",     160, 0, 160, 1,"Photon");
  hists->fill1DHist(ph.SigmaIEtaIEta(), "ph_SigmaIEtaIEta",";SigmaIEtaIEta", 100, 0, 0.05, 1,"Photon");
  hists->fill1DHist(ph.SigmaIPhiIPhi(), "ph_SigmaIPhiIPhi",";SigmaIPhiIPhi", 100, 0, 0.04, 1,"Photon");
  hists->fill1DHist(ph.ConversionVeto(),"ph_ConversionVeto",";ConversionVeto", 3, 0, 3,   1,"Photon");
  //hists->fill1DHist(ph., "ph_",";", 3, 0, 3, 1, "Photon");
  //hists->fill1DHist(ph., "ph_",";", 3, 0, 3, 1, "Photon");
  //hists->fill1DHist(ph., "ph_",";", 3, 0, 3, 1, "Photon");

}

void egamma::MakeElectronPlots(TCElectron el, string dir)
{
  //cout<<"Making electron plots in dir="<<dir<<endl;


  hists->fill1DHist(el.InverseEnergyMomentumDiff(), dir+"_el_fabsEPDiff",";|1/E - 1/p|", 100, 0, 0.1, 1, dir);

  hists->fill1DHist(el.R9(),    dir+"_el_R9",   ";R9",     100, 0, 1, 1,dir);
  hists->fill1DHist(el.MvaID(), dir+"_el_mvaID",";mva ID", 100, -1,1, 1,dir);
  hists->fill1DHist(el.FBrem(), dir+"_el_fbrem","; fbrem", 100, 0, 1, 1,dir);
  hists->fill1DHist(el.HadOverEm(),  dir+"_el_HadOverEm", ";HadOverEm",  100, 0, 0.1,     1,dir);
  hists->fill1DHist(el.SCDeltaPhi(), dir+"_el_SCDeltaPhi",";SCDeltaPhi", 100, -0.2, 0.2,  1,dir);
  hists->fill1DHist(el.SCDeltaEta(), dir+"_el_SCDeltaEta",";SCDeltaEta", 100, -0.02, 0.02,1,dir);

  hists->fill1DHist(el.PtError()/el.Pt(),     dir+"_el_ptErrorOverPt",";ptErrorOverPt", 100, 0, 1,   1,dir);
  hists->fill1DHist(el.SigmaIEtaIEta(),       dir+"_el_SigmaIEtaIEta",";SigmaIEtaIEta", 100, 0, 0.1, 1,dir);
  hists->fill1DHist(el.SigmaIPhiIPhi(),       dir+"_el_SigmaIPhiIPhi",";SigmaIPhiIPhi", 100, 0, 0.1, 1,dir);
  hists->fill1DHist(el.NormalizedChi2(),      dir+"_el_gsfChi2",      ";gsfChi2",       100, 0, 3,   1,dir);
  hists->fill1DHist(el.NormalizedChi2Kf(),    dir+"_el_kfChi2",       ";kfChi2",        100, 0, 3,   1,dir);
  hists->fill1DHist(el.DeltaEtaSeedCluster(), dir+"_el_dEtaAtCalo",   ";dEtaAtCalo",    100, -0.05, 0.05, 1,dir);
  hists->fill1DHist(el.DeltaPhiSeedCluster(), dir+"_el_dPhiAtCalo",   ";dPhiAtCalo",    100, -0.05, 0.05, 1,dir);
  hists->fill1DHist(el.PreShowerOverRaw(),    dir+"_el_preShowerOverRaw",";preShowerOverRaw", 100, 0, 0.3,1,dir);
  hists->fill1DHist(el.SCEtaWidth(),          dir+"_el_SCEtaWidth",   ";SCEtaWidth",    100, 0, 0.04, 1,dir);
  hists->fill1DHist(el.SCPhiWidth(),          dir+"_el_SCPhiWidth",   ";SCPhiWidth",    100, 0, 0.1,  1,dir);
  hists->fill1DHist(1 - el.E1x5()/el.E5x5(),  dir+"_el_ome1x5oe5x5",  ";1 - e1x5/e5x5", 100, 0, 1,    1,dir);
  hists->fill1DHist(el.EoP(),                 dir+"_el_EoP",    ";EoP",     100, 0, 6,   1,dir);
  hists->fill1DHist(el.EoPout(),              dir+"_el_EoPout", ";EoPout",  100, 0, 6,   1,dir);
  hists->fill1DHist(el.IP3d(),                dir+"_el_ip3d",   ";ip3d",    100, 0, 0.04,1,dir);
  hists->fill1DHist(el.IP3dSig(),             dir+"_el_ip3dSig",";ip3dSig", 100, 0, 3,   1,dir);
  hists->fill1DHist(el.NumberOfValidHits(),   dir+"_el_NumberOfValidHits", ";NumberOfValidHits",   25, 0,25,   1,dir);
  hists->fill1DHist(el.TrackerLayersWithMeasurement(), dir+"_el_TrackerLayersWithMeasurement",
		    ";TrackerLayersWithMeasurement", 20, 0,20,1,dir);

  hists->fill1DHist(el.PassConversionVeto(), dir+"_el_passConv",";pass conv veto", 3, 0, 3, 1, dir);
  hists->fill1DHist(el.ConversionMissHits(), dir+"_el_convMissHits",";conv miss hits", 10, 0, 10, 1, dir);

  hists->fill1DHist(el.ConversionDist(),   dir+"_el_convDist",  ";conv Dist",   200, 0, 50, 1, dir);
  hists->fill1DHist(el.ConversionDcot(),   dir+"_el_convDcot",  ";conv Dcot",   200, -5, 5, 1, dir);
  hists->fill1DHist(el.ConversionRadius(), dir+"_el_convRadius",";conv Radius", 200, 0, 50, 1, dir);

  hists->fill1DHist(el.GetTracks().size(),   dir+"_el_nTracks", ";nTracks",   10, 0,10,   1,dir+"_tracks");

  if (el.GetTracks().size()>=2)

    {
      Float_t mll = (el.GetTracks()[0]+el.GetTracks()[1]).M();
      hists->fill1DHist(mll, dir+"_el_mll_full", ";mll",  200, 0,120, 1,dir);
      hists->fill1DHist(mll, dir+"_el_mll_low",  ";mll",  200, 0,20,  1,dir);
      hists->fill1DHist(mll, dir+"_el_mll_jpsi", ";mll",  200, 1, 6,  1,dir);

    }
  for (Int_t i=0; i<el.GetTracks().size() && i<3; i++){

    string subdir = dir+"_tracks";
    hists->fill1DHist(el.GetTracks()[i].Pt(), subdir+Form("_el_track_%i",i)+"_Pt","Pt", 200, 0,100, 1,subdir);
    hists->fill1DHist(el.GetTracks()[i].PtError()/el.GetTracks()[i].Pt(), subdir+Form("_el_track_%i",i)+"_PtError","Pt Error /Pt", 200, 0,100, 1,subdir);
    hists->fill1DHist(el.GetTracks()[i].NormalizedChi2(), subdir+Form("_el_track_%i",i)+"_NormalizedChi2","NormalizedChi2", 100, 0,3, 1,subdir);

    TCTrack::ConversionInfo inf = el.GetTracks()[i].GetConversionInfo();
    hists->fill1DHist(inf.isValid, subdir+Form("_el_conv_%i",i)+"_isValid","isValid", 3, 0,3, 1,subdir);
    if (inf.isValid){
      hists->fill1DHist(inf.nHitsMax, subdir+Form("_el_conv_%i",i)+"_nHitsMax","nHitsMax", 20, 0,20, 1,subdir);
      hists->fill1DHist(inf.vtxProb,  subdir+Form("_el_conv_%i",i)+"_vtxProb", "vtxProb",  200,0, 1, 1,subdir);
      hists->fill1DHist(inf.lxyPV,    subdir+Form("_el_conv_%i",i)+"_lxyPV",   "lxyPV",    200,0,50, 1,subdir);
      hists->fill1DHist(inf.lxyBS,    subdir+Form("_el_conv_%i",i)+"_lxyBS",   "lxyBS",    200,0,50, 1,subdir);
    }

  }

  Double_t eleISO = (el.PfIsoCharged() + TMath::Max(0.0, (Double_t)(el.PfIsoPhoton() + el.PfIsoNeutral() - rho25Factor*el.EffArea())))/el.Pt();

  hists->fill1DHist(el.PfIsoCharged(), dir+"_el_iso_Charged",";pfChIso_R04",  100, 0, 30, 1, dir);
  hists->fill1DHist(el.PfIsoPhoton(),  dir+"_el_iso_Photon", ";pfPhoIso_R04", 100, 0, 20, 1, dir);
  hists->fill1DHist(el.PfIsoNeutral(), dir+"_el_iso_Neutral",";pfNeuIso_R04", 100, 0, 50, 1, dir);
  hists->fill1DHist(eleISO, dir+"_el_iso",";iso", 100, 0, 1, 1, dir);
  //hists->fill1DHist(el.IdMap(""), dir+"_el_",";", 3, 0, 3, 1, dir);
  //hists->fill1DHist(el., dir+"_el_",";", 3, 0, 3, 1, dir);

}

TCGenParticle * egamma::GetPrimaryAncestor(TCGenParticle *p)
{
  TCGenParticle *a = p;
  while (a->Mother())
    a = a->Mother();
  return a;
}


void egamma::DiscoverGeneology(TCGenParticle *p, ULong64_t ev)
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



void egamma::MuonDump(TCMuon mu, TVector3 *pv)
{
  Float_t muISO = (mu.PfIsoCharged() + TMath::Max(0.0, mu.PfIsoPhoton() + mu.PfIsoNeutral()  - 0.5*mu.IsoMap("pfPUPt_R04")))/mu.Pt();

  cout  << runNumber << " " << eventNumber << " " << mu.Pt()
	<< " " << mu.Eta() << " " << mu.IsGLB() << " " << mu.IsPF() << " "<<mu.IsTRK()<<endl;
  cout  << mu.NormalizedChi2_tracker() << "  " << muISO <<  " " << mu.Dxy(pv) << " " << mu.Dz(pv)
	<< "\n  good =" << mu.IsGood()
	<< "\n  trk =" << mu.TrackLayersWithMeasurement() << "  pix=" <<mu.PixelLayersWithMeasurement()
	<< "\n " << mu.NormalizedChi2() << " " << mu.NumberOfValidMuonHits() << " " << mu.NumberOfMatchedStations()

	<< " " << mu.NumberOfValidPixelHits() << " " << mu.TrackLayersWithMeasurement()
	<< endl;
}


void egamma::PhotonDump(TCPhoton pho, phIdAndIsoCuts cuts)
{
  //Float_t phoISO = 0;
  float chIsoCor,nhIsoCor,phIsoCor;
  float tmpEta = pho.SCEta();

  CalculatePhotonIso(&pho, chIsoCor,nhIsoCor,phIsoCor);

  cout  <<"\n"<< runNumber << " " << eventNumber << " " << pho.Pt() << " " << pho.Eta() << " " << pho.SCEta() <<endl;
  cout<<"conv:"<<pho.ConversionVeto()<<"  hoem:"<< pho.HadOverEm() <<"  sieie:"<<pho.SigmaIEtaIEta()<<endl;
  if (fabs(tmpEta)  < 1.4442)
    cout<<"cut conv:"<<cuts.PassedEleSafeVeto[0]<<"  hoem:"<< cuts.HadOverEm[0]<<"  sieie:"<<cuts.sigmaIetaIeta[0]<<endl;
  else if (fabs(tmpEta)  > 1.5566)
    cout<<"cut conv:"<<cuts.PassedEleSafeVeto[1]<<"  hoem:"<< cuts.HadOverEm[1]<<"  sieie:"<<cuts.sigmaIetaIeta[1]<<endl;
  cout<<"iso"<<endl;
  cout  <<"iso: ch:" << chIsoCor<< "  nh:" << nhIsoCor <<  " ph:"<<phIsoCor<<endl;
  if (fabs(tmpEta)  < 1.4442)
    cout<<"cuts ch="<<cuts.chIso03[0] <<"  nhso="<<cuts.nhIso03[0] + 0.040*pho.Pt()<<"  phiso="<<cuts.phIso03[0] + 0.005*pho.Pt()<<endl;
  else if (fabs(tmpEta)  > 1.5566)
    cout<<"cuts ch="<<cuts.chIso03[1] <<"  "<<cuts.nhIso03[1] + 0.040*pho.Pt()<<"  "<<cuts.phIso03[1] + 0.005*pho.Pt()<<endl;

}

void egamma::FillElecronMVATree(TCElectron el){

  mva_SCPhiWidth    = el.SCPhiWidth();
  mva_SCEtaWidth    = el.SCEtaWidth();
  mva_SigmaIEtaIEta = el.SigmaIEtaIEta();
  mva_SigmaIPhiIPhi = el.SigmaIPhiIPhi();

  mva_SCEta     = el.SCEta();
  mva_R9        = el.R9();
  mva_fbrem     = el.FBrem();
  mva_SCdPhi    = el.SCDeltaPhi();
  mva_SCdEta    = el.SCDeltaEta();
  mva_HadOverEm = el.HadOverEm();

  mva_EoP         = el.EoP();

  _mvaTree->Fill();
}
