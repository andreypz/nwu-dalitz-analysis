#define zgamma_cxx
#include "zgamma.h"

const UInt_t ntrig = 7;
string myTriggers[ntrig] = {
  "HLT_IsoMu24_v",
  "HLT_IsoMu24_eta2p1_v",
  "HLT_Mu13_Mu8_v",
  "HLT_Mu17_Mu8_v",
  "HLT_Mu17_TkMu8_v",
  "HLT_Mu22_TkMu8_v",
  "HLT_Mu22_Photon22_CaloIdL_v"
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
Float_t cut_l2pt  = 7, cut_l2pt_low = 4;
const Float_t cut_iso_mu1 = 0.4, cut_iso_mu2=0;
Float_t cut_gammapt = 25;
Float_t mllMax = 20;

Float_t global_Mll = 0;
Bool_t applyPhosphor = 0;
Int_t  dal=0, nodal=0;
Bool_t makeFitTree = 1;
const UInt_t hisEVTS[] = {9331};
Int_t evSize = sizeof(hisEVTS)/sizeof(int);

bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());}

void zgamma::Begin(TTree * tree)
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
  weighter = new WeightUtils(sample, period, selection, 0);
  roch = new rochcor2012();

  histoFile = new TFile(Form("a_%s_higgsHistograms.root", selection.c_str()), "RECREATE");

  histoFile->mkdir("fitTree", "fitTree");

  hists = new HistManager(histoFile);

  //phoCorrector.reset(new ammagz::PhosphorCorrectionFunctor("../data/PHOSPHOR_NUMBERS_EXPFIT_ERRORS.txt", true));

  for (Int_t n1=0; n1<nC; n1++){
    nEvents[n1]=0;
    for (UInt_t n2=0; n2<ntrig; n2++)
      nEventsTrig[n1][n2]=0;
  }

  //Cuts

  //Tight
  elIdAndIsoCutsTight.ptErrorOverPt[0] = 9999.;
  elIdAndIsoCutsTight.dPhiIn[0] = 0.06;
  elIdAndIsoCutsTight.dEtaIn[0] = 0.004;
  elIdAndIsoCutsTight.sigmaIetaIeta[0] = 0.01;
  elIdAndIsoCutsTight.HadOverEm[0] = 0.12;
  elIdAndIsoCutsTight.fabsEPDiff[0] = 0.05;
  elIdAndIsoCutsTight.dxy[0] = 0.02;
  elIdAndIsoCutsTight.dz[0] = 0.1;
  elIdAndIsoCutsTight.pfIso04[0] = 0.15;

  elIdAndIsoCutsTight.ptErrorOverPt[1] = 9999.;
  elIdAndIsoCutsTight.dPhiIn[1] = 0.03;
  elIdAndIsoCutsTight.dEtaIn[1] = 0.007;
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

  histoFile->cd("fitTree");
  if (makeFitTree){
    _fitTree = new TTree("fitTree", "Mllg tree for fitting");
    _fitTree->Branch("m_llg",  &fit_m_llg, "m_llg/D");
    _fitTree->Branch("m_ll",   &fit_m_ll,  "m_ll/D");
    _fitTree->Branch("ph_eta", &fit_phEta, "m_phEta/D");
    _fitTree->Branch("ph_pt",  &fit_phPt,  "m_phPt/D");
    _fitTree->Branch("di_pt",  &fit_diPt,  "m_diPt/D");
    _fitTree->Branch("isLowPt",&fit_isLowPt,"isLowPt/O");
    _fitTree->Branch("weight", &fit_weight,"weight/D");
  }
  histoFile->cd();

}

void zgamma::SlaveBegin(TTree * /*tree*/)
{   TString option = GetOption(); }

Bool_t zgamma::Process(Long64_t entry)
{
  GetEntry(entry);
  CountEvents(0);
  FillHistoCounts(0, 1);
  if (nEvents[0] % (int)5e4 == 0) cout<<nEvents[8]<<" events passed of "<<nEvents[0]<<" checked!"<<endl;

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
  if (trigger=="double-mu")
    {
      myTrigger = "HLT_Mu13_Mu8_v";
      //myTrigger = "HLT_Mu17_Mu8_v";
      //cut_l1pt = 15;
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
      cut_l2pt = 7;
      cut_gammapt = 23;
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
      
      //GEN MUONS
      //This GEN selection is only valid for a sample with LHE from MCFM
      if (abs(thisParticle->GetPDGId()) == 13 && thisParticle->GetStatus()==1) {
        if (zgamma::GetPrimaryAncestor(thisParticle)->GetPDGId()==A)
          {
            hists->fill1DHist(thisParticle->M(), "gen_mu_mass",";gen mu mass",  200, -1,1, 1,"GEN");
            gen_mu.push_back(*thisParticle);
          }
        //this is for the Madgraph case, where there are events with no Higgs in LHE file
        //:(while there are always should be):
        else if (h==0 && thisParticle->Mother() //good boys should have a mother
                 && zgamma::GetPrimaryAncestor(thisParticle)->GetPDGId() == thisParticle->GetPDGId())
          {
            //cout<<"   --- this one --->>>"<<endl;
            //zgamma::DiscoverGeneology(thisParticle, eventNumber);
            //cout<<"   <<--- this one"<<endl;
            hists->fill1DHist(thisParticle->M(), "gen_mu_mass",";gen mu mass",  200, -1,1, 1,"GEN");
            gen_mu.push_back(*thisParticle);
          }
      }


      //GEN ELECTRONS
      if (abs(thisParticle->GetPDGId()) == 11 && thisParticle->GetStatus()==1) {
        if (zgamma::GetPrimaryAncestor(thisParticle)->GetPDGId()==A)
          gen_el.push_back(*thisParticle);
        else if (h==0 && thisParticle->Mother()
                 && zgamma::GetPrimaryAncestor(thisParticle)->GetPDGId() == thisParticle->GetPDGId())
          gen_el.push_back(*thisParticle);
      }


      //PHOTON from the Higgs
      if (sample=="dalitz" && thisParticle->GetPDGId()==22 && thisParticle->GetStatus()==1
          //&& zgamma::GetPrimaryAncestor(thisParticle)->GetPDGId()==A
          && thisParticle->Mother() &&  abs(thisParticle->Mother()->GetPDGId())!=13 && abs(thisParticle->Mother()->GetPDGId())!=11)
        {
          gen_gamma = *thisParticle;
          ph++;
        }
      else if (sample=="DY" &&
               thisParticle->GetPDGId()==22 && thisParticle->GetStatus()==1 && zgamma::GetPrimaryAncestor(thisParticle)->GetPDGId()==23)
        gen_gamma = *thisParticle;



      // **** Find FSR PHOTONS *** //

      if (thisParticle->GetPDGId()==22 //&& thisParticle->GetStatus()==1
          && thisParticle->Mother()
          //&& zgamma::GetPrimaryAncestor(thisParticle)->GetPDGId()==A
          && selection=="mu" && abs(thisParticle->Mother()->GetPDGId())==13
          )
        {
          hists->fill1DHist(thisParticle->Pt(), "gen_fsr_pho_pt",";gen fsr photon pt",  200, 0,50, 1,"GEN");
          // if (zgamma::GetPrimaryAncestor(thisParticle)->GetPDGId()!=A)
          //zgamma::DiscoverGeneology(thisParticle, eventNumber);
         
          Float_t preFSR_pt = thisParticle->Mother()->Pt();
          for (int j = 0; j < genParticles->GetSize(); ++j) {
            TCGenParticle* fsr = (TCGenParticle*) genParticles->At(j);
            if(fsr->GetPDGId()  == thisParticle->Mother()->GetPDGId()
               && fsr->Mother() == thisParticle->Mother()){
              
              //cout<<"   0--- this one --->>>"<<endl;
              //zgamma::DiscoverGeneology(fsr, eventNumber);
              //cout<<"   <<---0 this one"<<endl;
              
              hists->fill1DHist((preFSR_pt - fsr->Pt())/preFSR_pt, "gen_fsr_mu_dPt",";gen mu (Pt_preFSR - Pt_afterFSR)/pt",  200, -1,1, 1,"GEN");

              if(thisParticle->Pt() > 1)
                hists->fill1DHist((preFSR_pt - fsr->Pt())/preFSR_pt, "gen_fsr_mu_dPt_photonPt1",";gen mu (Pt_preFSR - Pt_afterFSR)/pt",  200, -1,1, 1,"GEN");
              fsr_mu_count++;
            }
          }
        }//end of fsr photons loop
    }// end of genparticles loop

    if (ph!=1 && sample=="dalitz")
      Abort(Form("ev #%i NONONO There has to be exactly one photon from the Higgs! \n \t\t but there is %i", 
                 (int)eventNumber, (int)ph));
    //if (h!=1 && sample=="dalitz")
    //return kTRUE;
    //Abort(Form(" NONONO There has to be exactly one Higgs! \n \t\t but there is %i", h));

    sort(gen_el.begin(), gen_el.end(), P4SortCondition);
    sort(gen_mu.begin(), gen_mu.end(), P4SortCondition);

    hists->fill1DHist(fsr_mu_count, "gen_mu_fsrcount",";gen mu FSR",  4, 0,4, 1,"GEN");
    hists->fill1DHist(fsr_el_count, "gen_el_fsrcount",";gen el FSR",  4, 0,4, 1,"GEN");


    if(selection=="mu"){
      if (gen_mu.size()!=2) return kTRUE;
      //Abort(Form("ev #%i NONONO There has to be exactly 2 muons from the Higgs! \n \t\t but there are %i", eventNumber, (int)gen_mu.size()));
      

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


      //cout<<"\t\t event = "<<eventNumber<<"  "<<endl;
      //cout<<" charge = "<<gen_l1.Charge()<<" pt = "<<gen_l1.Pt()<<" eta="<<gen_l1.Eta()<<" phi="<<gen_l1.Phi()<<endl;
      //cout<<" charge = "<<gen_l2.Charge()<<" pt = "<<gen_l2.Pt()<<" eta="<<gen_l2.Eta()<<" phi="<<gen_l2.Phi()<<endl;

    }


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

    hists->fill1DHist(genMll,"gen_Mll_0",";gen_Mll",100,0,mllMax, 1,"eff");
    hists->fill1DHist(gendR, "gen_dR_0", ";gen_dR", 50,0,0.3,1,"eff");
    hists->fillProfile(gendR, genMll, "gen_Mll_vs_dR", ";dR(l1,l2);M(l1,l2)", 100, 0, 0.5, 0, 20, 1,"eff");

    
    //Acceptance study (do not erase!)
    if (sample =="dalitz" && gen_gamma.Pt()>cut_gammapt && fabs(gen_gamma.Eta())<2.5){
      hists->fill1DHist(genMll, "gen_Mll_acc_gamma",";gen_Mll",100,0,mllMax, 1,"eff");
      hists->fill1DHist(gendR,  "gen_dR_acc_gamma", ";gen_dR", 50,0,0.3,1,"eff");
      
      if(sample=="dalitz" &&
         gen_lPt1.Pt()>cut_l1pt && fabs(gen_lPt1.Eta())<2.4 &&
         gen_lPt2.Pt()>cut_l2pt_low && fabs(gen_lPt2.Eta())<2.4
         ){
        hists->fill1DHist(genMll,  "gen_Mll_acc_lept",  ";gen_Mll",100,0,mllMax, 1,"eff");
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

  FillHistoCounts(2, eventWeight);
  CountEvents(2);


  nVtx = primaryVtx->GetSize();
  // --- Primary vertex ---//

  if (primaryVtx->GetSize()<1) return kTRUE;
  TCPrimaryVtx *mainPrimaryVertex = 0, *secondPrimaryVertex = 0;
  mainPrimaryVertex   = (TCPrimaryVtx*)(primaryVtx->At(0));
  TVector3* pvPosition;// = new TVector3();

  nDofVtx1 = mainPrimaryVertex->NDof();
  if(primaryVtx->GetSize() >1){
    secondPrimaryVertex = (TCPrimaryVtx*)(primaryVtx->At(1));
    nDofVtx2 = secondPrimaryVertex->NDof();
  }
  pvPosition = mainPrimaryVertex;

  vector<TCElectron> electrons0, electrons, electrons_dalitz, fake_electrons0;
  vector<TCMuon> muons0, muons;
  vector<TCPhoton> photons0, photonsTight, photonsHZG, fake_photons;
  TCPhysObject l1,l2, lPt1, lPt2;
  TCPhoton gamma0, gamma, ufoton, ufelectron;


  for (Int_t i = 0; i < recoPhotons->GetSize(); ++i) {
    //cout<<"new photon!!!!!!!"<<endl;
    TCPhoton* thisPhoton = (TCPhoton*) recoPhotons->At(i);
    //cout<<"in photon loop event = "<<eventNumber
    //  <<"\n pt="<<thisPhoton->Pt()<<" eta="<<thisPhoton->Eta()<<" phi="<<thisPhoton->Phi()<<" px="<<thisPhoton->Px()<<endl;

    if (isRealData || !makeGen || (makeGen &&gen_lPt1.DeltaR(*thisPhoton) < 0.1 && thisPhoton->Pt()>20))
      fake_photons.push_back(*thisPhoton);

    if (!(fabs(thisPhoton->Eta()) < 2.5)) continue;
    //if (!(fabs(thisPhoton->SCEta()) < 2.5)) continue;


    if(thisPhoton->Pt() > 15){
      if (isRealData || !makeGen || (sample=="dalitz" && makeGen && gen_gamma.DeltaR(*thisPhoton) < 0.1) )
        photons0.push_back(*thisPhoton);


      if(PassPhotonIdAndIso(thisPhoton, phIdAndIsoCutsHZG, pvPosition)
         //&& (isRealData || !makeGen || (sample=="dalitz" && makeGen && gen_gamma.DeltaR(*thisPhoton) < 0.2) )
         )
        {    

          if (applyPhosphor)
            {
              Float_t corrPhoPt = thisPhoton->Pt();
              
              if (isRealData){
                //corrPhoPt = phoCorrector->GetCorrEtData(thisPhoton->R9(), 2012, thisPhoton->Pt(), thisPhoton->Eta());
                thisPhoton->SetPtEtaPhiM(corrPhoPt,thisPhoton->Eta(),thisPhoton->Phi(),0.0);
              }
              else
                {
                  TCGenParticle goodGenPhoton;
                  Float_t testDr = 9999;

                  //PhotonR9Corrector(*thisPhoton);
                  
                  for (int j = 0; j < genParticles->GetSize(); ++j) {
                    TCGenParticle* thisParticle = (TCGenParticle*) genParticles->At(j);
                    if (thisParticle->GetPDGId()==22 && thisParticle->Mother() &&  thisParticle->GetStatus()==1
                        && (thisParticle->Mother()->GetPDGId()==25 || thisParticle->Mother()->GetPDGId()==22)){
                      
                      if(thisPhoton->DeltaR(*thisParticle)<testDr){
                        goodGenPhoton = *thisParticle;
                        testDr = thisPhoton->DeltaR(*thisParticle);
                      }
                    }
                  }
                  //if (testDr < 0.2)
                  //corrPhoPt = phoCorrector->GetCorrEtMC(thisPhoton->R9(), 2012, thisPhoton->Pt(), thisPhoton->Eta(), goodGenPhoton.E());
                  thisPhoton->SetPtEtaPhiM(corrPhoPt,thisPhoton->Eta(),thisPhoton->Phi(),0.0);
                }
            }
          photonsHZG.push_back(*thisPhoton);
        }
      
      if(PassPhotonIdAndIso(thisPhoton, phIdAndIsoCutsTight, pvPosition)
         //&& (isRealData || !makeGen || (sample=="dalitz" && makeGen && gen_gamma.DeltaR(*thisPhoton) < 0.2) )
         )
        photonsTight.push_back(*thisPhoton);
    }
  }


  sort(photons0.begin(),     photons0.end(),     P4SortCondition);
  sort(photonsTight.begin(), photonsTight.end(), P4SortCondition);
  sort(photonsHZG.begin(),   photonsHZG.end(),   P4SortCondition);
  sort(fake_photons.begin(), fake_photons.end(), P4SortCondition);

  //if (fake_photons.size()>1)
  // Abort("Nah, this is too much.");
  if (fake_photons.size()>=1)
    ufoton = fake_photons[0];


  //if (photonsTight.size()<1) return kTRUE; 
  //gamma = photonsTight[0];

  if (photonsHZG.size()<1) return kTRUE;
  gamma = photonsHZG[0];

  if (!isRealData)
    eventWeight *= weighter->PhotonSF(gamma);

  //if (photons0.size()<1) return kTRUE;
  //gamma0 = photons0[0];
  //if (gamma0.Pt() < cut_gammapt) return kTRUE;
 
  for (Int_t i = 0; i <  recoElectrons->GetSize(); ++i) {
    TCElectron* thisElec = (TCElectron*) recoElectrons->At(i);
    if (!(fabs(thisElec->Eta()) < 2.5)) continue;
    if (thisElec->Pt() > 10.0){
      if(isRealData || !makeGen || ( selection=="el" && makeGen && (gen_l1.DeltaR(*thisElec) < 0.1 || gen_l2.DeltaR(*thisElec)<0.1)  && thisElec->Pt()>10))
        {
          if (thisElec->MvaID() > -0.20)
            electrons0.push_back(*thisElec);
        }
      if (
          //PassElectronIdAndIsoMVA(thisElec)
          PassElectronIdAndIso(thisElec, elIdAndIsoCutsTight, pvPosition)
          )
        electrons.push_back(*thisElec);


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

     
      Float_t qter = 1;
      //cout<<"DBG before rochcor"<<endl;
      if (isRealData)
        roch->momcor_data(*thisMuon, thisMuon->Charge(), 0, qter); 
      else
        roch->momcor_mc(*thisMuon,  thisMuon->Charge(), 0, qter ); 

      //cout<<"DBG after rochcor"<<endl;
      
      muons.push_back(*thisMuon);
    }
  }

  sort(muons.begin(), muons.end(), P4SortCondition);


  hists->fill1DHist(muons.size(),     Form("size_mu_cut%i",  1),";Number of muons",    5,0,5, 1, "Muons");
  //hists->fill1DHist(electrons.size(), Form("size_el_cut%i",  1),";Number of electrons",5,0,5, 1, "Electrons");
  //hists->fill1DHist(electrons0.size(),Form("size_el0_cut%i", 1),";Number of electrons",5,0,5, 1, "Electrons");
  //hists->fill1DHist(photons.size(),   Form("size_ph_cut%i",  1),";Number of photons",  5,0,5, 1, "Photon");
  //hists->fill1DHist(photons0.size(),  Form("size_ph0_cut%i", 1),";Number of photons",  5,0,5, 1, "Photon");


  if(!isRealData && makeGen){
    hists->fill1DHist(genMll,  "gen_Mll_reco_gamma",  ";gen_Mll",100,0,mllMax, 1,"eff");
    hists->fill1DHist(gendR,   "gen_dR_reco_gamma",   ";gen_dR", 50,0,0.3,1,"eff");
  }


  if (checkTrigger){
    triggerSelector->SelectTrigger(myTrigger, triggerStatus, hltPrescale, isFound, triggerPass, prescale);
    if (!triggerPass) return kTRUE;
  }

  if (gamma.Pt() < cut_gammapt) return kTRUE;
  if (fabs(gamma.Eta()) > 1.444) return kTRUE;
  //if (fabs(gamma.SCEta()) > 1.444) return kTRUE;

  if(!isRealData && makeGen){
    hists->fill1DHist(genMll,  "gen_Mll_reco_gamma_iso",  ";gen_Mll",100,0,mllMax, 1,"eff");
    hists->fill1DHist(gendR,   "gen_dR_reco_gamma_iso",   ";gen_dR", 50,0,0.3,1,"eff");
  }

  FillHistoCounts(3, eventWeight);
  CountEvents(3);

  if (muons.size()<2) return kTRUE;
  
  Float_t tmpMll = 99999;
  if (muons.size()==2){
    lPt1 = muons[0];
    lPt2 = muons[1];
  }
  else 
    for (UInt_t m1=0; m1 < muons.size(); m1++){
      for (UInt_t m2=m1; m2 < muons.size(); m2++){
        if (muons[m1].Charge()*muons[m2].Charge()==1) continue; //same charge - don't care
        
        if( (muons[m1]+muons[m2]).M() < tmpMll)
          {
            lPt1 = muons[m1];
            lPt2 = muons[m2];
            tmpMll = (muons[m1]+muons[m2]).M();
          }
      }
    }
  
  if (lPt1.Charge()==1 && lPt2.Charge()==-1){
    l1 = lPt1;
    l2 = lPt2;
  }
  else if (lPt1.Charge()==-1 && lPt2.Charge()==1){
    l1 = lPt2;
    l2 = lPt1;
  }
  else{
    //cout<<"ch1 = "<<muons[0].Charge()<<"  ch2 = "<<muons[1].Charge()<<endl;
    return kTRUE;//Abort("reco * They are the same charge!");
  }
    
  Double_t Mll = (l1+l2).M();

  hists->fill1DHist(Mll,  Form("diLep_mass_low_cut%i", 3), ";M(ll)", 50, 0,20,  1, "");
  hists->fill1DHist(Mll,  Form("diLep_mass_high_cut%i", 3),";M(ll)", 50, 0,120, 1, "");

  
  //if (lPt1.Pt() > 20 && zgamma::CalculateMuonIso(&muons[0]) > 0.4)
  //return kTRUE;
  //if (lPt2.Pt() > 20 && zgamma::CalculateMuonIso(&muons[1]) > 0.4)
  //return kTRUE;

  //if (zgamma::CalculateMuonIso(&muons[1]) > 0.4)
  //return kTRUE;
  
  MakePhotonPlots(gamma);
  
  if (lPt1.Pt() < cut_l1pt || lPt2.Pt() < cut_l2pt_low)   return kTRUE;


  //if (!isfinite(eventWeight))
  //cout<<"nan/inf "<<eventWeight<<endl;
  //    if (fabs(eventWeight)>50 || fabs(eventWeight)<0.00001)

  if (!isRealData){
    eventWeight *= weighter->MuonSF(lPt1);
    eventWeight *= weighter->MuonSF(lPt2);
  }

  if (!isfinite(eventWeight))
    cout<<"nan after muon sf  "<<eventWeight<<endl;

  if(!isRealData && makeGen){
    hists->fill1DHist(genMll,  "gen_Mll_two_lep_reco", ";gen_Mll",100,0,mllMax, 1,"eff");
    hists->fill1DHist(gendR,   "gen_dR_two_lep_reco",  ";gen_dR", 50,0,0.3,1,"eff");
  }

  FillHistosFull(4, eventWeight, l1, l2, lPt1, lPt2, gamma, "");
  FillHistoCounts(4, eventWeight);
  CountEvents(4);


  //for(Int_t ev=0; ev<evSize;ev++){
  //if (eventNumber==hisEVTS[ev]){cout<<eventNumber<<" Found an event after lept selection "<<endl;
  //  cout<<"Mll = "<<Mll<<endl;
  //  break;}
  //}

  if (Mll > 20) return kTRUE;

  fout<<runNumber<<" "<<eventNumber<<endl;
  
  
  FillHistosFull(5, eventWeight, l1, l2, lPt1, lPt2, gamma, "");
  FillHistoCounts(5, eventWeight);
  CountEvents(5);
  
  if (lPt1.DeltaR(gamma)<1.0 || lPt2.DeltaR(gamma)<1.0) return kTRUE;

  FillHistosFull(6, eventWeight, l1, l2, lPt1, lPt2, gamma, "");
  FillHistoCounts(6, eventWeight);
  CountEvents(6);
  //if (fabs(gamma.Eta())<1.444)
  //FillHistosFull(6, eventWeight, l1, l2, lPt1, lPt2, gamma, "EB");
  //else
  //FillHistosFull(6, eventWeight, l1, l2, lPt1, lPt2, gamma, "EE");

  global_Mll = Mll;
  Float_t Mllg = (l1+l2+gamma).M();

  MakeMuonPlots(muons[0], pvPosition);
  MakeMuonPlots(muons[1], pvPosition);
  


  if (Mll>2.9 && Mll<3.3){ //jpsi window
    FillHistosFull(12, eventWeight, l1, l2, lPt1, lPt2, gamma, "jpsi");
    FillHistoCounts(12, eventWeight);
    CountEvents(12);
    //return kTRUE;
  }


  if ( (Mll>2.9 && Mll<3.3) || (Mll>9.3 && Mll<9.7)) return kTRUE; //jpsi and upsilon removeal
  FillHistosFull(7, eventWeight, l1, l2, lPt1, lPt2, gamma, "");
  FillHistoCounts(7, eventWeight);
  CountEvents(7);

  for (UInt_t i =0; i<ntrig; i++){
    triggerSelector->SelectTrigger(myTriggers[i], triggerStatus, hltPrescale, isFound, triggerPass, prescale);
    if(triggerPass) nEventsTrig[7][i]++;

    if (!isFound)
      cout<<"TRG **** Warning ***\n The trigger name "<<myTriggers[i]<<" is not in the list of trigger names"<<endl;
  }


  if(!isRealData && makeGen && sample=="dalitz"){
    hists->fill1DHist(gen_l1.DeltaR(l1),"reco_gen_l1_deltaR","reco_gen_l1_deltaR",100,0,5, 1,"");
    hists->fill1DHist(gen_l2.DeltaR(l2),"reco_gen_l2_deltaR","reco_gen_l2_deltaR",100,0,5, 1,"");

    hists->fill2DHist(gen_l1.DeltaR(l1), gen_l2.DeltaR(l2),      "reco_gen_2D_ll_deltaR",     "reco_gen_2D_ll_deltaR",     100,0,5, 100,0,5, 1,"");
    hists->fill2DHist(gen_l1.DeltaR(l1), gen_gamma.DeltaR(gamma),"reco_gen_2D_l1gamma_deltaR","reco_gen_2D_l1gamma_deltaR",100,0,5, 100,0,5, 1,"");

    hists->fill1DHist(gen_gamma.DeltaR(gamma),"reco_gen_gamma_deltaR","reco_gen_gamma_deltaR",100,0,5, 1,"");
  }



  if (Mllg>122 && Mllg<128){
    FillHistosFull(13, eventWeight, l1, l2, lPt1, lPt2, gamma, "");
    FillHistoCounts(13, eventWeight);
    CountEvents(13);
  }

  if(!isRealData && makeGen){
    for (Int_t i = 1; i<=6; i++){
      if ((l1+l2).Pt()>20+i*5)
        hists->fill1DHist(genMll,  Form("gen_Mll_%i",i),";gen_Mll", 100,0,mllMax, 1,"eff");
      if (gamma.Pt()>20+i*5)
        hists->fill1DHist(genMll,  Form("gen_Mll_%i",6+i),";gen_Mll", 100,0,mllMax, 1,"eff");
      if (gamma.Pt()/Mllg> 0.15+0.05*i)
        hists->fill1DHist(genMll,  Form("gen_Mll_%i",12+i),";gen_Mll", 100,0,mllMax, 1,"eff");
      if ((l1+l2).Pt()/Mllg> 0.15+0.05*i)
        hists->fill1DHist(genMll,  Form("gen_Mll_%i",18+i),";gen_Mll", 100,0,mllMax, 1,"eff");
      if ((l1+l2).Pt()/Mllg> 0.15+0.05*i   && gamma.Pt()/Mllg> 0.15+0.05*i)
        hists->fill1DHist(genMll,  Form("gen_Mll_%i",24+i),";gen_Mll", 100,0,mllMax, 1,"eff");
    }
  }

  
  
  if ((l1+l2).Pt()/Mllg < 0.30 || gamma.Pt()/Mllg < 0.30)
    //if ((l1+l2).Pt() < 40 || gamma.Pt() < 40)
    return kTRUE;

  FillHistosFull(8, eventWeight, l1, l2, lPt1, lPt2, gamma, "");
  FillHistoCounts(8, eventWeight);
  CountEvents(8);
  

  if (checkTrigger){
    triggerSelector->SelectTrigger(myTrigger, triggerStatus, hltPrescale, isFound, triggerPass, prescale);
    if (!triggerPass) return kTRUE;
    //  else Abort("Event trigger should mu or el");

  }

  FillHistosFull(9, eventWeight, l1, l2, lPt1, lPt2, gamma, "");
  FillHistoCounts(9, eventWeight);
  CountEvents(9);

  //fout<<" nEvt = "<<nEvents[0]<<" : Run/lumi/event = "<<runNumber<<"/"<<lumiSection<<"/"<<eventNumber<<endl;


  fit_m_llg   = Mllg;
  fit_m_ll    = Mll;
  fit_weight  = eventWeight;
  fit_phEta   = gamma.SCEta();
  fit_phPt   = gamma.Pt();
  fit_diPt   = (l1+l2).Pt();
  fit_isLowPt = false;
      
  if (lPt2.Pt() < cut_l2pt)
    fit_isLowPt = true;
  
  _fitTree->Fill();



  if (Mllg>110 && Mllg<170){
    FillHistosFull(10, eventWeight, l1, l2, lPt1, lPt2, gamma, "");
    FillHistoCounts(10, eventWeight);
    CountEvents(10);  

  }

  if (Mllg<122 || Mllg>128) return kTRUE;

  FillHistosFull(11, eventWeight, l1, l2, lPt1, lPt2, gamma, "");
  FillHistoCounts(11, eventWeight);
  CountEvents(11);



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

void zgamma::SlaveTerminate() {}

void zgamma::Terminate()
{
  cout<<" ** FOR PAS **"<<endl;
  cout<<"| CUT DESCRIPTION            |\t"<< "\t|"<<endl;
  cout<<"| 0: initial                 |\t"<< nEvents[0]  <<"\t|"<<float(nEvents[0])/nEvents[0]<<"\t|"<<endl;
  cout<<"| 1: n/a                     |\t"<< nEvents[1]  <<"\t|"<<float(nEvents[1])/nEvents[0]<<"\t|"<<endl;
  cout<<"| 2: n/a                     |\t"<< nEvents[2]  <<"\t|"<<float(nEvents[2])/nEvents[1]<<"\t|"<<endl;
  cout<<"| 3: trigger & reco gamma iso|\t"<< nEvents[3]  <<"\t|"<<float(nEvents[3])/nEvents[2]<<"\t|"<<endl;
  cout<<"| 4: reco muons, and pt      |\t"<< nEvents[4]  <<"\t|"<<float(nEvents[4])/nEvents[3]<<"\t|"<<endl;
  cout<<"| 5: Mll < 20                |\t"<< nEvents[5]  <<"\t|"<<float(nEvents[5])/nEvents[4]<<"\t|"<<endl;
  cout<<"| 6:  dR (l,g) > 1           |\t"<< nEvents[6]  <<"\t|"<<float(nEvents[6])/nEvents[5]<<"\t|"<<endl;
  cout<<"| 7: jpsi/Upsi veto          |\t"<< nEvents[7]  <<"\t|"<<float(nEvents[7])/nEvents[6]<<"\t|"<<endl;
  cout<<"| 8: pT/Mllg > 0.30          |\t"<< nEvents[8]  <<"\t|"<<float(nEvents[8])/nEvents[7]<<"\t|"<<endl;
  cout<<"| 9:  trigger 2-check        |\t"<< nEvents[9]  <<"\t|"<<float(nEvents[9])/nEvents[8]<<"\t|"<<endl;
  cout<<"| 10: 110<mllg<170           |\t"<< nEvents[10] <<"\t|"<<float(nEvents[10])/nEvents[8]<<"\t|"<<endl;
  cout<<"| 11: 122<mllg<128           |\t"<< nEvents[11] <<"\t|"<<float(nEvents[11])/nEvents[8]<<"\t|"<<endl;
  cout<<"| 12: in Jpsi before #7      |\t"<< nEvents[12] <<"\t|"<<float(nEvents[12])/nEvents[6]<<"\t|"<<endl;
  cout<<"| 13: before pT/Mllg in [122,128]|\t"<< nEvents[13] <<"\t|"<<float(nEvents[13])/nEvents[7]<<"\t|"<<endl;
  
  /*
  cout<<"| CUT DESCRIPTION             |\t"<< "\t|"<<endl;
  cout<<"| 0: initial          |\t"<< nEvents[0]  <<"\t|"<<float(nEvents[0])/nEvents[0]<<"\t|"<<endl;
  cout<<"| 1: gen matching     |\t"<< nEvents[1]  <<"\t|"<<float(nEvents[1])/nEvents[0]<<"\t|"<<endl;
  cout<<"| 2: al in acc        |\t"<< nEvents[2]  <<"\t|"<<float(nEvents[2])/nEvents[1]<<"\t|"<<endl;
  cout<<"| 3: reco gamma iso   |\t"<< nEvents[3]  <<"\t|"<<float(nEvents[3])/nEvents[2]<<"\t|"<<endl;
  cout<<"| 4: reco lep         |\t"<< nEvents[4]  <<"\t|"<<float(nEvents[4])/nEvents[3]<<"\t|"<<endl;
  cout<<"| 5: mll <20          |\t"<< nEvents[5]  <<"\t|"<<float(nEvents[5])/nEvents[4]<<"\t|"<<endl;
  cout<<"| 6: dR(l,g)>1        |\t"<< nEvents[6]  <<"\t|"<<float(nEvents[6])/nEvents[5]<<"\t|"<<endl;
  cout<<"| 7: pT/Mllg > 0.30  |\t"<< nEvents[7]  <<"\t|"<<float(nEvents[7])/nEvents[6]<<"\t|"<<endl;
  cout<<"| 8: no jpsi/Ups      |\t"<< nEvents[8]  <<"\t|"<<float(nEvents[8])/nEvents[7]<<"\t|"<<endl;
  cout<<"| 9:  trigger         |\t"<< nEvents[9]  <<"\t|"<<float(nEvents[9])/nEvents[7]<<"\t|"<<endl;
  cout<<"| 10: pTll+pTg>95     |\t"<< nEvents[10] <<"\t|"<<float(nEvents[10])/nEvents[7]<<"\t|"<<endl;
  cout<<"| 11: 100<mH<150      |\t"<< nEvents[11] <<"\t|"<<float(nEvents[11])/nEvents[7]<<"\t|"<<endl;
  cout<<"| 12: jpsi            |\t"<< nEvents[12] <<"\t|"<<float(nEvents[12])/nEvents[7]<<"\t|"<<endl;
  */  
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
  //if(period=="2012")
    //eleISO = (lep->IsoMap("pfChIso_R04") + TMath::Max(0.0, (Double_t)(lep->IsoMap("pfPhoIso_R04") + lep->IsoMap("pfNeuIso_R04") - rho25Factor*lep->IsoMap("EffArea_R04"))))/lep->Pt();
    //else if(period=="2011")
    //eleISO = (lep->IsoMap("SumPt_R04") + TMath::Max(0.0, (Double_t)(lep->IsoMap("EmIso_R04") + lep->IsoMap("HadIso_R04") - rhoFactor*TMath::Pi()*0.09)))/lep->Pt();
    //else cout<<"Ele iso - Period is not defined!"<<period<<endl;

  return eleISO;
}



bool zgamma::PassElectronIdAndIso(TCElectron *lep, elIdAndIsoCuts cuts, TVector3 *pv)
{
  bool pass = false;

  //Float_t eleISO = CalculateElectronIso(lep);

  if(((fabs(lep->Eta()) < 1.442
       && lep->PtError()/lep->Pt() < cuts.ptErrorOverPt[0]
       && lep->SigmaIEtaIEta() < cuts.sigmaIetaIeta[0]
       && fabs(lep->SCDeltaPhi()) < cuts.dPhiIn[0]
       && fabs(lep->SCDeltaEta()) < cuts.dEtaIn[0]
       && lep->HadOverEm() < cuts.HadOverEm[0]
       && fabs(lep->Dxy(pv)) < cuts.dxy[0]
       && fabs(lep->Dz(pv)) < cuts.dz[0]
       //&& eleISO < cuts.pfIso04[0]
       )||
      (fabs(lep->Eta()) > 1.556
       && lep->PtError()/lep->Pt() < cuts.ptErrorOverPt[1]
       && lep->SigmaIEtaIEta() < cuts.sigmaIetaIeta[1]
       && fabs(lep->SCDeltaPhi()) < cuts.dPhiIn[1]
       && fabs(lep->SCDeltaEta()) < cuts.dEtaIn[1]
       && lep->HadOverEm() < cuts.HadOverEm[1]
       && fabs(lep->Dxy(pv)) < cuts.dxy[1]
       && fabs(lep->Dz(pv)) < cuts.dz[1]
       //&& eleISO < cuts.pfIso04[1]
       ))
     //&& lep->ConversionVeto()
     ) pass = true;

  //pass = true;

  //cout<<"ele pass? "<<pass<<endl;
  return pass;
}


bool zgamma::PassElectronIdAndIsoMVA(TCElectron *lep)
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

float zgamma::CalculateMuonIso(TCMuon *lep)
{
  float muISO = (lep->PfIsoCharged() + TMath::Max(0.0, lep->PfIsoPhoton() + lep->PfIsoNeutral() - 0.5*lep->PfIsoPU()))/lep->Pt();

  //  Float_t muISO = (lep->IsoMap("pfChargedHadronPt_R04") +
  //TMath::Max(0.0, lep->IsoMap("pfNeutralHadronEt_R04") + lep->IsoMap("pfPhotonEt_R04") - 0.5*lep->IsoMap("pfPUPt_R04")))/lep->Pt();

  return muISO;
}


bool zgamma::PassMuonIdAndIso(TCMuon *lep, muIdAndIsoCuts cuts, TVector3 *pv)
{
  bool pass = false;

  //Float_t muISO =  CalculateMuonIso(lep);

  if(1
     && lep->IsPF() 
     && ((lep->IsTRK() && lep->NumberOfMatches() > 0) || lep->IsGLB())
     && fabs(lep->Dxy(pv))     < cuts.dxy
     && fabs(lep->Dz(pv))      < cuts.dz
     && (lep->Pt() < 20 || zgamma::CalculateMuonIso(lep) < 0.4)
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
  chIsoCor = ph->PfIsoCharged()-rhoFactor*chEA;
  nhIsoCor = ph->PfIsoNeutral()-rhoFactor*nhEA;
  phIsoCor = ph->PfIsoPhoton()-rhoFactor*phEA;
  //cout<<chIsoCor<<"  "<<nhIsoCor<<"  "<<phIsoCor<<endl;

}

bool zgamma::PassPhotonIdAndIso(TCPhoton *ph, phIdAndIsoCuts cuts, TVector3 *pv)
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


void zgamma::FillHistoCounts(Int_t num, Double_t weight)
{
  hists->fill1DHist(num, "evt_byCut",";cut;weighted events", nC+1, -1,nC, weight, "Counts");
  hists->fill1DHist(num, "evt_byCut_raw", ";cut;events",     nC+1, -1,nC, 1, "Counts");
}

void zgamma::FillHistosFull(Int_t num, Double_t weight,
                            TCPhysObject l1, TCPhysObject l2, TCPhysObject lPt1, TCPhysObject lPt2, TCPhysObject gamma, string dir)
{

  char * d = new char [dir.length()+1];
  std::strcpy (d, dir.c_str());

  TLorentzVector diLep = l1+l2;
  TLorentzVector tri = diLep + gamma;

  TLorentzVector lPt1Gamma = lPt1+gamma;
  TLorentzVector lPt2Gamma = lPt2+gamma;

  TLorentzVector gammaCM = gamma;
  TLorentzVector diLepCM = diLep;
  TVector3    b1  = -tri.BoostVector();
  gammaCM.Boost(b1);
  diLepCM.Boost(b1);


  hists->fill1DHist(weight,     Form("weight_%s_cut%i",         d, num),";weight",  200, 0,3, 1, dir);

  hists->fill1DHist(tri.M(),     Form("tri_mass_%s_cut%i",         d, num),";M(ll#gamma)",  100, 0,200,  weight, dir);
  hists->fill1DHist(tri.M(),     Form("tri_mass80_%s_cut%i",       d, num),";M(ll#gamma)",  100,80,200,  weight, dir);
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

  hists->fill1DHist(diLep.Pt()/tri.M(),  Form("diLep_ptOverMllg_%s_cut%i",   d, num),
                    ";di-Lepton p_{T}/M_{ll#gamma}",100, 0,1, weight, dir);
  hists->fill1DHist(gamma.Pt()/tri.M(), Form("gamma_ptOverMllg_%s_cut%i",   d, num),
                    ";Photon p_{T}/M_{ll#gamma}",  100, 0,1, weight, dir);
  hists->fill1DHist(gammaCM.E(),Form("gamma_Ecom_%s_cut%i", d, num),";E_{#gamma} in CoM",  50, 0,120,  weight, dir);
  hists->fill1DHist(gamma.E(),  Form("gamma_E_%s_cut%i",    d, num),";Photon Energy", 50, 0,200, weight, dir);
  hists->fill1DHist(gamma.Pt(), Form("gamma_pt_%s_cut%i",   d, num),";Photon Pt",  50, 0,120, weight, dir);
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
                    ";q_{T}^{ll} (GeV); p_{T}^{#gamma} (GeV)",      50, 0,100, 50,0,100, weight, dir);
  
  hists->fill2DHist(diLep.Pt()/tri.M(),gamma.Pt()/tri.M(),Form("h2D_diLep_vs_gamma_ptOverMllg_%s_cut%i",d, num),
                    ";q_{T}^{ll}/M_{ll#gamma}; p_{T}^{#gamma}/M_{ll#gamma}",  100, 0,1, 100,0,1, weight, dir);

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

    hists->fill1DHist(co1, Form("co1_cut%i", num), ";cos_lp",  100,-1,1, 1,"Angles");
    hists->fill1DHist(co2, Form("co2_cut%i", num), ";cos_lm,", 100,-1,1, 1,"Angles");
    hists->fill1DHist(co3, Form("co3_cut%i", num), ";cosTheta",100,-1,1, 1,"Angles");
    hists->fill1DHist(phi, Form("phi_cut%i", num), ";phi lp",  100, -TMath::Pi(), TMath::Pi(), 1,"Angles");
  }

  //hists->fill1DHist(nVtxTotal, Form("vtx_nPV_tot_%s_cut%i",   d, num), "vtx_nPV_tot",    40, 0,40, weight, dir);
  hists->fill1DHist(nVtx,      Form("vtx_nPV_raw_%s_cut%i",   d, num), "Raw;nPV",    40, 0,40,      1, dir);
  hists->fill1DHist(nVtx,      Form("vtx_nPV_weight_%s_cut%i",d, num), "Re-weighted;nPV", 40, 0,40, weight, dir);

  hists->fill1DHist(nDofVtx1, Form("vtx_ndof1_%s_cut%i",  d, num), ";First vertex ndof", 50, 0,200, weight, dir);
  hists->fill1DHist(nDofVtx2, Form("vtx_ndof2_%s_cut%i",  d, num), ";Second vertex ndof", 50, 0,200, weight, dir);


}

void zgamma::MakeMuonPlots(TCMuon mu, TVector3 *pv)
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
  hists->fill1DHist(mu.PtError()/mu.Pt(), "mu_ptErrorOverPt", ";ptErrorOverPt", 50, 0,0.2, 1, "Muons");

  Float_t muIso = CalculateMuonIso(&mu);
  hists->fill1DHist(muIso, "mu_iso", ";rel isolation, pu corr", 50, 0,0.7, 1, "Muons");
  hists->fill2DHist(muIso, global_Mll, "mu_iso_vs_Mll", ";pfIsolationR04;M(l_{1},l_{2})", 50, 0,0.7, 50, 0,20, 1, "Muons");
  //Float_t iso2 = (mu.IsoMap("pfChargedHadronPt_R04") + mu.IsoMap("pfNeutralHadronEt_R04") + mu.IsoMap("pfPhotonEt_R04"))/mu.Pt();
  //hists->fill1DHist(iso2, "mu_iso_official", ";pfIsolationR04", 50, 0,0.7, 1, "Muons");


}

void zgamma::MakePhotonPlots(TCPhoton ph)
{

  hists->fill1DHist(ph.R9(),            "ph_R9",           ";R9",        100,    0, 1,    1,"Photon");
  hists->fill1DHist(ph.E2x5()/ph.R9(),  "ph_E2OverE9",     ";E2OverE9",  100,    0, 0.5,  1,"Photon");
  hists->fill1DHist(ph.TrackVeto(),     "ph_trackVeto",    ";track veto",  3,    0, 3,    1,"Photon");
  hists->fill1DHist(ph.HadOverEm(),     "ph_HadOverEm",    ";HadOverEm", 100,    0, 0.05, 1,"Photon");
  //hists->fill1DHist(ph.NormChi2(),    "ph_NormChi2",     ";NormChi2",  100,    0, 0.01, 1,"Photon");
  hists->fill1DHist(ph.SCDeltaPhi(),    "ph_SCDeltaPhi",   ";SCDeltaPhi", 100, -0.2, 0.2, 1,"Photon");
  hists->fill1DHist(ph.SCDeltaEta(),    "ph_SCDeltaEta",   ";SCDeltaEta", 100,-0.02, 0.02,1,"Photon");
  hists->fill1DHist(ph.GetNCrystals(),  "ph_nCrystals",    ";nCrystals",     160, 0, 160, 1,"Photon");
  hists->fill1DHist(ph.SigmaIEtaIEta(), "ph_SigmaIEtaIEta",";SigmaIEtaIEta", 100, 0, 0.1, 1,"Photon");
  hists->fill1DHist(ph.SigmaIPhiIPhi(), "ph_SigmaIPhiIPhi",";SigmaIPhiIPhi", 100, 0, 0.1, 1,"Photon");
  hists->fill1DHist(ph.ConversionVeto(),"ph_ConversionVeto",";ConversionVeto", 3, 0, 3,   1,"Photon");
  //hists->fill1DHist(ph., "ph_",";", 3, 0, 3, 1, "Photon");
  //hists->fill1DHist(ph., "ph_",";", 3, 0, 3, 1, "Photon");
  //hists->fill1DHist(ph., "ph_",";", 3, 0, 3, 1, "Photon");

}


TCGenParticle * zgamma::GetPrimaryAncestor(TCGenParticle *p)
{
  TCGenParticle *a = p;
  while (a->Mother())
    a = a->Mother();
  return a;
}


void zgamma::DiscoverGeneology(TCGenParticle *p, ULong64_t ev)
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



void zgamma::MuonDump(TCMuon mu, TVector3 *pv)
{
  Float_t muISO = zgamma::CalculateMuonIso(&mu);
  //Float_t muISO = (mu.IsoMap("pfChargedHadronPt_R04") +
  //               TMath::Max(0.0, mu.IsoMap("pfNeutralHadronEt_R04") + mu.IsoMap("pfPhotonEt_R04") - 0.5*mu.IsoMap("pfPUPt_R04")))/mu.Pt();


  cout  << runNumber << " " << eventNumber << "  pt =" << mu.Pt()
        << " eta=" << mu.Eta() << " GLB=" << mu.IsGLB() << "  isPF=" << mu.IsPF() << " isTRK="<<mu.IsTRK()<<endl;
  cout  << "chi2_tr="<< mu.NormalizedChi2_tracker() << "  iso=" << muISO <<  "  dxy=" << mu.Dxy(pv) << " dz=" << mu.Dz(pv)
        << "\n  good =" << mu.IsGood()
    //<< "\n  trk =" << mu.TrackLayersWithMeasurement() << "  pix=" <<mu.PixelLayersWithMeasurement()
    //   << "\n " << mu.NormalizedChi2() << " " << mu.NumberOfValidMuonHits() << " " << mu.NumberOfMatchedStations()
    //  << " " << mu.NumberOfValidPixelHits() << " " << mu.TrackLayersWithMeasurement()
        << endl;
}


void zgamma::PhotonDump(TCPhoton pho, phIdAndIsoCuts cuts)
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

void zgamma::PhotonR9Corrector(TCPhoton& ph){
  //old R9 correction
  float R9Cor;
  R9Cor = ph.R9();
  
  if (fabs(ph.SCEta()) < 1.479)
    R9Cor = ph.R9()*1.0045 + 0.0010;
  else
    R9Cor = ph.R9()*1.0086 - 0.0007;
    
  ph.SetR9(R9Cor);
}

