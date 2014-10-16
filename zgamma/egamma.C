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

string  myTrigger = "";
Float_t cut_l1pt  = 20;
Float_t cut_l2pt  = 5;
Float_t cut_gammapt = 25;
Float_t mllMax = 25;

Bool_t checkTrigger = kTRUE;
//Bool_t applyPhosphor = 0;
Bool_t doRochCorr = 1;
Bool_t doZee = 0;
Bool_t makeApzTree = 1;

const UInt_t hisEVTS[] = {635,907};
Int_t evSize = sizeof(hisEVTS)/sizeof(int);

bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());}
bool SCESortCondition(const TCEGamma& p1, const TCEGamma& p2) {return (p1.SCEnergy() > p2.SCEnergy());}
void egamma::Begin(TTree * tree)
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

  if (selection == "zee")
    {
      doZee = true;
      trigger = "ee";
    }

  //cout<<sample<<selection<<trigger<<gen<<endl;

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

  histoFile = new TFile(Form("a_%s_higgsHistograms.root", selection.c_str()), "RECREATE");
  histoFile->mkdir("apzTree", "apzTree");
  hists    = new HistManager(histoFile);
  HM       = new HistMaker(hists);
  ObjID    = new ObjectID();
  weighter = new WeightUtils(sample, period, selection, 0);
  roch     = new rochcor2012();
  triggerSelector = new TriggerSelector("", period, *triggerNames);
  myRandom = new TRandom3();

  for (Int_t n1=0; n1<nC; n1++){
    nEvents[n1]=0;
    for (UInt_t n2=0; n2<ntrig; n2++)
      nEventsTrig[n1][n2]=0;
  }

  //Angles class:
  ang = new ZGAngles();

  fcuts.open("./out_cutlist.txt", ofstream::out);
  fout.open("./out_synch_.txt",ofstream::out);
  fout.precision(3); fout.setf(ios::fixed, ios::floatfield);
  cout.precision(5); //fout.setf(ios::fixed, ios::floatfield);

  histoFile->cd("apzTree");
  if (makeApzTree){
    _apzTree = new TTree("apzTree", "A tree for studying new particles");
    _apzTree->Branch("weight",&apz_w,     "weight/D");
    _apzTree->Branch("dr1234",&apz_dr1234,"dr1234/D");
    _apzTree->Branch("dr12",  &apz_dr12,  "dr12/D");
    _apzTree->Branch("dr13",  &apz_dr13,  "dr13/D");
    _apzTree->Branch("dr23",  &apz_dr23,  "dr23/D");
    _apzTree->Branch("dr34",  &apz_dr34,  "dr34/D");
    _apzTree->Branch("pt12",  &apz_pt12,  "pt12/D");
    _apzTree->Branch("pt34",  &apz_pt34,  "pt34/D");
    _apzTree->Branch("m12",   &apz_m12,   "m12/D");
    _apzTree->Branch("m123",  &apz_m123,  "m123/D");
    _apzTree->Branch("m34",   &apz_m34,   "m34/D");
    _apzTree->Branch("m4l",   &apz_m4l,   "m4l/D");
    _apzTree->Branch("pt1",   &apz_pt1,   "pt1/D");
    _apzTree->Branch("pt2",   &apz_pt2,   "pt2/D");
    _apzTree->Branch("pt3",   &apz_pt3,   "pt3/D");
    _apzTree->Branch("pt4",   &apz_pt4,   "pt4/D");

    _apzTree->Branch("eta1234",&apz_eta1234,"eta1234/D");
    _apzTree->Branch("eta12",  &apz_eta12,  "eta12/D");
    _apzTree->Branch("eta34",  &apz_eta34,  "eta12/D");
    _apzTree->Branch("eta1",   &apz_eta1,   "eta1/D");
    _apzTree->Branch("eta2",   &apz_eta2,   "eta2/D");
    _apzTree->Branch("eta3",   &apz_eta3,   "eta3/D");
    _apzTree->Branch("eta4",   &apz_eta4,   "eta4/D");

    _apzTree->Branch("run",  &runNumber,  "run/i");
    _apzTree->Branch("event",&eventNumber,"event/l");
  }

  //----------------------//
  // Selecting triggers   //
  //----------------------//

  if (trigger=="pho")
    {
      myTrigger = "HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v";
      //cut_l1pt = 20;
      //cut_l2pt = 1;
      //cut_gammapt = 38;
    }
  else if (trigger=="pho26")
    {
      myTrigger = "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v";
      //cut_l1pt = 23;
      //cut_l2pt = 4;
      //cut_gammapt = 25;
    }
  else if (trigger=="ee")
    {
      myTrigger = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v";
      cut_l1pt = 25;
      cut_l2pt = 7;
      cut_gammapt = 23;
    }
  else
    checkTrigger = kFALSE;
  //cout<<"Other options are not supported!"<<endl;

}

void egamma::SlaveBegin(TTree * /*tree*/)
{   TString option = GetOption(); }

Bool_t egamma::Process(Long64_t entry)
{
  GetEntry(entry);
  CountEvents(0, "Ntuple events",fcuts);
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
  ObjID->SetEventInfo(isRealData, runNumber, eventNumber, rhoFactor);
  float phoMvaScore;
  float DAleMvaScore;

  //cout<<eventNumber<<"  ev rho = "<<rhoFactor<<endl;
  hists->fill1DHist(nPUVerticesTrue, "nPUVerticesTrue",";nPUVerticesTrue",  100, 0,100, 1,"GEN");
  hists->fill1DHist(nPUVertices,     "nPUVertices",";nPUVertices",          100, 0,100, 1,"GEN");

  nVtx = primaryVtx->GetSize();
  HM->Reset(rhoFactor, nVtx);

  // --- Primary vertex ---//
  if (nVtx<1) return kTRUE;
  TCPrimaryVtx *mainPrimaryVertex = 0; mainPrimaryVertex   = (TCPrimaryVtx*)(primaryVtx->At(0));
  TVector3* pvPosition; pvPosition = mainPrimaryVertex;

  HM->SetVtx(*mainPrimaryVertex);


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
	if (ObjID->GetPrimaryAncestor(thisParticle)->GetPDGId()==A)
	  {
	    hists->fill1DHist(thisParticle->M(), "gen_mu_mass",";gen mu mass",  200, -1,1, 1,"GEN");
	    gen_mu.push_back(*thisParticle);
	  }
	//this is for the Madgraph case, where there are events with no Higgs in LHE file
	//:(while there are always should be):
	else if (h==0 && thisParticle->Mother() //good boys should have a mother
		 && ObjID->GetPrimaryAncestor(thisParticle)->GetPDGId() == thisParticle->GetPDGId())
	  {
	    //cout<<"   --- this one --->>>"<<endl;
	    //ObjID->DiscoverGeneology(thisParticle);
	    //cout<<"   <<--- this one"<<endl;
	    hists->fill1DHist(thisParticle->M(), "gen_mu_mass",";gen mu mass",  200, -1,1, 1,"GEN");
	    gen_mu.push_back(*thisParticle);
	  }
      }


      //GEN ELECTRONS
      if (abs(thisParticle->GetPDGId()) == 11 && thisParticle->GetStatus()==1) {
	if (ObjID->GetPrimaryAncestor(thisParticle)->GetPDGId()==A)
	  gen_el.push_back(*thisParticle);
	else if (h==0 && thisParticle->Mother()
		 && ObjID->GetPrimaryAncestor(thisParticle)->GetPDGId() == thisParticle->GetPDGId())
	  gen_el.push_back(*thisParticle);
      }


      //PHOTON from the Higgs
      if (sample=="dalitz" && thisParticle->GetPDGId()==22 && thisParticle->GetStatus()==1
	  //&& ObjID->GetPrimaryAncestor(thisParticle)->GetPDGId()==A
	  && thisParticle->Mother() &&  abs(thisParticle->Mother()->GetPDGId())!=13 && abs(thisParticle->Mother()->GetPDGId())!=11)
	{
	  gen_gamma = *thisParticle;
	  ph++;
	}
      else if (sample=="DY" &&
	       thisParticle->GetPDGId()==22 && thisParticle->GetStatus()==1 &&
	       ObjID->GetPrimaryAncestor(thisParticle)->GetPDGId()==23)
	gen_gamma = *thisParticle;


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


    if(selection=="el"){
      if (gen_el.size()!=2) return kTRUE;
      //Abort(Form("ev #%i NONONO There has to be exactly 2 muons from the Higgs! \n \t\t but there are %i", eventNumber, (int)gen_el.size()));

      gen_lPt1 = gen_el[0];
      gen_lPt2 = gen_el[1];
      // If there are more muons - it's confusing, so let's just ignore them for now.
      // (they may come from ZH -> ZZgamma process)

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

    }


    if (gen_el.size()!=2 && gen_el.size()!=2)
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
    CountEvents(1,"Pass Gen acceptance",fcuts);

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
	 gen_lPt2.Pt()>cut_l2pt && fabs(gen_lPt2.Eta())<2.4
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
  CountEvents(2, "VTX and after GEN", fcuts);

  //if (checkTrigger){
  //triggerSelector->SelectTrigger(myTrigger, triggerStatus, hltPrescale, isFound, triggerPass, prescale);
  //if (!triggerPass) return kTRUE;
  //}

  if (checkTrigger){
    Bool_t pass = false;
    for (Int_t i=0; i<6; i++){
      triggerSelector->SelectTrigger(pho26OR[i], triggerStatus, hltPrescale, isFound, triggerPass, prescale);
      if (triggerPass) {pass = true; break;}
    }
    if (!pass) return kTRUE;
    // else Abort("Event trigger should mu or el");
  }

  FillHistoCounts(3, eventWeight);
  CountEvents(3,"Pass Trigger",fcuts);

  vector<TCElectron> electrons0, electrons, DAlectrons, fake_electrons0;
  vector<TCMuon> muons0, muons;
  vector<TCPhoton> photons0, photonsTight, photonsHZG, photonsMVA, fake_photons;
  TCPhysObject l1,l2, lPt1, lPt2;
  TCPhoton gamma0, gamma, ufoton, ufelectron;
  TCElectron DALE;


  for (Int_t i = 0; i < recoPhotons->GetSize(); ++i) {
    TCPhoton* thisPhoton = (TCPhoton*) recoPhotons->At(i);
    if (fabs(thisPhoton->SCEta()) > 1.5) continue;

    if(thisPhoton->Pt() > 15){
      if (isRealData || !makeGen || (sample=="dalitz" && makeGen && gen_gamma.DeltaR(*thisPhoton) < 0.1) )
	photons0.push_back(*thisPhoton);
    }
  }

  sort(photons0.begin(), photons0.end(), P4SortCondition);


  for (Int_t i = 0; i < recoElectrons->GetSize(); ++i) {
    TCElectron* thisElec = (TCElectron*) recoElectrons->At(i);
    if (!(fabs(thisElec->Eta()) < 2.5)) continue;
    if (thisElec->Pt() > 15.0) {
      //thisElec->DeltaR(gamma) > 0.2 ) {

      if(isRealData || !makeGen ||
	 ( selection=="el" && makeGen && (gen_l1.DeltaR(*thisElec) < 0.1 || gen_l2.DeltaR(*thisElec)<0.1)  && thisElec->Pt()>10))
	{
	  if (ObjID->PassElectronIdAndIsoMVA(*thisElec)){
	    thisElec->SetIdMap("mvaScore", 0);
	    electrons.push_back(*thisElec);
	  }
	}
      if (ObjID->PassElectronIdAndIso(*thisElec, pvPosition, "Loose") &&
	  thisElec->Pt() > 15.0){
	thisElec->SetIdMap("mvaScore", 0);
	electrons0.push_back(*thisElec);
      }

      TCPhoton match;
      Bool_t isMatched = 0;
      for(UInt_t p=0; p<photons0.size(); p++)
	{
	  if (fabs(thisElec->SCEta()-photons0[p].SCEta()) < 0.0001){
	    match = photons0[p];
	    thisElec->SetSCRawEnergy(match.SCRawEnergy());
	    isMatched = 1;
	    //cout<<" Ele matched to photon! \t ev="<<nEvents[0]-1<<endl;
	    //cout<<"PHo: SCEne = "<<match.SCEnergy()<<" pt = "<<match.Pt()<<"  SCEta = "<<match.SCEta()<<"  SCPhi = "<<match.SCPhi()<<endl;
	    //cout<<"Ele: SCEne = "<<thisElec->SCEnergy()<<" pt = "<<thisElec->Pt()<<"   SCEta = "<<thisElec->SCEta()<<"  SCPhi = "<<thisElec->SCPhi()<<endl;
	    break;
	  }
	}

      if (thisElec->Pt() > 30 && fabs(thisElec->Eta()) < 1.4442 &&
	  ObjID->PassDalitzEleID(*thisElec, pvPosition, "MVA", DAleMvaScore) ){
	thisElec->SetIdMap("mvaScore", DAleMvaScore);

	//cout<<"reco ele: "<<i<<endl;
	if (isMatched && ObjID->HggPreselection(match))
	  DAlectrons.push_back(*thisElec);
      }

      bool idpass = ObjID->PassDalitzEleID(*thisElec, pvPosition, "MVA", DAleMvaScore);
      thisElec->SetIdMap("mvaScore", DAleMvaScore);

      /*
      for(Int_t ev=0; ev<evSize;ev++){
	if (nEvents[0]-1 == hisEVTS[ev]){cout<<nEvents[0]-1<<" evNumber = "<<eventNumber<<endl;
	  cout<<"gamma: pt = "<<gamma.Pt()<<"  SCEta = "<<gamma.SCEta()<<"  SCPhi = "<<gamma.SCPhi()<<endl;
	  cout<<"\t R9 = "<<gamma.R9()<<"  hadoverEm = "<<gamma.HadOverEm()<<endl;
	  cout<<"Ele  : pt = "<<thisElec->Pt()<<"   SCEta = "<<thisElec->SCEta()<<"  SCPhi = "<<thisElec->SCPhi()<<"  MVA = "<<thisElec->IdMap("mvaScore")<<endl;
	  cout<<"\t R9 = "<<thisElec->R9()<<"  hadoverEm = "<<thisElec->HadOverEm()<<endl;
	  cout<<"\t SCEnergy: "<<thisElec->SCEnergy()<<"   SCRaw energy: "<<thisElec->SCRawEnergy()<<endl;
	  cout<<"\t matched to pho: "<<isMatched<<"  pass Hgg pre: "<<ObjID->HggPreselection(match)<<endl;
	  vector<TCEGamma> sc0 = thisElec->BaseSC();
	  sort(sc0.begin(),  sc0.end(),  SCESortCondition);
	  const vector<TCElectron::Track> trk0 = thisElec->GetTracks();

	  cout<<" n tracks="<<trk0.size() <<"  N SC="<<sc0.size()<<endl;

	  for (int s = 0; s <sc0.size(); s++)
	    cout<<"sc "<<s<<"  SCEnergy: "<<sc0[s].SCEnergy()<<"  eta: "<<sc0[s].SCEta()<<" phi: "<<sc0[s].SCPhi()
		<<"  sieie: "<<sc0[s].SigmaIEtaIEta()<<"  e25/e55: "<<sc0[s].E2x5()/sc0[s].E5x5()<<endl;
	  //cout<<"gsf trk1, pt ="<<trk0[0].Pt()<<"  eta= "<<trk0[0].Eta()<<endl;
	  //cout<<"gsf trk2, pt ="<<trk0[1].Pt()<<"  eta= "<<trk0[1].Eta()<<endl;
	  //cout<<"SCRawEn = "<<thisElec->SCEnergy()<<"  trk1+trk2 Pt ="<<(trk[0]+trk[1]).Pt()
	  //  <<" EoverPt = "<<thisElec->SCEnergy()/(trk[0]+trk[1]).Pt()<<endl;
	  //cout<<"dzy = "<<thisElec->Dxy(pvPosition)<<"  dz = "<<thisElec->Dz(pvPosition)<<endl;
	  //thisElec->Dump();
	  //ObjID->ElectronDump(*thisElec);

	  break;}
      }
      */

    }

  }

  sort(electrons0.begin(), electrons0.end(), P4SortCondition);
  sort(electrons.begin(),  electrons.end(),  P4SortCondition);
  sort(DAlectrons.begin(), DAlectrons.end(), P4SortCondition);

  if (electrons0.size()>0)
    HM->MakeElectronPlots(electrons0[0], "AnElectron");

  if (DAlectrons.size()<1)
    return kTRUE;

  DALE = DAlectrons[0];
  if (DALE.Pt()<30)
    return kTRUE;

  HM->SetDalectron(DALE);
  FillHistoCounts(4, eventWeight);
  CountEvents(4,"DALectron selection",fcuts);

  vector<TCElectron::Track> trk = DALE.GetTracks();
  sort(trk.begin(), trk.end(), P4SortCondition);

  float xx[6] = {-1,-1,-1,-1,-1,-1};
  TCTrack::ConversionInfo cI = trk[0].GetConversionInfo();
  //if (cI.vtxProb>1e-6 && cI.lxyBS>2 && cI.nHitsMax>1)
  //return kTRUE;
  xx[0] = cI.vtxProb;
  xx[1] = cI.lxyBS;
  xx[2] = cI.nHitsMax;

  cI = trk[1].GetConversionInfo();
  //if (cI.vtxProb>1e-6 && cI.lxyBS>2 && cI.nHitsMax>1)
  //return kTRUE;
  xx[3] = cI.vtxProb;
  xx[4] = cI.lxyBS;
  xx[5] = cI.nHitsMax;


  l1 = trk[0];
  l2 = trk[1];

  if (l1.Pt() > l2.Pt()){
    lPt1 = l1; lPt2 = l2;}
  else {lPt1 = l2; lPt2 = l1;}

  HM->SetLeptons(lPt1, lPt2);

  if (DALE.BaseSC()[0].Pt()<DALE.BaseSC()[1].Pt())
    cout<<"Warinng SC are  not ordered in pt!"<<endl;

  if (!l1.IdMap("GSF_PassConv")) return kTRUE;
  if (!l2.IdMap("GSF_PassConv")) return kTRUE;

  //if (l1.IdMap("GSF_MissHits")!=0) return kTRUE;
  //if (l2.IdMap("GSF_MissHits")!=0) return kTRUE;

  FillHistoCounts(5, eventWeight);
  CountEvents(5,"Conv. Veto etc",fcuts);

  if (lPt1.Pt() + lPt2.Pt() < 44) return kTRUE;
  //if (lPt1.Pt() < cut_l1pt || lPt2.Pt() < cut_l2pt)   return kTRUE;

  FillHistoCounts(6, eventWeight);
  CountEvents(6,"pt1+pt2 > 44",fcuts);
  HM->MakeElectronPlots(DALE, "DalitzEle-cut6");

  Double_t scale = DALE.SCEnergy()/DALE.E();
  DALE.SetPxPyPzE(scale*DALE.Px(), scale*DALE.Py(),scale*DALE.Pz(), DALE.SCEnergy());


  for (Int_t i = 0; i < recoPhotons->GetSize(); ++i) {
    //cout<<"new photon!!!!!!!"<<endl;
    TCPhoton* thisPhoton = (TCPhoton*) recoPhotons->At(i);
    //cout<<"in photon loop event = "<<eventNumber
    //  <<"\n pt="<<thisPhoton->Pt()<<" eta="<<thisPhoton->Eta()<<" phi="<<thisPhoton->Phi()<<" px="<<thisPhoton->Px()<<endl;

    if (isRealData || !makeGen || (makeGen && gen_lPt1.DeltaR(*thisPhoton) < 0.1 && thisPhoton->Pt()>20))
      fake_photons.push_back(*thisPhoton);

    //if (!(fabs(thisPhoton->Eta()) < 2.5)) continue;
    if (fabs(thisPhoton->SCEta()) > 1.4442) continue;
    if (DALE.DeltaR(*thisPhoton) < 0.2) continue;

    if(thisPhoton->Pt() > 15){
      //for(Int_t ev=0; ev<evSize;ev++){
      //if (nEvents[0]-1 == hisEVTS[ev]){cout<<nEvents[0]-1<<" evNumber="<<eventNumber<<endl;
      //  cout<<"gamma: pt = "<<thisPhoton->Pt()<<"  SCEta = "<<thisPhoton->SCEta()<<"  SCPhi = "<<thisPhoton->SCPhi()<<endl;
      //  break;}
      //}

      if(ObjID->PassPhotonIdAndIso(*thisPhoton, "CutBased-MediumWP", phoMvaScore)
	 //&& (isRealData || !makeGen || (sample=="dalitz" && makeGen && gen_gamma.DeltaR(*thisPhoton) < 0.2) )
	 )
	{
	  Double_t corrPhoEnSC   = correctedPhotonEnergy(thisPhoton->SCEnergy(), thisPhoton->SCEta(), thisPhoton->R9(), runNumber,
							 0, "Moriond2014", !isRealData, myRandom);
	  Double_t corrPhoEnReco = correctedPhotonEnergy(thisPhoton->E(), thisPhoton->SCEta(), thisPhoton->R9(), runNumber,
							 0, "Moriond2014", !isRealData, myRandom);

	  HM->MakePhotonEnergyCorrPlots(*thisPhoton, corrPhoEnReco, corrPhoEnSC);

	  //cout<<runNumber<<"  uncor en="<< thisPhoton->E()<<"  cor en = "<<corrPhoEn<<endl;
	  Double_t scale = 1;
	  //Double_t scale = corrPhoEnReco/thisPhoton->E();
	  //Double_t scale = corrPhoEnSC/thisPhoton->E();

	  thisPhoton->SetXYZM(scale*thisPhoton->Px(), scale*thisPhoton->Py(),scale*thisPhoton->Pz(),0);

	  //cout<<runNumber<<"  uncor en="<< thisPhoton->E()<<"  cor en = "<<corrPhoEn<<endl;
	  //thisPhoton->SetPtEtaPhiM(corrPhoEnSC/cosh(SCEta()));

	  thisPhoton->SetIdMap("mvaScore", phoMvaScore);
	  photonsHZG.push_back(*thisPhoton);
	}

      if(ObjID->PassPhotonIdAndIso(*thisPhoton, "CutBased-TightWP", phoMvaScore)){
	thisPhoton->SetIdMap("mvaScore", phoMvaScore);
	photonsTight.push_back(*thisPhoton);
      }

      if(ObjID->PassPhotonIdAndIso(*thisPhoton, "MVA", phoMvaScore)){
	thisPhoton->SetIdMap("mvaScore", phoMvaScore);
	photonsMVA.push_back(*thisPhoton);
      }
    }
  }

  sort(photonsTight.begin(), photonsTight.end(), P4SortCondition);
  sort(photonsHZG.begin(),   photonsHZG.end(),   P4SortCondition);
  sort(photonsMVA.begin(),   photonsMVA.end(),   P4SortCondition);
  sort(fake_photons.begin(), fake_photons.end(), P4SortCondition);

  //if (fake_photons.size()>1)
  // Abort("Nah, this is too much.");
  if (fake_photons.size()>=1)
    ufoton = fake_photons[0];

  //if (photonsTight.size()<1) return kTRUE;
  //gamma = photonsTight[0];

  if (photonsHZG.size()<1) return kTRUE;
  gamma = photonsHZG[0];

  HM->MakePhotonPlots(gamma);
  HM->SetGamma(gamma);

  //if (photonsMVA.size()<1) return kTRUE;
  //gamma = photonsMVA[0];

  if (!isRealData)
    eventWeight *= weighter->PhotonSF(gamma);

  //if (photons0.size()<1) return kTRUE;
  //gamma0 = photons0[0];
  //if (gamma0.Pt() < cut_gammapt) return kTRUE;

  if (gamma.Pt() < cut_gammapt) return kTRUE;
  if (fabs(gamma.SCEta()) > 1.4442) return kTRUE;


  if(!isRealData && makeGen){
    hists->fill1DHist(genMll,  "gen_Mll_reco_gamma",  ";gen_Mll",100,0,mllMax, 1,"eff");
    hists->fill1DHist(gendR,   "gen_dR_reco_gamma",   ";gen_dR", 50, 0,0.3,    1,"eff");
  }

  if(!isRealData && makeGen){
    hists->fill1DHist(genMll,  "gen_Mll_reco_gamma_iso",  ";gen_Mll",100,0,mllMax, 1,"eff");
    hists->fill1DHist(gendR,   "gen_dR_reco_gamma_iso",   ";gen_dR", 50, 0,0.3,    1,"eff");
  }

  HM->FillHistosFull(7, eventWeight);
  FillHistoCounts(7, eventWeight);
  CountEvents(7,"Photon pt>30; |eta|<1.4442",fcuts);

  FillHistoCounts(8, eventWeight);
  CountEvents(8,"Eight",fcuts);



  if (doZee){
    if (electrons.size()<2 || photons0.size()<2) return kTRUE;
    //cout<<"doing zee.  Enough photons and electrons!"<<endl;

    if (checkTrigger){
      triggerSelector->SelectTrigger(myTrigger, triggerStatus, hltPrescale, isFound, triggerPass, prescale);
      if (!triggerPass) return kTRUE;
    }

    //cout<<"passed the trigger"<<endl;
    Int_t match1 = -1, match2 = -1;
    for (UInt_t p=0; p<photons0.size(); p++)
      {
	if (photons0[p].Pt() < 40 || fabs(photons0[p].SCEta())>1.444) continue;
	for (UInt_t e=0; e<electrons.size(); e++)
	  {
	    if (photons0[p].DeltaR(electrons[e]) < 0.2)
	      {
		//cout<<"matched! photon "<<p<<"  to electron "<<e<<endl;
		if (match1==-1) match1=p;
		else {match2=p; break;}
	      }
	  }
      }

    if (match1!=-1 && match2!=-1)
      HM->MakeZeePlots(photons0[match1], photons0[match2]);

    return kTRUE;
  }


  for (Int_t i = 0; i < recoMuons->GetSize(); ++ i) {
    TCMuon* thisMuon = (TCMuon*) recoMuons->At(i);
    //cout<<"event = "<<eventNumber
    //  <<"\n pt="<<thisMuon->Pt()<<" eta="<<thisMuon->Eta()<<" phi="<<thisMuon->Phi()<<" px="<<thisMuon->Px()<<endl;

    if (!(fabs(thisMuon->Eta()) < 2.4)) continue;
    if(thisMuon->Pt() > 4 && ObjID->PassMuonId(*thisMuon, pvPosition, "Soft") ) {
      if (doRochCorr){
	Float_t qter = 1;
	if (isRealData)	roch->momcor_data(*thisMuon, thisMuon->Charge(), 0, qter);
	else     	roch->momcor_mc(  *thisMuon, thisMuon->Charge(), 0, qter );
      }
      muons.push_back(*thisMuon);
    }
  }

  sort(muons.begin(), muons.end(), P4SortCondition);

  hists->fill1DHist(muons.size(),     Form("size_mu_cut%i",  1),";Number of muons",    5,0,5, 1, "N");
  hists->fill1DHist(electrons.size(), Form("size_el_cut%i",  1),";Number of electrons (HZZ)",   5,0,5, 1, "N");
  hists->fill1DHist(electrons0.size(),Form("size_el0_cut%i", 1),";Number of electrons (Loose)", 5,0,5, 1, "N");
  hists->fill1DHist(DAlectrons.size(),Form("size_DALE_cut%i",1),";Number of electrons (Dalitz)",5,0,5, 1, "N");
  hists->fill1DHist(photonsHZG.size(),   Form("size_phHZG_cut%i",  1),";Number of photons (HZG)",   5,0,5, 1, "N");
  hists->fill1DHist(photonsTight.size(), Form("size_phTight_cut%i",1),";Number of photons (Tight)", 5,0,5, 1, "N");
  hists->fill1DHist(photonsMVA.size(),   Form("size_phMVA_cut%i",  1),";Number of photons (MVA)",   5,0,5, 1, "N");



  //fout<<eventNumber<<endl;

  /*
  if (electrons.size()<2) return kTRUE;

  Float_t tmpMll = 99999;
  if (electrons.size()==2){
    lPt1 = electrons[0];
    lPt2 = electrons[1];
  }
  else
    for (UInt_t m1=0; m1 < electrons.size(); m1++){
      for (UInt_t m2=m1; m2 < electrons.size(); m2++){
	if (electrons[m1].Charge()*electrons[m2].Charge()==1) continue; //same charge - don't care

	if( (electrons[m1]+electrons[m2]).M() < tmpMll)
	  {
	    lPt1 = electrons[m1];
	    lPt2 = electrons[m2];
	    tmpMll = (electrons[m1]+electrons[m2]).M();
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
    //cout<<"ch1 = "<<electrons[0].Charge()<<"  ch2 = "<<electrons[1].Charge()<<endl;
    return kTRUE;//Abort("reco * They are the same charge!");
  }

  Double_t Mll = (l1+l2).M();

  hists->fill1DHist(Mll,  Form("diLep_mass_low_cut%i", 3), ";M(ll)", 50, 0,20,  1, "");
  hists->fill1DHist(Mll,  Form("diLep_mass_high_cut%i", 3),";M(ll)", 50, 0,120, 1, "");

  */



  /*
    //if (DAlectrons.size()<1)
    return kTRUE;
  if (electrons0[0].Pt()<30)
    return kTRUE;

  if (electrons0[0].GetTracks().size()<2) return kTRUE;

  l1 = electrons0[0].GetTracks()[0];
  l2 = electrons0[0].GetTracks()[1];

  HM->SetLeptons(l1, l2);

  HM->FillHistosFull(5, eventWeight);
  FillHistoCounts(5, eventWeight);
  CountEvents(5,"",fcuts);

  //if (fabs(DAlectrons[0].Eta()) < 1.44)
  //HM->MakeElectronPlots(DAlectrons[0], "DalitzEleEB");
  //else
  //HM->MakeElectronPlots(DAlectrons[0], "DalitzEleEE");
  //cout<<" EEEEE "<<DAlectrons[0].Pt()<<endl;

  if (!electrons0[0].PassConversionVeto())
    return kTRUE;
  //else if (lPt1.DeltaR(gamma)<1.0 || lPt2.DeltaR(gamma)<1.0) return kTRUE;

  */



  /*
  for(Int_t ev=0; ev<evSize;ev++){
    if (nEvents[0]-1 == hisEVTS[ev]){cout<<nEvents[0]-1<<" evNumber="<<eventNumber<<endl;
      cout<<"gamma: pt = "<<gamma.Pt()<<"  SCEta = "<<gamma.SCEta()<<"  SCPhi = "<<gamma.SCPhi()<<endl;
      //cout<<"R9 = "<<gamma.R9()<<"  hadoverEm = "<<gamma.HadOverEm()<<endl;
      cout<<"Ele  : pt = "<<DALE.Pt()<<"   SCEta = "<<DALE.SCEta()<<"  SCPhi = "<<DALE.SCPhi()<<"  MVA = "<<DALE.IdMap("mvaScore")<<endl;
      cout<<"rho = "<<rhoFactor<<endl;
      //cout<<"Track 1: pt="<<DALE.GetTracks()[0].Pt()<<endl;
      //cout<<"Track 2: pt="<<DALE.GetTracks()[1].Pt()<<endl;
      //const vector<TCElectron::Track> trk = DALE.GetTracks();
      //cout<<"pt2/pt1="<<trk[1].Pt()/trk[0].Pt()
      //  <<"   dR = "<<trk[1].DeltaR(trk[0])<<endl;

      cout<<"gsf trk1, pt ="<<trk[0].Pt()<<"  eta= "<<trk[0].Eta()<<endl;
      cout<<"gsf trk2, pt ="<<trk[1].Pt()<<"  eta= "<<trk[1].Eta()<<endl;
      //cout<<"SCEn = "<<DALE.SCEnergy()<<"  trk1+trk2 Pt ="<<(trk[0]+trk[1]).Pt()
      //	  <<" EoverPt = "<<DALE.SCEnergy()/(trk[0]+trk[1]).Pt()<<endl;
      //cout<<"dzy = "<<DALE.Dxy(pvPosition)<<"  dz = "<<DALE.Dz(pvPosition)<<endl;
      //DALE.Dump();

      cout<<" trk1 prob = "<<xx[0]<<" lxyBS = "<<xx[1]<<" nHits = "<<xx[2]<< endl;
      cout<<" trk2 prob = "<<xx[3]<<" lxyBS = "<<xx[4]<<" nHits = "<<xx[5]<< endl;
      //const vector<TCEGamma> sc = DALE.BaseSC();
      //sc[0].Dump();
      //sc[1].Dump();
      //sc[2].Dump();
      //cout<<"sc0 "<<sc[0].Pt()<<"sc1 "<<sc[1].Pt()<<"sc2 "<<sc[2].Pt()<<endl;


      break;}
  }
  */

  //fout<<runNumber<<" "<<eventNumber<<" "<<DALE.IdMap("mvaScore")<<endl;

  Float_t MdalG = (DALE+gamma).M();
  Float_t Mllg  = (l1+l2+gamma).M();
  Float_t Mll   = (l1+l2).M();

  if (makeApzTree){
    apz_w = eventWeight;
    apz_pt1 = lPt1.Pt();
    apz_pt2 = lPt2.Pt();
    apz_pt3 = gamma.Pt();
    apz_pt4 = -1;
    apz_pt12 = DALE.Pt();
    apz_pt34 = gamma.Pt();

    apz_dr1234 = DALE.DeltaR(gamma);
    apz_dr12 = lPt1.DeltaR(lPt2);
    apz_dr13 = lPt1.DeltaR(gamma);
    apz_dr23 = lPt2.DeltaR(gamma);
    apz_dr34 = -1;

    apz_eta1 = lPt1.Eta();
    apz_eta2 = lPt2.Eta();
    apz_eta3 = gamma.SCEta();
    apz_eta4 = -999;

    apz_eta12 = DALE.Eta();
    apz_eta34 = gamma.SCEta();
    apz_eta1234 = (DALE+gamma).Eta();

    apz_m12 = (lPt1+lPt2).M();
    apz_m34 = gamma.M();

    apz_m123 = (DALE+gamma).M();
    //apz_m123 = (lPt1+lPt2+gamma).M();

    apz_m4l = 0;

    _apzTree->Fill();

  }



  if (MdalG<110 || MdalG>170) return kTRUE;
  //if (Mllg<70)
  HM->FillHistosFull(9, eventWeight);
  FillHistoCounts(9, eventWeight);
  CountEvents(9,"110 < m(e,g) < 170",fcuts);
  HM->MakeElectronPlots(DALE, "DalitzEle-cut9");


  if (electrons.size() >=2) {
    //this is second category, where two good electrons are reconstructed;
    //for now, lets just ignore it

    hists->fill1DHist(DALE.DeltaR(electrons[0]), Form("00_dR_DALE-El1_cut%i", 9),";dR(Dale, el1)", 100, 0,4, 1, "");
    hists->fill1DHist(DALE.DeltaR(electrons[1]), Form("00_dR_DALE-El2_cut%i", 9),";dR(Dale, el2)", 100, 0,4, 1, "");
    hists->fill1DHist(electrons[0].DeltaR(electrons[1]), Form("00_dR_Ele1-El2_cut%i", 9),";dR(el1, el2)", 100, 0,4, 1, "");
    //return kTRUE;
  }


  if (Mll > 20) return kTRUE;
  //if (Mll > mllMax) return kTRUE;
  HM->FillHistosFull(10, eventWeight);
  FillHistoCounts(10, eventWeight);
  CountEvents(10,"m(ll) < 20 GeV",fcuts);

  if (DALE.DeltaR(gamma)<1.0)
    //  if (MdalG<80)
    return kTRUE;

  HM->FillHistosFull(11, eventWeight);
  FillHistoCounts(11, eventWeight);
  CountEvents(11,"dR(g,e) > 1",fcuts);
  HM->MakeElectronPlots(DALE, "DalitzEle-cut11");

  //if (DALE.Pt()/Mllg < 0.30 || gamma.Pt()/Mllg < 0.30)
  if (DALE.Pt()/MdalG < 0.30 || gamma.Pt()/MdalG < 0.30)
    //if ((l1+l2).Pt() < 40 || gamma.Pt() < 40)
    return kTRUE;

  fout<<nEvents[0]-1<<endl;

  HM->FillHistosFull(12, eventWeight);
  FillHistoCounts(12, eventWeight);
  CountEvents(12,"DALE.Pt/m(e,g) > 0.3 and g.Pt/M(e,g) > 0.3",fcuts);
  HM->MakeElectronPlots(DALE, "DalitzEle-cut12");

  //fout<<"Run/lumi/event = "<<runNumber<<"/"<<lumiSection<<"/"<<eventNumber<<endl;
  //fout<<"MvaID = "<<DALE.IdMap("mvaScore")<<"  "<<DALE<<endl;
  //fout<<" SC1 = "<<DALE.BaseSC()[0]<<endl;
  //fout<<" SC2 = "<<DALE.BaseSC()[1]<<endl;


  if (!DALE.PassConversionVeto())
    return kTRUE;

  //if (DALE.ConversionMissHits()!=0)
  //return kTRUE;

  FillHistoCounts(13, eventWeight);
  CountEvents(13,"Conv gsf Ele",fcuts);
  HM->MakeElectronPlots(DALE, "DalitzEle-cut13");

  cI = trk[0].GetConversionInfo();
  if (cI.vtxProb>1e-6 && cI.lxyBS>2 && cI.nHitsMax>1)
    return kTRUE;

  HM->FillHistosFull(14, eventWeight);
  FillHistoCounts(14, eventWeight);
  CountEvents(14,"Conv gsf trk 1",fcuts);

  cI = trk[1].GetConversionInfo();
  if (cI.vtxProb>1e-6 && cI.lxyBS>2 && cI.nHitsMax>1)
    return kTRUE;

  HM->FillHistosFull(15, eventWeight);
  FillHistoCounts(15, eventWeight);
  CountEvents(15,"Conv gsf trk 2",fcuts);
  HM->MakeElectronPlots(DALE, "DalitzEle-cut15");


  if (MdalG<122 || MdalG>128) return kTRUE;
  //if (Mllg<122 || Mllg>128) return kTRUE;

  HM->FillHistosFull(16, eventWeight);
  FillHistoCounts(16, eventWeight);
  CountEvents(16,"122 < m(e,g) < 128",fcuts);
  HM->MakeElectronPlots(DALE, "DalitzEle-cut16");



  /*
  for (UInt_t i =0; i<ntrig; i++){
    triggerSelector->SelectTrigger(myTriggers[i], triggerStatus, hltPrescale, isFound, triggerPass, prescale);
    if(triggerPass) nEventsTrig[7][i]++;

    if (!isFound)
    cout<<"TRG **** Warning ***\n The trigger name "<<myTriggers[i]<<" is not in the list of trigger names"<<endl;
    }
  */


  return kTRUE;

}

void egamma::SlaveTerminate() {}

void egamma::Terminate()
{

  cout<<" ***  Terminating... **"<<endl;
  fout.close();
  fcuts.close();

  string allCuts[nC];
  string line;
  ifstream myfile("./out_cutlist.txt");
  if (myfile.is_open()){
    while (! myfile.eof() ){
      getline (myfile,line);
      Int_t n = atoi(line.substr(0,2).c_str());
      allCuts[n] = line;
      //cout << n<<"  "<<allCuts[n]<< endl;
    }
    myfile.close();
  }

  cout<<"n |"<<setw(45)<<" CUT DESCRIPTION \t|"<<" events \t"<< "eff\t|"<<endl;
  for (Int_t n=0; n<nC; n++){
    if (n==0)
      cout<<  "0 |"<<setw(45)<<allCuts[0]<<"\t |"<< nEvents[0]  <<"\t|"<<float(nEvents[0])/nEvents[0]<<"\t|"<<endl;
    else
      cout<<n<<" |"<<setw(45)<<allCuts[n]<<"\t |"<< nEvents[n]  <<"\t|"<<float(nEvents[n])/nEvents[n-1]<<"\t|"<<endl;
  }


  /*
  cout<<"\n\n | Trigger efficiency 3              |\t"<<endl;
  for(UInt_t n=0; n<ntrig; n++){
    UInt_t N=7;
    cout<<n<<"  "<<float(nEventsTrig[N][n])/nEvents[N]<<"   "<<myTriggers[n]<<endl;;
  }
  */

  hists->fill1DHist(-1, "evt_byCut",";cut #;weighted events", nC+1,-1,nC, totEvents, "Counts");
  hists->fill1DHist(-1, "evt_byCut_raw", ";cut #;events",     nC+1,-1,nC, totEvents, "Counts");

  histoFile->cd();

  histoFile->Write();
  histoFile->Close();
  cout<<" ** End of Analyzer **"<<endl;
}


void egamma::CountEvents(Int_t num, string cutName, ofstream& s)
{
  if (nEvents[num]==0)
    s<<num<<" "<<cutName<<endl;
  nEvents[num]++;
}


void egamma::FillHistoCounts(Int_t num, Double_t weight)
{
  hists->fill1DHist(num, "evt_byCut",";cut;weighted events", nC+1, -1,nC, weight, "Counts");
  hists->fill1DHist(num, "evt_byCut_raw", ";cut;events",     nC+1, -1,nC, 1, "Counts");
}
