#define zgamma_cxx
#include "zgamma.h"

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

string  myTrigger = "";
Float_t cut_l1pt  = 23;
Float_t cut_l2pt  = 7, cut_l2pt_low = 4;
const Float_t cut_iso_mu1 = 0.4, cut_iso_mu2=0;
Float_t cut_gammapt = 25;
Float_t mllMax = 25;

Bool_t checkTrigger = kTRUE;
Float_t global_Mll = 0;
//Bool_t applyPhosphor = 0;
Bool_t doRochCorr = 1;
Bool_t doZee = 0;
Bool_t makeApzTree = 1;
Bool_t accCut = 0;  //Acceptance study

const UInt_t hisEVTS[] = {9331};
Int_t evSize = sizeof(hisEVTS)/sizeof(int);

bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());}

void zgamma::Begin(TTree * tree)
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

  histoFile->cd();
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
    _apzTree->Branch("isVBF",  &apz_vbf,    "vbf/O");

    _apzTree->Branch("run",  &runNumber,  "run/i");
    _apzTree->Branch("event",&eventNumber,"event/l");
  }

  //----------------------//
  // Selecting triggers   //
  //----------------------//

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

}

void zgamma::SlaveBegin(TTree * /*tree*/)
{   TString option = GetOption(); }

Bool_t zgamma::Process(Long64_t entry)
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

  //ZGAnlgles:
  double gen_co1,gen_co2,gen_phi,gen_co3;

  Float_t genMll  = 0;
  Float_t gendR   = 0;
  Float_t genMllg = 0;


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
      if (thisParticle->GetPDGId()==23)
	hists->fill1DHist(thisParticle->M(), "gen_Z_mass",";m_{Z}^{gen} (GeV)",  200, 0,200, 1,"GEN");

      if (thisParticle->GetPDGId()==25 && (thisParticle->GetStatus()==3 || thisParticle->GetStatus()==22)){
	gen_higgs = *thisParticle;
	h++;
      }

      //GEN MUONS
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

      //if (thisParticle->GetPDGId()==22 && thisParticle->Pt()>2)
      //	ObjID->DiscoverGeneology(thisParticle);

      //PHOTON from the Higgs
      if (thisParticle->GetPDGId()==22 && (thisParticle->GetStatus()==1 || thisParticle->GetStatus()==23)
	  && thisParticle->Mother() &&  abs(thisParticle->Mother()->GetPDGId())!=11 && abs(thisParticle->Mother()->GetPDGId())!=13
	  && abs(thisParticle->Mother()->GetPDGId())!=15)
	{
	  //cout<<"  DBG gen particles. Ancestorr = "<<ObjID->GetPrimaryAncestor(thisParticle)->GetPDGId()<<endl;
	  if (ObjID->GetPrimaryAncestor(thisParticle)->GetPDGId()==A || thisParticle->Mother()->GetPDGId()==A){
	    gen_gamma = *thisParticle;
	    ph++;
	  }
	  else if (h==0 && ObjID->GetPrimaryAncestor(thisParticle)->GetPDGId() == thisParticle->GetPDGId()) {
	    gen_gamma = *thisParticle;
	    ph++;
	  }

	}
      else if (sample=="DY" &&
	       thisParticle->GetPDGId()==22 && thisParticle->GetStatus()==1 && ObjID->GetPrimaryAncestor(thisParticle)->GetPDGId()==23)
	gen_gamma = *thisParticle;



      // **** Find FSR PHOTONS *** //

      if (thisParticle->GetPDGId()==22 //&& thisParticle->GetStatus()==1
	  && thisParticle->Mother()
	  //&& ObjID->GetPrimaryAncestor(thisParticle)->GetPDGId()==A
	  && selection=="mugamma" && abs(thisParticle->Mother()->GetPDGId())==13
	  )
	{
	  hists->fill1DHist(thisParticle->Pt(), "gen_fsr_pho_pt",";gen fsr photon pt",  200, 0,30, 1,"GEN");
	  // if (ObjID->GetPrimaryAncestor(thisParticle)->GetPDGId()!=A)

	  //if (thisParticle->Pt() > 30)
	  //ObjID->DiscoverGeneology(thisParticle);

	  Float_t preFSR_pt = thisParticle->Mother()->Pt();
	  for (int j = 0; j < genParticles->GetSize(); ++j) {
	    TCGenParticle* fsr = (TCGenParticle*) genParticles->At(j);
	    if(fsr->GetPDGId()  == thisParticle->Mother()->GetPDGId()
	       && fsr->Mother() == thisParticle->Mother()){

	      //cout<<"   0--- this one --->>>"<<endl;
	      //ObjID->DiscoverGeneology(fsr);
	      //cout<<"   <<---0 this one"<<endl;

	      hists->fill1DHist((preFSR_pt - fsr->Pt())/preFSR_pt, "gen_fsr_mu_dPt",";gen mu (Pt_preFSR - Pt_afterFSR)/pt",  200, -0.4,1, 1,"GEN");

	      if(thisParticle->Pt() > 1)
		hists->fill1DHist((preFSR_pt - fsr->Pt())/preFSR_pt, "gen_fsr_mu_dPt_photonPt1",";gen mu (Pt_preFSR - Pt_afterFSR)/pt",  200, -0.4,1, 1,"GEN");
	      fsr_mu_count++;
	    }
	  }
	}//end of fsr photons loop
    }// end of genparticles loop


    if (ph!=1 && (sample=="dalitz" || sample =="HZG"))
      Abort(Form("ev #%i NONONO There has to be exactly one photon from the Higgs! \n \t\t but there is %i",
		 (int)eventNumber, (int)ph));
    //if (h!=1 && sample=="dalitz")
    //return kTRUE;
    //Abort(Form(" NONONO There has to be exactly one Higgs! \n \t\t but there is %i", h));

    sort(gen_el.begin(), gen_el.end(), P4SortCondition);
    sort(gen_mu.begin(), gen_mu.end(), P4SortCondition);

    hists->fill1DHist(fsr_mu_count, "gen_mu_fsrcount",";gen mu FSR",  4, 0,4, 1,"GEN");
    hists->fill1DHist(fsr_el_count, "gen_el_fsrcount",";gen el FSR",  4, 0,4, 1,"GEN");


    if(selection=="mugamma"){
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


    if (sample=="dalitz" || sample=="HZG"){
      ang->GetAngles(gen_l1, gen_l2, gen_gamma, gen_co1,gen_co2,gen_phi,gen_co3);
      //cout<<eventNumber<<" gen Angles: c1= "<<co1<<"  c2="<<co2<<"   gen_phi="<<gen_phi<<"   coTh="<<co3<<endl;
      hists->fill1DHist(gen_co1, "gen_co1",";gen cos_lp",  100,-1,1, 1,"GEN");
      hists->fill1DHist(gen_co2, "gen_co2",";gen cos_lm",  100,-1,1, 1,"GEN");
      hists->fill1DHist(gen_co3, "gen_co3",";gen cosTheta",100,-1,1, 1,"GEN");
      hists->fill1DHist(gen_phi, "gen_phi",";gen phi lp",  100, -TMath::Pi(), TMath::Pi(), 1,"GEN");
    }

    FillHistoCounts(1, eventWeight);
    CountEvents(1,"Pass Gen acceptance",fcuts);

    gendR   = gen_l1.DeltaR(gen_l2);
    genMll  = (gen_l1+gen_l2).M();
    genMllg = (gen_l1+gen_l2+gen_gamma).M();

    hists->fill1DHist(genMll,"gen_Mll_low",  ";m_{#mu#mu}",100,0, 50, 1,"GEN");
    hists->fill1DHist(genMll,"gen_Mll_full", ";m_{#mu#mu}",200,0,130, 1,"GEN");
    hists->fill1DHist(genMll,"gen_Mll_inter",";m_{#mu#mu}",100,20,80, 1,"GEN");

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
      else if (accCut)
	return kTRUE;
    }
    else if (accCut)
      return kTRUE;

    if (sample == "dalitz" || sample=="HZG"){
      hists->fill1DHist(gen_l1.DeltaR(gen_gamma),"gen_l1gamma_deltaR",";#Delta R(#mu_{1}, #gamma)",100,0,5, 1,"GEN");
      hists->fill1DHist(gen_l2.DeltaR(gen_gamma),"gen_l2gamma_deltaR",";#Delta R(#mu_{2}, #gamma)",100,0,5, 1,"GEN");

      hists->fill1DHist(gen_higgs.Pt(),  "gen_higgs_pt",  ";p_{T}^{H} (GeV)",   100, 0,150,  1, "GEN");
      hists->fill1DHist(gen_higgs.M(),   "gen_higgs_mass",";m_{H} (GeV)",    100, 50,180, 1, "GEN");
      hists->fill1DHist(gen_gamma.Pt(),  "gen_gamma_pt",  ";p_{T}^{#gamma}",  50, 0,150,  1, "GEN");
      hists->fill1DHist(gen_gamma.Eta(), "gen_gamma_eta", ";#eta^{#gamma}",   50, -3,3,   1, "GEN");

      //hists->fill1DHist(gen_gamma_st1.Pt(), "gen_gamma_st1_pt","Pt of gamma",   50, 0,150,  1, "");
      //hists->fill1DHist(gen_gamma_st1.Eta(),"gen_gamma_st1_eta","Eta of gamma", 50, -3,3,  1, "");
    }

  }

  FillHistoCounts(2, eventWeight);
  CountEvents(2, "Acceptance", fcuts);

  vector<TCElectron> electrons0, electrons;
  vector<TCMuon> muons0, muons;
  vector<TCPhoton> photons0, photonsTight, photonsHZG, photonsMVA, fake_photons;
  vector<TCJet> jets;
  TCPhysObject l1,l2, lPt1, lPt2;
  TCPhoton gamma0, gamma;
  TCJet jet1, jet2;

  //------------------------//
  //--- PHOTON SELECTION ---//
  //-----------------------//

  for (Int_t i = 0; i < recoPhotons->GetSize(); ++i) {
    //cout<<"new photon!!!!!!!"<<endl;
    TCPhoton* thisPhoton = (TCPhoton*) recoPhotons->At(i);
    //cout<<"in photon loop event = "<<eventNumber
    //  <<"\n pt="<<thisPhoton->Pt()<<" eta="<<thisPhoton->Eta()<<" phi="<<thisPhoton->Phi()<<" px="<<thisPhoton->Px()<<endl;

    //if (isRealData || !makeGen || (makeGen &&gen_lPt1.DeltaR(*thisPhoton) < 0.1 && thisPhoton->Pt()>20))
    //fake_photons.push_back(*thisPhoton);

    if (!(fabs(thisPhoton->SCEta()) < 2.5)) continue;


    if(thisPhoton->Pt() > 15){
      //if (isRealData || !makeGen || (sample=="dalitz" && makeGen && gen_gamma.DeltaR(*thisPhoton) < 0.1) )
      photons0.push_back(*thisPhoton);

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
	  //Double_t scale = 1;
	  Double_t scale = corrPhoEnReco/thisPhoton->E(); //Yes - use this one for mumugamma
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

  sort(photons0.begin(),     photons0.end(),     P4SortCondition);
  sort(photonsTight.begin(), photonsTight.end(), P4SortCondition);
  sort(photonsHZG.begin(),   photonsHZG.end(),   P4SortCondition);
  sort(photonsMVA.begin(),   photonsMVA.end(),   P4SortCondition);
  //sort(fake_photons.begin(), fake_photons.end(), P4SortCondition);

  //if (fake_photons.size()>1)
  // Abort("Nah, this is too much.");

  //if (photonsTight.size()<1) return kTRUE;
  //gamma = photonsTight[0];

  if (photonsHZG.size()<1) return kTRUE;
  gamma = photonsHZG[0];

  //if (photonsMVA.size()<1) return kTRUE;
  //gamma = photonsMVA[0];

  if (!isRealData)
    eventWeight *= weighter->PhotonSF(gamma);

  //if (photons0.size()<1) return kTRUE;
  //gamma0 = photons0[0];
  //if (gamma0.Pt() < cut_gammapt) return kTRUE;

  //------------------------//
  //--- ELECTRON SELECTION ---//
  //-----------------------//
  for (Int_t i = 0; i <  recoElectrons->GetSize(); ++i) {
    TCElectron* thisElec = (TCElectron*) recoElectrons->At(i);
    if (!(fabs(thisElec->Eta()) < 2.5)) continue;
    if (thisElec->Pt() > 10.0){
      if(isRealData || !makeGen || ( selection=="elgamma" && makeGen && (gen_l1.DeltaR(*thisElec) < 0.1 || gen_l2.DeltaR(*thisElec)<0.1)  && thisElec->Pt()>10))
	{
	  if (thisElec->MvaID() > -0.50 && thisElec->Pt()>20)
	    electrons.push_back(*thisElec);
	}
      if (
	  //PassElectronIdAndIsoMVA(thisElec)
	  ObjID->PassElectronIdAndIso(*thisElec, pvPosition, "Loose")
	  )
	electrons0.push_back(*thisElec);

    }
  }
  sort(electrons0.begin(), electrons0.end(), P4SortCondition);
  sort(electrons.begin(),  electrons.end(),  P4SortCondition);

  if (doZee){
    //----------------------//
    //------ Z->ee STUDY ---//
    //----------------------//

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
	if (photons0[p].Pt() < 40 || fabs(photons0[p].SCEta())>1.4442) continue;
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


  //------------------------//
  //--- MUON SELECTION ---//
  //-----------------------//

  for (Int_t i = 0; i < recoMuons->GetSize(); ++ i) {
    TCMuon* thisMuon = (TCMuon*) recoMuons->At(i);
    //cout<<"event = "<<eventNumber
    //  <<"\n pt="<<thisMuon->Pt()<<" eta="<<thisMuon->Eta()<<" phi="<<thisMuon->Phi()<<" px="<<thisMuon->Px()<<endl;

    if (!(fabs(thisMuon->Eta()) < 2.4)) continue;

    if(thisMuon->Pt() > 4)
      for(Int_t ev=0; ev<evSize;ev++){
	if (eventNumber==hisEVTS[ev]){cout<<"   mu>4 cut ---> "<<eventNumber<<" <--- Found an event after cut "<<endl;
	  ObjID->MuonDump(*thisMuon, pvPosition);
	  break;}}

    if(thisMuon->Pt() > 4 && ObjID->PassMuonId(*thisMuon, pvPosition, "Soft") ) {
      if (doRochCorr){
	Float_t qter = 1;
	//cout<<"DBG before rochcor"<<endl;
	if (isRealData)	  roch->momcor_data(*thisMuon, thisMuon->Charge(), 0, qter);
	else         	  roch->momcor_mc(*thisMuon,  thisMuon->Charge(), 0, qter );

	//cout<<"DBG after rochcor"<<endl;
      }

      muons.push_back(*thisMuon);
    }
  }

  sort(muons.begin(), muons.end(), P4SortCondition);


  if(!isRealData && makeGen){
    hists->fill1DHist(genMll,  "gen_Mll_reco_gamma",  ";gen_Mll",100,0,mllMax, 1,"eff");
    hists->fill1DHist(gendR,   "gen_dR_reco_gamma",   ";gen_dR", 50, 0,0.3,    1,"eff");
  }


  if (checkTrigger){
    triggerSelector->SelectTrigger(myTrigger, triggerStatus, hltPrescale, isFound, triggerPass, prescale);
    if (!triggerPass) return kTRUE;
  }

  if (gamma.Pt() < cut_gammapt) return kTRUE;
  //if (fabs(gamma.SCEta()) > 1.4442) return kTRUE;

  if(!isRealData && makeGen){
    hists->fill1DHist(genMll,  "gen_Mll_reco_gamma_iso",  ";gen_Mll",100,0,mllMax, 1,"eff");
    hists->fill1DHist(gendR,   "gen_dR_reco_gamma_iso",   ";gen_dR", 50, 0,0.3,    1,"eff");
  }

  FillHistoCounts(3, eventWeight);
  CountEvents(3,"Pass Trigger and Photon reco selection",fcuts);
  HM->MakeNPlots(3,muons.size(), electrons.size(), electrons0.size(), photonsHZG.size(),
		 photonsTight.size(), photonsMVA.size(), jets.size(), eventWeight);

  if (muons.size()<2) return kTRUE;

  //------------------------//
  //------ PICK THE Muons --//
  //------------------------//
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

  //------------------------//
  //------JET SELECTION ---//
  //-----------------------//

  for (Int_t i = 0; i < recoJets->GetSize(); ++i) {
    TCJet* thisJ = (TCJet*) recoJets->At(i);
    if (thisJ->Pt() < 30 || fabs(thisJ->Eta())>4.7) continue;
    if (thisJ->DeltaR(gamma) < 0.4) continue;
    if (thisJ->DeltaR(lPt1) < 0.3 || thisJ->DeltaR(lPt2) < 0.3) continue;

    if (ObjID->PassJetID(*thisJ, nVtx))
      jets.push_back(*thisJ);
  }

  sort(jets.begin(), jets.end(), P4SortCondition);

  if (jets.size()>=1){
    jet1 = jets[0];
    HM->SetJet1(jet1);
  }

  Bool_t isVBF = 0;

  if (jets.size()>=2){
    jet2 = jets[1];
    HM->SetJet2(jet2);

    //--------------------------//
    //----- VBF SELECTION ------//
    //--------------------------//
    if (fabs(jet1.Eta() - jet2.Eta()) > 3.
	&& (jet1+jet2).M() > 450
	&& fabs(HM->Zeppenfeld((lPt1+lPt2+gamma),jet1,jet2)) < 4
	//&& (jet1+jet2).DeltaPhi(lPt1+lPt2+gamma)) > 0.1
	)
    isVBF = 1;
  }

  //-------------------------//
  //------MISSING ENERGY ---//
  //------------------------//

  HM->SetMet(*recoMET);
  //cout<<"MET: "<<recoMET->Mod()<<endl;

  //---------------------------------//
  //------ KINEMATIC SELECTION ------//
  //---------------------------------//

  Double_t Mll  = (l1+l2).M();
  Double_t Mllg = (l1+l2+gamma).M();
  Double_t dR   = l1.DeltaR(l2);

  hists->fill1DHist(Mll,  Form("diLep_mass_low_cut%i",  3),";M(ll)", 50, 0,20,  1, "");
  hists->fill1DHist(Mll,  Form("diLep_mass_high_cut%i", 3),";M(ll)", 50, 0,120, 1, "");


  if (ObjID->CalculateMuonIso(muons[0]) > 0.4)
    return kTRUE;
  //if (lPt2.Pt() > 20 && ObjID->CalculateMuonIso(&muons[1]) > 0.4)
  //return kTRUE;

  HM->MakePhotonPlots(gamma);

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
    hists->fill1DHist(gendR,   "gen_dR_two_lep_reco",  ";gen_dR", 50, 0,0.3,    1,"eff");


    double co1,co2,phi,co3;

    if (sample=="dalitz" || sample=="HZG"){
      ang->GetAngles(l1, l2, gamma, co1,co2,phi,co3);
      //cout<<eventNumber<<" gen Angles: c1= "<<co1<<"  c2="<<co2<<"   phi="<<phi<<"   coTh="<<co3<<endl;
      hists->fill1DHist((gen_co1 - co1)/gen_co1, "res_co1",";Resolution on cos(#vec{l-},#vec{Z})",  100, -0.2,0.2, 1,"GEN");
      hists->fill1DHist((gen_co2 - co2)/gen_co2, "res_co2",";Resolution on cos(#vec{l+},#vec{Z})",  100, -0.2,0.2, 1,"GEN");
      hists->fill1DHist((gen_co3 - co3)/gen_co3, "res_co3",";Resolution on cos(#Theta)",  100, -0.2,0.2, 1,"GEN");
      hists->fill1DHist((gen_phi - phi)/gen_phi, "res_phi",";Resolution on #phi(l+)",     100, -0.2,0.2, 1,"GEN");
    }


    if (fabs(gamma.SCEta()) < 1.4442){
      hists->fill1DHist((genMll - Mll)/genMll, "res_EB_Mll",
			";(m_{#mu#mu}^{gen} - m_{#mu#mu}^{reco})/m_{#mu#mu}^{gen}",                  100, -0.2,0.2, 1,"GEN");
      hists->fill1DHist((genMllg - Mllg)/genMllg, "res_EB_Mllg",
			";(m_{#mu#mu#gamma}^{gen} - m_{#mu#mu#gamma}^{reco})/m_{#mu#mu#gamma}^{gen}",100, -0.2,0.2, 1,"GEN");
      hists->fill1DHist((gendR - dR)/gendR,    "res_EB_dR",
			";(dR_{#mu#mu}^{gen} - dR_{#mu#mu}^{reco})/dR_{#mu#mu}^{gen}",               100, -0.2,0.2, 1,"GEN");

      float mllBins[9] = {0,0.2,0.5,1,2,4,9,20,50};
      hists->fill2DHistUnevenBins(genMll,Mll , "res_2D_EB_Mll", ";m_{#mu#mu}^{gen};m_{#mu#mu}^{reco}", 8,mllBins, 8,mllBins, 1,"GEN");

      for(int m=1; m<8; m++){
	if (Mll > mllBins[m-1] && Mll < mllBins[m]){

	  hists->fill1DHist((genMll - Mll)/genMll, Form("res_EB_Mll_%.1fto%.1f",mllBins[m-1],mllBins[m]),
			    ";(m_{#mu#mu}^{gen} - m_{#mu#mu}^{reco})/m_{#mu#mu}^{gen}", 100, -0.1,0.1, 1,"GEN");
	  hists->fill1DHist((genMllg - Mllg)/genMllg, Form("res_EB_Mllg_%.1fto%.1f",mllBins[m-1],mllBins[m]),
			    ";(m_{#mu#mu#gamma}^{gen} - m_{#mu#mu#gamma}^{reco})/m_{#mu#mu#gamma}^{gen}",100, -0.1,0.1, 1,"GEN");
	}
      }
    }
    else if  (fabs(gamma.SCEta()) > 1.566){
      hists->fill1DHist((genMll - Mll)/genMll, "res_EE_Mll",
			";(m_{#mu#mu}^{gen} - m_{#mu#mu}^{reco})/m_{#mu#mu}^{gen}",                  100, -0.2,0.2, 1,"GEN");
      hists->fill1DHist((genMllg - Mllg)/genMllg, "res_EE_Mllg",
			";(m_{#mu#mu#gamma}^{gen} - m_{#mu#mu#gamma}^{reco})/m_{#mu#mu#gamma}^{gen}",100, -0.2,0.2, 1,"GEN");
      hists->fill1DHist((gendR - dR)/gendR,    "res_EE_dR",
			";(dR_{#mu#mu}^{gen} - dR_{#mu#mu}^{reco})/dR_{#mu#mu}^{gen}",               100, -0.2,0.2, 1,"GEN");


    }
  }

  HM->SetLeptons(lPt1, lPt2);
  HM->SetGamma(gamma);


  HM->FillHistosFull(4, eventWeight);
  FillHistoCounts(4, eventWeight);
  CountEvents(4, "Two muons selected", fcuts);
  HM->MakeNPlots(4,muons.size(), electrons.size(), electrons0.size(), photonsHZG.size(),
		 photonsTight.size(), photonsMVA.size(), jets.size(), eventWeight);


  if (Mll > 50)
    return kTRUE;

  if (makeApzTree){
    apz_w = eventWeight;
    apz_pt1 = lPt1.Pt();
    apz_pt2 = lPt2.Pt();
    apz_pt3 = gamma.Pt();
    apz_pt4 = -1;
    apz_pt12 = (lPt1+lPt2).Pt();
    apz_pt34 = gamma.Pt();

    apz_dr1234 = (lPt1+lPt2).DeltaR(gamma);
    apz_dr12 = lPt1.DeltaR(lPt2);
    apz_dr13 = lPt1.DeltaR(gamma);
    apz_dr23 = lPt2.DeltaR(gamma);
    apz_dr34 = -1;

    apz_eta1 = lPt1.Eta();
    apz_eta2 = lPt2.Eta();
    apz_eta3 = gamma.SCEta();
    apz_eta4 = -999;

    apz_eta12 = (lPt1+lPt2).Eta();
    apz_eta34 = gamma.SCEta();
    apz_eta1234 = (lPt1+lPt2+gamma).Eta();

    apz_m12 = (lPt1+lPt2).M();
    apz_m34 = gamma.M();

    apz_m123 = (lPt1+lPt2+gamma).M();

    apz_m4l = 0;
    apz_vbf = isVBF;

    _apzTree->Fill();

  }

  //for(Int_t ev=0; ev<evSize;ev++){
  //if (eventNumber==hisEVTS[ev]){cout<<eventNumber<<" Found an event after lept selection "<<endl;
  //  cout<<"Mll = "<<Mll<<endl;
  //  break;}
  //}

  if (Mllg<110 || Mllg>170)
    return kTRUE;

  HM->FillHistosFull(5, eventWeight);
  FillHistoCounts(5, eventWeight);
  CountEvents(5, "110 < m(llg) < 170; m(ll) < 50", fcuts);

  if (lPt1.DeltaR(gamma)<1.0 || lPt2.DeltaR(gamma)<1.0) return kTRUE;
  if ( (Mll>2.9 && Mll<3.3) || (Mll>9.3 && Mll<9.7)) return kTRUE; //jpsi and upsilon removeal

  HM->FillHistosFull(6, eventWeight);
  FillHistoCounts(6, eventWeight);
  CountEvents(6, "dR(l,g) > 1; removed J/Psi, Ups", fcuts);
  HM->MakeNPlots(6,muons.size(), electrons.size(), electrons0.size(), photonsHZG.size(),
		 photonsTight.size(), photonsMVA.size(), jets.size(), eventWeight);

  if (Mll < 20) {
    HM->FillHistosFull(7, eventWeight);
    FillHistoCounts(7, eventWeight);
    CountEvents(7, "m(ll) < 20 GeV", fcuts);
    HM->MakeNPlots(7,muons.size(), electrons.size(), electrons0.size(), photonsHZG.size(),
		   photonsTight.size(), photonsMVA.size(), jets.size(), eventWeight);

    HM->MakeMuonPlots(muons[0]);
    HM->MakeMuonPlots(muons[1]);

    if (fabs(gamma.SCEta()) < 1.4442){

      HM->FillHistosFull(8, eventWeight);
      FillHistoCounts(8, eventWeight);
      CountEvents(8, "gamma eta < 1.4442", fcuts);

      if ((l1+l2).Pt()/Mllg > 0.30 && gamma.Pt()/Mllg > 0.30) {
	HM->FillHistosFull(9, eventWeight);
	FillHistoCounts(9, eventWeight);
	CountEvents(9, "pT/m(llg) > 0.3", fcuts);
	HM->MakeNPlots(9,muons.size(), electrons.size(), electrons0.size(), photonsHZG.size(),
		       photonsTight.size(), photonsMVA.size(), jets.size(), eventWeight);

	HM->MakePhotonPlots(gamma, "Pho-after");
	HM->MakeMuonPlots(muons[0],"Mu-after");
	HM->MakeMuonPlots(muons[1],"Mu-after");

	if (Mllg>122 && Mllg<128){
	  HM->FillHistosFull(10, eventWeight);
	  FillHistoCounts(10, eventWeight);
	  CountEvents(10, "122 GeV < m(llg) < 128 GeV", fcuts);
	}

	if (isVBF){
	  HM->FillHistosFull(11, eventWeight);
	  FillHistoCounts(11, eventWeight);
	  CountEvents(11, "VBF ID", fcuts);
	}

	if ((lPt1.Pt()+lPt2.Pt()) > 40){
	  HM->FillHistosFull(12, eventWeight);
	  FillHistoCounts(12, eventWeight);
	  CountEvents(12, "ptl1+ptl2 > 40 GeV", fcuts);
	}


      }
    }
    else if (fabs(gamma.SCEta()) > 1.566) {
      HM->FillHistosFull(15, eventWeight);
      FillHistoCounts(15, eventWeight);
      CountEvents(15, "gamma eta > 1.566", fcuts);

      if ((l1+l2).Pt()/Mllg > 0.30 && gamma.Pt()/Mllg > 0.30) {
	HM->FillHistosFull(16, eventWeight);
	FillHistoCounts(16, eventWeight);
	CountEvents(16, "pT/m(llg) > 0.3", fcuts);

	if (Mllg>122 && Mllg<128){
	  HM->FillHistosFull(17, eventWeight);
	  FillHistoCounts(17, eventWeight);
	  CountEvents(17, "122 GeV < m(llg) < 128 GeV", fcuts);
	}
      }
    }
  }
  else if (Mll < 50 && fabs(gamma.SCEta()) < 1.4442){

    HM->FillHistosFull(18, eventWeight);
    FillHistoCounts(18, eventWeight);
    CountEvents(18, "20 GeV < m(ll) < 50 GeV; EB only", fcuts);

    if ((l1+l2).Pt()/Mllg > 0.30 && gamma.Pt()/Mllg > 0.30) {
      HM->FillHistosFull(19, eventWeight);
      FillHistoCounts(19, eventWeight);
      CountEvents(19, "pT/m(llg) > 0.3", fcuts);

      if (Mllg>122 && Mllg<128){
	HM->FillHistosFull(20, eventWeight);
	FillHistoCounts(20, eventWeight);
	CountEvents(20, "122 GeV < m(llg) < 128 GeV", fcuts);
      }
    }
  }


  for (UInt_t i=0; i<ntrig; i++){
    triggerSelector->SelectTrigger(myTriggers[i], triggerStatus, hltPrescale, isFound, triggerPass, prescale);
    if(triggerPass) nEventsTrig[6][i]++;

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




  if(!isRealData && makeGen){
    for (Int_t i = 1; i<=6; i++){
      if ((l1+l2).Pt()>20+i*5)
	hists->fill1DHist(genMll,  Form("gen_Mll_%i",i),";gen_Mll", 100,0,mllMax,    1,"eff");
      if (gamma.Pt()>20+i*5)
	hists->fill1DHist(genMll,  Form("gen_Mll_%i",6+i),";gen_Mll", 100,0,mllMax,  1,"eff");
      if (gamma.Pt()/Mllg> 0.15+0.05*i)
	hists->fill1DHist(genMll,  Form("gen_Mll_%i",12+i),";gen_Mll", 100,0,mllMax, 1,"eff");
      if ((l1+l2).Pt()/Mllg> 0.15+0.05*i)
	hists->fill1DHist(genMll,  Form("gen_Mll_%i",18+i),";gen_Mll", 100,0,mllMax, 1,"eff");
      if ((l1+l2).Pt()/Mllg> 0.15+0.05*i   && gamma.Pt()/Mllg> 0.15+0.05*i)
	hists->fill1DHist(genMll,  Form("gen_Mll_%i",24+i),";gen_Mll", 100,0,mllMax, 1,"eff");
    }
  }




  //if (checkTrigger){
  //triggerSelector->SelectTrigger(myTrigger, triggerStatus, hltPrescale, isFound, triggerPass, prescale);
  //if (!triggerPass) return kTRUE;
  ////  else Abort("Event trigger should mu or el");
  //}

  if (NoiseFilters_isScraping ||
      NoiseFilters_isNoiseHcalHBHE || NoiseFilters_isNoiseHcalLaser
      )
    return kTRUE;

  if (NoiseFilters_isNoiseEcalTP || NoiseFilters_isNoiseEcalBE || NoiseFilters_isNoiseEEBadSc
      )
    return kTRUE;

  if (NoiseFilters_isCSCTightHalo// || NoiseFilters_isCSCLooseHalo
      )
    return kTRUE;

  if (NoiseFilters_isNoiseTracking ||
      !NoiseFilters_isNoisetrkPOG1 || !NoiseFilters_isNoisetrkPOG2 || !NoiseFilters_isNoisetrkPOG3
      )
    return kTRUE;

  FillHistoCounts(21, eventWeight);
  CountEvents(21, "noise filters, after cuts #6",fcuts);


  /*
  if (lPt1.Phi() > -1.8 && lPt1.Phi() < -1.6){
    HM->FillHistosFull(18, eventWeight);
    fout<<" nEvt = "<<nEvents[0]<<" : Run/lumi/event = "<<runNumber<<"/"<<lumiSection<<"/"<<eventNumber<<endl;
    fout<<"lPt1:   "<<&lPt1<<endl;
    fout<<"dilep phi = "<<(lPt1+lPt2).Phi()<<endl;
  }
  */

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

  cout<<" ** FOR PAS **"<<endl;
  cout<<"n |"<<setw(45)<<" CUT DESCRIPTION \t|"<<" events \t"<< "eff\t|"<<endl;
  for (Int_t n=0; n<nC; n++){
    if (n==0)
      cout<<  "0 |"<<setw(45)<<allCuts[0]<<"\t |"<< nEvents[0]  <<"\t|"<<float(nEvents[0])/nEvents[0]<<"\t|"<<endl;
    else
      cout<<n<<" |"<<setw(45)<<allCuts[n]<<"\t |"<< nEvents[n]  <<"\t|"<<float(nEvents[n])/nEvents[n-1]<<"\t|"<<endl;
  }


  cout<<"\n\n | Trigger efficiency 3              |\t"<<endl;
  for(UInt_t n=0; n<ntrig; n++){
    UInt_t N=6;
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


void zgamma::CountEvents(Int_t num, string cutName, ofstream& s)
{
  if (nEvents[num]==0)
    s<<num<<" "<<cutName<<endl;
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


void zgamma::FillHistoCounts(Int_t num, Double_t weight)
{
  hists->fill1DHist(num, "evt_byCut",";cut;weighted events", nC+1, -1,nC, weight, "Counts");
  hists->fill1DHist(num, "evt_byCut_raw", ";cut;events",     nC+1, -1,nC, 1, "Counts");
}
