#define analyzer_higgs_cxx

#include "analyzer_higgs.h"

using namespace std;

Float_t cut_vz = 24, cut_vd0 = 2, cut_vndof = 4;  //PV filter cuts
Float_t cut_Mll_low = 76.2, cut_Mll_high = 106.2, cut_pTll = 25.; //Di-lepton quantities
 
UInt_t nEvents[20], nEventsPassHcal[20], nEventsPassEcal[20];
UInt_t nEventsPassNoiseFilter[6][20]; //1-Hcal, 2-Ecal, 3- Scraping, 4-CSCTight, 5-CSCLoose,   0-passed all

Bool_t isElectronSample = kFALSE, isMuonSample = kTRUE;

UInt_t verboseLvl = 1;

string sel, sample;
Float_t CS, nEv;  Int_t lColor, fColor;

Float_t ct_pfMet, ct_projMet, ct_MT, ct_qT; //variables for cutTree

//Variables to Fill into histos. Have to be global
Float_t qT, diEta, Mll, Mll_EB, Mll_EE, Mll_EX;
Float_t MET, MET_phi, MET1, MET1_phi, projMET, puCorrMET, puSigMET;
Float_t METqt, MET1qt, projMETqt, puCorrMETqt;
Float_t MT, MTZ; 
Int_t nJets, nJetsB, jetCount;
Float_t dPhiClos1,dPhiClos2, dPhiLead1, dPhiLead2;
 

void analyzer_higgs::Begin(TTree * /*tree*/)
{ 
  TString option = GetOption();
  TH1::SetDefaultSumw2(kTRUE);

  histoFile = new TFile("hhhh.root", "RECREATE"); 
  //histoFile = new TFile(Form("hhistograms_%s.root", sample.Data()), "RECREATE"); 
  for(Int_t n=0; n<nC; n++)
    {
      met0_phi[n]      = new TH1F(Form("met0_phi_%i",n), "met0_phi", 40, -TMath::Pi(), TMath::Pi()); 
      met0_et[n]       = new TH1F(Form("met0_et_%i",n), "met0_et", 40, 0,400);
      met0_over_qt[n]  = new TH1F(Form("met0_over_qt_%i",n), "met0_over_qt", 40, 0,4);

      met0_et_ovQt[n]  = new TH2F(Form("met0_et_ovQt_%i",n), "pfMET vs Met/qt", 40, 0,400, 40, 0,4);
      met1_et_ovQt[n]  = new TH2F(Form("met1_et_ovQt_%i",n), "MET1 vs Met/qt", 40, 0,400, 40, 0,4);
      met2_et_ovQt[n]  = new TH2F(Form("met2_et_ovQt_%i",n), "pfMet noise vs Met/qt", 40, 0,400, 40, 0,4);
      met3_et_ovQt[n]  = new TH2F(Form("met3_et_ovQt_%i",n), "projMET vs Met/qt", 40, 0,400, 40, 0,4);
      met4_et_ovQt[n]  = new TH2F(Form("met4_et_ovQt_%i",n), "puCorrMET vs Met/qt", 40, 0,400, 40, 0,4);
    
      met1_phi[n]     = new TH1F(Form("met1_phi_%i",n), "met1_phi", 40, -TMath::Pi(), TMath::Pi()); 
      met1_et[n]      = new TH1F(Form("met1_et_%i",n), "met1_et", 40, 0,400);
      met1_over_qt[n] = new TH1F(Form("met1_over_qt_%i",n), "met1_over_qt", 40, 0,4);

      met2_phi[n]     = new TH1F(Form("met2_phi_%i",n), "met2_phi", 40, -TMath::Pi(), TMath::Pi()); 
      met2_et[n]      = new TH1F(Form("met2_et_%i",n), "met2_et", 40, 0,400);
      met2_over_qt[n] = new TH1F(Form("met2_over_qt_%i",n), "met2_over_qt", 40, 0,4);

      met3_phi[n]     = new TH1F(Form("met3_phi_%i",n), "met3_phi", 40, -TMath::Pi(), TMath::Pi()); 
      met3_et[n]      = new TH1F(Form("met3_et_%i",n), "met3_et", 40, 0,400);
      met3_over_qt[n] = new TH1F(Form("met3_over_qt_%i",n), "met2_over_qt", 40, 0,4);

      met4_phi[n]     = new TH1F(Form("met4_phi_%i",n), "met4_phi", 40, -TMath::Pi(), TMath::Pi()); 
      met4_et[n]      = new TH1F(Form("met4_et_%i",n), "met4_et", 42, -20,400);
      met4_over_qt[n] = new TH1F(Form("met4_over_qt_%i",n), "met4_over_qt", 42, -0.2,4);
      met4_puSig[n]   = new TH1F(Form("met4_sig_%i",n), "met4_sig", 42, -0.4,8);

      met2_dPhiLeadJet1[n] =  new TH1F(Form("met2_dPhiLeadJet1_%i",n), "delta phi to lead jets", 40, 0, TMath::Pi());
      met2_dPhiLeadJet2[n] =  new TH1F(Form("met2_dPhiLeadJet2_%i",n), "delta phi to lead jets", 40, 0, TMath::Pi());
      met2_dPhiClosJet1[n] =  new TH1F(Form("met2_dPhiClosJet1_%i",n), "delta phi to closest jets", 40, 0, TMath::Pi());
      met2_dPhiClosJet2[n] =  new TH1F(Form("met2_dPhiClosJet2_%i",n), "delta phi to closest jets", 40, 0, TMath::Pi());


      mu1_phi[n] = new TH1F(Form("mu1_phi_%i",n), "mu1_phi", 40, -TMath::Pi(), TMath::Pi()); 
      mu1_eta[n] = new TH1F(Form("mu1_eta_%i",n), "mu1_eta", 40, -2.6, 2.6); 
      mu1_pt[n]  = new TH1F(Form("mu1_pt_%i",n), "mu1_pt", 40, 0, 200); 

      mu2_phi[n] = new TH1F(Form("mu2_phi_%i",n), "mu2_phi", 40, -TMath::Pi(), TMath::Pi()); 
      mu2_eta[n] = new TH1F(Form("mu2_eta_%i",n), "mu2_eta", 40, -2.6, 2.6); 
      mu2_pt[n]  = new TH1F(Form("mu2_pt_%i",n), "mu2_pt", 40, 0, 200); 
      btag_hp[n] = new TH1F(Form("btag_hp_%i",n), "btag_hp", 40, -10, 10); 

      mtZ[n]     = new TH1F(Form("mtZ_%i",n), "MTZ pf", 50, 100, 600); 
      mt0[n]     = new TH1F(Form("mt0_%i",n), "MT pf", 50, 100, 600); 
      mt1[n]     = new TH1F(Form("mt1_%i",n), "MT type1", 50, 100, 600); 
      mt2[n]     = new TH1F(Form("mt2_%i",n), "MT noise", 50, 100, 600); 
      mt3[n]     = new TH1F(Form("mt3_%i",n), "MT proj", 50, 100, 600); 
      mt4[n]     = new TH1F(Form("mt4_%i",n), "MT pu corr", 50, 100, 600);

      mtZ_met2[n]   = new TH2F(Form("mtZ_met2_%i",n), "MTZ pf vs pfMet noise", 50, 100, 600, 40, 0,400); 
      mt2_met2[n]   = new TH2F(Form("mt2_met2_%i",n), "MT pf vs pfMet noise", 50, 100, 600, 40, 0,400); 

      mtZ_met3[n]   = new TH2F(Form("mtZ_met3_%i",n), "MTZ pf vs projMet noise", 50, 100, 600, 40, 0,400); 
      mt2_met3[n]   = new TH2F(Form("mt2_met3_%i",n), "MT pf vs projfMet noise", 50, 100, 600, 40, 0,400); 
 
      di_qt[n]      = new TH1F(Form("di_qt_%i",n), "di-lepton qt", 40, 0,400);  
      di_eta[n]     = new TH1F(Form("di_eta_%i",n), "di-lepton Eta", 40, -2.6,2.6);  
      di_mass[n]    = new TH1F(Form("di_mass_%i",n), "di-lepton Mass", 50, 70,120);  
      di_mass_EB[n] = new TH1F(Form("di_mass_EB_%i",n), "di-lepton Mass in Ecal Barrel", 50, 70,120);  
      di_mass_EE[n] = new TH1F(Form("di_mass_EE_%i",n), "di-lepton Mass in Ecal Endcap", 50, 70,120);  
      di_mass_EX[n] = new TH1F(Form("di_mass_EX_%i",n), "di-lepton Mass in Ecal E/B mix", 50, 70,120);  

      jet_N[n]      = new TH1F(Form("jet_N_%i",n), "Number of jets", 20,0,20);
      jet_dRlep1[n] = new TH1F(Form("jet_dRlep1_%i",n), "dR(jet, lepton1)", 50, 0,1);
      jet_dRlep2[n] = new TH1F(Form("jet_dRlep2_%i",n), "dR(jet, lepton2)", 50, 0,1);
      jet_b_N[n]    = new TH1F(Form("jet_b_N_%i",n), "Number of b-jets", 10,0,10);

      jet_pt[n]   = new TH1F(Form("jet_pt_%i",n), "Pt of jets, eta<2.4", 40,0,200);
      jet_b_pt[n] = new TH1F(Form("jet_b_pt_%i",n), "Pt of b-jets, eta<2.4", 40,0,200);

    }
  
  cutTree = new TTree("cutTree","Tree for cuts");
  cutTree->Branch("ct_pfMet",&ct_pfMet, "ct_pfMet/F");
  cutTree->Branch("ct_projMet",&ct_projMet, "ct_projMet/F");
  cutTree->Branch("ct_MT",&ct_MT, "ct_MT/F");
  cutTree->Branch("ct_qT",&ct_qT, "ct_qT/F");
  
  ffout.open("./events_printout_andrey_final.txt",ofstream::out);
  ncout.open("./counts_for_tex.txt",ofstream::out);
  ifstream params;
  params.open("./params.txt"); 
  params>>sel;
  if (sel == "muon")      { ncout<<"Mu selection!"<<endl; isMuonSample = kTRUE;  isElectronSample = kFALSE;}
  if (sel == "electron")  { ncout<<"El selection!"<<endl; isMuonSample = kFALSE; isElectronSample = kTRUE;}
  ncout<<endl;
  params>>sample;  params>>nEv;  params>>CS;  params>>lColor;  params>>fColor;
  ncout<<"sample: "<<sample<<"  cs: "<<CS<<"  events: "<<nEv<<" colors: "<<lColor<<"  "<<fColor<<endl;
  ncout<<endl;
  
  ffout.precision(4); ffout.setf(ios::fixed, ios::floatfield);
  for(Int_t i=0; i<nC; i++)
    {
      fout[i].open(Form("./events_printout_andrey_%i.txt",i),ofstream::out);
      fout[i].precision(3); fout[i].setf(ios::fixed, ios::floatfield);
      fout[i]<<
	"******\t*********\t*****\t********\t************\t*********\t********\t*****\t*\n"<<
	"* Run \t* Event  \t* LS \t*  m(ll)\t*   MT(llvv)\t* PFMET  \t* PT(ll)\t* rho\t*\n"<<
	"******\t*********\t*****\t********\t************\t*********\t********\t*****\t*\n";
      nout[i].open(Form("./events_filtered_AllFilters_%i.txt",i),ofstream::out);
      nout[i].precision(3); fout[i].setf(ios::fixed, ios::floatfield);
      nout[i]<<
	"******\t*********\t*********\t******************************************\n"<<
	"* Run \t* Event  \t*  LS   \t*  PFMET * Hcal * Ecal * scrap * CSCTight * \n"<<
	"******\t*********\t*********\t******************************************\n";
    }

}

void analyzer_higgs::SlaveBegin(TTree * /*tree*/)
{ TString option = GetOption(); }

Bool_t analyzer_higgs::Process(Long64_t entry)
{
  GetEntry(entry); 
  //if ( nEvents[0]>250) return kTRUE; 
  CountEvents(0);


  //cout<<"run: "<runNumber<<" \t event: "<<eventNumber<<" \t lumi: "<<lumiSection<<" \t bx: "<<bunchCross<<endl;
  //cout<<"rho:  "<<rhoFactor<<endl;

  //---------------------------------------------------
  //--- Primary vertex filter -------------------------
  //-------------------------------------------------
  
  Int_t PVind[2] = {-1, -1}; //indecies for two best PVs
  Bool_t vertexFilter = kFALSE;
  Int_t nVtx = 0;
  for (Int_t vv=0; vv<primaryVtxDAWithBS->GetSize(); vv++)
    {
      TCPrimaryVtx* pVtx = (TCPrimaryVtx*) primaryVtxDAWithBS->At(vv);
      //if(!pVtx->isValid()) cout<<"PV is not valid"<<endl;
      //if(pVtx->isFake()) cout<<"PV is fake"<<endl;

      if( !pVtx->isFake() && 
	  pVtx->NDof() > cut_vndof   &&                    //4
	  fabs(pVtx->Position().z()) <= cut_vz &&     //15
	  fabs(pVtx->Position().Perp()) <= cut_vd0   //Not Mag()!; 2
	  )
	{
	  vertexFilter = kTRUE;

	  if (PVind[0]==-1) PVind[0] = vv; //Set first main PV
	  if (PVind[1]==-1 && PVind[0]!=vv) PVind[1] = vv; //Second PV
	  nVtx++;
	}
      //if (pVtx->NDof()==0) PVndof0++;
      
      //      if (!pVtx->isValid()) 

      //h1_PV_chi2OverNdof -> Fill(pVtx->Chi2()/pVtx->NDof());
      //h1_PV_Ndof -> Fill(pVtx->NDof());
      //h1_PV_d0 -> Fill(pVtx->Position().Perp());
    }

  TCPrimaryVtx *mainPrimaryVertex = 0, *secondPrimaryVertex = 0;
  if (PVind[0]!=-1) mainPrimaryVertex   = (TCPrimaryVtx*)(primaryVtxDAWithBS->At(PVind[0]));
  if (PVind[1]!=-1) secondPrimaryVertex = (TCPrimaryVtx*)(primaryVtxDAWithBS->At(PVind[1]));
  
  // Apply the PV filters here! -------------
  if (!vertexFilter) return kTRUE;  

  Bool_t passPreSelection = kFALSE, passZpeak = kFALSE;
  // Bool_t passMuon = kFALSE, passElectron = kFALSE;


  //Lepton variables:
  //  Float_t 
  Mll=0, Mll_EE=0, Mll_EB=0, Mll_EX=0;//, pTll=0;
  //--------------//
  //--- Muons------//
  //-------------//
  TLorentzVector diLepton(0.,0.,0.,0.), Lepton1(0.,0.,0.,0.), Lepton2(0.,0.,0.,0.);
 
  TCMuon *myMuon1=0, *myMuon2=0, *myMuon3=0, *myMuon4=0;
  Int_t muIndex[5] = {-1,-1,-1,-1,-1};

  Int_t muCount=0, muCountLoose = 0;
  for (Int_t i = 0; i < recoMuons->GetSize(); ++i) 
    {    
      TCMuon* thisMuon   = (TCMuon*) recoMuons->At(i);     
      TVector3 *muVertex =  new TVector3(thisMuon->Vtx());
      TVector3 *muMomentum = new TVector3( thisMuon->p4().Vect());
      
      double dxyVtx = dxy(mainPrimaryVertex, muVertex, muMomentum);
      double dzVtx  = dz(mainPrimaryVertex, muVertex, muMomentum);
      
      if (fabs(thisMuon->p4().Eta()) < 2.4 && 
	  thisMuon->p4().Pt()>10     &&
	  thisMuon->isGLB()  && 
	  fabs(dxyVtx)<0.02  && 
	  fabs(dzVtx)<0.1    &&
	  thisMuon->normChi2()     < 10 &&
	  thisMuon->nTRKHits()     > 10 &&
	  thisMuon->nPXLHits()     > 0  &&
	  thisMuon->nValidMuHits() > 0  &&
	  thisMuon->nMatchSeg()    > 1  &&
	  (thisMuon->trkIso() + thisMuon->hadIso() + thisMuon->emIso() - rhoFactor*TMath::Pi()*0.09)/thisMuon->pt() < 0.15
	  
	  ) 
	{
	  //cout<<i<<" pt: "<<thisMuon->p4().Pt() <<" eta: "<<thisMuon->p4().Eta()<<endl;
	  
	  if (muIndex[1]==-1) muIndex[1] = i;                                                                  
	  if (muIndex[2]==-1 && muIndex[1]!=i)   muIndex[2] = i;
	  if (muIndex[3]==-1 && muIndex[1]!=i && muIndex[2]!=i)   muIndex[3] = i;
	  if (muIndex[4]==-1 && muIndex[1]!=i && muIndex[2]!=i && muIndex[3]!=i) muIndex[4] = i;
	  
	  muCount++;
	}

      if (fabs(thisMuon->p4().Eta()) < 2.4 &&
          thisMuon->p4().Pt()>10     &&
          thisMuon->isGLB()  &&
	  fabs(dxyVtx)<0.05  && 
	  fabs(dzVtx)<0.2    &&
	  thisMuon->normChi2()     < 200 &&
	  thisMuon->nTRKHits()     > 2   &&
	  //thisMuon->nPXLHits()     > 0   &&
	  thisMuon->nValidMuHits() > 0   
	  // thisMuon->nMatchSeg()    > 0   
	  )  muCountLoose++;
      
    }
  
  if (muCount>=2) 
    {
      myMuon1 = (TCMuon*)recoMuons->At(muIndex[1]); 
      myMuon2 = (TCMuon*)recoMuons->At(muIndex[2]);
      //swap the muons to have first one has larger pt
      if (myMuon1->pt() < myMuon2->pt()) 
	{
	  myMuon1 = (TCMuon*)recoMuons->At(muIndex[2]); 
	  myMuon2 = (TCMuon*)recoMuons->At(muIndex[1]);
	}
      
      if (muCount>2) myMuon3 = (TCMuon*)recoMuons->At(muIndex[3]); 
      if (muCount>3) myMuon4 = (TCMuon*)recoMuons->At(muIndex[4]); 
      //if (muCount>2) cout<<"more then 2 Muons!  "<<muCount<<endl;
      
      Lepton1 = (TLorentzVector)(myMuon1->p4());
      Lepton2 = (TLorentzVector)(myMuon2->p4());
      //if (myMuon1->charge()*myMuon2->charge()!=-1) cout<<"Not opposite sign: "<<myMuon1->charge()<<"  "<<myMuon2->charge()<<endl;
      //cout<<"Muons indeces: "<<muIndex[1]<<"  "<<muIndex[2]<<endl;
      //cout<< muIndex[1]<<"  "<<myMuon1->p4().Pt()<<" "<<myMuon1->p4().Eta()<<endl;
      //cout<< muIndex[2]<<"  "<<myMuon2->p4().Pt()<<" "<<myMuon2->p4().Eta()<<endl;
	
    }
  

  //-------------------//
  //--- Electrons------//
  //-------------------//
  TCElectron *myElectron1=0, *myElectron2=0;

  Int_t eleIndex[3] = {-1,-1,-1};

  //if (recoElectrons->GetSize()>1)  cout<<"N ele: "<<  recoElectrons->GetSize()<<endl;
  Int_t eleCount=0, eleCountLoose=0;
  for (Int_t i = 0; i < recoElectrons->GetSize(); ++i) 
    {    
      TCElectron* thisElectron = (TCElectron*) recoElectrons->At(i);     
      TVector3* eleVertex = new TVector3(thisElectron->Vtx());
      TVector3* eleMomentum = new TVector3(thisElectron->p4().Vect());
      
      double dxyVtx = dxy(mainPrimaryVertex, eleVertex, eleMomentum);
      double dzVtx  = dz(mainPrimaryVertex, eleVertex, eleMomentum);
      
      Float_t  eleISO = (thisElectron->trkIso() + thisElectron->emIso() + thisElectron->hadIso() - rhoFactor*TMath::Pi()*0.09)/thisElectron->pt(); 
      Float_t  eleISObarrel = 0;
      eleISObarrel = (thisElectron->trkIso() + TMath::Max(0.0, thisElectron->emIso()-1.0)+ thisElectron->hadIso() - rhoFactor*TMath::Pi()*0.09)/thisElectron->pt(); 
      
      if (fabs(thisElectron->eta())   < 2.5 && 
	  thisElectron->p4().Pt()     > 10  && 
	  fabs(dxyVtx) < 0.02               && 
	  fabs(dzVtx)  < 0.1                &&
	  (
	   (fabs(thisElectron->eta()) < 1.4442       && 
	    eleISObarrel < 0.15                      &&
	    thisElectron->sigmaIetaIeta() < 0.01     &&
	    fabs(thisElectron->dPhiSuperCluster()) < 0.06  &&
	    fabs(thisElectron->dEtaSuperCluster()) < 0.004 &&
	    thisElectron->hadOverEm() < 0.04      
	    )  ||
	   (fabs(thisElectron->eta()) >  1.566       && 
	    eleISO < 0.15                            &&
	    thisElectron->sigmaIetaIeta()    < 0.03  &&
	    fabs(thisElectron->dPhiSuperCluster()) < 0.03  &&
	    fabs(thisElectron->dEtaSuperCluster()) < 0.007 &&
	    thisElectron->hadOverEm() < 0.025
	    )
	   )
	  ) 
	{
	  //cout<<i<<" pt: "<<thisElectron->p4().Pt() <<" eta: "<<thisElectron->p4().Eta()<<endl;
	  if (eleIndex[1]==-1) eleIndex[1] = i;                                                                  
	  if (eleIndex[2]==-1 && eleIndex[1]!=i) eleIndex[2] = i;
	  
	  eleCount++;
	}
	  
      if (fabs(thisElectron->eta())   < 2.5 && 
	  thisElectron->p4().Pt()     > 10  && 
	  fabs(dxyVtx) < 0.05               && 
	  fabs(dzVtx)  < 0.2                &&
	  thisElectron->sigmaIetaIeta()          < 0.03 &&
	  fabs(thisElectron->dPhiSuperCluster()) < 0.8  &&
	  fabs(thisElectron->dEtaSuperCluster()) < 0.01 &&
	  thisElectron->hadOverEm() < 0.15
	  )  eleCountLoose++;
      
    }
  if (eleCount==2) 
    {
      myElectron1 = (TCElectron*) recoElectrons->At(eleIndex[1]); 
      myElectron2 = (TCElectron*) recoElectrons->At(eleIndex[2]); 
      if (myElectron1->pt() < myElectron2->pt())
	{
	  myElectron1 = (TCElectron*) recoElectrons->At(eleIndex[2]); 
	  myElectron2 = (TCElectron*) recoElectrons->At(eleIndex[1]); 
	}
      Lepton1 = (TLorentzVector)(myElectron1->p4());
      Lepton2 = (TLorentzVector)(myElectron2->p4());

      //cout<<"Electrons indeces: "<<eleIndex[1]<<"  "<<eleIndex[2]<<endl;
      //cout<< eleIndex[1]<<"  "<<myElectron1->p4().Pt()<<" "<<myElectron1->p4().Eta()<<endl;
      //cout<< eleIndex[2]<<"  "<<myElectron2->p4().Pt()<<" "<<myElectron2->p4().Eta()<<endl;
      
    }

  Float_t pTll =0;

  if (isMuonSample && muCount==2)
    {
      Mll  = (myMuon1->p4()+myMuon2->p4()).M();
      pTll = (myMuon1->p4()+myMuon2->p4()).Pt();
      diLepton = (TLorentzVector)(myMuon1->p4()+myMuon2->p4());
      //if((diLepton.M()-Mll)>0.01 || (diLepton.Pt()-pTll)>0.01) cout<<" Warning!  Dilepton is wrong"
      
      //if (eleCount==0 && myMuon1->charge()*myMuon2->charge()==-1 && Mll>40 && Mll<500)	 
      if (eleCountLoose==0 && muCountLoose==2 && myMuon1->charge()*myMuon2->charge()==-1 && Mll>40 && Mll<500)	 
	{
	  passPreSelection= kTRUE;
	  if (Mll>cut_Mll_low && Mll<cut_Mll_high) passZpeak = kTRUE;
	}
    }
  if (isElectronSample && eleCount==2)
    {
      Mll  = (myElectron1->p4()+myElectron2->p4()).M();
      pTll = (myElectron1->p4()+myElectron2->p4()).Pt();
     
      diLepton = (TLorentzVector)(myElectron1->p4()+myElectron2->p4());
      //if (muCount==0 && myElectron1->charge()*myElectron2->charge()==-1 && Mll>40 && Mll<500)	 
      if (muCountLoose==0 && eleCountLoose==2 && myElectron1->charge()*myElectron2->charge()==-1 && Mll>40 && Mll<500)	 
	{
	  //cout<<"preselection?!"<<endl;
	  passPreSelection= kTRUE;
	  if (Mll>cut_Mll_low && Mll<cut_Mll_high) passZpeak = kTRUE;
	}

    }

  if (recoMET->GetSize()>1) cout<<"warning: to many mets"<<endl;
  
  TCMET* pfMet = 0, *pfMetCorr = 0;

  if (recoMET->GetSize()==1)
    {
      pfMet = (TCMET*) recoMET->At(0);     
      MET = pfMet->Met();
      MET_phi = pfMet->Phi();
    }

  if (corrMET->GetSize()==1)
    {
      pfMetCorr = (TCMET*) corrMET->At(0);     
      MET1 = pfMetCorr->Met();
      MET1_phi = pfMetCorr->Phi();
    }


  qT = diLepton.Pt();

  Double_t weight = 1;
  if(sample=="MC_ZZtoAny") weight *= weightZZ(qT);


  projMET     = projectedMET(MET1, MET1_phi, Lepton1, Lepton2);
  projMETqt   = projMET/qT;
  puCorrMET   = pileUpCorrectedMET(projMET, nVtx);
  puCorrMETqt = puCorrMET/qT;
  Float_t sigma0      = 7.343;
  puSigMET    = puCorrMET/sigma0;
  
  //cout<<"pfMET: "<<MET<<"  pfType1: "<<MET1<<"  proj1: "<<projMET<<"  puCorrMET: "<< puCorrMET<<"   signif: "<<puSigMET<<"  nVtx: "<<nVtx<<endl;
  TVector2 ZptXY(0.0,0.0), METXY(0.0,0.0);  
  ZptXY.Set(diLepton.X(), diLepton.Y());
  METXY.Set(MET*cos(MET_phi), MET*sin(MET_phi));
  MT  = higgsMT(pTll, Mll, MET, ZptXY, METXY);
  MTZ = higgsMT(pTll, 91.1876, MET, ZptXY, METXY);
  //if (muCount==2) cout<<"MT: "<<MT<<"    "<<ZptXY.X()<<" "<<ZptXY.Y()<<"    "<<METXY.X()<<" "<<METXY.Y()<<endl;

  METqt  = MET/qT;
  MET1qt = MET1/qT;


   FillHistosNoise(0, weight);
   if(!isNoiseHcal && !isDeadEcalCluster  && !isScraping && !isCSCTightHalo) FillHistos(0, weight);
  else
    if(verboseLvl>1) 
      nout[0]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<MET<<"\t* "
	     <<isNoiseHcal<<"\t* "<<isDeadEcalCluster<<"\t* "<<isScraping<<"\t* "<<isCSCTightHalo<<"\t* "<<endl;  
  
  if (!passPreSelection) return kTRUE;
  
  // Only make sense after preselection (two leptons):
  if(fabs(Lepton1.Eta()) < 1.444 && fabs(Lepton2.Eta())<1.444) Mll_EB = diLepton.M();
  else  if(fabs(Lepton1.Eta()) > 1.566 && fabs(Lepton2.Eta())>1.566)  Mll_EE = diLepton.M();
  else  Mll_EX = diLepton.M();

  diEta = diLepton.Eta();

  Bool_t antiB  = kFALSE;
  jetCount=0, nJets = 0, nJetsB = 0;

  Float_t ptJ1 = 0, ptJ2=0;
  dPhiLead1 = 110,  dPhiLead2 = 110; 
  dPhiClos1 = 110,  dPhiClos2 = 110; 
  for (Int_t i = 0; i < recoJets->GetSize(); ++i) 
    {    
      TCJet* thisJet = (TCJet*) recoJets->At(i);     
      if(thisJet->NumConstit()    > 1
	 && thisJet->NumChPart()  > 0 
	 && thisJet->ChHadFrac()  > 0
	 && thisJet->ChEmFrac()   < 1
	 && (thisJet->NeuHadFrac() + thisJet->NeuEmFrac()) < 0.9
	 && thisJet->P4().DeltaR(Lepton1) > 0.4
	 && thisJet->P4().DeltaR(Lepton2) > 0.4
	 ){
	
	if (thisJet->P4(7).Pt() > 30 && fabs(thisJet->P4(7).Eta()) <= 2.4) 
	  {
	    //do anti-b tag
	    if (thisJet->BDiscrTrkCountHiEff() > 2)
	      {
		antiB = kTRUE;
		nJetsB++;
	      }
	    nJets++;
	  }
	
	Float_t ptTemp = thisJet->P4(7).Pt();
	if(ptTemp>ptJ1)
	  {
	    ptJ1 = ptTemp;
	    dPhiLead1 = fabs(TVector2::Phi_mpi_pi(thisJet->P4(7).Phi() - MET_phi));
	  }
	
	if(ptTemp>ptJ2 && fabs(thisJet->P4(7).Eta()) <= 2.4)
	  {
	    ptJ2 = ptTemp;
	    dPhiLead2 = fabs(TVector2::Phi_mpi_pi(thisJet->P4(7).Phi() - MET_phi));
	  }
	
	Float_t dphiTemp = fabs(TVector2::Phi_mpi_pi(thisJet->P4(7).Phi() - MET_phi));
	
	if(dphiTemp<dPhiClos1)
	  dPhiClos1 = dphiTemp;
	
	if(dphiTemp<dPhiClos2 && fabs(thisJet->P4(7).Eta()) <= 2.4)
	  dPhiClos2 = dphiTemp;
	
	
	jetCount++;
      }
    }

  for(Int_t ii=0; ii<nC; ii++)
    if(ii!=8 && ii!=4) jet_b_pt[ii] -> Fill(66, weight);



  FillHistosNoise(1, weight);
  CountEvents(1);

  if(!isNoiseHcal && !isDeadEcalCluster  && !isScraping && !isCSCTightHalo) 
    FillHistos(1, weight);
  else
    if(verboseLvl>0) 
      nout[1]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<MET<<"\t* "
	     <<isNoiseHcal<<"\t* "<<isDeadEcalCluster<<"\t* "<<isScraping<<"\t* "<<isCSCTightHalo<<"\t* "<<endl;  
  
  if (!passZpeak) return kTRUE;
  CountEvents(2);
  FillHistosNoise(2, weight);

  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo) 
    FillHistos(2, weight);
     
  if (isMuonSample)
    {
      if (myMuon1->pt() < 20) return kTRUE;
      CountEvents(3);

      if (myMuon2->pt() < 20) return kTRUE;
      CountEvents(4);
    }

  if (isElectronSample)
    {
      if (myElectron1->pt() < 20) return kTRUE;
      CountEvents(3);
  
      if (myElectron2->pt() < 20) return kTRUE;
      CountEvents(4);
    }
     

  FillHistosNoise(4, weight);
  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo)
    FillHistos(4, weight);  

  if(verboseLvl>1) 
    fout[4]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<Mll<<"\t* "<<MT<<"\t* "<<MET<<"\t* "<<pTll<<"\t* "<<rhoFactor<<"\t*"<<endl;  


  for (Int_t i = 0; i < recoJets->GetSize(); ++i) 
    {    
      TCJet* thisJet = (TCJet*) recoJets->At(i);     
      if(thisJet->NumConstit()    > 1
	 && thisJet->NumChPart()  > 0 
	 && thisJet->ChHadFrac()  > 0
	 && thisJet->ChEmFrac()   < 1
	 && (thisJet->NeuHadFrac() + thisJet->NeuEmFrac()) < 0.9
	 ){
	
	if (fabs(thisJet->P4(7).Eta()) <= 2.4 && thisJet->BDiscrTrkCountHiEff() > 2) 
	  jet_b_pt[4] ->Fill(thisJet->P4(7).Pt(), weight);

	jet_dRlep1[4] -> Fill(thisJet->P4(7).DeltaR(Lepton1), weight);
	jet_dRlep2[4] -> Fill(thisJet->P4(7).DeltaR(Lepton2), weight);

      }
  }
  
  if (projMET>70)
    {
      jet_b_N[9]   -> Fill(nJetsB, weight);
      if(MT>327 && MT<441)
	jet_b_N[10]   -> Fill(nJetsB, weight);

     for (Int_t i = 0; i < recoJets->GetSize(); ++i) 
	{    
	  TCJet* thisJet2 = (TCJet*) recoJets->At(i);     
	  if (fabs(thisJet2->P4(7).Eta()) <= 2.4 && thisJet2->BDiscrTrkCountHiEff() > 2) 
	    jet_b_pt[10] ->Fill(thisJet2->P4(7).Pt(), weight);
	}
    }


  if(pTll<cut_pTll) return kTRUE;
  if (antiB) return kTRUE;
  if(jetCount>0  && dPhiClos1 < 0.3) return kTRUE;

  CountEvents(5);
  FillHistosNoise(5, weight);

  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo) 
    {
      FillHistos(5, weight);

      if(projMET<70) {
	//Reversing Met cut
	FillHistos(11, weight);
      }
    }


  if(verboseLvl>1)   
    fout[5]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<Mll<<"\t* "<<MT<<"\t* "<<MET<<"\t* "<<pTll<<"\t* "<<rhoFactor<<"\t*"<<muCountLoose<<"\t* "<<eleCountLoose<<endl;  

  // if (pTll<50) return kTRUE;
 
  if (projMET<50) return kTRUE;
 
  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo)
    {
      ct_pfMet = MET; ct_projMet=projMET; ct_MT = MT; ct_qT = qT;
      cutTree -> Fill();
    }


  CountEvents(6);
  FillHistosNoise(6, weight);
  
  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo) 
    FillHistos(6, weight);
  else
    if(verboseLvl>0) 
      nout[6]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<MET<<"\t* "
	     <<isNoiseHcal<<"\t* "<<isDeadEcalCluster<<"\t* "<<isScraping<<"\t* "<<isCSCTightHalo<<"\t* "<<endl;  


  if(verboseLvl>0) 
    fout[6]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<Mll<<"\t* "<<MT<<"\t* "<<MET<<"\t* "<<pTll<<"\t* "<<rhoFactor<<"\t* "<<muCountLoose<<"\t* "<<eleCountLoose<<endl;  

  //  cout<<"dbg 322"<<endl;
  if (projMET<112) return kTRUE;
  if (MT<292) return kTRUE;
  
  //  if(projMETqt <0.6 || projMETqt>1.8) return kTRUE;
  
  FillHistosNoise(7, weight);
  CountEvents(7);
  
  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo) 
    FillHistos(7, weight);
  else
    if(verboseLvl>0) 
      nout[7]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<MET<<"\t* "
	     <<isNoiseHcal<<"\t* "<<isDeadEcalCluster<<"\t* "<<isScraping<<"\t* "<<isCSCTightHalo<<"\t* "<<endl;  
  

    fout[7]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<Mll<<"\t* "<<MT<<"\t* "<<MET<<"\t* "<<pTll<<"\t* "<<rhoFactor<<"\t* "<<muCountLoose<<"\t* "<<eleCountLoose<<endl;  

    CountEvents(8);
    CountEvents(9);
    CountEvents(10);

  //  if(verboseLvl>1)   fout[9]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<Mll<<"\t* "<<pTll<<"\t* "<<MET<<"\t* "<<pTll<<"\t* "<<rhoFactor<<"\t*"<<endl;  


    if(isMuonSample && verboseLvl>0)
    { 
    ffout<<"## Run: "<<runNumber<<"  event: "<<eventNumber<<"  lumiSection: "<<lumiSection<<" ##"<<endl;  
    ffout<<"Vertex: x="<<mainPrimaryVertex->Position().x()<<" y="<<mainPrimaryVertex->Position().y()<<" z="<<mainPrimaryVertex->Position().z()<<
    "    d0: "<<mainPrimaryVertex->Position().Perp()<<"  ndof: "<<mainPrimaryVertex->NDof()<<endl;
    ffout<<"   Muon1. q="<<myMuon1->charge()<<" pt: "<<myMuon1->p4().Pt()<<"   eta: "<<myMuon1->eta()<<"   phi: "<<myMuon1->phi()<<endl;
    ffout<<"   Muon2. q="<<myMuon2->charge()<<" pt: "<<myMuon2->p4().Pt()<<"   eta: "<<myMuon2->eta()<<"   phi: "<<myMuon2->phi()<<endl;
    ffout<<"   M(mumu): "<<Mll<<"  pT(mumu): "<<pTll<<"      MET: "<<MET<<endl;
    }


  return kTRUE;
}

void analyzer_higgs::SlaveTerminate(){}

void analyzer_higgs::Terminate()
{ 
  cout<<" Begin Terminate"<<endl;
  cout<<endl;

  for(Int_t i=0; i<nC; i++) {nout[i].close();  fout[i].close();}
  ffout.close();

  ncout.precision(3); cout.setf(ios::fixed, ios::floatfield);

  ncout<<"0 Total:  \n| "<<nEvents[0]<<" |\t"<<(Float_t)nEvents[0]/nEvents[0]<<" |"<<(Float_t)nEvents[0]/nEvents[0]<<" |"<<endl;
  ncout<<"1 presel       \n| "<<nEvents[1]<<" |\t"<<(Float_t)nEvents[1]/nEvents[0]<<" |"<<(Float_t)nEvents[1]/nEvents[0]<<" |"<<endl;
  ncout<<"2 Mll   \n| "<<nEvents[2]<<" |\t"<<(Float_t)nEvents[2]/nEvents[0]<<" |"<<(Float_t)nEvents[2]/nEvents[1]<<" |"<<endl;
  ncout<<"3 pt1     \n| "<<nEvents[3]<<" |\t"<<(Float_t)nEvents[3]/nEvents[0]<<" |"<<(Float_t)nEvents[3]/nEvents[2]<<" |"<<endl;
  ncout<<"4 pt2   \n| "<<nEvents[4]<<" |\t"<<(Float_t)nEvents[4]/nEvents[0]<<" |"<<(Float_t)nEvents[4]/nEvents[3]<<" |"<<endl;
  ncout<<"5 anti-b  \n| "<<nEvents[5]<<" |\t"<<(Float_t)nEvents[5]/nEvents[0]<<" |"<<(Float_t)nEvents[5]/nEvents[4]<<" |"<<endl;
  ncout<<"6 Met    \n| "<<nEvents[6]<<" |\t"<<(Float_t)nEvents[6]/nEvents[0]<<" |"<<(Float_t)nEvents[6]/nEvents[5]<<" |"<<endl;
  //ncout<<"7        \n| "<<nEvents[7]<<" |\t"<<(Float_t)nEvents[7]/nEvents[0]<<" |"<<(Float_t)nEvents[7]/nEvents[6]<<" |"<<endl;
  //ncout<<"8        \n| "<<nEvents[8]<<" |\t"<<(Float_t)nEvents[8]/nEvents[0]<<" |"<<(Float_t)nEvents[8]/nEvents[7]<<" |"<<endl;
  // ncout<<"9      \n| "<<nEvents[9]<<" |\t"<<(Float_t)nEvents[9]/nEvents[0]<<" |"<<(Float_t)nEvents[9]/nEvents[8]<<" |"<<endl;
  // ncout<<"10     \n| "<<nEvents[10]<<" |\t"<<(Float_t)nEvents[10]/nEvents[0]<<" |"<<(Float_t)nEvents[10]/nEvents[9]<<" |"<<endl;
  
  ncout<<endl;
  ncout<<" & total   &  passHcal  & passEcal  &  Scraping  & CSCTight & CSCLoose & & passAll & frac \\\\ \n \\hline"<<endl;
  for (Int_t i=0; i<nC; i++)
    {
      ncout<<" & "<<nEvents[i]<<
	" & "<<nEventsPassNoiseFilter[1][i]<<"\t& "<<nEventsPassNoiseFilter[2][i]<<"\t&  "<<nEventsPassNoiseFilter[3][i]<<
	" & "<<nEventsPassNoiseFilter[4][i]<<"\t& "<<nEventsPassNoiseFilter[5][i]<<"\t& &"<<nEventsPassNoiseFilter[0][i];
      if(nEvents[i]!=0)  ncout<<" & "<<(Float_t)nEventsPassNoiseFilter[0][i]/(Float_t)nEvents[i];
      else               ncout<<" & n/a ";
      ncout<<"   \\\\ \n \\hline"<<endl;
    }
  ncout.close();

  cout<<"End of Job. Let's write the files."<<endl;

  Float_t nEv1 = getNevents(sample.c_str(), met1_et[0]);

  histoFile->Write();
  scaleAndColor(sample.c_str(), CS, nEv1, 1000.0, lColor, fColor); //1fb

  histoFile->Close();

  cout<<"Wrote and Closed"<<endl;

}


double analyzer_higgs::dxy(TCPrimaryVtx *primVtx, TVector3* trackVtx, TVector3* trackMom)
{
  //Calculating track dxy parameter wrt primary vertex
  //d0 = - dxy
  Float_t vx = trackVtx->x(), vy = trackVtx->y();
  Float_t px = trackMom->Px(), py = trackMom->Py(), pt = trackMom->Pt();
  Float_t pvx = primVtx->Position().x(), pvy = primVtx->Position().y();
  
  Double_t ret =  (-(vx-pvx)*py + (vy-pvy)*px)/pt;
  return ret;
}


double analyzer_higgs::dz(TCPrimaryVtx *primVtx, TVector3* trackVtx, TVector3* trackMom)
{
  //Calculating track dz parameter wrt primary vertex
  Float_t vx = trackVtx->x(), vy = trackVtx->y(), vz = trackVtx->z();
  Float_t px = trackMom->Px(), py = trackMom->Py();
  Float_t pz = trackMom->Pz(), pt = trackMom->Pt();
  Float_t pvx = primVtx->Position().x(), pvy = primVtx->Position().y(), pvz = primVtx->Position().z();
  
  Float_t ret =  (vz-pvz)-((vx-pvx)*px +(vy-pvy)*py)/pt*(pz/pt);
  return ret;
}

double analyzer_higgs::higgsMT(Float_t f_pTll, Float_t f_Mll, Float_t f_MET, TVector2 f_ZptXY, TVector2 f_METXY)
{
  double f_MT = sqrt( pow(sqrt(f_pTll*f_pTll + f_Mll*f_Mll) + sqrt(f_MET*f_MET + f_Mll*f_Mll), 2) - (f_ZptXY + f_METXY).Mod2());
  return f_MT;
}

double analyzer_higgs::projectedMET(Float_t f_MET, Float_t f_MET_phi, TLorentzVector f_Lepton1, TLorentzVector f_Lepton2)
{
  Float_t deltaPhiMin =   TMath::Min( fabs(TVector2::Phi_mpi_pi(f_Lepton1.Phi() - f_MET_phi)), fabs(TVector2::Phi_mpi_pi(f_Lepton2.Phi() - f_MET_phi)));
  //cout<<deltaPhiMin<<endl;
  Float_t f_projMET = 0;
  if (deltaPhiMin > TMath::Pi()/2) f_projMET = f_MET;
  else f_projMET = f_MET*sin(deltaPhiMin);
  return f_projMET;
}

double analyzer_higgs::pileUpCorrectedMET(Float_t f_projMET, Int_t nVtx)
{
  //See talk by Marionneau Matthieu at DibosonMeeting on May 11, 2011
  Float_t sigma0  = 7.343;//, sigma0_err = 0.111;
  Float_t sigmaPU = 3.9;
  Float_t a = 0.856; //slope
  Float_t sigmaT = sqrt(sigma0*sigma0 + (nVtx-1)*sigmaPU*sigmaPU);
  Float_t projMET2 = f_projMET - a*(nVtx-1);
  Float_t f_puCorrMET = projMET2*sigma0/sigmaT;
  return f_puCorrMET;
}



//void analyzer_higgs::scale(string sample1, TFile *file, TFile * outFile, Float_t cs, Float_t nTotEv, Float_t lumi)
void analyzer_higgs::scaleAndColor(TString sample1, Float_t cs, Float_t nTotEv, Float_t lumi, Int_t line, Int_t fill)
{
  Float_t scaleFactor      = lumi*cs/nTotEv;
  Bool_t isData = kFALSE; 
  if(sample1.Contains("2011A_")) {scaleFactor = 1.0; isData = kTRUE;}
  if(sample1.Contains("MC_ggH")) scaleFactor*=10;

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

float analyzer_higgs::getNevents(string sample1, TH1F* h)
{

  Int_t nEv2 = 0;
  if(sample1=="MC_ZllG") nEv2 = 165000.;
  else if(sample1=="MC_Wjets") nEv2 = h    -> GetEntries();
  else if(sample1=="MC_tt2l2nu2b") nEv2 = h       -> GetEntries();
  else if(sample1=="MC_tSchannel") nEv2 = h       -> GetEntries();
  else if(sample1=="MC_tTchannel") nEv2 = h       -> GetEntries();
  else if(sample1=="MC_tWchannel") nEv2 = h       -> GetEntries();
  else if(sample1=="MC_ZZtoAny") nEv2 = h   -> GetEntries();
  else if(sample1=="MC_DYmumu") nEv2 = 1984154.;
  else if(sample1=="MC_DYee") nEv2 = 1992276.;
  else if(sample1=="MC_DYtautau") nEv2 = 1995369.;
  else if(sample1=="MC_WW") nEv2 = 110000.;
  else if(sample1=="MC_WZ") nEv2 = 110000.;
  else if(sample1=="MC_Zbb0") nEv2 = 347393.;
  else if(sample1=="MC_Zbb1") nEv2 = 2025811.;
  else if(sample1=="MC_Zbb2") nEv2 = 10842.;
  else if(sample1=="MC_Zbb3") nEv2 = 10898.;
  else if(sample1=="MC_Zcc0") nEv2 = 438067.;
  else if(sample1=="MC_Zcc1") nEv2 = 184365.;
  else if(sample1=="MC_Zcc2") nEv2 = 10735.;
  else if(sample1=="MC_Zcc3") nEv2 = 10177.;
  else if(sample1=="MC_ggH200") nEv2 = 109274;
  else if(sample1=="MC_ggH400") nEv2 = h   -> GetEntries();
  else if(sample1=="MC_ggH600") nEv2 = 109255;
  else if(sample1=="MC_ggH200_WW2l2nu") nEv2 = 109988;
  else if(sample1=="MC_ggH400_WW2l2nu") nEv2 = 109584;
  else if(sample1=="MC_ggH600_WW2l2nu") nEv2 = 109974;
  else if(sample1=="MC_ggH200_WW2t2nu") nEv2 = 65995;
  else if(sample1=="MC_ggH400_WW2t2nu") nEv2 = 65991;
  else if(sample1=="MC_ggH600_WW2t2nu") nEv2 = 65997;
  return nEv2;
}


float analyzer_higgs::weightZZ(Float_t x)
{
  Float_t a = 1.108, b = 0.002429, c = -1.655e-06;
  Float_t w = a + b*x + c*x*x;  
  return w;
}

inline void analyzer_higgs::FillHistosNoise(Int_t num, Double_t weight){
  met0_over_qt[num] -> Fill(METqt, weight);  
  met0_et[num]      -> Fill(MET, weight);
  met0_et_ovQt[num] -> Fill(MET, METqt, weight);
  met0_phi[num]     -> Fill(MET_phi, weight);

  met1_over_qt[num] -> Fill(MET1qt, weight);
  met1_et[num]      -> Fill(MET1, weight);
  met1_et_ovQt[num] -> Fill(MET1, MET1qt, weight);
  met1_phi[num]     -> Fill(MET1_phi, weight);

  mt0[num]          -> Fill(MT, weight);
}

inline void analyzer_higgs::FillHistos(Int_t num, Double_t weight){
      met2_over_qt[num] -> Fill(METqt, weight);
      met2_et[num]      -> Fill(MET, weight);
      met2_et_ovQt[num] -> Fill(MET, METqt, weight);
      met2_phi[num]     -> Fill(MET_phi, weight);

      met3_over_qt[num]  -> Fill(projMETqt, weight);
      met3_et[num]       -> Fill(projMET, weight);
      met3_et_ovQt[num] -> Fill(projMET, projMETqt, weight);

      met4_over_qt[num] -> Fill(puCorrMETqt, weight);
      met4_et[num]      -> Fill(puCorrMET, weight);
      met4_et_ovQt[num] -> Fill(puCorrMET, puCorrMETqt, weight);
      met4_puSig[num]   -> Fill(puSigMET, weight);

      mt2[num]      -> Fill(MT, weight);
      mtZ[num]      -> Fill(MTZ, weight);
      mt2_met2[num] -> Fill(MT, MET, weight);
      mtZ_met2[num] -> Fill(MTZ, MET, weight);
      mt2_met3[num] -> Fill(MT, projMET, weight);
      mtZ_met3[num] -> Fill(MTZ, projMET, weight);

      di_qt[num]    -> Fill(qT, weight);
      di_eta[num] -> Fill(diEta);
      di_mass[num]  -> Fill(Mll, weight);
      if(Mll_EB!=0) di_mass_EB[num]  -> Fill(Mll_EB, weight);
      if(Mll_EE!=0) di_mass_EE[num]  -> Fill(Mll_EE, weight);
      if(Mll_EX!=0) di_mass_EX[num]  -> Fill(Mll_EX, weight);

      jet_N[num]    -> Fill(nJets, weight);
      jet_b_N[num]  -> Fill(nJetsB, weight);

      if(jetCount>0)  
	{
	  met2_dPhiClosJet1[num] -> Fill(dPhiClos1, weight); 
	  met2_dPhiClosJet2[num] -> Fill(dPhiClos2, weight); 
	  met2_dPhiLeadJet1[num] -> Fill(dPhiLead1, weight); 
	  met2_dPhiLeadJet2[num] -> Fill(dPhiLead2, weight); 
	}

}

inline void analyzer_higgs::CountEvents(Int_t num)
{
  nEvents[num]++;
  if(!isNoiseHcal)        nEventsPassNoiseFilter[1][num]++; 
  if(!isDeadEcalCluster)  nEventsPassNoiseFilter[2][num]++; 
  if(!isScraping)         nEventsPassNoiseFilter[3][num]++;
  if(!isCSCTightHalo)     nEventsPassNoiseFilter[4][num]++;
  if(!isCSCLooseHalo)     nEventsPassNoiseFilter[5][num]++;
  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo) 
      nEventsPassNoiseFilter[0][num]++; 
}

/*
void analyzer_higgs::scale(string sample1, TH1 * h, Float_t cs, Float_t nTotEv, Float_t lumi)
{

  Float_t scaleFactor      = lumi*cs/nTotEv;
  if(sample1=="2011A_May10_DoubleMu" ||
     sample1=="2011A_May10_DoubleElectron" ||
     sample1=="2011A_HZZSkim_DoubleMu" ||
     sample1=="2011A_HZZSkim_DoubleElectron") scaleFactor = 1.0;
  if(sample1=="MC_ggH200" ||
     sample1=="MC_ggH400" ||
     sample1=="MC_ggH600" ) scaleFactor*=10;
     

  // cout<<sample1<<"  nEv: "<<nTotEv<<"   cs: "<<cs<<"   Rescaling all the histogram "<<h->GetName()<<" by: "<<scaleFactor<<endl;
  h -> Scale(scaleFactor);
}
*/

