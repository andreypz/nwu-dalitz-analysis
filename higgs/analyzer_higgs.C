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
Float_t CS, nEv;
Int_t lColor, fColor;

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

      mtZ_met2[n]     = new TH2F(Form("mtZ_met2_%i",n), "MTZ pf vs pfMet noise", 50, 100, 600, 40, 0,400); 
      mt2_met2[n]     = new TH2F(Form("mt2_met2_%i",n), "MT pf vs pfMet noise", 50, 100, 600, 40, 0,400); 

      mtZ_met3[n]     = new TH2F(Form("mtZ_met3_%i",n), "MTZ pf vs projMet noise", 50, 100, 600, 40, 0,400); 
      mt2_met3[n]     = new TH2F(Form("mt2_met3_%i",n), "MT pf vs projfMet noise", 50, 100, 600, 40, 0,400); 
 
      di_qt[n]     = new TH1F(Form("di_qt_%i",n), "di-lepton qt", 40, 0,400);  
      di_mass[n]   = new TH1F(Form("di_mass_%i",n), "di-lepton Mass", 50, 70,120);  
      di_mass_EB[n]   = new TH1F(Form("di_mass_EB_%i",n), "di-lepton Mass in Ecal Barrel", 50, 70,120);  
      di_mass_EE[n]   = new TH1F(Form("di_mass_EE_%i",n), "di-lepton Mass in Ecal Endcap", 50, 70,120);  
      di_mass_EX[n]   = new TH1F(Form("di_mass_EX_%i",n), "di-lepton Mass in Ecal E/B mix", 50, 70,120);  

      jet_N[n]  = new TH1F(Form("jet_N_%i",n), "Number of jets", 20,0,20);
      jet_b_N[n] = new TH1F(Form("jet_b_N_%i",n), "Number of b-jets", 10,0,10);

      jet_pt[n] = new TH1F(Form("jet_pt_%i",n), "Pt of jets, eta<2.4", 40,0,200);
      jet_b_pt[n] = new TH1F(Form("jet_b_pt_%i",n), "Pt of b-jets, eta<2.4", 40,0,200);


    }


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
  nEvents[0]++; 
  if(!isNoiseHcal)        nEventsPassNoiseFilter[1][0]++; 
  if(!isDeadEcalCluster)  nEventsPassNoiseFilter[2][0]++; 
  if(!isScraping)         nEventsPassNoiseFilter[3][0]++;
  if(!isCSCTightHalo)     nEventsPassNoiseFilter[4][0]++;
  if(!isCSCLooseHalo)     nEventsPassNoiseFilter[5][0]++;
  if(!isNoiseHcal && !isDeadEcalCluster&& !isScraping && !isCSCTightHalo) 
    nEventsPassNoiseFilter[0][0]++; 
    
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
  Float_t Mll=0, Mll_EE=0, Mll_EB=0, Mll_EX=0, pTll=0;
  //--------------//
  //--- Muons------//
  //-------------//
  TLorentzVector diLepton(0.,0.,0.,0.), Lepton1(0.,0.,0.,0.), Lepton2(0.,0.,0.,0.);
 
  TCMuon *myMuon1=0, *myMuon2=0, *myMuon3=0, *myMuon4=0;
  Int_t muIndex[5] = {-1,-1,-1,-1,-1};

  Int_t muCount=0, muCountLoose = 0;
  for (Int_t i = 0; i < recoMuons->GetSize(); ++i) 
    {    
      TCMuon* thisMuon = (TCMuon*) recoMuons->At(i);     
      TVector3 *muVertex =  new TVector3(thisMuon->Vtx());
      TVector3 *muMomentum = new TVector3( thisMuon->p4().Vect());
      
      double dxyVtx = dxy(mainPrimaryVertex, muVertex, muMomentum);
      double dzVtx = dz(mainPrimaryVertex, muVertex, muMomentum);
      
      if (fabs(thisMuon->p4().Eta()) < 2.4 && 
	  thisMuon->p4().Pt()>10     &&
	  thisMuon->isGLB()  && 
	  fabs(dxyVtx)<0.02  && 
	  fabs(dzVtx)<0.1    &&
	  thisMuon->normChi2()<10    &&
	  thisMuon->nTRKHits()>10     &&
	  thisMuon->nPXLHits()>0     &&
	  thisMuon->nValidMuHits()>0 &&
	  thisMuon->nMatchSeg()>1    &&
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
          fabs(dzVtx)<0.2   
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
	  fabs(dzVtx)  < 0.2                
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


  if (isMuonSample && muCount==2)
    {
      Mll  = (myMuon1->p4()+myMuon2->p4()).M();
      pTll = (myMuon1->p4()+myMuon2->p4()).Pt();
      diLepton = (TLorentzVector)(myMuon1->p4()+myMuon2->p4());
      //if((diLepton.M()-Mll)>0.01 || (diLepton.Pt()-pTll)>0.01) cout<<" Warning!  Dilepton is wrong"
      //							   <<diLepton.M()<<"  "<<Mll<<endl;
      //cout<<"dbg  "<<muCount<<endl;
      if (eleCount==0 && myMuon1->charge()*myMuon2->charge()==-1 && Mll>40 && Mll<500)	 
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
      if (muCount==0 && myElectron1->charge()*myElectron2->charge()==-1 && Mll>40 && Mll<500)	 
	{
	  //cout<<"preselection?!"<<endl;
	  passPreSelection= kTRUE;
	  if (Mll>cut_Mll_low && Mll<cut_Mll_high) passZpeak = kTRUE;
	}

    }


  Float_t MET=0, MET_phi=0;
  Float_t MET1=0, MET1_phi=0;
  Float_t MT = 0, MTZ = 0;
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


  Float_t qT = diLepton.Pt();

  Float_t projMET     = projectedMET(MET1, MET1_phi, Lepton1, Lepton2);
  Float_t projMETqt   = projMET/qT;
  Float_t puCorrMET   = pileUpCorrectedMET(projMET, nVtx);
  Float_t puCorrMETqt = puCorrMET/qT;
  Float_t sigma0      = 7.343;
  Float_t puSigMET    = puCorrMET/sigma0;

  //cout<<"pfMET: "<<MET<<"  pfType1: "<<MET1<<"  proj1: "<<projMET<<"  puCorrMET: "<< puCorrMET<<"   signif: "<<puSigMET<<"  nVtx: "<<nVtx<<endl;
  TVector2 ZptXY(0.0,0.0), METXY(0.0,0.0);  
  ZptXY.Set(diLepton.X(), diLepton.Y());
  METXY.Set(MET*cos(MET_phi), MET*sin(MET_phi));
  MT = higgsMT(pTll, Mll, MET, ZptXY, METXY);
  MTZ = higgsMT(pTll, 91.1876, MET, ZptXY, METXY);
  //if (muCount==2) cout<<"MT: "<<MT<<"    "<<ZptXY.X()<<" "<<ZptXY.Y()<<"    "<<METXY.X()<<" "<<METXY.Y()<<endl;

  Float_t METqt = MET/qT;
  Float_t MET1qt = MET1/qT;

  Bool_t antiB  = kFALSE;
  Int_t jetCount=0, nJets = 0, nJetsB = 0;

  Float_t ptJ1 = 0, ptJ2=0;
  Float_t  dPhiLead1 = 110,  dPhiLead2 = 110; 
  Float_t  dPhiClos1 = 110,  dPhiClos2 = 110; 
  for (Int_t i = 0; i < recoJets->GetSize(); ++i) 
    {    
      TCJet* thisJet = (TCJet*) recoJets->At(i);     
      if (thisJet->P4(7).Pt() > 30 && fabs(thisJet->P4(7).Eta()) <= 2.4) 
	{
	  //do anti-b tag
          if (thisJet->BDiscrTrkCountHiEff() > 2)
	    {
	      antiB = kTRUE;//return kTRUE; 
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

  for(Int_t ii=0; ii<nC; ii++)
    {
      //Fill out 
      if(ii!=8 && ii!=4) jet_b_pt[ii] ->Fill(66);
    }

  met0_over_qt[0] -> Fill(METqt);
  met0_et[0]      -> Fill(MET);
  met0_et_ovQt[0] -> Fill(MET, METqt);
  met0_phi[0]     -> Fill(MET_phi);

  met1_over_qt[0] -> Fill(MET1qt);
  met1_et[0]      -> Fill(MET1);
  met1_et_ovQt[0] -> Fill(MET1, MET1qt);
  met1_phi[0]     -> Fill(MET1_phi);

  mt0[0] -> Fill(MT);

  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo) 
    {
      di_qt[0] -> Fill(qT);
      di_mass[0]  -> Fill(Mll);
      if(Mll_EB!=0) di_mass_EB[0]  -> Fill(Mll_EB);
      if(Mll_EE!=0) di_mass_EE[0]  -> Fill(Mll_EE);
      if(Mll_EX!=0) di_mass_EX[0]  -> Fill(Mll_EX);
      
      met2_over_qt[0] -> Fill(METqt);
      met2_et[0]      -> Fill(MET);
      met2_et_ovQt[0] -> Fill(MET, METqt);
      met2_phi[0]     -> Fill(MET_phi);

      met3_over_qt[0] -> Fill(projMETqt);
      met3_et[0]      -> Fill(projMET);
      met3_et_ovQt[0] -> Fill(projMET, projMETqt);

      met4_over_qt[0] -> Fill(puCorrMETqt);
      met4_et[0]      -> Fill(puCorrMET);
      met4_et_ovQt[0] -> Fill(puCorrMET, puCorrMETqt);
      met4_puSig[0]   -> Fill(puSigMET);

      mt2[0]      -> Fill(MT);
      mtZ[0]      -> Fill(MTZ);
      mt2_met2[0] -> Fill(MT, MET);
      mtZ_met2[0] -> Fill(MTZ, MET);
      mt2_met3[0] -> Fill(MT, projMET);
      mtZ_met3[0] -> Fill(MTZ, projMET);

      jet_N[0]    -> Fill(nJets);
      jet_b_N[0]   -> Fill(nJetsB);

      if(jetCount>0)  
	{
	  met2_dPhiClosJet1[0] -> Fill(dPhiClos1); 
	  met2_dPhiClosJet2[0] -> Fill(dPhiClos2); 
	  met2_dPhiLeadJet1[0] -> Fill(dPhiLead1); 
	  met2_dPhiLeadJet2[0] -> Fill(dPhiLead2); 
	}


    }
  else
    if(verboseLvl>1) 
      nout[0]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<MET<<"\t* "
	     <<isNoiseHcal<<"\t* "<<isDeadEcalCluster<<"\t* "<<isScraping<<"\t* "<<isCSCTightHalo<<"\t* "<<endl;  

  if (!passPreSelection) return kTRUE;

  // Only make sense after preselection (two leptons):
  if(fabs(Lepton1.Eta()) < 1.444 && fabs(Lepton2.Eta())<1.444) Mll_EB = diLepton.M();
  else  if(fabs(Lepton1.Eta()) > 1.566 && fabs(Lepton2.Eta())>1.566)  Mll_EE = diLepton.M();
  else  Mll_EX = diLepton.M();
  
  nEvents[1]++;
  if(!isNoiseHcal)        nEventsPassNoiseFilter[1][1]++; 
  if(!isDeadEcalCluster)  nEventsPassNoiseFilter[2][1]++; 
  if(!isScraping)         nEventsPassNoiseFilter[3][1]++;
  if(!isCSCTightHalo)     nEventsPassNoiseFilter[4][1]++;
  if(!isCSCLooseHalo)     nEventsPassNoiseFilter[5][1]++;
  if(!isNoiseHcal && !isDeadEcalCluster  && !isScraping && !isCSCTightHalo) 
    {
      nEventsPassNoiseFilter[0][1]++; 
      met2_over_qt[1] -> Fill(METqt);
      met2_et[1]      -> Fill(MET);
      met2_et_ovQt[1] -> Fill(MET, METqt);
      met2_phi[1]     -> Fill(MET_phi);

      met3_over_qt[1]  -> Fill(projMETqt);
      met3_et[1]       -> Fill(projMET);
      met3_et_ovQt[1] -> Fill(projMET, projMETqt);

      met4_over_qt[1] -> Fill(puCorrMETqt);
      met4_et[1]      -> Fill(puCorrMET);
      met4_et_ovQt[1] -> Fill(puCorrMET, puCorrMETqt);
      met4_puSig[1]   -> Fill(puSigMET);

      mt2[1]      -> Fill(MT);
      mtZ[1]      -> Fill(MTZ);
      mt2_met2[1] -> Fill(MT, MET);
      mtZ_met2[1] -> Fill(MTZ, MET);
      mt2_met3[1] -> Fill(MT, projMET);
      mtZ_met3[1] -> Fill(MTZ, projMET);

      di_qt[1]    -> Fill(qT);
      di_mass[1]  -> Fill(Mll);
      if(Mll_EB!=0) di_mass_EB[1]  -> Fill(Mll_EB);
      if(Mll_EE!=0) di_mass_EE[1]  -> Fill(Mll_EE);
      if(Mll_EX!=0) di_mass_EX[1]  -> Fill(Mll_EX);

      jet_N[1]    -> Fill(nJets);
      jet_b_N[1]  -> Fill(nJetsB);

      if(jetCount>0)  
	{
	  met2_dPhiClosJet1[1] -> Fill(dPhiClos1); 
	  met2_dPhiClosJet2[1] -> Fill(dPhiClos2); 
	  met2_dPhiLeadJet1[1] -> Fill(dPhiLead1); 
	  met2_dPhiLeadJet2[1] -> Fill(dPhiLead2); 
	}
    }
  else
    if(verboseLvl>0) 
      nout[1]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<MET<<"\t* "
	     <<isNoiseHcal<<"\t* "<<isDeadEcalCluster<<"\t* "<<isScraping<<"\t* "<<isCSCTightHalo<<"\t* "<<endl;  

  met0_over_qt[1] -> Fill(METqt);
  met0_et[1]      -> Fill(MET);
  met0_et_ovQt[1] -> Fill(MET, METqt);
  met0_phi[1]     -> Fill(MET_phi);

  met1_over_qt[1] -> Fill(MET1qt);
  met1_et[1]      -> Fill(MET1);
  met1_et_ovQt[1] -> Fill(MET1, MET1qt);
  met1_phi[1]     -> Fill(MET1_phi);

  mt0[1] -> Fill(MT);
    
  
if (!passZpeak) return kTRUE;
  nEvents[2]++;
  if(!isNoiseHcal)        nEventsPassNoiseFilter[1][2]++; 
  if(!isDeadEcalCluster)  nEventsPassNoiseFilter[2][2]++; 
  if(!isScraping)         nEventsPassNoiseFilter[3][2]++;
  if(!isCSCTightHalo)     nEventsPassNoiseFilter[4][2]++;
  if(!isCSCLooseHalo)     nEventsPassNoiseFilter[5][2]++;
  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo) 
    {
      nEventsPassNoiseFilter[0][2]++; 
      met2_over_qt[2] -> Fill(METqt);
      met2_et[2]      -> Fill(MET);
      met2_et_ovQt[2] -> Fill(MET, METqt);
      met2_phi[2]     -> Fill(MET_phi);

      met3_over_qt[2] -> Fill(projMETqt);
      met3_et[2]      -> Fill(projMET);
      met3_et_ovQt[2] -> Fill(projMET, projMETqt);
      
      met4_over_qt[2] -> Fill(puCorrMETqt);
      met4_et[2]      -> Fill(puCorrMET);
      met4_et_ovQt[2] -> Fill(puCorrMET, puCorrMETqt);
      met4_puSig[2]     -> Fill(puSigMET);

      mt2[2]      -> Fill(MT);
      mtZ[2]      -> Fill(MTZ);
      mt2_met2[2] -> Fill(MT, MET);
      mtZ_met2[2] -> Fill(MTZ, MET);
      mt2_met3[2] -> Fill(MT, projMET);
      mtZ_met3[2] -> Fill(MTZ, projMET);

      di_qt[2] -> Fill(qT);
      di_mass[2]  -> Fill(Mll);
      if(Mll_EB!=0) di_mass_EB[2]  -> Fill(Mll_EB);
      if(Mll_EE!=0) di_mass_EE[2]  -> Fill(Mll_EE);
      if(Mll_EX!=0) di_mass_EX[2]  -> Fill(Mll_EX);

      jet_N[2]    -> Fill(nJets);
      jet_b_N[2]   -> Fill(nJetsB);


      if(jetCount>0)  
	{
	  met2_dPhiClosJet1[2] -> Fill(dPhiClos1); 
	  met2_dPhiClosJet2[2] -> Fill(dPhiClos2); 
	  met2_dPhiLeadJet1[2] -> Fill(dPhiLead1); 
	  met2_dPhiLeadJet2[2] -> Fill(dPhiLead2); 
	}

    }

  met0_over_qt[2] -> Fill(METqt);
  met0_et[2]      -> Fill(MET);
  met0_et_ovQt[2] -> Fill(MET, METqt);
  met0_phi[2]     -> Fill(MET_phi);

  met1_over_qt[2] -> Fill(MET1qt);
  met1_et[2]      -> Fill(MET1);
  met1_et_ovQt[2] -> Fill(MET1, MET1qt);
  met1_phi[2]     -> Fill(MET1_phi);

  mt0[2] -> Fill(MT);
 
  if (isMuonSample)
    {
      mu1_eta[2] -> Fill(myMuon1->eta());
      mu1_phi[2] -> Fill(myMuon1->phi());
      mu1_pt[2]  -> Fill(myMuon1->pt());

      mu2_eta[2] -> Fill(myMuon2->eta());
      mu2_phi[2] -> Fill(myMuon2->phi());
      mu2_pt[2]  -> Fill(myMuon2->pt());

      if (myMuon1->pt() < 20) return kTRUE;
      nEvents[3]++;
      if(!isNoiseHcal)        nEventsPassNoiseFilter[1][3]++; 
      if(!isDeadEcalCluster)  nEventsPassNoiseFilter[2][3]++; 
      if(!isScraping)         nEventsPassNoiseFilter[3][3]++;
      if(!isCSCTightHalo)     nEventsPassNoiseFilter[4][3]++;
      if(!isCSCLooseHalo)     nEventsPassNoiseFilter[5][3]++;
      if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo) 
	nEventsPassNoiseFilter[0][3]++; 

      if (myMuon2->pt() < 20) return kTRUE;
      nEvents[4]++;
      if(!isNoiseHcal)        nEventsPassNoiseFilter[1][4]++; 
      if(!isDeadEcalCluster)  nEventsPassNoiseFilter[2][4]++; 
      if(!isScraping)         nEventsPassNoiseFilter[3][4]++;
      if(!isCSCTightHalo)     nEventsPassNoiseFilter[4][4]++;
      if(!isCSCLooseHalo)     nEventsPassNoiseFilter[5][4]++;
      if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo) 
	nEventsPassNoiseFilter[0][4]++; 

    }

  if (isElectronSample)
    {
      if (myElectron1->pt() < 20) return kTRUE;
      nEvents[3]++;
      if(!isNoiseHcal)        nEventsPassNoiseFilter[1][3]++; 
      if(!isDeadEcalCluster)  nEventsPassNoiseFilter[2][3]++; 
      if(!isScraping)         nEventsPassNoiseFilter[3][3]++;
      if(!isCSCTightHalo)     nEventsPassNoiseFilter[4][3]++;
      if(!isCSCLooseHalo)     nEventsPassNoiseFilter[5][3]++;
      if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo) 
	nEventsPassNoiseFilter[0][3]++; 
  
      if (myElectron2->pt() < 20) return kTRUE;
      nEvents[4]++;
      if(!isNoiseHcal)        nEventsPassNoiseFilter[1][4]++; 
      if(!isDeadEcalCluster)  nEventsPassNoiseFilter[2][4]++; 
      if(!isScraping)         nEventsPassNoiseFilter[3][4]++;
      if(!isCSCTightHalo)     nEventsPassNoiseFilter[4][4]++;
      if(!isCSCLooseHalo)     nEventsPassNoiseFilter[5][4]++;
      if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo) 
	nEventsPassNoiseFilter[0][4]++; 
  
    }
     
  if (isMuonSample)
    {
      mu1_eta[4] -> Fill(myMuon1->eta());
      mu1_phi[4] -> Fill(myMuon1->phi());
      mu1_pt[4]  -> Fill(myMuon1->pt());

      mu2_eta[4] -> Fill(myMuon2->eta());
      mu2_phi[4] -> Fill(myMuon2->phi());
      mu2_pt[4]  -> Fill(myMuon2->pt());

    }

  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo)
    { 
      met2_over_qt[4] -> Fill(METqt);
      met2_et[4]      -> Fill(MET);
      met2_et_ovQt[4] -> Fill(MET, METqt);
      met2_phi[4]     -> Fill(MET_phi);

      met3_over_qt[4] -> Fill(projMETqt);
      met3_et[4]      -> Fill(projMET);
      met3_et_ovQt[4] -> Fill(projMET, projMETqt);
      
      met4_over_qt[4] -> Fill(puCorrMETqt);
      met4_et[4]      -> Fill(puCorrMET);
      met4_et_ovQt[4] -> Fill(puCorrMET, puCorrMETqt);
      met4_puSig[4]   -> Fill(puSigMET);

      mt2[4]      -> Fill(MT);
      mtZ[4]      -> Fill(MTZ);
      mt2_met2[4] -> Fill(MT, MET);
      mtZ_met2[4] -> Fill(MTZ, MET);
      mt2_met3[4] -> Fill(MT, projMET);
      mtZ_met3[4] -> Fill(MTZ, projMET);

      di_qt[4]    -> Fill(qT);
      di_mass[4]  -> Fill(Mll);
      if(Mll_EB!=0) di_mass_EB[4]  -> Fill(Mll_EB);
      if(Mll_EE!=0) di_mass_EE[4]  -> Fill(Mll_EE);
      if(Mll_EX!=0) di_mass_EX[4]  -> Fill(Mll_EX);

      jet_N[4]    -> Fill(nJets);
      jet_b_N[4]  -> Fill(nJetsB);

      if(jetCount>0)  
	{
	  met2_dPhiClosJet1[4] -> Fill(dPhiClos1); 
	  met2_dPhiClosJet2[4] -> Fill(dPhiClos2); 
	  met2_dPhiLeadJet1[4] -> Fill(dPhiLead1); 
	  met2_dPhiLeadJet2[4] -> Fill(dPhiLead2); 
	}

    }
  
  met0_over_qt[4] -> Fill(METqt);
  met0_et[4]      -> Fill(MET);
  met0_et_ovQt[4] -> Fill(MET, METqt);
  met0_phi[4]     -> Fill(MET_phi);

  met1_over_qt[4] -> Fill(MET1qt);
  met1_et[4]      -> Fill(MET1);
  met1_et_ovQt[4] -> Fill(MET1, MET1qt);
  met1_phi[4]     -> Fill(MET1_phi);

  mt0[4] -> Fill(MT);

  if(verboseLvl>1) 
    fout[4]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<Mll<<"\t* "<<MT<<"\t* "<<MET<<"\t* "<<pTll<<"\t* "<<rhoFactor<<"\t*"<<endl;  


  for (Int_t i = 0; i < recoJets->GetSize(); ++i) 
    {    
      TCJet* thisJet = (TCJet*) recoJets->At(i);     
      if (fabs(thisJet->P4(7).Eta()) <= 2.4 && thisJet->BDiscrTrkCountHiEff() > 2) 
	jet_b_pt[4] ->Fill(thisJet->P4(7).Pt());
    }
  
  if (projMET>70)
    {
      jet_b_N[7]   -> Fill(nJetsB);
      if(MT>327 && MT<441)
	jet_b_N[8]   -> Fill(nJetsB);

     for (Int_t i = 0; i < recoJets->GetSize(); ++i) 
	{    
	  TCJet* thisJet2 = (TCJet*) recoJets->At(i);     
	  if (fabs(thisJet2->P4(7).Eta()) <= 2.4 && thisJet2->BDiscrTrkCountHiEff() > 2) 
	    jet_b_pt[8] ->Fill(thisJet2->P4(7).Pt());
	}
    }



  if(antiB  &&  pTll<cut_pTll) return kTRUE;


  nEvents[5]++;

  if(!isNoiseHcal)        nEventsPassNoiseFilter[1][5]++; 
  if(!isDeadEcalCluster)  nEventsPassNoiseFilter[2][5]++; 
  if(!isScraping)         nEventsPassNoiseFilter[3][5]++;
  if(!isCSCTightHalo)     nEventsPassNoiseFilter[4][5]++;
  if(!isCSCLooseHalo)     nEventsPassNoiseFilter[5][5]++;
  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo) 
    {
      nEventsPassNoiseFilter[0][5]++; 
      met2_over_qt[5] -> Fill(METqt);
      met2_et[5]      -> Fill(MET);
      met2_et_ovQt[5] -> Fill(MET, METqt);
      met2_phi[5]     -> Fill(MET_phi);
    
      met3_over_qt[5] -> Fill(projMETqt);
      met3_et[5]      -> Fill(projMET);
      met3_et_ovQt[5] -> Fill(projMET, projMETqt);
 
      met4_over_qt[5] -> Fill(puCorrMETqt);
      met4_et[5]      -> Fill(puCorrMET);
      met4_et_ovQt[5] -> Fill(puCorrMET, puCorrMETqt);
      met4_puSig[5]   -> Fill(puSigMET);

      mt2[5]      -> Fill(MT);
      mtZ[5]      -> Fill(MTZ);
      mt2_met2[5] -> Fill(MT, MET);
      mtZ_met2[5] -> Fill(MTZ, MET);
      mt2_met3[5] -> Fill(MT, projMET);
      mtZ_met3[5] -> Fill(MTZ, projMET);

      di_qt[5]    -> Fill(qT);
      di_mass[5]  -> Fill(Mll);
      if(Mll_EB!=0) di_mass_EB[5]  -> Fill(Mll_EB);
      if(Mll_EE!=0) di_mass_EE[5]  -> Fill(Mll_EE);
      if(Mll_EX!=0) di_mass_EX[5]  -> Fill(Mll_EX);

      jet_N[5]    -> Fill(nJets);
      jet_b_N[5]  -> Fill(nJetsB);

      if(jetCount>0)  
	{
	  met2_dPhiClosJet1[5] -> Fill(dPhiClos1); 
	  met2_dPhiClosJet2[5] -> Fill(dPhiClos2); 
	  met2_dPhiLeadJet1[5] -> Fill(dPhiLead1); 
	  met2_dPhiLeadJet2[5] -> Fill(dPhiLead2); 
	}

      if(projMET<70)
	{
	  //Reversing Met cut
	  met2_over_qt[9] -> Fill(METqt);
	  met2_et[9]      -> Fill(MET);
	  met2_et_ovQt[9] -> Fill(MET, METqt);
	  met2_phi[9]     -> Fill(MET_phi);
	  
	  met3_over_qt[9] -> Fill(projMETqt);
	  met3_et[9]      -> Fill(projMET);
	  met3_et_ovQt[9] -> Fill(projMET, projMETqt);
	  
	  met4_over_qt[9] -> Fill(puCorrMETqt);
	  met4_et[9]      -> Fill(puCorrMET);
	  met4_et_ovQt[9] -> Fill(puCorrMET, puCorrMETqt);
	  met4_puSig[9]   -> Fill(puSigMET);
	  
	  mt2[9]      -> Fill(MT);
	  mtZ[9]      -> Fill(MTZ);
	  mt2_met2[9] -> Fill(MT, MET);
	  mtZ_met2[9] -> Fill(MTZ, MET);
	  mt2_met3[9] -> Fill(MT, projMET);
	  mtZ_met3[9] -> Fill(MTZ, projMET);
	  
	  di_qt[9]    -> Fill(qT);
	  di_mass[9]  -> Fill(Mll);
	  if(Mll_EB!=0) di_mass_EB[9]  -> Fill(Mll_EB);
	  if(Mll_EE!=0) di_mass_EE[9]  -> Fill(Mll_EE);
	  if(Mll_EX!=0) di_mass_EX[9]  -> Fill(Mll_EX);

	  jet_N[9]    -> Fill(nJets);
	  jet_b_N[9]  -> Fill(nJetsB);
	  
	  if(jetCount>0)  
	    {
	      met2_dPhiClosJet1[9] -> Fill(dPhiClos1); 
	      met2_dPhiClosJet2[9] -> Fill(dPhiClos2); 
	      met2_dPhiLeadJet1[9] -> Fill(dPhiLead1); 
	      met2_dPhiLeadJet2[9] -> Fill(dPhiLead2); 
	    }
	}
    }

  met0_over_qt[5] -> Fill(METqt);
  met0_et[5]      -> Fill(MET);
  met0_et_ovQt[5] -> Fill(MET, METqt);
  met0_phi[5]     -> Fill(MET_phi);

  met1_over_qt[5] -> Fill(MET1qt);
  met1_et[5]      -> Fill(MET1);
  met1_et_ovQt[5] -> Fill(MET1, MET1qt);
  met1_phi[5]     -> Fill(MET1_phi);

  mt0[5] -> Fill(MT);

  if (isMuonSample)
    {
      mu1_eta[5] -> Fill(myMuon1->eta());
      mu1_phi[5] -> Fill(myMuon1->phi());
      mu1_pt[5]  -> Fill(myMuon1->pt());

      mu2_eta[5] -> Fill(myMuon2->eta());
      mu2_phi[5] -> Fill(myMuon2->phi());
      mu2_pt[5]  -> Fill(myMuon2->pt());

    }
  

  if(verboseLvl>1)   
    fout[5]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<Mll<<"\t* "<<MT<<"\t* "<<MET<<"\t* "<<pTll<<"\t* "<<rhoFactor<<"\t*"<<endl;  

  // if (pTll<50) return kTRUE;
 


  if (projMET<70) return kTRUE;
 
 nEvents[6]++;

  if(!isNoiseHcal)        nEventsPassNoiseFilter[1][6]++; 
  if(!isDeadEcalCluster)  nEventsPassNoiseFilter[2][6]++; 
  if(!isScraping)         nEventsPassNoiseFilter[3][6]++;
  if(!isCSCTightHalo)     nEventsPassNoiseFilter[4][6]++;
  if(!isCSCLooseHalo)     nEventsPassNoiseFilter[5][6]++;
  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo) 
    {
      nEventsPassNoiseFilter[0][6]++; 
      met2_over_qt[6] -> Fill(METqt);
      met2_et[6]      -> Fill(MET);
      met2_et_ovQt[6] -> Fill(MET, METqt);
      met2_phi[6]     -> Fill(MET_phi);
      
      met3_over_qt[6] -> Fill(projMETqt);
      met3_et[6]      -> Fill(projMET);
      met3_et_ovQt[6] -> Fill(projMET, projMETqt);

      met4_over_qt[6] -> Fill(puCorrMETqt);
      met4_et[6]      -> Fill(puCorrMET);
      met4_et_ovQt[6] -> Fill(puCorrMET, puCorrMETqt);
      met4_puSig[6]   -> Fill(puSigMET);

      mt2[6]      -> Fill(MT);
      mtZ[6]      -> Fill(MTZ);
      mt2_met2[6] -> Fill(MT, MET);
      mtZ_met2[6] -> Fill(MTZ, MET);
      mt2_met3[6] -> Fill(MT, projMET);
      mtZ_met3[6] -> Fill(MTZ, projMET);

      di_qt[6]    -> Fill(qT);
      di_mass[6]  -> Fill(Mll);
      if(Mll_EB!=0) di_mass_EB[6]  -> Fill(Mll_EB);
      if(Mll_EE!=0) di_mass_EE[6]  -> Fill(Mll_EE);
      if(Mll_EX!=0) di_mass_EX[6]  -> Fill(Mll_EX);

      jet_N[6]    -> Fill(nJets);
      jet_b_N[6]  -> Fill(nJetsB);

      if(jetCount>0)  
	{
	  met2_dPhiClosJet1[6] -> Fill(dPhiClos1); 
	  met2_dPhiClosJet2[6] -> Fill(dPhiClos2); 
	  met2_dPhiLeadJet1[6] -> Fill(dPhiLead1); 
	  met2_dPhiLeadJet2[6] -> Fill(dPhiLead2); 
	}


    }
  else
    if(verboseLvl>0) 
      nout[6]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<MET<<"\t* "
	     <<isNoiseHcal<<"\t* "<<isDeadEcalCluster<<"\t* "<<isScraping<<"\t* "<<isCSCTightHalo<<"\t* "<<endl;  
  
  met0_over_qt[6] -> Fill(METqt);
  met0_et[6]      -> Fill(MET);
  met0_et_ovQt[6] -> Fill(MET, METqt);
  met0_phi[6]     -> Fill(MET_phi);

  met1_over_qt[6] -> Fill(MET1qt);
  met1_et[6]      -> Fill(MET1);
  met1_et_ovQt[6] -> Fill(MET1, MET1qt);
  met1_phi[6]     -> Fill(MET1_phi);

  mt0[6] -> Fill(MT);

  if (isMuonSample)
    {
      mu1_eta[6] -> Fill(myMuon1->eta());
      mu1_phi[6] -> Fill(myMuon1->phi());
      mu1_pt[6]  -> Fill(myMuon1->pt());

      mu2_eta[6] -> Fill(myMuon2->eta());
      mu2_phi[6] -> Fill(myMuon2->phi());
      mu2_pt[6]  -> Fill(myMuon2->pt());

    }

  if(verboseLvl>0) 
    fout[6]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<Mll<<"\t* "<<MT<<"\t* "<<MET<<"\t* "<<pTll<<"\t* "<<rhoFactor<<"\t* "<<muCountLoose<<"\t* "<<eleCountLoose<<endl;  


 nEvents[7]++;
  if(!isNoiseHcal)        nEventsPassNoiseFilter[1][7]++; 
  if(!isDeadEcalCluster)  nEventsPassNoiseFilter[2][7]++; 
  if(!isScraping)         nEventsPassNoiseFilter[3][7]++;
  if(!isCSCTightHalo)     nEventsPassNoiseFilter[4][7]++;
  if(!isCSCLooseHalo)     nEventsPassNoiseFilter[5][7]++;
  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo) 
    {
      nEventsPassNoiseFilter[0][7]++; 
      met2_over_qt[7] -> Fill(METqt);
      met2_et[7]      -> Fill(MET);
      met2_et_ovQt[7] -> Fill(MET, METqt);
      met2_phi[7]     -> Fill(MET_phi);

      met3_over_qt[7] -> Fill(projMETqt);
      met3_et[7]      -> Fill(projMET);
      met3_et_ovQt[7] -> Fill(projMET, projMETqt);
      
      met4_over_qt[7] -> Fill(puCorrMETqt);
      met4_et[7]      -> Fill(puCorrMET);
      met4_et_ovQt[7] -> Fill(puCorrMET, puCorrMETqt);
      met4_puSig[7]   -> Fill(puSigMET);
    }
  else
    if(verboseLvl>0) 
      nout[7]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<MET<<"\t* "
	     <<isNoiseHcal<<"\t* "<<isDeadEcalCluster<<"\t* "<<isScraping<<"\t* "<<isCSCTightHalo<<"\t* "<<endl;  
  
  met0_over_qt[7] -> Fill(METqt);
  met0_et[7]      -> Fill(MET);
  met0_et_ovQt[7] -> Fill(MET, METqt);
  met0_phi[7]     -> Fill(MET_phi);

  met1_over_qt[7] -> Fill(MET1qt);
  met1_et[7]      -> Fill(MET1);
  met1_et_ovQt[7] -> Fill(MET1, MET1qt);
  met1_phi[7]     -> Fill(MET1_phi);

  mt0[7] -> Fill(MT);

  if (isMuonSample)
    {
      mu1_eta[7] -> Fill(myMuon1->eta());
      mu1_phi[7] -> Fill(myMuon1->phi());
      mu1_pt[7]  -> Fill(myMuon1->pt());

      mu2_eta[7] -> Fill(myMuon2->eta());
      mu2_phi[7] -> Fill(myMuon2->phi());
      mu2_pt[7]  -> Fill(myMuon2->pt());

    }

  if(verboseLvl>0) 
    fout[7]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<Mll<<"\t* "<<MT<<"\t* "<<MET<<"\t* "<<pTll<<"\t* "<<rhoFactor<<"\t*"<<endl;  

  nEvents[8]++;

  nEvents[9]++;
  nEvents[10]++;

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

  //TFile *temp = new TFile("temp.root", "RECREATE"); 
  //temp->Write();
  //   gDirectory -> GetListOfKeys() -> Print();

  histoFile->Write();
  scaleAndColor(sample.c_str(), CS, nEv1, 1000.0, lColor, fColor); //1fb
  //  scaleAndColor(sample.c_str(), CS, nEv1, 191.0, lColor, fColor);
  //scale(sample.c_str(), CS, nEv1, 191.0);
  
  // cout<<"met0_et_0 integral: "<<met0_et[0]->Integral()<<endl;  //histoFile->Write();

  //histoFile->Print("-d");
  histoFile->Close();
  //temp ->Close();
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

double analyzer_higgs::higgsMT(Float_t pTll, Float_t Mll, Float_t MET, TVector2 ZptXY, TVector2 METXY)
{
  double MT = sqrt( pow(sqrt(pTll*pTll + Mll*Mll) + sqrt(MET*MET + Mll*Mll), 2) - (ZptXY + METXY).Mod2());
  return MT;
}

double analyzer_higgs::projectedMET(Float_t MET, Float_t MET_phi, TLorentzVector Lepton1, TLorentzVector Lepton2)
{
  Float_t deltaPhiMin =   TMath::Min( fabs(TVector2::Phi_mpi_pi(Lepton1.Phi() - MET_phi)), fabs(TVector2::Phi_mpi_pi(Lepton2.Phi() - MET_phi)));
  //cout<<deltaPhiMin<<endl;
  Float_t projMET = 0;
  if (deltaPhiMin > TMath::Pi()/2) projMET = MET;
  else projMET = MET*sin(deltaPhiMin);
  return projMET;
}

double analyzer_higgs::pileUpCorrectedMET(Float_t projMET, Int_t nVtx)
{
  //See talk by Marionneau Matthieu at DibosonMeeting on May 11, 2011
  Float_t sigma0  = 7.343;//, sigma0_err = 0.111;
  Float_t sigmaPU = 3.9;
  Float_t a = 0.856; //slope
  Float_t sigmaT = sqrt(sigma0*sigma0 + (nVtx-1)*sigmaPU*sigmaPU);
  Float_t projMET2 = projMET - a*(nVtx-1);
  Float_t puCorrMET = projMET2*sigma0/sigmaT;
  return puCorrMET;
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

/*
void analyzer_higgs::CountEvents(UInt_t nEv[], UInt_t nEvPassNoiseFilter[])
{
  nEv[0]++; 
  if(!isNoiseHcal)  nEvPassNoiseFilter[1][0]++; 
  if(!isDeadEcalCluster)  nEvPassNoiseFilter[2][0]++; 
  if(!isScraping)  nEvPassNoiseFilter[3][0]++;
  if(!isCSCTightHalo)  nEvPassNoiseFilter[4][0]++;
  if(!isCSCLooseHalo)  nEvPassNoiseFilter[5][0]++;
  if(!isNoiseHcal && !isDeadEcalCluster&& !isScraping && !isCSCTightHalo)  nEvPassNoiseFilter[0][0]++; 
}
*/
