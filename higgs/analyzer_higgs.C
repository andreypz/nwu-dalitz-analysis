#define analyzer_higgs_cxx

#include "analyzer_higgs.h"
#include <TH2.h>
#include <TStyle.h>
#include <TString.h>
#include <TMath.h>

using namespace std;

Float_t cut_vz = 24, cut_vd0 = 2, cut_vndof = 4;  //PV filter cuts
Float_t cut_Mll_low = 76.2, cut_Mll_high = 106.2;
 
UInt_t nEvents[20], nEventsPassHcal[20], nEventsPassEcal[20];
UInt_t nEventsPassNoiseFilter[6][20]; //1-Hcal, 2-Ecal, 3- Scraping, 4-CSCTight, 5-CSCLoose,   0-passed all

Bool_t isElectronSample = kFALSE, isMuonSample = kTRUE;

UInt_t verboseLvl = 1;
//TString sample("none");

void analyzer_higgs::Begin(TTree * /*tree*/)
{ 
  TString option = GetOption();
  histoFile = new TFile("hhhh.root", "RECREATE"); 
  //histoFile = new TFile(Form("hhistograms_%s.root", sample.Data()), "RECREATE"); 
  for(Int_t n=0; n<nCuts; n++)
    {
      met_phi[n]      = new TH1F(Form("met_phi_%i",n), "met_phi", 50, -TMath::Pi(), TMath::Pi()); 
      met_et[n]       = new TH1F(Form("met_et_%i",n), "met_et", 50, 0,350);
      met_over_qt[n]  = new TH1F(Form("met_over_qt_%i",n), "met_over_qt", 50, 0,4);
    
      met1_phi[n]     = new TH1F(Form("met1_phi_%i",n), "met1_phi", 50, -TMath::Pi(), TMath::Pi()); 
      met1_et[n]      = new TH1F(Form("met1_et_%i",n), "met1_et", 50, 0,350);
      met1_over_qt[n] = new TH1F(Form("met1_over_qt_%i",n), "met1_over_qt", 50, 0,4);

      met2_phi[n]     = new TH1F(Form("met2_phi_%i",n), "met2_phi", 50, -TMath::Pi(), TMath::Pi()); 
      met2_et[n]      = new TH1F(Form("met2_et_%i",n), "met2_et", 50, 0,350);
      met2_over_qt[n] = new TH1F(Form("met2_over_qt_%i",n), "met2_over_qt", 50, 0,4);

      met3_phi[n]     = new TH1F(Form("met3_phi_%i",n), "met3_phi", 50, -TMath::Pi(), TMath::Pi()); 
      met3_et[n]      = new TH1F(Form("met3_et_%i",n), "met3_et", 50, 0,350);
      met3_over_qt[n] = new TH1F(Form("met3_over_qt_%i",n), "met2_over_qt", 50, 0,4);

      met4_phi[n]     = new TH1F(Form("met4_phi_%i",n), "met4_phi", 50, -TMath::Pi(), TMath::Pi()); 
      met4_et[n]      = new TH1F(Form("met4_et_%i",n), "met4_et", 52, -14,350);
      met4_over_qt[n] = new TH1F(Form("met4_over_qt_%i",n), "met4_over_qt", 52, -0.16,4);
      met4_sig[n]     = new TH1F(Form("met4_sig_%i",n), "met4_sig", 52, -0.32,8);

      mu1_phi[n] = new TH1F(Form("mu1_phi_%i",n), "mu1_phi", 50, -TMath::Pi(), TMath::Pi()); 
      mu1_eta[n] = new TH1F(Form("mu1_eta_%i",n), "mu1_eta", 50, -2.6, 2.6); 
      mu1_pt[n]  = new TH1F(Form("mu1_pt_%i",n), "mu1_pt", 50, 0, 200); 

      mu2_phi[n] = new TH1F(Form("mu2_phi_%i",n), "mu2_phi", 50, -TMath::Pi(), TMath::Pi()); 
      mu2_eta[n] = new TH1F(Form("mu2_eta_%i",n), "mu2_eta", 50, -2.6, 2.6); 
      mu2_pt[n]  = new TH1F(Form("mu2_pt_%i",n), "mu2_pt", 50, 0, 200); 
      btag_hp[n] = new TH1F(Form("btag_hp_%i",n), "btag_hp", 50, -10, 10); 
    }

  ffout.open("./events_printout_andrey_final.txt",ofstream::out);
  ncout.open("./counts_for_tex.txt",ofstream::out);
  string sel;
  ifstream params;
  params.open("./params.txt"); params>>sel;
  if (sel == "muon")      { ncout<<"Mu selection!"<<endl; isMuonSample = kTRUE;  isElectronSample = kFALSE;}
  if (sel == "electron")  { ncout<<"El selection!"<<endl; isMuonSample = kFALSE; isElectronSample = kTRUE;}
  ncout<<endl;

 
  ffout.precision(4); ffout.setf(ios::fixed, ios::floatfield);
  for(Int_t i=0; i<nCuts; i++)
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
  //  if ( nEvents[0]>500) return kTRUE; 
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
  Float_t Mll=0, pTll=0;
  //--------------//
  //--- Muons------//
  //-------------//
  TLorentzVector diLepton(0.,0.,0.,0.), Lepton1(0.,0.,0.,0.), Lepton2(0.,0.,0.,0.);
 
  TCMuon *myMuon1=0, *myMuon2=0, *myMuon3=0, *myMuon4=0;
  Int_t muIndex[5] = {-1,-1,-1,-1,-1};

  Int_t muCount=0;
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
  Int_t eleCount=0;
  for (Int_t i = 0; i < recoElectrons->GetSize(); ++i) 
    {    
      TCElectron* thisElectron = (TCElectron*) recoElectrons->At(i);     
      TVector3* eleVertex = new TVector3(thisElectron->Vtx());
      TVector3* eleMomentum = new TVector3(thisElectron->p4().Vect());
      
      double dxyVtx = dxy(mainPrimaryVertex, eleVertex, eleMomentum);
      double dzVtx  = dz(mainPrimaryVertex, eleVertex, eleMomentum);
      
      Float_t  eleISO = (thisElectron->trkIso() + thisElectron->emIso() + thisElectron->hadIso() - rhoFactor*TMath::Pi()*0.09)/thisElectron->pt(); 
      Float_t  eleISObarrel = 0;
      //if((thisElectron->emIso()-1)>0) eleISObarrel = eleISO;
      //else eleISObarrel = (thisElectron->trkIso() + thisElectron->hadIso() - rhoFactor*TMath::Pi()*0.09)/thisElectron->pt(); 
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
      if (eleCount==0 && myMuon1->charge()*myMuon2->charge()==-1 && Mll>20 && Mll<500)	 
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
      if (muCount==0 && myElectron1->charge()*myElectron2->charge()==-1 && Mll>20 && Mll<500)	 
	{
	  //cout<<"preselection?!"<<endl;
	  passPreSelection= kTRUE;
	  if (Mll>cut_Mll_low && Mll<cut_Mll_high) passZpeak = kTRUE;
	}

    }

  Float_t MET=0, MET_phi=0;
  Float_t MET1=0, MET1_phi=0;
  Float_t MT = 0;
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
  //if (muCount==2) cout<<"MT: "<<MT<<"    "<<ZptXY.X()<<" "<<ZptXY.Y()<<"    "<<METXY.X()<<" "<<METXY.Y()<<endl;

  Float_t METqt = MET/qT;
  met_over_qt[0] -> Fill(METqt);

  met_et[0] -> Fill(MET);
  met_phi[0] -> Fill(MET_phi);

  met1_et[0] -> Fill(MET1);
  met1_phi[0] -> Fill(MET1_phi);

  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo) 
    {
      met2_et[0]      -> Fill(MET);
      met2_phi[0]     -> Fill(MET_phi);

      met3_over_qt[0] -> Fill(projMETqt);
      met3_et[0]      -> Fill(projMET);

      met4_over_qt[0] -> Fill(puCorrMETqt);
      met4_et[0]      -> Fill(puCorrMET);
      met4_sig[0]     -> Fill(puSigMET);
    }
  else
    if(verboseLvl>1) 
      nout[0]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<MET<<"\t* "
	     <<isNoiseHcal<<"\t* "<<isDeadEcalCluster<<"\t* "<<isScraping<<"\t* "<<isCSCTightHalo<<"\t* "<<endl;  

  if (!passPreSelection) return kTRUE;
  
  nEvents[1]++;
  if(!isNoiseHcal)        nEventsPassNoiseFilter[1][1]++; 
  if(!isDeadEcalCluster)  nEventsPassNoiseFilter[2][1]++; 
  if(!isScraping)         nEventsPassNoiseFilter[3][1]++;
  if(!isCSCTightHalo)     nEventsPassNoiseFilter[4][1]++;
  if(!isCSCLooseHalo)     nEventsPassNoiseFilter[5][1]++;
  if(!isNoiseHcal && !isDeadEcalCluster  && !isScraping && !isCSCTightHalo) 
    {
      nEventsPassNoiseFilter[0][1]++; 
      met2_et[1]         -> Fill(MET);
      met2_phi[1]        -> Fill(MET_phi);

      met3_over_qt[1] -> Fill(projMETqt);
      met3_et[1]         -> Fill(projMET);

      met4_over_qt[1] -> Fill(puCorrMETqt);
      met4_et[1]         -> Fill(puCorrMET);
      met4_sig[1]        -> Fill(puSigMET);
    }
  else
    if(verboseLvl>0) 
      nout[1]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<MET<<"\t* "
	     <<isNoiseHcal<<"\t* "<<isDeadEcalCluster<<"\t* "<<isScraping<<"\t* "<<isCSCTightHalo<<"\t* "<<endl;  
  

  met_et[1]  -> Fill(MET);
  met_phi[1] -> Fill(MET_phi);

  met1_et[1]  -> Fill(MET1);
  met1_phi[1] -> Fill(MET1_phi);

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
      met2_et[2]      -> Fill(MET);
      met2_phi[2]     -> Fill(MET_phi);
      met3_over_qt[2] -> Fill(projMETqt);
      met3_et[2]      -> Fill(projMET);

      met4_over_qt[2] -> Fill(puCorrMETqt);
      met4_et[2]      -> Fill(puCorrMET);
      met4_sig[2]     -> Fill(puSigMET);
    }


  met_et[2]  -> Fill(MET);
  met_phi[2] -> Fill(MET_phi);

  met1_et[2]  -> Fill(MET1);
  met1_phi[2] -> Fill(MET1_phi);

  if (isMuonSample)
    {
      mu1_eta[2] -> Fill(myMuon1->eta());
      mu1_phi[2] -> Fill(myMuon1->phi());
      mu1_pt[2]  -> Fill(myMuon1->pt());

      mu2_eta[2] -> Fill(myMuon2->eta());
      mu2_phi[2] -> Fill(myMuon2->phi());
      mu2_pt[2]  -> Fill(myMuon2->pt());

      if (myMuon1->pt() < 30) return kTRUE;
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
      if (myElectron1->pt() < 30) return kTRUE;
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
  met_et[4]  -> Fill(MET);
  met_phi[4] -> Fill(MET_phi);

  met1_et[4]  -> Fill(MET1);
  met1_phi[4] -> Fill(MET1_phi);

  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo)
    { 
      met2_et[4]         -> Fill(MET);
      met2_phi[4]        -> Fill(MET_phi);
      met3_over_qt[4] -> Fill(projMETqt);
      met3_et[4]         -> Fill(projMET);
      
      met4_over_qt[4] -> Fill(puCorrMETqt);
      met4_et[4]         -> Fill(puCorrMET);
      met4_sig[4]        -> Fill(puSigMET);
    }

  if(verboseLvl>1) 
    fout[4]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<Mll<<"\t* "<<MT<<"\t* "<<MET<<"\t* "<<pTll<<"\t* "<<rhoFactor<<"\t*"<<endl;  

  Int_t jetCount=0;
  for (Int_t i = 0; i < recoJets->GetSize(); ++i) 
    {    
      TCJet* thisJet = (TCJet*) recoJets->At(i);     
      if (fabs(thisJet->P4(7).Eta()) <= 2.5) 
	{
	  //do anti-b tag
          if (thisJet->P4(7).Pt() > 30 && thisJet->BDiscrTrkCountHiEff() > 2) return kTRUE; 
	}
      jetCount++;
    }
  
  nEvents[5]++;

  if(!isNoiseHcal)        nEventsPassNoiseFilter[1][5]++; 
  if(!isDeadEcalCluster)  nEventsPassNoiseFilter[2][5]++; 
  if(!isScraping)         nEventsPassNoiseFilter[3][5]++;
  if(!isCSCTightHalo)     nEventsPassNoiseFilter[4][5]++;
  if(!isCSCLooseHalo)     nEventsPassNoiseFilter[5][5]++;
  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo) 
    {
      nEventsPassNoiseFilter[0][5]++; 
      met2_et[5]         -> Fill(MET);
      met2_phi[5]        -> Fill(MET_phi);
      met3_over_qt[5] -> Fill(projMETqt);
      met3_et[5]         -> Fill(projMET);
      
      met4_over_qt[5] -> Fill(puCorrMETqt);
      met4_et[5]         -> Fill(puCorrMET);
      met4_sig[5]        -> Fill(puSigMET);
    }

  met_et[5]      -> Fill(MET);
  met_phi[5]     -> Fill(MET_phi);
  met_over_qt[5] -> Fill(METqt);

  met1_et[5]      -> Fill(MET1);
  met1_phi[5]     -> Fill(MET1_phi);
  met1_over_qt[5] -> Fill(METqt);

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

  if (pTll<50) return kTRUE;
  nEvents[6]++;

  if(!isNoiseHcal)        nEventsPassNoiseFilter[1][6]++; 
  if(!isDeadEcalCluster)  nEventsPassNoiseFilter[2][6]++; 
  if(!isScraping)         nEventsPassNoiseFilter[3][6]++;
  if(!isCSCTightHalo)     nEventsPassNoiseFilter[4][6]++;
  if(!isCSCLooseHalo)     nEventsPassNoiseFilter[5][6]++;
  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo) 
    {
      nEventsPassNoiseFilter[0][6]++; 
      met2_et[6]         -> Fill(MET);
      met2_phi[6]        -> Fill(MET_phi);
      met3_over_qt[6] -> Fill(projMETqt);
      met3_et[6]         -> Fill(projMET);
      
      met4_over_qt[6] -> Fill(puCorrMETqt);
      met4_et[6]         -> Fill(puCorrMET);
      met4_sig[6]        -> Fill(puSigMET);
    }
  else
    if(verboseLvl>0) 
      nout[6]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<MET<<"\t* "
	     <<isNoiseHcal<<"\t* "<<isDeadEcalCluster<<"\t* "<<isScraping<<"\t* "<<isCSCTightHalo<<"\t* "<<endl;  
  

  met_et[6]      -> Fill(MET);
  met_phi[6]     -> Fill(MET_phi);
  met_over_qt[6] -> Fill(METqt);

  met1_et[6]      -> Fill(MET1);
  met1_phi[6]     -> Fill(MET1_phi);
  met1_over_qt[6] -> Fill(METqt);

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
    fout[6]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<Mll<<"\t* "<<MT<<"\t* "<<MET<<"\t* "<<pTll<<"\t* "<<rhoFactor<<"\t*"<<endl;  

  if (MET<80) return kTRUE;
  nEvents[7]++;
  if(!isNoiseHcal)        nEventsPassNoiseFilter[1][7]++; 
  if(!isDeadEcalCluster)  nEventsPassNoiseFilter[2][7]++; 
  if(!isScraping)         nEventsPassNoiseFilter[3][7]++;
  if(!isCSCTightHalo)     nEventsPassNoiseFilter[4][7]++;
  if(!isCSCLooseHalo)     nEventsPassNoiseFilter[5][7]++;
  if(!isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo) 
    {
      nEventsPassNoiseFilter[0][7]++; 
      met2_et[7]         -> Fill(MET);
      met2_phi[7]        -> Fill(MET_phi);
      met3_over_qt[7] -> Fill(projMETqt);
      met3_et[7]         -> Fill(projMET);
      
      met4_over_qt[7] -> Fill(puCorrMETqt);
      met4_et[7]         -> Fill(puCorrMET);
      met4_sig[7]        -> Fill(puSigMET);
    }
  else
    if(verboseLvl>0) 
      nout[7]<<"* "<<runNumber<<"\t* "<<eventNumber<<"  \t* "<<lumiSection<<"\t* "<<MET<<"\t* "
	     <<isNoiseHcal<<"\t* "<<isDeadEcalCluster<<"\t* "<<isScraping<<"\t* "<<isCSCTightHalo<<"\t* "<<endl;  
  
  met_et[7]  -> Fill(MET);
  met_phi[7] -> Fill(MET_phi);

  met1_et[7]  -> Fill(MET1);
  met1_phi[7] -> Fill(MET1_phi);

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

  for(Int_t i=0; i<nCuts; i++) {nout[i].close();  fout[i].close();}
  ffout.close();

  ncout.precision(3); cout.setf(ios::fixed, ios::floatfield);

  ncout<<"0 Total:  \n| "<<nEvents[0]<<" |\t"<<(Float_t)nEvents[0]/nEvents[0]<<" |"<<(Float_t)nEvents[0]/nEvents[0]<<" |"<<endl;
  ncout<<"1 presel       \n| "<<nEvents[1]<<" |\t"<<(Float_t)nEvents[1]/nEvents[0]<<" |"<<(Float_t)nEvents[1]/nEvents[0]<<" |"<<endl;
  ncout<<"2 Mll   \n| "<<nEvents[2]<<" |\t"<<(Float_t)nEvents[2]/nEvents[0]<<" |"<<(Float_t)nEvents[2]/nEvents[1]<<" |"<<endl;
  ncout<<"3 pt1     \n| "<<nEvents[3]<<" |\t"<<(Float_t)nEvents[3]/nEvents[0]<<" |"<<(Float_t)nEvents[3]/nEvents[2]<<" |"<<endl;
  ncout<<"4 pt2   \n| "<<nEvents[4]<<" |\t"<<(Float_t)nEvents[4]/nEvents[0]<<" |"<<(Float_t)nEvents[4]/nEvents[3]<<" |"<<endl;
  ncout<<"5 anti-b  \n| "<<nEvents[5]<<" |\t"<<(Float_t)nEvents[5]/nEvents[0]<<" |"<<(Float_t)nEvents[5]/nEvents[4]<<" |"<<endl;
  ncout<<"6 pTll  \n| "<<nEvents[6]<<" |\t"<<(Float_t)nEvents[6]/nEvents[0]<<" |"<<(Float_t)nEvents[6]/nEvents[5]<<" |"<<endl;
  ncout<<"7 pfMet  \n| "<<nEvents[7]<<" |\t"<<(Float_t)nEvents[7]/nEvents[0]<<" |"<<(Float_t)nEvents[7]/nEvents[6]<<" |"<<endl;
  ncout<<"8      \n| "<<nEvents[8]<<" |\t"<<(Float_t)nEvents[8]/nEvents[0]<<" |"<<(Float_t)nEvents[8]/nEvents[7]<<" |"<<endl;
  // ncout<<"9      \n| "<<nEvents[9]<<" |\t"<<(Float_t)nEvents[9]/nEvents[0]<<" |"<<(Float_t)nEvents[9]/nEvents[8]<<" |"<<endl;
  // ncout<<"10     \n| "<<nEvents[10]<<" |\t"<<(Float_t)nEvents[10]/nEvents[0]<<" |"<<(Float_t)nEvents[10]/nEvents[9]<<" |"<<endl;
  
  ncout<<endl;
  ncout<<" & total   &  passHcal  & passEcal  &  Scraping  & CSCTight & CSCLoose & & passAll & frac \\\\ \n \\hline"<<endl;
  for (Int_t i=0; i<nCuts; i++)
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
  histoFile->cd();
  histoFile->Write();
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

double analyzer_higgs::higgsMT(Float_t pTll, Float_t Mll, Float_t MET, TVector2 ZptXY, TVector2 METXY)
{
  double MT = sqrt( pow(sqrt(pTll*pTll + Mll*Mll) + sqrt(MET*MET + Mll*Mll), 2) - (ZptXY + METXY).Mod2());
  return MT;
}

double analyzer_higgs::projectedMET(Float_t MET, Float_t MET_phi, TLorentzVector Lepton1, TLorentzVector Lepton2)
{
  //Float_t deltaPhiMin = TMath::Min(deltaPhi(Lepton1.Phi(), MET_phi), deltaPhi(Lepton2.Phi(), MET_phi));
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
