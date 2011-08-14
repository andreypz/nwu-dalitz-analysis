#define analyzer_higgs_cxx

#include "analyzer_higgs.h"
#include <TH2.h>
#include <TStyle.h>
#include <TString.h>
#include <TMath.h>
#include "TF1.h"

#include <algorithm>
using namespace std;

Float_t cut_vz = 24, cut_vd0 = 2, cut_vndof = 4;  //PV filter cuts
Float_t cut_Mll_low = 76.2, cut_Mll_high = 106.2;

UInt_t nEvents[20], nEventsPassHcal[20], nEventsPassEcal[20];
UInt_t nEventsPassNoiseFilter[6][20]; //1-Hcal, 2-Ecal, 3- Scraping, 4-CSCTight, 5-CSCLoose,   0-passed all

Bool_t isElectronSample = kFALSE, isMuonSample = kFALSE, isPhotonSample = kTRUE;

UInt_t verboseLvl = 1;
//TString sample("none");


// sorting of photons
bool higherPhotonPt(TCPhoton* p1, TCPhoton* p2) {
	return p1->pt() > p2->pt();
}

Long64_t eventCounter = 0;


void analyzer_higgs::Begin(TTree * /*tree*/)
{ 
	// read weight histograms for the photons
  //	TFile* weightFile = new TFile("histo/photonWeight_allr1_noBgSub_110713.root");
        TFile* weightFile = new TFile("histo/photonWeight_110714-1.root");
	h1_photonWeightEe = (TH1F*) weightFile->Get("h1_eeScale")->Clone();
	h1_photonWeightMumu = (TH1F*) weightFile->Get("h1_mumuScale")->Clone();
	h1_photonWeightLl = (TH1F*) weightFile->Get("h1_llScale")->Clone(); 
	//	weightFile->Close();
	h1_photonMassFromEe = (TH1F*) weightFile->Get("h1_photonMassFromEe")->Clone();
	h1_photonMassFromMumu = (TH1F*) weightFile->Get("h1_photonMassFromMumu")->Clone();


	TString option = GetOption();
	histoFile = new TFile("hhhhnew2.root", "RECREATE"); 
	histoFile->cd();

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

	h_hlt_fired = new TH1F("hlt_fired","hlt_fired",50, -0.5, 49.5);
	h_hlt_fired_prescale = new TH2F("hlt_fired_prescale","hlt_fired_prescale",50, -0.5, 49.5, 10001, -0.5, 10000.5);
	h_hlt_fired2 = new TH1F("hlt_fired_1gamma","hlt_fired_1gamma",50, -0.5, 49.5);
	h_hlt_fired_prescale2 = new TH2F("hlt_fired_prescale_1gamma","hlt_fired_prescale_1gamma",50, -0.5, 49.5, 10001, -0.5, 10000.5);


	// all events
	h_pho_pt = new TH1F("photon_pt","photon_pt", 40, 0., 400.);
	h_pho_eta = new TH1F("photon_eta", "photon_eta", 50, -5, 5);
	h_pho_phi  = new TH1F("photon_phi", "photon_phi", 50, -TMath::Pi(), TMath::Pi());
	h_pho_mass = new TH1F("photon_mass","photon_mass", 100, 0., 300.);
	h_phow_pt = new TH1F("photon_weighted_pt","photon_weighted_pt",40, 0., 400.);
	h_phow_eta = new TH1F("photon_weighted_eta", "photon_weighted_eta", 50, -5, 5);
	h_phow_phi  = new TH1F("photon_weighted_phi", "photon_weighted_phi", 50, -TMath::Pi(), TMath::Pi());
	h_phow_mass = new TH1F("photon_weighted_mass","photon_weighted_mass",100, 0., 300.);


	// split sample in two based on odd/even entries
	h_pho1_pt = new TH1F("photon1_pt","photon1_pt", 40, 0., 400.);
	h_pho1_eta = new TH1F("photon1_eta", "photon1_eta", 50, -5, 5);
	h_pho1_phi = new TH1F("photon1_phi", "photon1_phi", 50, -TMath::Pi(), TMath::Pi());
	h_pho1_mass = new TH1F("photon1_mass","photon1_mass", 100, 0., 300.);

	h_phow1_pt = new TH1F("photon1_weighted_pt","photon1_weighted_pt",40, 0., 400.);
	h_phow1_eta = new TH1F("photon1_weighted_eta", "photon1_weighted_eta", 50, -5, 5);
	h_phow1_phi = new TH1F("photon1_weighted_phi", "photon1_weighted_phi", 50, -TMath::Pi(), TMath::Pi());
	h_phow1_mass = new TH1F("photon1_weighted_mass","photon1_weighted_mass",100, 0., 300.);


	h_pho2_pt = new TH1F("photon2_pt","photon2_pt", 40, 0., 400.);
	h_pho2_eta = new TH1F("photon2_eta", "photon2_eta", 50, -5, 5);
	h_pho2_phi = new TH1F("photon2_phi", "photon2_phi", 50, -TMath::Pi(), TMath::Pi());
	h_pho2_mass = new TH1F("photon2_mass","photon2_mass", 100, 0., 300.);

	h_phow2_pt = new TH1F("photon2_weighted_pt","photon2_weighted_pt",40, 0., 400.);
	h_phow2_eta = new TH1F("photon2_weighted_eta", "photon2_weighted_eta", 50, -5, 5);
	h_phow2_phi = new TH1F("photon2_weighted_phi", "photon2_weighted_phi", 50, -TMath::Pi(), TMath::Pi());
	h_phow2_mass = new TH1F("photon2_weighted_mass","photon2_weighted_mass",100, 0., 300.);


	for (int i=0; i<12; ++i) {
		h_phow1_mt[i] = new TH1F(Form("photon1_weighted_mt_%i",i) ,Form("photon1_weighted_mt_%i",i), 50, 100., 600.);
		h_phow1_met[i] = new TH1F(Form("photon1_weighted_met_%i",i),Form("photon1_weighted_met_%i",i), 40, 0., 400.);
		h_phow1_metOverQt[i] = new TH1F(Form("photon1_weighted_metOverQt_%i",i),Form("photon1_weighted_metOverQt_%i",i), 40, 0., 4.);
		h_phow1_nJets[i] = new TH1F(Form("photon1_weighted_nJets_%i",i),Form("photon1_weighted_nJets_%i",i), 20, 0., 20.);
		h_phow1_projMet[i] = new TH1F(Form("photon1_weighted_projMet_%i",i),Form("photon1_weighted_projMet_%i",i), 40, 0., 400.);
		h_phow1_projMetOverQt[i] = new TH1F(Form("photon1_weighted_projMetOverQt_%i",i),Form("photon1_weighted_projMetOverQt_%i",i), 40, 0., 4.);

		h_phow2_mt[i] = new TH1F(Form("photon2_weighted_mt_%i",i) ,Form("photon2_weighted_mt_%i",i), 50, 100., 600.);
		h_phow2_met[i] = new TH1F(Form("photon2_weighted_met_%i",i),Form("photon2_weighted_met_%i",i), 40, 0., 400.);
		h_phow2_metOverQt[i] = new TH1F(Form("photon2_weighted_metOverQt_%i",i),Form("photon2_weighted_metOverQt_%i",i), 40, 0., 4.);
		h_phow2_nJets[i] = new TH1F(Form("photon2_weighted_nJets_%i",i),Form("photon2_weighted_nJets_%i",i), 20, 0., 20.);
		h_phow2_projMet[i] = new TH1F(Form("photon2_weighted_projMet_%i",i),Form("photon2_weighted_projMet_%i",i), 40, 0., 400.);
		h_phow2_projMetOverQt[i] = new TH1F(Form("photon2_weighted_projMetOverQt_%i",i),Form("photon2_weighted_projMetOverQt_%i",i), 40, 0., 4.);

		h_phow1_mt[i]->Sumw2();
		h_phow1_met[i]->Sumw2();
		h_phow1_metOverQt[i]->Sumw2();
		h_phow1_nJets[i]->Sumw2();
		h_phow1_projMet[i]->Sumw2();
		h_phow1_projMetOverQt[i]->Sumw2();

		h_phow2_mt[i]->Sumw2();
		h_phow2_met[i]->Sumw2();
		h_phow2_metOverQt[i]->Sumw2();
		h_phow2_nJets[i]->Sumw2();
		h_phow2_projMet[i]->Sumw2();
		h_phow2_projMetOverQt[i]->Sumw2();	
	}


	h_pho_pt->Sumw2();
	h_pho_eta->Sumw2();
	h_pho_phi->Sumw2();
	h_pho_mass->Sumw2();
	h_phow_pt->Sumw2();
	h_phow_eta->Sumw2();
	h_phow_phi->Sumw2();
	h_phow_mass->Sumw2();

	h_pho1_pt->Sumw2();
	h_pho1_eta->Sumw2();
	h_pho1_phi->Sumw2();
	h_pho1_mass->Sumw2();
	h_phow1_pt->Sumw2();
	h_phow1_eta->Sumw2();
	h_phow1_phi->Sumw2();
	h_phow1_mass->Sumw2();

	h_pho2_pt->Sumw2();
	h_pho2_eta->Sumw2();
	h_pho2_phi->Sumw2();
	h_pho2_mass->Sumw2();
	h_phow2_pt->Sumw2();
	h_phow2_eta->Sumw2();
	h_phow2_phi->Sumw2();
	h_phow2_mass->Sumw2();




	ffout.open("./events_printout_andrey_final.txt",ofstream::out);
	ncout.open("./counts_for_tex.txt",ofstream::out);
	string sel;
	ifstream params;
	params.open("./params.txt"); params>>sel;
	if (sel == "muon")      { ncout<<"Mu selection!"<<endl; isMuonSample = kTRUE;  isElectronSample = kFALSE; isPhotonSample = kFALSE;}
	if (sel == "electron")  { ncout<<"El selection!"<<endl; isMuonSample = kFALSE; isElectronSample = kTRUE; isPhotonSample = kFALSE;}
	if (sel == "photon") { ncout<<"Photon selection!"<<endl; isMuonSample = kFALSE; isElectronSample = kFALSE; isPhotonSample = kTRUE;}

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

	++eventCounter;

	if (eventCounter%100000 == 0) {
		cout << "Entry: " << eventCounter << "    Run: " << runNumber << "  Event: " << eventNumber << endl;
	}
// 	if (entry>10000) return kFALSE;


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


	//-------------------//
	//--- Trigger Selection------//
	//-------------------//

	Bool_t passSinglePhotonHLT = kFALSE;  // for use in gamma+jet(s) selection

	//  Bool_t = passHLT = kFALSE;

	// 	  cout << "trigSize: " << hltTriggerResults->GetSize() << endl;

	for (Int_t i = 0; i < hltTriggerResults->GetSize(); i++)
	{
		TCTrigger* thisTrigger = (TCTrigger*) hltTriggerResults->At(i);

		//cout << "ind: " << i << endl;
		//cout << "valid :" << thisTrigger->isValid() << endl;
		//cout << "fired :" << thisTrigger->hasFired() << endl << endl;

		Int_t trighasfired = 0;
		Int_t trigprescale = 0;
		if (thisTrigger->isValid() > 0)
		{
			trighasfired = thisTrigger->hasFired();
			trigprescale = thisTrigger->prescale();
			if (trighasfired > 0)
			{
				// plot triggers 1D
				h_hlt_fired->Fill(i);
				// plot 2D fired vs. prescale
				h_hlt_fired_prescale->Fill(i,trigprescale);
				// make here the selection of a particular trigger
				//	      if (i==1) passHLT = kTRUE;
			}
		}
	}
	//  if (!passHLT) return kTRUE;

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


	//-------------------//
	//--- Photons------//
	//-------------------//

	vector<TCPhoton*> goodPhotons; // pass all cuts

	TCPhoton *myPhoton1=0, *myPhoton2=0;

	Int_t phoIndex[3] = {-1,-1,-1};

	//if (recoPhotons->GetSize()>1)  cout<<"N ele: "<<  recoPhotons->GetSize()<<endl;
	Int_t phoCount=0;
	for (Int_t i = 0; i < recoPhotons->GetSize(); ++i)
	{
		TCPhoton* thisPhoton = (TCPhoton*) recoPhotons->At(i);

		Float_t photonPT = thisPhoton->p4().Pt();
		if ( (thisPhoton->emIso() < (4.2 + 0.006*photonPT)) &&
			(thisPhoton->hadIso() < (2.2 + 0.0025*photonPT)) &&
			(thisPhoton->trkIso() < (2.0 + 0.001*photonPT)) &&
			(thisPhoton->hadOverEm() < 0.05 ) &&
			(thisPhoton->sigmaIetaIeta() > 0.001 ) &&
			(thisPhoton->sigmaIetaIeta() < 0.013 ) &&
			(thisPhoton->etaSupercluster() < 1.4442 ) &&
			(thisPhoton->trackVeto() == 0 )
			)
		{
			goodPhotons.push_back(thisPhoton);
			if (phoIndex[1]==-1) phoIndex[1] = i;  // Radek's stuff
			if (phoIndex[2]==-1 && phoIndex[1]!=i) phoIndex[2] = i;

			phoCount++;
		}
	}
	if (phoCount>0)
	{
		myPhoton1 = (TCPhoton*) recoPhotons->At(phoIndex[1]);

		for (Int_t i = 0; i < hltTriggerResults->GetSize(); i++)
		{
			TCTrigger* thisTrigger = (TCTrigger*) hltTriggerResults->At(i);
			Int_t trighasfired = 0;
			Int_t trigprescale = 0;
			if (thisTrigger->isValid())
			{
				trighasfired = thisTrigger->hasFired();
				trigprescale = thisTrigger->prescale();
				if (trighasfired > 0)
				{
					// plot triggers 1D
					h_hlt_fired2->Fill(i);
					// plot 2D fired vs. prescale
					h_hlt_fired_prescale2->Fill(i,trigprescale);
				}
			}
		}
	}

	// sort the photons for use in the analysis
	sort(goodPhotons.begin(), goodPhotons.end(), higherPhotonPt);


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
			//	  passPreSelection= kTRUE;
			//	  if (Mll>cut_Mll_low && Mll<cut_Mll_high) passZpeak = kTRUE;
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
			//	  passPreSelection= kTRUE;
			//	  if (Mll>cut_Mll_low && Mll<cut_Mll_high) passZpeak = kTRUE;
		}

	}

	if (isPhotonSample)
	{
		passPreSelection = kTRUE;
		passZpeak = kTRUE;

		// check the HLT: select only photons on the single trigger
		for (Int_t i = 0; i < hltTriggerResults->GetSize(); i++) {

			// ignore all triggers exept single photon
			// this is from radek's indexing:   
			// HLT_Photon20_CaloIdVL_IsoL_v1-6     0-  5
			// HLT_Photon30_CaloIdVL_v1-6         60- 65
			// HLT_Photon30_CaloIdVL_IsoL_v1-6    66- 71
			// HLT_Photon50_CaloIdVL_IsoL_v1-6    84- 89
			// HLT_Photon75_CaloIdVL_v1-6        120-125
			// HLT_Photon75_CaloIdVL_IsoL_v1-6   126-131
		        // HLT_Photon90_CaloIdVL_v1-6        138-143
   		        // HLT_Photon90_CaloIdVL_IsoL_v1-6   144-149

			if (!(i>=0 && i<=5) &&
				!(i>=60 && i<=71) &&
				!(i>=84 && i<=89) &&
				!(i>=120 && i<=131) &&
			        !(i>=138 && i<=149)
				)
				continue;

			TCTrigger* thisTrigger = (TCTrigger*) hltTriggerResults->At(i);

			if (thisTrigger->isValid() && thisTrigger->hasFired()) {
				passSinglePhotonHLT = kTRUE;
//				cout << "passSinglePhotonHLT" << endl;
				break;
			}

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



	// define some local cuts that do not impact the general analysis flow for other 
	// final states (quick fix to selection etc)
	// this loop is done elsewhere, but the info is not accessible from here
	// get also number of jets (AA)

	// Note: when counting the jets one should exclude all e,mu, etc...
	// However, we consider photon events with e,mu=0 and the structure of the code
	// does not allow consistent general implementation.
	// In addition, one should apply some jet ID... which I'm not sure this analysis does. AA)

	Bool_t passBJetVeto =  kTRUE;
	Bool_t passDPhiJetMet = kTRUE;
	Int_t numJets = 0;  // in specified eta region abovce a threshold
	for (Int_t i = 0; i < recoJets->GetSize(); ++i)  {    
		TCJet* thisJet = (TCJet*) recoJets->At(i);     

		Bool_t matchedToPhoton = kFALSE;
		for (vector<TCPhoton*>::iterator it = goodPhotons.begin(); it!=goodPhotons.end(); ++it) {
			if ( thisJet->P4().DeltaR((*it)->p4())<0.4 ) {
				matchedToPhoton = kTRUE;
				break;
			}
		}

		if (thisJet->P4(7).Pt() > 20)  { // still using uncorrected met....
			if (fabs(TVector2::Phi_mpi_pi(MET_phi - thisJet->P4(7).Phi())) < 0.3) passDPhiJetMet = kFALSE;
		}

		if (!matchedToPhoton && fabs(thisJet->P4(7).Eta()) <= 2.4 && thisJet->P4(7).Pt() > 30)  {
			++numJets;
			//do anti-b tag, thresholds for jet counting and b tagging are the same at the moment
			if (thisJet->BDiscrTrkCountHiEff() > 2) {
				passBJetVeto =  kFALSE; 
			}
		}

	} // loop over jets


	if (isPhotonSample && phoCount>0 && 
		(passSinglePhotonHLT || !isRealData)
		&& eleCount==0 && muCount == 0  // to be consistent with the extra electron/muon veto 
		&& !isNoiseHcal && !isDeadEcalCluster && !isScraping && !isCSCTightHalo)  // kill noise events -> they should have been applied globally!!!
        //		&& passBJetVeto) // seems to be a requirement before making the spectrum... split histogremas later for with/withot this veto
	{

		// assigned mass to the photon according to
		// the selected ll distribution

		TCPhoton* leadPhoton = *goodPhotons.begin();

		Float_t photonPt = leadPhoton->p4().Pt();
		Float_t photonEta = leadPhoton->p4().Eta();
		Float_t photonPhi = leadPhoton->p4().Phi();

		Float_t photonMassFromEe = h1_photonMassFromEe->GetRandom(); 
		Float_t photonMassFromMumu = h1_photonMassFromMumu->GetRandom(); 


		// this is the Mt definition from the original code (see below) - just replaced the mass.
		// There does not seem to be a need to redefine the energy of the p4 of the photoen
		// to account for the mass as it is not used anywhere yet...

		// these have to be moved to the "if" blocks for ee/mum -> keep them here till the functional form is confirmed
		// and synchronised with Andrey's definitions
		Float_t Mt_ee = TMath::Power( (sqrt(TMath::Power(leadPhoton->p4().Pt(),2) + TMath::Power(photonMassFromEe,2) ) 
			+ sqrt(TMath::Power(MET,2) + TMath::Power(photonMassFromEe,2)) ) , 2);
		Mt_ee = TMath::Sqrt(Mt_ee);

		Float_t Mt_mumu = TMath::Power( (sqrt(TMath::Power(leadPhoton->p4().Pt(),2) + TMath::Power(photonMassFromMumu,2) ) 
			+ sqrt(TMath::Power(MET,2) + TMath::Power(photonMassFromMumu,2)) ) , 2);
		Mt_mumu = TMath::Sqrt(Mt_mumu);


		// note that Andrey used pfMet for MET, and pfMet Corr for MET1 - for now we are comparing without corrections...
		// MET and MET_phi are previously assigned..

		Float_t metOverQt = MET/photonPt;

		// new definitions of "projected met"... I'm not happy neither with this, nor the previous definition, but people want it...

		Float_t newProjMet = MET;
		Float_t deltaPhiPhoMet = TVector2::Phi_mpi_pi(MET_phi - photonPhi);
		if (fabs(deltaPhiPhoMet) < TMath::PiOver2()) newProjMet *= TMath::Sin(deltaPhiPhoMet);

		Float_t newProjMetOverQt = newProjMet/photonPt;


		// cumulative cuts
		Bool_t passCuts[12] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};

		if (photonPt>25 && passDPhiJetMet) passCuts[0] =  kTRUE;
		if (passBJetVeto && passCuts[0]) passCuts[1] =  kTRUE; 
		if (MET>70 && passCuts[1])  passCuts[2] =  kTRUE;
		if (newProjMet>70 && passCuts[1])  passCuts[3] =  kTRUE; // this is not on top of the previous cut, but instead of

		// new cuts as per 15.07.11 (Andrey's request)
		if (MET>50 && passCuts[1]) passCuts[4] = kTRUE;
                if (newProjMet>50 && passCuts[1]) passCuts[5] = kTRUE;
		// cuts 6 and 7 see subsamples

		Double_t y2d = MET/photonPt;
		Double_t y2d_up = 0.0053333 * MET + 0.36666;
		Double_t y2d_down = 0.0046666 * MET + 0.033333;
                Double_t y2dp = newProjMet/photonPt;
                Double_t y2dp_up = 0.0053333 * newProjMet + 0.36666;
                Double_t y2dp_down = 0.0046666 * newProjMet + 0.033333;
		if (passCuts[1] && (MET>100 && MET<400 && y2d>y2d_down && y2d<y2d_up)) passCuts[8] = kTRUE;
                if (passCuts[1] && (newProjMet>100 && newProjMet<400 && y2dp>y2dp_down && y2dp<y2dp_up)) passCuts[9] = kTRUE;

		if (eventCounter%3 == 0) {  // this will be the normalization sample

			// plot photon pT and MET, MT
			h_pho_pt->Fill(photonPt);
			// Transverse Higgs mass
			//      Float_t Mt = TMath::Power( (sqrt(TMath::Power(leadPhoton->p4().Pt(),2) + TMath::Power(leadPhoton->p4().M(),2) ) 
			//				  + sqrt(TMath::Power(MET,2) + TMath::Power(leadPhoton->p4().M(),2)) ) , 2)
			//	- ();
			//      h_pho_Mt->Fill();
//			h_pho_mass->Fill(photonMass);
			h_pho_phi->Fill(photonPhi);
			h_pho_eta->Fill(photonEta);
		}

		// for split-sample for e/mu for now (as far as the wighted distributions are concerned)

		else if (eventCounter%3 == 1) {

		        if (MET>112 && Mt_ee > 292 && passCuts[1]) passCuts[6] = kTRUE;
		        if (newProjMet>112 && Mt_ee > 292 && passCuts[1]) passCuts[7] = kTRUE;

			h_pho1_pt->Fill(photonPt);
			h_pho1_mass->Fill(photonMassFromEe);
			h_pho1_phi->Fill(photonPhi);
			h_pho1_eta->Fill(photonEta);

			Float_t weightEe = photonWeightEe(photonPt);

			h_phow1_pt->Fill(photonPt, weightEe);
			h_phow1_mass->Fill(photonMassFromEe, weightEe);
			h_phow1_phi->Fill(photonPhi, weightEe);
			h_phow1_eta->Fill(photonEta, weightEe);

			for (int c=0; c<12; ++c) {
				if (passCuts[c]) {
					h_phow1_met[c]->Fill(MET, weightEe);
					h_phow1_mt[c]->Fill(Mt_ee, weightEe);
					h_phow1_metOverQt[c]->Fill(metOverQt, weightEe);
					h_phow1_nJets[c]->Fill(numJets, weightEe);
					h_phow1_projMet[c]->Fill(newProjMet, weightEe);
					h_phow1_projMetOverQt[c]->Fill(newProjMetOverQt, weightEe);
				}
			}
		}
		else {
 		        if (MET>112 && Mt_mumu > 292 && passCuts[1]) passCuts[6] = kTRUE;
		        if (newProjMet>112 && Mt_mumu > 292 && passCuts[1]) passCuts[7] = kTRUE;

			h_pho2_pt->Fill(photonPt);
			h_pho2_mass->Fill(photonMassFromMumu);
			h_pho2_phi->Fill(photonPhi);
			h_pho2_eta->Fill(photonEta);

			Float_t weightMumu = photonWeightMumu(photonPt);

			h_phow2_pt->Fill(photonPt, weightMumu);
			h_phow2_mass->Fill(photonMassFromMumu, weightMumu);
			h_phow2_phi->Fill(photonPhi, weightMumu);
			h_phow2_eta->Fill(photonEta, weightMumu);

			for (int c=0; c<12; ++c) {
				if (passCuts[c]) {
					h_phow2_met[c]->Fill(MET, weightMumu);
					h_phow2_mt[c]->Fill(Mt_mumu, weightMumu);
					h_phow2_metOverQt[c]->Fill(metOverQt, weightMumu);
					h_phow2_nJets[c]->Fill(numJets, weightMumu);
					h_phow2_projMet[c]->Fill(newProjMet, weightMumu);
					h_phow2_projMetOverQt[c]->Fill(newProjMetOverQt, weightMumu);
				}
			}

		}

	} // if photonSample


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
		if (fabs(thisJet->P4(7).Eta()) <= 2.4) 
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

Float_t analyzer_higgs::photonWeightMumu(Float_t pt) {
	// the binning of the histograms are hardcoded for now

	Float_t X_MIN = 0;
	Float_t X_MAX = 400;
	Float_t N_BINS = 40;

	Float_t BIN_WIDTH = (X_MAX-X_MIN)/N_BINS;

	Float_t weight = 0;
	Int_t bin = 0;

	if  (pt<X_MAX) {
		weight = h1_photonWeightMumu->GetBinContent(Int_t(pt/BIN_WIDTH) + 1);
	}
	else {
		// we will put sonme extrapolation here...
		// possibly the low stat bins will also use some functions
		weight = h1_photonWeightMumu->GetBinContent(N_BINS);
	}
	return weight;
}

Float_t analyzer_higgs::photonWeightEe(Float_t pt) {
	// the binning of the histograms are hardcoded for now

	Float_t X_MIN = 0;
	Float_t X_MAX = 400;
	Float_t N_BINS = 40;

	Float_t BIN_WIDTH = (X_MAX-X_MIN)/N_BINS;

	Float_t weight = 0;
	Int_t binNum = Int_t(pt/BIN_WIDTH) + 1;

	//cout << "pt: " << pt << endl;
	//cout << "Bin num: " << Int_t(pt/BIN_WIDTH) + 1 << endl;

	if  (pt<X_MAX) {
		weight = h1_photonWeightEe->GetBinContent(binNum);
	}
	else {
		// we will put sonme extrapolation here...
		// possibly the low stat bins will also use some functions
		weight = h1_photonWeightEe->GetBinContent(N_BINS);
	}

	//	cout << "weight: " << weight << cout;
	return weight;
}

