#include "WeightUtils.h"

WeightUtils::WeightUtils(string sampleName, string dataPeriod, string selection, bool isRealData)
{
    _sampleName = sampleName;
    _dataPeriod = dataPeriod;
    _selection  = selection;
    _isRealData = isRealData;

    Initialize();
    rnGen   = new TRandom3();
    _puFile = new TFile("../data/puReweight.root", "OPEN");

    TString MC("S10");
    //TString MC("RD1");
    _phFileID = new TFile("../data/Photon_ID_CSEV_SF_Jan22rereco_Full2012_"+MC+"_MC_V01.root", "OPEN");
    //_muFileID  = new TFile("../data/MuonEfficiencies_Run2012ReReco_53X.root", "OPEN");
    _muFileISO = new TFile("../data/MuonEfficiencies_ISO_Run_2012ReReco_53X.root", "OPEN");

    h2_photonIDSF   = (TH2D*)_phFileID->Get("PhotonIDSF_MediumWP_Jan22rereco_Full2012_"+MC+"_MC_V01");
    h2_photonCSEVSF = (TH2D*)_phFileID->Get("PhotonCSEVSF_MediumWP_Jan22rereco_Full2012_"+MC+"_MC_V01");
    
    //h2_photonIDSF->Print("all");
    //h2_photonCSEVSF->Print("all");
    //for (Int_t i = 1; i < 9; i++)
    //cout<<i<<" bin "<<h2_photonIDSF->GetXaxis()->GetBinLowEdge(i)<<endl;

    //  Bin edges are:
    // X(Pt): 0 - 10 - 15 - 20 - 30 - 50 - 1000
    // Y(Eta): 0 - 0.8 - 1.444 - 1.566 - 2.0 - 2.5 - 3 - 3.5

    //gr_muonIDSF_eta09 = (TGraph*)_muFileID->Get("DATA_over_MC_Loose_pt_abseta<0.9");
    //gr_muonIDSF_eta12 = (TGraph*)_muFileID->Get("DATA_over_MC_Loose_pt_abseta0.9-1.2");
    //gr_muonIDSF_eta21 = (TGraph*)_muFileID->Get("DATA_over_MC_Loose_pt_abseta1.2-2.1");
    //gr_muonIDSF_eta24 = (TGraph*)_muFileID->Get("DATA_over_MC_Loose_pt_abseta2.1-2.4");
    
    //Because the scale factors for muons are stupidly saved in TGraphs, we need to know 
    //the binning beforehand. And it is (for pt bins):
    //10-20-25-30-35-40-50-60-90-140-300
    // This works somehow:
    //_muonIDSF_eta09 = gr_muonIDSF_eta09->GetY();
    //_muonIDSF_eta12 = gr_muonIDSF_eta12->GetY();
    //_muonIDSF_eta21 = gr_muonIDSF_eta21->GetY();
    //_muonIDSF_eta24 = gr_muonIDSF_eta24->GetY();

    //gr_muonIDSF_eta12->Print("all");
    //for (Int_t i = 0; i < 10; i++)
    // cout<<i<<" bin "<<_muonIDSF_eta12[i]<<endl;
    
    /*
    gr_muonIDSF_eta09->Print("all");
    for (Int_t i = 0; i < 10; i++)
      cout<<i<<" bin "<<_muonIDSF_eta09[i]<<endl;

      gr_muonIDSF_eta12->Print("all");
      for (Int_t i = 0; i < 10; i++)
      cout<<i<<" bin "<<_muonIDSF_eta12[i]<<endl;
    
    gr_muonIDSF_eta21->Print("all");
    for (Int_t i = 0; i < 10; i++)
      cout<<i<<" bin "<<_muonIDSF_eta21[i]<<endl;
    
      gr_muonIDSF_eta24->Print("all");
      for (Int_t i = 0; i < 10; i++)
      cout<<i<<" bin "<<_muonIDSF_eta24[i]<<endl;
    */

    gr_muonISOSF_eta09 = (TGraph*)_muFileISO->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Loose_pt_abseta<0.9");
    gr_muonISOSF_eta12 = (TGraph*)_muFileISO->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Loose_pt_abseta0.9-1.2");
    gr_muonISOSF_eta21 = (TGraph*)_muFileISO->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Loose_pt_abseta1.2-2.1");
    gr_muonISOSF_eta24 = (TGraph*)_muFileISO->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Loose_pt_abseta2.1-2.4");

    _muonISOSF_eta09 = gr_muonISOSF_eta09->GetY();
    _muonISOSF_eta12 = gr_muonISOSF_eta12->GetY();
    _muonISOSF_eta21 = gr_muonISOSF_eta21->GetY();
    _muonISOSF_eta24 = gr_muonISOSF_eta24->GetY();

    //gr_muonISOSF_eta24->Print("all");
    //for (Int_t i = 0; i < 10; i++)
    // cout<<i<<" bin "<<_muonISOSF_eta24[i]<<endl;

    
    // Photon weights
    //h1_eGammaPt     = (TH1D*)_inFile->Get("h1_eGammaPtWeightCombined");
    //h1_muGammaPt    = (TH1D*)_inFile->Get("h1_muGammaPtWeightCombined");
    //h1_eGammaPV     = (TH1D*)_inFile->Get("h1_eGammaPVWeightCombined");
    //h1_muGammaPV    = (TH1D*)_inFile->Get("h1_muGammaPVWeightCombined");
    //h1_eGammaMass   = (TH1D*)_inFile->Get("h1_diElectronMassCombined");
    //h1_muGammaMass  = (TH1D*)_inFile->Get("h1_diMuonMassCombined");

    // PU weights
    h1_puReweight2011  = (TH1D*)_puFile->Get("h1_PU2011");
    h1_puReweight2012  = (TH1D*)_puFile->Get("h1_PU2012");
  
    // Recoil weights
    //h1_recoilLongMuon       = (TH1D*)_inFile->Get("h1_muonRecoilTransWeight");
    //h1_recoilTransMuon      = (TH1D*)_inFile->Get("h1_muonRecoilTransWeight");
    //h1_recoilLongElectron   = (TH1D*)_inFile->Get("h1_electronRecoilTransWeight");
    //h1_recoilTransElectron  = (TH1D*)_inFile->Get("h1_electronRecoilTransWeight");

    // higgs pt weights
    for (int i = 0; i < 10; ++i) {
      int higgsMass =0;
      if (i==0)
	higgsMass = 125;
      else
        higgsMass = 200+ 50*(i-1);
      
      TFile *higgsFile = new TFile(Form("../data/Kfactors_%i_AllScales.root", higgsMass), "OPEN");
      TH1D *h1_tmp = (TH1D*)higgsFile->GetDirectory("kfactors")->Get(Form("kfact_mh%i_ren%i_fac%i", higgsMass, higgsMass, higgsMass));
      
      h1_Higgs[i] = (TH1D*)h1_tmp->Clone();
      //higgsFile->Close();
      
    }
}

void WeightUtils::Initialize()
{
  _puWeight = 1.;
  _zzWeight = 1.;
  _glugluWeight = 1.;
  _vbfWeight = 1.;
  _recoWeight = 1.;
  _triggerWeight = 1.;
  _gammaPtWeight = 1.;
  _gammaPVWeight = 1.;
  _gammaJetWeight = 1.;
  _recoilLongWeight = 1.;
  _recoilTransWeight = 1.;
}

void WeightUtils::SetDataBit(bool isRealData)
{
    _isRealData = isRealData;
}

void WeightUtils::SetDataPeriod(string dataPeriod)
{
    _dataPeriod = dataPeriod;
}

void WeightUtils::SetSampleName(string sampleName)
{
    _sampleName = sampleName;
}

void WeightUtils::SetSelection(string selection)
{
    _selection = selection;
}

float WeightUtils::GetTotalWeight(float nPV, int nJets, TLorentzVector l1, TLorentzVector l2)
{
    Initialize();
    float weight = 1.;

    weight *= PUWeight(nPV);

    if (!_isRealData) {
      weight *= RecoWeight(l1, l2);
      //if (_sampleName.compare(0,2,"ZZ") == 0 && _dataPeriod!="2012") weight *= ZZWeight(l1, l2);
    } else {
      //weight *= GammaWeight(nPV, nJets, l1);
    }
    return weight;
}


float WeightUtils::PhotonSF(TLorentzVector ph)
{
  Float_t wID = 1;
  Float_t wCS = 1;

  Int_t myBinPt  = h2_photonIDSF->GetXaxis()->FindBin(ph.Pt());
  Int_t myBinEta = h2_photonIDSF->GetYaxis()->FindBin(fabs(ph.Eta()));
  //cout<<"xbin(pt) = "<<myBinPt<<"  ybin(eta)="<<myBinEta<<endl;
  wID =  h2_photonIDSF->GetBinContent(myBinPt, myBinEta);
  myBinPt  = h2_photonCSEVSF->GetXaxis()->FindBin(ph.Pt());
  myBinEta = h2_photonCSEVSF->GetYaxis()->FindBin(fabs(ph.Eta()));
  wCS =  h2_photonCSEVSF->GetBinContent(myBinPt, myBinEta);
  
  //cout<<"Photon pt = "<<ph.Pt()<<"  eta="<<ph.Eta()<<"   IDSF ="<<wID<<"  CSCF="<<wCS<<endl;
  if (!isfinite(wID*wCS))
    cout<<"Photon pt = "<<ph.Pt()<<"  eta="<<ph.Eta()<<"   IDSF ="<<wID<<"  CSCF="<<wCS<<endl;
  return wID*wCS;
}

float WeightUtils::MuonSF(TLorentzVector mu)
{
  Float_t wID  = 1;
  Float_t wISO = 1;

  //cout<<"mu pt="<<mu.Pt()<<"  eta="<<mu.Eta()<<"  wID="<<wID<<"  wISO="<<wISO<<endl;
  //Double_t * mySFID;
  Double_t * mySFISO;
  // ** We do not apply muID SF for dalitz analysis! **
  //10-20-25-30-35-40-50-60-90-140-300
  if (fabs(mu.Eta())<0.9){
    //mySFID  = _muonIDSF_eta09;
    mySFISO = _muonISOSF_eta09;
  } 
  else if (fabs(mu.Eta())<1.2){
    //mySFID  = _muonIDSF_eta12; 
    mySFISO = _muonISOSF_eta12;
  }
  else if (fabs(mu.Eta())<2.1){
    //mySFID  = _muonIDSF_eta21;
    mySFISO = _muonISOSF_eta21;
  }
  else if (fabs(mu.Eta())<2.4){
    //mySFID  = _muonIDSF_eta24; 
    mySFISO = _muonISOSF_eta24;
  }
  else return 1.0;

  //if (mu.Pt() > 10){
    //wID  = mySFID[0];
    //wISO = mySFISO[0]; only apply to pt>20 muons
  //}
  if (mu.Pt() > 20){
    //wID  = mySFID[1];
    wISO = mySFISO[1];
  }
  if (mu.Pt() > 25){
    //wID  = mySFID[2];
    wISO = mySFISO[2];
  }
  if (mu.Pt() > 30){
    //wID  = mySFID[3];
    wISO = mySFISO[3];
  }
  if (mu.Pt() > 35){
    //wID  = mySFID[4];
    wISO = mySFISO[4];
  }
  if (mu.Pt() > 40){
    //wID  = mySFID[5];
    wISO = mySFISO[5];
  }
  if (mu.Pt() > 50){
    //wID  = mySFID[6];
    wISO = mySFISO[6];
  }
  if (mu.Pt() > 60){
    //wID  = mySFID[7];
    wISO = mySFISO[7];
  }
  if (mu.Pt() > 90){
    //wID  = mySFID[8];
    wISO = mySFISO[8];
  }
  if (mu.Pt() > 140){
    //wID  = mySFID[9];
    wISO = mySFISO[9];
  }

  //cout<<"mu pt="<<mu.Pt()<<"  eta="<<mu.Eta()<<"  wID="<<wID<<"  wISO="<<wISO<<endl;
  if (!isfinite(wID*wISO))
    cout<<"mu pt="<<mu.Pt()<<"  eta="<<mu.Eta()<<"  wID="<<wID<<"  wISO="<<wISO<<endl;
  return wID*wISO;
}

  float WeightUtils::PUWeight(float nPUtrue)
{
  if (_isRealData) return 1;

  if (nPUtrue < 60 && (_dataPeriod == "2011" || _dataPeriod == "2011A" || _dataPeriod == "2011B"))
    {
      Int_t myBin = h1_puReweight2011->FindBin(nPUtrue);
      _puWeight   = h1_puReweight2011->GetBinContent(myBin);
      //cout<<"period = "<<_dataPeriod<<"  Pu weight "<<_puWeight<<"  bin = "<<myBin<<endl;
    }
  else if  (nPUtrue < 60 &&  _dataPeriod == "2012") 
    {
      Int_t myBin = h1_puReweight2012->FindBin(nPUtrue);
      _puWeight   = h1_puReweight2012->GetBinContent(myBin);
      //_puWeight = 1;

      //cout<<"period = "<<_dataPeriod<<"  Pu weight "<<_puWeight<<"  bin = "<<myBin<<endl;
    } 
  else
    _puWeight = 1;
  //cout  <<" PU  weight = "<<_puWeight<<endl;
  return _puWeight;
}

float WeightUtils::RecoWeight(TLorentzVector l1, TLorentzVector l2)
{
    if (_selection == "muon") {
        _triggerWeight = GetMuTriggerEff(l1)*GetMuTriggerEff(l2);
        _recoWeight    = 1.;
    }
    if (_selection == "electron") {
        _triggerWeight = 0.995;
        _recoWeight    = GetElectronEff(l1)*GetElectronEff(l2);
    }
    return _triggerWeight*_recoWeight;
}

float WeightUtils::ZZWeight(TLorentzVector l1, TLorentzVector l2)
{
    float zPt = (l1 + l2).Pt();
    _zzWeight = 0.12 + (1.108 + 0.002429*zPt - (1.655e-6)*pow(zPt, 2));  
    return _zzWeight;
}


float WeightUtils::HiggsMassLineShapeWeight(float mass, float m_gen) 
{
  double _w;
  /*
  double mh, gh, mt, m, _w;
  double decay_width;
  int BWflag;

  if(m_gen == 130){    decay_width =   0.00487; 
  }else if(m_gen == 140){    decay_width =   0.00812; 
  }else if(m_gen == 150){    decay_width =   0.01730; 
  }else if(m_gen == 160){    decay_width =   0.08290; 
  }else if(m_gen == 170){    decay_width =   0.38000; 
  }else if(m_gen == 180){    decay_width =   0.63000; 
  }else if(m_gen == 190){    decay_width =   1.04000; 
  }else if(m_gen == 200){    decay_width =   1.43000; 
  }else if(m_gen == 250){    decay_width =   4.04000; 
  }else if(m_gen == 300){    decay_width =   8.43000; 
  }else if(m_gen == 350){    decay_width =  15.20000; 
  }else if(m_gen == 400){    decay_width =  29.20000; 
  }else if(m_gen == 450){    decay_width =  46.95000; 
  }else if(m_gen == 500){    decay_width =  68.00000;
  }else if(m_gen == 550){    decay_width =  93.15000;
  }else if(m_gen == 600){    decay_width = 123.00000;
  }

  mh = m_gen;
  gh = decay_width; 
  mt = 172.5;
  m  = mass;
  BWflag = 0;
  _w  =  -1;  // invalid initial value                                                                                        


  pwhg_cphto_reweight_(&mh, &gh, &mt, &BWflag, &m, &_w);
*/

/* Usage of  pwhg_cphto_reweight_() function:
c     INPUT
c     mh : Higgs boson mass (used in the POWHEG BOX generation)
c     gh : Higgs boson width (used in the POWHEG BOX generation)
c     mt : top quark mass
c     BWflag : 0    if the sample to reweight was produced with fixed Higgs width
c              1    if the sample to reweight was produced with running Higgs 
c                   width (this is the default in the POWHEG BOX)
c     m : virtuality of the produced Higgs boson resonance
c     OUTPUT
c     w : the reweighting factor 
*/

//  cout <<"mass= "<<mh<< "  w = " << _w << endl;
  _w=-1;
  return _w;
}


float WeightUtils::GluGluHiggsWeight(float higgsPt, int higgsMass) 
{
  int iMass = 0;
  if (higgsMass==125)
    iMass = 0;
  else
    iMass = 1+(higgsMass - 200)/50;
    _glugluWeight = h1_Higgs[iMass]->GetBinContent(h1_Higgs[iMass]->FindBin(higgsPt));
    return _glugluWeight;
}


//This looks wrong - we dont reweight like that
//float WeightUtils::VBFHiggsWeight(float genMass, int higgsMass)
//{
//    _vbfWeight = float(higgsMass)/genMass;
//    return _vbfWeight;
//}

/*
float WeightUtils::GammaWeight(int nPV, int nJets, TLorentzVector p1)
{
    if (_selection == "eGamma") {
        _gammaPtWeight = h1_eGammaPt->GetBinContent(h1_eGammaPt->FindBin(p1.Pt()));

        if (h1_eGammaPV->GetBinContent(1) != 1) {
	  _gammaPVWeight = h1_eGammaPV->GetBinContent(nPV);
	  if (nJets == 0) _gammaJetWeight = 1.0;
	  if (nJets == 1) _gammaJetWeight = 1.0;
	}
    }
    
    if (_selection == "muGamma") {
      _gammaPtWeight = h1_muGammaPt->GetBinContent(h1_muGammaPt->FindBin(p1.Pt()));
      
      if (h1_muGammaPV->GetBinContent(1) != 1) {
	_gammaPVWeight = h1_muGammaPV->GetBinContent(nPV);
	if (nJets == 0) _gammaJetWeight = 1.0;
	if (nJets == 1) _gammaJetWeight = 1.0;
      }
    
    }
    return _gammaPtWeight*_gammaPVWeight*_gammaJetWeight;
}
*/

float WeightUtils::GetMuTriggerEff(TLorentzVector l1) const
{
    //float doubleMu7Scale[4][4] = {
    //    {0.98, 0.974, 0.978, 0.978},
    //    {0.947, 0.947, 0.948, 0.947},
    //    {0.959, 0.958, 0.958, 0.956},
    //    {0.926, 0.925, 0.913, 0.921}
    //};

    float mu13_mu8Scale[4][4] = {
        {0.9753, 0.9751, 0.9746, 0.9745},
        {0.9627, 0.9568, 0.9563, 0.9575},
        {0.9490, 0.9482, 0.9503, 0.9462},
        {0.8594, 0.8614, 0.867, 0.8529}
    };

    int ptBin = 0;
    int etaBin = 0;
    float binningPt[]  = {20., 30., 40., 50., 99999.};
    float binningEta[] = {0., 0.8, 1.2, 2.1, 2.4};
    float weight = 1.;

    for (int i = 0; i < 4; ++i) {
        if (fabs(l1.Eta()) > binningEta[i] && fabs(l1.Eta()) <= binningEta[i+1]) {
            etaBin = i;
            break;
        }
    }
    for (int i = 0; i < 4; ++i) {
        if (l1.Pt() > binningPt[i] && l1.Pt() <= binningPt[i+1]) {
            ptBin = i;
            break;
        }
    }
    weight = mu13_mu8Scale[etaBin][ptBin];
    return weight;
}

float WeightUtils::GetElectronEff(TLorentzVector l1) const
{
    float eleScale[2][5] = {
        {1.0, 0.938, 0.971, 0.980, 0.987},
        {1.0, 0.967, 0.97, 0.989, 0.989}
    };

    int ptBin = 0;
    int etaBin = 0;
    float binning[] = {10., 20., 30., 40., 50., 9999.};
    float weight = 1.;

    if (fabs(l1.Eta()) < 1.479) {
        etaBin = 0;
    } else {
        etaBin = 1;
    }

    for (int i = 0; i < 5; ++i) {
        if (l1.Pt() > binning[i] && l1.Pt() <= binning [i+1]) {
            ptBin = i;
            break;
        }
    }
    weight = eleScale[etaBin][ptBin];
    return weight;
}

float WeightUtils::GetPhotonMass() const
{
    float photonMass = 91.2;
    //if (_selection == "eGamma") {
      //photonMass = h1_eGammaMass->GetRandom();
    //}
    //if (_selection == "muGamma") {
      //photonMass = h1_muGammaMass->GetRandom();
    //}
    return photonMass;
}


// Not used
//float WeightUtils::RecoilWeight(TLorentzVector recoilP4, TLorentzVector ZP4)/
//{
    //if (_selection == "muon") {
    //    _recoilLongWeight = h1_recoilLongMuon->GetBinContent()
    //}
//    return -1; //_recoilWeight;
//}
