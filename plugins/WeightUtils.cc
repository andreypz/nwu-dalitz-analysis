#include "WeightUtils.h"

WeightUtils::WeightUtils() {}

WeightUtils::~WeightUtils() {}

WeightUtils::WeightUtils(string sampleName, string dataPeriod, string selection, bool isRealData)
{
    _sampleName = sampleName;
    _dataPeriod = dataPeriod;
    _selection  = selection;
    _isRealData = isRealData;

    _puWeight = 1.;
    _zzWeight = 1.;
    _glugluWeight = 1.;
    _vbfWeight = 1.;
    _recoWeight = 1.;
    _triggerWeight = 1.;
    _gammaPtWeight = 1.;
    _gammaPVWeight = 1.;
    _gammaJetWeight = 1.;

    _inFile = new TFile("../data/Reweight.root", "OPEN");

    // Photon weights
    h1_eGammaPt     = (TH1D*)_inFile->Get("h1_eGammaPtWeightCombined");
    h1_muGammaPt    = (TH1D*)_inFile->Get("h1_muGammaPtWeightCombined");
    h1_eGammaPV     = (TH1D*)_inFile->Get("h1_eGammaPVWeightCombined");
    h1_muGammaPV    = (TH1D*)_inFile->Get("h1_muGammaPVWeightCombined");
    h1_eGammaMass   = (TH1D*)_inFile->Get("h1_diElectronMassCombined");
    h1_muGammaMass  = (TH1D*)_inFile->Get("h1_diMuonMassCombined");

    // PU weights
    h1_puReweight2011A  = (TH1D*)_inFile->Get("h1_PU2011A");
    h1_puReweight2011B  = (TH1D*)_inFile->Get("h1_PU2011B");

    // Higgs pt weights
    for (int i = 0; i < 8; ++i) {
        int higgsMass = 250+50*i;
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
    _recoWeight = 1.;
    _triggerWeight = 1.;
    _gammaPtWeight = 1.;
    _gammaPVWeight = 1.;
    _gammaJetWeight = 1.;
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

float WeightUtils::GetTotalWeight(int nPV, int nJets, TLorentzVector l1, TLorentzVector l2)
{
    Initialize();
    float weight = 1.;

    if (!_isRealData) {
        weight *= PUWeight(nPV);
        weight *= RecoWeight(l1, l2);
        if (_sampleName.compare(0,2,"ZZ") == 0) weight *= ZZWeight(l1, l2);
    } else {
        weight *= GammaWeight(nPV, nJets, l1);
    }
    return weight;
}

float WeightUtils::PUWeight(int nPU)
{
    if (nPU < 24 && _dataPeriod == "2011A"){
        _puWeight = h1_puReweight2011A->GetBinContent(nPU+1); 
    } else if (nPU < 35 && _dataPeriod == "2011B") {
        _puWeight = h1_puReweight2011B->GetBinContent(nPU+1); 
    } else _puWeight = 0;

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

float WeightUtils::GluGluHiggsWeight(float higgsPt, int higgsMass) 
{
    int iMass = (higgsMass - 250)/50;
    _glugluWeight = h1_Higgs[iMass]->GetBinContent(h1_Higgs[iMass]->FindBin(higgsPt));
    return _glugluWeight;
}

float WeightUtils::VBFHiggsWeight(float genMass, int higgsMass)
{
    _vbfWeight = float(higgsMass)/genMass;
    return _vbfWeight;
}

float WeightUtils::GammaWeight(int nPV, int nJets, TLorentzVector p1)
{
    if (_selection == "eGamma") {
      _gammaPtWeight = h1_eGammaPt->GetBinContent(h1_eGammaPt->FindBin(p1.Pt()));
      
      if (h1_eGammaPV->GetBinContent(1) != 1) {
	_gammaPVWeight = h1_eGammaPV->GetBinContent(nPV);
	if (_dataPeriod == "2011B") {
	  if (nJets == 0) _gammaJetWeight = 1.063;
	  if (nJets == 1) _gammaJetWeight = 1.0;
	} else if (_dataPeriod == "2011A") {
	  if (nJets == 0) _gammaJetWeight = 1.099;
	  if (nJets == 1) _gammaJetWeight = 0.99;
	}
      }
    }
    
    if (_selection == "muGamma") {
      _gammaPtWeight = h1_muGammaPt->GetBinContent(h1_muGammaPt->FindBin(p1.Pt()));
      
      if (h1_muGammaPV->GetBinContent(1) != 1) {
	_gammaPVWeight = h1_muGammaPV->GetBinContent(nPV);
	if (_dataPeriod == "2011B") {
	  if (nJets == 0) _gammaJetWeight = 1.149;
	  if (nJets == 1) _gammaJetWeight = 1.025;
	} else if (_dataPeriod == "2011A") {
	  if (nJets == 0) _gammaJetWeight = 1.072;
	  if (nJets == 1) _gammaJetWeight = 1.016;
	}
      }
    }
    return _gammaPtWeight*_gammaPVWeight*_gammaJetWeight;
}

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
  if (_selection == "eGamma") {
    photonMass = h1_eGammaMass->GetRandom();
  }
  if (_selection == "muGamma") {
    photonMass = h1_muGammaMass->GetRandom();
  }
  return photonMass;
}
