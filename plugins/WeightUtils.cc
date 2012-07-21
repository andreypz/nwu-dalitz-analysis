#include "WeightUtils.h"

WeightUtils::WeightUtils(string sampleName, string dataPeriod, string selection, bool isRealData)
{
    _sampleName = sampleName;
    _dataPeriod = dataPeriod;
    _selection  = selection;
    _isRealData = isRealData;

    Initialize();
    rnGen   = new TRandom3();
    _inFile = new TFile("../data/Reweight.root", "OPEN");

    // Photon weights
    h1_eGammaPt     = (TH1D*)_inFile->Get("h1_eGammaPtWeightCombined");
    h1_muGammaPt    = (TH1D*)_inFile->Get("h1_muGammaPtWeightCombined");
    h1_eGammaPV     = (TH1D*)_inFile->Get("h1_eGammaPVWeightCombined");
    h1_muGammaPV    = (TH1D*)_inFile->Get("h1_muGammaPVWeightCombined");
    h1_eGammaMass   = (TH1D*)_inFile->Get("h1_diElectronMassCombined");
    h1_muGammaMass  = (TH1D*)_inFile->Get("h1_diMuonMassCombined");

    // PU weights
    h1_puReweight2011  = (TH1D*)_inFile->Get("h1_PU2011");
  
    // Recoil weights
    //h1_recoilLongMuon       = (TH1D*)_inFile->Get("h1_muonRecoilTransWeight");
    //h1_recoilTransMuon      = (TH1D*)_inFile->Get("h1_muonRecoilTransWeight");
    //h1_recoilLongElectron   = (TH1D*)_inFile->Get("h1_electronRecoilTransWeight");
    //h1_recoilTransElectron  = (TH1D*)_inFile->Get("h1_electronRecoilTransWeight");

    // higgs pt weights
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
  if (_dataPeriod == "2011" || _dataPeriod == "2011A" || _dataPeriod == "2011B")
    _puWeight = h1_puReweight2011->GetBinContent(nPU+1);
  else 
    _puWeight = 0;
  
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

  return _w;
}

float WeightUtils::GluGluHiggsWeight(float higgsPt, int higgsMass) 
{
    int iMass = (higgsMass - 250)/50;
    _glugluWeight = h1_Higgs[iMass]->GetBinContent(h1_Higgs[iMass]->FindBin(higgsPt));
    return _glugluWeight;
}


//This looks wrong - we dont reweight like that
//float WeightUtils::VBFHiggsWeight(float genMass, int higgsMass)
//{
//    _vbfWeight = float(higgsMass)/genMass;
//    return _vbfWeight;
//}

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


// Not used
//float WeightUtils::RecoilWeight(TLorentzVector recoilP4, TLorentzVector ZP4)/
//{
    //if (_selection == "muon") {
    //    _recoilLongWeight = h1_recoilLongMuon->GetBinContent()
    //}
//    return -1; //_recoilWeight;
//}
