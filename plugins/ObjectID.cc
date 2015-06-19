#include "ObjectID.h"

bool SCESortCondition(const TCEGamma& p1, const TCEGamma& p2) {return (p1.SCEnergy() > p2.SCEnergy());}
bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());}

float _EAPho[7][3] = {
  {0.012,  0.030,   0.148}, //         eta < 1.0
  {0.010,  0.057,   0.130}, // 1.0   < eta < 1.479
  {0.014,  0.039,   0.112}, // 1.479 < eta < 2.0
  {0.012,  0.015,   0.216}, // 2.0   < eta < 2.2
  {0.016,  0.024,   0.262}, // 2.2   < eta < 2.3
  {0.020,  0.039,   0.260}, // 2.3   < eta < 2.4
  {0.012,  0.072,   0.266}  // 2.4   < eta
};
//                EB   EE
float _phoMvaCuts[2] = {0.08, 0.08};
float _DAleMvaCut    = 0.12;

ObjectID::ObjectID():
  _isRealData(0),
  _runNumber(0),
  _eventNumber(0),
  _rhoFactor(0)
{
  //Cuts
  //electrons are two types: Barrel/Endcap
  //Loose
  elIdAndIsoCutsLoose.ptErrorOverPt[0] = 99999;
  elIdAndIsoCutsLoose.dPhiIn[0]        = 0.1;
  elIdAndIsoCutsLoose.dEtaIn[0]        = 0.01;
  elIdAndIsoCutsLoose.sigmaIetaIeta[0] = 0.02;
  elIdAndIsoCutsLoose.HadOverEm[0]     = 0.3;
  elIdAndIsoCutsLoose.fabsEPDiff[0]    = 99999;
  elIdAndIsoCutsLoose.dxy[0]           = 0.5;
  elIdAndIsoCutsLoose.dz[0]            = 0.5;
  elIdAndIsoCutsLoose.pfIso04[0]       = 0.3;

  elIdAndIsoCutsLoose.ptErrorOverPt[1] = 99999;
  elIdAndIsoCutsLoose.dPhiIn[1]        = 0.1;
  elIdAndIsoCutsLoose.dEtaIn[1]        = 0.03;
  elIdAndIsoCutsLoose.sigmaIetaIeta[1] = 0.1;
  elIdAndIsoCutsLoose.HadOverEm[1]     = 0.3;
  elIdAndIsoCutsLoose.fabsEPDiff[1]    = 99999;
  elIdAndIsoCutsLoose.dxy[1]           = 0.5;
  elIdAndIsoCutsLoose.dz[1]            = 0.5;
  elIdAndIsoCutsLoose.pfIso04[1]       = 0.3;

  //Tight
  elIdAndIsoCutsTight.ptErrorOverPt[0] = 9999.;
  elIdAndIsoCutsTight.dPhiIn[0]        = 0.06;
  elIdAndIsoCutsTight.dEtaIn[0]        = 0.004;
  elIdAndIsoCutsTight.sigmaIetaIeta[0] = 0.01;
  elIdAndIsoCutsTight.HadOverEm[0]     = 0.12;
  elIdAndIsoCutsTight.fabsEPDiff[0]    = 0.05;
  elIdAndIsoCutsTight.dxy[0] = 0.02;
  elIdAndIsoCutsTight.dz[0]  = 0.1;
  elIdAndIsoCutsTight.pfIso04[0] = 0.15;

  elIdAndIsoCutsTight.ptErrorOverPt[1] = 9999.;
  elIdAndIsoCutsTight.dPhiIn[1]        = 0.03;
  elIdAndIsoCutsTight.dEtaIn[1]        = 0.007;
  elIdAndIsoCutsTight.sigmaIetaIeta[1] = 0.03;
  elIdAndIsoCutsTight.HadOverEm[1]     = 0.10;
  elIdAndIsoCutsTight.fabsEPDiff[1]    = 0.05;
  elIdAndIsoCutsTight.dxy[1] = 0.02;
  elIdAndIsoCutsTight.dz[1]  = 0.1;
  elIdAndIsoCutsTight.pfIso04[1] = 0.15;

  //Soft muon id (specific for dalitz)
  muIdAndIsoCutsSoft.TrackLayersWithMeasurement = 5;
  muIdAndIsoCutsSoft.PixelLayersWithMeasurement = 1;
  muIdAndIsoCutsSoft.NormalizedChi2_tracker = 3;
  muIdAndIsoCutsSoft.dxy             = 0.2; //3
  muIdAndIsoCutsSoft.dz              = 0.5; //30
  muIdAndIsoCutsSoft.pfIso04         = 0.4;

  //muon id cuts
  muIdAndIsoCutsTight.ptErrorOverPt = 9999.;
  muIdAndIsoCutsTight.TrackLayersWithMeasurement = 5;
  muIdAndIsoCutsTight.NumberOfValidMuonHits      = 0;
  muIdAndIsoCutsTight.NumberOfValidTrackerHits   = 10;
  muIdAndIsoCutsTight.NumberOfValidPixelHits     = 0;
  muIdAndIsoCutsTight.NumberOfMatches            = 1;
  muIdAndIsoCutsTight.NumberOfMatchedStations    = 1;
  muIdAndIsoCutsTight.NormalizedChi2 = 10;
  muIdAndIsoCutsTight.dxy = 0.2;
  muIdAndIsoCutsTight.dz  = 0.5;
  muIdAndIsoCutsTight.pfIso04 = 0.4;

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

  jetIdCutsVBF.betaStarC[0] = 0.2;
  jetIdCutsVBF.dR2Mean[0] = 0.06;
  jetIdCutsVBF.betaStarC[1] = 0.3;
  jetIdCutsVBF.dR2Mean[1] = 0.05;
  jetIdCutsVBF.dR2Mean[2] = 0.05;
  jetIdCutsVBF.dR2Mean[3] = 0.055;

}

ObjectID::~ObjectID(){}

void ObjectID::SetEventInfo(Bool_t is, UInt_t run, ULong64_t evt, Float_t rho) {
  _isRealData  = is;
  _runNumber   = run;
  _eventNumber = evt;
  _rhoFactor   = rho;

}

void ObjectID::CalculatePhotonIso(const TCPhoton& ph, float& chIsoCor, float& nhIsoCor, float& phIsoCor){
  float chEA,nhEA,phEA,tmpEta;

  //cout<<"ObjID.  rhoFactor = "<<_rhoFactor<<endl;
  tmpEta = ph.SCEta();

  if (fabs(tmpEta) < 1.0){
    chEA = _EAPho[0][0];
    nhEA = _EAPho[0][1];
    phEA = _EAPho[0][2];
  }else if (fabs(tmpEta) < 1.479){
    chEA = _EAPho[1][0];
    nhEA = _EAPho[1][1];
    phEA = _EAPho[1][2];
  }else if (fabs(tmpEta) < 2.0){
    chEA = _EAPho[2][0];
    nhEA = _EAPho[2][1];
    phEA = _EAPho[2][2];
  }else if (fabs(tmpEta) < 2.2){
    chEA = _EAPho[3][0];
    nhEA = _EAPho[3][1];
    phEA = _EAPho[3][2];
  }else if (fabs(tmpEta) < 2.3){
    chEA = _EAPho[4][0];
    nhEA = _EAPho[4][1];
    phEA = _EAPho[4][2];
  }else if (fabs(tmpEta) < 2.4){
    chEA = _EAPho[5][0];
    nhEA = _EAPho[5][1];
    phEA = _EAPho[5][2];
  }else{
    chEA = _EAPho[6][0];
    nhEA = _EAPho[6][1];
    phEA = _EAPho[6][2];
  }

  chIsoCor = ph.PfIsoCharged() - _rhoFactor*chEA;
  nhIsoCor = ph.PfIsoNeutral() - _rhoFactor*nhEA;
  phIsoCor = ph.PfIsoPhoton()  - _rhoFactor*phEA;
  //cout<<chIsoCor<<"  "<<nhIsoCor<<"  "<<phIsoCor<<endl;

}

bool ObjectID::PassPhotonIdAndIso(const TCPhoton& ph, TString n, float& mvaScore)
{
  mvaScore = -99;
  phIdAndIsoCuts cuts;

  if (n=="CutBased-MediumWP")
    cuts = phIdAndIsoCutsHZG;
  else if (n=="CutBased-MediumWP-El"){
    if (!HggPreselection(ph)) return false;
    cuts = phIdAndIsoCutsHZG;
  }
  else if (n=="CutBased-TightWP")
    cuts = phIdAndIsoCutsTight;
  else if (n=="MVA")
    return PassPhotonMVA(ph, mvaScore);
  else  {
    cout<<"PhotonID Warning: no such cut or method"<<n<<endl;
    return false;
  }

  bool pass = false;
  //Float_t phoISO = 0;
  float tmpEta = ph.SCEta();

  float chIsoCor=0, nhIsoCor=0, phIsoCor=0;
  CalculatePhotonIso(ph, chIsoCor,nhIsoCor,phIsoCor);

  if(
     (fabs(tmpEta)  < 1.4442
      && ph.ConversionVeto()       == cuts.PassedEleSafeVeto[0]
      && ph.HadOverEm()             < cuts.HadOverEm[0]
      && ph.SigmaIEtaIEta()         < cuts.sigmaIetaIeta[0]
      && max((double)chIsoCor,0.)   < cuts.chIso03[0]
      && max((double)nhIsoCor,0.)   < cuts.nhIso03[0] + 0.040*ph.Pt()
      && max((double)phIsoCor,0.)   < cuts.phIso03[0] + 0.005*ph.Pt()
      ) ||
     (fabs(tmpEta)  > 1.566
      && ph.ConversionVeto()       == cuts.PassedEleSafeVeto[1]
      && ph.HadOverEm()             < cuts.HadOverEm[1]
      && ph.SigmaIEtaIEta()         < cuts.sigmaIetaIeta[1]
      && max((double)chIsoCor,0.)   < cuts.chIso03[1]
      && max((double)nhIsoCor,0.)   < cuts.nhIso03[1] + 0.040*ph.Pt()
      && max((double)phIsoCor,0.)   < cuts.phIso03[1] + 0.005*ph.Pt()
      )
     ) pass = true;

  return pass;
}

/*
  void ObjectID::PhotonR9Corrector(const TCPhoton& ph){
  //old R9 correction
  float R9Cor;
  R9Cor = ph.R9();

  if (fabs(ph.SCEta()) < 1.479)
  R9Cor = ph.R9()*1.0045 + 0.0010;
  else
  R9Cor = ph.R9()*1.0086 - 0.0007;

  //ph.SetR9(R9Cor);
  }
*/

float ObjectID::CalculateMuonIso(const TCMuon& lep)
{
  float muISO = (lep.PfIsoCharged() + TMath::Max(0.0, lep.PfIsoPhoton() + lep.PfIsoNeutral() - 0.5*lep.PfIsoPU()))/lep.Pt();
  return muISO;
}

//if (ObjID->CalculateMuonIso(muons[1]) > 1.4)


bool ObjectID::PassMuonIso(const TCMuon& lep)
{
  float iso = ObjectID::CalculateMuonIso(lep);
  if (iso<0.12)
    return true;
  else
    return false;
}

bool ObjectID::PassMuonId(const TCMuon& lep, TVector3 *pv, TString n)
{
  muIdAndIsoCuts cuts;
  if (n=="Soft")
    cuts = muIdAndIsoCutsSoft;
  else if (n=="Tight")
    cuts = muIdAndIsoCutsTight;
  else {
    cout<<"MuonID Warning: no such cut "<<n<<endl;
    return false;
  }

  bool pass = false;


  if(n=="Soft"
     && lep.IsPF()
     && ((lep.IsTRK() && lep.NumberOfMatches() > 0) || lep.IsGLB())
     && fabs(lep.Dxy(pv))     < cuts.dxy
     && fabs(lep.Dz(pv))      < cuts.dz
     )
    pass = true;

  else if(n == "Tight"
          && lep.IsPF() && lep.IsGLB()
          && lep.NormalizedChi2()             < cuts.NormalizedChi2
          && lep.NumberOfValidMuonHits()      > cuts.NumberOfValidMuonHits
          && lep.NumberOfMatchedStations()    > cuts.NumberOfMatchedStations
          && lep.NumberOfValidPixelHits()     > cuts.NumberOfValidPixelHits
          && lep.TrackLayersWithMeasurement() > cuts.TrackLayersWithMeasurement
          && fabs(lep.Dxy(pv)) < cuts.dxy
          && fabs(lep.Dz(pv))  < cuts.dz
          )
    pass = true;
  return pass;
}


bool ObjectID::PassMuonIdAndIso(const TCMuon& lep, TVector3 *pv, TString n)
{
  if (PassMuonId(lep, pv, n) && PassMuonIso(lep))
    return true;
  else
    return false;
}

float ObjectID::CalculateElectronIso(const TCElectron& lep)
{
  float eleISO = 0;
  //if(period=="2012")
  //eleISO = (lep.IsoMap("pfChIso_R04") + TMath::Max(0.0, (Double_t)(lep.IsoMap("pfPhoIso_R04") + lep.IsoMap("pfNeuIso_R04") - rho25Factor*lep.IsoMap("EffArea_R04"))))/lep.Pt();
  //else if(period=="2011")
  //eleISO = (lep.IsoMap("SumPt_R04") + TMath::Max(0.0, (Double_t)(lep.IsoMap("EmIso_R04") + lep.IsoMap("HadIso_R04") - rhoFactor*TMath::Pi()*0.09)))/lep.Pt();
  //else cout<<"Ele iso - Period is not defined!"<<period<<endl;

  return eleISO;
}



bool ObjectID::PassDalitzEleID(const TCElectron& el, TVector3 *pv, TString n, float& mvaScore)
{
  mvaScore = -99;
  if (n=="MVA")
    {
      //bool mvaPass = false;

      //if (!HggPreselection(el)) return false;

      vector<TCElectron::Track> trk = el.GetTracks();
      vector<TCEGamma> sc = el.BaseSC();

      //if (_eventNumber == 131810 || _eventNumber == 131886) {
      //cout<<_eventNumber<<" DBG DAlitz MVA ID "<<trk.size() <<"  "<<sc.size()<<endl;
      //sc[0].Dump();
      //sc[1].Dump();
      //sc[2].Dump();
      //}

      if (trk.size() < 2 || sc.size()<2)
        return false;

      sort(sc.begin(),  sc.end(),  SCESortCondition);
      sort(trk.begin(), trk.end(), P4SortCondition);


      if (fabs(el.Dxy(pv)) > 0.02 || fabs(el.Dz(pv)) > 0.1)
        return false;

      // classification variables
      static Float_t rho2012_, r9_, HoverE_, EoverPt_;
      static Float_t eleBCS25_1_, eleBCS25_2_, eleGSFPt2over1_, eleGSFdR_;// eleGSFPt1_, eleGSFPt2_;
      static Float_t SCEtaWidth_, SCPhiWidth_, dEtaAtVtx_, dPhiAtVtx_;
      static Float_t EcalEnPin_, recoSCEta_, eleBCSieie2_;
      static Float_t eleBCS15_2_, eleBCSieie1_;
      // MVA classifiers for 0=ECAL barrel and 1=ECAL endcaps
      static TMVA::Reader* tmvaDalitz = NULL;

      if (!tmvaDalitz) {
        tmvaDalitz = new TMVA::Reader("!Color:Silent");

        tmvaDalitz->AddVariable("rho2012",    &rho2012_);
        tmvaDalitz->AddVariable("eleBCS25_1", &eleBCS25_1_);
        tmvaDalitz->AddVariable("eleBCS25_2", &eleBCS25_2_);
        tmvaDalitz->AddVariable("eleGSFPt2over1",  &eleGSFPt2over1_);
        //tmvaDalitz->AddVariable("eleGSFPt1",  &eleGSFPt1_);
        //tmvaDalitz->AddVariable("eleGSFPt2",  &eleGSFPt2_);
        tmvaDalitz->AddVariable("eleGSFdR",   &eleGSFdR_);
        tmvaDalitz->AddVariable("r9", &r9_);
        tmvaDalitz->AddVariable("SCEtaWidth", &SCEtaWidth_);
        tmvaDalitz->AddVariable("SCPhiWidth", &SCPhiWidth_);
        tmvaDalitz->AddVariable("HoverE",     &HoverE_);
        tmvaDalitz->AddVariable("EoverPt",    &EoverPt_);
        tmvaDalitz->AddVariable("dEtaAtVtx",  &dEtaAtVtx_);
        tmvaDalitz->AddVariable("dPhiAtVtx",  &dPhiAtVtx_);
        tmvaDalitz->AddVariable("EcalEnPin",  &EcalEnPin_);
        tmvaDalitz->AddVariable("recoSCEta",  &recoSCEta_);
        tmvaDalitz->AddVariable("eleBCSieie2",&eleBCSieie2_);

        // Spectators
        tmvaDalitz->AddSpectator("eleBCS15_2", &eleBCS15_2_);
        tmvaDalitz->AddSpectator("eleBCSieie1",&eleBCSieie1_);

        tmvaDalitz->BookMVA("BDT", "../data/ElectronDalitzMVAID.xml");
      }

      r9_ = el.R9();
      rho2012_    = _rhoFactor;
      eleBCS25_1_ = sc[0].E2x5()/sc[0].E5x5();
      eleBCS25_2_ = sc[1].E2x5()/sc[1].E5x5();
      eleGSFPt2over1_  = trk[1].Pt()/trk[0].Pt();
      //eleGSFPt1_  = el.GetTracks()[0].Pt();
      //eleGSFPt2_  = el.GetTracks()[1].Pt();
      eleGSFdR_   = trk[0].DeltaR(trk[1]);
      SCEtaWidth_ = el.SCEtaWidth();
      SCPhiWidth_ = el.SCPhiWidth();
      HoverE_     = el.HadOverEm();
      EoverPt_    = el.SCRawEnergy()/(trk[0]+trk[1]).Pt();
      //EoverPt_    = el.SCEnergy()/(trk[0]+trk[1]).Pt();
      dEtaAtVtx_  = el.SCDeltaEta();
      dPhiAtVtx_  = el.SCDeltaPhi();
      EcalEnPin_  = el.InverseEnergyMomentumDiff(); //| 1/E - 1/Pin |
      recoSCEta_  = el.SCEta();
      eleBCSieie1_= sc[0].SigmaIEtaIEta();
      eleBCSieie2_= sc[1].SigmaIEtaIEta();
      eleBCS15_2_ = sc[1].E1x5()/sc[1].E5x5();

      mvaScore = tmvaDalitz->EvaluateMVA("BDT");

      //cout<<"\t *** DBG MVA ID score: "<<mvaScore<<"  gsf1 pt="<<trk[0].Pt()<<endl;
      if (mvaScore > _DAleMvaCut)
        return true;
      else
        return false;

    }
  else {
    cout<<"EleID Warning: no such cut "<<n<<endl;
    return false;
  }


}

bool ObjectID::PassElectronMVAPreSel(const TCElectron& el)
{
  bool pass = false;
  if (
      el.IdMap("gsf_numberOfLostHits") == 0
      && (el.IdMap("dr03TkSumPt")) / el.Pt() < 0.2
      && (el.IdMap("dr03EcalRecHitSumEt")) /el.Pt() < 0.2
      && (el.IdMap("dr03HcalTowerSumEt")) / el.Pt() < 0.2
      ) {
    if (fabs(el.Eta()) < 1.479) {
      if (el.SigmaIEtaIEta()< 0.014 && el.IdMap("hadronicOverEm") < 0.15)
        pass = true;
    } else { //endcap
      if (el.SigmaIEtaIEta()< 0.035 && el.IdMap("hadronicOverEm") < 0.10)
        pass = true;
    }
  }
  return pass;
}


bool ObjectID::PassElectronIdAndIsoMVA(const TCElectron& lep)
{
  bool pass = false;

  Float_t eleISO = CalculateElectronIso(lep);
  if (!PassElectronMVAPreSel(lep)) return false;

  if(
     eleISO < 0.15
     && lep.ConversionMissHits()==0
     && lep.PassConversionVeto() == true
     &&
     (
      (lep.Pt()>10 && lep.Pt()<20 &&
       (
        (fabs(lep.Eta()) < 0.8
         && lep.MvaID() > 0.00)
        ||
        (fabs(lep.Eta()) >  0.8 && fabs(lep.Eta()) < 1.479
         && lep.MvaID() > 0.10)
        ||
        (fabs(lep.Eta()) >  1.479 && fabs(lep.Eta()) < 2.5
         && lep.MvaID() > 0.62)
        )
       )
      ||
      (lep.Pt()>20 &&
       (
        (fabs(lep.Eta()) < 0.8
         && lep.MvaID() > 0.94)
        ||
        (fabs(lep.Eta()) >  0.8 && fabs(lep.Eta()) < 1.479
         && lep.MvaID() > 0.85)
        ||
        (fabs(lep.Eta()) >  1.479 && fabs(lep.Eta()) < 2.5
         && lep.MvaID() > 0.92)
        )
       )
      )
     )
    pass = true;

  return pass;
}


bool ObjectID::PassElectronIdAndIso(const TCElectron& lep, TVector3 *pv, TString n)
{
  elIdAndIsoCuts cuts;
  if (n=="Loose")
    cuts = elIdAndIsoCutsLoose;
  else if (n=="Tight")
    cuts = elIdAndIsoCutsTight;
  else {
    cout<<"ElectronID Warning: no such cut "<<n<<endl;
    return false;
  }

  bool pass = false;

  //Float_t eleISO = CalculateElectronIso(lep);

  if(((fabs(lep.Eta()) < 1.4442
       && lep.PtError()/lep.Pt() < cuts.ptErrorOverPt[0]
       && lep.SigmaIEtaIEta()    < cuts.sigmaIetaIeta[0]
       && fabs(lep.SCDeltaPhi()) < cuts.dPhiIn[0]
       && fabs(lep.SCDeltaEta()) < cuts.dEtaIn[0]
       && lep.HadOverEm()   < cuts.HadOverEm[0]
       && fabs(lep.Dxy(pv)) < cuts.dxy[0]
       && fabs(lep.Dz(pv))  < cuts.dz[0]
       //&& eleISO < cuts.pfIso04[0]
       )||
      (fabs(lep.Eta()) > 1.566
       && lep.PtError()/lep.Pt() < cuts.ptErrorOverPt[1]
       && lep.SigmaIEtaIEta()    < cuts.sigmaIetaIeta[1]
       && fabs(lep.SCDeltaPhi()) < cuts.dPhiIn[1]
       && fabs(lep.SCDeltaEta()) < cuts.dEtaIn[1]
       && lep.HadOverEm()   < cuts.HadOverEm[1]
       && fabs(lep.Dxy(pv)) < cuts.dxy[1]
       && fabs(lep.Dz(pv))  < cuts.dz[1]
       //&& eleISO < cuts.pfIso04[1]
       ))
     //&& lep.ConversionVeto()
     ) pass = true;

  //pass = true;

  //cout<<"ele pass? "<<pass<<endl;
  return pass;
}


TCGenParticle * ObjectID::GetPrimaryAncestor(TCGenParticle *p)
{
  TCGenParticle *a = p;
  while (a->Mother())
    a = a->Mother();
  return a;
}


bool ObjectID::isFSR(TCGenParticle *p)
{
  TCGenParticle *a = p;
  bool fsr = false;
  while (a->Mother())
    {
      if (a->Mother()->GetPDGId() == 22)
        return true;
      else
        a = a->Mother();
    }
  return fsr;
}


void ObjectID::DiscoverGeneology(TCGenParticle *p)
{
  cout<<"---->> event = "<<_eventNumber<<"   "<<p->GetPDGId()<<"  st = "<<p->GetStatus()
      <<"\n primary ancestor = "<<GetPrimaryAncestor(p)->GetPDGId()
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


bool ObjectID::HggPreselection(const TCPhoton& ph)
{
  /* Preselection of photon candidates for the identification with MVA.
   *
   * Returns true if candidate passes the preselection or false otherwise.
   */

  /*if (_eventNumber == 154245 || _eventNumber == 154548){
    cout<<"\t   Hgg preselection DBG for evt: "<<_eventNumber<<endl;

    cout<<"R9 = "<<ph.R9()<<"  h/e = "<<ph.HadOverEm()<<"   sieie ="<<ph.SigmaIEtaIEta()<<"  convVeto: "<<ph.ConversionVeto()<<endl;
    cout<<"Hcal iso:  "<< ph.IdMap("HadIso_R03") - 0.005 * ph.Et()<<"  ET="<<ph.Et()<<"  pt="<<ph.Pt()<<endl;
    cout<<"trk  iso:  "<< ph.IdMap("TrkIso_R03") - 0.002 * ph.Et()<<"  ET="<<ph.Et()<<"  pt="<<ph.Pt()<<endl;
    cout<<"CiC: "<<ph.CiCPF4chgpfIso02()[0]<<endl;
    }*/

  if (ph.Et() < 15 || ph.Et() > 200) return false;
  //if (ph.ConversionVeto() == 0) return false;

  if (fabs(ph.SCEta()) < 1.479) {
    // ECAL barrel
    if (ph.SigmaIEtaIEta() > 0.014) return false;
    if (ph.R9() > 0.9) {
      if (ph.HadOverEm() > 0.082) return false;
    } else {
      if (ph.HadOverEm() > 0.075) return false;
    }
  } else {
    // ECAL endcaps
    if (ph.SigmaIEtaIEta() > 0.034) return false;
    if (ph.HadOverEm() > 0.075) return false;
  }

  // isolation
  if (ph.R9() > 0.9) {
    if (ph.IdMap("HadIso_R03") - 0.005 * ph.Et() > 50) return false;
    if (ph.IdMap("TrkIso_R03") - 0.002 * ph.Et() > 50) return false;
    if (ph.CiCPF4chgpfIso02()[0] > 4) return false;
  } else {
    if (ph.IdMap("HadIso_R03") - 0.005 * ph.Et() > 4) return false;
    if (ph.IdMap("TrkIso_R03") - 0.002 * ph.Et() > 4) return false;
    if (ph.CiCPF4chgpfIso02()[0] > 4) return false;
  }

  // preselection successful
  return true;
}


bool ObjectID::HggPreselection(const TCElectron& el)
{

  if (el.Et() < 15 || el.Et() > 200) return false;
  //if (el.ConversionVeto() == 0) return false;

  // ECAL barrel
  if (fabs(el.SCEta()) < 1.479) {
    if (el.SigmaIEtaIEta() > 0.014) return false;
    if (el.R9() > 0.9) {
      if (el.HadOverEm() > 0.082) return false;
    } else {
      if (el.HadOverEm() > 0.075) return false;
    }
    // ECAL endcaps
  } else {
    if (el.SigmaIEtaIEta() > 0.034) return false;
    if (el.HadOverEm() > 0.075) return false;
  }

  // isolation
  if (el.R9() > 0.9) {
    if (el.IdMap("dr03HcalTowerSumEt") - 0.005 * el.Et() > 50) return false;
    if (el.IdMap("dr03TkSumPt") - 0.002 * el.Et() > 50) return false;
    //if (el.CiCPF4chgpfIso02()[0] > 4) return false;
  } else {
    if (el.IdMap("dr03HcalTowerSumEt") - 0.005 * el.Et() > 4) return false;
    if (el.IdMap("dr03TkSumPt") - 0.002 * el.Et() > 4) return false;
    //if (el.CiCPF4chgpfIso02()[0] > 4) return false;
  }

  // preselection successful
  return true;
}


bool ObjectID::PassPhotonMVA(const TCPhoton& ph, float& mvaScore){

  mvaScore = -99;

  bool mvaPass = false;
  //if (!HggPreselection(ph)) return false;

  // classification variables
  static Float_t phoEt_, phoEta_, phoPhi_, phoR9_;
  static Float_t phoSigmaIEtaIEta_, phoSigmaIEtaIPhi_;
  static Float_t phoS13_, phoS4_, phoS25_, phoSCEta_, phoSCRawE_;
  static Float_t phoSCEtaWidth_, phoSCPhiWidth_, rho2012_;
  static Float_t phoPFPhoIso_, phoPFChIso_, phoPFChIsoWorst_;
  static Float_t phoESEnToRawE_, phoESEffSigmaRR_x_;

  // MVA classifiers for 0=ECAL barrel and 1=ECAL endcaps
  static TMVA::Reader* tmvaReader[2] = {NULL, NULL};

  // 0=ECAL barrel or 1=ECAL endcaps
  int iBE = (fabs(ph.SCEta()) < 1.479) ? 0 : 1;

  // one-time MVA initialization
  if (!tmvaReader[iBE]) {
    tmvaReader[iBE] = new TMVA::Reader("!Color:Silent");

    // add classification variables
    // Orger Matters!!
    tmvaReader[iBE]->AddVariable("phoPhi", &phoPhi_);
    tmvaReader[iBE]->AddVariable("phoR9",  &phoR9_);
    tmvaReader[iBE]->AddVariable("phoSigmaIEtaIEta", &phoSigmaIEtaIEta_);
    tmvaReader[iBE]->AddVariable("phoSigmaIEtaIPhi", &phoSigmaIEtaIPhi_);
    tmvaReader[iBE]->AddVariable("s13",       &phoS13_);
    tmvaReader[iBE]->AddVariable("s4ratio",   &phoS4_);
    tmvaReader[iBE]->AddVariable("s25",       &phoS25_);
    tmvaReader[iBE]->AddVariable("phoSCEta",  &phoSCEta_);
    tmvaReader[iBE]->AddVariable("phoSCRawE", &phoSCRawE_);
    tmvaReader[iBE]->AddVariable("phoSCEtaWidth", &phoSCEtaWidth_);
    tmvaReader[iBE]->AddVariable("phoSCPhiWidth", &phoSCPhiWidth_);

    if (iBE == 1) {
      tmvaReader[iBE]->AddVariable("phoESEn/phoSCRawE", &phoESEnToRawE_);
      tmvaReader[iBE]->AddVariable("phoESEffSigmaRR",   &phoESEffSigmaRR_x_);
    }

    tmvaReader[iBE]->AddVariable("rho2012",     &rho2012_);
    tmvaReader[iBE]->AddVariable("phoPFPhoIso", &phoPFPhoIso_);
    tmvaReader[iBE]->AddVariable("phoPFChIso",  &phoPFChIso_);

    tmvaReader[iBE]->AddVariable("phoPFChIsoWorst", &phoPFChIsoWorst_);

    tmvaReader[iBE]->AddSpectator("phoEt",  &phoEt_);
    tmvaReader[iBE]->AddSpectator("phoEta", &phoEta_);

    // weight files
    if (iBE == 0)
      tmvaReader[0]->BookMVA("BDT", "../data/PhotonMVAID_EB_weights.xml");
    else
      tmvaReader[1]->BookMVA("BDT", "../data/PhotonMVAID_EE_weights.xml");

  } // one-time initialization

  // set MVA variables
  phoSigmaIEtaIEta_ = ph.SigmaIEtaIEta();
  phoSigmaIEtaIPhi_ = ph.SigmaIEtaIPhi();
  phoS4_  = ph.E2x2()/ph.E5x5();
  phoS13_ = ph.E1x3()/ph.E5x5();
  phoS25_ = ph.E2x5Max()/ph.E5x5();
  phoEt_  = ph.Pt();
  phoEta_ = ph.Eta();
  phoPhi_ = ph.Phi();
  phoR9_  = ph.R9();
  phoSCEta_  = ph.SCEta();
  phoSCRawE_ = ph.SCRawEnergy();
  phoSCEtaWidth_ = ph.SCEtaWidth();
  phoSCPhiWidth_ = ph.SCPhiWidth();
  rho2012_       = _rhoFactor;
  phoPFPhoIso_   = ph.PfIsoPhoton();
  phoPFChIso_    = ph.PfIsoCharged();
  phoESEnToRawE_ = ph.PreShowerOverRaw();
  phoESEffSigmaRR_x_= ph.ESEffSigmaRR()[0];

  // evaluate largest isolation value
  phoPFChIsoWorst_ = 0;
  for (size_t k = 0; k < ph.CiCPF4chgpfIso03().size(); ++k){
    if (phoPFChIsoWorst_ < ph.CiCPF4chgpfIso03()[k]) phoPFChIsoWorst_ = ph.CiCPF4chgpfIso03()[k];
  }

  mvaScore = tmvaReader[iBE]->EvaluateMVA("BDT");
  if (mvaScore > _phoMvaCuts[iBE])
    mvaPass = true;


  return mvaPass;

}

void ObjectID::MuonDump(const TCMuon& mu, TVector3 *pv)
{
  Float_t muISO = CalculateMuonIso(mu);
  //Float_t muISO = (mu.IsoMap("pfChargedHadronPt_R04") +
  //               TMath::Max(0.0, mu.IsoMap("pfNeutralHadronEt_R04") + mu.IsoMap("pfPhotonEt_R04") - 0.5*mu.IsoMap("pfPUPt_R04")))/mu.Pt();


  cout  << _runNumber << " " << _eventNumber << "  pt =" << mu.Pt()
        << " eta=" << mu.Eta() << " GLB=" << mu.IsGLB() << "  isPF=" << mu.IsPF() << " isTRK="<<mu.IsTRK()<<endl;
  cout  << "chi2_tr="<< mu.NormalizedChi2_tracker() << "  iso=" << muISO <<  "  dxy=" << mu.Dxy(pv) << " dz=" << mu.Dz(pv)
        << "\n  good =" << mu.IsGood()
    //<< "\n  trk =" << mu.TrackLayersWithMeasurement() << "  pix=" <<mu.PixelLayersWithMeasurement()
    //   << "\n " << mu.NormalizedChi2() << " " << mu.NumberOfValidMuonHits() << " " << mu.NumberOfMatchedStations()
    //  << " " << mu.NumberOfValidPixelHits() << " " << mu.TrackLayersWithMeasurement()
        << endl;
}

void ObjectID::PhotonDump(const TCPhoton& pho, phIdAndIsoCuts cuts)
{
  //Float_t phoISO = 0;
  float chIsoCor,nhIsoCor,phIsoCor;
  float tmpEta = pho.SCEta();

  CalculatePhotonIso(pho, chIsoCor,nhIsoCor,phIsoCor);

  cout  << _runNumber << " " << _eventNumber << " " << pho.Pt() << " " << pho.Eta() << " " << pho.SCEta() <<endl;
  cout<<"conv:"<<pho.ConversionVeto()<<"  hoem:"<< pho.HadOverEm() <<"  sieie:"<<pho.SigmaIEtaIEta()<<endl;
  cout  <<"\n iso: ch:" << chIsoCor<< "  nh:" << nhIsoCor <<  " ph:"<<phIsoCor<<endl;
  if (fabs(tmpEta)  < 1.4442){
    cout<<"\n cuts: "<<cuts.chIso03[0] <<"  "<<cuts.nhIso03[0] + 0.040*pho.Pt()<<"  "<<cuts.phIso03[0] + 0.005*pho.Pt()<<endl;
    cout<<"conv:"<<cuts.PassedEleSafeVeto[0]<<"  hoem:"<< cuts.HadOverEm[0]<<"  sieie:"<<cuts.sigmaIetaIeta[0]<<endl;
  }
  else if (fabs(tmpEta)  > 1.566){
    cout<<"\n cuts: "<<cuts.chIso03[1] <<"  "<<cuts.nhIso03[1] + 0.040*pho.Pt()<<"  "<<cuts.phIso03[1] + 0.005*pho.Pt()<<endl;
    cout<<"conv:"<<cuts.PassedEleSafeVeto[1]<<"  hoem:"<< cuts.HadOverEm[1]<<"  sieie:"<<cuts.sigmaIetaIeta[1]<<endl;
  }

}

void ObjectID::ElectronDump(const TCElectron& el)
{

  vector<TCElectron::Track> trk = el.GetTracks();
  vector<TCEGamma> sc = el.BaseSC();

  if (trk.size() < 2 || sc.size()<2)
    return;

  sort(sc.begin(),  sc.end(),  SCESortCondition);
  sort(trk.begin(), trk.end(), P4SortCondition);

  cout<<"rho         "<<_rhoFactor<<endl;
  cout<<"eleBCS25_1: "<< sc[0].E2x5()/sc[0].E5x5()<<endl;
  cout<<"eleBCS25_2: "<<sc[1].E2x5()/sc[1].E5x5()<<endl;
  cout<<"eleGSFPt2over1: "<<trk[1].Pt()/trk[0].Pt()<<endl;
  cout<<"eleGSFPt2:  "<<trk[1].Pt()<<endl;
  cout<<"eleGSFPt1:  "<<trk[0].Pt()<<endl;
  cout<<"eleGSFdR:   "<<trk[0].DeltaR(trk[1])<<endl;
  cout<<"r9:         "<<el.R9()<<endl;
  cout<<"SCEtaWidth: "<<el.SCEtaWidth()<<endl;
  cout<<"SCPhiWidth: "<<el.SCPhiWidth()<<endl;
  cout<<"HoverE:     "<<el.HadOverEm()<<endl;
  cout<<"EoverPt:    "<<el.SCRawEnergy()/(trk[0]+trk[1]).Pt()<<endl;
  cout<<"E, ERaw Pt: "<<el.SCEnergy() <<" "<<el.SCRawEnergy()<<" "<<(trk[0]+trk[1]).Pt()<<endl;
  cout<<"dEtaAtVtx   "<<el.SCDeltaEta()<<endl;
  cout<<"dPhiAtVtx:  "<<el.SCDeltaPhi()<<endl;
  cout<<"EcalEnPin:  "<<el.InverseEnergyMomentumDiff()<<endl;
  cout<<"recoSCEta:  "<<el.SCEta()<<endl;
  cout<<"eleBCS15_2: "<<sc[1].E1x5()/sc[1].E5x5()<<endl;
  cout<<"eleBCSieie1:"<<sc[0].SigmaIEtaIEta()<<endl;
  cout<<"eleBCSieie2:"<<sc[1].SigmaIEtaIEta()<<endl;
  cout<<"=== BDT: "<<el.IdMap("mvaScore")<<endl;

}

bool ObjectID::PassJetID(const TCJet& jet, int nVtx){
  jetIdCuts cuts = jetIdCutsVBF;
  bool idPass = false;
  if (fabs(jet.Eta()) < 2.5){
    if(
       jet.BetaStarClassic()/log(nVtx-0.64) < cuts.betaStarC[0]
       && jet.DR2Mean() < cuts.dR2Mean[0]
       ) idPass = true;
  }else if (fabs(jet.Eta()) < 2.75){
    if(
       jet.BetaStarClassic()/log(nVtx-0.64) < cuts.betaStarC[1]
       && jet.DR2Mean() < cuts.dR2Mean[1]
       ) idPass = true;
  }else if (fabs(jet.Eta()) < 3.0){
    if(
       jet.DR2Mean() < cuts.dR2Mean[2]
       ) idPass = true;
  }else{
    if(
       jet.DR2Mean() < cuts.dR2Mean[3]
       ) idPass = true;
  }
  return idPass;
}

void ObjectID::JetDump(const TCJet& jet)
{
  cout << "pt\t:" << jet.Pt() << endl;
  cout << "eta\t:" << jet.Eta() << endl;
  cout << "number of constituents\t:" << jet.NumConstit() << endl;
  cout << "neutral hadronic\t:" << jet.NeuHadFrac() << endl;
  cout << "neutral em fraction\t:" << jet.NeuEmFrac() << endl;
  cout << "charged hadronic\t:" << jet.ChHadFrac() << endl;
  cout << "number of charged\t:" << jet.NumChPart() << endl;
  cout << "charged EM fraction\t:" << jet.ChEmFrac() << endl;
  //cout << "CSV b-jet discriminator\t:" << jet.BDiscriminatorMap("CSV") << endl;
}
