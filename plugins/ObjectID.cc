#include "ObjectID.h"

float _EAPho[7][3] = {
  {0.012,  0.030,   0.148}, //         eta < 1.0
  {0.010,  0.057,   0.130}, // 1.0   < eta < 1.479
  {0.014,  0.039,   0.112}, // 1.479 < eta < 2.0
  {0.012,  0.015,   0.216}, // 2.0   < eta < 2.2
  {0.016,  0.024,   0.262}, // 2.2   < eta < 2.3
  {0.020,  0.039,   0.260}, // 2.3   < eta < 2.4
  {0.012,  0.072,   0.266}  // 2.4   < eta
};

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
  elIdAndIsoCutsTight.dPhiIn[0] = 0.06;
  elIdAndIsoCutsTight.dEtaIn[0] = 0.004;
  elIdAndIsoCutsTight.sigmaIetaIeta[0] = 0.01;
  elIdAndIsoCutsTight.HadOverEm[0] = 0.12;
  elIdAndIsoCutsTight.fabsEPDiff[0] = 0.05;
  elIdAndIsoCutsTight.dxy[0] = 0.02;
  elIdAndIsoCutsTight.dz[0] = 0.1;
  elIdAndIsoCutsTight.pfIso04[0] = 0.15;

  elIdAndIsoCutsTight.ptErrorOverPt[1] = 9999.;
  elIdAndIsoCutsTight.dPhiIn[1] = 0.03;
  elIdAndIsoCutsTight.dEtaIn[1] = 0.007;
  elIdAndIsoCutsTight.sigmaIetaIeta[1] = 0.03;
  elIdAndIsoCutsTight.HadOverEm[1] = 0.10;
  elIdAndIsoCutsTight.fabsEPDiff[1] = 0.05;
  elIdAndIsoCutsTight.dxy[1] = 0.02;
  elIdAndIsoCutsTight.dz[1] = 0.1;
  elIdAndIsoCutsTight.pfIso04[1] = 0.15;

  //Soft muon id (specific for dalitz)
  muIdAndIsoCutsSoft.TrackLayersWithMeasurement = 5;
  muIdAndIsoCutsSoft.PixelLayersWithMeasurement = 1;
  muIdAndIsoCutsSoft.NormalizedChi2_tracker = 3;
  muIdAndIsoCutsSoft.dxy             = 0.2; //3
  muIdAndIsoCutsSoft.dz              = 0.5; //30
  muIdAndIsoCutsSoft.pfIso04         = 0.4;

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

bool ObjectID::PassPhotonIdAndIso(const TCPhoton& ph, TString n)
{
  phIdAndIsoCuts cuts;

  if (n=="CutBased-MediumWP")
    cuts = phIdAndIsoCutsHZG;
  else if (n=="CutBased-TightWP")
    cuts = phIdAndIsoCutsTight;
  else  {
    cout<<"PhotonID Warning: no such cut "<<n<<endl;
    return false;
  }

  bool pass = false;
  //Float_t phoISO = 0;
  float tmpEta = ph.SCEta();

  float chIsoCor=0, nhIsoCor=0, phIsoCor=0;
  CalculatePhotonIso(ph, chIsoCor,nhIsoCor,phIsoCor);

  if(
     (fabs(tmpEta)  < 1.442
      && ph.ConversionVeto()       == cuts.PassedEleSafeVeto[0]
      && ph.HadOverEm()             < cuts.HadOverEm[0]
      && ph.SigmaIEtaIEta()         < cuts.sigmaIetaIeta[0]
      && max((double)chIsoCor,0.)    < cuts.chIso03[0]
      && max((double)nhIsoCor,0.)    < cuts.nhIso03[0] + 0.040*ph.Pt()
      && max((double)phIsoCor,0.)    < cuts.phIso03[0] + 0.005*ph.Pt()
      ) ||
     (fabs(tmpEta)  > 1.566
      && ph.ConversionVeto()       == cuts.PassedEleSafeVeto[1]
      && ph.HadOverEm()             < cuts.HadOverEm[1]
      && ph.SigmaIEtaIEta()         < cuts.sigmaIetaIeta[1]

      && max((double)chIsoCor,0.)    < cuts.chIso03[1]
      && max((double)nhIsoCor,0.)    < cuts.nhIso03[1] + 0.040*ph.Pt()
      && max((double)phIsoCor,0.)    < cuts.phIso03[1] + 0.005*ph.Pt()
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

  //  Float_t muISO = (lep.IsoMap("pfChargedHadronPt_R04") +
  //TMath::Max(0.0, lep.IsoMap("pfNeutralHadronEt_R04") + lep.IsoMap("pfPhotonEt_R04") - 0.5*lep.IsoMap("pfPUPt_R04")))/lep.Pt();

  return muISO;
}

bool ObjectID::PassMuonIdAndIso(const TCMuon& lep, TVector3 *pv, TString n)
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

  //Float_t muISO =  CalculateMuonIso(lep);

  if(1
     && lep.IsPF()
     && ((lep.IsTRK() && lep.NumberOfMatches() > 0) || lep.IsGLB())
     && fabs(lep.Dxy(pv))     < cuts.dxy
     && fabs(lep.Dz(pv))      < cuts.dz
     //&& (lep.Pt() < 20 || zgamma::CalculateMuonIso(lep) < 0.4)
     )
    pass = true;
  /*
  if(1
     && lep.IsPF() && lep.IsTRK()
     && lep.NormalizedChi2_tracker()  < cuts.NormalizedChi2_tracker
     && fabs(lep.Dxy(pv))     < cuts.dxy
     && fabs(lep.Dz(pv))      < cuts.dz
     )
    pass = true;
  */

  return pass;
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





bool ObjectID::PassElectronIdAndIsoMVA(const TCElectron& lep)
{
  bool pass = false;

  Float_t eleISO = CalculateElectronIso(lep);

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

  if(((fabs(lep.Eta()) < 1.442
       && lep.PtError()/lep.Pt() < cuts.ptErrorOverPt[0]
       && lep.SigmaIEtaIEta() < cuts.sigmaIetaIeta[0]
       && fabs(lep.SCDeltaPhi()) < cuts.dPhiIn[0]
       && fabs(lep.SCDeltaEta()) < cuts.dEtaIn[0]
       && lep.HadOverEm() < cuts.HadOverEm[0]
       && fabs(lep.Dxy(pv)) < cuts.dxy[0]
       && fabs(lep.Dz(pv)) < cuts.dz[0]
       //&& eleISO < cuts.pfIso04[0]
       )||
      (fabs(lep.Eta()) > 1.556
       && lep.PtError()/lep.Pt() < cuts.ptErrorOverPt[1]
       && lep.SigmaIEtaIEta() < cuts.sigmaIetaIeta[1]
       && fabs(lep.SCDeltaPhi()) < cuts.dPhiIn[1]
       && fabs(lep.SCDeltaEta()) < cuts.dEtaIn[1]
       && lep.HadOverEm() < cuts.HadOverEm[1]
       && fabs(lep.Dxy(pv)) < cuts.dxy[1]
       && fabs(lep.Dz(pv)) < cuts.dz[1]
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


void ObjectID::DiscoverGeneology(TCGenParticle *p)
{
  cout<<"---->> event = "<<_eventNumber<<"   "<<p->GetPDGId()<<"  st = "<<p->GetStatus()
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
  if (fabs(tmpEta)  < 1.442){
    cout<<"\n cuts: "<<cuts.chIso03[0] <<"  "<<cuts.nhIso03[0] + 0.040*pho.Pt()<<"  "<<cuts.phIso03[0] + 0.005*pho.Pt()<<endl;
    cout<<"conv:"<<cuts.PassedEleSafeVeto[0]<<"  hoem:"<< cuts.HadOverEm[0]<<"  sieie:"<<cuts.sigmaIetaIeta[0]<<endl;
  }
  else if (fabs(tmpEta)  > 1.566){
    cout<<"\n cuts: "<<cuts.chIso03[1] <<"  "<<cuts.nhIso03[1] + 0.040*pho.Pt()<<"  "<<cuts.phIso03[1] + 0.005*pho.Pt()<<endl;
    cout<<"conv:"<<cuts.PassedEleSafeVeto[1]<<"  hoem:"<< cuts.HadOverEm[1]<<"  sieie:"<<cuts.sigmaIetaIeta[1]<<endl;
  }

}
