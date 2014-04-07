#ifndef _ObjectID_H
#define _ObjectID_H


#include "TLorentzVector.h"
#include "../interface/TCPhoton.h"
#include "../interface/TCEGamma.h"
#include "../interface/TCElectron.h"
#include "../interface/TCPhoton.h"
#include "../interface/TCMuon.h"
#include "../interface/TCGenParticle.h"

struct muIdAndIsoCuts{
  Bool_t IsPF;
  Bool_t IsGLB;
  Float_t ptErrorOverPt;
  Int_t TrackLayersWithMeasurement;
  Int_t PixelLayersWithMeasurement;
  Int_t NumberOfValidMuonHits;
  Int_t NumberOfValidTrackerHits;
  Int_t NumberOfValidPixelHits;
  Int_t NumberOfMatches;
  Int_t NumberOfMatchedStations;
  Float_t NormalizedChi2;
  Float_t NormalizedChi2_tracker;
  Float_t dxy;
  Float_t dz;
  Float_t chIso04;
  Float_t nhIso04;
  Float_t phIso04;
  Float_t pfIso04;
};

struct elIdAndIsoCuts{
  //broken into [0] barrel and [1] endcap
  Float_t ptErrorOverPt[2];
  Float_t dEtaIn[2];
  Float_t dPhiIn[2];
  Float_t sigmaIetaIeta[2];
  Float_t HadOverEm[2];
  Float_t dxy[2];
  Float_t dz[2];
  Float_t fabsEPDiff[2];
  Int_t ConversionMissHits[2];
  Int_t PassedConversionProb[2];
  Float_t pfIso04[2];
};

struct phIdAndIsoCuts{
  //broken into [0] barrel and [1] endcap
  Int_t PassedEleSafeVeto[2];
  Float_t sigmaIetaIeta[2];
  Float_t HadOverEm[2];
  float chIso03[2];
  float nhIso03[2];
  float phIso03[2];
};

muIdAndIsoCuts muIdAndIsoCutsTight, muIdAndIsoCutsLoose, muIdAndIsoCutsSoft;
elIdAndIsoCuts elIdAndIsoCutsTight, elIdAndIsoCutsLoose;
phIdAndIsoCuts phIdAndIsoCutsHZG,   phIdAndIsoCutsTight, phIdAndIsoCutsLoose;



class ObjectID {
 public:
  ObjectID();
  virtual ~ObjectID();

  virtual void CalculatePhotonIso(TCPhoton *ph, float& chIsoCor, float& nhIsoCor, float& phIsoCor);
  bool PassPhotonIdAndIso(TCPhoton *ph, phIdAndIsoCuts cuts, TVector3 *pv);
  virtual void SetEventInfo(Bool_t, UInt_t, ULong64_t, Float_t );
  virtual void PhotonR9Corrector(TCPhoton *ph);

  virtual float CalculateMuonIso(TCMuon *l);
  virtual float CalculateElectronIso(TCElectron *l);
  virtual bool PassMuonIdAndIso(TCMuon *l, muIdAndIsoCuts c, TVector3 *pv);
  virtual bool PassElectronIdAndIso(TCElectron *l, elIdAndIsoCuts c, TVector3 *pv);
  virtual bool PassElectronIdAndIsoMVA(TCElectron *l);

   TCGenParticle * GetPrimaryAncestor(TCGenParticle *p);
   virtual void DiscoverGeneology(TCGenParticle *p);
   virtual void MuonDump(TCMuon mu, TVector3 *pv);
   virtual void PhotonDump(TCPhoton pho,  phIdAndIsoCuts c);

 private:
   Bool_t    _isRealData;
   UInt_t    _runNumber;
   ULong64_t _eventNumber;
   Float_t   _rhoFactor;

};


#endif
