#ifndef _NUPHOTON_H
#define	_NUPHOTON_H

#include "TObject.h"
#include "TLorentzVector.h"

class TCPhoton : public TObject {
private:
    TLorentzVector _p4;
    TVector3 _vtx;
    //    TVector3 _assocPV;
    int _charge;
    float _normChi2;
    float _emIso; // use
    float _hadIso; // use
    float _trkIso; // use
    float _hadOverEm; // use
    float _dPhiSuperCluster;
    float _dEtaSuperCluster;
    float _sigmaIetaIeta; // use
    float _sigmaiphiiphi, _e2overe9, _etasupercluster;
    bool _trackveto;

public:
    TCPhoton();
    virtual ~TCPhoton();

    // "get" methods -----------

    TLorentzVector p4() const;
    float pt() const;
    TVector3 Vtx() const;
    float eta() const;
    float phi() const;
    int charge() const;
    float normChi2() const;
    float emIso() const;
    float hadIso() const;
    float trkIso() const;
    float hadOverEm() const;
    float dPhiSuperCluster() const;
    float dEtaSuperCluster() const;
    float sigmaIetaIeta() const;

    float sigmaIphiIphi() const {return _sigmaiphiiphi;}
    float e2overe9() const {return _e2overe9;}
    float etaSupercluster() const {return _etasupercluster;}
    bool trackVeto() const {return _trackveto;}

    //    TVector3 AssocVtx() const;

    // "set" methods ---------
    void Setp4(TLorentzVector p4);
    void Setp4(float px, float py, float pz, float e);
    void SetVtx(float vx, float vy, float vz);
    //  void SetAssocVtx(float vx, float vy, float vz);

    void SetCharge(int c);
    void SetNormChi2(float c);
    void SetEMIso(float e);
    void SetHADIso(float h);
    void SetTRKIso(float t);
 
    void SetHadOverEm(float h);
    void SetDPhiSuperCluster(float dp);
    void SetDEtaSuperCluster(float de);
    void SetSigmaIetaIeta(float sieie);

    void SetSigmaIphiIphi(float pw) {_sigmaiphiiphi = pw;}
    void Sete2overe9(float e2e9) {_e2overe9 = e2e9;}
    void SetEtaSupercluster(float es) {_etasupercluster = es;}
    void SetTrackVeto(bool tv) {_trackveto = tv;}

    ClassDef(TCPhoton, 1);

};

#endif	/* _NUPHOTON_H */


