#ifndef _NUELECTRON_H
#define	_NUELECTRON_H

#include "TObject.h"
#include "TLorentzVector.h"

class TCElectron : public TObject {
private:
    TLorentzVector _p4;
    TVector3 _vtx;
    //    TVector3 _assocPV;
    int _charge;
    float _normChi2;
    float _emIso;
    float _hadIso;
    float _trkIso;
    float _hadOverEm;
    float _dPhiSuperCluster;
    float _dEtaSuperCluster;
    float _sigmaIetaIeta;

public:
    TCElectron();
    virtual ~TCElectron();

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

    ClassDef(TCElectron, 1);

};

#endif	/* _NUELECTRON_H */


