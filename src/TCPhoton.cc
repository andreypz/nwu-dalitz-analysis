#include "TCPhoton.h"
#include<iostream>

TCPhoton::TCPhoton() {
}

TCPhoton::~TCPhoton() {
}



// "get" methods -------------------------------------

TLorentzVector TCPhoton::p4() const {
   return _p4;
}

float TCPhoton::pt() const {
   return _p4.Pt();
}

TVector3 TCPhoton::Vtx() const {
   return _vtx;
}

//TVector3 TCPhoton::AssocVtx() {
//   return _assocPV;
//}



float TCPhoton::eta() const {
 return _p4.Eta();
}

float TCPhoton::phi() const {
  return _p4.Phi();
}

int TCPhoton::charge() const {
   return _charge;
}

float TCPhoton::normChi2() const {
   return _normChi2;
}

float TCPhoton::emIso() const {
   return _emIso;
}

float TCPhoton::hadIso() const {
   return _hadIso;
}

float TCPhoton::trkIso() const {
   return _trkIso;
}

float TCPhoton::hadOverEm() const {
  return _hadOverEm;
}
float TCPhoton::dPhiSuperCluster() const {
  return _dPhiSuperCluster;
}
float TCPhoton::dEtaSuperCluster() const {
  return _dEtaSuperCluster;
}
float TCPhoton::sigmaIetaIeta() const {
  return _sigmaIetaIeta;
}


// "set" methods ---------------------------------------------

void TCPhoton::Setp4(TLorentzVector p4) {
  _p4 = p4;
}

void TCPhoton::Setp4(float px, float py, float pz, float e) {
  TLorentzVector p4(px, py, pz, e);
  _p4 = p4;
}

void TCPhoton::SetVtx(float vx, float vy, float vz) {
  TVector3 v3(vx, vy, vz);
  _vtx = v3;
}

//void TCPhoton::SetAssocVtx(float vx, float vy, float vz) {
//   TVector3 v3(vx, vy, vz);
//   _assocPV = v3;
//}

void TCPhoton::SetCharge(int c){
  _charge = c;
}
//void TCPhoton::Setdxy(float d){
//   _dxy = d;
//}

void TCPhoton::SetNormChi2(float c){
  _normChi2 = c;
}
void TCPhoton::SetEMIso(float e){
  _emIso = e;
}
void TCPhoton::SetHADIso(float h){
  _hadIso = h;
}
void TCPhoton::SetTRKIso(float t){
  _trkIso = t;
}
void TCPhoton::SetHadOverEm(float he){
  _hadOverEm = he;
}
void TCPhoton::SetDPhiSuperCluster(float dp){
  _dPhiSuperCluster = dp;
}
void TCPhoton::SetDEtaSuperCluster(float de){
  _dEtaSuperCluster = de;
}
void TCPhoton::SetSigmaIetaIeta(float sieie){
  _sigmaIetaIeta = sieie;
}


