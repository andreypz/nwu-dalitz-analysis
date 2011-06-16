#include "TCElectron.h"
#include<iostream>

TCElectron::TCElectron() {
}

TCElectron::~TCElectron() {
}



// "get" methods -------------------------------------

TLorentzVector TCElectron::p4() const {
   return _p4;
}

float TCElectron::pt() const {
   return _p4.Pt();
}

TVector3 TCElectron::Vtx() const {
   return _vtx;
}

//TVector3 TCElectron::AssocVtx() {
//   return _assocPV;
//}



float TCElectron::eta() const {
 return _p4.Eta();
}

float TCElectron::phi() const {
  return _p4.Phi();
}

int TCElectron::charge() const {
   return _charge;
}

float TCElectron::normChi2() const {
   return _normChi2;
}

float TCElectron::emIso() const {
   return _emIso;
}

float TCElectron::hadIso() const {
   return _hadIso;
}

float TCElectron::trkIso() const {
   return _trkIso;
}

float TCElectron::hadOverEm() const {
  return _hadOverEm;
}
float TCElectron::dPhiSuperCluster() const {
  return _dPhiSuperCluster;
}
float TCElectron::dEtaSuperCluster() const {
  return _dEtaSuperCluster;
}
float TCElectron::sigmaIetaIeta() const {
  return _sigmaIetaIeta;
}


// "set" methods ---------------------------------------------

void TCElectron::Setp4(TLorentzVector p4) {
  _p4 = p4;
}

void TCElectron::Setp4(float px, float py, float pz, float e) {
  TLorentzVector p4(px, py, pz, e);
  _p4 = p4;
}

void TCElectron::SetVtx(float vx, float vy, float vz) {
  TVector3 v3(vx, vy, vz);
  _vtx = v3;
}

//void TCElectron::SetAssocVtx(float vx, float vy, float vz) {
//   TVector3 v3(vx, vy, vz);
//   _assocPV = v3;
//}

void TCElectron::SetCharge(int c){
  _charge = c;
}
//void TCElectron::Setdxy(float d){
//   _dxy = d;
//}

void TCElectron::SetNormChi2(float c){
  _normChi2 = c;
}
void TCElectron::SetEMIso(float e){
  _emIso = e;
}
void TCElectron::SetHADIso(float h){
  _hadIso = h;
}
void TCElectron::SetTRKIso(float t){
  _trkIso = t;
}
void TCElectron::SetHadOverEm(float he){
  _hadOverEm = he;
}
void TCElectron::SetDPhiSuperCluster(float dp){
  _dPhiSuperCluster = dp;
}
void TCElectron::SetDEtaSuperCluster(float de){
  _dEtaSuperCluster = de;
}
void TCElectron::SetSigmaIetaIeta(float sieie){
  _sigmaIetaIeta = sieie;
}


