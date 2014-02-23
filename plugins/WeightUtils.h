/*
   Utilities for retrieving weights for PU,etc.
 */

#ifndef _WeightUtils_H
#define _WeightUtils_H

#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH2D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TLorentzVector.h"

using namespace std;

//extern"C" {
//  void pwhg_cphto_reweight_(double *mh, double *gh, double *mt, int *BWflag, double *m, double *w);
//}

class WeightUtils: public TObject {
 public:
  WeightUtils() {};
  virtual ~WeightUtils() {};
  WeightUtils(string sampleName, string dataPeriod, string selection, bool isRealData);
  void  Initialize();
  void  SetDataBit(bool isRealData);
  void  SetDataPeriod(string dataPeriod);
  void  SetSampleName(string sampleName);
  void  SetSelection(string selection);
  
  //Used in zgamma
  float PhotonSF(TLorentzVector ph);
  float MuonSF(TLorentzVector mu);
  float PUWeight(float);

  float RecoWeight(TLorentzVector l1, TLorentzVector l2);
  //float GammaWeight(int nPV, int nJets, TLorentzVector p1);
  float ZZWeight(TLorentzVector l1, TLorentzVector l2);
  float GluGluHiggsWeight(float higgsPt, int higgsMass);
  //float VBFHiggsWeight(float genMass, int higgsMass);
  //float RecoilWeight(TLorentzVector recoilP4, TLorentzVector ZP4);
  float GetTotalWeight(float nPV, int nJets, TLorentzVector l1, TLorentzVector l2);
  
  float HiggsMassLineShapeWeight(float, float);
  float GetElectronEff(TLorentzVector l1) const;
  float GetMuTriggerEff(TLorentzVector l1) const;
  float GetPhotonMass() const;
  
  float GetGammaPtWeight() {return _gammaPtWeight;}
  
  ClassDef(WeightUtils, 0);
  
 private:
  //input parameters
  string _dataPeriod;
  string _sampleName;
  string _selection;
  bool   _isRealData;
  
  //sources
  TFile *_phFile;
  TFile *_muFileID;
  TFile *_muFileISO;
  TFile *_puFile;
 
  TH2D * h2_photonIDSF;
  TH2D * h2_photonCSEVSF;

  Double_t *_muonIDSF_eta09;
  Double_t *_muonIDSF_eta12;
  Double_t *_muonIDSF_eta21;
  Double_t *_muonIDSF_eta24;
  TGraph * gr_muonIDSF_eta09;
  TGraph * gr_muonIDSF_eta12;
  TGraph * gr_muonIDSF_eta21;
  TGraph * gr_muonIDSF_eta24;

  Double_t *_muonISOSF_eta09;
  Double_t *_muonISOSF_eta12;
  Double_t *_muonISOSF_eta21;
  Double_t *_muonISOSF_eta24;
  TGraph * gr_muonISOSF_eta09;
  TGraph * gr_muonISOSF_eta12;
  TGraph * gr_muonISOSF_eta21;
  TGraph * gr_muonISOSF_eta24;
  //TGraph * gr_muonISOSF;
 //TH1D  *h1_eGammaPt;   
  //TH1D  *h1_muGammaPt;  
  //TH1D  *h1_eGammaPV;   
  //TH1D  *h1_muGammaPV;  
  //TH1D  *h1_eGammaMass; 
  //TH1D  *h1_muGammaMass;
  TH1D  *h1_puReweight2011;
  TH1D  *h1_puReweight2012;
  
  TH1D  *h1_Higgs[8];
  
  //weights
  float _photonIDWeight;
  float _puWeight;
  float _zzWeight;
  float _glugluWeight;
  float _vbfWeight;
  float _recoWeight;
  float _triggerWeight;
  float _gammaPtWeight;
  float _gammaPVWeight;
  float _gammaJetWeight;
  float _recoilLongWeight;
  float _recoilTransWeight;
  
  //Misc
  TRandom3* rnGen;
};

#endif

#if !defined(__CINT__)
ClassImp(WeightUtils);
#endif
