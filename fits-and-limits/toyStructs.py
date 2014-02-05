#!/usr/bin/env python
import sys
sys.argv.append('-b')
from ROOT import *

def makeToyStucts():
  print 'making TOYDATA..'
  gROOT.ProcessLine(
  'struct TOYDATA{\
    Int_t totalData;\
    Int_t sigWindowData;\
  };')

  print 'making BERN4...'
  gROOT.ProcessLine(
  'struct BERN4{\
    Double_t yieldBkg;\
    Double_t yieldBkgErr;\
    Double_t yieldSig;\
    Double_t yieldSigErr;\
    Double_t paramP1;\
    Double_t paramP1Err;\
    Double_t paramP2;\
    Double_t paramP2Err;\
    Double_t paramP3;\
    Double_t paramP3Err;\
    Double_t paramP4;\
    Double_t paramP4Err;\
    Double_t edm;\
    Double_t minNll;\
    Int_t statusAll;\
    Int_t statusMIGRAD;\
    Int_t statusHESSE;\
    Int_t covQual;\
    Int_t numInvalidNLL;\
  };')

  print 'making BERN5...'
  gROOT.ProcessLine(
  'struct BERN5{\
    Double_t yieldBkg;\
    Double_t yieldBkgErr;\
    Double_t yieldSig;\
    Double_t yieldSigErr;\
    Double_t paramP1;\
    Double_t paramP1Err;\
    Double_t paramP2;\
    Double_t paramP2Err;\
    Double_t paramP3;\
    Double_t paramP3Err;\
    Double_t paramP4;\
    Double_t paramP4Err;\
    Double_t paramP5;\
    Double_t paramP5Err;\
    Double_t edm;\
    Double_t minNll;\
    Int_t statusAll;\
    Int_t statusMIGRAD;\
    Int_t statusHESSE;\
    Int_t covQual;\
    Int_t numInvalidNLL;\
  };')

  print 'making BERN6...'
  gROOT.ProcessLine(
  'struct BERN6{\
    Double_t yieldBkg;\
    Double_t yieldBkgErr;\
    Double_t yieldSig;\
    Double_t yieldSigErr;\
    Double_t paramP1;\
    Double_t paramP1Err;\
    Double_t paramP2;\
    Double_t paramP2Err;\
    Double_t paramP3;\
    Double_t paramP3Err;\
    Double_t paramP4;\
    Double_t paramP4Err;\
    Double_t paramP5;\
    Double_t paramP5Err;\
    Double_t paramP6;\
    Double_t paramP6Err;\
    Double_t edm;\
    Double_t minNll;\
    Int_t statusAll;\
    Int_t statusMIGRAD;\
    Int_t statusHESSE;\
    Int_t covQual;\
    Int_t numInvalidNLL;\
  };')



  print 'making GAUSSBERN3...'
  gROOT.ProcessLine(
  'struct GAUSSBERN3{\
    Double_t yieldBkg;\
    Double_t yieldBkgErr;\
    Double_t yieldSig;\
    Double_t yieldSigErr;\
    Double_t paramSigma;\
    Double_t paramSigmaErr;\
    Double_t paramStep;\
    Double_t paramStepErr;\
    Double_t paramP1;\
    Double_t paramP1Err;\
    Double_t paramP2;\
    Double_t paramP2Err;\
    Double_t paramP3;\
    Double_t paramP3Err;\
    Double_t edm;\
    Double_t minNll;\
    Int_t statusAll;\
    Int_t statusMIGRAD;\
    Int_t statusHESSE;\
    Int_t covQual;\
    Int_t numInvalidNLL;\
  };')

  print 'making GAUSSBERN4...'
  gROOT.ProcessLine(
  'struct GAUSSBERN4{\
    Double_t yieldBkg;\
    Double_t yieldBkgErr;\
    Double_t yieldSig;\
    Double_t yieldSigErr;\
    Double_t paramSigma;\
    Double_t paramSigmaErr;\
    Double_t paramStep;\
    Double_t paramStepErr;\
    Double_t paramP1;\
    Double_t paramP1Err;\
    Double_t paramP2;\
    Double_t paramP2Err;\
    Double_t paramP3;\
    Double_t paramP3Err;\
    Double_t paramP4;\
    Double_t paramP4Err;\
    Double_t edm;\
    Double_t minNll;\
    Int_t statusAll;\
    Int_t statusMIGRAD;\
    Int_t statusHESSE;\
    Int_t covQual;\
    Int_t numInvalidNLL;\
  };')

  print 'making GAUSSBERN5...'
  gROOT.ProcessLine(
  'struct GAUSSBERN5{\
    Double_t yieldBkg;\
    Double_t yieldBkgErr;\
    Double_t yieldSig;\
    Double_t yieldSigErr;\
    Double_t paramSigma;\
    Double_t paramSigmaErr;\
    Double_t paramStep;\
    Double_t paramStepErr;\
    Double_t paramP1;\
    Double_t paramP1Err;\
    Double_t paramP2;\
    Double_t paramP2Err;\
    Double_t paramP3;\
    Double_t paramP3Err;\
    Double_t paramP4;\
    Double_t paramP4Err;\
    Double_t paramP5;\
    Double_t paramP5Err;\
    Double_t edm;\
    Double_t minNll;\
    Int_t statusAll;\
    Int_t statusMIGRAD;\
    Int_t statusHESSE;\
    Int_t covQual;\
    Int_t numInvalidNLL;\
  };')


  print 'making GEN...'
  gROOT.ProcessLine(
  'struct GEN{\
    Double_t yieldBkg;\
    Double_t yieldBkgErr;\
    Double_t yieldSig;\
    Double_t yieldSigErr;\
    Double_t edm;\
    Double_t minNll;\
    Int_t statusAll;\
    Int_t statusMIGRAD;\
    Int_t statusHESSE;\
    Int_t covQual;\
    Int_t numInvalidNLL;\
  };')



  print 'making TRUTH'
  gROOT.ProcessLine(
  'struct TRUTH{\
    Double_t yieldBkg;\
    Double_t yieldBkgErr;\
    Double_t yieldSig;\
    Double_t yieldSigErr;\
  };')

