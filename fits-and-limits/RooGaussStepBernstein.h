/*****************************************************************************
* Project: RooFit *
* Package: RooFitModels *
* File: $Id: RooStepBernstein.h 28259 2009-04-16 16:21:16Z wouter $
* Authors: *
* Kyle Cranmer (L. Gray for step function and gaussian convoluion piece)
* *
* *
* Redistribution and use in source and binary forms, *
* with or without modification, are permitted according to the terms *
* listed in LICENSE (http://roofit.sourceforge.net/license.txt) *
*****************************************************************************/
#ifndef ROO_GAUSSSTEPBERNSTEIN
#define ROO_GAUSSSTEPBERNSTEIN

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooRealVar;
class RooArgList ;

class RooGaussStepBernstein : public RooAbsPdf {
public:

  RooGaussStepBernstein() ;
  RooGaussStepBernstein(const char *name, const char *title,
      RooAbsReal& _x, RooAbsReal& _mean,
      RooAbsReal& _sigma, RooAbsReal& _stepVal,
      const RooArgList& _coefList) ;

  RooGaussStepBernstein(const RooGaussStepBernstein& other,
const char* name = 0);
  virtual TObject* clone(const char* newname) const
  { return new RooGaussStepBernstein(*this, newname); }
  inline virtual ~RooGaussStepBernstein() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

private:

  RooRealProxy _x;
  RooRealProxy _mean,_sigma,_stepVal;
  RooListProxy _coefList ;

  Double_t evaluate() const;

  // Bernstein polynomial PDF with step function convoluted with gaussian
  ClassDef(RooGaussStepBernstein,1)
};

#endif
