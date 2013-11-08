//written by L. Gray

#include "Math/Math.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/DistFunc.h"
#include "RooRealVar.h"
#ifndef __CINT__
#include "RooCFunction1Binding.h"
#include "RooCFunction3Binding.h"
#endif
#include "RooTFnBinding.h"


double beta_pdf_expand(double x, double a, double b) {
  //return std::exp( ROOT::Math::lgamma(a + b) - ROOT::Math::lgamma(a) - ROOT::Math::lgamma(b) + std::log((x-xlow)/range) * (a -1.) + ROOT::Math::log1p(-((x-xlow)/range) ) * (b - 1.) ); 
  return std::exp( ROOT::Math::lgamma(a + b) - ROOT::Math::lgamma(a) - ROOT::Math::lgamma(b) + std::log((x-108)/52.) * (a -1.) + ROOT::Math::log1p(-((x-108)/52.) ) * (b - 1.) ); 
}

RooAbsPdf* makeBetaPdf(const char* name,
           RooRealVar& x, 
           RooRealVar& alpha, 
           RooRealVar& beta) {
  using namespace RooFit;
  //return bindPdf(name,ROOT::Math::beta_pdf,x,alpha,beta);
  return bindPdf(name,beta_pdf_expand,x,alpha,beta);
}

