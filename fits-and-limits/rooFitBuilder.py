#!/usr/bin/env python
import sys
sys.argv.append('-b')
from ROOT import gSystem
gSystem.Load("libRooFit")
from ROOT import *
import numpy as np

gSystem.SetIncludePath( "-I$ROOFITSYS/include/" );
gROOT.ProcessLine('.x RooStepBernstein.cxx+')
gROOT.ProcessLine('.x RooGaussStepBernstein.cxx+')
#gROOT.ProcessLine('.x HZGRooPdfs.cxx++')

def BuildBetaAndBern(year,lepton,cat,mzg,rangeName,frac = 0.1, fracLow = 0, fracHigh = 0.9):
  suffix = '_'.join([year,lepton,'cat'+cat])
  beta = BuildBetaFunc(year,lepton,cat+'BB',mzg,rangeName)
  bern = BuildBern4(year,lepton,cat+'BB',mzg)
  fracVar = RooRealVar('fracBB_'+suffix,'fracBB_'+suffix,frac,fracLow,fracHigh)
  bbArgs = RooArgList(beta,bern)
  fracArg = RooArgList(fracVar)

  BB = RooAddPdf('BB_'+suffix,'BB_'+suffix,bbArgs,fracArg,True)
  #BB = RooFFTConvPdf('BB_'+suffix,'BB_'+suffix,mzg,beta,bern)
  SetOwnership(beta,0)
  SetOwnership(bern,0)
  SetOwnership(fracVar,0)
  return BB

def BuildGaussAndBern(year,lepton,cat,mzg,rangeName,frac = 0.1, fracLow = 0, fracHigh = 0.9):
  suffix = '_'.join([year,lepton,'cat'+cat])
  gauss = BuildRooGaussian(year,lepton,cat+'GB',mzg)
  bern = BuildBern3(year,lepton,cat+'GB',mzg)
  fracVar = RooRealVar('fracGB_'+suffix,'fracGB_'+suffix,frac,fracLow,fracHigh)
  gbArgs = RooArgList(gauss,bern)
  fracArg = RooArgList(fracVar)

  GB = RooAddPdf('GB_'+suffix,'GB_'+suffix,gbArgs,fracArg,True)
  #BB = RooFFTConvPdf('BB_'+suffix,'BB_'+suffix,mzg,beta,bern)
  SetOwnership(gauss,0)
  SetOwnership(bern,0)
  SetOwnership(fracVar,0)
  return GB


def BuildBetaFunc(year,lepton,cat,mzg,rangeName,alpha = 2, alphaLow = 1, alphaHigh = 10, beta = 5, betaLow = 1, betaHigh = 10):
  suffix = '_'.join([year,lepton,'cat'+cat])
  alphaVar = RooRealVar('alphaBetaFunc_'+suffix,'alphaBetaFunc_'+suffix, alpha, alphaLow, alphaHigh)
  betaVar = RooRealVar('betaBetaFunc_'+suffix,'betaBetaFunc_'+suffix, beta, betaLow, betaHigh)
  xLow = mzg.getMin(rangeName)
  xRange = mzg.getMax(rangeName) - mzg.getMin(rangeName)
  gROOT.ProcessLine('.L betaWrapper.cxx+')
  from ROOT import makeBetaPdf
  BetaFunc = makeBetaPdf('BetaFunc_'+suffix,mzg,alphaVar,betaVar)

  SetOwnership(alphaVar,0)
  SetOwnership(betaVar,0)
  return BetaFunc

def BuildKumaraswamy(year,lepton,cat,mzg,rangeName,alpha = 2, alphaLow = 1, alphaHigh = 10, beta = 5, betaLow = 1, betaHigh = 10):
  suffix = '_'.join([year,lepton,'cat'+cat])
  alphaVar = RooRealVar('alphaKumaraswamy_'+suffix,'alphaKumaraswamy_'+suffix, alpha, alphaLow, alphaHigh)
  betaVar = RooRealVar('betaKumaraswamy_'+suffix,'betaKumaraswamy_'+suffix, beta, betaLow, betaHigh)
  xLow = RooRealVar('xLowKumaraswamy_'+suffix, 'xLowKumaraswamy_'+suffix, mzg.getMin(rangeName))
  xRange = RooRealVar('xRangeKumaraswamy_'+suffix, 'xReangeKumaraswamy_'+suffix,mzg.getMax(rangeName) - mzg.getMin(rangeName))
  Kumaraswamy = RooGenericPdf('Kumaraswamy_'+suffix, 'Kumaraswamy_'+suffix, '(@0>@3)*(@0<(@3+@4))*@1*@2*((@0-@3)/@4)**(@1-1)*(1-((@0-@3)/@4)**@1)**(@2)', RooArgList(mzg,alphaVar,betaVar,xLow,xRange))

  SetOwnership(alphaVar,0)
  SetOwnership(betaVar,0)
  SetOwnership(xLow,0)
  SetOwnership(xRange,0)
  return Kumaraswamy

def BuildGaussExp(year,lepton,cat,mzg,mean = 120, meanLow = 90, meanHigh = 150, sigma = 1, sigmaLow = 0.01, sigmaHigh = 10, tau = 5, tauLow = 0, tauHigh = 50):
  suffix = '_'.join([year,lepton,'cat'+cat])
  meanVar = RooRealVar('meanGaussExp_'+suffix,'meanGaussExp_'+suffix, mean, meanLow, meanHigh)
  sigmaVar = RooRealVar('sigmaGaussExp_'+suffix,'sigmaGaussExp_'+suffix,sigma,sigmaLow,sigmaHigh)
  tauVar = RooRealVar('tauGaussExp_'+suffix,'tauGaussExp_'+suffix,tau,tauLow,tauHigh)

  turnOn = RooGaussModel('turnOnGaussExp_'+suffix,'turnOnGaussExp_'+suffix,mzg,meanVar,sigmaVar)
  GaussExp = RooDecay('GaussExp_'+suffix,'GaussExp_'+suffix,mzg,tauVar,turnOn,RooDecay.SingleSided)

  SetOwnership(meanVar,0)
  SetOwnership(sigmaVar,0)
  SetOwnership(tauVar,0)
  SetOwnership(turnOn,0)
  return GaussExp
  #sprintf(sbuffer1, "GaussExpFullExt_cat%i",i+1);
  #GaussExpFullExt[i] = new RooExtendPdf(sbuffer1,sbuffer1,*GaussExpFull[i],*nGaussExp[i]);
  #nGaussExp[i] = new RooRealVar(sbuffer1,sbuffer1,ntcat[i]->GetEntries(),0,3*ntcat[i]->GetEntries());


def BuildGaussPow(year,lepton,cat,mzg,mean = 0, sigma = 1, sigmaLow = 0.01, sigmaHigh = 10, alpha = 115, alphaLow = 50, alphaHigh = 200,beta = 2.7, betaLow = 0, betaHigh = 20):
  suffix = '_'.join([year,lepton,'cat'+cat])
  meanVar = RooRealVar('meanGaussPow_'+suffix,'meanGaussPow_'+suffix, mean)
  sigmaVar = RooRealVar('sigmaGaussPow_'+suffix,'sigmaGaussPow_'+suffix,sigma,sigmaLow,sigmaHigh)
  alphaVar = RooRealVar('alphaGaussPow_'+suffix,'alphaGaussPow_'+suffix,alpha,alphaLow,alphaHigh)
  betaVar = RooRealVar('betaGaussPow_'+suffix,'betaGaussPow_'+suffix,beta,betaLow,betaHigh)

  turnOn = RooGaussModel('turnOnGaussPow_'+suffix,'turnOnGaussPow_'+suffix,mzg,meanVar,sigmaVar)
  tail = RooGenericPdf('tailGaussPow_'+suffix,'tailGaussPow_'+suffix,'1e-20 + (@0 > @1)*((@0)^(-@2))',RooArgList(mzg,alphaVar,betaVar))
  GaussPow = RooFFTConvPdf('GaussPow_'+suffix,'GaussPow_'+suffix, mzg, tail, turnOn)

  SetOwnership(meanVar,0)
  SetOwnership(sigmaVar,0)
  SetOwnership(alphaVar,0)
  SetOwnership(betaVar,0)
  SetOwnership(turnOn,0)
  SetOwnership(tail,0)
  return GaussPow

def BuildSechExp(year,lepton,cat,mzg,mean = 0, sigma = 5, sigmaLow = 0.01, sigmaHigh = 20, tau = 35, tauLow = 0, tauHigh = 100, alpha = 105, alphaLow = 50, alphaHigh = 200):
  suffix = '_'.join([year,lepton,'cat'+cat])
  meanVar = RooRealVar('meanSechExp_'+suffix,'meanSechExp_'+suffix, mean)
  sigmaVar = RooRealVar('sigmaSechExp_'+suffix,'sigmaSechExp_'+suffix,sigma,sigmaLow,sigmaHigh)
  tauVar = RooRealVar('tauSechExp_'+suffix,'tauSechExp_'+suffix,tau,tauLow,tauHigh)
  alphaVar = RooRealVar('alphaSechExp_'+suffix,'alphaSechExp_'+suffix,alpha,alphaLow,alphaHigh)

  turnOn  = RooGenericPdf('turnOnSechExp_'+suffix,'turnOnSechExp_'+suffix, 'exp(-(@0-@1)/@2)/(@2*(1.0+exp(-(@0-@1)/@2))**2)',RooArgList(mzg,meanVar,sigmaVar))
  tail    = RooGenericPdf('tailSechExp_'+suffix,'tailSechExp_'+suffix,'1e-20 + (@0 > @1)*(exp(-@0/@2))',RooArgList(mzg,alphaVar,tauVar))
  SechExp = RooFFTConvPdf('SechExp_'+suffix,'SechExp_'+suffix,mzg,tail,turnOn)

  SetOwnership(meanVar,0)
  SetOwnership(sigmaVar,0)
  SetOwnership(alphaVar,0)
  SetOwnership(tauVar,0)
  SetOwnership(turnOn,0)
  SetOwnership(tail,0)
  return SechExp

def BuildSechPow(year,lepton,cat,mzg,mean = 0, sigma = 4, sigmaLow = 0.01, sigmaHigh = 20, alpha = 107, alphaLow = 50, alphaHigh = 200, beta = 5, betaLow = 0, betaHigh = 20):
  suffix = '_'.join([year,lepton,'cat'+cat])
  meanVar = RooRealVar('meanSechPow_'+suffix,'meanSechPow_'+suffix, mean)
  sigmaVar = RooRealVar('sigmaSechPow_'+suffix,'sigmaSechPow_'+suffix,sigma,sigmaLow,sigmaHigh)
  alphaVar = RooRealVar('alphaSechPow_'+suffix,'alphaSechPow_'+suffix,alpha,alphaLow,alphaHigh)
  betaVar = RooRealVar('betaSechPow_'+suffix,'betaSechPow_'+suffix,beta,betaLow,betaHigh)

  turnOn = RooGenericPdf('turnOnSechPow_'+suffix,'turnOnSechPow_'+suffix,'exp(-(@0-@1)/@2)/(@2*(1.0+exp(-(@0-@1)/@2))**2)',RooArgList(mzg,meanVar,sigmaVar))
  tail = RooGenericPdf('tailSechPow_'+suffix,'tailSechPow_'+suffix,'1e-20 + (@0 > @1)*((@0)^(-@2))',RooArgList(mzg,alphaVar,betaVar))
  SechPow = RooFFTConvPdf('SechPow_'+suffix,'SechPow_'+suffix, mzg, tail, turnOn)

  SetOwnership(meanVar,0)
  SetOwnership(sigmaVar,0)
  SetOwnership(alphaVar,0)
  SetOwnership(betaVar,0)
  SetOwnership(turnOn,0)
  SetOwnership(tail,0)
  return SechPow

def BuildGaussStepBern3(year,lepton,cat,mzg,mean = 0, sigma = 4, sigmaLow = 0.01, sigmaHigh = 20, step = 115, stepLow = 100, stepHigh = 130,
    p0 = 15, p1 = 0.3, p1Low = -1e-6, p1High = 900,p2 = 0.3, p2Low = -1e-6, p2High = 900,p3 = 0.3, p3Low = -1e-6, p3High = 900):
  suffix = '_'.join([year,lepton,'cat'+cat])
  meanVar = RooRealVar('meanGaussBern3_'+suffix,'meanGaussBern3_'+suffix, mean)
  sigmaVar = RooRealVar('sigmaGaussBern3_'+suffix,'sigmaGaussBern3_'+suffix,sigma,sigmaLow,sigmaHigh)
  stepVar = RooRealVar('stepGaussBern3_'+suffix,'stepGaussBern3_'+suffix,step,stepLow,stepHigh)
  p0Var = RooRealVar('p0GaussBern3_'+suffix,'p0GaussBern3_'+suffix, p0)
  p1Var = RooRealVar('p1GaussBern3_'+suffix,'p1GaussBern3_'+suffix,p1,p1Low,p1High)
  p2Var = RooRealVar('p2GaussBern3_'+suffix,'p2GaussBern3_'+suffix,p2,p2Low,p2High)
  p3Var = RooRealVar('p3GaussBern3_'+suffix,'p3GaussBern3_'+suffix,p3,p3Low,p3High)

  pArgs = RooArgList(p0Var,p1Var,p2Var,p3Var)
  GaussBern3 = RooGaussStepBernstein('GaussBern3_'+suffix,'GaussBern3_'+suffix,mzg,meanVar,sigmaVar,stepVar,pArgs)

  SetOwnership(meanVar,0)
  SetOwnership(sigmaVar,0)
  SetOwnership(stepVar,0)
  SetOwnership(p0Var,0)
  SetOwnership(p1Var,0)
  SetOwnership(p2Var,0)
  SetOwnership(p3Var,0)
  return GaussBern3

def BuildGaussStepBern4(year,lepton,cat,mzg,mean = 0, sigma = 4, sigmaLow = 0.001, sigmaHigh = 60, step = 115, stepLow = 100, stepHigh = 130,
    p0 = 15, p1 = 0.4, p1Low = -1e-6, p1High = 900,p2 = 0.4, p2Low = -1e-6, p2High = 900,p3 = 0.4, p3Low = -1e-6, p3High = 900, p4 = 0.4, p4Low = -1e-6, p4High = 900):
  #def BuildGaussStepBern4(year,lepton,cat,mzg,mean = 0, sigma = 4, sigmaLow = 0.001, sigmaHigh = 60, step = 115, stepLow = 100, stepHigh = 130,
    #p0 = 15, p1 = 0.4, p1Low = -1e2, p1High = 900,p2 = 0.4, p2Low = -1e2, p2High = 900,p3 = 0.4, p3Low = -1e2, p3High = 900, p4 = 0.4, p4Low = -1e2, p4High = 900):
  suffix = '_'.join([year,lepton,'cat'+cat])
  meanVar = RooRealVar('meanGaussBern4_'+suffix,'meanGaussBern4_'+suffix, mean)
  sigmaVar = RooRealVar('sigmaGaussBern4_'+suffix,'sigmaGaussBern4_'+suffix,sigma,sigmaLow,sigmaHigh)
  stepVar = RooRealVar('stepGaussBern4_'+suffix,'stepGaussBern4_'+suffix,step,stepLow,stepHigh)
  p0Var = RooRealVar('p0GaussBern4_'+suffix,'p0GaussBern4_'+suffix, p0)
  p1Var = RooRealVar('p1GaussBern4_'+suffix,'p1GaussBern4_'+suffix,p1,p1Low,p1High)
  p2Var = RooRealVar('p2GaussBern4_'+suffix,'p2GaussBern4_'+suffix,p2,p2Low,p2High)
  p3Var = RooRealVar('p3GaussBern4_'+suffix,'p3GaussBern4_'+suffix,p3,p3Low,p3High)
  p4Var = RooRealVar('p4GaussBern4_'+suffix,'p4GaussBern4_'+suffix,p4,p4Low,p4High)

  pArgs = RooArgList(p0Var,p1Var,p2Var,p3Var,p4Var)
  GaussBern4 = RooGaussStepBernstein('GaussBern4_'+suffix,'GaussBern4_'+suffix,mzg,meanVar,sigmaVar,stepVar,pArgs)

  SetOwnership(meanVar,0)
  SetOwnership(sigmaVar,0)
  SetOwnership(stepVar,0)
  SetOwnership(p0Var,0)
  SetOwnership(p1Var,0)
  SetOwnership(p2Var,0)
  SetOwnership(p3Var,0)
  SetOwnership(p4Var,0)
  return GaussBern4

def BuildGaussStepBern5(year,lepton,cat,mzg,mean = 0, sigma = 5, sigmaLow = 0.001, sigmaHigh = 60, step = 115, stepLow = 100, stepHigh = 130,
    p0 = 15, p1 = 0.5, p1Low = -1e-6, p1High = 900,p2 = 0.5, p2Low = -1e-6, p2High = 900,p3 = 0.5, p3Low = -1e-6, p3High = 900, p4 = 0.5, p4Low = -1e-6, p4High = 900, p5 = 0.5, p5Low = -1e-6, p5High = 900):
  #def BuildGaussStepBern5(year,lepton,cat,mzg,mean = 0, sigma = 5, sigmaLow = 0.001, sigmaHigh = 60, step = 115, stepLow = 100, stepHigh = 130,
   # p0 = 15, p1 = 0.5, p1Low = -1e2, p1High = 900,p2 = 0.5, p2Low = -1e2, p2High = 900,p3 = 0.5, p3Low = -1e2, p3High = 900, p4 = 0.5, p4Low = -1e2, p4High = 900, p5 = 0.5, p5Low = -1e2, p5High = 900):
  suffix = '_'.join([year,lepton,'cat'+cat])
  meanVar = RooRealVar('meanGaussBern5_'+suffix,'meanGaussBern5_'+suffix, mean)
  sigmaVar = RooRealVar('sigmaGaussBern5_'+suffix,'sigmaGaussBern5_'+suffix,sigma,sigmaLow,sigmaHigh)
  stepVar = RooRealVar('stepGaussBern5_'+suffix,'stepGaussBern5_'+suffix,step,stepLow,stepHigh)
  p0Var = RooRealVar('p0GaussBern5_'+suffix,'p0GaussBern5_'+suffix, p0)
  p1Var = RooRealVar('p1GaussBern5_'+suffix,'p1GaussBern5_'+suffix,p1,p1Low,p1High)
  p2Var = RooRealVar('p2GaussBern5_'+suffix,'p2GaussBern5_'+suffix,p2,p2Low,p2High)
  p3Var = RooRealVar('p3GaussBern5_'+suffix,'p3GaussBern5_'+suffix,p3,p3Low,p3High)
  p4Var = RooRealVar('p4GaussBern5_'+suffix,'p4GaussBern5_'+suffix,p4,p4Low,p4High)
  p5Var = RooRealVar('p5GaussBern5_'+suffix,'p5GaussBern5_'+suffix,p5,p5Low,p5High)

  pArgs = RooArgList(p0Var,p1Var,p2Var,p3Var,p4Var,p5Var)
  GaussBern5 = RooGaussStepBernstein('GaussBern5_'+suffix,'GaussBern5_'+suffix,mzg,meanVar,sigmaVar,stepVar,pArgs)

  SetOwnership(meanVar,0)
  SetOwnership(sigmaVar,0)
  SetOwnership(stepVar,0)
  SetOwnership(p0Var,0)
  SetOwnership(p1Var,0)
  SetOwnership(p2Var,0)
  SetOwnership(p3Var,0)
  SetOwnership(p4Var,0)
  SetOwnership(p5Var,0)
  return GaussBern5

#def BuildGaussStepBern6(year,lepton,cat,mzg,mean = 0, sigma = 5, sigmaLow = 0.001, sigmaHigh = 60, step = 115, stepLow = 100, stepHigh = 130,
#    p0 = 15, p1 = 0.5, p1Low = -1e-6, p1High = 900,p2 = 0.5, p2Low = -1e-6, p2High = 900,p3 = 0.5, p3Low = -1e-6, p3High = 900,
#    p4 = 0.5, p4Low = -1e-6, p4High = 900, p5 = 0.5, p5Low = -1e-6, p5High = 900, p6 = 0.5, p6Low = -1e-6, p6High = 900):
def BuildGaussStepBern6(year,lepton,cat,mzg,mean = 0, sigma = 5, sigmaLow = 0.001, sigmaHigh = 60, step = 115, stepLow = 100, stepHigh = 130,
    p0 = 15, p1 = 0.5, p1Low = -1e2, p1High = 900,p2 = 0.5, p2Low = -1e2, p2High = 900,p3 = 0.5, p3Low = -1e2, p3High = 900,
    p4 = 0.5, p4Low = -1e2, p4High = 900, p5 = 0.5, p5Low = -1e2, p5High = 900, p6 = 0.5, p6Low = -1e2, p6High = 900):
  suffix = '_'.join([year,lepton,'cat'+cat])
  meanVar = RooRealVar('meanGaussBern6_'+suffix,'meanGaussBern6_'+suffix, mean)
  sigmaVar = RooRealVar('sigmaGaussBern6_'+suffix,'sigmaGaussBern6_'+suffix,sigma,sigmaLow,sigmaHigh)
  stepVar = RooRealVar('stepGaussBern6_'+suffix,'stepGaussBern6_'+suffix,step,stepLow,stepHigh)
  p0Var = RooRealVar('p0GaussBern6_'+suffix,'p0GaussBern6_'+suffix, p0)
  p1Var = RooRealVar('p1GaussBern6_'+suffix,'p1GaussBern6_'+suffix,p1,p1Low,p1High)
  p2Var = RooRealVar('p2GaussBern6_'+suffix,'p2GaussBern6_'+suffix,p2,p2Low,p2High)
  p3Var = RooRealVar('p3GaussBern6_'+suffix,'p3GaussBern6_'+suffix,p3,p3Low,p3High)
  p4Var = RooRealVar('p4GaussBern6_'+suffix,'p4GaussBern6_'+suffix,p4,p4Low,p4High)
  p5Var = RooRealVar('p5GaussBern6_'+suffix,'p5GaussBern6_'+suffix,p5,p5Low,p5High)
  p6Var = RooRealVar('p6GaussBern6_'+suffix,'p6GaussBern6_'+suffix,p6,p6Low,p6High)

  pArgs = RooArgList(p0Var,p1Var,p2Var,p3Var,p4Var,p5Var,p6Var)
  GaussBern6 = RooGaussStepBernstein('GaussBern6_'+suffix,'GaussBern6_'+suffix,mzg,meanVar,sigmaVar,stepVar,pArgs)

  SetOwnership(meanVar,0)
  SetOwnership(sigmaVar,0)
  SetOwnership(stepVar,0)
  SetOwnership(p0Var,0)
  SetOwnership(p1Var,0)
  SetOwnership(p2Var,0)
  SetOwnership(p3Var,0)
  SetOwnership(p4Var,0)
  SetOwnership(p5Var,0)
  SetOwnership(p6Var,0)
  return GaussBern6

def BuildSechStepBern3(year,lepton,cat,mzg,mean = 0, sigma = 3, sigmaLow = 0.01, sigmaHigh = 20, step = 0.1, stepLow = 0, stepHigh = 10,
    p0 = 15, p1 = 0.3, p1Low = -1e-6, p1High = 900,p2 = 0.3, p2Low = -1e-6, p2High = 900,p3 = 0.3, p3Low = -1e-6, p3High = 900):
  suffix = '_'.join([year,lepton,'cat'+cat])
  meanVar = RooRealVar('meanSechBern3_'+suffix,'meanSechBern3_'+suffix, mean)
  sigmaVar = RooRealVar('sigmaSechBern3_'+suffix,'sigmaSechBern3_'+suffix,sigma,sigmaLow,sigmaHigh)
  stepVar = RooRealVar('stepSechBern3_'+suffix,'stepSechBern3_'+suffix,step,stepLow,stepHigh)
  p0Var = RooRealVar('p0SechBern3_'+suffix,'p0SechBern3_'+suffix, p0)
  p1Var = RooRealVar('p1SechBern3_'+suffix,'p1SechBern3_'+suffix,p1,p1Low,p1High)
  p2Var = RooRealVar('p2SechBern3_'+suffix,'p2SechBern3_'+suffix,p2,p2Low,p2High)
  p3Var = RooRealVar('p3SechBern3_'+suffix,'p3SechBern3_'+suffix,p3,p3Low,p3High)

  pArgs = RooArgList(p0Var,p1Var,p2Var,p3Var)
  turnOn = RooGenericPdf('turnOnSechBern3_'+suffix,'turnOnSechBern3_'+suffix,'exp(-(@0-@1)/@2)/(@2*(1.0+exp(-(@0-@1)/@2))**2)',RooArgList(mzg,meanVar,sigmaVar))
  tail = RooStepBernstein('tailSechBern3_'+suffix,'tailSechBern3_'+suffix,mzg,stepVar,pArgs)
  SechBern3 = RooFFTConvPdf('SechBern3_'+suffix,'SechBern3_'+suffix,mzg,tail,turnOn)

  SetOwnership(meanVar,0)
  SetOwnership(sigmaVar,0)
  SetOwnership(stepVar,0)
  SetOwnership(p0Var,0)
  SetOwnership(p1Var,0)
  SetOwnership(p2Var,0)
  SetOwnership(p3Var,0)
  SetOwnership(turnOn,0)
  SetOwnership(tail,0)
  return SechBern3

def BuildSechStepBern4(year,lepton,cat,mzg,mean = 0, sigma = 3, sigmaLow = 0.01, sigmaHigh = 20, step = 0.1, stepLow = 0, stepHigh = 10,
    p0 = 15, p1 = 0.4, p1Low = -1e-6, p1High = 900,p2 = 0.4, p2Low = -1e-6, p2High = 900,p3 = 0.4, p3Low = -1e-6, p3High = 900, p4 = 0.4, p4Low = -1e-6, p4High = 900):
  suffix = '_'.join([year,lepton,'cat'+cat])
  meanVar = RooRealVar('meanSechBern4_'+suffix,'meanSechBern4_'+suffix, mean)
  sigmaVar = RooRealVar('sigmaSechBern4_'+suffix,'sigmaSechBern4_'+suffix,sigma,sigmaLow,sigmaHigh)
  stepVar = RooRealVar('stepSechBern4_'+suffix,'stepSechBern4_'+suffix,step,stepLow,stepHigh)
  p0Var = RooRealVar('p0SechBern4_'+suffix,'p0SechBern4_'+suffix, p0)
  p1Var = RooRealVar('p1SechBern4_'+suffix,'p1SechBern4_'+suffix,p1,p1Low,p1High)
  p2Var = RooRealVar('p2SechBern4_'+suffix,'p2SechBern4_'+suffix,p2,p2Low,p2High)
  p3Var = RooRealVar('p3SechBern4_'+suffix,'p3SechBern4_'+suffix,p3,p3Low,p3High)
  p4Var = RooRealVar('p4SechBern4_'+suffix,'p4SechBern4_'+suffix,p4,p4Low,p4High)

  pArgs = RooArgList(p0Var,p1Var,p2Var,p3Var,p4Var)
  turnOn = RooGenericPdf('turnOnSechBern4_'+suffix,'turnOnSechBern4_'+suffix,'exp(-(@0-@1)/@2)/(@2*(1.0+exp(-(@0-@1)/@2))**2)',RooArgList(mzg,meanVar,sigmaVar))
  tail = RooStepBernstein('tailSechBern4_'+suffix,'tailSechBern4_'+suffix,mzg,stepVar,pArgs)
  SechBern4 = RooFFTConvPdf('SechBern4_'+suffix,'SechBern4_'+suffix,mzg,tail,turnOn)

  SetOwnership(meanVar,0)
  SetOwnership(sigmaVar,0)
  SetOwnership(stepVar,0)
  SetOwnership(p0Var,0)
  SetOwnership(p1Var,0)
  SetOwnership(p2Var,0)
  SetOwnership(p3Var,0)
  SetOwnership(p4Var,0)
  SetOwnership(turnOn,0)
  SetOwnership(tail,0)
  return SechBern4

def BuildSechStepBern5(year,lepton,cat,mzg,mean = 0, sigma = 5, sigmaLow = 0.001, sigmaHigh = 60, step = 0.1, stepLow = 0, stepHigh = 10,
    p0 = 15, p1 = 0.5, p1Low = -1e-6, p1High = 900,p2 = 0.5, p2Low = -1e-6, p2High = 900,p3 = 0.5, p3Low = -1e-6, p3High = 900, p4 = 0.5, p4Low = -1e-6, p4High = 900, p5 = 0.5, p5Low = -1e-6, p5High = 900):
  suffix = '_'.join([year,lepton,'cat'+cat])
  meanVar = RooRealVar('meanSechBern5_'+suffix,'meanSechBern5_'+suffix, mean)
  sigmaVar = RooRealVar('sigmaSechBern5_'+suffix,'sigmaSechBern5_'+suffix,sigma,sigmaLow,sigmaHigh)
  stepVar = RooRealVar('stepSechBern5_'+suffix,'stepSechBern5_'+suffix,step,stepLow,stepHigh)
  p0Var = RooRealVar('p0SechBern5_'+suffix,'p0SechBern5_'+suffix, p0)
  p1Var = RooRealVar('p1SechBern5_'+suffix,'p1SechBern5_'+suffix,p1,p1Low,p1High)
  p2Var = RooRealVar('p2SechBern5_'+suffix,'p2SechBern5_'+suffix,p2,p2Low,p2High)
  p3Var = RooRealVar('p3SechBern5_'+suffix,'p3SechBern5_'+suffix,p3,p3Low,p3High)
  p4Var = RooRealVar('p4SechBern5_'+suffix,'p4SechBern5_'+suffix,p4,p4Low,p4High)
  p5Var = RooRealVar('p5SechBern5_'+suffix,'p5SechBern5_'+suffix,p5,p5Low,p5High)

  pArgs = RooArgList(p0Var,p1Var,p2Var,p3Var,p4Var,p5Var)
  turnOn = RooGenericPdf('turnOnSechBern5_'+suffix,'turnOnSechBern5_'+suffix,'exp(-(@0-@1)/@2)/(@2*(1.0+exp(-(@0-@1)/@2))**2)',RooArgList(mzg,meanVar,sigmaVar))
  tail = RooStepBernstein('tailSechBern5_'+suffix,'tailSechBern5_'+suffix,mzg,stepVar,pArgs)
  SechBern5 = RooFFTConvPdf('SechBern5_'+suffix,'SechBern5_'+suffix,mzg,tail,turnOn)

  SetOwnership(meanVar,0)
  SetOwnership(sigmaVar,0)
  SetOwnership(stepVar,0)
  SetOwnership(p0Var,0)
  SetOwnership(p1Var,0)
  SetOwnership(p2Var,0)
  SetOwnership(p3Var,0)
  SetOwnership(p4Var,0)
  SetOwnership(p5Var,0)
  SetOwnership(turnOn,0)
  SetOwnership(tail,0)
  return SechBern5

def BuildExp(year,lepton,cat,mzg,tau = 1, tauLow = -50, tauHigh = 50):
  suffix = '_'.join([year,lepton,'cat'+cat])
  tauVar = RooRealVar('tauExp_'+suffix,'tauExp_'+suffix,tau,tauLow,tauHigh)
  Exp = RooExponential('Exp_'+suffix,'Exp_'+suffix,mzg,tauVar)
  SetOwnership(tauVar,0)
  return Exp

def BuildPow(year,lepton,cat,mzg,alpha = 115, alphaLow = 50, alphaHigh = 200, beta = 3, betaLow = -20, betaHigh = 20):
  suffix = '_'.join([year,lepton,'cat'+cat])
  alphaVar = RooRealVar('alphaPow_'+suffix,'alphaPow_'+suffix,alpha,alphaLow,alphaHigh)
  betaVar = RooRealVar('betaPow_'+suffix,'betaPow_'+suffix,beta,betaLow,betaHigh)
  Pow = RooGenericPdf('Pow_'+suffix,'Pow_'+suffix,'1e-20 + (@1)*((@0)^(-@2))',RooArgList(mzg,alphaVar,betaVar))
  SetOwnership(alphaVar,0)
  SetOwnership(betaVar,0)
  return Pow

def BuildBern2(year,lepton,cat,mzg,p0 = 1 ,p1 = 5, p1Low = -1e-6, p1High = 30, p2 = 5, p2Low = -1e-6, p2High = 30):
  suffix = '_'.join([year,lepton,'cat'+cat])
  p0Var = RooRealVar('p0Bern2_'+suffix, 'p0Bern2_'+suffix,p0)
  p1Var = RooRealVar('p1Bern2_'+suffix, 'p1Bern2_'+suffix,p1,p1Low,p1High)
  p2Var = RooRealVar('p2Bern2_'+suffix, 'p2Bern2_'+suffix,p2,p2Low,p2High)
  Bern2 = RooBernstein('Bern2_'+suffix,'Bern2_'+suffix,mzg,RooArgList(p0Var,p1Var,p2Var))
  SetOwnership(p0Var,0)
  SetOwnership(p1Var,0)
  SetOwnership(p2Var,0)
  return Bern2

def BuildBern3(year,lepton,cat,mzg,p0 = 1 ,p1 = 5, p1Low = -1e-6, p1High = 30, p2 = 5, p2Low = -1e-6, p2High = 30, p3 = 5, p3Low = -1e-6, p3High = 30):
  suffix = '_'.join([year,lepton,'cat'+cat])
  p0Var = RooRealVar('p0Bern3_'+suffix, 'p0Bern3_'+suffix,p0)
  p1Var = RooRealVar('p1Bern3_'+suffix, 'p1Bern3_'+suffix,p1,p1Low,p1High)
  p2Var = RooRealVar('p2Bern3_'+suffix, 'p2Bern3_'+suffix,p2,p2Low,p2High)
  p3Var = RooRealVar('p3Bern3_'+suffix, 'p3Bern3_'+suffix,p3,p3Low,p3High)
  Bern3 = RooBernstein('Bern3_'+suffix,'Bern3_'+suffix,mzg,RooArgList(p0Var,p1Var,p2Var, p3Var))
  SetOwnership(p0Var,0)
  SetOwnership(p1Var,0)
  SetOwnership(p2Var,0)
  SetOwnership(p3Var,0)
  return Bern3

def BuildBern4(year,lepton,cat,mzg,p0 = 1 ,p1 = 5, p1Low = -1e-6, p1High = 30, p2 = 5, p2Low = -1e-6, p2High = 30, p3 = 5, p3Low = -1e-6, p3High = 30, p4 = 5, p4Low = -1e-6, p4High = 30):
  suffix = '_'.join([year,lepton,'cat'+cat])
  p0Var = RooRealVar('p0Bern4_'+suffix, 'p0Bern4_'+suffix,p0)
  p1Var = RooRealVar('p1Bern4_'+suffix, 'p1Bern4_'+suffix,p1,p1Low,p1High)
  p2Var = RooRealVar('p2Bern4_'+suffix, 'p2Bern4_'+suffix,p2,p2Low,p2High)
  p3Var = RooRealVar('p3Bern4_'+suffix, 'p3Bern4_'+suffix,p3,p3Low,p3High)
  p4Var = RooRealVar('p4Bern4_'+suffix, 'p4Bern4_'+suffix,p4,p4Low,p4High)
  Bern4 = RooBernstein('Bern4_'+suffix,'Bern4_'+suffix,mzg,RooArgList(p0Var,p1Var,p2Var, p3Var, p4Var))
  SetOwnership(p0Var,0)
  SetOwnership(p1Var,0)
  SetOwnership(p2Var,0)
  SetOwnership(p3Var,0)
  SetOwnership(p4Var,0)
  return Bern4

def BuildBern5(year,lepton,cat,mzg,p0 = 1 ,p1 = 5, p1Low = -1e-6, p1High = 30, p2 = 5, p2Low = -1e-6, p2High = 30, p3 = 5, p3Low = -1e-6, p3High = 30, p4 = 5, p4Low = -1e-6, p4High = 30, p5 = 5, p5Low = -1e-6, p5High = 30):
  suffix = '_'.join([year,lepton,'cat'+cat])
  p0Var = RooRealVar('p0Bern5_'+suffix, 'p0Bern5_'+suffix,p0)
  p1Var = RooRealVar('p1Bern5_'+suffix, 'p1Bern5_'+suffix,p1,p1Low,p1High)
  p2Var = RooRealVar('p2Bern5_'+suffix, 'p2Bern5_'+suffix,p2,p2Low,p2High)
  p3Var = RooRealVar('p3Bern5_'+suffix, 'p3Bern5_'+suffix,p3,p3Low,p3High)
  p4Var = RooRealVar('p4Bern5_'+suffix, 'p4Bern5_'+suffix,p4,p4Low,p4High)
  p5Var = RooRealVar('p5Bern5_'+suffix, 'p5Bern5_'+suffix,p5,p5Low,p5High)
  Bern5 = RooBernstein('Bern5_'+suffix,'Bern5_'+suffix,mzg,RooArgList(p0Var,p1Var,p2Var, p3Var, p4Var, p5Var))
  SetOwnership(p0Var,0)
  SetOwnership(p1Var,0)
  SetOwnership(p2Var,0)
  SetOwnership(p3Var,0)
  SetOwnership(p4Var,0)
  SetOwnership(p5Var,0)
  return Bern5


def BuildBern6(year,lepton,cat,mzg,p0 = 1 ,p1 = 5, p1Low = -1e-6, p1High = 30, p2 = 5, p2Low = -1e-6, p2High = 30, p3 = 5, p3Low = -1e-6, p3High = 30, p4 = 5, p4Low = -1e-6, p4High = 30, p5 = 5, p5Low = -1e-6, p5High = 30, p6 = 5, p6Low = -1e-6, p6High = 30):
  suffix = '_'.join([year,lepton,'cat'+cat])
  p0Var = RooRealVar('p0Bern6_'+suffix, 'p0Bern6_'+suffix,p0)
  p1Var = RooRealVar('p1Bern6_'+suffix, 'p1Bern6_'+suffix,p1,p1Low,p1High)
  p2Var = RooRealVar('p2Bern6_'+suffix, 'p2Bern6_'+suffix,p2,p2Low,p2High)
  p3Var = RooRealVar('p3Bern6_'+suffix, 'p3Bern6_'+suffix,p3,p3Low,p3High)
  p4Var = RooRealVar('p4Bern6_'+suffix, 'p4Bern6_'+suffix,p4,p4Low,p4High)
  p5Var = RooRealVar('p5Bern6_'+suffix, 'p5Bern6_'+suffix,p5,p5Low,p5High)
  p6Var = RooRealVar('p6Bern6_'+suffix, 'p6Bern6_'+suffix,p6,p6Low,p5High)
  Bern6 = RooBernstein('Bern6_'+suffix,'Bern6_'+suffix,mzg,RooArgList(p0Var,p1Var,p2Var, p3Var, p4Var, p5Var, p6Var))
  SetOwnership(p0Var,0)
  SetOwnership(p1Var,0)
  SetOwnership(p2Var,0)
  SetOwnership(p3Var,0)
  SetOwnership(p4Var,0)
  SetOwnership(p5Var,0)
  SetOwnership(p6Var,0)
  return Bern6


def BuildRooGaussian(year,lepton,cat,mzg, mean = 125,meanLow = 100, meanHigh = 150, sigma = 1.5, sigmaLow = 0.3, sigmaHigh = 70):
  suffix = '_'.join([year,lepton,'cat'+cat])
  meanVar = RooRealVar('mean_'+suffix,'mean_'+suffix, mean, meanLow, meanHigh)
  sigmaVar = RooRealVar('sigma_'+suffix,'sigma_'+suffix,sigma,sigmaLow,sigmaHigh)
  gauss = RooGaussian('gauss_'+suffix,'gauss_'+suffix,mzg,meanVar,sigmaVar)
  SetOwnership(meanVar,0)
  SetOwnership(sigmaVar,0)
  return gauss

def BuildCrystalBallGauss(year,lepton,cat,sig,mass,piece,mzg, mean = 125,meanLow = -1, meanHigh = -1, sigmaCB = 1.5, sigmaCBLow = 0.3, sigmaCBHigh = 20, alpha = 1, alphaLow = 0.5, alphaHigh = 10, n = 4, nLow = 0.5, nHigh = 50, sigmaG = 2, sigmaGLow = 0.3, sigmaGHigh = 20, frac = 0.1, fracLow = 0.0, fracHigh = 1.0):
  suffix = '_'.join([year,lepton,'cat'+cat,sig,mass,piece])
  if meanLow==-1:  meanLow  = mean-5
  if meanHigh==-1: meanHigh = mean+5
  meanVar    = RooRealVar('meanCBG_'   +suffix,'meanCBG_'   +suffix, mean, meanLow, meanHigh)
  sigmaCBVar = RooRealVar('sigmaCBCBG_'+suffix,'sigmaCBCBG_'+suffix, sigmaCB,sigmaCBLow,sigmaCBHigh)
  alphaVar   = RooRealVar('alphaCBG_'  +suffix,'alphaCBG_'  +suffix, alpha,alphaLow,alphaHigh)
  nVar       = RooRealVar('nCBG_'      +suffix,'nCBG_'      +suffix, n,nLow,nHigh)
  sigmaGVar  = RooRealVar('sigmaGCBG_' +suffix,'sigmaGCBG_' +suffix, sigmaG,sigmaGLow,sigmaGHigh)
  fracVar    = RooRealVar('fracCBG_'   +suffix,'fracCBG_'   +suffix, frac,fracLow,fracHigh)

  crystal = RooCBShape('crystalCBG_'+suffix,'crystalCBG_'+suffix,mzg,meanVar,sigmaCBVar,alphaVar,nVar)
  gauss = RooGaussian('gaussCBG_'+suffix,'gaussCBG_'+suffix,mzg,meanVar,sigmaGVar)
  cbArgs = RooArgList(gauss,crystal)
  fracArg = RooArgList(fracVar)
  CBG = RooAddPdf('CBG_'+suffix,'CBG_'+suffix,cbArgs,fracArg,True)

  SetOwnership(meanVar,0)
  SetOwnership(sigmaCBVar,0)
  SetOwnership(alphaVar,0)
  SetOwnership(nVar,0)
  SetOwnership(sigmaGVar,0)
  SetOwnership(fracVar,0)
  SetOwnership(crystal,0)
  SetOwnership(gauss,0)
  paramList = [meanVar,sigmaCBVar,alphaVar,nVar,sigmaGVar,fracVar]
  return CBG, paramList
  return CBG
