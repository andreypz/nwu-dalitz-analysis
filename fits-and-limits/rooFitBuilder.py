#!/usr/bin/env python
import sys
from ROOT import gSystem
gSystem.Load("libRooFit")
from ROOT import *
import numpy as np

gROOT.SetBatch()
gSystem.SetIncludePath( "-I$ROOFITSYS/include/" );
#gROOT.ProcessLine('.x RooStepBernstein.cxx+')
#gROOT.ProcessLine('.x RooGaussStepBernstein.cxx+')
#gROOT.ProcessLine('.x HZGRooPdfs.cxx++')

class FitBuilder:

  def __init__(self,mzg,tev,lepton,cat,sig=None,mass=None):
    if sig == None:
      self.suffix = '_'.join([tev,lepton,'cat'+cat])
    else:
      self.suffix = '_'.join([tev,lepton,'cat'+cat,sig,mass])
      self.sig = sig
      self.mass = mass

    self.mzg = mzg
    self.tev = tev
    self.lepton = lepton
    self.cat = cat
    self.BuildDict = {'BB':self.BuildBetaAndBern,
        'GB': self.BuildGaussAndBern,
        'BetaFunc': self.BuildBetaFunc,
        'Kumaraswamy': self.BuildKumaraswamy,
        'GaussExp': self.BuildGaussExp,
        'GaussPow': self.BuildGaussPow,
        'SechExp': self.BuildSechExp,
        'SechPow': self.BuildSechPow,
        'GaussBern3': self.BuildGaussStepBern3,
        'GaussBern4': self.BuildGaussStepBern4,
        'GaussBern5': self.BuildGaussStepBern5,
        'GaussBern6': self.BuildGaussStepBern6,
        'SechBern3': self.BuildSechStepBern3,
        'SechBern4': self.BuildSechStepBern4,
        'SechBern5': self.BuildSechStepBern5,
        'Exp': self.BuildExp,
        'Pow': self.BuildPow,
        'PowDecay': self.BuildPowDecay,
        'PowLog': self.BuildPowLog,
        'Exp2': self.BuildExp2,
        'ExpSum': self.BuildExpSum,
        'Laurent': self.BuildLaurent,
        'LaurentHZG': self.BuildLaurentHZG,
        'Bern2': self.BuildBern2,
        'Bern3': self.BuildBern3,
        'Bern4': self.BuildBern4,
        'Bern5': self.BuildBern5,
        'CBG': self.BuildCrystalBallGauss,
        'CBGM': self.BuildCrystalBallGaussM,
        'TripG': self.BuildTripleGauss}

    self.FitColorDict = {
      'GaussBern3':kViolet,
      'SechBern3':kMagenta,
      'GaussExp':kBlue,
      'GaussPow':kCyan,
      'SechExp':kRed,
      'SechPow':kYellow,
      'GaussBern4':kPink,
      'GaussBern5':kGray,
      'SechBern4':kBlack,
      'SechBern5':kGreen,
      'Exp':kBlue,
      'Exp2':kOrange,
      'ExpSum':kGreen+2,
      'Laurent':kMagenta,
      'LaurentHZG':kMagenta,
      'Pow':kCyan,
      'Bern2':kViolet,
      'Bern3':kPink,
      'Bern4':kGray,
      'Bern5':kGreen,
      'PowDecay':kYellow+1,
      'PowLog':kRed+1,
    }

    self.FitNdofDict = {
      'GaussBern3':5,
      'SechBern3':5,
      'GaussExp':3,
      'GaussPow':3,
      'SechExp':3,
      'SechPow':3,
      'GaussBern4':6,
      'GaussBern5':7,
      'SechBern4':6,
      'SechBern5':7,
      'Exp':1,
      'Exp2':2,
      'ExpSum':2,
      'Laurent':0,
      'LaurentHZG':1,
      'Pow':2,
      'Bern2':2,
      'Bern3':3,
      'Bern4':4,
      'Bern5':5,
      'PowDecay':2,
      'PowLog':3
    }



  def Build(self, funcName, **kargs):
    return self.BuildDict[funcName](**kargs)

  def BuildBetaAndBern(self,rangeName,frac = 0.1, fracLow = 0, fracHigh = 0.9):

    beta = BuildBetaFunc(tev,lepton,cat+'BB',self.mzg,rangeName)
    bern = BuildBern4(tev,lepton,cat+'BB',self.mzg)
    fracVar = RooRealVar('fracBB_'+self.suffix,'fracBB_'+self.suffix,frac,fracLow,fracHigh)
    bbArgs = RooArgList(beta,bern)
    fracArg = RooArgList(fracVar)

    BB = RooAddPdf('BB_'+self.suffix,'BB_'+self.suffix,bbArgs,fracArg,True)
    #BB = RooFFTConvPdf('BB_'+self.suffix,'BB_'+self.suffix,self.mzg,beta,bern)
    SetOwnership(beta,0)
    SetOwnership(bern,0)
    SetOwnership(fracVar,0)
    return BB

  def BuildGaussAndBern(self,rangeName,frac = 0.1, fracLow = 0, fracHigh = 0.9):

    gauss = BuildRooGaussian(tev,lepton,cat+'GB',self.mzg)
    bern = BuildBern3(tev,lepton,cat+'GB',self.mzg)
    fracVar = RooRealVar('fracGB_'+self.suffix,'fracGB_'+self.suffix,frac,fracLow,fracHigh)
    gbArgs = RooArgList(gauss,bern)
    fracArg = RooArgList(fracVar)

    GB = RooAddPdf('GB_'+self.suffix,'GB_'+self.suffix,gbArgs,fracArg,True)
    #BB = RooFFTConvPdf('BB_'+self.suffix,'BB_'+self.suffix,self.mzg,beta,bern)
    SetOwnership(gauss,0)
    SetOwnership(bern,0)
    SetOwnership(fracVar,0)
    return GB


  def BuildBetaFunc(self,rangeName,alpha = 2, alphaLow = 1, alphaHigh = 10, beta = 5, betaLow = 1, betaHigh = 10):

    alphaVar = RooRealVar('alphaBetaFunc_'+self.suffix,'alphaBetaFunc_'+self.suffix, alpha, alphaLow, alphaHigh)
    betaVar = RooRealVar('betaBetaFunc_'+self.suffix,'betaBetaFunc_'+self.suffix, beta, betaLow, betaHigh)
    xLow = self.mzg.getMin(rangeName)
    xRange = self.mzg.getMax(rangeName) - self.mzg.getMin(rangeName)
    gROOT.ProcessLine('.L betaWrapper.cxx+')
    from ROOT import makeBetaPdf
    BetaFunc = makeBetaPdf('BetaFunc_'+self.suffix,self.mzg,alphaVar,betaVar)

    SetOwnership(alphaVar,0)
    SetOwnership(betaVar,0)
    return BetaFunc

  def BuildKumaraswamy(self,rangeName,alpha = 2, alphaLow = 1, alphaHigh = 10, beta = 5, betaLow = 1, betaHigh = 10):

    alphaVar = RooRealVar('alphaKumaraswamy_'+self.suffix,'alphaKumaraswamy_'+self.suffix, alpha, alphaLow, alphaHigh)
    betaVar = RooRealVar('betaKumaraswamy_'+self.suffix,'betaKumaraswamy_'+self.suffix, beta, betaLow, betaHigh)
    xLow = RooRealVar('xLowKumaraswamy_'+self.suffix, 'xLowKumaraswamy_'+self.suffix, self.mzg.getMin(rangeName))
    xRange = RooRealVar('xRangeKumaraswamy_'+self.suffix, 'xReangeKumaraswamy_'+self.suffix,self.mzg.getMax(rangeName) - self.mzg.getMin(rangeName))
    Kumaraswamy = RooGenericPdf('Kumaraswamy_'+self.suffix, 'Kumaraswamy_'+self.suffix, '(@0>@3)*(@0<(@3+@4))*@1*@2*((@0-@3)/@4)**(@1-1)*(1-((@0-@3)/@4)**@1)**(@2)', RooArgList(self.mzg,alphaVar,betaVar,xLow,xRange))

    SetOwnership(alphaVar,0)
    SetOwnership(betaVar,0)
    SetOwnership(xLow,0)
    SetOwnership(xRange,0)
    return Kumaraswamy

  def BuildGaussExp(self,mean = 120, meanLow = 90, meanHigh = 150, sigma = 1, sigmaLow = 0.01, sigmaHigh = 10, tau = 5, tauLow = 0, tauHigh = 50):

    meanVar = RooRealVar('meanGaussExp_'+self.suffix,'meanGaussExp_'+self.suffix, mean, meanLow, meanHigh)
    sigmaVar = RooRealVar('sigmaGaussExp_'+self.suffix,'sigmaGaussExp_'+self.suffix,sigma,sigmaLow,sigmaHigh)
    tauVar = RooRealVar('tauGaussExp_'+self.suffix,'tauGaussExp_'+self.suffix,tau,tauLow,tauHigh)

    turnOn = RooGaussModel('turnOnGaussExp_'+self.suffix,'turnOnGaussExp_'+self.suffix,self.mzg,meanVar,sigmaVar)
    GaussExp = RooDecay('GaussExp_'+self.suffix,'GaussExp_'+self.suffix,self.mzg,tauVar,turnOn,RooDecay.SingleSided)

    SetOwnership(meanVar,0)
    SetOwnership(sigmaVar,0)
    SetOwnership(tauVar,0)
    SetOwnership(turnOn,0)
    return GaussExp
    #sprintf(sbuffer1, "GaussExpFullExt_cat%i",i+1);
    #GaussExpFullExt[i] = new RooExtendPdf(sbuffer1,sbuffer1,*GaussExpFull[i],*nGaussExp[i]);
    #nGaussExp[i] = new RooRealVar(sbuffer1,sbuffer1,ntcat[i]->GetEntries(),0,3*ntcat[i]->GetEntries());


  def BuildGaussPow(self,mean = 0, sigma = 2, sigmaLow = 0.01, sigmaHigh = 10, alpha = 105 , alphaLow = 50, alphaHigh = 200,beta = 6, betaLow = 0, betaHigh = 20):

    meanVar = RooRealVar('meanGaussPow_'+self.suffix,'meanGaussPow_'+self.suffix, mean)
    sigmaVar = RooRealVar('sigmaGaussPow_'+self.suffix,'sigmaGaussPow_'+self.suffix,sigma,sigmaLow,sigmaHigh)
    alphaVar = RooRealVar('alphaGaussPow_'+self.suffix,'alphaGaussPow_'+self.suffix,alpha,alphaLow,alphaHigh)
    betaVar = RooRealVar('betaGaussPow_'+self.suffix,'betaGaussPow_'+self.suffix,beta,betaLow,betaHigh)

    turnOn = RooGaussModel('turnOnGaussPow_'+self.suffix,'turnOnGaussPow_'+self.suffix,self.mzg,meanVar,sigmaVar)
    tail = RooGenericPdf('tailGaussPow_'+self.suffix,'tailGaussPow_'+self.suffix,'1e-20 + (@0 > @1)*((@0)^(-@2))',RooArgList(self.mzg,alphaVar,betaVar))
    GaussPow = RooFFTConvPdf('GaussPow_'+self.suffix,'GaussPow_'+self.suffix, self.mzg, tail, turnOn)

    SetOwnership(meanVar,0)
    SetOwnership(sigmaVar,0)
    SetOwnership(alphaVar,0)
    SetOwnership(betaVar,0)
    SetOwnership(turnOn,0)
    SetOwnership(tail,0)
    return GaussPow

  def BuildSechExp(self,mean = 0, sigma = 5, sigmaLow = 0.01, sigmaHigh = 20, tau = 35, tauLow = 0, tauHigh = 100, alpha = 105, alphaLow = 50, alphaHigh = 200):

    meanVar = RooRealVar('meanSechExp_'+self.suffix,'meanSechExp_'+self.suffix, mean)
    sigmaVar = RooRealVar('sigmaSechExp_'+self.suffix,'sigmaSechExp_'+self.suffix,sigma,sigmaLow,sigmaHigh)
    tauVar = RooRealVar('tauSechExp_'+self.suffix,'tauSechExp_'+self.suffix,tau,tauLow,tauHigh)
    alphaVar = RooRealVar('alphaSechExp_'+self.suffix,'alphaSechExp_'+self.suffix,alpha,alphaLow,alphaHigh)

    turnOn  = RooGenericPdf('turnOnSechExp_'+self.suffix,'turnOnSechExp_'+self.suffix, 'exp(-(@0-@1)/@2)/(@2*(1.0+exp(-(@0-@1)/@2))**2)',RooArgList(self.mzg,meanVar,sigmaVar))
    tail    = RooGenericPdf('tailSechExp_'+self.suffix,'tailSechExp_'+self.suffix,'1e-20 + (@0 > @1)*(exp(-@0/@2))',RooArgList(self.mzg,alphaVar,tauVar))
    SechExp = RooFFTConvPdf('SechExp_'+self.suffix,'SechExp_'+self.suffix,self.mzg,tail,turnOn)

    SetOwnership(meanVar,0)
    SetOwnership(sigmaVar,0)
    SetOwnership(alphaVar,0)
    SetOwnership(tauVar,0)
    SetOwnership(turnOn,0)
    SetOwnership(tail,0)
    return SechExp

  def BuildSechPow(self,mean = 0, sigma = 4, sigmaLow = 0.01, sigmaHigh = 20, alpha = 107, alphaLow = 50, alphaHigh = 200, beta = 5, betaLow = 0, betaHigh = 20):

    meanVar = RooRealVar('meanSechPow_'+self.suffix,'meanSechPow_'+self.suffix, mean)
    sigmaVar = RooRealVar('sigmaSechPow_'+self.suffix,'sigmaSechPow_'+self.suffix,sigma,sigmaLow,sigmaHigh)
    alphaVar = RooRealVar('alphaSechPow_'+self.suffix,'alphaSechPow_'+self.suffix,alpha,alphaLow,alphaHigh)
    betaVar = RooRealVar('betaSechPow_'+self.suffix,'betaSechPow_'+self.suffix,beta,betaLow,betaHigh)

    turnOn = RooGenericPdf('turnOnSechPow_'+self.suffix,'turnOnSechPow_'+self.suffix,'exp(-(@0-@1)/@2)/(@2*(1.0+exp(-(@0-@1)/@2))**2)',RooArgList(self.mzg,meanVar,sigmaVar))
    tail = RooGenericPdf('tailSechPow_'+self.suffix,'tailSechPow_'+self.suffix,'1e-20 + (@0 > @1)*((@0)^(-@2))',RooArgList(self.mzg,alphaVar,betaVar))
    SechPow = RooFFTConvPdf('SechPow_'+self.suffix,'SechPow_'+self.suffix, self.mzg, tail, turnOn)

    SetOwnership(meanVar,0)
    SetOwnership(sigmaVar,0)
    SetOwnership(alphaVar,0)
    SetOwnership(betaVar,0)
    SetOwnership(turnOn,0)
    SetOwnership(tail,0)
    return SechPow

  def BuildGaussStepBern3(self,mean = 0, sigma = 3, sigmaLow = 0.01, sigmaHigh = 20, step = 110, stepLow = 100, stepHigh = 130,
      p0 = 15, p1 = 0.3, p1Low = -1e-6, p1High = 900,p2 = 0.3, p2Low = -1e-6, p2High = 900,p3 = 0.3, p3Low = -1e-6, p3High = 900):

    meanVar = RooRealVar('meanGaussBern3_'+self.suffix,'meanGaussBern3_'+self.suffix, mean)
    sigmaVar = RooRealVar('sigmaGaussBern3_'+self.suffix,'sigmaGaussBern3_'+self.suffix,sigma,sigmaLow,sigmaHigh)
    stepVar = RooRealVar('stepGaussBern3_'+self.suffix,'stepGaussBern3_'+self.suffix,step,stepLow,stepHigh)
    p0Var = RooRealVar('p0GaussBern3_'+self.suffix,'p0GaussBern3_'+self.suffix, p0)
    p1Var = RooRealVar('p1GaussBern3_'+self.suffix,'p1GaussBern3_'+self.suffix,p1,p1Low,p1High)
    p2Var = RooRealVar('p2GaussBern3_'+self.suffix,'p2GaussBern3_'+self.suffix,p2,p2Low,p2High)
    p3Var = RooRealVar('p3GaussBern3_'+self.suffix,'p3GaussBern3_'+self.suffix,p3,p3Low,p3High)

    pArgs = RooArgList(p0Var,p1Var,p2Var,p3Var)
    GaussBern3 = RooGaussStepBernstein('GaussBern3_'+self.suffix,'GaussBern3_'+self.suffix,self.mzg,meanVar,sigmaVar,stepVar,pArgs)

    returnArgs = [meanVar,sigmaVar,stepVar,p0Var,p1Var,p2Var,p3Var]
    SetOwnership(meanVar,0)
    SetOwnership(sigmaVar,0)
    SetOwnership(stepVar,0)
    SetOwnership(p0Var,0)
    SetOwnership(p1Var,0)
    SetOwnership(p2Var,0)
    SetOwnership(p3Var,0)
    return GaussBern3,returnArgs

  def BuildGaussStepBern3Const(self,mean = 0, sigma = 4, step = 115,
      p0 = 15, p1 = 0.3, p2 = 0.3, p3 = 0.3):

    meanVar = RooRealVar('meanGaussBern3_'+self.suffix,'meanGaussBern3_'+self.suffix, mean)
    sigmaVar = RooRealVar('sigmaGaussBern3_'+self.suffix,'sigmaGaussBern3_'+self.suffix,sigma)
    stepVar = RooRealVar('stepGaussBern3_'+self.suffix,'stepGaussBern3_'+self.suffix,step)
    p0Var = RooRealVar('p0GaussBern3_'+self.suffix,'p0GaussBern3_'+self.suffix,p0)
    p1Var = RooRealVar('p1GaussBern3_'+self.suffix,'p1GaussBern3_'+self.suffix,p1)
    p2Var = RooRealVar('p2GaussBern3_'+self.suffix,'p2GaussBern3_'+self.suffix,p2)
    p3Var = RooRealVar('p3GaussBern3_'+self.suffix,'p3GaussBern3_'+self.suffix,p3)

    pArgs = RooArgList(p0Var,p1Var,p2Var,p3Var)
    GaussBern3 = RooGaussStepBernstein('GaussBern3_'+self.suffix,'GaussBern3_'+self.suffix,self.mzg,meanVar,sigmaVar,stepVar,pArgs)

    returnArgs = [meanVar,sigmaVar,stepVar,p0Var,p1Var,p2Var,p3Var]
    SetOwnership(meanVar,0)
    SetOwnership(sigmaVar,0)
    SetOwnership(stepVar,0)
    SetOwnership(p0Var,0)
    SetOwnership(p1Var,0)
    SetOwnership(p2Var,0)
    SetOwnership(p3Var,0)
    return GaussBern3,returnArgs

  def BuildGaussStepBern4(self,mean = 0, sigma = 4, sigmaLow = 0.001, sigmaHigh = 60, step = 115, stepLow = 100, stepHigh = 130,
      p0 = 15, p1 = 0.4, p1Low = -1e-6, p1High = 900,p2 = 0.4, p2Low = -1e-6, p2High = 900,p3 = 0.4, p3Low = -1e-6, p3High = 900, p4 = 0.4, p4Low = -1e-6, p4High = 900):
    #def BuildGaussStepBern4(self,mean = 0, sigma = 4, sigmaLow = 0.001, sigmaHigh = 60, step = 115, stepLow = 100, stepHigh = 130,
      #p0 = 15, p1 = 0.4, p1Low = -1e2, p1High = 900,p2 = 0.4, p2Low = -1e2, p2High = 900,p3 = 0.4, p3Low = -1e2, p3High = 900, p4 = 0.4, p4Low = -1e2, p4High = 900):

    meanVar = RooRealVar('meanGaussBern4_'+self.suffix,'meanGaussBern4_'+self.suffix, mean)
    sigmaVar = RooRealVar('sigmaGaussBern4_'+self.suffix,'sigmaGaussBern4_'+self.suffix,sigma,sigmaLow,sigmaHigh)
    stepVar = RooRealVar('stepGaussBern4_'+self.suffix,'stepGaussBern4_'+self.suffix,step,stepLow,stepHigh)
    p0Var = RooRealVar('p0GaussBern4_'+self.suffix,'p0GaussBern4_'+self.suffix, p0)
    p1Var = RooRealVar('p1GaussBern4_'+self.suffix,'p1GaussBern4_'+self.suffix,p1,p1Low,p1High)
    p2Var = RooRealVar('p2GaussBern4_'+self.suffix,'p2GaussBern4_'+self.suffix,p2,p2Low,p2High)
    p3Var = RooRealVar('p3GaussBern4_'+self.suffix,'p3GaussBern4_'+self.suffix,p3,p3Low,p3High)
    p4Var = RooRealVar('p4GaussBern4_'+self.suffix,'p4GaussBern4_'+self.suffix,p4,p4Low,p4High)

    pArgs = RooArgList(p0Var,p1Var,p2Var,p3Var,p4Var)
    GaussBern4 = RooGaussStepBernstein('GaussBern4_'+self.suffix,'GaussBern4_'+self.suffix,self.mzg,meanVar,sigmaVar,stepVar,pArgs)

    SetOwnership(meanVar,0)
    SetOwnership(sigmaVar,0)
    SetOwnership(stepVar,0)
    SetOwnership(p0Var,0)
    SetOwnership(p1Var,0)
    SetOwnership(p2Var,0)
    SetOwnership(p3Var,0)
    SetOwnership(p4Var,0)
    return GaussBern4

  def BuildGaussStepBern5(self,mean = 0, sigma = 5, sigmaLow = 0.001, sigmaHigh = 60, step = 115, stepLow = 100, stepHigh = 130,
      p0 = 15, p1 = 0.5, p1Low = -1e-6, p1High = 900,p2 = 0.5, p2Low = -1e-6, p2High = 900,p3 = 0.5, p3Low = -1e-6, p3High = 900, p4 = 0.5, p4Low = -1e-6, p4High = 900, p5 = 0.5, p5Low = -1e-6, p5High = 900):
    #def BuildGaussStepBern5(self,mean = 0, sigma = 5, sigmaLow = 0.001, sigmaHigh = 60, step = 115, stepLow = 100, stepHigh = 130,
     # p0 = 15, p1 = 0.5, p1Low = -1e2, p1High = 900,p2 = 0.5, p2Low = -1e2, p2High = 900,p3 = 0.5, p3Low = -1e2, p3High = 900, p4 = 0.5, p4Low = -1e2, p4High = 900, p5 = 0.5, p5Low = -1e2, p5High = 900):

    meanVar = RooRealVar('meanGaussBern5_'+self.suffix,'meanGaussBern5_'+self.suffix, mean)
    sigmaVar = RooRealVar('sigmaGaussBern5_'+self.suffix,'sigmaGaussBern5_'+self.suffix,sigma,sigmaLow,sigmaHigh)
    stepVar = RooRealVar('stepGaussBern5_'+self.suffix,'stepGaussBern5_'+self.suffix,step,stepLow,stepHigh)
    p0Var = RooRealVar('p0GaussBern5_'+self.suffix,'p0GaussBern5_'+self.suffix, p0)
    p1Var = RooRealVar('p1GaussBern5_'+self.suffix,'p1GaussBern5_'+self.suffix,p1,p1Low,p1High)
    p2Var = RooRealVar('p2GaussBern5_'+self.suffix,'p2GaussBern5_'+self.suffix,p2,p2Low,p2High)
    p3Var = RooRealVar('p3GaussBern5_'+self.suffix,'p3GaussBern5_'+self.suffix,p3,p3Low,p3High)
    p4Var = RooRealVar('p4GaussBern5_'+self.suffix,'p4GaussBern5_'+self.suffix,p4,p4Low,p4High)
    p5Var = RooRealVar('p5GaussBern5_'+self.suffix,'p5GaussBern5_'+self.suffix,p5,p5Low,p5High)

    pArgs = RooArgList(p0Var,p1Var,p2Var,p3Var,p4Var,p5Var)
    GaussBern5 = RooGaussStepBernstein('GaussBern5_'+self.suffix,'GaussBern5_'+self.suffix,self.mzg,meanVar,sigmaVar,stepVar,pArgs)

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

#def BuildGaussStepBern6(self,mean = 0, sigma = 5, sigmaLow = 0.001, sigmaHigh = 60, step = 115, stepLow = 100, stepHigh = 130,
#    p0 = 15, p1 = 0.5, p1Low = -1e-6, p1High = 900,p2 = 0.5, p2Low = -1e-6, p2High = 900,p3 = 0.5, p3Low = -1e-6, p3High = 900,
#    p4 = 0.5, p4Low = -1e-6, p4High = 900, p5 = 0.5, p5Low = -1e-6, p5High = 900, p6 = 0.5, p6Low = -1e-6, p6High = 900):
  def BuildGaussStepBern6(self,mean = 0, sigma = 5, sigmaLow = 0.001, sigmaHigh = 60, step = 115, stepLow = 100, stepHigh = 130,
      p0 = 15, p1 = 0.5, p1Low = -1e2, p1High = 900,p2 = 0.5, p2Low = -1e2, p2High = 900,p3 = 0.5, p3Low = -1e2, p3High = 900,
      p4 = 0.5, p4Low = -1e2, p4High = 900, p5 = 0.5, p5Low = -1e2, p5High = 900, p6 = 0.5, p6Low = -1e2, p6High = 900):

    meanVar = RooRealVar('meanGaussBern6_'+self.suffix,'meanGaussBern6_'+self.suffix, mean)
    sigmaVar = RooRealVar('sigmaGaussBern6_'+self.suffix,'sigmaGaussBern6_'+self.suffix,sigma,sigmaLow,sigmaHigh)
    stepVar = RooRealVar('stepGaussBern6_'+self.suffix,'stepGaussBern6_'+self.suffix,step,stepLow,stepHigh)
    p0Var = RooRealVar('p0GaussBern6_'+self.suffix,'p0GaussBern6_'+self.suffix, p0)
    p1Var = RooRealVar('p1GaussBern6_'+self.suffix,'p1GaussBern6_'+self.suffix,p1,p1Low,p1High)
    p2Var = RooRealVar('p2GaussBern6_'+self.suffix,'p2GaussBern6_'+self.suffix,p2,p2Low,p2High)
    p3Var = RooRealVar('p3GaussBern6_'+self.suffix,'p3GaussBern6_'+self.suffix,p3,p3Low,p3High)
    p4Var = RooRealVar('p4GaussBern6_'+self.suffix,'p4GaussBern6_'+self.suffix,p4,p4Low,p4High)
    p5Var = RooRealVar('p5GaussBern6_'+self.suffix,'p5GaussBern6_'+self.suffix,p5,p5Low,p5High)
    p6Var = RooRealVar('p6GaussBern6_'+self.suffix,'p6GaussBern6_'+self.suffix,p6,p6Low,p6High)

    pArgs = RooArgList(p0Var,p1Var,p2Var,p3Var,p4Var,p5Var,p6Var)
    GaussBern6 = RooGaussStepBernstein('GaussBern6_'+self.suffix,'GaussBern6_'+self.suffix,self.mzg,meanVar,sigmaVar,stepVar,pArgs)

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

  def BuildSechStepBern3(self,mean = 0, sigma = 10, sigmaLow = 0.01, sigmaHigh = 20, step = 0.1, stepLow = 0, stepHigh = 10,
      p0 = 15, p1 = 0.3, p1Low = -1e-6, p1High = 900,p2 = 0.3, p2Low = -1e-6, p2High = 900,p3 = 0.3, p3Low = -1e-6, p3High = 900):

    meanVar = RooRealVar('meanSechBern3_'+self.suffix,'meanSechBern3_'+self.suffix, mean)
    sigmaVar = RooRealVar('sigmaSechBern3_'+self.suffix,'sigmaSechBern3_'+self.suffix,sigma,sigmaLow,sigmaHigh)
    stepVar = RooRealVar('stepSechBern3_'+self.suffix,'stepSechBern3_'+self.suffix,step,stepLow,stepHigh)
    p0Var = RooRealVar('p0SechBern3_'+self.suffix,'p0SechBern3_'+self.suffix, p0)
    p1Var = RooRealVar('p1SechBern3_'+self.suffix,'p1SechBern3_'+self.suffix,p1,p1Low,p1High)
    p2Var = RooRealVar('p2SechBern3_'+self.suffix,'p2SechBern3_'+self.suffix,p2,p2Low,p2High)
    p3Var = RooRealVar('p3SechBern3_'+self.suffix,'p3SechBern3_'+self.suffix,p3,p3Low,p3High)

    pArgs = RooArgList(p0Var,p1Var,p2Var,p3Var)
    turnOn = RooGenericPdf('turnOnSechBern3_'+self.suffix,'turnOnSechBern3_'+self.suffix,'exp(-(@0-@1)/@2)/(@2*(1.0+exp(-(@0-@1)/@2))**2)',RooArgList(self.mzg,meanVar,sigmaVar))
    tail = RooStepBernstein('tailSechBern3_'+self.suffix,'tailSechBern3_'+self.suffix,self.mzg,stepVar,pArgs)
    SechBern3 = RooFFTConvPdf('SechBern3_'+self.suffix,'SechBern3_'+self.suffix,self.mzg,tail,turnOn)

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

  def BuildSechStepBern4(self,mean = 0, sigma = 3, sigmaLow = 0.01, sigmaHigh = 20, step = 0.1, stepLow = 0, stepHigh = 10,
      p0 = 15, p1 = 0.4, p1Low = -1e-6, p1High = 900,p2 = 0.4, p2Low = -1e-6, p2High = 900,p3 = 0.4, p3Low = -1e-6, p3High = 900, p4 = 0.4, p4Low = -1e-6, p4High = 900):

    meanVar = RooRealVar('meanSechBern4_'+self.suffix,'meanSechBern4_'+self.suffix, mean)
    sigmaVar = RooRealVar('sigmaSechBern4_'+self.suffix,'sigmaSechBern4_'+self.suffix,sigma,sigmaLow,sigmaHigh)
    stepVar = RooRealVar('stepSechBern4_'+self.suffix,'stepSechBern4_'+self.suffix,step,stepLow,stepHigh)
    p0Var = RooRealVar('p0SechBern4_'+self.suffix,'p0SechBern4_'+self.suffix, p0)
    p1Var = RooRealVar('p1SechBern4_'+self.suffix,'p1SechBern4_'+self.suffix,p1,p1Low,p1High)
    p2Var = RooRealVar('p2SechBern4_'+self.suffix,'p2SechBern4_'+self.suffix,p2,p2Low,p2High)
    p3Var = RooRealVar('p3SechBern4_'+self.suffix,'p3SechBern4_'+self.suffix,p3,p3Low,p3High)
    p4Var = RooRealVar('p4SechBern4_'+self.suffix,'p4SechBern4_'+self.suffix,p4,p4Low,p4High)

    pArgs = RooArgList(p0Var,p1Var,p2Var,p3Var,p4Var)
    turnOn = RooGenericPdf('turnOnSechBern4_'+self.suffix,'turnOnSechBern4_'+self.suffix,'exp(-(@0-@1)/@2)/(@2*(1.0+exp(-(@0-@1)/@2))**2)',RooArgList(self.mzg,meanVar,sigmaVar))
    tail = RooStepBernstein('tailSechBern4_'+self.suffix,'tailSechBern4_'+self.suffix,self.mzg,stepVar,pArgs)
    SechBern4 = RooFFTConvPdf('SechBern4_'+self.suffix,'SechBern4_'+self.suffix,self.mzg,tail,turnOn)

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

  def BuildSechStepBern5(self,mean = 0, sigma = 5, sigmaLow = 0.001, sigmaHigh = 60, step = 0.1, stepLow = 0, stepHigh = 10,
      p0 = 15, p1 = 0.5, p1Low = -1e-6, p1High = 900,p2 = 0.5, p2Low = -1e-6, p2High = 900,p3 = 0.5, p3Low = -1e-6, p3High = 900, p4 = 0.5, p4Low = -1e-6, p4High = 900, p5 = 0.5, p5Low = -1e-6, p5High = 900):

    meanVar = RooRealVar('meanSechBern5_'+self.suffix,'meanSechBern5_'+self.suffix, mean)
    sigmaVar = RooRealVar('sigmaSechBern5_'+self.suffix,'sigmaSechBern5_'+self.suffix,sigma,sigmaLow,sigmaHigh)
    stepVar = RooRealVar('stepSechBern5_'+self.suffix,'stepSechBern5_'+self.suffix,step,stepLow,stepHigh)
    p0Var = RooRealVar('p0SechBern5_'+self.suffix,'p0SechBern5_'+self.suffix, p0)
    p1Var = RooRealVar('p1SechBern5_'+self.suffix,'p1SechBern5_'+self.suffix,p1,p1Low,p1High)
    p2Var = RooRealVar('p2SechBern5_'+self.suffix,'p2SechBern5_'+self.suffix,p2,p2Low,p2High)
    p3Var = RooRealVar('p3SechBern5_'+self.suffix,'p3SechBern5_'+self.suffix,p3,p3Low,p3High)
    p4Var = RooRealVar('p4SechBern5_'+self.suffix,'p4SechBern5_'+self.suffix,p4,p4Low,p4High)
    p5Var = RooRealVar('p5SechBern5_'+self.suffix,'p5SechBern5_'+self.suffix,p5,p5Low,p5High)

    pArgs = RooArgList(p0Var,p1Var,p2Var,p3Var,p4Var,p5Var)
    turnOn = RooGenericPdf('turnOnSechBern5_'+self.suffix,'turnOnSechBern5_'+self.suffix,'exp(-(@0-@1)/@2)/(@2*(1.0+exp(-(@0-@1)/@2))**2)',RooArgList(self.mzg,meanVar,sigmaVar))
    tail = RooStepBernstein('tailSechBern5_'+self.suffix,'tailSechBern5_'+self.suffix,self.mzg,stepVar,pArgs)
    SechBern5 = RooFFTConvPdf('SechBern5_'+self.suffix,'SechBern5_'+self.suffix,self.mzg,tail,turnOn)

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

  def BuildExp(self,tau = 1, tauLow = -50, tauHigh = 50):

    tauVar = RooRealVar('tauExp_'+self.suffix,'tauExp_'+self.suffix,tau,tauLow,tauHigh)
    Exp = RooExponential('Exp_'+self.suffix,'Exp_'+self.suffix,self.mzg,tauVar)
    SetOwnership(tauVar,0)
    return Exp

  def BuildPow(self,alpha = 115, alphaLow = -20, alphaHigh = 200, beta = 3, betaLow = -20, betaHigh = 20):

    alphaVar = RooRealVar('alphaPow_'+self.suffix,'alphaPow_'+self.suffix,alpha,alphaLow,alphaHigh)
    betaVar = RooRealVar('betaPow_'+self.suffix,'betaPow_'+self.suffix,beta,betaLow,betaHigh)
    Pow = RooGenericPdf('Pow_'+self.suffix,'Pow_'+self.suffix,'1e-20 + (@1)*((@0)^(-@2))',RooArgList(self.mzg,alphaVar,betaVar))
    SetOwnership(alphaVar,0)
    SetOwnership(betaVar,0)
    return Pow

  def BuildPowDecay(self,p1 = 1, p1Low = -20, p1High = 20, p2 = 1, p2Low = -20, p2High = 20):

    p1Var = RooRealVar('p1PowDecay_'+self.suffix,'p1PowDecay_'+self.suffix,p1,p1Low,p1High)
    p2Var = RooRealVar('p2PowDecay_'+self.suffix,'p2PowDecay_'+self.suffix,p2,p2Low,p2High)
    PowDecay = RooGenericPdf('PowDecay_'+self.suffix,'PowDecay_'+self.suffix,'exp(-@1)*((@0)^(-@2))',RooArgList(self.mzg,p1Var,p2Var))
    SetOwnership(p1Var,0)
    SetOwnership(p2Var,0)
    return PowDecay

  def BuildPowLog(self,p1 = 2, p1Low = -20, p1High = 20, p2 = -1, p2Low = -20, p2High = 20, p3 = 0.5, p3Low = -20, p3High = 20):

    p1Var = RooRealVar('p1PowLog_'+self.suffix,'p1PowLog_'+self.suffix,p1,p1Low,p1High)
    p2Var = RooRealVar('p2PowLog_'+self.suffix,'p2PowLog_'+self.suffix,p2,p2Low,p2High)
    p3Var = RooRealVar('p3PowLog_'+self.suffix,'p3PowLog_'+self.suffix,p3,p3Low,p3High)
    PowLog = RooGenericPdf('PowLog_'+self.suffix,'PowLog_'+self.suffix,'(1-@0)^(@1)/(@0)^(@2+@3*log(@0))',RooArgList(self.mzg,p1Var,p2Var,p3Var))
    SetOwnership(p1Var,0)
    SetOwnership(p2Var,0)
    SetOwnership(p3Var,0)
    return PowLog

  def BuildExp2(self,p1 = 1, p1Low = -20, p1High = 20, p2 = 1, p2Low = -20, p2High = 20):

    p1Var = RooRealVar('p1Exp2_'+self.suffix,'p1Exp2_'+self.suffix,p1,p1Low,p1High)
    p2Var = RooRealVar('p2Exp2_'+self.suffix,'p2Exp2_'+self.suffix,p2,p2Low,p2High)
    Exp2 = RooGenericPdf('Exp2_'+self.suffix,'Exp2_'+self.suffix,'exp((-@0)/(@1+@2*@0))',RooArgList(self.mzg,p1Var,p2Var))
    SetOwnership(p1Var,0)
    SetOwnership(p2Var,0)
    return Exp2

  def BuildExpSum(self,p1 = 1.3, p1Low = -20, p1High = 20, p2 = 0.1, p2Low = -20, p2High = 20, p3 = 0.01, p3Low = -20, p3High = 20):

    p1Var = RooRealVar('p1ExpSum_'+self.suffix,'p1ExpSum_'+self.suffix,p1,p1Low,p1High)
    p2Var = RooRealVar('p2ExpSum_'+self.suffix,'p2ExpSum_'+self.suffix,p2,p2Low,p2High)
    p3Var = RooRealVar('p3ExpSum_'+self.suffix,'p3ExpSum_'+self.suffix,p3,p3Low,p3High)
    ExpSum = RooGenericPdf('ExpSum_'+self.suffix,'ExpSum_'+self.suffix,'(1-@1)*exp(-@2*@0)+@1*exp(-@3*@0)',RooArgList(self.mzg,p1Var,p2Var,p3Var))
    SetOwnership(p1Var,0)
    SetOwnership(p2Var,0)
    SetOwnership(p3Var,0)
    return ExpSum

  def BuildLaurentHZG(self,p1 = 1, p1Low = -20, p1High = 20):

    p1Var   = RooRealVar('p1Laurent_'+ self.suffix,'p1Laurent_'+self.suffix,p1,p1Low,p1High)
    Laurent = RooGenericPdf('Laurent_'+self.suffix,'Laurent_'+self.suffix,'(1-@1)*@0^(-4)+@1*@0^(-5)',RooArgList(self.mzg,p1Var))
    SetOwnership(p1Var,0)
    return Laurent

  def BuildLaurent(self, order = 1):
    # This is an implementation similar to h2gglobe.
    # Zeroth order is power -4

    if order>=1:
      p1Var   = RooRealVar('p1Laurent_'+ self.suffix,'p1Laurent_'+self.suffix, 0.25/order, 0.00001,0.99999)
      Laurent = RooGenericPdf('Laurent_'+self.suffix,'Laurent_'+self.suffix,'(1-@1)*@0^(-4)+@1*@0^(-5)',RooArgList(self.mzg,p1Var))

      SetOwnership(p1Var,0)
    if order==2:
      p2Var   = RooRealVar('p2Laurent_'+ self.suffix,'p2Laurent_'+self.suffix, 0.25/order, 0.00001,0.99999)
      Laurent = RooGenericPdf('Laurent_'+self.suffix,'Laurent_'+self.suffix,'(1-@1-@2)*@0^(-3)+@1*@0^(-4)+@2*@0^(-5)',RooArgList(self.mzg,p1Var,p2Var))
      SetOwnership(p2Var,0)
    if order>2:
      print '\n \t\t Error only orders 1 and 2 are supported for Laurent polynomials! \n'
      sys.exit(0)

    self.FitNdofDict['Laurent'] = order
    return Laurent

  def BuildBern2(self, p1 = 0.1, p1Low = 1e-3, p1High = 1, p2 = 0.1, p2Low = 1e-3, p2High = 1):
    #p0Var = RooRealVar('p0Bern2_'+self.suffix, 'p0Bern2_'+self.suffix,p0)
    p1Var = RooRealVar('p1Bern2_'+self.suffix, 'p1Bern2_'+self.suffix,p1,p1Low,p1High)
    p2Var = RooRealVar('p2Bern2_'+self.suffix, 'p2Bern2_'+self.suffix,p2,p2Low,p2High)
    Bern2 = RooBernstein('Bern2_'+self.suffix,'Bern2_'+self.suffix,
                         self.mzg,RooArgList(RooFit.RooConst(1.0),p1Var,p2Var))
    #SetOwnership(p0Var,0)
    SetOwnership(p1Var,0)
    SetOwnership(p2Var,0)
    return Bern2

  def BuildBern3(self, p1 = 0.1, p1Low =1e-3, p1High = 1, p2 = 0.1, p2Low =1e-3, p2High = 1,
                 p3 = 0.1, p3Low =1e-3, p3High = 1):
    #p0Var = RooRealVar('p0Bern3_'+self.suffix, 'p0Bern3_'+self.suffix,p0)
    p1Var = RooRealVar('p1Bern3_'+self.suffix, 'p1Bern3_'+self.suffix,p1,p1Low,p1High)
    p2Var = RooRealVar('p2Bern3_'+self.suffix, 'p2Bern3_'+self.suffix,p2,p2Low,p2High)
    p3Var = RooRealVar('p3Bern3_'+self.suffix, 'p3Bern3_'+self.suffix,p3,p3Low,p3High)
    Bern3 = RooBernstein('Bern3_'+self.suffix,'Bern3_'+self.suffix,
                         self.mzg,RooArgList(RooFit.RooConst(1.0),p1Var,p2Var, p3Var))
    #SetOwnership(p0Var,0)
    SetOwnership(p1Var,0)
    SetOwnership(p2Var,0)
    SetOwnership(p3Var,0)
    return Bern3

  def BuildBern4(self, p1 = 0.1, p1Low = 1e-3, p1High = 1, p2 = 0.1, p2Low = 1e-3, p2High = 1,
                 p3 = 0.1, p3Low = 1e-3, p3High = 1, p4 = 0.1, p4Low = 1e-3, p4High = 1):
    #p0Var = RooRealVar('p0Bern4_'+self.suffix, 'p0Bern4_'+self.suffix,p0)
    p1Var = RooRealVar('p1Bern4_'+self.suffix, 'p1Bern4_'+self.suffix,p1,p1Low,p1High)
    p2Var = RooRealVar('p2Bern4_'+self.suffix, 'p2Bern4_'+self.suffix,p2,p2Low,p2High)
    p3Var = RooRealVar('p3Bern4_'+self.suffix, 'p3Bern4_'+self.suffix,p3,p3Low,p3High)
    p4Var = RooRealVar('p4Bern4_'+self.suffix, 'p4Bern4_'+self.suffix,p4,p4Low,p4High)
    Bern4 = RooBernstein('Bern4_'+self.suffix,'Bern4_'+self.suffix,
                         self.mzg,RooArgList(RooFit.RooConst(1.0),p1Var,p2Var, p3Var, p4Var))
    #SetOwnership(p0Var,0)
    SetOwnership(p1Var,0)
    SetOwnership(p2Var,0)
    SetOwnership(p3Var,0)
    SetOwnership(p4Var,0)
    return Bern4

  def BuildBern5(self,p0 = 1 ,p1 = 0.1, p1Low = 1e-3, p1High = 1, p2 = 0.1, p2Low = 1e-3, p2High = 1,
                 p3 = 0.1, p3Low = 1e-3, p3High = 1, p4 = 0.1, p4Low = 1e-3, p4High = 1, p5 = 0.1, p5Low = 1e-3, p5High = 1):

    p0Var = RooRealVar('p0Bern5_'+self.suffix, 'p0Bern5_'+self.suffix,p0)
    p1Var = RooRealVar('p1Bern5_'+self.suffix, 'p1Bern5_'+self.suffix,p1,p1Low,p1High)
    p2Var = RooRealVar('p2Bern5_'+self.suffix, 'p2Bern5_'+self.suffix,p2,p2Low,p2High)
    p3Var = RooRealVar('p3Bern5_'+self.suffix, 'p3Bern5_'+self.suffix,p3,p3Low,p3High)
    p4Var = RooRealVar('p4Bern5_'+self.suffix, 'p4Bern5_'+self.suffix,p4,p4Low,p4High)
    p5Var = RooRealVar('p5Bern5_'+self.suffix, 'p5Bern5_'+self.suffix,p5,p5Low,p5High)
    Bern5 = RooBernstein('Bern5_'+self.suffix,'Bern5_'+self.suffix,self.mzg,RooArgList(p0Var,p1Var,p2Var, p3Var, p4Var, p5Var))
    SetOwnership(p0Var,0)
    SetOwnership(p1Var,0)
    SetOwnership(p2Var,0)
    SetOwnership(p3Var,0)
    SetOwnership(p4Var,0)
    SetOwnership(p5Var,0)
    return Bern5

  def BuildPoly3(self,p0 = 1 ,p1 = 5, p1Low = -3000, p1High = 3000, p2 = 5, p2Low = -3000, p2High = 3000, p3 = 5, p3Low = -3000, p3High = 3000):
    p0Var = RooRealVar('p0Poly3_'+self.suffix, 'p0Poly3_'+self.suffix,p0)
    p1Var = RooRealVar('p1Poly3_'+self.suffix, 'p1Poly3_'+self.suffix,p1,p1Low,p1High)
    p2Var = RooRealVar('p2Poly3_'+self.suffix, 'p2Poly3_'+self.suffix,p2,p2Low,p2High)
    p3Var = RooRealVar('p3Poly3_'+self.suffix, 'p3Poly3_'+self.suffix,p3,p3Low,p3High)
    Poly3 = RooPolynomial('Poly3_'+self.suffix,'Poly3_'+self.suffix,self.mzg,RooArgList(p0Var,p1Var,p2Var, p3Var))
    SetOwnership(p0Var,0)
    SetOwnership(p1Var,0)
    SetOwnership(p2Var,0)
    SetOwnership(p3Var,0)
    return Poly3

  def BuildRooGaussian(self, mean = 125,meanLow = 100, meanHigh = 150, sigma = 1.5, sigmaLow = 0.3, sigmaHigh = 70):

    meanVar = RooRealVar('mean_'+self.suffix,'mean_'+self.suffix, mean, meanLow, meanHigh)
    sigmaVar = RooRealVar('sigma_'+self.suffix,'sigma_'+self.suffix,sigma,sigmaLow,sigmaHigh)
    gauss = RooGaussian('gauss_'+self.suffix,'gauss_'+self.suffix,self.mzg,meanVar,sigmaVar)
    SetOwnership(meanVar,0)
    SetOwnership(sigmaVar,0)
    return gauss

  def BuildCrystalBallGauss(self, piece, mean = 125, meanG = -1, meanGLow = -1, meanGHigh = -1,
                            meanCB = -1, meanCBLow = -1, meanCBHigh = -1,
                            sigmaCB = 1.5, sigmaCBLow = 0.3, sigmaCBHigh = 20, alpha = 1, alphaLow = 0.5, alphaHigh = 10,
                            n = 4, nLow = 0.5, nHigh = 50, sigmaG = 2, sigmaGLow = 0.3, sigmaGHigh = 20,
                            frac = 0.1, fracLow = 0.0, fracHigh = 1.0):

    suffix = self.suffix+'_'+piece
    meanG = meanCB = mean

    if meanGLow   == -1: meanGLow   = meanG-5
    if meanGHigh  == -1: meanGHigh  = meanG+5
    if meanCBLow  == -1: meanCBLow  = meanCB-5
    if meanCBHigh == -1: meanCBHigh = meanCB+5
    meanGVar   = RooRealVar('meanGCBG_'  +suffix,'meanGCBG_'  +suffix, meanG, meanGLow, meanGHigh)
    meanCBVar  = RooRealVar('meanCBCBG_' +suffix,'meanCBCBG_' +suffix, meanCB, meanCBLow, meanCBHigh)
    sigmaCBVar = RooRealVar('sigmaCBCBG_'+suffix,'sigmaCBCBG_'+suffix, sigmaCB,sigmaCBLow,sigmaCBHigh)
    alphaVar   = RooRealVar('alphaCBG_'  +suffix,'alphaCBG_'  +suffix, alpha,alphaLow,alphaHigh)
    nVar       = RooRealVar('nCBG_'      +suffix,'nCBG_'      +suffix, n,nLow,nHigh)
    sigmaGVar  = RooRealVar('sigmaGCBG_' +suffix,'sigmaGCBG_' +suffix, sigmaG,sigmaGLow,sigmaGHigh)
    fracVar    = RooRealVar('fracCBG_'   +suffix,'fracCBG_'   +suffix, frac,fracLow,fracHigh)

    crystal = RooCBShape('crystalCBG_'+suffix,'crystalCBG_'+suffix,self.mzg,meanCBVar,sigmaCBVar,alphaVar,nVar)
    gauss   = RooGaussian('gaussCBG_' +suffix,'gaussCBG_'  +suffix,self.mzg,meanGVar, sigmaGVar)
    cbArgs  = RooArgList(gauss,crystal)
    fracArg = RooArgList(fracVar)
    CBG = RooAddPdf('CBG_'+suffix,'CBG_'+suffix,cbArgs,fracArg,True)

    SetOwnership(meanGVar,0)
    SetOwnership(meanCBVar,0)
    SetOwnership(sigmaCBVar,0)
    SetOwnership(alphaVar,0)
    SetOwnership(nVar,0)
    SetOwnership(sigmaGVar,0)
    SetOwnership(fracVar,0)
    SetOwnership(crystal,0)
    SetOwnership(gauss,0)
    paramList = [meanGVar,meanCBVar,sigmaCBVar,alphaVar,nVar,sigmaGVar,fracVar]
    return CBG, paramList



  def BuildCrystalBallGaussM(self, piece, mean = 125, meanG = -1, meanGLow = -1, meanGHigh = -1,
                             meanCB = -1, meanCBLow = -1, meanCBHigh = -1,
                             sigmaCB = 1.5, sigmaCBLow = 0.3, sigmaCBHigh = 7, alpha = 1, alphaLow = 0.5, alphaHigh = 10,
                             n = 4, nLow = 0.5, nHigh = 50, sigmaG = 2, sigmaGLow = 0.3, sigmaGHigh = 7,
                             frac = 0.1, fracLow = 0.05, fracHigh = 0.4):

    suffix = self.suffix+'_'+piece
    meanG = meanCB = mean

    if meanGLow  == -1: meanGLow  = meanG-5
    if meanGHigh == -1: meanGHigh = meanG+5
    meanVar    = RooRealVar('meanCBG_'   +suffix,'meanCBG_'   +suffix, meanG, meanGLow, meanGHigh)
    sigmaCBVar = RooRealVar('sigmaCBCBG_'+suffix,'sigmaCBCBG_'+suffix, sigmaCB,sigmaCBLow,sigmaCBHigh)
    alphaVar   = RooRealVar('alphaCBG_'  +suffix,'alphaCBG_'  +suffix, alpha,alphaLow,alphaHigh)
    nVar       = RooRealVar('nCBG_'      +suffix,'nCBG_'      +suffix, n,nLow,nHigh)
    sigmaGVar  = RooRealVar('sigmaGCBG_' +suffix,'sigmaGCBG_' +suffix, sigmaG,sigmaGLow,sigmaGHigh)
    fracVar    = RooRealVar('fracCBG_'   +suffix,'fracCBG_'   +suffix, frac,fracLow,fracHigh)

    crystal = RooCBShape('crystalCBG_'+suffix,'crystalCBG_'+suffix,self.mzg,meanVar,sigmaCBVar,alphaVar,nVar)
    gauss   = RooGaussian('gaussCBG_' +suffix,'gaussCBG_'  +suffix,self.mzg,meanVar, sigmaGVar)
    cbArgs  = RooArgList(gauss,crystal)
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

  def BuildTripleGauss(self, piece, mean = 125, mean1 = -1, mean1Low = -1, mean1High = -1, sigma1 = 2, sigma1Low = 1, sigma1High = 8,
                       delta21 = 0, delta21Low = -2, delta21High = 2, s21 = 3, s21Low = 1, s21High = 30,
                       delta31 = 0, delta31Low = -2, delta31High = 2, s32 = 3, s32Low = 1, s32High = 30,
                       frac23 = 0.9, frac23Low = 0, frac23High = 1, frac123 = 0.9, frac123Low = 0, frac123High = 1):

    suffix = self.suffix+'_'+piece
    mean1 = mean
    if mean1Low is -1: mean1Low = mean1-1
    if mean1High is -1: mean1High = mean1+1

    mean1Var = RooRealVar('mean1TripG_'+suffix,'mean1TripG_'+suffix, mean1, mean1Low, mean1High)
    sigma1Var = RooRealVar('sigma1TripG_'+suffix,'sigma1TripG_'+suffix, sigma1, sigma1Low, sigma1High)
    delta21Var = RooRealVar('delta21TripG_'+suffix,'delta21TripG_'+suffix, delta21, delta21Low, delta21High)
    s21Var = RooRealVar('s21TripG_'+suffix,'s21TripG_'+suffix,s21,s21Low,s21High)
    delta31Var = RooRealVar('delta31TripG_'+suffix,'delta31TripG_'+suffix,delta31,delta31Low,delta31High)
    s32Var = RooRealVar('s32TripG_'+suffix,'s32TripG_'+suffix,s32,s32Low,s32High)
    frac23Var = RooRealVar('frac23TripG_'+suffix,'frac23TripG_'+suffix,frac23,frac23Low,frac23High)
    frac123Var = RooRealVar('frac123TripG_'+suffix,'frac123TripG_'+suffix,frac123,frac123Low,frac123High)

    mean2Var = RooFormulaVar('mean2TripG_'+suffix,'@0 + @1', RooArgList(mean1Var, delta21Var))
    sigma2Var = RooFormulaVar('sigma2TripG_'+suffix,'@0 * @1', RooArgList(sigma1Var, s21Var))
    mean3Var = RooFormulaVar('mean3TripG_'+suffix,'@0 + @1', RooArgList(mean1Var, delta31Var))
    sigma3Var = RooFormulaVar('sigma3TripG_'+suffix,'@0 * @1', RooArgList(sigma2Var, s32Var))

    gauss1 = RooGaussian('gauss1TripG_'+suffix,'gauss1TripG_'+suffix,self.mzg,mean1Var,sigma1Var)
    gauss2 = RooGaussian('gauss2TripG_'+suffix,'gauss2TripG_'+suffix,self.mzg,mean2Var,sigma2Var)
    gauss3 = RooGaussian('gauss3TripG_'+suffix,'gauss3TripG_'+suffix,self.mzg,mean3Var,sigma3Var)
    #gaussArgs = RooArgList(gauss1,gauss2,gauss3)
    #fracArgs = RooArgList(frac1Var, frac2Var)
    pdf23 = RooAddPdf('pdf23TripG_'+suffix,'pdf12TripG_'+suffix, gauss2, gauss3, frac23Var)
    TripG = RooAddPdf('TripG_'+suffix,'TripG_'+suffix, gauss1, pdf23, frac123Var)

    SetOwnership(mean1Var,0)
    SetOwnership(sigma1Var,0)
    SetOwnership(mean2Var,0)
    SetOwnership(sigma2Var,0)
    SetOwnership(mean3Var,0)
    SetOwnership(sigma3Var,0)
    SetOwnership(frac23Var,0)
    SetOwnership(frac123Var,0)
    SetOwnership(gauss1,0)
    SetOwnership(gauss2,0)
    SetOwnership(gauss3,0)
    SetOwnership(pdf23,0)
    SetOwnership(delta21Var,0)
    SetOwnership(delta31Var,0)
    SetOwnership(s21Var,0)
    SetOwnership(s32Var,0)
    paramList = [mean1Var, sigma1Var, frac23Var, frac123Var, delta21Var, s21Var, s32Var, delta31Var]
    return TripG, paramList
