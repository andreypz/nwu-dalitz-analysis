#!/usr/bin/env python

import sys,os
sys.path.append("../zgamma")
import utils as u

from ROOT import *
gROOT.SetBatch()
gROOT.LoadMacro("integralError.C")
gROOT.ProcessLine(".L ../tdrstyle.C")
gROOT.LoadMacro("../CMS_lumi.C")
setTDRStyle()

print len(sys.argv), sys.argv
verbose = 0
doExt   = 0

import ConfigParser as cp
cf = cp.ConfigParser()
cf.read('config.cfg')
yearList    = [a.strip() for a in (cf.get("fits","yearList")).split(',')]
leptonList  = [a.strip() for a in (cf.get("fits","leptonList")).split(',')]
catList     = [a.strip() for a in (cf.get("fits","catList")).split(',')]
sigNameList = [a.strip() for a in (cf.get("fits","sigNameList")).split(',')]
doBlind     = int(cf.get("fits","blind"))
hjp = 0
if 'hjp' in sigNameList:
  hjp = 1
  leptonList = ['mu']

mllBins = u.mllBins()

subdir = cf.get("path","ver")
plotBase = cf.get("path","base")+"/fits-"+subdir
u.createDir(plotBase+'/BestFits')

rooWsFile = TFile(subdir+'/testRooFitOut_Dalitz.root')
myWs    = rooWsFile.Get('ws')
card_ws = RooWorkspace('ws_card')
card_ws.autoImportClassCode(True)

c = TCanvas("c","c",0,0,500,400)
c.cd()
mzg = myWs.var('CMS_hzg_mass')
mzg.setRange('signal',120,130)
mzg.setRange('r1',110,120)
mzg.setRange('r2',130,170)

bkgModel = {'hjp': 'Bern2',
            'HR9': 'Bern4',
            'LR9': 'Bern4',
            'VBF': 'Bern4',
            'EE':  'Bern4'}
# #######################################
# prep the background and data card    #
# we're going to the extend the bg pdf #
# and rename the parameters to work    #
# with the higgs combination tool      #
# #######################################


def BackgroundNameFixer(fitName, year, lep, cat, ws, Ext=True):
  dataName      = '_'.join(['data',      lep,year,'cat'+cat])
  dataNameNew   = '_'.join(['data','obs',lep,year,'cat'+cat])
  if Ext:
    fitExtName = '_'.join(['bkgTmp',lep,year,'cat'+cat])
  else:
    fitExtName = '_'.join([bkgModel[cat],year,lep,'cat'+cat])

  fitExtNameNew = '_'.join(['bkg',lep,year,'cat'+cat])

  BernNames = ['Bern2','Bern3','Bern4','Bern5','Bern6']
  for n in BernNames:
    if n in fitName:
      print "renaming " + fitName
      suffix = '_'.join([year,lep,'cat'+cat])
      if Ext: normName  = 'norm'+n+'_'+suffix
      #p0Name = 'p0'+n+'_'+suffix
      p1Name = 'p1'+n+'_'+suffix
      p2Name = 'p2'+n+'_'+suffix

      if Ext: print "Normname from renaming:", normName

      if Ext: normNameNew  = '_'.join(['bkg',lep,year,'cat'+cat,'norm'])
      #p0NameNew = '_'.join(['bkg','p0',lep,year,'cat'+cat])
      p1NameNew = '_'.join(['bkg','p1',lep,year,'cat'+cat])
      p2NameNew = '_'.join(['bkg','p2',lep,year,'cat'+cat])

      if Ext: ws.factory(normNameNew+'[{0},{1},{2}]'.format(ws.function(normName).getVal(),
                                                            ws.function(normName).getMin(), ws.function(normName).getMax()))
      #ws.factory(p0NameNew+'[{0}]'.format(ws.function(p0Name).getVal()))
      ws.factory(p1NameNew+'[{0},{1},{2}]'.format(ws.function(p1Name).getVal(),
                                                  ws.function(p1Name).getMin(),ws.function(p1Name).getMax()))
      ws.factory(p2NameNew+'[{0},{1},{2}]'.format(ws.function(p2Name).getVal(),
                                                  ws.function(p2Name).getMin(),ws.function(p2Name).getMax()))

      if n in ['Bern3','Bern4','Bern5','Bern6']:
        p3Name = 'p3'+n+'_'+suffix
        p3NameNew = '_'.join(['bkg','p3',lep,year,'cat'+cat])
        ws.factory(p3NameNew+'[{0},{1},{2}]'.format(ws.function(p3Name).getVal(),
                                                    ws.function(p3Name).getMin(),ws.function(p3Name).getMax()))
      if n in ['Bern4','Bern5','Bern6']:
        p4Name = 'p4'+n+'_'+suffix
        p4NameNew = '_'.join(['bkg','p4',lep,year,'cat'+cat])
        ws.factory(p4NameNew+'[{0},{1},{2}]'.format(ws.function(p4Name).getVal(),
                                                    ws.function(p4Name).getMin(),ws.function(p4Name).getMax()))
      if n in ['Bern5','Bern6']:
        p5Name = 'p5'+n+'_'+suffix
        p5NameNew = '_'.join(['bkg','p5',lep,year,'cat'+cat])
        ws.factory(p5NameNew+'[{0},{1},{2}]'.format(ws.function(p5Name).getVal(),
                                                    ws.function(p5Name).getMin(),ws.function(p5Name).getMax()))
      if n in ['Bern6']:
        p6Name = 'p6'+n+'_'+suffix
        p6NameNew = '_'.join(['bkg','p6',lep,year,'cat'+cat])
        ws.factory(p6NameNew+'[{0},{1},{2}]'.format(ws.function(p6Name).getVal(),
                                                    ws.function(p6Name).getMin(),ws.function(p6Name).getMax()))
      if n=='Bern2':
        if Ext:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','+normName+'='+normNameNew+','
                     +p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+')')
        else:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','
                     +p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+')')
      if n=='Bern3':
        if Ext:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','+normName+'='+normNameNew+','
                     +p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+')')
        else:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','
                     +p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+')')
      if n=='Bern4':
        if Ext:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','+normName+'='+normNameNew+','
                     +p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+')')
        else:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','
                     +p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+')')

      elif n =='Bern5':
        if Ext:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','+normName+'='+normNameNew+','
                     +p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+','+p5Name+'='+p5NameNew+')')
        else:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','
                     +p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+','+p5Name+'='+p5NameNew+')')
      elif n =='Bern6':
        if Ext:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','+normName+'='+normNameNew+','
                     +p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+','+p5Name+'='+p5NameNew+p6Name+'='+p6NameNew+')')
        else:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','
                     +p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+','+p5Name+'='+p5NameNew+p6Name+'='+p6NameNew+')')

  BernNames = ['GaussBern4','GaussBern5']
  for n in BernNames:
    if n in fitName:
      suffix = '_'.join([year,lep,'cat'+cat])
      if Ext: normName  = 'norm'+n+'_'+suffix
      meanName  = 'mean'+n+'_'+suffix
      sigmaName = 'sigma'+n+'_'+suffix
      stepName  = 'step'+n+'_'+suffix
      p0Name = 'p0'+n+'_'+suffix
      p1Name = 'p1'+n+'_'+suffix
      p2Name = 'p2'+n+'_'+suffix
      p3Name = 'p3'+n+'_'+suffix
      p4Name = 'p4'+n+'_'+suffix
      if n=="GaussBern5":
        p5Name = 'p5'+n+'_'+suffix
        if Ext: normNameNew  = '_'.join(['bkg',lep,year,'cat'+cat,'norm'])
        meanNameNew  = '_'.join(['bkg','mean', lep,year,'cat'+cat])
        sigmaNameNew = '_'.join(['bkg','sigma',lep,year,'cat'+cat])
        stepNameNew  = '_'.join(['bkg','step', lep,year,'cat'+cat])
        p0NameNew = '_'.join(['bkg','p0',lep,year,'cat'+cat])
        p1NameNew = '_'.join(['bkg','p1',lep,year,'cat'+cat])
        p2NameNew = '_'.join(['bkg','p2',lep,year,'cat'+cat])
        p3NameNew = '_'.join(['bkg','p3',lep,year,'cat'+cat])
        p4NameNew = '_'.join(['bkg','p4',lep,year,'cat'+cat])
      if n=="GaussBern5":
        p5NameNew = '_'.join(['bkg','p5',lep,year,'cat'+cat])

        if Ext: ws.factory(normNameNew+'[{0},{1},{2}]'.format(ws.function(normName).getVal(),
                                                      ws.function(normName).getMin(), ws.function(normName).getMax()))
        ws.factory(meanNameNew+'[{0}]'.format(ws.function(meanName).getVal()))
        ws.factory(sigmaNameNew+'[{0},{1},{2}]'.format(ws.function(sigmaName).getVal(),
                                                       ws.function(sigmaName).getMin(),ws.function(sigmaName).getMax()))
        ws.factory(stepNameNew+'[{0},{1},{2}]'.format(ws.function(stepName).getVal(),
                                                      ws.function(stepName).getMin(),ws.function(stepName).getMax()))
        ws.factory(p0NameNew+'[{0}]'.format(ws.function(p0Name).getVal()))
        ws.factory(p1NameNew+'[{0},{1},{2}]'.format(ws.function(p1Name).getVal(),
                                                    ws.function(p1Name).getMin(),ws.function(p1Name).getMax()))
        ws.factory(p2NameNew+'[{0},{1},{2}]'.format(ws.function(p2Name).getVal(),
                                                    ws.function(p2Name).getMin(),ws.function(p2Name).getMax()))
        ws.factory(p3NameNew+'[{0},{1},{2}]'.format(ws.function(p3Name).getVal(),
                                                    ws.function(p3Name).getMin(),ws.function(p3Name).getMax()))
        ws.factory(p4NameNew+'[{0},{1},{2}]'.format(ws.function(p4Name).getVal(),
                                                    ws.function(p4Name).getMin(),ws.function(p4Name).getMax()))
      if n=="GaussBern5":
        ws.factory(p5NameNew+'[{0},{1},{2}]'.format(ws.function(p5Name).getVal(),
                                                    ws.function(p5Name).getMin(),ws.function(p4Name).getMax()))

        if Ext:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','+meanName+'='+meanNameNew+','
                     +sigmaName+'='+sigmaNameNew+','+stepName+'='+stepNameNew+','+normName+'='+normNameNew+','
                     +p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+','+p5Name+'='+p5NameNew+')')
        else:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','+meanName+'='+meanNameNew+','
                     +sigmaName+'='+sigmaNameNew+','+stepName+'='+stepNameNew+','
                     +p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+','+p5Name+'='+p5NameNew+')')
      elif n=="GaussBern4":
        if Ext:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','+meanName+'='+meanNameNew+','
                     +sigmaName+'='+sigmaNameNew+','+stepName+'='+stepNameNew+','+normName+'='+normNameNew+','
                     +p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+')')
        else:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','+meanName+'='+meanNameNew+','
                     +sigmaName+'='+sigmaNameNew+','+stepName+'='+stepNameNew+','
                     +p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+')')
  print '\n \t\t ** END OF RENAMING  ***\n'

#myWs.Print()

def doBandsFit(onesigma, twosigma, hmass, cpdf, nomcurve, datanorm, plot, year, lep):
  print '\n \t \t *** starting bands \n'
  nlim = RooRealVar("nlim","", 0, 0,100)
  print 'total steps needed:', plot.GetXaxis().GetNbins()
  oldhi = oldlo = 9999

  cpdf.Print()
  datanorm.Print()
  nomcurve.Print()

  raw_input('Enter')

  for i in range(1,plot.GetXaxis().GetNbins()+1):
    #r = TRandom(i)
    lowedge = plot.GetXaxis().GetBinLowEdge(i)
    upedge  = plot.GetXaxis().GetBinUpEdge(i)
    center  = plot.GetXaxis().GetBinCenter(i)
    nombkg  = nomcurve.interpolate(center)
    print 'bin number:', i, 'nombkg=', nombkg
    nlim.setVal(nombkg)

    #hmass.removeRange()
    #nlim.removeRange()
    nlim.setRange(nombkg*0.3, nombkg*3.0)
    hmass.setRange("errRange",lowedge,upedge)
    #hmass.Print()
    epdf = RooExtendPdf("epdf","",cpdf,nlim,"errRange")

    cl_one = 1.0 - 2.0*(RooStats.SignificanceToPValue(1.0))
    cl_two = 1.0 - 2.0*(RooStats.SignificanceToPValue(2.0))
    #print 'cl_one', cl_one, 'cl_two', cl_two
    onesigma.SetPoint(i-1,center,nombkg)
    tempi = i

    #nll = epdf.createNLL(datanorm)
    nll = epdf.createNLL(datanorm, RooFit.Extended())
    #nll.Print()
    minim = RooMinimizer(nll)
    minim.setErrorLevel(0.5*pow(ROOT.Math.normal_quantile(1-0.5*(1-cl_one),1.0),2)) #0.5 is because qmu is -2*NLL
    minim.setStrategy(2)
    #minim.setPrintLevel(-1)
    minim.migrad()
    minim.hesse()
    minim.minos(RooArgSet(nlim))
    onelo = -nlim.getErrorLo()
    onehi =  nlim.getErrorHi()
    val = nlim.getVal()

    onesigma.SetPointError(i-1,0.,0.,onelo,onehi)
    if fabs(onelo)<0.01:
      onesigma.SetPointError(i-1,0.,0.,onehi,onehi)
      onelo=onehi
    if fabs(onehi)<0.01:
      onesigma.SetPointError(i-1,0.,0.,onelo,onelo)
      onehi=onelo
    if(fabs(onelo) <0.01 and fabs(onehi)<0.01):
      onesigma.SetPointError(i-1,0.,0.,nlim.getError(),nlim.getError())
      onelo=nlim.getError()
      onehi=nlim.getError()

    print 'val=', val, 'one errHi', onehi, 'one errLo', onelo

    #minim.setErrorLevel(0.5*pow(ROOT.Math.normal_quantile(1-0.5*(1-cltwo),1.0),2)) #0.5 is because qmu is -2*NLL
    # eventually if cl = 0.95 this is the usual 1.92!
    twosigma.SetPoint(i-1,center,nombkg)
    twolo = 1.92*onelo
    twohi = 1.92*onehi
    twosigma.SetPointError(i-1,0.,0.,twolo,twohi)

    print 'two errHi', twohi, 'two errLo', twolo
  onesigma.Print("V")
  twosigma.Print("V")

  raw_input('Enter ')


if __name__=="__main__":
  for year in yearList:
    for lep in leptonList:
      for cat in catList:
        if lep=='el' and cat!='EB': continue

        dataName = '_'.join(['data',lep,year,'cat'+cat])
        suffix   = '_'.join([year,lep,'cat'+cat])
        print cat, dataName, suffix

        fitName  = '_'.join([bkgModel[cat],year,lep,'cat'+cat])
        normName = 'norm'+bkgModel[cat]+'_'+suffix

        if hjp: fs125 = TFile(subdir+'/s125-hjp.root','open')
        else:   fs125 = TFile(subdir+'/s125-'+lep+'.root','open')
        fs125.Print()
        fs125.ls()
        hsig = []
        if lep=='mu': factor=10
        if lep=='el': factor=10
        if hjp:       factor=500
        sigNameList   = [a.strip() for a in (cf.get("fits","sigNameList")).split(',')]
        for prod in sigNameList:
          if lep=='el' and prod=='v': continue
          if cat in ['m1','m2','m3','m4','m5','m6','m7']:
            if prod!='gg' or lep=='el': continue

          histName  = '_'.join(['sig',prod,lep,year,'cat'+cat,'M125','_CMS_hzg_mass'])
          hsig.append(fs125.Get(histName))
          print prod, lep, hsig
          hsig[-1].Print("all")
          hsig[-1].Scale(factor)

        # print hsig
        if not hjp:
          for i in range(1,len(hsig)):
            print i, 'hsig, scaled events:', hsig[i].Integral()
            hsig[0].Add(hsig[i])
            print "total signal now = ", hsig[0].Integral()

        print fitName, dataName
        data = myWs.data(dataName)
        fit  = myWs.pdf(fitName)

        ###### Extend the fit (give it a normalization parameter)
        print dataName
        sumEntriesBkg = data.sumEntries()
        sumEntriesSig = data.sumEntries('1','signal')

        if verbose:
          print sumEntriesBkg, sumEntriesSig
          raw_input("sumEntriesBkg and sumEntriesSig")

        dataYieldName = '_'.join(['data','yield',lep,year,'cat'+cat])
        dataYield     = RooRealVar(dataYieldName,dataYieldName,sumEntriesBkg)
        norm          = RooRealVar(normName,normName,sumEntriesBkg,sumEntriesBkg*0.25,sumEntriesBkg*1.75)

        fitExtName    = '_'.join(['bkgTmp',lep,year,'cat'+cat])
        fit_ext       = RooExtendPdf(fitExtName,fitExtName, fit, norm)

        fit_result = fit_ext.fitTo(data, RooFit.Range('DalitzRegion'), RooFit.Save())
        #fit_result = fit.fitTo(data, RooFit.Range('DalitzRegion'), RooFit.Save())

        if verbose:
          """ Calculating integrals and their errors!"""
          l0 = RooArgSet(mzg)
          # pars = fit.getParameters(l0)
          pars = fit_ext.getParameters(l0)
          # b =  integralError(RooArgSet(mzg), fit_ext, fit_result)
          # print b
          covmat = fit_result.covarianceMatrix()
          covmat.Print()
          pars.Print()
          parsiter = pars.createIterator()
          print 'First:', pars.first().GetName(), pars.first().getVal()
          var=parsiter.Next()
          while var!=None:
            print '%s: %.3f  +/ %.3f'%(var.GetName(), var.getVal(), var.getError())
            var=parsiter.Next()

          l1, l2 = RooArgList(mzg), RooArgList(pars)
          # f1 = fit.asTF(l1, l2)
          f1 = fit_ext.asTF(l1, l2)
          m1_f = mzg.getMin('DalitzRegion')
          m2_f = mzg.getMax('DalitzRegion')
          m1_s = mzg.getMin('signal')
          m2_s = mzg.getMax('signal')
          integ_norm = f1.Integral(m1_f,m2_f)
          integ_full = norm.getVal()*f1.Integral(m1_f,m2_f)/integ_norm
          integ_sig  = norm.getVal()*f1.Integral(m1_s,m2_s)/integ_norm

          params = f1.GetParameters()
          # parlist = [params[a] for a in xrange(f1.GetNpar())]

          dinteg_f = norm.getVal()*f1.IntegralError(m1_f, m2_f, params, covmat.GetMatrixArray()) / integ_norm
          dinteg_s = norm.getVal()*f1.IntegralError(m1_s, m2_s, params, covmat.GetMatrixArray()) / integ_norm
          # dinteg = f1.IntegralError(110, 170, params, covmat.GetMatrixArray())
          print 'Npar = ', f1.GetNpar()
          for p in xrange(f1.GetNpar()):
            print f1.GetParName(p), f1.GetParameter(p)

          print 'Integral full = %.3f  +/- %.3f'%(integ_full, dinteg_f)
          print 'Signal region = %.3f  +/- %.3f'%(integ_sig,  dinteg_s)


          #for p in parsiter:
          #  print 'Parameter', p.GetName(), p.getVal()
          print 'norm.getVal() = %.3f, norm.getError() = %.3f' %(norm.getVal(), norm.getError())
          #print 'Integral in signal region:', ifit.getVal()
          #print 'ratio? = ', norm.getVal()/ifit.getVal()
          raw_input("Enter to continue")

        if hjp:
          myBinning = 20
          binWidth  = 2.
        else:
          myBinning = 60
          binWidth  = 1.

        testFrame = mzg.frame(RooFit.Range('DalitzRegion'))
        if doBlind:
          data.plotOn(testFrame, RooFit.Binning(myBinning), RooFit.Name('data1'), RooFit.CutRange('r1'))
          data.plotOn(testFrame, RooFit.Binning(myBinning), RooFit.Name('data2'), RooFit.CutRange('r2'))
        else:
          data.plotOn(testFrame, RooFit.Binning(myBinning), RooFit.Invisible(), RooFit.Name('data0'))


        fit.plotOn(testFrame, RooFit.Name(bkgModel[cat]+"2sigma"),
                   RooFit.VisualizeError(fit_result,2), RooFit.FillColor(kCyan-10),RooFit.LineColor(kCyan-10))
        fit.plotOn(testFrame, RooFit.Name(bkgModel[cat]+"1sigma"),
                   RooFit.VisualizeError(fit_result,1), RooFit.FillColor(kCyan-6), RooFit.LineColor(kCyan-6))
        #fit.plotOn(testFrame, RooFit.Name(bkgModel[cat]), RooFit.LineColor(kBlue), RooFit.LineWidth(2))
        fit.plotOn(testFrame, RooFit.Name(bkgModel[cat]), RooFit.LineColor(kBlue), RooFit.LineWidth(2), RooFit.FillColor(kCyan-6))
        #fit.paramOn(testFrame, RooFit.Layout(0.30,0.99,0.9))
        #fit.statOn(testFrame)


        if verbose:
          #print 'have unit norm?? ', sigP.haveUnitNorm()
          if doBlind:
            print "\n\t WARNING: your Chi2 would not make sence when Blinded!"
          chi2 = testFrame.chiSquare(bkgModel[cat],'data')
          for a in xrange(6):
            print 'nDof = ', a, testFrame.chiSquare(bkgModel[cat],'data0',a)
            print testFrame.chiSquare(a)
            #print "Figuring out norms of PDFs",sigP.getVal(), sigP.analyticalIntegral()
          raw_input("Enter to continue ")


        '''
        onesigma = TGraphAsymmErrors()
        twosigma = TGraphAsymmErrors()
        tmpCurve = RooCurve(testFrame.findObject(bkgModel[cat]))
        doBandsFit(onesigma, twosigma, mzg, fit, tmpCurve, data, testFrame, year, lep)
        twosigma.SetLineColor(kBlack)
        twosigma.SetFillColor(kCyan-10)
        onesigma.SetLineColor(kBlack)
        onesigma.SetFillColor(kCyan-6)
        twosigma.Draw("L3 same")
        onesigma.Draw("L3 same")
        '''

        testFrame.SetMinimum(0)
        if doBlind:
          data.plotOn(testFrame, RooFit.Binning(myBinning), RooFit.Name('data1'), RooFit.CutRange('r1'))
          data.plotOn(testFrame, RooFit.Binning(myBinning), RooFit.Name('data2'), RooFit.CutRange('r2'))
          testFrame.SetMinimum(0.1)
        else:
          data.plotOn(testFrame, RooFit.Binning(myBinning), RooFit.XErrorSize(0), RooFit.Name('data'))

        if hjp:
          testFrame.SetMaximum(12)
        elif cat in ['0']:
          testFrame.SetMaximum(82)
        elif cat in ['m1','m2','m3','m4','m5','m6']:
          testFrame.SetMaximum(22)
        elif cat in ['m7']:
          testFrame.SetMaximum(42)
        elif lep=='el':
          testFrame.SetMaximum(35)
        elif cat=='VBF':
          testFrame.SetMaximum(12)
        else:
          testFrame.SetMaximum(72)

        testFrame.Draw()
        hsig[0].SetAxisRange(115,135,"X")
        hsig[0].SetLineColor(kRed+1)
        hsig[0].SetLineWidth(2)
        hsig[0].Draw('same hist')


        if lep=='mu':
          testFrame.SetTitle(";m_{#mu#mu#gamma} (GeV);Events/"+str(binWidth)+" GeV")
        elif lep=='el':
          testFrame.SetTitle(";m_{ee#gamma} (GeV);Events/"+str(binWidth)+" GeV")

        if hjp:
          leg   = TLegend(0.40,0.62,0.91,0.88)
          leg2  = TLegend(0.51,0.63,0.81,0.89)
        else:
          leg   = TLegend(0.47,0.62,0.92,0.88)
          leg2  = TLegend(0.58,0.63,0.88,0.89)

        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg2.SetFillStyle(0)
        leg2.SetBorderSize(0)


        leg.AddEntry(testFrame.findObject('data'),'Data','p')
        leg.AddEntry(testFrame.findObject(bkgModel[cat]),"Background model",'l')
        leg.AddEntry(0,'','')
        if hjp:
          leg.AddEntry(hsig[0],str(factor)+'x SM H#rightarrow(J/#Psi)#gamma #rightarrow #mu#mu#gamma','f')
        else:
          if lep=='el':
            leg.AddEntry(hsig[0], str(factor) + 'x SM H #rightarrow #gamma*#gamma #rightarrow ee#gamma','f')
          else:
            leg.AddEntry(hsig[0], str(factor) + 'x SM H #rightarrow #gamma*#gamma #rightarrow #mu#mu#gamma','f')
        #leg.AddEntry(testFrame.findObject(bkgModel[cat]+'1sigma'),"Background Model",'le')
        #if not hjp:
        leg.SetTextFont(42)
        print '\n \t \t Legend font: ', gStyle.GetLegendFont()

        leg.SetTextSize(0.045)
        leg.Draw()


        leg2.SetNColumns(2)
        leg2.AddEntry(0,'','')
        leg2.AddEntry(0,'','')
        leg2.AddEntry(0,'','')
        leg2.AddEntry(0,'','')
        leg2.AddEntry(testFrame.findObject(bkgModel[cat]+'1sigma'),'#pm 1 #sigma','f')
        leg2.AddEntry(testFrame.findObject(bkgModel[cat]+'2sigma'),"#pm 2 #sigma",'f')
        leg2.AddEntry(0,'','')
        leg2.AddEntry(0,'','')
        leg2.SetTextSize(0.037)

        leg2.Draw()


        #proc = 'Z#rightarrow J/#Psi#gamma#rightarrow#mu#mu#gamma'
        """
        lat = TLatex()
        lat.SetNDC()
        lat.SetTextFont(42)
        lat.SetTextSize(0.035)
        if hjp:
          lat.DrawLatex(0.18,0.95, 'H #rightarrow J/#Psi#gamma #rightarrow #mu#mu#gamma')
        else:
          if lep=='el':
            lat.DrawLatex(0.18,0.95, 'H #rightarrow #gamma*#gamma #rightarrow ee#gamma')
          else:
            lat.DrawLatex(0.18,0.95, 'H #rightarrow #gamma*#gamma #rightarrow #mu#mu#gamma')
        """
        CMS_lumi(c, 4, 11,"")

        gPad.RedrawAxis()
        for e in ['.png', '.pdf']:
          if hjp: sel = 'hjp'
          else: sel = lep
          c.SaveAs(plotBase+'/BestFits/'+'_'.join(['best_fit',year,sel,'cat'+cat])+e)

        ###### Import the fit and data, and rename them to the card convention
        dataNameNew = '_'.join(['data','obs',lep,year,'cat'+cat])

        getattr(card_ws,'import')(data,RooFit.Rename(dataNameNew))
        if doExt:
          getattr(card_ws,'import')(fit_ext)
        else:
          getattr(card_ws,'import')(fit)
          normNameFixed = '_'.join(['bkg',lep,year,'cat'+cat,'norm'])
          norm.SetName(normNameFixed)
          getattr(card_ws,'import')(norm)

        getattr(card_ws,'import')(dataYield)
        card_ws.commitTransaction()

        print '\t\t fit_ext'
        fit_ext.Print()
        print '\t\t fit'
        fit.Print()
        print 'printing the WS before renaming!'
        card_ws.Print()

        #print normName
        BackgroundNameFixer(fitName, year,lep,cat,card_ws, doExt)

        print 'Now print it After renaming'
        card_ws.Print()


        print "\n * The end * \n"

  card_ws.writeToFile(subdir+'/testCardBackground_Dalitz.root')

  comments = ['Limits and fits']
  defaultPage = 'BestFits'

  plot_types =[]
  dirlist = os.listdir(plotBase)
  for d in dirlist:
    if os.path.isdir(plotBase+"/"+d):
      plot_types.append(d)
  # ht.makeHTML("h &rarr; dalitz decay plots", plotBase, plot_types, comments, defaultPage, doYields=False)
