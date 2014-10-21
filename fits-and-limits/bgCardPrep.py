#!/usr/bin/env python
import sys,os
from ROOT import *
gROOT.SetBatch()
sys.path.append("../zgamma")
sys.path.append("../scripts")
gROOT.LoadMacro("integralError.C")
from array import *
import utils as u
import makeHTML as ht
import ConfigParser as cp
gROOT.ProcessLine(".L ../tdrstyle.C")
setTDRStyle()
print len(sys.argv), sys.argv
verbose = 0
doExt   = 0

cf = cp.ConfigParser()
cf.read('config.cfg')
subdir = cf.get("path","ver")
yearList    = [a.strip() for a in (cf.get("fits","yearList")).split(',')]
leptonList  = [a.strip() for a in (cf.get("fits","leptonList")).split(',')]
catList     = [a.strip() for a in (cf.get("fits","catList")).split(',')]
sigNameList = [a.strip() for a in (cf.get("fits","sigNameList")).split(',')]
doBlind     = int(cf.get("fits","blind"))
hjp = 0
if 'hjp' in sigNameList:  hjp = 1

mllBins = u.mllBins()

plotBase = cf.get("path","htmlbase")+"/html/zgamma/dalitz/fits-"+subdir
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

if hjp:
  bkgModel = 'Bern2'
else:
  bkgModel = 'Bern4'

# #######################################
# prep the background and data card    #
# we're going to the extend the bg pdf #
# and rename the parameters to work    #
# with the higgs combination tool      #
# #######################################


def BackgroundNameFixer(fitName, year, lepton, cat, ws, Ext=True):
  dataName      = '_'.join(['data',      lepton,year,'cat'+cat])
  dataNameNew   = '_'.join(['data','obs',lepton,year,'cat'+cat])
  if Ext:
    fitExtName = '_'.join(['bkgTmp',lepton,year,'cat'+cat])
  else:
    fitExtName = '_'.join([bkgModel,year,lepton,'cat'+cat])

  fitExtNameNew = '_'.join(['bkg',lepton,year,'cat'+cat])

  BernNames = ['Bern2','Bern3','Bern4','Bern5','Bern6']
  for n in BernNames:
    if n in fitName:
      print "renaming " + fitName
      suffix = '_'.join([year,lepton,'cat'+cat])
      if Ext: normName  = 'norm'+n+'_'+suffix
      p0Name = 'p0'+n+'_'+suffix
      p1Name = 'p1'+n+'_'+suffix
      p2Name = 'p2'+n+'_'+suffix

      if Ext: print "Normname from renaming:", normName

      if Ext: normNameNew  = '_'.join(['bkg',lepton,year,'cat'+cat,'norm'])
      p0NameNew = '_'.join(['bkg','p0',lepton,year,'cat'+cat])
      p1NameNew = '_'.join(['bkg','p1',lepton,year,'cat'+cat])
      p2NameNew = '_'.join(['bkg','p2',lepton,year,'cat'+cat])

      if Ext: ws.factory(normNameNew+'[{0},{1},{2}]'.format(ws.function(normName).getVal(),
                                                            ws.function(normName).getMin(), ws.function(normName).getMax()))
      ws.factory(p0NameNew+'[{0}]'.format(ws.function(p0Name).getVal()))
      ws.factory(p1NameNew+'[{0},{1},{2}]'.format(ws.function(p1Name).getVal(),
                                                  ws.function(p1Name).getMin(),ws.function(p1Name).getMax()))
      ws.factory(p2NameNew+'[{0},{1},{2}]'.format(ws.function(p2Name).getVal(),
                                                  ws.function(p2Name).getMin(),ws.function(p2Name).getMax()))

      if n in ['Bern3','Bern4','Bern5','Bern6']:
        p3Name = 'p3'+n+'_'+suffix
        p3NameNew = '_'.join(['bkg','p3',lepton,year,'cat'+cat])
        ws.factory(p3NameNew+'[{0},{1},{2}]'.format(ws.function(p3Name).getVal(),
                                                    ws.function(p3Name).getMin(),ws.function(p3Name).getMax()))
      if n in ['Bern4','Bern5','Bern6']:
        p4Name = 'p4'+n+'_'+suffix
        p4NameNew = '_'.join(['bkg','p4',lepton,year,'cat'+cat])
        ws.factory(p4NameNew+'[{0},{1},{2}]'.format(ws.function(p4Name).getVal(),
                                                    ws.function(p4Name).getMin(),ws.function(p4Name).getMax()))
      if n in ['Bern5','Bern6']:
        p5Name = 'p5'+n+'_'+suffix
        p5NameNew = '_'.join(['bkg','p5',lepton,year,'cat'+cat])
        ws.factory(p5NameNew+'[{0},{1},{2}]'.format(ws.function(p5Name).getVal(),
                                                    ws.function(p5Name).getMin(),ws.function(p5Name).getMax()))
      if n in ['Bern6']:
        p6Name = 'p6'+n+'_'+suffix
        p6NameNew = '_'.join(['bkg','p6',lepton,year,'cat'+cat])
        ws.factory(p6NameNew+'[{0},{1},{2}]'.format(ws.function(p6Name).getVal(),
                                                    ws.function(p6Name).getMin(),ws.function(p6Name).getMax()))
      if n=='Bern2':
        if Ext:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','+normName+'='+normNameNew+','
                     +p0Name+'='+p0NameNew+','+p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+')')
        else:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','
                     +p0Name+'='+p0NameNew+','+p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+')')
      if n=='Bern3':
        if Ext:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','+normName+'='+normNameNew+','
                     +p0Name+'='+p0NameNew+','+p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+')')
        else:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','
                     +p0Name+'='+p0NameNew+','+p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+')')
      if n=='Bern4':
        if Ext:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','+normName+'='+normNameNew+','
                     +p0Name+'='+p0NameNew+','+p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+')')
        else:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','
                     +p0Name+'='+p0NameNew+','+p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+')')

      elif n =='Bern5':
        if Ext:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','+normName+'='+normNameNew+','
                     +p0Name+'='+p0NameNew+','+p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+','+p5Name+'='+p5NameNew+')')
        else:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','
                     +p0Name+'='+p0NameNew+','+p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+','+p5Name+'='+p5NameNew+')')
      elif n =='Bern6':
        if Ext:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','+normName+'='+normNameNew+','
                     +p0Name+'='+p0NameNew+','+p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+','+p5Name+'='+p5NameNew+p6Name+'='+p6NameNew+')')
        else:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','
                     +p0Name+'='+p0NameNew+','+p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+','+p5Name+'='+p5NameNew+p6Name+'='+p6NameNew+')')

  BernNames = ['GaussBern4','GaussBern5']
  for n in BernNames:
    if n in fitName:
      suffix = '_'.join([year,lepton,'cat'+cat])
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
        if Ext: normNameNew  = '_'.join(['bkg',lepton,year,'cat'+cat,'norm'])
        meanNameNew  = '_'.join(['bkg','mean', lepton,year,'cat'+cat])
        sigmaNameNew = '_'.join(['bkg','sigma',lepton,year,'cat'+cat])
        stepNameNew  = '_'.join(['bkg','step', lepton,year,'cat'+cat])
        p0NameNew = '_'.join(['bkg','p0',lepton,year,'cat'+cat])
        p1NameNew = '_'.join(['bkg','p1',lepton,year,'cat'+cat])
        p2NameNew = '_'.join(['bkg','p2',lepton,year,'cat'+cat])
        p3NameNew = '_'.join(['bkg','p3',lepton,year,'cat'+cat])
        p4NameNew = '_'.join(['bkg','p4',lepton,year,'cat'+cat])
      if n=="GaussBern5":
        p5NameNew = '_'.join(['bkg','p5',lepton,year,'cat'+cat])

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
                     +p0Name+'='+p0NameNew+','+p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+','+p5Name+'='+p5NameNew+')')
        else:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','+meanName+'='+meanNameNew+','
                     +sigmaName+'='+sigmaNameNew+','+stepName+'='+stepNameNew+','
                     +p0Name+'='+p0NameNew+','+p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+','+p5Name+'='+p5NameNew+')')
      elif n=="GaussBern4":
        if Ext:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','+meanName+'='+meanNameNew+','
                     +sigmaName+'='+sigmaNameNew+','+stepName+'='+stepNameNew+','+normName+'='+normNameNew+','
                     +p0Name+'='+p0NameNew+','+p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+')')
        else:
          ws.factory('EDIT::'+fitExtNameNew+'('+fitExtName+','+meanName+'='+meanNameNew+','
                     +sigmaName+'='+sigmaNameNew+','+stepName+'='+stepNameNew+','
                     +p0Name+'='+p0NameNew+','+p1Name+'='+p1NameNew+','+p2Name+'='+p2NameNew+','
                     +p3Name+'='+p3NameNew+','+p4Name+'='+p4NameNew+')')


#myWs.Print()

def doBandsFit(onesigma, twosigma, hmass, cpdf, nomcurve, datanorm, plot, year, lepton):
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
    for lepton in leptonList:
      for cat in catList:
        dataName = '_'.join(['data',lepton,year,'cat'+cat])
        suffix   = '_'.join([year,lepton,'cat'+cat])
        print cat, dataName, suffix

        fitName  = '_'.join([bkgModel,year,lepton,'cat'+cat])
        normName = 'norm'+bkgModel+'_'+suffix

        hPath    = cf.get("path","base")+"/batch_output/zgamma/8TeV/"+subdir
        if hjp:
          sigFile_hjp  = TFile(hPath+"/jp-mugamma_"+year+"/hhhh_HiggsToJPsiGamma_1.root", "OPEN")
          fsig = [sigFile_hjp]
        else:
          if lepton == 'mu': tag = 'mugamma'
          if lepton == 'el': tag = 'elgamma'
          sigFile_gg   = TFile(hPath+"/"+tag+"_"+year+"/hhhh_ggH-mad125_1.root", "OPEN")
          sigFile_vbf  = TFile(hPath+"/"+tag+"_"+year+"/hhhh_vbf-mad125_1.root", "OPEN")
          sigFile_vh   = TFile(hPath+"/"+tag+"_"+year+"/hhhh_vh-mad125_1.root",  "OPEN")
          if lepton == 'mu':
            fsig = [sigFile_gg, sigFile_vbf, sigFile_vh]
          if lepton == 'el':
            fsig = [sigFile_gg]

        hsig = []
        for i,f in enumerate(fsig):
          Nev = f.Get("Counts/evt_byCut").GetBinContent(2)
          if hjp:
            cro = u.getCS("HtoJPsiGamma")/100
          else:
            if i==0:
              cro = u.getCS("ggH-125", mySel=lepton)
            elif i==1:
              cro = u.getCS("vbfH-125",mySel=lepton)
            elif i==2:
              cro = u.getCS("vH-125",  mySel=lepton)

          lumi  = u.getLumi("2012")
          scale = float(lumi*cro)/Nev

          print i, 'File:', str(f.GetName())[-50:], cro, Nev, scale

          if cat=='0':
            hsig.append(f.Get("Main/00_tri_mass_Main_cut9"))
          elif cat=='EB':
            if hjp: hsig.append(f.Get("Main/00_tri_mass_Main_cut10"))
            else:   hsig.append(f.Get("Main/00_tri_mass_Main_cut9"))
          elif cat=='EE':
            hsig.append(f.Get("Main/00_tri_mass_Main_cut16"))
          elif cat=='mll50':
            hsig.append(f.Get("Main/00_tri_mass_Main_cut19"))
          else:
            hsig.append(f.Get("Main/00_tri_mass_Main_cut19"))

          adjust = 1
          if cat in ['m1','m2','m3','m4','m5','m6','m7']:
            adjust=mllBins[int(cat[1])][1]/mllBins[7][1]
          if hjp: factor = 500
          else:   factor = 10
          print len(hsig), hsig, cat
          hsig[-1].Scale(factor*adjust*scale)

        #hsig[-1].Print("all")
        #print hsig
        if not hjp and lepton=='mu':
          print 'hsig0, scaled events:', hsig[0].Integral(), hsig[1].Integral(), hsig[2].Integral()
          #print 'hsigs:', hsig[0], hsig[1], hsig[2]
          hsig[0].Add(hsig[1])
          hsig[0].Add(hsig[2])

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

        dataYieldName = '_'.join(['data','yield',lepton,year,'cat'+cat])
        dataYield     = RooRealVar(dataYieldName,dataYieldName,sumEntriesBkg)
        norm          = RooRealVar(normName,normName,sumEntriesBkg,sumEntriesBkg*0.25,sumEntriesBkg*1.75)

        fitExtName    = '_'.join(['bkgTmp',lepton,year,'cat'+cat])
        fit_ext       = RooExtendPdf(fitExtName,fitExtName, fit, norm)

        fit_result = fit_ext.fitTo(data, RooFit.Range('DalitzRegion'), RooFit.Save(),  RooFit.Optimize(0))
        #fit_result = fit_ext.fitTo(data, RooFit.Range('DalitzRegion'), RooFit.Save())
        #fit_result = RooFitResult(fit_ext.fitTo(data, RooFit.Range('DalitzRegion'), RooFit.Save()))
        #ifit = fit_ext.createIntegral(RooArgSet(mzg), RooFit.Range('signal'))
        #ifit = fit_ext.createIntegral(RooArgSet(mzg), RooFit.NormSet(RooArgSet(mzg)), RooFit.Range('DalitzRegion'))

        l0 = RooArgSet(mzg)
        pars = fit_ext.getParameters(l0)
        #b =  integralError(RooArgSet(mzg), fit_ext, fit_result)
        #print b

        l1, l2 = RooArgList(mzg), RooArgList(pars)
        f1 = fit_ext.asTF(l1, l2)
        #f1 = fit_ext.asTF(RooArgList(mzg), RooArgList(pars))
        integ_full = f1.Integral(110,170)
        integ  = norm.getVal()*f1.Integral(110,170)/integ_full
        covmat = fit_result.covarianceMatrix()
        covmat.Print()
        cormat =  fit_result.correlationMatrix()
        cormat.Print()

        params = f1.GetParameters()

        parlist = [params[a] for a in xrange(f1.GetNpar())]

        #dinteg = norm.getVal()*f1.IntegralError(110, 170, params, covmat.GetMatrixArray()) / integ_full
        #dinteg = f1.IntegralError(110, 170)
        dinteg = f1.IntegralError(110, 170, params, covmat.GetMatrixArray())
        print 'Npar = ', f1.GetNpar()
        print 'Integral full   = %.3f  +/- %.3f'%(integ_full, 0)
        print 'Integral region = %.3f  +/- %.3f'%(integ, dinteg)


        if verbose:
          pars.Print()
          parsiter = pars.createIterator()
          print 'First:', pars.first().GetName(), pars.first().getVal()
          var=parsiter.Next()
          while var!=None:
            print '%s: %.3f  +/ %.3f'%(var.GetName(), var.getVal(), var.getError())
            var=parsiter.Next()

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
          myBinning = 30
          binWidth  = 2.

        testFrame = mzg.frame(RooFit.Range('DalitzRegion'))
        if doBlind:
          data.plotOn(testFrame, RooFit.Binning(myBinning), RooFit.Name('data'), RooFit.CutRange('r1'))
          data.plotOn(testFrame, RooFit.Binning(myBinning), RooFit.Name('data'), RooFit.CutRange('r2'))
        else:
          data.plotOn(testFrame, RooFit.Binning(myBinning), RooFit.Name('data'))


        fit.plotOn(testFrame, RooFit.Name(bkgModel+"2sigma"),
                   RooFit.VisualizeError(fit_result,2), RooFit.FillColor(kCyan-10),RooFit.LineColor(kBlack))
        fit.plotOn(testFrame, RooFit.Name(bkgModel+"1sigma"),
                   RooFit.VisualizeError(fit_result,1), RooFit.FillColor(kCyan-6), RooFit.LineColor(kBlack))
        #fit.plotOn(testFrame, RooFit.Name(bkgModel), RooFit.LineColor(kBlue), RooFit.LineWidth(2))
        fit.plotOn(testFrame, RooFit.Name(bkgModel), RooFit.LineColor(kBlue), RooFit.LineWidth(2))
        #fit.paramOn(testFrame, RooFit.Layout(0.30,0.99,0.9))
        #fit.statOn(testFrame)


        if verbose:
          #print 'have unit norm?? ', sigP.haveUnitNorm()
          chi2 = testFrame.chiSquare(bkgModel,'data')
          chi2_4 = testFrame.chiSquare(4)
          print ' chiSquare=', chi2, chi2_4
          # print "Figuring out norms of PDFs",sigP.getVal(), sigP.analyticalIntegral()
          #raw_input("pdf norm / chi2  ")


        '''
        onesigma = TGraphAsymmErrors()
        twosigma = TGraphAsymmErrors()
        tmpCurve = RooCurve(testFrame.findObject(bkgModel))
        doBandsFit(onesigma, twosigma, mzg, fit, tmpCurve, data, testFrame, year, lepton)
        twosigma.SetLineColor(kBlack)
        twosigma.SetFillColor(kCyan-10)
        onesigma.SetLineColor(kBlack)
        onesigma.SetFillColor(kCyan-6)
        twosigma.Draw("L3 same")
        onesigma.Draw("L3 same")
        '''

        if doBlind:
          data.plotOn(testFrame, RooFit.Binning(myBinning), RooFit.Name('data'), RooFit.CutRange('r1'))
          data.plotOn(testFrame, RooFit.Binning(myBinning), RooFit.Name('data'), RooFit.CutRange('r2'))
          testFrame.SetMinimum(0.1)
        else:
          data.plotOn(testFrame, RooFit.Binning(myBinning), RooFit.Name('data'))

        if hjp:
          testFrame.SetMaximum(12)
        elif cat in ['0']:
          testFrame.SetMaximum(82)
        elif cat in ['m1','m2','m3','m4','m5','m6']:
          testFrame.SetMaximum(22)
        elif cat in ['m7']:
          testFrame.SetMaximum(42)
        else:
          testFrame.SetMaximum(65)

        testFrame.Draw()
        hsig[0].SetAxisRange(115,135,"X")
        hsig[0].SetLineColor(kRed+1)
        hsig[0].SetLineWidth(2)
        hsig[0].Draw('same hist')


        testFrame.SetTitle(";m_{#mu#mu#gamma} (GeV);Events/"+str(binWidth)+" GeV")

        if hjp: leg  = TLegend(0.45,0.62,0.91,0.87)
        else:   leg  = TLegend(0.51,0.62,0.91,0.87)
        leg.SetFillColor(0)
        leg.SetBorderSize(1)
        leg.AddEntry(hsig[0],'Expected signal x'+str(factor),'f')
        #leg.AddEntry(testFrame.findObject(bkgModel+'1sigma'),"Background Model",'f')
        leg.AddEntry(testFrame.findObject(bkgModel),"Background Model",'le')
        #if not hjp: leg.AddEntry(0,'','')
        leg.AddEntry(data,'Data','lep')
        leg.SetTextSize(0.045)
        leg.Draw()

        leg2  = TLegend(0.55,0.62,0.91,0.7)
        leg2.SetNColumns(2)
        leg2.SetFillColor(0)
        leg2.SetBorderSize(0)
        leg2.AddEntry(testFrame.findObject(bkgModel+'1sigma'),"#pm 1 #sigma",'f')
        leg2.AddEntry(testFrame.findObject(bkgModel+'2sigma'),"#pm 2 #sigma",'f')
        leg2.SetTextSize(0.045)
        #leg2.Draw()

        #proc = 'Z#rightarrow J/#Psi#gamma#rightarrow#mu#mu#gamma'
        lat = TLatex()
        lat.SetNDC()
        lat.SetTextSize(0.035)
        if hjp:
          lat.DrawLatex(0.18,0.95, 'H #rightarrow J/#Psi#gamma#rightarrow#mu#mu#gamma')
        else:
          if lepton=='el':
            lat.DrawLatex(0.18,0.95, 'H #rightarrow#gamma*#gamma#rightarrow ee#gamma')
          else:
            lat.DrawLatex(0.18,0.95, 'H #rightarrow#gamma*#gamma#rightarrow#mu#mu#gamma')
        CMS_lumi(c, 2, 11)

        gPad.RedrawAxis()
        for e in ['.png', '.pdf']:
          c.SaveAs(plotBase+'/BestFits/'+'_'.join(['best_fit',year,lepton,'cat'+cat])+e)

        ###### Import the fit and data, and rename them to the card convention
        dataNameNew = '_'.join(['data','obs',lepton,year,'cat'+cat])

        getattr(card_ws,'import')(data,RooFit.Rename(dataNameNew))
        if doExt:
          getattr(card_ws,'import')(fit_ext)
        else:
          getattr(card_ws,'import')(fit)
          normNameFixed = '_'.join(['bkg',lepton,year,'cat'+cat,'norm'])
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
        BackgroundNameFixer(fitName, year,lepton,cat,card_ws, doExt)

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
  ht.makeHTML("h &rarr; dalitz decay plots", plotBase, plot_types, comments, defaultPage, doYields=False)
