#!/usr/bin/env python
import sys
from ROOT import *
gROOT.SetBatch()
#from rooFitBuilder import *
sys.path.append("../zgamma")
import utils as u
import ConfigParser as cp
gROOT.ProcessLine(".L ~/tdrstyle.C")
setTDRStyle()

print len(sys.argv), sys.argv

verbose = 0
doExt   = 0

cf = cp.ConfigParser()
cf.read('config.cfg')
subdir = cf.get("fits","ver")  
yearList   = [a.strip() for a in (cf.get("fits","yearList")).split(',')]
leptonList = [a.strip() for a in (cf.get("fits","leptonList")).split(',')]
catList    = [a.strip() for a in (cf.get("fits","catList")).split(',')]
doBlind    = int(cf.get("fits","blind"))

plotBase = "/uscms_data/d2/andreypz/html/zgamma/dalitz/fits-"+subdir
u.createDir(plotBase)
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
    fitExtName    = '_'.join(['bkgTmp',lepton,year,'cat'+cat])
  else:
    fitExtName = '_'.join(['Bern4',year,lepton,'cat'+cat])

  fitExtNameNew = '_'.join(['bkg',lepton,year,'cat'+cat])

  BernNames = ['Bern3','Bern4','Bern5','Bern6']
  for n in BernNames:
    if n in fitName:
      print "renaming " + fitName
      suffix = '_'.join([year,lepton,'cat'+cat])
      if Ext: normName  = 'norm'+n+'_'+suffix
      p0Name = 'p0'+n+'_'+suffix
      p1Name = 'p1'+n+'_'+suffix
      p2Name = 'p2'+n+'_'+suffix
      p3Name = 'p3'+n+'_'+suffix
      p4Name = 'p4'+n+'_'+suffix

      if Ext: print "Normname from renaming:", normName

      if Ext: normNameNew  = '_'.join(['bkg',lepton,year,'cat'+cat,'norm'])
      p0NameNew = '_'.join(['bkg','p0',lepton,year,'cat'+cat])
      p1NameNew = '_'.join(['bkg','p1',lepton,year,'cat'+cat])
      p2NameNew = '_'.join(['bkg','p2',lepton,year,'cat'+cat])
      p3NameNew = '_'.join(['bkg','p3',lepton,year,'cat'+cat])

      if Ext: ws.factory(normNameNew+'[{0},{1},{2}]'.format(ws.function(normName).getVal(),
                                                            ws.function(normName).getMin(), ws.function(normName).getMax()))
      ws.factory(p0NameNew+'[{0}]'.format(ws.function(p0Name).getVal()))
      ws.factory(p1NameNew+'[{0},{1},{2}]'.format(ws.function(p1Name).getVal(),
                                                  ws.function(p1Name).getMin(),ws.function(p1Name).getMax()))
      ws.factory(p2NameNew+'[{0},{1},{2}]'.format(ws.function(p2Name).getVal(),
                                                  ws.function(p2Name).getMin(),ws.function(p2Name).getMax()))
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

for year in yearList:
  for lepton in leptonList:
    for cat in catList:
      dataName = '_'.join(['data',lepton,year,'cat'+cat])
      suffix   = '_'.join([year,lepton,'cat'+cat])
      print dataName, suffix

      fitName  = '_'.join([bkgModel,year,lepton,'cat'+cat])
      normName = 'norm'+bkgModel+'_'+suffix      

      hPath    = "/eos/uscms/store/user/andreypz/batch_output/zgamma/8TeV/"+subdir
      sigFileMAD  = TFile(hPath+"/mugamma_"+year+"/hhhh_ggH-mad125_1.root", "OPEN")

      Nev = sigFileMAD.Get("Counts/evt_byCut").GetBinContent(2)
      cro = u.getCS("ggH-125",mySel='mu')
      lumi = u.getLumi("2012")
      scale = float(lumi*cro)/Nev
      
      hc7 = sigFileMAD.Get("tri_mass__cut10").Clone()
      hc7.Scale(10*scale)
                        
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
      fit_ext       = RooExtendPdf(fitExtName,fitExtName, fit,norm)

      if verbose:
        print norm.getVal(), norm.getError()
        raw_input("norm.getVal(), norm.getError()")
      
      
      fit_ext.fitTo(data,RooFit.Range('DalitzRegion'))
      myBinning = 30
      binWidth = 2.
      

      testFrame = mzg.frame(RooFit.Range('DalitzRegion'))
      if doBlind:
        data.plotOn(testFrame, RooFit.Binning(myBinning), RooFit.Name('data'), RooFit.CutRange('r1'))
        data.plotOn(testFrame, RooFit.Binning(myBinning), RooFit.Name('data'), RooFit.CutRange('r2'))
      else:
        data.plotOn(testFrame, RooFit.Binning(myBinning), RooFit.Name('data'))

      fit_ext.plotOn(testFrame, RooFit.Name(bkgModel),  RooFit.LineColor(kBlue))
      
      sigDSName = '_'.join(['ds_sig','gg',lepton,year,'cat'+cat,'M125'])
      sigName = 'pdf_sig_mu_2012_cat0_M125'
      myWs.Print()
      sigP  = myWs.pdf(sigName)
      if doBlind:
        sigP.plotOn(testFrame,  RooFit.Name('signal'), RooFit.LineColor(kRed+2),RooFit.Range(115,135),
                    RooFit.Normalization(15.e-02))
      else:
        sigP.plotOn(testFrame,  RooFit.Name('signal'), RooFit.LineColor(kRed+2),RooFit.Range(115,135),
                    RooFit.Normalization(10* 3.33/binWidth/151))

      #sig_dh = RooDataHist("sig_dh","sig histogram",RooArgList(mzg), hc7);
      #sig_dh.plotOn(testFrame,  RooFit.Name('signal2'), RooFit.LineColor(kRed+2),
      #              RooFit.Range(115,135), RooFit.MarkerStyle(0), RooFit.LineStyle(1))
                    
      #hc7.plotOn(testFrame, RooFit.Name('signal2'))
      hc7.Draw('same hist')
  

      if verbose:
        print 'have unit norm?? ', sigP.haveUnitNorm()
        chi2 = testFrame.chiSquare(bkgModel,'data')
        chi2_4 = testFrame.chiSquare(4)
        print ' chiSquare=', chi2, chi2_4
        # print "Figuring out norms of PDFs",sigP.getVal(), sigP.analyticalIntegral()
        raw_input("pdf norm / chi2  ")
        
      testFrame.SetMaximum(50)
      testFrame.Draw()
      testFrame.SetTitle(";m_{H} (GeV);Events/"+str(binWidth)+" GeV")

      leg  = TLegend(0.65,0.7,0.87,0.87)
      leg.SetFillColor(0)
      leg.SetShadowColor(0)
      leg.SetBorderSize(1)
      leg.AddEntry(testFrame.findObject(bkgModel),bkgModel,'l')
      leg.AddEntry(testFrame.findObject('signal'),'10x ggH','l')
      leg.Draw()

      prelim = TLatex(0.15,0.95, "CMS Preliminary, #sqrt{s} = 8 TeV")
      prelim.SetNDC();
      prelim.SetTextSize(0.03)
      prelim.Draw()
              
      c.SaveAs(plotBase+'/'+'_'.join(['best_fit',year,lepton,'cat'+cat])+'.png')

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

      fit_ext.Print()
      fit.Print()
      print 'printing the card'
      card_ws.Print()

      #print normName
      BackgroundNameFixer(fitName, year,lepton,cat,card_ws, doExt)

      card_ws.Print()
            
      print "\n * The end * \n"

card_ws.writeToFile(subdir+'/testCardBackground_Dalitz.root')
