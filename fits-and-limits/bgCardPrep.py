#!/usr/bin/env python
import sys
from ROOT import *
gROOT.SetBatch()
from rooFitBuilder import *
sys.path.append("../zgamma")
import utils as u

gROOT.ProcessLine(".L ~/tdrstyle.C")
setTDRStyle()

print len(sys.argv), sys.argv
if len(sys.argv) != 3:
  sys.exit(0)
subdir = sys.argv[1]

leptonList = ['mu']
yearList   = ['2012']
catList    = ['0']


plotBase = "/uscms_data/d2/andreypz/html/zgamma/dalitz/fits-"+subdir+"/bgCardPrep/"
u.createDir(plotBase)
rooWsFile = TFile(subdir+'/testRooFitOut_Dalitz.root')
myWs    = rooWsFile.Get('ws')
card_ws = RooWorkspace('ws_card')
card_ws.autoImportClassCode(True)

c = TCanvas("c","c",0,0,500,400)
c.cd()
mzg = myWs.var('CMS_hzg_mass')
mzg.setRange('signal',120,130)


# #######################################
# prep the background and data card    #
# we're going to the extend the bg pdf #
# and rename the parameters to work    #
# with the higgs combination tool      #
# #######################################

#myWs.Print()

for year in yearList:
  for lepton in leptonList:
    for cat in catList:
      dataName = '_'.join(['data',lepton,year,'cat'+cat])
      suffix   = '_'.join([year,lepton,'cat'+cat])
      print dataName, suffix

      #fitName  = '_'.join(['GaussExp',year,lepton,'cat'+cat])
      #normName = 'normGaussExp_'+suffix

      fitName  = '_'.join(['GaussBern4',year,lepton,'cat'+cat])
      normName = 'normGaussBern4_'+suffix

      # possibly have different fits for separate categories
      if cat is 'a':
        fitName  = '_'.join(['GaussBern6',year,lepton,'cat'+cat])
        normName = 'normGaussBern6_'+suffix


      print fitName, dataName
      data = myWs.data(dataName)
      fit  = myWs.pdf(fitName)

      ###### Extend the fit (give it a normalization parameter)
      print dataName
      sumEntriesBkg = data.sumEntries()
      sumEntriesSig = data.sumEntries('1','signal')
      print sumEntriesBkg, sumEntriesSig

      raw_input()

      dataYieldName = '_'.join(['data','yield',lepton,year,'cat'+cat])
      dataYield     = RooRealVar(dataYieldName,dataYieldName,sumEntriesBkg)
      norm          = RooRealVar(normName,normName,sumEntriesBkg,sumEntriesBkg*0.25,sumEntriesBkg*1.75)

      fitExtName    = '_'.join(['bkgTmp',lepton,year,'cat'+cat])
      fit_ext       = RooExtendPdf(fitExtName,fitExtName, fit,norm)


      fit_ext.fitTo(data,RooFit.Range('dalitzRegion'))

      testFrame = mzg.frame()
      data.plotOn(testFrame, RooFit.Binning(50))
      fit_ext.plotOn(testFrame)
      testFrame.Draw()
      testFrame.SetTitle(";m_{H};unbinned fit pdf")
      c.SaveAs(plotBase+'_'.join(['best_fit',year,lepton,'cat'+cat])+'.png')

      ###### Import the fit and data, and rename them to the card convention
      dataNameNew = '_'.join(['data','obs',lepton,year,'cat'+cat])

      getattr(card_ws,'import')(data,RooFit.Rename(dataNameNew))
      getattr(card_ws,'import')(fit_ext)
      getattr(card_ws,'import')(dataYield)
      card_ws.commitTransaction()
      fit_ext.Print()
      card_ws.Print()

      BackgroundNameFixer(year,lepton,cat,card_ws)

      print "\n * The end * \n"

card_ws.writeToFile(subdir+'/testCardBackground_Dalitz.root')
