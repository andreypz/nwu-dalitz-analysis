#!/usr/bin/env python
import sys, os
import numpy as np
#import pdb
from rooFitBuilder import *
from ROOT import *
gROOT.SetBatch()
sys.path.append("../zgamma")
import utils as u
class AutoVivification(dict):
  """Implementation of perl's autovivification feature."""
  def __getitem__(self, item):
    try:
      return dict.__getitem__(self, item)
    except KeyError:
      value = self[item] = type(self)()
      return value

import ConfigParser as cp
cf = cp.ConfigParser()
cf.read('config.cfg')
massList        = [a.strip()[0:3] for a in (cf.get("fits","massList")).split(',')]
yearList        = [a.strip() for a in (cf.get("fits","yearList")).split(',')]
leptonList      = [a.strip() for a in (cf.get("fits","leptonList")).split(',')]
catList         = [a.strip() for a in (cf.get("fits","catList")).split(',')]
sigNameList     = [a.strip() for a in (cf.get("fits","sigNameList")).split(',')]

colors = {}
for f,col in cf.items("colors"):
  colors[f] = int(col)

#print colors
#print yearList, leptonList, catList, massList, sigNameList
#print "Begin, INPUT"
#raw_input()

#gSystem.SetIncludePath( "-I$ROOFITSYS/include/" );
#gROOT.ProcessLine(".L ~/tdrstyle.C")
#setTDRStyle()
TH1.SetDefaultSumw2(kTRUE)
gStyle.SetOptTitle(0)

debugPlots = True
verbose    = 0
rootrace   = False

EBetaCut = 1.444
ptMllgCut = 0.30

# OK listen, we're gunna need to do some bullshit just to get a uniform RooRealVar name for our data objects.
# Because I named branches different than the trees, there's a lot of tree loops here
# So we'll loop through the Branch, set mzg to the event value (if it's in range), and add that to our RooDataSet.
# This way, we make a RooDataSet that uses our 'CMS_hzg_mass' variable.

TH1.SetDefaultSumw2(kTRUE)
if rootrace: RooTrace.active(kTRUE)

def LumiXSWeighter(mH, prod, sel, Nev=None):
  cro = float(u.conf.get(prod+"H-"+str(mH), "cs-"+sel))
  if Nev==None:
    Nev = int(u.conf.get(prod+"H-"+str(mH), "Nev-"+sel))
  sc = float(u.lumi*cro)/Nev
  print 'M=',mH, 'cs=',cro, 'Nev=', Nev, 'lumi=',u.lumi, "scale", sc
  return sc


def doInitialFits(subdir):
  print 'loading up the files'

  plotBase = cf.get("path","htmlbase")+'/html/zgamma/dalitz/fits-'+subdir+'/init/'
  u.createDir(plotBase)
  basePath1 = cf.get("path","base")+'/batch_output/zgamma/8TeV/'+subdir+'/'
  basePath2 = cf.get("path","base")+'/zgamma/v34-ele2/'

  if not os.path.exists(basePath1):
    print basePath1, "does not exist!"
    sys.exit(0)
  if not os.path.exists(basePath1):
    print basePath2, "does not exist!"
    sys.exit(0)

  #musel = "mumu"
  musel = "mugamma"

  dataDict   = {'mu2012':TFile(basePath1+'m_Data_'+musel+'_2012.root','r'),
                'el2012':TFile(basePath2+'m_Data_electron_2012.root','r')}

  signalDict = {'gg_mu2012_M120':TFile(basePath1+musel+'_2012/hhhh_ggH-mad120_1.root','r'),
                'gg_mu2012_M125':TFile(basePath1+musel+'_2012/hhhh_ggH-mad125_1.root','r'),
                'gg_mu2012_M130':TFile(basePath1+musel+'_2012/hhhh_ggH-mad130_1.root','r'),
                'gg_mu2012_M135':TFile(basePath1+musel+'_2012/hhhh_ggH-mad135_1.root','r'),
                'gg_mu2012_M140':TFile(basePath1+musel+'_2012/hhhh_ggH-mad140_1.root','r'),
                'gg_mu2012_M145':TFile(basePath1+musel+'_2012/hhhh_ggH-mad145_1.root','r'),
                'gg_mu2012_M150':TFile(basePath1+musel+'_2012/hhhh_ggH-mad150_1.root','r'),

                'vbf_mu2012_M120':TFile(basePath1+musel+'_2012/hhhh_vbf-mad120_1.root','r'),
                'vbf_mu2012_M125':TFile(basePath1+musel+'_2012/hhhh_vbf-mad125_1.root','r'),
                'vbf_mu2012_M130':TFile(basePath1+musel+'_2012/hhhh_vbf-mad130_1.root','r'),
                'vbf_mu2012_M135':TFile(basePath1+musel+'_2012/hhhh_vbf-mad135_1.root','r'),
                'vbf_mu2012_M140':TFile(basePath1+musel+'_2012/hhhh_vbf-mad140_1.root','r'),
                'vbf_mu2012_M145':TFile(basePath1+musel+'_2012/hhhh_vbf-mad145_1.root','r'),
                'vbf_mu2012_M150':TFile(basePath1+musel+'_2012/hhhh_vbf-mad150_1.root','r'),

                'v_mu2012_M120':TFile(basePath1+musel+'_2012/hhhh_vh-mad120_1.root','r'),
                'v_mu2012_M125':TFile(basePath1+musel+'_2012/hhhh_vh-mad125_1.root','r'),
                'v_mu2012_M130':TFile(basePath1+musel+'_2012/hhhh_vh-mad130_1.root','r'),
                'v_mu2012_M135':TFile(basePath1+musel+'_2012/hhhh_vh-mad135_1.root','r'),
                'v_mu2012_M140':TFile(basePath1+musel+'_2012/hhhh_vh-mad140_1.root','r'),
                'v_mu2012_M145':TFile(basePath1+musel+'_2012/hhhh_vh-mad145_1.root','r'),
                'v_mu2012_M150':TFile(basePath1+musel+'_2012/hhhh_vh-mad150_1.root','r')

                #'gg_el2012_M125':TFile(basePath2+'electron_2012/hhhh_dal-mad125_1.root','r')
                }

  treeName = 'fitTree/fitTree'


  lowCutOff  = 110
  highCutOff = 170
  binning = 30
  binWidth = 2

  weight  = RooRealVar('Weight','Weight',0,100)
  mzg  = RooRealVar('CMS_hzg_mass','CMS_hzg_mass', lowCutOff,highCutOff)
  mzg.setRange('FullRegion',   lowCutOff, highCutOff)
  mzg.setRange('DalitzRegion', lowCutOff, highCutOff)
  mzg.setRange('SignalRegion', 122, 128)
  mzg.setBins(50000,'cache')

  c = TCanvas("c","c",0,0,500,400)
  c.cd()

  ws = RooWorkspace("ws")

  yi_da0  = AutoVivification()
  yi_da1  = AutoVivification()
  yi_sig0 = AutoVivification()
  yi_sig1 = AutoVivification()
  mean_sig  = AutoVivification()
  sigma_sig = AutoVivification()
  # ###################################
  # start loop over all year/lep/cat #
  # ###################################
  for year in yearList:
    for lepton in leptonList:
      for cat in catList:
        if rootrace:
          print "rootrace"
          RooTrace.dump()
          raw_input()

        signalList    = []
        signalListDH  = []
        signalListPDF = []
        if verbose: print 'top of loop',year,lepton,cat

        # ##################################################
        # set up the signal histograms and the mzg ranges #
        # ##################################################

        for prod in sigNameList:
          signalListDS = []
          for mass in massList:
            # store the unbinned signals for CB fitting
            signalTree = signalDict[prod+"_"+lepton+year+"_M"+mass].Get(treeName)
            #signalTree.Print()
            sigName = '_'.join(['ds_sig',prod,lepton,year,'cat'+cat,'M'+mass])

            sig_argSW = RooArgSet(mzg,weight)
            sig_ds    = RooDataSet(sigName,sigName,sig_argSW,'Weight')

            Nev = u.getTotalEvents(signalDict[prod+"_"+lepton+year+"_M"+mass])
            lumiWeight = LumiXSWeighter(int(mass), prod, lepton, Nev)
            print signalTree, "A signal tree", prod, '  categ=',cat, "mass =", mass
            for i in signalTree:
              if   cat=="EB" and fabs(i.ph_eta)>EBetaCut: continue
              elif cat=="EE" and fabs(i.ph_eta)<EBetaCut: continue

              if fabs(i.ph_eta) > EBetaCut: continue
              if i.ph_pt/i.m_llg < ptMllgCut or i.di_pt/i.m_llg < ptMllgCut: continue

              if i.m_llg> lowCutOff and i.m_llg<highCutOff:
                mzg.setVal(i.m_llg)

                sigWeight = lumiWeight*i.weight
                #sigWeight = lumiWeight
                sig_ds.add(sig_argSW, sigWeight)
                #sig_argSW.Print()

            #raw_input()
            signalListDS.append(sig_ds)
            getattr(ws,'import')(signalListDS[-1])
            #signalTree.ResetBranchAddresses()

            # do some histogramming for gg signal for bias study
            # we don't need or use unbinned signal or complicated fits
            # but this is mostly for compatibility, we may change to unbinned
            # during a future iteration

            histName  = '_'.join(['sig',prod,lepton,year,'cat'+cat,'M'+mass])
            signalList.append(TH1F(histName, histName, 100, lowCutOff, highCutOff))
            signalList[-1].SetLineColor(kRed)
            signalTree = signalDict[prod+"_"+lepton+year+"_M"+mass].Get(treeName)


            if verbose: print 'signal mass loop', mass
            print 'in gg   INPUT ', cat, prod,mass
            # raw_input()
            
            if verbose:
              print histName
              signalTree.Print()
              print
                
            if cat=="0":
              signalTree.Draw('m_llg>>'+histName,'weight')
            elif cat=="EB":
              signalTree.Draw('m_llg>>'+histName,"weight*(fabs(ph_eta)<"+str(EBetaCut)+")")
            elif cat=="EE":
              signalTree.Draw('m_llg>>'+histName,"weight*(fabs(ph_eta)>"+str(EBetaCut)+")")
              
            if signalList[-1].Integral()!=0:
              signalList[-1].Scale(1/signalList[-1].Integral())
            signalList[-1].Smooth(2)

            if prod=='gg':
              # range is +/- 1 RMS centered around signal peak
              rangeLow = signalList[-1].GetMean()-1.0*signalList[-1].GetRMS()
              rangeHi  = signalList[-1].GetMean()+1.0*signalList[-1].GetRMS()
              rangeName = '_'.join(['range',lepton,year,'cat'+cat,'M'+mass])
              mzg.setRange(rangeName,rangeLow,rangeHi)

              mzg_argL = RooArgList(mzg)
              mzg_argS = RooArgSet(mzg)
              signalListDH.append(RooDataHist('dh_'+histName, 'dh_' +histName,mzg_argL,signalList[-1]))
              signalListPDF.append(RooHistPdf('pdf_'+histName,'pdf_'+histName,mzg_argS,signalListDH[-1],2))
              getattr(ws,'import')(signalListPDF[-1])

              if verbose: print '\n\n ** finished one mass -->>', mass

            yi_sig0[year][mass][lepton][cat][prod] = sig_ds.sumEntries()
            yi_sig1[year][mass][lepton][cat][prod] = sig_ds.sumEntries('1','SignalRegion')

            #mean_sig[year][mass][lepton][cat][prod]  = [signalList[-1].mean(mzg), signalList[-1].GetMeanError()]
            #sigma_sig[year][mass][lepton][cat][prod] = [signalList[-1].GetRMS(), signalList[-1].GetRMSError()]
            mean_sig[year][mass][lepton][cat][prod]  = [signalList[-1].GetMean(), signalList[-1].GetMeanError()]
            sigma_sig[year][mass][lepton][cat][prod] = [signalList[-1].GetRMS(), signalList[-1].GetRMSError()]

            # mean_dh[year][mass][lepton][cat]  = [signalListDH[-1].GetMean(), signalListDH[-1].GetMeanError()]
            # sigma_dh[year][mass][lepton][cat] = [signalListDH[-1].GetRMS(), signalListDH[-1].GetRMSError()]


          print '\n\n Now make some ploOns!\n'
          testFrame = mzg.frame()
          for i,signal in enumerate(signalListPDF):
            signalListDH[i].plotOn(testFrame)
            signal.plotOn(testFrame)
            # signal.paramOn(testFrame)
            # signal.statOn(testFrame)
            # mzg_argS.writeToFile('myargset.txt')
            print i,signal
            testFrame.Draw()
          # raw_input("All DONE. INPUT NEEDED")
          c.SaveAs(plotBase+'_'.join(['signals',prod,year,lepton,'cat'+cat])+'.png')

          testFrame = mzg.frame()
          for signal in signalListDS:
            signal.plotOn(testFrame, RooFit.DrawOption('pl'))
            testFrame.Draw()
          c.SaveAs(plotBase+'_'.join(['ds','sig',prod,year,lepton,'cat'+cat])+'.png')

        del signalTree


        # ###############
        # get the data #
        # ###############
        if verbose: print 'starting data section'

        dataName = '_'.join(['data',lepton,year,'cat'+cat])
        dataTree = dataDict[lepton+year].Get(treeName)

        data_argS = RooArgSet(mzg)
        data_ds   = RooDataSet(dataName,dataName,data_argS)

        for i in dataTree:
          if cat=="EB" and fabs(i.ph_eta)>EBetaCut: continue
          elif cat=="EE" and fabs(i.ph_eta)<EBetaCut: continue

          if fabs(i.ph_eta)>EBetaCut: continue
          if i.ph_pt/i.m_llg < ptMllgCut or i.di_pt/i.m_llg < ptMllgCut: continue

          if i.m_llg> lowCutOff and i.m_llg<highCutOff:
            mzg.setVal(i.m_llg)
            data_ds.add(data_argS)

        #dataTree.ResetBranchAddresses()

        if verbose:
          print dataName
          data_ds.Print()
          print
        if debugPlots:
          testFrame = mzg.frame()
          data_ds.plotOn(testFrame,RooFit.Binning(binning))

          testFrame.Draw()
          c.SaveAs(plotBase+'_'.join(['data',year,lepton,'cat'+cat,'M'+mass])+'.png')

        getattr(ws,'import')(data_ds)



        # ############
        # make fits #
        # ############
        if verbose: 'starting fits'

        if cat!='5':
          #GaussExp = BuildGaussExp(year, lepton, cat, mzg)
          #SechExp  = BuildSechExp(year, lepton, cat, mzg)
          #GaussBern3 = BuildGaussStepBern3(year, lepton, cat, mzg)
          #GaussBern4 = BuildGaussStepBern4(year, lepton, cat, mzg)
          #GaussBern5 = BuildGaussStepBern5(year, lepton, cat, mzg)
          #GaussBern6 = BuildGaussStepBern6(year, lepton, cat, mzg)
          Exp        = BuildExp(year, lepton, cat, mzg)
          Pow        = BuildPow(year, lepton, cat, mzg)
          Bern2      = BuildBern2(year, lepton, cat, mzg)
          Bern3      = BuildBern3(year, lepton, cat, mzg)
          Bern4      = BuildBern4(year, lepton, cat, mzg)
          Bern5      = BuildBern5(year, lepton, cat, mzg)
          Bern6      = BuildBern6(year, lepton, cat, mzg)

          #SechPow  = BuildSechPow(year, lepton, cat, mzg)
          #GaussPow = BuildGaussPow(year, lepton, cat, mzg)

          #GaussBern4  = BuildGaussStepBern4(year, lepton, cat, mzg, step = 105, stepLow = 100, stepHigh = 150, sigma = 2.5)
          #GaussBern5  = BuildGaussStepBern5(year, lepton, cat, mzg, step = 105, stepLow = 100, stepHigh = 150, sigma = 2.5)
          #GaussBern6  = BuildGaussStepBern6(year, lepton, cat, mzg, step = 105, stepLow = 90, stepHigh = 190, sigma = 2.5)
          #SechBern3   = BuildSechStepBern3(year,  lepton, cat, mzg)
          #if lepton == 'mu' and cat == '3': SechBern4 = BuildSechStepBern4(year, lepton, cat, mzg,sigma=2)
          ##else: SechBern4 = BuildSechStepBern4(year, lepton, cat, mzg)
          #if lepton == 'mu' and cat == '3': SechBern5 = BuildSechStepBern5(year, lepton, cat, mzg,sigma=2)
          #else: SechBern5 = BuildSechStepBern5(year, lepton, cat, mzg)

          #gauss = BuildRooGaussian(year, lepton, cat, mzg)
          #BetaFunc = BuildBetaFunc(year, lepton, cat, mzg, 'DalitzRegion')
          #Kumaraswamy = BuildKumaraswamy(year, lepton, cat, mzg, 'DalitzRegion')
                    #BB = BuildBetaAndBern(year, lepton, cat, mzg, 'DalitzRegion')
          #GB = BuildGaussAndBern(year, lepton, cat, mzg, 'DalitzRegion')

          if verbose:
            Exp.Print()
            #GaussExp.Print()
            #GaussPow.Print()
            #SechExp.Print()
            #SechPow.Print()
            #GaussBern3.Print()
            #GaussBern4.Print()
            #GaussBern5.Print()
            #GaussBern6.Print()

          #GaussExp.fitTo(data_ds,RooFit.Range('FullRegion'))
          #SechExp.fitTo(data_ds,RooFit.Range('FullRegion'))
          #GaussBern3.fitTo(data_ds,RooFit.Range('FullRegion'))
          #GaussBern4.fitTo(data_ds,RooFit.Range('FullRegion'))
          #GaussBern5.fitTo(data_ds,RooFit.Range('FullRegion'))
          #GaussBern6.fitTo(data_ds,RooFit.Range('FullRegion'))
          Exp.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          Pow.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          Bern2.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          Bern3.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          Bern4.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          Bern5.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          Bern6.fitTo(data_ds,RooFit.Range('DalitzRegion'))

          #GaussPow.fitTo(data_ds,RooFit.Range('FullRegion'))
          #SechPow.fitTo(data_ds,RooFit.Range('FullRegion'))


          #SechBern3.fitTo(data_ds,RooFit.Range('FullRegion'))
          #SechBern4.fitTo(data_ds,RooFit.Range('FullRegion'))
          #SechBern5.fitTo(data_ds,RooFit.Range('FullRegion'))
          #gauss.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          #BetaFunc.fitTo(data_ds,RooFit.Range('FullRegion'))
          #Kumaraswamy.fitTo(data_ds,RooFit.Range('FullRegion'))

          #BB.fitTo(data_ds,RooFit.Range('FullRegion'))
          #GB.fitTo(data_ds,RooFit.Range('FullRegion'))

          if debugPlots:
            leg  = TLegend(0.7,0.7,1.0,1.0)
            leg.SetFillColor(0)
            leg.SetShadowColor(0)
            leg.SetBorderSize(1)
            leg.SetHeader('_'.join(['test','fits',year,lepton,'cat'+cat]))
            testFrame = mzg.frame(RooFit.Range('DalitzRegion'))
            data_ds.plotOn(testFrame, RooFit.Binning(binning))
            #GaussExp.plotOn(testFrame,RooFit.Name('GaussExp'))
            #SechExp.plotOn(testFrame,RooFit.LineColor(kRed),     RooFit.Name('SechExp'))
            #GaussBern3.plotOn(testFrame,RooFit.LineColor(kGreen),  RooFit.Name('GaussBern3'))
            #GaussBern4.plotOn(testFrame,RooFit.LineColor(kPink),   RooFit.Name('GaussBern4'))
            #GaussBern5.plotOn(testFrame,RooFit.LineColor(kOrange), RooFit.Name('GaussBern5'))
            #GaussBern6.plotOn(testFrame,RooFit.LineColor(kPink),   RooFit.Name('GaussBern6'))
            Exp.plotOn(testFrame,RooFit.LineColor(colors['exp']), RooFit.Name('Exp'))
            Exp.plotOn(testFrame,RooFit.LineColor(colors['pow']), RooFit.Name('Pow'))
            Bern2.plotOn(testFrame,RooFit.LineColor(colors['bern2']), RooFit.Name('Bern2'))
            Bern3.plotOn(testFrame,RooFit.LineColor(colors['bern3']), RooFit.Name('Bern3'))
            Bern4.plotOn(testFrame,RooFit.LineColor(colors['bern4']), RooFit.Name('Bern4'))
            Bern5.plotOn(testFrame,RooFit.LineColor(colors['bern5']), RooFit.Name('Bern5'))
            Bern6.plotOn(testFrame,RooFit.LineColor(colors['bern6']), RooFit.Name('Bern6'))

            #GaussPow.plotOn(testFrame,RooFit.LineColor(kCyan),  RooFit.Name('GaussPow'))
            #SechPow.plotOn(testFrame,RooFit.LineColor(kOrange),  RooFit.Name('SechPow'))
            #SechBern3.plotOn(testFrame,RooFit.LineColor(kMagenta))
            #SechBern4.plotOn(testFrame,RooFit.LineColor(kBlack))
            #SechBern5.plotOn(testFrame,RooFit.LineColor(kGreen))
            #gauss.plotOn(testFrame,RooFit.LineColor(kBlue), RooFit.Name('Gauss'))
            #BetaFunc.plotOn(testFrame,RooFit.LineColor(kBlack), RooFit.Name('Beta'))
            #Kumaraswamy.plotOn(testFrame,RooFit.LineColor(kCyan), RooFit.Name('Kumaraswamy'))

            #BB.plotOn(testFrame,RooFit.LineColor(kViolet), RooFit.Name('Beta+Bern4'))
            #GB.plotOn(testFrame,RooFit.LineColor(kGreen), RooFit.Name('GaussBern3'))
            testFrame.Draw()
            testFrame.SetTitle(";m_{#mu#mu#gamma} (GeV);Events/"+str(binWidth)+" GeV")
            leg.AddEntry(testFrame.findObject('Exp'),'Exp','l')
            leg.AddEntry(testFrame.findObject('Pow'),'Pow','l')
            leg.AddEntry(testFrame.findObject('Bern2'),'Bern2','l')
            leg.AddEntry(testFrame.findObject('Bern3'),'Bern3','l')
            leg.AddEntry(testFrame.findObject('Bern4'),'Bern4','l')
            leg.AddEntry(testFrame.findObject('Bern5'),'Bern5','l')
            leg.AddEntry(testFrame.findObject('Bern6'),'Bern6','l')

            #leg.AddEntry(testFrame.findObject('Beta'),'Beta','l')
            #leg.AddEntry(testFrame.findObject('Kumaraswamy'),'Kumaraswamy','l')
            #leg.AddEntry(testFrame.findObject('Beta+Bern4'),'Beta+Bern4','l')

            #leg.AddEntry(testFrame.findObject('GaussPow'),'GaussPow','l')
            #leg.AddEntry(testFrame.findObject('SechPow'), 'SechPow', 'l')

            #leg.AddEntry(testFrame.findObject('GaussExp'),'GaussExp','l')
            #leg.AddEntry(testFrame.findObject('SechExp'), 'SechExp', 'l')
            #leg.AddEntry(testFrame.findObject('GaussBern3'),'GaussBern3','l')
            #leg.AddEntry(testFrame.findObject('GaussBern4'),'GaussBern4','l')
            #leg.AddEntry(testFrame.findObject('GaussBern5'),'GaussBern5','l')
            #leg.AddEntry(testFrame.findObject('GaussBern6'),'GaussBern6','l')
            leg.Draw()
            c.Print(plotBase+'_'.join(['fits',year,lepton,'cat'+cat])+'.png')


          #raw_input()
          #getattr(ws,'import')(GaussExp)
          #getattr(ws,'import')(SechExp)

          getattr(ws,'import')(Exp)
          getattr(ws,'import')(Pow)
          getattr(ws,'import')(Bern2)
          getattr(ws,'import')(Bern3)
          getattr(ws,'import')(Bern4)
          getattr(ws,'import')(Bern5)
          getattr(ws,'import')(Bern6)
          #getattr(ws,'import')(GaussBern3)
          #getattr(ws,'import')(GaussBern4)
          #getattr(ws,'import')(GaussBern5)
          #getattr(ws,'import')(GaussBern6)

          #getattr(ws,'import')(SechBern3)
          #getattr(ws,'import')(SechBern4)
          #getattr(ws,'import')(SechBern5)
          #getattr(ws,'import')(GB)
          #getattr(ws,'import')(GaussPow)
          #getattr(ws,'import')(SechPow)

        else:
          print cat


        ws.commitTransaction()


        yi_da0[year][lepton][cat] = data_ds.sumEntries()
        yi_da1[year][lepton][cat] = data_ds.sumEntries('1','SignalRegion')




  u.createDir(subdir)
  ws.writeToFile(subdir+'/testRooFitOut_Dalitz.root')

  ws.Print()

  print '*** Some yields from data'
  print 'Full range:', yi_da0
  print 'in 122-128:', yi_da1

  print '*** Some yields from ggH'
  print 'Full range:', yi_sig0
  print 'in 122-128:', yi_sig1


  print "\n ** You should also notice that EBeta cut was: ", EBetaCut
  print "And that the range was from ", lowCutOff, 'to', highCutOff

  print "\n Mean from dataset: ", mean_sig
  print "\n Sigma DS:", sigma_sig

  #print "\n Mean from Datahist:", mean_gg0_dh
  #print "\n Sigma DH:", sigma_gg0_dh


  for year in yearList:
    for lepton in leptonList:
      for cat in catList:
        t_mean = []
        t_sigma = []
        for mass in massList:
          l_mean = []
          l_sigma = []
          l_mean.append(mass)
          l_sigma.append(mass)
          for prod in sigNameList:
            print year, lepton, cat, prod
            print mean_sig[year][mass][lepton][cat][prod]
            print sigma_sig[year][mass][lepton][cat][prod]
            #a,b = mean_sig[year][mass][lepton][cat][prod]
            #l.append("%.2f &pm; %.2f"%(a,b))
            l_mean.append( "%.3f &pm; %.3f"%(mean_sig[year][mass][lepton][cat][prod][0], mean_sig[year][mass][lepton][cat][prod][1]))
            l_sigma.append("%.3f &pm; %.3f"%(sigma_sig[year][mass][lepton][cat][prod][0],sigma_sig[year][mass][lepton][cat][prod][1]))
            #print l
          t_mean.append(l_mean)
          t_sigma.append(l_sigma)
        u.makeTable(t_mean,  "mean",  opt="twiki")
        u.makeTable(t_sigma, "sigma", opt="twiki")

  print '\n \t we did it!\t'



if __name__=="__main__":

  print len(sys.argv)
  print sys.argv
  if len(sys.argv) != 3:
    sys.exit()

  s = sys.argv[1]
  if 'vv/' in s: s = s[3:].rstrip('/')
  print s
  cf.set("path","ver", s)
  if '/tthome' in os.getcwd():
      cf.set("path","base", '/tthome/andrey')
      cf.set("path","htmlbase", '/tthome/andrey')
  else:
      cf.set("path","base", '/eos/uscms/store/user/andreypz')
      cf.set("path","htmlbase", '/uscms_data/d2/andreypz/')

  cf.set("colors","bern2",kCyan+1)
  cf.set("colors","bern3",kOrange)
  cf.set("colors","bern4",kGreen+1)
  cf.set("colors","bern5",kBlue+1)
  cf.set("colors","bern6",kRed+1)

  cf.set("colors","gaussbern3",kSpring-1)
  cf.set("colors","gaussbern4",kGreen-1)
  cf.set("colors","gaussbern5",kBlue-1)
  cf.set("colors","gaussbern6",kRed-1)

  cf.set("colors","gaussexp",kGray)
  cf.set("colors","gausspow",kGray+1)

  with open(r'config.cfg', 'wb') as configfile:
    cf.write(configfile)

  doInitialFits(s)
  print "done"
