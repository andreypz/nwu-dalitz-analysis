#!/usr/bin/env python
import sys, os
import numpy as np
#import pdb
from rooFitBuilder import *
from ROOT import *
gROOT.SetBatch()
sys.path.append("../zgamma")
import utils as u

gSystem.SetIncludePath( "-I$ROOFITSYS/include/" );
gROOT.ProcessLine(".L ~/tdrstyle.C")
setTDRStyle()
TH1.SetDefaultSumw2(kTRUE)
gStyle.SetOptTitle(0)

debugPlots = True
verbose    = 0
rootrace   = False

# OK listen, we're gunna need to do some bullshit just to get a uniform RooRealVar name for our data objects.
# Because I named branches different than the trees, there's a lot of tree loops here
# So we'll loop through the Branch, set mzg to the event value (if it's in range), and add that to our RooDataSet.
# This way, we make a RooDataSet that uses our 'CMS_hzg_mass' variable.

TH1.SetDefaultSumw2(kTRUE)
if rootrace: RooTrace.active(kTRUE)

def LumiXSWeighter(mH, sel):
  #Mad:
  #cro = 0.00091
  cro = float(u.conf.get("ggH-"+str(mH), "cs-"+sel))
  Nev = int(u.conf.get("ggH-"+str(mH), "Nev-"+sel))
  
  sc = float(u.lumi*cro)/Nev
  print 'M=',mH, 'cs=',cro, 'Nev=', Nev, 'lumi=',u.lumi, "scale", sc
  return sc


def doInitialFits(subdir):
  print 'loading up the files'
    
  plotBase = '/uscms_data/d2/andreypz/html/zgamma/dalitz/fits-'+subdir+'/init/'
  u.createDir(plotBase)
  #basePath1 = '/uscms_data/d2/andreypz/zgamma/'+subdir+'/'
  basePath1 = '/eos/uscms/store/user/andreypz/batch_output/zgamma/8TeV/'+subdir+'/'
  #basePath2 = '/eos/uscms/store/user/andreypz/batch_output/zgamma/8TeV/v30/'
  basePath2 = '/uscms_data/d2/andreypz/zgamma/v34-ele2/'
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

  signalDict = {'mu2012_M120':TFile(basePath1+musel+'_2012/hhhh_dal-mad120_1.root','r'),
                'mu2012_M125':TFile(basePath1+musel+'_2012/hhhh_dal-mad125_1.root','r'),
                'mu2012_M130':TFile(basePath1+musel+'_2012/hhhh_dal-mad130_1.root','r'),
                'mu2012_M135':TFile(basePath1+musel+'_2012/hhhh_dal-mad135_1.root','r'),
                'mu2012_M140':TFile(basePath1+musel+'_2012/hhhh_dal-mad140_1.root','r'),
                'mu2012_M145':TFile(basePath1+musel+'_2012/hhhh_dal-mad145_1.root','r'),
                'mu2012_M150':TFile(basePath1+musel+'_2012/hhhh_dal-mad150_1.root','r'),

                'el2012_M125':TFile(basePath2+'electron_2012/hhhh_dal-mad125_1.root','r')}

  treeName = 'fitTree/fitTree'

  leptonList  = ['mu']
  yearList    = ['2012']
  #catList     = ['0',"EB","EE","LowPt"]
  catList     = ['0']
  EBetaCut = 1.0 
  category    = {"0":0,"EB":1,"EE":2,"Low":3,"LowMll":4,"HighMll":5}
  #massList    = ['125']
  massList    = ['120','125','130','135','140','145','150']
  sigNameList = ['gg']

  
  weight  = RooRealVar('Weight','Weight',0,100)
  mzg  = RooRealVar('CMS_hzg_mass','CMS_hzg_mass',90,190)
  mzg.setRange('fullRegion', 80,200)
  mzg.setRange('DalitzRegion',  90,190)
  mzg.setBins(50000,'cache')

  c = TCanvas("c","c",0,0,500,400)
  c.cd()

  ws =RooWorkspace("ws")

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
          
        signalList = []
        signalListDH = []
        signalListPDF = []
        if verbose: print 'top of loop',year,lepton,cat
        
        # ##################################################
        # set up the signal histograms and the mzg ranges #
        # ##################################################

        for prod in sigNameList:
          signalListDS = []
          for mass in massList:
            # store the unbinned signals for CB fitting
            signalTree = signalDict[lepton+year+"_M"+mass].Get(treeName)
            #signalTree.Print()
            sigName = '_'.join(['ds_sig',prod,lepton,year,'cat'+cat,'M'+mass])

            sig_argSW = RooArgSet(mzg,weight)
            sig_ds    = RooDataSet(sigName,sigName,sig_argSW,'Weight')

            lumiWeight = LumiXSWeighter(int(mass), lepton)
            for i in signalTree:
              if cat=="EB" and fabs(i.ph_eta)>EBetaCut: continue
              elif cat=="EE" and fabs(i.ph_eta)<EBetaCut: continue
              #print "loop in signal tree", cat, category[cat], tmpType, "mass =", tmpSigMass[0] 

              if i.m_llg> 90 and i.m_llg<190:
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
            if prod is 'gg':
              if verbose: print 'signal mass loop', mass
              histName  = '_'.join(['sig',lepton,year,'cat'+cat,'M'+mass])
              rangeName = '_'.join(['range',lepton,year,'cat'+cat,'M'+mass])

              signalList.append(TH1F(histName, histName, 100, 90, 190))
              signalList[-1].SetLineColor(kRed)
              signalTree = signalDict[lepton+year+"_M"+mass].Get(treeName)
              
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

              # range is +/- 1 RMS centered around signal peak
              rangeLow = signalList[-1].GetMean()-1.0*signalList[-1].GetRMS()
              rangeHi  = signalList[-1].GetMean()+1.0*signalList[-1].GetRMS()
              mzg.setRange(rangeName,rangeLow,rangeHi)

              mzg_argL = RooArgList(mzg)
              mzg_argS = RooArgSet(mzg)
              signalListDH.append(RooDataHist('dh_'+histName, 'dh_' +histName,mzg_argL,signalList[-1]))
              signalListPDF.append(RooHistPdf('pdf_'+histName,'pdf_'+histName,mzg_argS,signalListDH[-1],2))
              getattr(ws,'import')(signalListPDF[-1])
              if verbose: print '\n\n ** finshed one mass -->>', mass

            if debugPlots and prod is 'gg':
              testFrame = mzg.frame()
              for i,signal in enumerate(signalListPDF):
                signalListDH[i].plotOn(testFrame)
                signal.plotOn(testFrame)
              testFrame.Draw()
              c.SaveAs(plotBase+'_'.join(['signals',prod,year,lepton,'cat'+cat,'M'+mass])+'.png')
              
            if debugPlots:
              testFrame = mzg.frame()
              for signal in signalListDS:
                signal.plotOn(testFrame, RooFit.DrawOption('pl'))
              testFrame.Draw()
              c.SaveAs(plotBase+'_'.join(['ds','sig',prod,year,lepton,'cat'+cat,'M'+mass])+'.png')
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
          
          if i.m_llg> 90 and i.m_llg<190:
            mzg.setVal(i.m_llg)
            data_ds.add(data_argS)

        #dataTree.ResetBranchAddresses()

        if verbose:
          print dataName
          data_ds.Print()
          print
        if debugPlots:
          testFrame = mzg.frame()
          data_ds.plotOn(testFrame,RooFit.Binning(50))
          testFrame.Draw()
          c.SaveAs(plotBase+'_'.join(['data',year,lepton,'cat'+cat,'M'+mass])+'.png')
          
        getattr(ws,'import')(data_ds)



        # ############
        # make fits #
        # ############
        if verbose: 'starting fits'
        
        if cat is not '5':
          GaussExp = BuildGaussExp(year, lepton, cat, mzg)
          GaussPow = BuildGaussPow(year, lepton, cat, mzg)
          SechExp  = BuildSechExp(year, lepton, cat, mzg)
          SechPow  = BuildSechPow(year, lepton, cat, mzg)
          #GaussBern3  = BuildGaussStepBern3(year, lepton, cat, mzg)
          GaussBern4 = BuildGaussStepBern4(year, lepton, cat, mzg)
          #GaussBern5 = BuildGaussStepBern5(year, lepton, cat, mzg)
          #GaussBern6 = BuildGaussStepBern6(year, lepton, cat, mzg)
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
          #Bern5 = BuildBern5(year, lepton, cat, mzg)
          #BB = BuildBetaAndBern(year, lepton, cat, mzg, 'DalitzRegion')
          GB = BuildGaussAndBern(year, lepton, cat, mzg, 'DalitzRegion')

          if verbose:
            GaussExp.Print()
            GaussPow.Print()
            SechExp.Print()
            SechPow.Print()
            #GaussBern3.Print()
            GaussBern4.Print()
            #GaussBern5.Print()
            #GaussBern6.Print()
            #SechBern3.Print()
            #SechBern4.Print()
            #SechBern5.Print()
            #BB.Print()
            GB.Print()
            
          GaussExp.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          GaussPow.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          SechExp.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          SechPow.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          #GaussBern3.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          GaussBern4.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          #GaussBern5.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          #GaussBern4.fitTo(data_ds,RooFit.Range('DalitzRegion'), RooFit.Strategy(1))
          #GaussBern5.fitTo(data_ds,RooFit.Range('DalitzRegion'), RooFit.Strategy(1))
          #GaussBern6.fitTo(data_ds,RooFit.Range('DalitzRegion'), RooFit.Strategy(1))
          #GaussBern6.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          #SechBern3.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          #SechBern4.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          #SechBern5.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          #gauss.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          #BetaFunc.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          #Kumaraswamy.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          #Bern5.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          #BB.fitTo(data_ds,RooFit.Range('DalitzRegion'))
          GB.fitTo(data_ds,RooFit.Range('DalitzRegion'))

          if debugPlots:
            leg  = TLegend(0.7,0.7,1.0,1.0)
            leg.SetFillColor(0)
            leg.SetShadowColor(0)
            leg.SetBorderSize(1)
            leg.SetHeader('_'.join(['test','fits',year,lepton,'cat'+cat]))
            testFrame = mzg.frame(90,190)
            data_ds.plotOn(testFrame,RooFit.Binning(50))
            GaussExp.plotOn(testFrame,RooFit.Name('GaussExp'))
            GaussPow.plotOn(testFrame,RooFit.LineColor(kCyan),  RooFit.Name('GaussPow'))
            SechExp.plotOn(testFrame,RooFit.LineColor(kRed),    RooFit.Name('SechExp'))
            SechPow.plotOn(testFrame,RooFit.LineColor(kOrange), RooFit.Name('SechPow'))
            #GaussBern3.plotOn(testFrame,RooFit.LineColor(kViolet))
            GaussBern4.plotOn(testFrame,RooFit.LineColor(kPink), RooFit.Name('GaussBern4'))
            #GaussBern5.plotOn(testFrame,RooFit.LineColor(kGray))
            #GaussBern6.plotOn(testFrame,RooFit.LineColor(kGreen+2), RooFit.Name('GaussBern6'))
            #SechBern3.plotOn(testFrame,RooFit.LineColor(kMagenta))
            #SechBern4.plotOn(testFrame,RooFit.LineColor(kBlack))
            #SechBern5.plotOn(testFrame,RooFit.LineColor(kGreen))
            #gauss.plotOn(testFrame,RooFit.LineColor(kBlue), RooFit.Name('Gauss'))
            #BetaFunc.plotOn(testFrame,RooFit.LineColor(kBlack), RooFit.Name('Beta'))
            #Kumaraswamy.plotOn(testFrame,RooFit.LineColor(kCyan), RooFit.Name('Kumaraswamy'))
            #Bern5.plotOn(testFrame,RooFit.LineColor(kRed), RooFit.Name('Bern5'))
            #BB.plotOn(testFrame,RooFit.LineColor(kViolet), RooFit.Name('Beta+Bern4'))
            GB.plotOn(testFrame,RooFit.LineColor(kGreen), RooFit.Name('GaussBern3'))
            testFrame.Draw()
            #leg.AddEntry(testFrame.findObject('Beta'),'Beta','l')
            #leg.AddEntry(testFrame.findObject('Kumaraswamy'),'Kumaraswamy','l')
            ##leg.AddEntry(testFrame.findObject('Bern5'),'Bern5','l')
            #leg.AddEntry(testFrame.findObject('Beta+Bern4'),'Beta+Bern4','l')
            leg.AddEntry(testFrame.findObject('GaussExp'),'GaussExp','l')
            leg.AddEntry(testFrame.findObject('GaussPow'),'GaussPow','l')
            leg.AddEntry(testFrame.findObject('SechExp'), 'SechExp', 'l')
            leg.AddEntry(testFrame.findObject('SechPow'), 'SechPow', 'l')
            leg.AddEntry(testFrame.findObject('GaussBern3'),'GaussBern3','l')
            leg.AddEntry(testFrame.findObject('GaussBern4'),'GaussBern4','l')
            #leg.AddEntry(testFrame.findObject('GaussBern6'),'GaussBern6','l')
            leg.Draw()
            c.Print(plotBase+'_'.join(['fits',year,lepton,'cat'+cat])+'.png')
 

          #raw_input()
          getattr(ws,'import')(GaussExp)
          getattr(ws,'import')(GB)
          getattr(ws,'import')(GaussPow)
          getattr(ws,'import')(SechExp)
          getattr(ws,'import')(SechPow)
          #getattr(ws,'import')(GaussBern3)
          getattr(ws,'import')(GaussBern4)
          #getattr(ws,'import')(GaussBern5)
          #getattr(ws,'import')(GaussBern6)
          #getattr(ws,'import')(SechBern3)
          #getattr(ws,'import')(SechBern4)
          #getattr(ws,'import')(SechBern5)

        else:
          print cat
          

        ws.commitTransaction()

  u.createDir(subdir)
  ws.writeToFile(subdir+'/testRooFitOut_Dalitz.root')

  ws.Print()

  print '\n \t we did it!\t'



if __name__=="__main__":

  print len(sys.argv)
  print sys.argv
  if len(sys.argv) != 3:
    sys.exit()
    
  s = sys.argv[1]
  print s
  doInitialFits(s)
  print "done"
