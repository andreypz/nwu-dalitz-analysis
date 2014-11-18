#!/usr/bin/env python
import sys, os
import numpy as np
#import pdb
from rooFitBuilder import *
from ROOT import *
gROOT.SetBatch()
sys.path.append("../zgamma")
gSystem.SetIncludePath( "-I$ROOFITSYS/include/" );
gROOT.ProcessLine(".L ../tdrstyle.C")
TH1.SetDefaultSumw2(kTRUE)
gStyle.SetOptTitle(0)
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

verbose    = 1
#massList   = ['125']
massList      = [a.strip()[0:3] for a in (cf.get("fits","massList")).split(',')]
yearList      = [a.strip() for a in (cf.get("fits","yearList")).split(',')]
leptonList    = [a.strip() for a in (cf.get("fits","leptonList")).split(',')]
catList       = [a.strip() for a in (cf.get("fits","catList")).split(',')]
sigNameList   = [a.strip() for a in (cf.get("fits","sigNameList")).split(',')]
mllBins = u.mllBins()

EBetaCut = 1.4442
EEetaCut = 1.566
ptMllgCut  = 0.3
ptGammaCut = 0
ptDiLepCut = 0
ptSumLep   = 0

lowCutOff  = 110
highCutOff = 170

doBlind    = int(cf.get("fits","blind"))
hjp=0
if 'hjp' in sigNameList: hjp = 1

colors = {}
for f,col in cf.items("colors"):
  colors[f] = eval(col)


debugPlots = 0
noSFweight = 0
rootrace   = 0
if rootrace: RooTrace.active(kTRUE)
# OK listen, we're gunna need to do some bullshit just to get a uniform RooRealVar name for our data objects.
# So we'll loop through the Branch, set mzg to the event value (if it's in range), and add that to our RooDataSet.
# This way, we make a RooDataSet that uses our 'CMS_hzg_mass' variable.

def LumiXSWeighter(mH, prod, sel, Nev=None):
  if prod == 'hjp':
    cro = u.getCS("HtoJPsiGamma")
  else:
    cro = float(u.conf.get(prod+"H-"+str(mH), "cs-"+sel))

  sc = float(u.lumi*cro)/Nev
  print 'M=',mH, 'cs=',cro, 'Nev=', Nev, 'lumi=',u.lumi, "scale", sc
  return sc


def LoopOverTree(myTree, cat, mzg, args, ds, lumiWeight):
  for i in myTree:
    if cat=='mll50':
      if i.m12<20 or i.m12>50: continue
      if fabs(i.eta3)>EBetaCut: continue  # only EB in 20 to 50 region
    elif cat in ['m1','m2','m3','m4','m5','m6','m7']:
      if fabs(i.eta3)>EBetaCut: continue
      if i.m12<mllBins[int(cat[1])-1][0] or i.m12>mllBins[int(cat[1])][0]: continue
    else:
      if i.m12>20: continue
      if i.pt1+i.pt2 < ptSumLep: continue
    if   cat=="EB" and fabs(i.eta3)>EBetaCut: continue
    elif cat=="EE" and fabs(i.eta3)<EEetaCut: continue

    if i.pt3/i.m123 < ptMllgCut or i.pt12/i.m123 < ptMllgCut: continue
    if i.pt3 < ptGammaCut or i.pt12 < ptDiLepCut: continue

    if i.dr13<1 or i.dr23<1: continue

    if hjp: #select j/psi
      if i.m12<2.9 or i.m12>3.3: continue
    else: #veto j/psi and upsions;
      if i.m12>2.9 and i.m12<3.3: continue
      if i.m12>9.3 and i.m12<9.7: continue

    if i.m123> lowCutOff and i.m123<highCutOff:
      mzg.setVal(i.m123)
      if lumiWeight!=None:
        Weight = lumiWeight*i.weight
        ds.add(args, Weight)
      else:
        ds.add(args)
  # sig_argSW.Print()


def doInitialFits(subdir):
  print '\t Loading up the files'

  plotBase1 = cf.get("path","htmlbase")+'/html/zgamma/dalitz/fits-'+subdir+'/init'
  basePath  = cf.get("path","base")+'/batch_output/zgamma/8TeV/'+subdir+'/'
  if not os.path.exists(basePath):
    print basePath1, "does not exist!"
    sys.exit(0)

  print "\t ==== All the input files are initialized into these dictionaries ===="""
  dataDict = {}
  signalDict = {}
  if hjp:
    signalDict = {'hjp_mu2012_M125':TFile(basePath+'jp-mugamma_2012/hhhh_HiggsToJPsiGamma_1.root','r')}
    dataDict['mu2012'] = TFile(basePath+'m_Data_jp-mugamma_2012.root','r')
  else:
    for y in yearList:
      for l in leptonList:
        if l == 'mu': tag = 'mugamma'
        if l == 'el': tag = 'elgamma'
        dataDict[l+y] = TFile(basePath+'m_Data_'+tag+'_'+y+'.root','r')
        for s in sigNameList:
          for m in massList:
            signalDict[s+'_'+l+y+'_M'+m] = TFile(basePath+tag+'_'+y+'/hhhh_'+s+'H-mad'+m+'_1.root','r')
  print "\t ===  All files are set === \n"

  treeName = 'apzTree/apzTree'

  binning = 30
  binWidth = 2

  weight  = RooRealVar('Weight','Weight',0,100)
  mzg  = RooRealVar('CMS_hzg_mass','CMS_hzg_mass', lowCutOff,highCutOff)
  mzg.setRange('FullRegion',   lowCutOff, highCutOff)
  mzg.setRange('DalitzRegion', lowCutOff, highCutOff)
  mzg.setRange('SignalRegion', 122, 128)
  mzg.setBins(50000,'cache')
  mzg.setRange('r1',110,120)
  mzg.setRange('r2',130,170)

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
        plotBase = plotBase1+'_'.join([year,lepton,'cat'+cat])+'/'
        u.createDir(plotBase)
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

            print signalTree, "A signal tree", prod, '  categ=',cat, "mass =", mass
            Nev = u.getTotalEvents(signalDict[prod+"_"+lepton+year+"_M"+mass])
            lumiWeight = LumiXSWeighter(int(mass), prod, lepton, Nev)
            print 'Nev = ', Nev
            LoopOverTree(signalTree, cat, mzg, sig_argSW, sig_ds, lumiWeight)

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


            if verbose:
              print 'signal mass loop', mass
              print 'in gg   INPUT ', cat, prod,mass
              # raw_input()

              print histName
              signalTree.Print()
              print

            if cat=="0":
              if noSFweight:
                signalTree.Draw('m123>>'+histName,'')
              else:
                signalTree.Draw('m123>>'+histName,'weight')

            elif cat=="EB":
              signalTree.Draw('m123>>'+histName,"weight*(fabs(eta3)<"+str(EBetaCut)+")")
            elif cat=="EE":
              signalTree.Draw('m123>>'+histName,"weight*(fabs(eta3)>"+str(EEetaCut)+")")

            if signalList[-1].Integral()!=0:
              signalList[-1].Scale(1/signalList[-1].Integral())
            signalList[-1].Smooth(2)

            if prod in ['gg','hjp'] :
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


          if debugPlots:
            print '\n\n Now make some plots!\n'
            testFrame = mzg.frame()
            for i,signal in enumerate(signalListPDF):
              sigName = 'sig-'+str(i)
              signalListDH[i].plotOn(testFrame)
              signal.plotOn(testFrame, RooFit.Name(sigName))
              # signal.paramOn(testFrame)
              # signal.statOn(testFrame)
              # print i,signal
              testFrame.SetTitle(';m_{#mu#mu#gamma} (GeV);Events')
              testFrame.SetMinimum(0.0001)
              testFrame.Draw()
            CMS_lumi(c, 2, 11)
            c.SaveAs(plotBase+'_'.join(['signals',prod,year,lepton,'cat'+cat])+'.png')


            # testFrame = mzg.frame()
            for signal in signalListDS:
              signal.plotOn(testFrame, RooFit.DrawOption('pl'))
              testFrame.Draw()
            CMS_lumi(c, 2, 11)
            c.SaveAs(plotBase+'_'.join(['ds_sig',year,lepton,'cat'+cat])+'.png')

        del signalTree


        # ###############
        # get the data #
        # ###############
        if verbose: print 'starting data section'

        dataName = '_'.join(['data',lepton,year,'cat'+cat])
        dataTree = dataDict[lepton+year].Get(treeName)

        data_argS = RooArgSet(mzg)
        data_ds   = RooDataSet(dataName,dataName,data_argS)

        LoopOverTree(dataTree, cat, mzg, data_argS, data_ds, None)

        if verbose:
          print dataName
          data_ds.Print()
          print

        if debugPlots:
          testFrame = mzg.frame()
          data_ds.plotOn(testFrame,RooFit.Binning(binning))

          testFrame.Draw()
          CMS_lumi(c, 2, 11)
          c.SaveAs(plotBase+'_'.join(['data',year,lepton,'cat'+cat])+'.png')

        getattr(ws,'import')(data_ds)



        # ############
        # make fits #
        # ############
        if verbose: '\t\t *** Starting BKG fits ***\n'

        fitBuilder = FitBuilder(mzg, year, lepton, cat)
        bgFitList = ['Exp','Pow','Bern2','Bern3','Bern4','Bern5','Laurent']

        testFrame = mzg.frame(RooFit.Range('DalitzRegion'))
        if doBlind:
          data_ds.plotOn(testFrame, RooFit.Binning(binning),RooFit.Name('data'),RooFit.CutRange('r1'))
          data_ds.plotOn(testFrame, RooFit.Binning(binning),RooFit.Name('data'),RooFit.CutRange('r2'))
        else:
          data_ds.plotOn(testFrame, RooFit.Binning(binning),RooFit.Name('data'))

        leg = TLegend(0.63,0.6,0.91,0.88)
        leg.SetFillColor(0)
        leg.SetShadowColor(0)
        leg.SetBorderSize(1)
        leg.SetHeader(', '.join([year,lepton,'cat: '+cat]))

        for fitName in bgFitList:
          ndof = fitBuilder.FitNdofDict[fitName]
          fit  = fitBuilder.Build(fitName)

          if verbose:
            fit.Print()

          fit.fitTo(data_ds,RooFit.Range('DalitzRegion'), RooFit.Strategy(1))
          fit.plotOn(testFrame, RooFit.LineColor(colors[fitName.lower()]), RooFit.Name(fitName))

          testFrame.Draw()
          chi2 = testFrame.chiSquare(fitName,'data',ndof)

          leg.AddEntry(testFrame.findObject(fitName),fitName+'   #chi2 = {0:.3f}'.format(chi2),'l')
          getattr(ws,'import')(fit)


        if doBlind:
          testFrame.SetMinimum(0.1)
        testFrame.SetTitle(";m_{#mu#mu#gamma} (GeV);Events/"+str(binWidth)+" GeV")
        leg.Draw()
        CMS_lumi(c, 2, 11)
        c.SaveAs(plotBase+'_'.join(['fits',year,lepton,'cat'+cat])+'.png')

        ws.commitTransaction()

        print '\t Commit transaction end \n'

        yi_da0[year][lepton][cat] = data_ds.sumEntries()
        yi_da1[year][lepton][cat] = data_ds.sumEntries('1','SignalRegion')

  u.createDir(subdir)
  ws.writeToFile(subdir+'/testRooFitOut_Dalitz.root')



  print '*** Some yields from data'
  print 'Full range:', yi_da0
  if not doBlind:
    print 'in 122-128:', yi_da1

  if verbose:
    ws.Print()

    print '*** Some yields from ggH'
    print 'Full range:', yi_sig0
    print 'in 122-128:', yi_sig1

    print "\n ** You should also notice that EBeta cut was: ", EBetaCut
    print "And that the range was from ", lowCutOff, 'to', highCutOff

    print "\n Mean from dataset: ", mean_sig
    print "\n Sigma DS:", sigma_sig

    # print "\n Mean from Datahist:", mean_gg0_dh
    # print "\n Sigma DH:", sigma_gg0_dh
    raw_input('Enter')


  if verbose:
    for year in yearList:
      for lepton in leptonList:
        for prod in sigNameList:
          for cat in catList:
            t_mean = []
            t_sigma = []
            for mass in massList:
              l_mean = []
              l_sigma = []
              l_mean.append(mass)
              l_sigma.append(mass)

              print year, lepton, cat, prod
              print mean_sig[year][mass][lepton][cat][prod]
              print sigma_sig[year][mass][lepton][cat][prod]
              # a,b = mean_sig[year][mass][lepton][cat][prod]
              # l.append("%.2f &pm; %.2f"%(a,b))
              l_mean.append("%.3f&pm;%.3f"%(mean_sig[year][mass][lepton][cat][prod][0], mean_sig[year][mass][lepton][cat][prod][1]))
              l_sigma.append("%.3f&pm;%.3f"%(sigma_sig[year][mass][lepton][cat][prod][0],sigma_sig[year][mass][lepton][cat][prod][1]))

              t_mean.append(l_mean)
              t_sigma.append(l_sigma)
            u.makeTable(t_mean,  "mean",  opt="twiki")
            u.makeTable(t_sigma, "sigma", opt="twiki")

  print '\n \t we did it!\t'



if __name__=="__main__":

  print len(sys.argv)
  print sys.argv
  if len(sys.argv) != 2:
    sys.exit()

  s = sys.argv[1]
  if 'vv/' in s: s = s[3:].rstrip('/')
  print 'Params=', s
  cf.set("path","ver", s)
  if '/tthome' in os.getcwd():
    cf.set("path","base", '/tthome/andrey')
    cf.set("path","htmlbase", '/tthome/andrey')
  else:
    cf.set("path","base", '/eos/uscms/store/user/andreypz')
    cf.set("path","htmlbase", '/uscms_data/d2/andreypz/')

  with open(r'config.cfg', 'wb') as configfile:
    cf.write(configfile)

  print '\t \t start \n'
  doInitialFits(s)
  print "\n \t \t done"
