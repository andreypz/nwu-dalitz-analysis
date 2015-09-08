#!/usr/bin/env python
import sys, os
import numpy as np
from rooFitBuilder import *
from ROOT import *
gROOT.SetBatch()
sys.path.append("../zgamma")
import utils as u
gSystem.SetIncludePath( "-I$ROOFITSYS/include/" );
gROOT.ProcessLine(".L ../tdrstyle.C")
TH1.SetDefaultSumw2(kTRUE)
gStyle.SetOptTitle(0)
from optparse import OptionParser
parser = OptionParser(usage="usage: %prog ver [options -v]")
parser.add_option("-v","--verbose", dest="verbose", action="store_true", default=False, help="Verbose mode (print out stuff)")
parser.add_option("--tw", dest="tw",   action="store_true", default=False, help="Use a tree from TW group (in el channel)")
parser.add_option("--dbg", dest="dbg", action="store_true", default=False, help="Debug option: make more debugging plots")

(opt, args) = parser.parse_args()

import ConfigParser as cp
cf = cp.ConfigParser()
cf.read('config.cfg')

verbose    = opt.verbose

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
#ptSumLep   = 44

lowCutOff  = 110
highCutOff = 170
doBlind    = int(cf.get("fits","blind"))
hjp=0
if 'hjp' in sigNameList: hjp = 1

if hjp:
  ptMllgCut  = 0.0
  ptGammaCut = 40.0
  ptDiLepCut = 40.0
  highCutOff = 150
  massList   = ['125']
  leptonList = ['mu']
deltaMllg = highCutOff-lowCutOff


colors = {}
for f,col in cf.items("colors"):
  colors[f] = eval(col)

debugPlots = opt.dbg
noSFweight = 0
rootrace   = 0

if rootrace: RooTrace.active(kTRUE)
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


def LoopOverTreeTW(myTree, cat, mzg, args, ds):
  print 'TW tree:', myTree
  hcount=0
  for i in myTree:
    if i.Meg> lowCutOff and i.Meg<highCutOff:
      mzg.setVal(i.Meg)
      Weight = i.mcwei*i.puwei
      ds.add(args, Weight)

      if Weight==1: #(This is data)
        if i.Meg>120 and i.Meg<130:
          hcount+=1
          # print hcount,i.run, i.event, i.Meg



def LoopOverTree(myTree, cat, mzg, args, ds, lumiWeight):
  print 'Looping my tree, ', myTree
  hcount=0
  for i in myTree:
    if cat=='mll50':
      if i.m12<20 or i.m12>50: continue
      if fabs(i.eta3)>EBetaCut: continue  # only EB in 20 to 50 region
    elif cat in ['m1','m2','m3','m4','m5','m6','m7']:
      if fabs(i.eta3)>EBetaCut: continue
      if i.m12<mllBins[int(cat[1])-1][0] or i.m12>mllBins[int(cat[1])][0]: continue
    else:
      if i.m12>20: continue
      # if i.pt1+i.pt2 < ptSumLep: continue
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
      if hjp and lumiWeight!=None:
        # The H->JPsi+g sample was produced with mH=125.8 GeV.
        # This is an offset to bring it to 125 GeV
        # Only applied to MC (lumiWeight=None)
        mzg.setVal(i.m123-0.8)
      else:
        mzg.setVal(i.m123)

      Weight = 1
      if lumiWeight!=None:
        Weight = lumiWeight*i.weight

      ds.add(args, Weight)

      if lumiWeight==None: #(This is data)
        if i.m123>120 and i.m123<130:
          hcount+=1
          # print hcount,i.run, i.event, i.m123


  # sig_argSW.Print()


def doInitialFits(subdir):
  print '\t Loading up the files'
  u.createDir(subdir)

  basePath  = cf.get("path","base")+'/batch_output/zgamma/8TeV/'
  plotBase1 = cf.get("path","htmlbase")+'/html/zgamma/dalitz/fits-'+subdir+'/init'

  print "\t ==== All the input files are initialized into these dictionaries ===="""
  dataDict   = {}
  signalDict = {}
  if hjp:
    signalDict = {'hjp_mu2012_M125':TFile(basePath+cf.get("path","ver-jp")+'/jp-mugamma_2012/hhhh_HiggsToJPsiGamma_1.root','r')}
    dataDict['mu2012'] = TFile(basePath+cf.get("path","ver-jp")+'/m_Data_jp-mugamma_2012.root','r')
  else:
    for y in yearList:
      for l in leptonList:
        if l == 'mu': tag = 'mugamma'
        if l == 'el': tag = 'elgamma'
        if l == 'ee': tag = 'eegamma'

        if not os.path.exists(basePath):
          print basePath, "does not exist!"
          sys.exit(0)

        if opt.tw and l=='el':
          dataDict[l+y] = TFile('dalitzTW4/data.root','r')
        else:
          dataDict[l+y] = TFile(basePath+cf.get("path","ver-"+l)+'/m_Data_'+tag+'_'+y+'.root','r')

        for s in sigNameList:
          for m in massList:
            if l=='el' and s=='v': continue
            if l=='el' and opt.tw:
              if s=='vbf':
                signalDict[s+'_'+l+y+'_M'+m] = TFile('dalitzTW/output_job_hzg_eeg_dalitz_VBFH_'+m+'.root','r')
              else:
                signalDict[s+'_'+l+y+'_M'+m] = TFile('dalitzTW/output_job_hzg_eeg_dalitz_'+s+'H_'+m+'.root','r')
            else:
              signalDict[s+'_'+l+y+'_M'+m] = TFile(basePath+cf.get("path","ver-"+l)+'/'+tag+'_'+y+'/hhhh_'+s+'H-mad'+m+'_1.root','r')


  print "\t ===  All files are set === \n"

  apzTreeName = 'apzTree/apzTree'
  twTreeName  = 't'

  binWidth = 0.2
  #binning = 30
  binning = int(float(deltaMllg)/binWidth)

  weight  = RooRealVar('Weight','Weight',0,100)
  mzg  = RooRealVar('CMS_hzg_mass','CMS_hzg_mass', lowCutOff,highCutOff)
  mzg.setRange('FullRegion',   lowCutOff, highCutOff)
  mzg.setRange('DalitzRegion', lowCutOff, highCutOff)
  mzg.setRange('SignalRegion', 120, 130)
  mzg.setRange('myFitRange', 115, 135)
  mzg.setBins(50000,'cache')
  mzg.setRange('r1',110,120)
  mzg.setRange('r2',130,170)

  #res  = RooRealVar('Mass_res','Mss_res', 0.8,1.2)

  c = TCanvas("c","c",0,0,500,400)
  c.cd()

  ws = RooWorkspace("ws")

  shist = u.AutoVivification()
  yi_da0  = u.AutoVivification()
  yi_da1  = u.AutoVivification()
  yi_sig0 = u.AutoVivification()
  yi_sig1 = u.AutoVivification()
  mean_sig  = u.AutoVivification()
  sigma_sig = u.AutoVivification()
  FWHM_sig = u.AutoVivification()
  # ###################################
  # start loop over all year/lep/cat #
  # ###################################
  for year in yearList:
    for lepton in leptonList:
      if hjp:
        fs125 = TFile(subdir+'/s125-hjp.root','recreate')
      else:
        fs125 = TFile(subdir+'/s125-'+lepton+'.root','recreate')

      if opt.tw and lepton=='el':
        treeName = twTreeName
      else:
        treeName = apzTreeName

      for cat in catList:
        if lepton=='el' and cat!='EB': continue
        plotBase = plotBase1+'_'.join([year,lepton,'cat'+cat])+'/'
        u.createDir(plotBase)
        if rootrace:
          print "rootrace"
          RooTrace.dump()
          raw_input()

        if verbose: print 'top of loop',year,lepton,cat

        # ##################################################
        # set up the signal histograms and the mzg ranges #
        # ##################################################

        for prod in sigNameList:
          if lepton=='el' and prod=='v': continue
          signalList    = []
          signalListDS = []
          signalListDH  = []
          signalListPDF = []

          for mass in massList:
            histName  = '_'.join(['sig',prod,lepton,year,'cat'+cat,'M'+mass])

            shist[year][lepton][cat][prod][mass] = mzg.createHistogram(histName,
                                                                       RooFit.Binning(binning))

            # store the unbinned signals for CB fitting
            signalTree = signalDict[prod+"_"+lepton+year+"_M"+mass].Get(treeName)

            #signalTree.Print()
            sigName = '_'.join(['ds_sig',prod,lepton,year,'cat'+cat,'M'+mass])

            sig_argSW = RooArgSet(mzg,weight)
            sig_ds    = RooDataSet(sigName,sigName,sig_argSW,'Weight')

            print signalTree, "A signal tree", prod, '  categ=',cat, "mass =", mass

            if opt.tw and lepton=='el':
              LoopOverTreeTW(signalTree, cat, mzg, sig_argSW, sig_ds)

            else:
              Nev = u.getTotalEvents(signalDict[prod+"_"+lepton+year+"_M"+mass])
              sel = lepton
              if sel=='ee': sel ='el'
              lumiWeight = LumiXSWeighter(int(mass), prod, sel, Nev)
              print 'Nev = ', Nev

              LoopOverTree(signalTree, cat, mzg, sig_argSW, sig_ds, lumiWeight)

            signalListDS.append(sig_ds)
            getattr(ws,'import')(signalListDS[-1])

            #raw_input()
            sig_ds.fillHistogram(shist[year][lepton][cat][prod][mass], RooArgList(mzg))
            if float(mass)==125:
              fs125.cd()
              shist[year][lepton][cat][prod][mass].Write()


            signalList.append(shist[year][lepton][cat][prod][mass])

            '''
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
            '''

            signalList[-1].Smooth(2)

            # do some histogramming for gg signal for bias study
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

            if verbose: print '\n\n\t ** Finished the mass -->>', mass, '\n\n'

            yi_sig0[year][mass][lepton][cat][prod] = sig_ds.sumEntries()
            yi_sig1[year][mass][lepton][cat][prod] = sig_ds.sumEntries('1','SignalRegion')

            reduced = signalListDS[-1].reduce("CMS_hzg_mass>"+str(float(mass)*(1-0.1))+"&&CMS_hzg_mass<"+str(float(mass)*(1+0.1)))

            mean_sig[year][mass][lepton][cat][prod]  = [reduced.mean(mzg),  signalList[-1].GetMeanError()]
            sigma_sig[year][mass][lepton][cat][prod] = [reduced.sigma(mzg), signalList[-1].GetRMSError()]
            #mean_sig[year][mass][lepton][cat][prod]  = [signalListDS[-1].mean(mzg), signalList[-1].GetMeanError()]
            #sigma_sig[year][mass][lepton][cat][prod] = [signalListDS[-1].sigma(mzg), signalList[-1].GetRMSError()]
            #mean_sig[year][mass][lepton][cat][prod]  = [signalList[-1].GetMean(), signalList[-1].GetMeanError()]
            #sigma_sig[year][mass][lepton][cat][prod] = [signalList[-1].GetRMS(), signalList[-1].GetRMSError()]

            bin1 = signalList[-1].FindFirstBinAbove(signalList[-1].GetMaximum()/2)
            bin2 = signalList[-1].FindLastBinAbove(signalList[-1].GetMaximum()/2)
            print bin1, bin2
            FWHM_sig[year][mass][lepton][cat][prod] = (signalList[-1].GetBinCenter(bin2) - signalList[-1].GetBinCenter(bin1))
            # mean_dh[year][mass][lepton][cat]  = [signalListDH[-1].GetMean(), signalListDH[-1].GetMeanError()]
            # sigma_dh[year][mass][lepton][cat] = [signalListDH[-1].GetRMS(), signalListDH[-1].GetRMSError()]


          if debugPlots:
            fitBuilder = FitBuilder(mzg, year, lepton, cat)
            print '\n\n Now make some plots!\n'
            testFrame = mzg.frame()
            SigFits = {'0':None} # Format 'mass': fit-function
            for i,signal in enumerate(signalListPDF):
              signalListDH[i].plotOn(testFrame)
              signal.plotOn(testFrame, RooFit.Name(sigName))
              if lepton=='mu': testFrame.SetTitle(';m_{#mu#mu#gamma} (GeV);Events')
              else:            testFrame.SetTitle(';m_{ee#gamma} (GeV);Events')
              testFrame.SetMinimum(0.0001)
              testFrame.Draw()
            CMS_lumi(c, 2, 11, 'Simulation')
            c.SaveAs(plotBase+'_'.join(['signals',prod,year,lepton,'cat'+cat])+'.png')


            # testFrame = mzg.frame()
            for signal in signalListDS:
              """
              sigName = '_'.join(['DBG-DH',prod,lepton,year,'cat'+cat,'M'+mass])
              SigFits[sigName],paramList = fitBuilder.Build('Gauss', mean = 125, meanLow = 120, meanHigh = 130, sigma = 0.2, sigmaLow = 0.05, sigmaHigh = 1)
              SigFits[sigName].fitTo(signalListDS[i], RooFit.Range('myFitRange'), RooFit.PrintLevel(-1))

              print paramList
              for param in paramList:
                print sigName, param.GetName(), param.getVal()
              raw_input("\n ---> Param List. hit enter to continue")
              SigFits[sigName].plotOn(testFrame)
              SigFits[sigName].paramOn(testFrame)
              """
              testFrame = mzg.frame()
              signal.plotOn(testFrame, RooFit.DrawOption('pl'))
              testFrame.Draw()
              myS125 = shist[year][lepton][cat][prod]['125'].Clone()
              myS125.Draw('same hist')
              myS125.SetLineColor(kGreen+1)

              f1 = TF1("f1", "gaus", 120, 130);
              f1.SetParameters(0,125, 0.2);
              f1.FixParameter(0, 0);
              f1.SetParLimits(1, 122, 127);
              myS125.Fit('f1','R')
              myS125.GetFunction("f1").Print()
              #raw_input("\n ---> Param List. hit enter to continue")

            CMS_lumi(c, 2, 11, 'Simulation')
            c.SaveAs(plotBase+'_'.join(['ds_sig',prod,year,lepton,'cat'+cat])+'.png')

        del signalTree


        # ###############
        # get the data #
        # ###############
        if verbose: print '\n\n \t ** Starting data selection. \n'

        dataName = '_'.join(['data',lepton,year,'cat'+cat])
        dataTree = dataDict[lepton+year].Get(treeName)
        data_argS = RooArgSet(mzg)
        data_argL = RooArgList(mzg)
        data_ds   = RooDataSet(dataName,dataName,data_argS)

        if opt.tw and lepton=='el':
          LoopOverTreeTW(dataTree, cat, mzg, data_argS, data_ds)
        else:
          LoopOverTree(dataTree, cat, mzg, data_argS, data_ds, None)

        # data_ds.Print('all')

        dahist = mzg.createHistogram('hist_'+dataName, RooFit.Binning(binning))
        data_ds.fillHistogram(dahist, RooArgList(mzg))
        data_dh   = RooDataHist('dh_'+dataName, 'dh_' +dataName, data_argL, dahist)
        #dahist.Print('all')
        # data_dh.Print()

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
        #if doBlind:
        #  data_ds.plotOn(testFrame, RooFit.Binning(binning),RooFit.Name('r1'),RooFit.CutRange('r1'))
        #  data_ds.plotOn(testFrame, RooFit.Binning(binning),RooFit.Name('r2'),RooFit.CutRange('r2'))
        #  data_ds.plotOn(testFrame, RooFit.Binning(binning),RooFit.Name('data'),RooFit.Invisible())
        #else:
        #  data_ds.plotOn(testFrame, RooFit.Binning(binning),RooFit.Name('data'))

        data_dh.plotOn(testFrame, RooFit.Binning(binning), RooFit.Name('data2'))

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

          fit.fitTo(data_dh,RooFit.Range('DalitzRegion'), RooFit.Strategy(1))
          #fit.fitTo(data_ds,RooFit.Range('DalitzRegion'), RooFit.Strategy(1))
          fit.plotOn(testFrame, RooFit.LineColor(colors[fitName.lower()]), RooFit.Name(fitName))

          getattr(ws,'import')(fit)

          testFrame.Draw()
          chi2 = testFrame.chiSquare(fitName,'data2',ndof)

          leg.AddEntry(testFrame.findObject(fitName),fitName+(10-len(fitName))*'  '+'#chi2 = {0:.2f}'.format(chi2),'l')


        if doBlind:
          testFrame.SetMinimum(0.1)
        if lepton=='mu':
          testFrame.SetTitle(";m_{#mu#mu#gamma} (GeV);Events/"+str(binWidth)+" GeV")
        else:
          testFrame.SetTitle(";m_{ee#gamma} (GeV);Events/"+str(binWidth)+" GeV")

        testFrame.Draw()

        # dahist.Draw('same hist')

        #shist[year][lepton][cat]['gg']['125'].Scale(10)
        #shist[year][lepton][cat]['gg']['125'].Draw('same hist')
        #shist[year][lepton][cat]['gg']['125'].SetLineColor(kRed+1)

        leg.Draw()
        CMS_lumi(c, 2, 11)
        c.SaveAs(plotBase+'_'.join(['fits',year,lepton,'cat'+cat])+'.png')

        ws.commitTransaction()

        print '\t Commit transaction end \n'

        yi_da0[year][lepton][cat] = data_ds.sumEntries()
        yi_da1[year][lepton][cat] = data_ds.sumEntries('1','SignalRegion')

      fs125.Close()

  ws.writeToFile(subdir+'/testRooFitOut_Dalitz.root')



  print '*** Some yields from data'
  print 'Full range:', yi_da0
  if not doBlind:
    print 'in SignalRegion:', yi_da1

  if verbose:
    ws.Print()

    print '*** Some yields from the Signal'
    print 'Full range:', yi_sig0
    print 'in SignalRegion:', yi_sig1

    print "\n ** You should also notice that EBeta cut was: ", EBetaCut
    print "And that the range was from ", lowCutOff, 'to', highCutOff

    print "\n Mean from dataset: ", mean_sig
    print "\n Sigma DS:", sigma_sig
    print "\n FWHM:", FWHM_sig

    # print "\n Mean from Datahist:", mean_gg0_dh
    # print "\n Sigma DH:", sigma_gg0_dh
    raw_input('Enter')


  if verbose:
    for year in yearList:
      for lepton in leptonList:
        for prod in sigNameList:
          if lepton=='el' and prod=='v': continue
          for cat in catList:
            t_mean = []
            t_sigma = []
            t_fwhm  = []
            for mass in massList:
              l_mean = []
              l_sigma = []
              l_fwhm = []
              l_mean.append(mass)
              l_sigma.append(mass)
              l_fwhm.append(mass)

              print year, lepton, cat, prod
              print mean_sig[year][mass][lepton][cat][prod]
              print sigma_sig[year][mass][lepton][cat][prod]
              print FWHM_sig[year][mass][lepton][cat][prod]
              # a,b = mean_sig[year][mass][lepton][cat][prod]
              # l.append("%.2f &pm; %.2f"%(a,b))
              l_mean.append("%.3f&pm;%.3f" %(mean_sig[year][mass][lepton][cat][prod][0],  mean_sig[year][mass][lepton][cat][prod][1]))
              l_sigma.append("%.3f&pm;%.3f"%(sigma_sig[year][mass][lepton][cat][prod][0],sigma_sig[year][mass][lepton][cat][prod][1]))
              l_fwhm.append("%.3f"% (FWHM_sig[year][mass][lepton][cat][prod]))
              t_mean.append(l_mean)
              t_sigma.append(l_sigma)
              t_fwhm.append(l_fwhm)
            u.makeTable(t_mean,  "mean",  opt="twiki")
            u.makeTable(t_sigma, "sigma", opt="twiki")
            u.makeTable(t_fwhm, "fwhm", opt="twiki")

  print '\n \t We did it!\t'



if __name__=="__main__":

  print len(sys.argv)
  print sys.argv
  if len(sys.argv) < 2:
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
