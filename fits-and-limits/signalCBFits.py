#!/usr/bin/env python
import sys
from ROOT import *
from rooFitBuilder import *
gROOT.SetBatch()
#gSystem.SetIncludePath( "-I$ROOFITSYS/include/" );
sys.path.append("../zgamma")
import utils as u
from initialFitProducer import AutoVivification
import ConfigParser as cp
cf = cp.ConfigParser()
cf.read('config.cfg')

gROOT.ProcessLine(".L ~/tdrstyle.C")
setTDRStyle()

# rounding function for interpolation
def roundTo5(x, base=5):
  return int(base * round(float(x)/base))

def SignalNameParamFixer(year,lepton,cat,sig,mass,ws):
  fitName    = '_'.join(['CBG',year,lepton,'cat'+cat,sig,mass,'Interp'])
  newFitName = '_'.join(['sig',sig,lepton,year,'cat'+cat])
  mean       = '_'.join(['meanCBG',   year,lepton,'cat'+cat,sig,mass, 'Interp'])
  sigmaCB    = '_'.join(['sigmaCBCBG',year,lepton,'cat'+cat,sig,mass, 'Interp'])
  sigmaG     = '_'.join(['sigmaGCBG', year,lepton,'cat'+cat,sig,mass, 'Interp'])
  meanNew    = '_'.join(['sig',sig,'mean',   lepton,year,'cat'+cat])
  sigmaCBNew = '_'.join(['sig',sig,'sigmaCB',lepton,year,'cat'+cat])
  sigmaGNew  = '_'.join(['sig',sig,'sigmaG', lepton,year,'cat'+cat])
  mShift     = '_'.join(['sig',sig,'mShift', lepton,year,'cat'+cat])
  sigmaShift = '_'.join(['sig',sig,'sigmaShift',lepton,year,'cat'+cat])
  ws.factory(mShift+'[1]')
  ws.factory(sigmaShift+'[1]')
  ws.factory('prod::'+meanNew+'('+mean+','+mShift+')')
  ws.factory('prod::'+sigmaCBNew+'('+sigmaCB+','+sigmaShift+')')
  ws.factory('prod::'+sigmaGNew+'('+sigmaG+','+sigmaShift+')')
  ws.factory('EDIT::'+newFitName+'('+fitName+','+mean+'='+meanNew+','+sigmaCB+'='+sigmaCBNew+','+sigmaG+'='+sigmaGNew+')')
  
  
def SignalFitMaker(lep, year, cat, subdir):
  #massList        = [str(a)+".0" for a in xrange(120,123)]
  massList   = ['%.1f'%(a) for a in u.drange(120,150,0.5)]
    
  #print massList
  #raw_input("Input raws  ")
      
  # reasd these from a config file:
  #massList        = [a.strip() for a in (cf.get("fits","massList-more")).split(',')]
  sigNameList     = [a.strip() for a in (cf.get("fits","sigNameList")).split(',')]

  #print massList
  #raw_input()
  
  plotBase = cf.get("path","htmlbase")+"/html/zgamma/dalitz/fits-"+subdir
  u.createDir(plotBase)
  rooWsFile = TFile(subdir+'/testRooFitOut_Dalitz.root')
  myWs = rooWsFile.Get('ws')
  #myWs.Print()

  RooRandom.randomGenerator().SetSeed(8675309)

  c = TCanvas("c","c",0,0,500,400)
  c.cd()
  mzg = myWs.var('CMS_hzg_mass')
  cardDict = AutoVivification()
  for mass in massList:
    cardDict[lep][year][cat][mass] = RooWorkspace('ws_card')

    print cardDict[lep][year][cat][mass]
    
    # we need crystal ball + gaussian fits for all mass points, and for all production methods.
    # we also need to interpolate the distributions for 0.5 mass bins, so we use some tricks
    # in order to create new fits out of the existing 5GeV steps
    
    # ##################
    # Start the Loop! #
    # ##################

  for prod in sigNameList:
    dsList   = []
    fitList  = []
    normList = []
    oldMassHi = oldMassLow = 0
    for massString in massList:

      # ############################################
      # Get the high and low mass references      #
      # If the point does not need interpolation, #
      # we just use it for high and low           #
      # ############################################

      mass = float(massString)
      if mass%5.0 == 0.0:
        massHi  = int(mass)
        massLow = int(mass)
      else:
        massRound = roundTo5(massString)
        if mass<massRound:
          massHi  = massRound
          massLow = massRound-5
        else:
          massHi  = massRound+5
          massLow = massRound
      ###### only calc the low and high points if they change
      if not(oldMassLow == massLow and oldMassHi == massHi):
        oldMassLow = massLow
        oldMassHi  = massHi

        ###### fit the low mass point
        #if massLow<=125:
        if massLow<=125:
          mzg.setRange('fitRegion1',110,int(massLow)+10)
        else:
          mzg.setRange('fitRegion1',int(massLow)-15,int(massLow)+10)
        sigNameLow = '_'.join(['ds','sig',prod,lep,year,'cat'+cat,'M'+str(massLow)])
        sig_ds_Low = myWs.data(sigNameLow)
        if massLow == massHi: dsList.append(sig_ds_Low)

        #print massLow
        #raw_input("Mass low")
        CBG_Low = BuildCrystalBallGauss(year,lep,cat,prod,str(massLow),'Low',mzg,mean = massLow)[0]

        CBG_Low.fitTo(sig_ds_Low, RooFit.Range('fitRegion1'), RooFit.SumW2Error(kTRUE), RooFit.Strategy(1), RooFit.NumCPU(4), RooFit.PrintLevel(-1))

        ###### fit the hi mass point
        #if massHi<=125:
        if massHi<=100:
          mzg.setRange('fitRegion2',115,int(massHi)+10)
        else:
          mzg.setRange('fitRegion2',int(massHi)-15,int(massHi)+10)
        sigNameHi = '_'.join(['ds','sig',prod,lep,year,'cat'+cat,'M'+str(massHi)])
        sig_ds_Hi = myWs.data(sigNameHi)

        CBG_Hi = BuildCrystalBallGauss(year,lep,cat,prod,str(massHi),'Hi',mzg,mean = massHi)[0]

        CBG_Hi.fitTo(sig_ds_Hi, RooFit.Range('fitRegion2'), RooFit.SumW2Error(kTRUE), RooFit.Strategy(1), RooFit.NumCPU(4), RooFit.PrintLevel(-1))

      ###### interpolate the two mass points
      massDiff = (massHi - mass)/5.
      #if mass<=125:
      if mass<=100:
        mzg.setRange('fitRegion_'+massString,115,mass+10)
      else:
        mzg.setRange('fitRegion_'+massString,mass-15,mass+15)
      beta = RooRealVar('beta','beta', 0.5, 0., 1.)
      if massHi == massLow:
        beta.setVal(1);
      else:
        beta.setVal(massDiff)

      interp_pdf = RooIntegralMorph('interp_pdf', 'interp_pdf', CBG_Low, CBG_Hi, mzg, beta)
      interp_ds  = interp_pdf.generate(RooArgSet(mzg), 10000)
      normList.append(sig_ds_Low.sumEntries()*massDiff+sig_ds_Hi.sumEntries()*(1-massDiff))
      yieldName = '_'.join(['sig',prod,'yield',lep,year,'cat'+cat])
      yieldVar = RooRealVar(yieldName,yieldName,sig_ds_Low.sumEntries()*massDiff+sig_ds_Hi.sumEntries()*(1-massDiff))


      sigNameInterp = '_'.join(['ds','sig',prod,lep,year,'cat'+cat,'M'+str(mass)])

      CBG_Interp,paramList = BuildCrystalBallGauss(year,lep,cat,prod,str(mass),'Interp',mzg,mean = mass)

      CBG_Interp.fitTo(interp_ds, RooFit.Range('fitRegion_'+massString), RooFit.SumW2Error(kTRUE), RooFit.Strategy(1), RooFit.NumCPU(4), RooFit.PrintLevel(-1))

      #print paramList
      #raw_input("Param List")
      
      for param in paramList:
        param.setConstant(True)
      fitList.append(CBG_Interp)

      print lep,year,cat,mass, cardDict[lep][year][cat][str(mass)]

      getattr(cardDict[lep][year][cat][str(mass)],'import')(CBG_Interp)
      getattr(cardDict[lep][year][cat][str(mass)],'import')(yieldVar)
      cardDict[lep][year][cat][str(mass)].commitTransaction()

      testFrame = mzg.frame(float(massList[0])-10,float(massList[-1])+10)
    for signal in dsList:
      signal.plotOn(testFrame)
      # print signal.mean(mzg)
      # raw_input("Mean of the fit ")

    for i,fit in enumerate(fitList):
      regionName = fit.GetName().split('_')[-1]
      # fit.plotOn(testFrame)
      # fit.plotOn(testFrame, RooFit.NormRange('fitRegion_'+regionName))
      fit.plotOn(testFrame, RooFit.Normalization(normList[i],RooAbsReal.NumEvent))
      fit.paramOn(testFrame)
      # testFrame.getAttText().SetTextSize(0.02)
      # testFrame.getAttText().SetTextColor(kRed)
      # fit.statOn(testFrame,RooFit.Layout(0.18,0.43,0.87))

    testFrame.SetMaximum(0.7)
    testFrame.Draw()
    testFrame.SetTitle(";m_{H} (GeV);fit pdf")
    c.SaveAs(plotBase+"/"+'_'.join(['sig','fit',prod,lep,year,'cat'+cat])+'.png')

  for prod in sigNameList:
    for mass in massList:
      SignalNameParamFixer(year,lep,cat,prod,mass,cardDict[lep][year][cat][mass])

      #print cardDict[lep][year][cat][mass]
      #raw_input("RAW INPUT")

    file = TFile("parampampam.root",'recreate')
    file.cd()
    
    #mass="121.0"
    #mean  = RooRealVar("mean","mean", 130,120,140) 
    #sigma = RooRealVar("sigma","sigma", 5, 0.1,30.0) 
    #fff = RooGaussian("fff","gauss model",mzg, mean, sigma) 

    #fff = (cardDict[lep][year][cat][mass]).pdf('_'.join(['sig',prod,lep,year,'cat'+cat]))    
    #fff.fitTo(sig_ds_Low, RooFit.Range(110, 150))
    #fff.fitTo(sig_ds_Low, RooFit.Range(110, 130), RooFit.SumW2Error(kTRUE), RooFit.Strategy(1), RooFit.NumCPU(4), RooFit.PrintLevel(-1))
    #testFrame = mzg.frame(110,150)
    #dsList[0].plotOn(testFrame)
    #fff.plotOn(testFrame)
    #fff.paramOn(testFrame, RooFit.Layout(0.18,0.43,0.67))
    #testFrame.Draw()
    ##c.SaveAs(plotBase+"/"+'_'.join(['refit_sig','fit',prod,lep,year,'cat'+cat])+'.png')
    #c.Write()  
        
      
  for mass in massList:
    fileName = subdir+'/'+'_'.join(['SignalOutput',lep,year,'cat'+cat,mass])
    cardDict[lep][year][cat][mass].writeToFile(fileName+'.root')
    #cardDict[lep][year][cat][mass].writeToFile('testCards/'+fileName+'.root')


#signal = myWs.data('ds_sig_gg_el_2012_cat4_M125')
#CBG = BuildCrystalBallGauss('2012','el','4','gg','125',mzg)
#CBG.fitTo(signal, RooFit.Range('fitRegion1'), RooFit.SumW2Error(kTRUE), RooFit.PrintLevel(-1))
#CBG.fitTo(signal, RooFit.SumW2Error(kTRUE), RooFit.PrintLevel(-1))

if __name__=="__main__":
  print len(sys.argv), sys.argv

  yearList        = [a.strip() for a in (cf.get("fits","yearList")).split(',')]
  leptonList      = [a.strip() for a in (cf.get("fits","leptonList")).split(',')]
  catList         = [a.strip() for a in (cf.get("fits","catList")).split(',')]

  ver = cf.get("fits","ver")

  for y in yearList:
    for c in catList:
      for l in leptonList:
        SignalFitMaker(l,y,c, ver)

  print "\n \t\t All finished! \n"
