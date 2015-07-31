#!/usr/bin/env python
import sys
from ROOT import *
from rooFitBuilder import *
gROOT.SetBatch()
gROOT.SetMacroPath(".:../plugins/")
gROOT.LoadMacro("HistManager.cc+")
sys.path.append("../zgamma")
import utils as u
#from apzFitProducer import AutoVivification
import ConfigParser as cp
cf = cp.ConfigParser()
cf.read('config.cfg')

massList0  = ['125']
massList   = ['125.0']
#massList0  = ['120','125','130','135','140','145','150']
#massList   = ['%.1f'%(a) for a in u.drange(120,150,1.0)]

myFunc = 'CBGM' #Crystall-Ball + Gauss (with same Mean)

lumi = u.getLumi()

# rounding function for interpolation
def roundTo5(x, base=5):
  return int(base * round(float(x)/base))

def SignalNameParamFixer(year,lepton,cat,sig,mass,ws, sameMean=True):
  # print 'DBG SignalNameParamFixer'
  fitName    = '_'.join(['CBG',year,lepton,'cat'+cat,sig,mass,'Interp'])
  newFitName = '_'.join(['sig',sig,lepton,year,'cat'+cat])
  mean       = '_'.join(['meanCBG',   year,lepton,'cat'+cat,sig,mass, 'Interp'])
  meanG      = '_'.join(['meanGCBG',  year,lepton,'cat'+cat,sig,mass, 'Interp'])
  meanCB     = '_'.join(['meanCBCBG', year,lepton,'cat'+cat,sig,mass, 'Interp'])
  sigmaCB    = '_'.join(['sigmaCBCBG',year,lepton,'cat'+cat,sig,mass, 'Interp'])
  sigmaG     = '_'.join(['sigmaGCBG', year,lepton,'cat'+cat,sig,mass, 'Interp'])
  meanNew    = '_'.join(['sig',sig,'mean',   lepton,year,'cat'+cat])
  meanGNew   = '_'.join(['sig',sig,'meanG',  lepton,year,'cat'+cat])
  meanCBNew  = '_'.join(['sig',sig,'meanCB', lepton,year,'cat'+cat])
  sigmaCBNew = '_'.join(['sig',sig,'sigmaCB',lepton,year,'cat'+cat])
  sigmaGNew  = '_'.join(['sig',sig,'sigmaG', lepton,year,'cat'+cat])
  mShift     = '_'.join(['sig',sig,'mShift', lepton,year,'cat'+cat])
  sigmaShift = '_'.join(['sig',sig,'sigmaShift',lepton,year,'cat'+cat])
  ws.factory(mShift+'[1]')
  ws.factory(sigmaShift+'[1]')

  if sameMean:
    ws.factory('prod::'+meanNew+'('+mean+','+mShift+')')
    ws.factory('prod::'+sigmaCBNew+'('+sigmaCB+','+sigmaShift+')')
    ws.factory('prod::'+sigmaGNew+'('+sigmaG+','+sigmaShift+')')
    ws.factory('EDIT::'+newFitName+'('+fitName+','+mean+'='+meanNew+','+sigmaCB+'='+sigmaCBNew+','+sigmaG+'='+sigmaGNew+')')
  else:
    ws.factory('prod::'+meanGNew+'('+meanG+','+mShift+')')
    ws.factory('prod::'+meanCBNew+'('+meanCB+','+mShift+')')
    ws.factory('prod::'+sigmaCBNew+'('+sigmaCB+','+sigmaShift+')')
    ws.factory('prod::'+sigmaGNew+'('+sigmaG+','+sigmaShift+')')
    ws.factory('EDIT::'+newFitName+'('+fitName+','+meanG+'='+meanGNew+','+meanCB+'='+meanCBNew+','
               +sigmaCB+'='+sigmaCBNew+','+sigmaG+'='+sigmaGNew+')')


def SignalFitMaker(lep, year, cat, subdir):
  f = TFile('../data/Dalitz_BR50.root','READ')
  g = f.Get('csbr_mu')
  xsBr = g.GetFunction('pol4')
  #g.Print('all')

  hfile = TFile('hists.root','recreate')
  hists = HistManager(hfile)
  #print massList

  u.set_palette()
  # reasd these from a config file:
  sigNameList = [a.strip() for a in (cf.get("fits","sigNameList")).split(',')]
  tev = u.yearToTeV(year)

  if 'hjp' in sigNameList:
    if lep!='mu': return
    massList0  = ['125']
    massList  = ['125.0']

  #print massList
  #raw_input()

  plotBase = cf.get("path","htmlbase")+"/html/zgamma/dalitz/fits-"+subdir+'/init' + '_'.join([year,lep,'cat'+cat])+'/'
  u.createDir(plotBase)

  rooWsFile = TFile(subdir+'/testRooFitOut_Dalitz.root')
  myWs = rooWsFile.Get('ws')
  #myWs.Print()

  RooRandom.randomGenerator().SetSeed(8675309)

  c = TCanvas("c","c",0,0,500,400)
  c.cd()
  mzg = myWs.var('CMS_hzg_mass')
  paraNames = {}
  cardDict = u.AutoVivification()

  mzg.setRange('statRange125',120,130)

  global massList0, massList
  #print massList, massList0
  masses0 = massList0
  masses1 = massList
  if cat in ['m1','m2','m3','m4','m5','m6','m7']:
    masses0 = ['125']
    masses1 = ['125.0']


  for mass in masses1:
    cardDict[lep][year][cat][mass] = RooWorkspace('ws_card')

    print cardDict[lep][year][cat][mass]

  # we need crystal ball + gaussian fits for all mass points, and for all production methods.
  # we also need to interpolate the distributions for 0.5 mass bins, so we use some tricks
  # in order to create new fits out of the existing 5GeV steps

  # ##################
  # Start the Loop! #
  # ##################

  for prod in sigNameList:
    if lep=='el' and prod=='v': continue
    if cat in ['m1','m2','m3','m4','m5','m6','m7']:
      if prod!='gg': continue
    # First of all make the fits to the existing MC samples
    SigFits = {'0':None} # Format 'mass': fit-function

    dsList   = []
    for m in masses0:
      fitBuilder = FitBuilder(mzg, year,lep,cat,sig=prod, mass=m+'.0')
      mzg.setRange('fitRegion1',int(m)-11,int(m)+11)

      sigName = '_'.join(['ds','sig',prod,lep,year,'cat'+cat,'M'+m])
      sig_ds = myWs.data(sigName)
      dsList.append(sig_ds)

      SigFits[m],paramList = fitBuilder.Build(myFunc, piece = 'Interp', mean = float(m))
      SigFits[m].fitTo(sig_ds, RooFit.Range('fitRegion1'), RooFit.SumW2Error(kTRUE),
                       RooFit.Strategy(1), RooFit.NumCPU(4), RooFit.PrintLevel(-1))

      # print paramList
      # raw_input("\n ---> Param List. hit enter to continue")

      paraNames[prod] = []
      for param in paramList:
        param.setConstant(True)
        print param.GetName(), param.getVal()

        if 'sigmaCBCBG' in param.GetName():
          n = '_'.join(['sigmaCBCBG',year,lep,'cat'+cat,prod])
          hists.fill1DHist(param.getVal(), n,";sigmaCB", 50, 0,5, 1, '')
        elif 'fracCBG' in param.GetName():
          n = '_'.join(['fracCBS',year,lep,'cat'+cat,prod])
          hists.fill1DHist(param.getVal(), n,";fracCBG", 50, 0,1, 1, '')
        elif 'mean' in param.GetName():
          n = '_'.join(['mean',year,lep,'cat'+cat,prod])
          hists.fill1DHist(param.getVal(), n,";mean", 160, 115,155, 1, '')
        elif 'sigmaGCBG' in param.GetName():
          n = '_'.join(['sigmaGCBS',year,lep,'cat'+cat,prod])
          hists.fill1DHist(param.getVal(), n,";sigmaGCBG", 50, 0,5, 1, '')
        else:
          n = None
        paraNames[prod].append(n)

      #r.Print()
      #a = r.floatParsFinal()
      #a[0].Print()
      #print a[0].GetName(), a[0].getVal()

      #raw_input('RooFitResults above.')

    fitList  = []
    normList = []
    oldMassHi = oldMassLow = 0
    for massString in masses1:

      fitBuilder = FitBuilder(mzg, year,lep,cat,sig=prod,mass=massString)
      #fitBuilder.__init__(mzg,year,lep,cat,sig=prod,mass=str(massLow))
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
        if not(oldMassLow == massLow and oldMassHi == massHi):
          # only calc the low and high points if they change
          oldMassLow = massLow
          oldMassHi  = massHi


      # print 'massLow and massHi=', massLow, massHi
      # raw_input("Mass low and Hi"+massString+"; hit enter")


      massDiff = (massHi - mass)/5.

      sigName = '_'.join(['ds','sig',prod,lep,year,'cat'+cat,'M'+str(massLow)])
      sig_ds_Low = myWs.data(sigName)
      sigName = '_'.join(['ds','sig',prod,lep,year,'cat'+cat,'M'+str(massHi)])
      sig_ds_Hi = myWs.data(sigName)

      if massHi==massLow:
        # No interpolateion needed, just take the existing fit
        SigFit_Interp = SigFits[str(mass)[:3]]
      else:
        # ##### interpolate the two mass points
        # print 'massdiff=',massDiff
        # raw_input('Massdiff...')

        mzg.setRange('fitRegion_'+massString,mass-11,mass+11)
        #mzg.setRange('fitRegion_'+massString,mass-5,mass+5)
        beta = RooRealVar('beta','beta', 0.5, 0., 1.)
        beta.setVal(massDiff)

        SigFit_Low = SigFits[str(massLow)]
        SigFit_Hi  = SigFits[str(massHi)]

        interp_pdf = RooIntegralMorph('interp_pdf', 'interp_pdf', SigFit_Low, SigFit_Hi, mzg, beta)
        interp_ds  = interp_pdf.generate(RooArgSet(mzg), 10000)

        fitBuilder.__init__(mzg,year,lep,cat,sig=prod,mass=str(mass))
        SigFit_Interp,paramList = fitBuilder.Build(myFunc, piece = 'Interp', mean = mass)

        SigFit_Interp.fitTo(interp_ds, RooFit.Range('fitRegion_'+massString), RooFit.SumW2Error(kTRUE),
                            RooFit.Strategy(1), RooFit.NumCPU(4), RooFit.PrintLevel(-1))


        # r.Print()
        # raw_input('RooFitResults above.')

        for param in paramList:
          param.setConstant(True)

          if 'sigmaCBCBG' in param.GetName():
            n = '_'.join(['sigmaCBCBG',year,lep,'cat'+cat,prod])
            hists.fill1DHist(param.getVal(), n,";sigmaCB", 50, 0.5,2.5, 1, '')
          elif 'fracCBG' in param.GetName():
            n = '_'.join(['fracCBS',year,lep,'cat'+cat,prod])
            hists.fill1DHist(param.getVal(), n,";fracCBG", 50, 0,0.5, 1, '')
          elif 'sigmaGCBG' in param.GetName():
            n = '_'.join(['sigmaGCBS',year,lep,'cat'+cat,prod])
            hists.fill1DHist(param.getVal(), n,";sigmaGCBG", 50, 0,5, 1, '')
          else:
            n = None


      normList.append(sig_ds_Low.sumEntries()*massDiff+sig_ds_Hi.sumEntries()*(1-massDiff))

      #acceff_low = (sig_ds_Low.sumEntries()/(xsBr(massLow)*lumi))
      #acceff_Hi  = (sig_ds_Hi.sumEntries()/(xsBr(massHi)*lumi))
      #yie_calc   = xsBr(mass)*lumi*(massDiff*acceff_low + (1-massDiff)*acceff_Hi)

      #print lep,year,cat,mass, cardDict[lep][year][cat][str(mass)]
      #print "\n\n Signal Yield?? = ", normList[-1], yie_calc
      #print "Low and Hi, sumentries:",  sig_ds_Low.sumEntries(), sig_ds_Hi.sumEntries()
      #print 'Low and Hi acc*eff:', acceff_low, acceff_Hi
      #if massLow!=massHi:
      #  raw_input('Yields for signals ')

      yieldName = '_'.join(['sig',prod,'yield',lep,year,'cat'+cat])
      #yieldVar  = RooRealVar(yieldName,yieldName, yie_calc)
      yieldVar  = RooRealVar(yieldName,yieldName, normList[-1])

      # sigNameInterp = '_'.join(['ds','sig',prod,lep,year,'cat'+cat,'M'+str(mass)])

      fitList.append(SigFit_Interp)

      getattr(cardDict[lep][year][cat][str(mass)],'import')(SigFit_Interp)
      getattr(cardDict[lep][year][cat][str(mass)],'import')(yieldVar)
      cardDict[lep][year][cat][str(mass)].commitTransaction()

    gStyle.SetOptStat(1)
    for n in paraNames[prod]:
      if n==None: continue
      s = hfile.Get(n)
      #s.Print()
      #raw_input('s print')
      s.Draw('hist')
      c.SaveAs(plotBase+'params_'+n+'.png')
    gStyle.SetOptStat(1)

    testFrame = mzg.frame(float(masses1[0])-10,float(masses1[-1])+10)

    for signal in dsList:
      signal.plotOn(testFrame, RooFit.MarkerStyle(6))
      # print signal.mean(mzg)
      # raw_input("Mean of the fit ")

    for i,fit in enumerate(fitList):
      regionName = fit.GetName().split('_')[-1]
      # fit.plotOn(testFrame)
      # fit.plotOn(testFrame, RooFit.NormRange('fitRegion_'+regionName))
      fit.plotOn(testFrame, RooFit.NormRange('fitRegion_'+regionName), RooFit.Normalization(normList[i],RooAbsReal.NumEvent),
                 RooFit.LineColor(TColor.GetColorPalette(i*10)))
      # fit.paramOn(testFrame)
      # testFrame.getAttText().SetTextSize(0.02)
      # testFrame.getAttText().SetTextColor(kRed)
    for i,signal in enumerate(dsList):
      signal.plotOn(testFrame, RooFit.MarkerStyle(20+i), RooFit.MarkerSize(1), RooFit.Binning(150))
      #signal.statOn(testFrame,RooFit.What("MR"),RooFit.Layout(0.60,0.83,0.87),RooFit.CutRange('statRange125'))
      print signal.mean(mzg), signal.sigma(mzg)
      # raw_input("Mean and sigma from dataset above")

    m = testFrame.GetMaximum()

    testFrame.SetMaximum(1.1*m)
    '''
    if prod == 'gg':
      testFrame.SetMaximum(1.1*m)
    elif prod =='vbf':
      testFrame.SetMaximum(0.05)
    elif prod =='v':
      testFrame.SetMaximum(0.03)
    '''
    testFrame.SetMinimum(1e-4)

    testFrame.Draw()
    if lep=="el":
      testFrame.SetTitle(";m_{ee#gamma} (GeV);Signal shape")
    else:
      testFrame.SetTitle(";m_{#mu#mu#gamma} (GeV);Signal shape")

    lat = TLatex()
    lat.SetNDC()
    lat.SetTextSize(0.045)
    lat.DrawLatex(0.66,0.85, 'Mean = %.1f'%dsList[0].mean(mzg))
    lat.DrawLatex(0.66,0.80, 'RMS  = %.2f'%dsList[0].sigma(mzg))

    CMS_lumi(c, 2, 11, "Simulation")
    c.SaveAs(plotBase+'_'.join(['sig','fit',prod,lep,year,'cat'+cat])+'.png')

  for prod in sigNameList:
    if lep=='el' and prod=='v': continue
    if cat in ['m1','m2','m3','m4','m5','m6','m7']:
      if prod!='gg': continue

    for mass in masses1:
      SignalNameParamFixer(year,lep,cat,prod,mass,cardDict[lep][year][cat][mass])

      #print cardDict[lep][year][cat][mass]
      #raw_input("RAW INPUT")

    #file = TFile("parampampam.root",'recreate')
    #file.cd()

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


  for mass in masses1:
    fileName = subdir+'/'+'_'.join(['SignalOutput',lep,year,'cat'+cat,mass])
    cardDict[lep][year][cat][mass].writeToFile(fileName+'.root')
    #cardDict[lep][year][cat][mass].writeToFile('testCards/'+fileName+'.root')


if __name__=="__main__":
  print len(sys.argv), sys.argv

  yearList   = [a.strip() for a in (cf.get("fits","yearList")).split(',')]
  leptonList = [a.strip() for a in (cf.get("fits","leptonList")).split(',')]
  catList    = [a.strip() for a in (cf.get("fits","catList")).split(',')]

  ver = cf.get("path","ver")

  for y in yearList:
    for c in catList:
      for l in leptonList:
        if l=='el' and c!='EB': continue
        SignalFitMaker(l,y,c, ver)

  print "\n \t\t All finished! \n"
