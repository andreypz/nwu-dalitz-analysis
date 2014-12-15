#!/usr/bin/env python
import sys,os
from ROOT import *
from collections import defaultdict
gROOT.SetBatch()
sys.path.append("../zgamma")
#gROOT.ProcessLine(".L ~/tdrstyle.C")
#setTDRStyle()
import utils as u
from optparse import OptionParser
parser2 = OptionParser(usage="usage: %prog [--opt]")
parser2.add_option("--copy", dest="copy",  default=None, help="Directory where the results from Batch are stored. Will copy them in a #better place.")
(opt, args) = parser2.parse_args()

#from apzFitProducer import AutoVivification

import ConfigParser as cp
cf = cp.ConfigParser()
cf.read('config.cfg')
s = cf.get("path","ver")
plotBase = cf.get("path","htmlbase")+'/html/zgamma/dalitz/fits-'+s+'/Pulls'
#plotBase = '/tthome/andrey/html/zgamma/bias-mar11/'
#toysDir = '../biasToys-2014-Sep23'
toysDir = '../biasToys-2014-Dec08'

sigMuDict = u.AutoVivification()
taMuDict  = u.AutoVivification()

#massList  = ['120','125','130']
massList     = [a.strip()[0:3] for a in (cf.get("fits","massList")).split(',')]
sigNameList  = [a.strip() for a in (cf.get("fits","sigNameList")).split(',')]
yearList     = [a.strip() for a in (cf.get("fits","yearList")).split(',')]
leptonList = ['mu']
catList    = ['EB']

genFuncList = ['Exp','Pow','Laurent']

def loopThruPulls():
  for year in yearList:
    for lep in leptonList:
      for mass in massList:
        for cat in catList:
          for genFunc in genFuncList:
            if opt.copy!=None:
              bdir = opt.copy
              toyFileName  = toysDir+'_'.join(['/biasToys',year,lep,'cat'+cat,genFunc,mass])+'.root'
              toyFileNameJ = bdir   +'_'.join(['/biasToys',year,lep,'cat'+cat,genFunc,mass])+'_job*.root'
              os.system(' '.join(['../scripts/ahadd.py', toyFileName, toyFileNameJ]))
            makePullPlots(year,lep,genFunc,cat,mass)

def makePullPlots(year='2012', lepton='mu', genFunc='Exp', cat='EB', mass='125'):

  #get the toy file and tree
  toyFileName  = toysDir+'_'.join(['/biasToys',year,lepton,'cat'+cat,genFunc,mass])+'.root'
  print toyFileName
  toyFile = TFile(toyFileName)
  toyTree = toyFile.Get('toys')

  #toyTree.Print()

  gStyle.SetOptStat(0)
  gStyle.SetTitleFontSize(2.7)
  gStyle.SetTitleH(0.06) # Set the height of the title box
  gStyle.SetTitleW(0.2)    # Set the width of the title box
  gStyle.SetTitleX(0.4)    # Set the position of the title box
  gStyle.SetTitleY(0.99)    # Set the position of the title box
  gStyle.SetLineWidth(2)
  #define the fit functions we'll be using today (and their associated plot colors)

  #turnOnList = ['Sech','Gauss']
  turnOnList = ['']
  tailList = ['Bern2','Bern3','Bern4','Bern5']
  #tailList.append('gen')
  colors ={}
  for f,col in cf.items("colors"): colors[f] = eval(col)

  #make a TCanvas and TLegend for each type of distribution

  distList = ['sigPull', 'bgPull', 'typeA']
  #distList = ['sigPull','bgPull','bgPull-2','bgPull-3','nSig','nBG','sigErr','bgErr','typeA']
  canList = []
  legList = []
  for dist in distList:
    canvasTmp = TCanvas(dist+'Can',dist+'Can',700,600)
    canList.append(canvasTmp)

    legendTmp = TLegend(0.66,0.70,0.98,0.95)
    legendTmp.SetFillColor(0)
    #legendTmp.SetFillStyle(0)
    legendTmp.SetTextSize(0.03)
    legList.append(legendTmp)

  #this is a dictionary where each key is the distribution string, and each value is the list of histograms associated with that dist
  histListDict = defaultdict(list)

  #start the main loop and make all the histos

  # build the fit func name
  for turnOn in turnOnList:
    for tail in tailList:
      fitFunc = turnOn+tail
      cutStr = ''
      cutStr = fitFunc+'.statusAll==0&&'+fitFunc+'.covQual>=1&&'+fitFunc+'.statusMIGRAD==0'
      if fitFunc in ['Bern3','Bern4','Bern5']:
        cutStr = cutStr+'&&'+fitFunc+'.yieldBkgErr>1.5&&'+fitFunc+'.yieldSigErr<30&&fabs('+fitFunc+'.yieldSig)<28'
        cutStr = cutStr+'&&'+fitFunc+'.paramP1Err<15&&'+fitFunc+'.paramP2Err<15&&'+fitFunc+'.paramP3Err<15'
      #go through all the distributions you want
      for i,dist in enumerate(distList):

        if dist in ['sigPull','bgPull','bgPull-2','bgPull-3']:
          tmpHist = TH1F(dist+'_'+fitFunc, dist+'_'+fitFunc, 100, -6, 6)
        elif dist in ['sigErr','bgErr']:
          tmpHist = TH1F(dist+'_'+fitFunc, dist+'_'+fitFunc, 100, 0, 20)
        elif dist in ['nSig']:
          tmpHist = TH1F(dist+'_'+fitFunc, dist+'_'+fitFunc, 100, -60, 60)
        elif dist in ['nBG']:
          tmpHist = TH1F(dist+'_'+fitFunc, dist+'_'+fitFunc, 100, 0, 150)
        elif dist == 'typeA':
          tmpHist = TH1F(dist+'_'+fitFunc, dist+'_'+fitFunc, 100, -10, 10)

        tmpHist.SetLineWidth(2)
        tmpHist.SetLineColor(colors[fitFunc.lower()])
        tmpHist.SetTitle(dist)
        canList[i].cd()

        #build the histos

        if dist == 'sigPull':
          toyTree.Draw('('+fitFunc+'.yieldSig/'+fitFunc+'.yieldSigErr)>>'+dist+'_'+fitFunc ,cutStr,'goff')
          tmpHist.SetTitle(mass+';nSig/#sigma(nSig);A.U.')

        elif dist == 'bgPull':
          toyTree.Draw('(('+fitFunc+'.yieldBkg - truth.yieldBkg)/'+fitFunc+'.yieldBkgErr)>>'+dist+'_'+fitFunc ,cutStr,'goff')
          tmpHist.SetTitle(mass+';(nBG-nTrue)/#sigma(nBG);A.U.')
        elif dist == 'bgPull-2':
          toyTree.Draw('(('+fitFunc+'.yieldBkg - truth.yieldBkg)/sqrt(truth.yieldBkg))>>'+dist+'_'+fitFunc ,cutStr,'goff')
          tmpHist.SetTitle(mass+';(nBG-nTrue)/#sqrt{nTrue};A.U.')
        elif dist == 'bgPull-3':
          toyTree.Draw('(('+fitFunc+'.yieldBkg - truth.yieldBkg)/sqrt('+fitFunc+'.yieldBkg))>>'+dist+'_'+fitFunc ,cutStr,'goff')
          tmpHist.SetTitle(mass+';(nBG-nTrue)/#sqrt{nBG};A.U.')

        elif dist == 'nSig':
          toyTree.Draw('('+fitFunc+'.yieldSig)>>'+dist+'_'+fitFunc ,cutStr,'goff')
          tmpHist.SetTitle(mass+';nSig;A.U.')

        elif dist == 'nBG':
          toyTree.Draw('('+fitFunc+'.yieldBkg)>>'+dist+'_'+fitFunc ,cutStr,'goff')
          tmpHist.SetTitle(mass+';nBG;A.U.')

        elif dist == 'sigErr':
          toyTree.Draw('('+fitFunc+'.yieldSigErr)>>'+dist+'_'+fitFunc ,cutStr,'goff')
          tmpHist.SetTitle(mass+';#sigma(nSig);A.U.')

        elif dist == 'bgErr':
          toyTree.Draw('('+fitFunc+'.yieldBkgErr)>>'+dist+'_'+fitFunc ,cutStr,'goff')
          tmpHist.SetTitle(mass+';#sigma(nBG);A.U.')

        elif dist == 'typeA':
          toyTree.Draw('('+fitFunc+'.yieldSig/'+fitFunc+'.yieldBkgErr)>>'+dist+'_'+fitFunc ,cutStr,'goff')
          tmpHist.SetTitle(mass+';nSig/#sigma(nBG);A.U.')

        tmpHist.GetYaxis().CenterTitle()

        print fitFunc,dist, tmpHist.GetEntries(), tmpHist.GetMean()
        if(tmpHist.Integral()>0): tmpHist.Scale(1./tmpHist.Integral())

        histListDict[dist].append(tmpHist)
        legList[i].AddEntry(tmpHist,fitFunc+': #mu={0:.2f}, #sigma={1:.2f}'.format(tmpHist.GetMean(), tmpHist.GetRMS()),'l')



  #make the plots


  newPath = plotBase+'-'.join(['',year,lepton,'cat'+cat,genFunc])

  if not os.path.isdir(newPath):
    os.makedirs(newPath)

  for i,dist in enumerate(distList):
    canList[i].cd()

    ymax = max(map(lambda x:x.GetMaximum(),histListDict[dist]))*1.1

    histListDict[dist][0].Draw()
    histListDict[dist][0].SetMaximum(ymax)
    #print  histListDict[dist][0]
    #raw_input()
    for j in range(1,len(histListDict[dist])):
      #if dist in ['sigPull','typeA'] and j==len(histListDict[dist])-1: continue
      histListDict[dist][j].Draw('same')
    legList[i].Draw('same')
    canList[i].SaveAs(newPath+'/'+dist+'_'+lepton+'_'+year+'_'+genFunc+'_cat'+cat+'_mH'+mass+'.png')


    mus = []
    if dist == 'typeA':
      for hist in histListDict[dist]:
        mus.append(hist.GetMean())
      taMuDict[year][lepton][cat][mass][genFunc]=mus
    if dist == 'sigPull':
      for hist in histListDict[dist]:
        mus.append(hist.GetMean())
      sigMuDict[year][lepton][cat][mass][genFunc]=mus

  print sigMuDict
  print taMuDict

if __name__=="__main__":
  #mass = options.mass
  #makePullPlots(mass=mass)
  loopThruPulls()

  table1 = []
  table2 = []
  for year in yearList:
    for lep in leptonList:
      for cat in catList:
        for mass in massList:
          line1 = [mass]
          line2 = [mass]
          for genFunc in genFuncList:
            line1.extend(sigMuDict[year][lep][cat][mass][genFunc][1:])
            line2.extend( taMuDict[year][lep][cat][mass][genFunc][1:])
            print genFunc, line1
            print genFunc, line2
          table1.append(line1)
          table2.append(line2)
        #print table1
        #print table2


  u.makeTable(table1,"pulls", "twiki", precision='%.2f')
  u.makeTable(table1,"pulls1", "tex",  precision='%.2f')
  u.makeTable(table2,"pulls2", "tex",  precision='%.2f')






