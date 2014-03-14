#!/usr/bin/env python
import sys,os
from ROOT import *
from collections import defaultdict
gROOT.SetBatch()
import ConfigParser as cp
cf = cp.ConfigParser()
cf.read('config.cfg')
s = cf.get("path","ver")

gROOT.ProcessLine(".L ~/tdrstyle.C")
setTDRStyle()
#plotBase = '/uscms_data/d2/andreypz/html/zgamma/dalitz/fits-'+s+'/toyFits/'
plotBase = '/tthome/andrey/html/zgamma/bias-mar11/'
toysDir = '../biasToys-mar11'

def loopThruPulls():
  #massList  = ['120','125','130']
  massList     = [a.strip()[0:3] for a in (cf.get("fits","massList")).split(',')]
  sigNameList  = [a.strip() for a in (cf.get("fits","sigNameList")).split(',')]
  yearList     = [a.strip() for a in (cf.get("fits","yearList")).split(',')]
  leptonList   = [a.strip() for a in (cf.get("fits","leptonList")).split(',')] 
  catList = ['0']

  genFuncList = ['Exp','Pow','Bern3']
  #genFuncList = ['GaussPow','GaussExp','SechPow','SechExp']
  for year in yearList:
    for lepton in leptonList:
      for genFunc in genFuncList:
        for cat in catList:
          for mass in massList:
            makePullPlots(year,lepton,genFunc,cat,mass)

def makePullPlots(year='2012', lepton='mu', genFunc='Exp', cat='0', mass='125'):

  #get the toy file and tree
  toyFileName = toysDir+'/toys_'+genFunc+'_m'+mass+'.root'
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
  for f,col in cf.items("colors"): colors[f] = int(col)

  #make a TCanvas and TLegend for each type of distribution

  #distList = ['sigPull', 'bgPull']
  distList = ['sigPull','bgPull','bgPull-2','bgPull-3','nSig','nBG','sigErr','bgErr','typeA']
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
          tmpHist.GetXaxis().SetTitle('nSig/#sigma(nSig)')

        elif dist == 'bgPull':
          toyTree.Draw('(('+fitFunc+'.yieldBkg - truth.yieldBkg)/'+fitFunc+'.yieldBkgErr)>>'+dist+'_'+fitFunc ,cutStr,'goff')
          tmpHist.GetXaxis().SetTitle('(nBG-nTrue)/#sigma(nBG)')
        elif dist == 'bgPull-2':
          toyTree.Draw('(('+fitFunc+'.yieldBkg - truth.yieldBkg)/sqrt(truth.yieldBkg))>>'+dist+'_'+fitFunc ,cutStr,'goff')
          tmpHist.GetXaxis().SetTitle('(nBG-nTrue)/#sqrt{nTrue}')
        elif dist == 'bgPull-3':
          toyTree.Draw('(('+fitFunc+'.yieldBkg - truth.yieldBkg)/sqrt('+fitFunc+'.yieldBkg))>>'+dist+'_'+fitFunc ,cutStr,'goff')
          tmpHist.GetXaxis().SetTitle('(nBG-nTrue)/#sqrt{nBG}')

        elif dist == 'nSig':
          toyTree.Draw('('+fitFunc+'.yieldSig)>>'+dist+'_'+fitFunc ,cutStr,'goff')
          tmpHist.GetXaxis().SetTitle('nSig')

        elif dist == 'nBG':
          toyTree.Draw('('+fitFunc+'.yieldBkg)>>'+dist+'_'+fitFunc ,cutStr,'goff')
          tmpHist.GetXaxis().SetTitle('nBG')

        elif dist == 'sigErr':
          toyTree.Draw('('+fitFunc+'.yieldSigErr)>>'+dist+'_'+fitFunc ,cutStr,'goff')
          tmpHist.GetXaxis().SetTitle('#sigma(nSig)')

        elif dist == 'bgErr':
          toyTree.Draw('('+fitFunc+'.yieldBkgErr)>>'+dist+'_'+fitFunc ,cutStr,'goff')
          tmpHist.GetXaxis().SetTitle('#sigma(nBG)')

        elif dist == 'typeA':
          toyTree.Draw('('+fitFunc+'.yieldSig/'+fitFunc+'.yieldBkgErr)>>'+dist+'_'+fitFunc ,cutStr,'goff')
          tmpHist.GetXaxis().SetTitle('nSig/#sigma(nBG)')

        tmpHist.GetYaxis().SetTitle('A.U.')
        tmpHist.GetYaxis().CenterTitle()

        print fitFunc,dist, tmpHist.GetEntries(), tmpHist.GetMean()
        if(tmpHist.Integral()>0): tmpHist.Scale(1./tmpHist.Integral())

        histListDict[dist].append(tmpHist)
        legList[i].AddEntry(tmpHist,fitFunc+': #mu={0:.2f}, #sigma={1:.2f}'.format(tmpHist.GetMean(), tmpHist.GetRMS()),'l')



  #make the plots

  newPath = plotBase+'/biasPulls/'+lepton+'_'+year+'/'+genFunc
  if not os.path.isdir(newPath):
    os.makedirs(newPath)

  #get a txt file ready for the latex
  f = open(plotBase+'/biasPulls/'+lepton+'_'+year+'/'+lepton+'_'+year+'_'+genFunc+'_cat'+cat+'_mH'+mass+'.txt', 'w')
  for i,dist in enumerate(distList):
    canList[i].cd()

    ymax = max(map(lambda x:x.GetMaximum(),histListDict[dist]))*1.1 #this is awesome if it works, python rocks

    histListDict[dist][0].Draw()
    histListDict[dist][0].SetMaximum(ymax)
    #print  histListDict[dist][0]
    #raw_input()
    for j in range(1,len(histListDict[dist])):
      #if dist in ['sigPull','typeA'] and j==len(histListDict[dist])-1: continue
      histListDict[dist][j].Draw('same')
    legList[i].Draw('same')
    canList[i].SaveAs(plotBase+'/biasPulls/'+lepton+'_'+year+'/'+genFunc+'/'+dist+'_'+lepton+'_'+year+'_'+genFunc+'_cat'+cat+'_mH'+mass+'.png')

    if dist == 'typeA':
      f.write('typeA:\n')
      for hist in histListDict[dist]:
        f.write(hist.GetName().strip('typeA_')+', ')
      f.write('\n')
      for hist in histListDict[dist]:
        f.write('{0:.2f}, '.format(hist.GetMean()))
      f.write('\n')
    if dist == 'bgPull':
      f.write('bgPull:\n')
      for hist in histListDict[dist]:
        f.write(hist.GetName().strip('bgPull_')+', ')
      f.write('\n')
      for hist in histListDict[dist]:
        f.write('{0:.2f}, '.format(hist.GetMean()))
      f.write('\n')
  f.close()


if __name__=="__main__":
  #mass = options.mass
  #makePullPlots(mass=mass)
  loopThruPulls()
  









