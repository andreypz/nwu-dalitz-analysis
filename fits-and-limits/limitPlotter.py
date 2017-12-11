#!/usr/bin/env python
import os,sys

import numpy as np
from ROOT import *
gROOT.SetBatch()
sys.path.append("../zgamma")
sys.path.append("../scripts")
import makeHTML as ht
import utils as u
from optparse import OptionParser
parser = OptionParser(usage="usage: %prog  [options]")
parser.add_option("--br", dest="br",  action="store_true", default=False, help="Do the limit on BR*cs instead of the mu")
parser.add_option("--mll",dest="mll", action="store_true", default=False, help="Do the limit on BR*cs in bins of mll")
parser.add_option("--divide",dest="divide", action="store_true", default=False, help="Divide by the binn-size (for --mll option)")
parser.add_option("--cat",dest="cat", default='HR9', help="Category: EB, EE, mll50, comb (for combination)")
parser.add_option("--lep",dest="lep", default='mu', help="Channel: mu or el")
parser.add_option("--path",dest="path", default=None, help="path where the limit results are stored")
(opt, args) = parser.parse_args()

import ConfigParser as cp
cf = cp.ConfigParser()
cf.read('config.cfg')

gStyle.SetOptTitle(0)

s = cf.get("path","ver")
if opt.path!=None:
  s = opt.path
lep = opt.lep
#s = "comboCardsMing/"
#s = "MingDataCards/Dalitz_electron/"

method = 'Asymptotic'
#method = 'ProfileLikelihood'
out = TFile(s+"/limit-data-cat"+opt.cat+"-"+lep+".root","recreate")

plotBase = cf.get("path","base")+'/fits-'+s

massList  = [float(a) for a in (cf.get("fits","massList-more")).split(',')]
#massList   = [float(a) for a in u.drange(120,150,1.0)]
doBlind   = int(cf.get("fits","blind"))
doObs = not doBlind
mllBins = u.mllBins()

if opt.cat!='comb': suff=opt.cat+'_'+lep
else:  suff='Combo'
# print massList
if opt.br:  combineOutDir = '/combineOut_xsbr/'
else:       combineOutDir = '/combineOut/'


c = TCanvas("c","c",0,0,500,400)
c.cd()

if __name__ == "__main__":

  xAxis = []
  xErr  = []
  obs = []
  exp = []
  exp1SigHi = []
  exp1SigLow = []
  exp2SigHi = []
  exp2SigLow = []
  SM = []
  fileListTmp = os.listdir(s+combineOutDir)
  #print fileListTmp
  fileList = filter(lambda fileName: 'higgsCombineTest.'+method in fileName,fileListTmp)
  #print fileList
  SMtot = 0.732

  if opt.mll:
    for mbin in ['m1','m2','m3','m4','m5','m6','m7']:
      #center = mllBins[int(mbin[1])-1][0]+ 0.5*(mllBins[int(mbin[1])][0]-mllBins[int(mbin[1])-1][0])
      r1 = mllBins[int(mbin[1])-1][0]
      r2 = mllBins[int(mbin[1])][0]
      center = 0.5*(r2+r1)
      print 'mbin and center',mbin, mllBins[int(mbin[1])], center

      thisFile = "higgsCombineTest.Asymptotic.mH125.0_cat_"+mbin+"_"+lep+".root"
      # print thisFile
      xAxis.append(float(center))
      xErr.append(float(0.5*(r2-r1)))
      dmScale = 1
      if opt.divide:
        dmScale = r2-r1

      SM.append(10*SMtot*mllBins[int(mbin[1])][1]/dmScale)
      print '10x SM prediction:', SM[-1]

      f = TFile(s+combineOutDir+thisFile,"open")
      # f.Print()
      tree = f.Get("limit")

      for i,l in enumerate(tree):
        # print i, l, l.limit
        if i==0: exp2SigLow.append(float(l.limit)/dmScale)
        if i==1: exp1SigLow.append(float(l.limit)/dmScale)
        if i==2: exp.append(float(l.limit)/dmScale)
        if i==3: exp1SigHi.append(float(l.limit)/dmScale)
        if i==4: exp2SigHi.append(float(l.limit)/dmScale)
        if i==5: obs.append(float(l.limit)/dmScale)

  else:
    for mass in massList:
      #print str(mass)
      #if str(mass)[4] == "0": # nasty trick
      #  mass = int(mass)
        #print 'we are her'
      filtName = str(mass)+'_cat_'+suff+'.root'
      print mass, filtName
      thisFile = filter(lambda fileName: filtName in fileName,fileList)[0]
      print thisFile
      xAxis.append(float(mass))
      f = TFile(s+combineOutDir+thisFile,"open")
      #f.Print()
      tree = f.Get("limit")
      #tree.Print('all')
      for i,l in enumerate(tree):
        #print i, l, l.limit
        if method =='ProfileLikelihood' and i==0:
          obs.append(float(l.limit))
          exp2SigLow.append(float(l.limit))
          exp1SigLow.append(float(l.limit))
          exp.append(float(l.limit))
          exp1SigHi.append(float(l.limit))
          exp2SigHi.append(float(l.limit))
          continue

        if i==0: exp2SigLow.append(float(l.limit))
        if i==1: exp1SigLow.append(float(l.limit))
        if i==2: exp.append(float(l.limit))
        if i==3: exp1SigHi.append(float(l.limit))
        if i==4: exp2SigHi.append(float(l.limit))
        if i==5: obs.append(float(l.limit))


  print 'masses:', xAxis
  print 'exp:\n', ["{0:0.2f}".format(i) for i in exp]
  if doObs:
    print 'obs:\n',["{0:0.2f}".format(i) for i in obs]
  #print exp2SigLow
  #print exp1SigLow
  #print exp1SigHi
  #print exp2SigHi

  exp2SigLowErr = [a-b for a,b in zip(exp,exp2SigLow)]
  exp1SigLowErr = [a-b for a,b in zip(exp,exp1SigLow)]
  exp2SigHiErr  = [b-a for a,b in zip(exp,exp2SigHi)]
  exp1SigHiErr  = [b-a for a,b in zip(exp,exp1SigHi)]

  print '-1 sigma:\n', ["{0:0.2f}".format(i) for i in exp1SigLowErr]
  print '+1 sigma:\n', ["{0:0.2f}".format(i) for i in exp1SigHiErr]

  print 'exp + 1 sigma up:\n',   ["{0:0.2f}".format(i) for i in exp1SigHi]
  print 'exp - 1 sigma down:\n', ["{0:0.2f}".format(i) for i in exp1SigLow]


  if opt.mll: print '\n 10x SM prediction \n',SM

  zeros_Array = np.zeros(len(xAxis),dtype = float)
  xAxis_Array = np.array(xAxis)
  if opt.mll:
    xErr_Array = np.array(xErr)
  else:
    xErr_Array = zeros_Array
  obs_Array = np.array(obs)
  exp_Array = np.array(exp)
  exp2SigLowErr_Array = np.array(exp2SigLowErr)
  exp1SigLowErr_Array = np.array(exp1SigLowErr)
  exp1SigHiErr_Array  = np.array(exp1SigHiErr)
  exp2SigHiErr_Array  = np.array(exp2SigHiErr)

  if opt.mll: SM_Array = np.array(SM)

  mg = TMultiGraph()
  mg.SetTitle('')

  nPoints  = len(xAxis)
  expected = TGraphAsymmErrors(nPoints,xAxis_Array,exp_Array,xErr_Array,xErr_Array,zeros_Array,zeros_Array)
  oneSigma = TGraphAsymmErrors(nPoints,xAxis_Array,exp_Array,xErr_Array,xErr_Array,exp1SigLowErr_Array,exp1SigHiErr_Array)
  twoSigma = TGraphAsymmErrors(nPoints,xAxis_Array,exp_Array,xErr_Array,xErr_Array,exp2SigLowErr_Array,exp2SigHiErr_Array)
  observed = TGraphAsymmErrors(nPoints,xAxis_Array,obs_Array,zeros_Array,zeros_Array,zeros_Array,zeros_Array)

  if opt.mll:
    standardM = TGraphAsymmErrors(nPoints,xAxis_Array,SM_Array,xErr_Array,xErr_Array,zeros_Array,zeros_Array)

  #expected.Print("all")
  if opt.mll:
    oneSigma.SetFillColor(kCyan-10)
    twoSigma.SetFillColor(kOrange-4)
  else:
    oneSigma.SetFillColor(kGreen)
    twoSigma.SetFillColor(kYellow)

  expected.SetMarkerColor(kBlack)
  #expected.SetMarkerStyle(kFullCircle)
  #expected.SetMarkerSize(1.5)
  expected.SetLineColor(kBlack)
  expected.SetLineWidth(2)
  expected.SetLineStyle(2)


  if opt.mll:
    standardM.SetFillColor(kRed+1)
    standardM.SetLineColor(kRed+1)
    standardM.SetLineWidth(2)
    standardM.SetMarkerColor(kRed+1)

    observed.SetMarkerColor(kBlue+2)
    observed.SetMarkerSize(1.2)
    observed.SetMarkerStyle(kFullCircle)

    mg.Add(twoSigma,'P2')
    mg.Add(oneSigma,'P2')
    mg.Add(expected,'P')
  else:
    observed.SetLineWidth(2)

    mg.Add(twoSigma)
    mg.Add(oneSigma)
    mg.Add(expected)

  if doObs:
    if opt.mll:
      mg.Add(observed,'P')
    else:
      mg.Add(observed)

  if opt.mll:
    mg.Add(standardM,'P')

  mg.SetMinimum(0)

  if opt.mll:
    mg.Draw('A')
    mg.GetXaxis().SetTitle('m_{#mu#mu} (GeV)')
    mg.GetXaxis().SetLimits(0, 80)
    c.SetLogx()
    #c.SetLogy()
  else:
    mg.Draw('AL3')
    mg.GetXaxis().SetTitle('m_{H} (GeV)')
    mg.GetXaxis().SetLimits(massList[0],massList[-1])
  if opt.br:
    #mg.GetXaxis().SetTitle('m_{H} (GeV)')
    mg.GetYaxis().SetTitle('#sigma(pp #rightarrow H) #times B(H #rightarrow #mu#mu#gamma)_{95% CL} (fb)')
    mg.SetMaximum(21)
    if opt.lep=='el':
      mg.GetYaxis().SetTitle('#sigma(pp #rightarrow H) #times B(H #rightarrow ee#gamma)_{95% CL} (fb)')
      mg.SetMaximum(52)
    if opt.mll and opt.lep=='mu':
      mg.SetMaximum(12)
      if opt.divide:
        mg.GetYaxis().SetTitle('#frac{1}{#Deltam} #sigma(pp #rightarrow H)#times B(H#rightarrow#mu#mu#gamma) (fb/GeV)')
        mg.SetMaximum(200)
        mg.SetMinimum(0.03)
        c.SetLogy()
  else:
    mg.GetYaxis().SetTitle('95% CL limit on #sigma/#sigma_{SM}')
    if opt.cat in ['EB']:
      mg.SetMaximum(42)
    elif opt.cat in ['comb']:
      mg.SetMaximum(9)
    elif opt.cat in ['mll50']:
      mg.SetMaximum(150)
    elif opt.cat in ['HR9','LR9','VBF']:
      mg.SetMaximum(30)
    elif opt.cat in ['EE']:
      mg.SetMaximum(50)
  gPad.RedrawAxis()

  lat = TLatex()
  lat.SetNDC()
  lat.SetTextSize(0.04)
  lat.SetTextFont(42)

  if opt.mll:
    lat.DrawLatex(0.18,0.95, 'H #rightarrow #gamma*#gamma #rightarrow #mu#mu#gamma')
    lat.SetTextSize(0.05)
    lat.DrawLatex(0.40,0.95,'m_{H} = 125 GeV')
    lat.SetTextSize(0.035)
  elif opt.br and opt.lep=='mu':
    lat.DrawLatex(0.20,0.75, 'H #rightarrow #gamma*#gamma #rightarrow #mu#mu#gamma')
  elif opt.br and opt.lep=='el':
    lat.DrawLatex(0.20,0.75, 'H #rightarrow #gamma*#gamma #rightarrow ee#gamma')
  elif opt.cat=='comb':
    lat.DrawLatex(0.20,0.72, '#splitline{H #rightarrow #gamma*#gamma #rightarrow ll#gamma}{(#mu#mu#gamma and ee#gamma)}')
  else:
    if lep=="mu":
      lat.DrawLatex(0.20,0.75, 'H #rightarrow #gamma*#gamma #rightarrow #mu#mu#gamma')
    elif lep=="el":
      lat.DrawLatex(0.20,0.75, 'H #rightarrow #gamma*#gamma #rightarrow ee#gamma')
  #CMS_lumi(c, 2, 11, "")

  leg = TLegend(0.50,0.68,0.75,0.88)
  if doObs:
    leg = TLegend(0.50,0.66,0.75,0.89)
    if opt.mll:
      leg.AddEntry(observed,"Observed", "p")
    else:
      leg.AddEntry(observed,"Observed", "l")
  leg.AddEntry(expected,"Expected", "l")
  leg.AddEntry(oneSigma,"Expected #pm 1#sigma", "f")
  leg.AddEntry(twoSigma,"Expected #pm 2#sigma", "f")

  #selection = TLatex()
  #if "mu" in s:
  #  selection = TLatex(0.70,0.95, "#mu#mu selection");
  #  selection.SetTextColor(kBlue-3);
  #elif "el" in s:
  #  selection = TLatex(0.70,0.95, "#font[32]{ee} selection");
  #  selection.SetTextColor(kGreen-3);
  #selection.SetNDC();
  #selection.Draw();

  leg.SetTextFont(42)
  leg.SetTextSize(0.04)
  leg.SetFillColor(kWhite)
  leg.Draw()

  if opt.br and not opt.mll:
    f = TFile('../data/Dalitz_BR20.root','READ')
    g = f.Get('csbr_'+opt.lep)
    for i in range(g.GetN()):
      g.GetY()[i] *= 10

    #fit = g.GetFunction('pol4')
    #fit.Print()
    g.SetLineColor(kRed+2)
    g.SetLineWidth(2)
    #g.SetMarkerStyle(0)

    g.Draw('same L')
    #print fit(125)
    #g.Print('all')
    gStyle.SetOptFit(0)
    leg.SetY1NDC(0.64)
    leg.AddEntry(g,"10 #times SM", "l")
  elif not opt.br:
    l = TLine(120, 1, 130, 1)
    l.SetLineColor(kRed+2)
    l.SetLineWidth(2)
    #l.SetLineStyle(21)
    l.Draw()
  elif opt.mll:
    leg.SetY1NDC(0.64)
    leg.AddEntry(standardM,"10 #times SM", "lp")


  for e in ['.png', '.pdf']:
    CMS_lumi(c, 4, 11, "")
    limitDir = plotBase+'/Limits'
    if opt.br:
      limitDir = plotBase+'/Limits-xsBr'

    u.createDir(limitDir)

    if opt.mll and opt.divide:
      c.SaveAs(limitDir+'/limits_Mll_dM'+e)
    elif opt.mll:
      c.SaveAs(limitDir+'/limits_Mll'+e)
    else:
      if opt.br:
        c.SaveAs(limitDir+'/limits_xsBR_cat_'+suff+e)
      else:
        c.SaveAs(limitDir+'/limits_cat_'+suff+e)

  out.cd()
  observed.Write("observed")
  expected.Write("expected")
  oneSigma.Write("oneSigma")



  comments = ['Limits and fits']
  defaultPage = 'BestFits'

  plot_types =[]
  dirlist = os.listdir(plotBase)
  for d in dirlist:
    if os.path.isdir(plotBase+"/"+d):
      plot_types.append(d)
  ht.makeHTML("h &rarr; dalitz decay plots", plotBase, plot_types, comments, defaultPage, doYields=False)

