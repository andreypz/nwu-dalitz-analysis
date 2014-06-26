#!/usr/bin/env python
import os,sys

import numpy as np
from ROOT import *
gROOT.SetBatch()
sys.path.append("../zgamma")
import utils as u
from optparse import OptionParser
parser = OptionParser(usage="usage: %prog  [options]")
parser.add_option("--br",dest="br", action="store_true", default=False, help="Do the limit on BR*cs instead of the mu")
(options, args) = parser.parse_args()

import ConfigParser as cp
cf = cp.ConfigParser()
cf.read('config.cfg')

gROOT.ProcessLine(".L ../tdrstyle.C")
setTDRStyle()
gStyle.SetOptTitle(0)

doObs = 1

s = cf.get("path","ver")
method = 'Asymptotic'
#method = 'ProfileLikelihood'
out = TFile(s+"/limit-data.root","recreate")

plotBase = cf.get("path","htmlbase")+'/html/zgamma/dalitz/fits-'+s
u.createDir(plotBase)

fullCombo = True
massList   = [float(a) for a in (cf.get("fits","massList-more")).split(',')]
#massList        = [float(a) for a in u.drange(120,150,5)]

print massList

c = TCanvas("c","c",0,0,500,400)
c.cd()

if fullCombo:
  xAxis = []
  obs = []
  exp = []
  exp1SigHi = []
  exp1SigLow = []
  exp2SigHi = []
  exp2SigLow = []
  fileListTmp = os.listdir(s)
  fileList = filter(lambda fileName: 'higgsCombineTest.'+method in fileName,fileListTmp)
  print fileList
  for mass in massList:
    #print str(mass)
    if str(mass)[4] == "0": # nasty trick
      mass = int(mass)
      #print 'we are her'
    #print mass
    thisFile = filter(lambda fileName: str(mass)+'.root' in fileName,fileList)[0]
    print thisFile
    xAxis.append(float(mass))
    f = TFile(s+"/"+thisFile,"open")
    #f.Print()
    tree = f.Get("limit")
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
  print 'exp:',exp
  if doObs:
    print 'obs:',obs
  #print exp2SigLow
  #print exp1SigLow
  #print exp1SigHi
  #print exp2SigHi

  exp2SigLowErr = [a-b for a,b in zip(exp,exp2SigLow)]
  exp1SigLowErr = [a-b for a,b in zip(exp,exp1SigLow)]
  exp2SigHiErr = [fabs(a-b) for a,b in zip(exp,exp2SigHi)]
  exp1SigHiErr = [fabs(a-b) for a,b in zip(exp,exp1SigHi)]

  print '2 sig low:',exp2SigLowErr
  print '1 sig low:',exp1SigLowErr
  print '1 sig hi:',exp1SigHiErr
  print '2 sig hi:',exp2SigHiErr

  xAxis_Array = np.array(xAxis)
  obs_Array = np.array(obs)
  exp_Array = np.array(exp)
  exp2SigLowErr_Array = np.array(exp2SigLowErr)
  exp1SigLowErr_Array = np.array(exp1SigLowErr)
  exp1SigHiErr_Array  = np.array(exp1SigHiErr)
  exp2SigHiErr_Array  = np.array(exp2SigHiErr)
  zeros_Array = np.zeros(len(xAxis),dtype = float)

  mg = TMultiGraph()
  mg.SetTitle('')

  nPoints = len(xAxis)
  expected = TGraphAsymmErrors(nPoints,xAxis_Array,exp_Array,zeros_Array,zeros_Array,zeros_Array,zeros_Array)
  oneSigma = TGraphAsymmErrors(nPoints,xAxis_Array,exp_Array,zeros_Array,zeros_Array,exp1SigLowErr_Array,exp1SigHiErr_Array)
  twoSigma = TGraphAsymmErrors(nPoints,xAxis_Array,exp_Array,zeros_Array,zeros_Array,exp2SigLowErr_Array,exp2SigHiErr_Array)
  observed = TGraphAsymmErrors(nPoints,xAxis_Array,obs_Array,zeros_Array,zeros_Array,zeros_Array,zeros_Array)

  #expected.Print("all")

  oneSigma.SetFillColor(kGreen)
  twoSigma.SetFillColor(kYellow)

  #expected.SetMarkerColor(kBlack)
  #expected.SetMarkerStyle(kFullCircle)
  #expected.SetMarkerSize(1.5)
  expected.SetLineColor(kBlack)
  expected.SetLineWidth(2)
  expected.SetLineStyle(2)

  observed.SetLineWidth(2)
  observed.SetMarkerStyle(kFullCircle)

  mg.Add(twoSigma)
  mg.Add(oneSigma)
  mg.Add(expected)
  if doObs:
    mg.Add(observed)


  mg.SetMinimum(0)

  mg.Draw('AL3')
  mg.GetXaxis().SetTitle('m_{H} (GeV)')
  mg.GetXaxis().SetLimits(massList[0],massList[-1]);
  if options.br:
    mg.GetYaxis().SetTitle('#sigma(gg #rightarrow a)#timesBR(a#rightarrow#mu#mu#gamma) (fb)')
    mg.SetMaximum(21)
  else:
    mg.GetYaxis().SetTitle('95% CL limit on #sigma#timesBR/#sigma_{SM}#times BR_{SM}')
    mg.SetMaximum(31)
  gPad.RedrawAxis()

  prelim = TLatex()
  prelim.SetNDC();
  prelim.SetTextSize(0.04)
  prelim.DrawLatex(0.20,0.95, "CMS Preliminary")
  prelim.DrawLatex(0.60,0.8, "H#rightarrow#gamma*#gamma#rightarrow#mu#mu#gamma")
  prelim.DrawLatex(0.60,0.74, "m_{#mu#mu} < 20 GeV")
  #  prelim.Draw()

  sqrt = TLatex(0.50,0.95, "#sqrt{s} = 8 Tev  #it{L_{int}}  =  19.7  fb^{-1}")
  sqrt.SetNDC();
  sqrt.SetTextSize(0.04);
  sqrt.Draw();

  leg = TLegend(0.25,0.7,0.50,0.88)
  if doObs:
    leg.AddEntry(observed,"Observed", "l")
  leg.AddEntry(expected,"Expected", "l")
  leg.AddEntry(oneSigma,"Expected #pm 1#sigma", "f")
  leg.AddEntry(twoSigma,"Expected #pm 2#sigma", "f")

  selection = TLatex()
  if "mu" in s:
    selection = TLatex(0.70,0.95, "#mu#mu selection");
    selection.SetTextColor(kBlue-3);
  elif "el" in s:
    selection = TLatex(0.70,0.95, "#font[32]{ee} selection");
    selection.SetTextColor(kGreen-3);

  #selection.SetNDC();
  #selection.Draw();

  leg.SetTextSize(0.04)
  leg.SetFillColor(kWhite)
  leg.Draw()
  if options.br:
    f = TFile('../data/Dalitz_BR20.root','READ')
    g = f.Get('csbr_mu')
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
    leg.AddEntry(g,"10#timesSM", "l")

  for e in ['.png', '.pdf']:
    c.SaveAs(plotBase+'/Limits'+e)
    
  out.cd()
  observed.Write("observed")
  expected.Write("expected")
  oneSigma.Write("oneSigma")
