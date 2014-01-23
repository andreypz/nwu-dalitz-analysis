#!/usr/bin/env python
import os,sys

import numpy as np
from ROOT import *
gROOT.SetBatch()
sys.path.append("../zgamma")
import utils as u

gROOT.ProcessLine(".L ~/tdrstyle.C")
setTDRStyle()
gStyle.SetOptTitle(0)

print len(sys.argv), sys.argv
if len(sys.argv) != 2:
  sys.exit()
s = sys.argv[1]
  
out = TFile(s+"/limit-data.root","recreate")

plotBase = '/uscms_data/d2/andreypz/html/zgamma/dalitz/fits-'+s
u.createDir(plotBase)

fullCombo = True
byParts = False

#massList = [125.0]
massList = [120.0,125.0,130.0,135.0,140.0,145.0,150.0]

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
  fileList = filter(lambda fileName: 'higgsCombineTest.Asymptotic' in fileName,fileListTmp)
  print fileList
  for mass in massList:
    thisFile = filter(lambda fileName: str(int(mass)) in fileName,fileList)[0]
    print thisFile
    xAxis.append(mass)
    f = TFile(s+"/"+thisFile,"open")
    #f.Print()
    tree = f.Get("limit")
    for i,l in enumerate(tree):
      #print i, l, l.limit
      if i==0: exp2SigLow.append(float(l.limit))
      if i==1: exp1SigLow.append(float(l.limit))
      if i==2: exp.append(float(l.limit))
      if i==3: exp1SigHi.append(float(l.limit))
      if i==4: exp2SigHi.append(float(l.limit))
      if i==5: obs.append(float(l.limit))    

  print 'masses:', xAxis
  print 'obs:',obs
  print 'exp:',exp
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

  oneSigma.SetFillColor(kGreen)

  twoSigma.SetFillColor(kYellow)

  expected.SetMarkerColor(kBlack)
  expected.SetMarkerStyle(kFullCircle)
  expected.SetMarkerSize(1.5)
  expected.SetLineColor(kBlack)
  expected.SetLineWidth(2)
  expected.SetLineStyle(2)

  observed.SetLineWidth(2)

  mg.Add(twoSigma)
  mg.Add(oneSigma)
  mg.Add(expected)
  #mg.Add(observed)

  mg.Draw('AL3')
  mg.GetXaxis().SetTitle('m_{H} (GeV)')
  mg.GetYaxis().SetTitle('95% CL limit on #sigma/#sigma_{SM}')
  mg.GetXaxis().SetLimits(massList[0],massList[-1]);
  gPad.RedrawAxis()

  prelim = TLatex(0.15,0.95, "CMS Preliminary")
  prelim.SetNDC();
  prelim.SetTextSize(0.05);
  prelim.Draw();

  sqrt = TLatex(0.25,0.85, "#sqrt{s} = 8 Tev; #it{L_{int}} = 19.7 fb^{-1}")
  sqrt.SetNDC();
  sqrt.SetTextSize(0.04);
  sqrt.Draw();
  
  
  selection = TLatex()
  if "mu" in s:
    selection = TLatex(0.70,0.95, "#mu#mu selection");
    selection.SetTextColor(kBlue-3);
  elif "el" in s:
    selection = TLatex(0.70,0.95, "#font[32]{ee} selection");
    selection.SetTextColor(kGreen-3);

  selection.Draw()
  selection.SetNDC();
  prelim.SetTextSize(0.03);
  selection.Draw();            

  c.SaveAs(plotBase+'Limits.png')


  out.cd()
  observed.Write("observed")
  expected.Write("expected")
  oneSigma.Write("oneSigma")
