#!/usr/bin/env python
import os,sys

import numpy as np
from ROOT import *
gROOT.SetBatch()
sys.path.append("../zgamma")
import utils as u
from optparse import OptionParser
parser = OptionParser(usage="usage: %prog  [options]")
parser.add_option("--br", dest="br",  action="store_true", default=False, help="Do the limit on BR*cs instead of the mu")
parser.add_option("--mll",dest="mll", action="store_true", default=False, help="Do the limit on BR*cs in bins of mll")
parser.add_option("--divide",dest="divide", action="store_true", default=False, help="Divide by the binn-size (for --mll option)")
parser.add_option("--cat",dest="cat", default='EB', help="Category: EB, EE, mll50, etc")
(opt, args) = parser.parse_args()

import ConfigParser as cp
cf = cp.ConfigParser()
cf.read('config.cfg')

gStyle.SetOptTitle(0)

s = cf.get("path","ver")
#s = "comboCardsMing/"
#s = "MingDataCards/Dalitz_electron/"

method = 'Asymptotic'
#method = 'ProfileLikelihood'
out = TFile(s+"/limit-data-cat"+opt.cat+".root","recreate")

plotBase = cf.get("path","htmlbase")+'/html/zgamma/dalitz/fits-'+s
u.createDir(plotBase+'/Limits')

massList  = [float(a) for a in (cf.get("fits","massList-more")).split(',')]
doBlind   = int(cf.get("fits","blind"))
doObs = not doBlind
mllBins = u.mllBins()

# print massList

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
  fileListTmp = os.listdir(s)
  fileList = filter(lambda fileName: 'higgsCombineTest.'+method in fileName,fileListTmp)
  #print fileList
  SMtot = 0.757

  if opt.mll:
    for mbin in ['m1','m2','m3','m4','m5','m6','m7']:
      #center = mllBins[int(mbin[1])-1][0]+ 0.5*(mllBins[int(mbin[1])][0]-mllBins[int(mbin[1])-1][0])
      r1 = mllBins[int(mbin[1])-1][0]
      r2 = mllBins[int(mbin[1])][0]
      center = 0.5*(r2+r1)
      print 'mbin and center',mbin, mllBins[int(mbin[1])], center

      thisFile = "higgsCombineTest.Asymptotic.mH125.0_cat"+mbin+".root"
      # print thisFile
      xAxis.append(float(center))
      xErr.append(float(0.5*(r2-r1)))
      dmScale = 1
      if opt.divide:
        dmScale = r2-r1

      SM.append(10*SMtot*mllBins[int(mbin[1])][1]/dmScale)
      print '10x SM prediction:', SM[-1]

      f = TFile(s+"/"+thisFile,"open")
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
      print mass
      thisFile = filter(lambda fileName: str(mass)+'_cat'+opt.cat+'.root' in fileName,fileList)[0]
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
  print '1 sig hi:', exp1SigHiErr
  print '2 sig hi:', exp2SigHiErr

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

  observed.SetLineWidth(2)
  observed.SetMarkerStyle(kFullCircle)

  if opt.mll:
    standardM.SetFillColor(kRed+1)
    standardM.SetLineColor(kRed+1)
    standardM.SetLineWidth(2)
    # standardM.SetMarkerStyle(24)

    mg.Add(twoSigma,'P2')
    mg.Add(oneSigma,'P2')
    mg.Add(expected,'P')
  else:
    mg.Add(twoSigma)
    mg.Add(oneSigma)
    mg.Add(expected)

  if doObs:
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
    mg.GetYaxis().SetTitle('#sigma(pp #rightarrow A)#timesBR(A#rightarrow#mu#mu#gamma) (fb)')
    mg.SetMaximum(21)
    if opt.mll:
      mg.SetMaximum(12)
      if opt.divide:
        mg.GetYaxis().SetTitle('#frac{1}{#Deltam} #sigma(pp#rightarrowA)#timesBR(A#rightarrow#mu#mu#gamma) (fb/GeV)')
        mg.SetMaximum(200)
        mg.SetMinimum(0.03)
        c.SetLogy()
  else:
    mg.GetYaxis().SetTitle('95% CL limit on #sigma#timesBR/#sigma_{SM}#times BR_{SM}')
    if opt.cat in ['EB','Combo']:
      mg.SetMaximum(31)
      if 'el' in s:
        mg.SetMaximum(42)
    elif opt.cat in ['EE','mll50']:
      mg.SetMaximum(150)
  gPad.RedrawAxis()

  lat = TLatex()
  lat.SetNDC()
  lat.SetTextSize(0.035)

  if opt.mll:
    lat.DrawLatex(0.18,0.95, 'A #rightarrow#gamma*#gamma#rightarrow #mu#mu#gamma')
    lat.SetTextSize(0.05)
    lat.DrawLatex(0.40,0.95,'m_{A} = 125 GeV')
    lat.SetTextSize(0.035)
  else:
    if "mu" in s:
      lat.DrawLatex(0.18,0.95, 'H #rightarrow#gamma*#gamma#rightarrow #mu#mu#gamma  CAT: '+opt.cat)
    elif "el" in s:
      lat.DrawLatex(0.18,0.95, 'H #rightarrow#gamma*#gamma#rightarrow #font[32]{ee}#gamma  CAT: '+opt.cat)
  CMS_lumi(c, 2, 11)

  leg = TLegend(0.50,0.68,0.75,0.88)
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

  if opt.br and not opt.mll:
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
    leg.AddEntry(g,"10 #times SM", "l")
  elif not opt.br:
    l = TLine(120, 1, 150, 1)
    l.SetLineColor(kRed+2)
    l.SetLineStyle(21)
    l.Draw()
  elif opt.mll:
    leg.SetY1NDC(0.64)
    leg.AddEntry(standardM,"10 #times SM Higgs", "lp")


  for e in ['.png', '.pdf']:
    CMS_lumi(c, 2, 11)
    if opt.mll:
      c.SaveAs(plotBase+'/Limits/Limits_Mll'+e)
    else:
      c.SaveAs(plotBase+'/Limits/Limits_cat'+opt.cat+e)

  out.cd()
  observed.Write("observed")
  expected.Write("expected")
  oneSigma.Write("oneSigma")
