#!/usr/bin/env python
import os,sys

from ROOT import *
gROOT.SetBatch()

gROOT.ProcessLine(".L ~/tdrstyle.C")
setTDRStyle()
gStyle.SetOptTitle(0)

print len(sys.argv), sys.argv
if len(sys.argv) < 2:
    sys.exit()

 
f0 = TFile("v35-mu-newiso/limit-data.root")
g0 = f0.Get("expected")
f1 = TFile("v35-mu-newiso-noSyst/limit-data.root")
g1 = f1.Get("expected")

g0.UseCurrentStyle()
g0.SetMarkerColor(kBlack)
g0.SetLineWidth(2)
#g0.SetLineStyle(2)

g1.UseCurrentStyle()
g1.SetLineColor(kRed+1)
g1.SetLineWidth(2)
#g1.SetLineStyle(2)

mg = TMultiGraph()
mg.SetTitle('')
mg.Add(g0)
mg.Add(g1)


leg = TLegend(0.40,0.7,0.70,0.90);
leg.AddEntry(g0,"Nominal", "l")
leg.AddEntry(g1,"No sytematics", "l")

f = {}
g = {}
for i,s in enumerate(sys.argv[1:]):
    print i,s

    f[s] = TFile(s+"/limit-data.root") 
    g[s] = f[s].Get("expected")
    g[s].UseCurrentStyle()
    g[s].SetLineColor(3+i)
    g[s].SetLineWidth(2)
    g[s].SetLineStyle(2)

    mg.Add(g[s])
    if "v35-mumu" in s:
        leg.AddEntry(g[s],"DoubleMu data","l")
    elif "v35-mu-oldiso" in s:
        leg.AddEntry(g[s],"Muon isolation","l")

        
mg.Draw('AL3')
mg.SetMinimum(0)
mg.SetMaximum(20)
mg.GetXaxis().SetTitle('m_{H} (GeV)')
mg.GetYaxis().SetTitle('95% CL limit on #sigma/#sigma_{SM}')
#mg.GetXaxis().SetLimits(massList[0],massList[-1]);
gPad.RedrawAxis()
gPad.Update()
leg.SetTextSize(0.04)
leg.SetFillColor(kWhite)
leg.Draw()

plotBase = '/uscms_data/d2/andreypz/html/zgamma/dalitz/'

c1.SaveAs(plotBase+'Limits-Compare.png')



#sys.path.append("../zgamma")
#import utils as u



