#!/usr/bin/env python
import os,sys

from ROOT import *
gROOT.SetBatch()

gROOT.ProcessLine(".L ~/tdrstyle.C")
setTDRStyle()
gStyle.SetOptTitle(0)

print len(sys.argv), sys.argv
if len(sys.argv) < 1:
    sys.exit()

nominal = 'v42-pre-app'
 
f0 = TFile(nominal+"/limit-data.root")
g0 = f0.Get("expected")
f1 = TFile(nominal+"-noSyst/limit-data.root")
g1 = f1.Get("expected")
f2 = TFile(nominal+"-br/limit-data.root")
g2 = f1.Get("expected")

g0.UseCurrentStyle()
g0.SetMarkerColor(kBlack)
g0.SetLineWidth(2)
#g0.SetLineStyle(2)

g1.UseCurrentStyle()
g1.SetLineColor(kRed+1)
g1.SetLineWidth(2)
#g1.SetLineStyle(2)

g2.UseCurrentStyle()
g2.SetLineColor(kGreen+1)
g2.SetLineWidth(2)

mg = TMultiGraph()
mg.SetTitle('')
mg.Add(g0)
mg.Add(g1)
#mg.Add(g2)


leg = TLegend(0.30,0.7,0.64,0.90);
leg.AddEntry(g0,"Nominal", "l")
leg.AddEntry(g1,"No sytematics", "l")
leg.AddEntry(g2,"2x syst on BR", "l")

f = {}
g = {}
for i,s in enumerate(sys.argv[1:]):
    print i,s

    f[s] = TFile(s+"/limit-data.root") 
    g[s] = f[s].Get("expected")
    g[s].UseCurrentStyle()
    g[s].SetLineWidth(2)
    g[s].SetLineStyle(2)

    mg.Add(g[s])
    if "v35-mumu" in s:
        leg.AddEntry(g[s],"DoubleMu data","l")
    elif "v35-mu-oldiso" in s:
        leg.AddEntry(g[s],"Muon isolation","l")
        g[s].SetLineColor(4)
    elif "v38-mugamma-tight" in s:
        leg.AddEntry(g[s],"Tight Photon ID","l")
        g[s].SetLineColor(7)
        
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



