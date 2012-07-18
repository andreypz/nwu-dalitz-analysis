#!/usr/bin/env python
import sys,os

from ROOT import *
import shutil
import datetime

import config as c
import utils as u

nargs = len(sys.argv)
print sys.argv[0], nargs

F0 = 6
FNOB = 7
FB = 8

plotSpecial = 1
plotMet = 0
plotJet = 0
plotLepton = 0
plotDiLepton = 0
plotMisc = 0



def makePlots(sel=1, dir="./", hPath="v00"):
    gROOT.ProcessLine(".L ../data/tdrstyle.C");
    setTDRStyle()
    
    gStyle.SetPadGridX(1)
    gStyle.SetPadGridY(1);
    gStyle.SetHistLineWidth(2)
    gStyle.SetLegendBorderSize(0);
    TH1.SetDefaultSumw2(kTRUE);
    
    intLumi = 1000.0

    if sel==1:
        thissel = "muon"
        intLumi = c.Params().getLumi(1)
        imgpath = dir+hPath+"/muon/" 
    if sel==2:
        thissel = "electron"
        intLumi = c.Params().getLumi(2)
        imgpath = dir+hPath+"/electron/"
    print "Making plots from", hPath, "with lumi = ", intLumi
    
    topbg = ["tW","tbarW"]
    bg   = ["WW", "WZ","ZZ", "ttbar","DYjets"]
    sig1 = ["ggHZZ125","ggHZZ200","ggHZZ250","ggHZZ300","ggHZZ350","ggHZZ400","ggHZZ450","ggHZZ500","ggHZZ550","ggHZZ600"]
    sig2 = ["ggHWW125","ggHWW200","ggHWW250","ggHWW300","ggHWW350","ggHWW400","ggHWW450","ggHWW500","ggHWW550","ggHWW600"]
    sig3 = ["VBFHZZ125","VBFHZZ200","VBFHZZ250","VBFHZZ300","VBFHZZ350","VBFHZZ400","VBFHZZ450","VBFHZZ500","VBFHZZ550","VBFHZZ600"]

    li_topbg = TList()
    li_bg    = TList()
    li_sig1  = TList()
    li_ov = TList()

    for a in topbg:
        #print "Creating TList", a
        f = TFile(hPath+"/"+thissel+"/hhhh_"+a+"_1.root", "OPEN")
        li_topbg.Add(f)        
    for a in bg:
        #print a
        if a=='ttbar':
            f = TFile(hPath+"/m_ttbar_"+thissel+".root", "OPEN")
        elif a=='DYjets':
            f = TFile(hPath+"/m_DYjets_"+thissel+".root", "OPEN")
        else:
            f = TFile(hPath+"/"+thissel+"/hhhh_"+a+"_1.root", "OPEN")
        #f.Print()
        #print "Creating TList", a
        li_bg.Add(f)
            
    for a in sig1:
        #print a
        f = TFile(hPath+"/"+thissel+"/hhhh_"+a+"_1.root", "OPEN")
        #f.Print()
        li_sig1.Add(f)

    f_Data =  TFile(hPath+"/m_Data_"+thissel+".root", "OPEN")

    li_ov.Add(f_Data)
    li_ov.Add(li_sig1.At(1))

    c1 = TCanvas("c1","small canvas",600,500);
    c2 = TCanvas("c2","big canvas",600,700);

    if plotSpecial:
        #u.drawMultiPlot(imgpath+"Special/sp01", "","pfMet1/q_{T}", "met1_over_qt_"+str(F0), 1, 0.1, 1e6, 0,2.9, c2, li_ov, li_topbg, li_bg, sel)

        c1.cd()
        for a in ["ggHZZ400","ggHZZ450","ggHZZ500","ggHZZ550","ggHZZ600"]:
            ind = sig1.index(a)
            #print a, ind
            test1 = li_sig1.At(ind).Get("Andrey/higgs_w_mass_"+str(F0)).Clone()
            u.handleOverflowBinsScaleAndColors(test1,a, 0)
            test2 = li_sig1.At(ind).Get("Andrey/higgs_mass_"+str(F0)).Clone()
            u.handleOverflowBinsScaleAndColors(test1,a, 0)

            #test1.Print()
            #test2.Print()
            test1.SetLineColor(kBlue+1);
            test2.SetLineColor(kRed+1);
            test1.Draw("hist");
            test2.Draw("same hist");
            test1.SetTitle(";M_{H}; normalized");

            leg05 = TLegend(0.65,0.73,0.85,0.89);
            leg05.SetTextSize(0.04);
            leg05.AddEntry(test1, "weighted", "f");
            leg05.AddEntry(test2, "no weight", "f");
            leg05.SetFillColor(kWhite);
            leg05.Draw();
            c1.SaveAs(imgpath+"Special/sp"+str(ind)+".png");
            test1.Delete()
            test2.Delete()
        

    if plotMet:
        u.drawMultiPlot(imgpath+"Met/m01", "","pfMet", "met1_et_"+str(F0), 1, 0.1, 1e6, 0,2.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Met/m02", "","pfMet Phi", "met1_phi_"+str(F0), 1, 0.1, 1e6, 0,2.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Met/m03", "","MT", "mt2_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Met/m04", "","pfMet long to di-lepton", "met1_lg_"+str(F0), 1, 0.1, 1e7, 0,2.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Met/m05", "","Long Recoil = -(pfMet+Z/G)", "met1_recoil_lg_"+str(F0), 1, 0.1, 1e6, 0,2.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Met/m06", "","projMET", "met4_et_"+str(F0), 1, 0.1, 1e6, 0,2.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Met/m07", "","ZprojMET", "met5_et_"+str(F0), 1, 0.1, 1e6, 0,2.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Met/m08", "","redMET1", "met6_et_"+str(F0), 1, 0.1, 1e6, 0,2.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Met/m09", "","redMET2", "met7_et_"+str(F0), 1, 0.1, 1e6, 0,2.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Met/m10", "","pfMet1/q_{T}", "met1_over_qt_"+str(F0), 1, 0.1, 1e6, 0,2.9, c2, li_ov, li_topbg, li_bg, sel)

    if plotJet:
        u.drawMultiPlot(imgpath+"Jet/j01", "","N jets", "jet_N_"+str(F0), 1, 0.1, 1e6, 0,2.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Jet/j02", "","#Delta#phi(MET, clos jet), p_{T}>30, |#eta|<4.8", "met1_dPhiClosJet1_"+str(F0), 1, 0.1, 1e6, 0,2.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Jet/j03", "No b-veto, MET>70","N b-jets", "jet_b_N_"+str(FNOB), 1, 0.1, 1e6, 0,2.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Jet/j04", "","pt of leading jet", "jet_pt_"+str(F0), 1, 0.1, 1e6, 0,2.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Jet/j05", "","eta of leading jet", "jet_eta_"+str(F0), 1, 0.1, 1e6, 0,2.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Jet/j06", "","dR(jet,lep 1)", "jet_dRlep1_"+str(F0), 1, 0.1, 1e6, 0,2.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Jet/j07", "","dR(jet,lep 2)", "jet_dRlep2_"+str(F0), 1, 0.1, 1e6, 0,2.9, c2, li_ov, li_topbg, li_bg, sel)

    if plotDiLepton:
        u.drawMultiPlot(imgpath+"diLepton/di01", "","M(ll)", "di_mass_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"diLepton/di02", "","Leptons in Barrel, |#eta|<1.444, M(ll)", "di_mass_EB_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"diLepton/di03", "","Leptons in Barrel, |#eta|>1.566, M(ll)", "di_mass_EE_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"diLepton/di04", "","Leptons in EB/EE, mixed, M(ll)", "di_mass_EX_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"diLepton/di05", "","q_{T} (di-lepton p_{T})", "di_qt_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"diLepton/di06", "","Di-lepton Eta", "di_eta_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"diLepton/di07", "","dPhi(Di-lep, Met)", "di_dPhiMet_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)

    if plotLepton:
        u.drawMultiPlot(imgpath+"Lepton/l01", "","Leading Lepton eta", "l1_eta_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Lepton/l02", "","Trailing Lepton eta", "l2_eta_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Lepton/l03", "","Leading Lepton phi", "l1_phi_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Lepton/l04", "","Trailing Lepton phi", "l2_phi_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Lepton/l05", "","Leading Lepton pt", "l1_pt_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Lepton/l06", "","Trailing Lepton pt", "l2_pt_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)

        u.drawMultiPlot(imgpath+"Lepton/l07", "","dPhi(lep1, lep2)", "l0_dPhi_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Lepton/l08", "","dEta(lep1, lep2)", "l0_dEta_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Lepton/l09", "","dR(lep1, lep2)", "l0_dR_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Lepton/l10", "","ratio p_{T2}/p_{T1}", "l0_ptRatio_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Lepton/l11", "","ratio p_{T2}/p_{T1}", "l0_dR_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
        
    if plotMisc:
        u.drawMultiPlot(imgpath+"Misc/mis01", "","evts cut by cut", "evt_byCut", 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Misc/mis02", "","nVtx total", "vtx_nPV_tot_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Misc/mis03", "","nVtx raw", "vtx_nPV_raw_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Misc/mis04", "","nVtx reweighted", "vtx_nPV_weight_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Misc/mis05", "","vtx 1 ndof", "vtx_ndof_1_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
        u.drawMultiPlot(imgpath+"Misc/mis06", "","vtx 2 ndof", "vtx_ndof_2_"+str(F0), 1, 0.1, 1e6, 0,1.9, c2, li_ov, li_topbg, li_bg, sel)
    
    #c2.cd()
    #hh = f_Data.Get('Andrey/di_qt_8')
    #hclone = hh.Clone()
    #hclone.Draw()
    
    #c2.SaveAs(imgpath+"/test.png")
