import sys,os

from ROOT import *
import datetime

import config as c
import utils as u

F0 = 6
FNOB =5
FB = 8
FMVA = 25

makeMvaTrees = 0
plotSpecial = 1
plotMVA = 0
plotMet = 0
plotJet = 0
plotLepton = 0
plotDiLepton = 0
plotMisc = 0

makePuWeights = 0

def makePlots(sel=1, dir="./", hPath="v00"):
    gROOT.ProcessLine(".L ../data/tdrstyle.C");
    setTDRStyle()
    
    gStyle.SetPadGridX(1)
    gStyle.SetPadGridY(1);
    gStyle.SetHistLineWidth(2)
    gStyle.SetLegendBorderSize(0);
    TH1.SetDefaultSumw2(kTRUE);
    
    if sel==1:
        thissel = "muon"
        imgpath = dir+hPath+"/muon/" 
    if sel==2:
        thissel = "electron"
        imgpath = dir+hPath+"/electron/"
    print "Making plots from", hPath
    
    #topbg = ["tW","tbarW"]
    bg   = ["tW","tbarW","ttbar","WW", "WZ","ZZ","DYjets"]
    sig1 = ["ggHZZ125","VBFHZZ200","ggHZZ300","ggHZZ350","ggHZZ400"]
    sig2 = ["VBFHZZ125","VBFHZZ200","VBFHZZ300","VBFHZZ350","VBFHZZ400"]
    sig3 = ["ggHWW125","ggHWW250","ggHWW300","ggHWW350","ggHWW400"]

    sig4 = [""]
    
    li_sig1  = {}
    li_sig2  = {}
    li_sig3  = {}
    li_sig4  = {}
    li_ov    = {}

    li_allbg = {}

    for a in bg:
        #print a
        if a=='ttbar':
            f = TFile(hPath+"/m_ttbar_"+thissel+".root", "OPEN")
        elif a=='DYjets':
            f = TFile(hPath+"/m_DYjets_"+thissel+".root", "OPEN")
        elif a=='vbfZ':
            f = TFile(hPath+"/m_vbfZ_"+thissel+".root", "OPEN")
        else:
            f = TFile(hPath+"/"+thissel+"/hhhh_"+a+"_1.root", "OPEN")
        li_allbg[a] = f
            
    for a in sig1:
        #print a
        f = TFile(hPath+"/"+thissel+"/hhhh_"+a+"_1.root", "OPEN")
        #f.Print()
        li_sig1[a] =f
        if a in ["ggHZZ125"]:
            li_ov[a] = f
    for a in sig2:
        #print a
        f = TFile(hPath+"/"+thissel+"/hhhh_"+a+"_1.root", "OPEN")
        #f.Print()
        li_sig2[a] =f
        if a in ["VBFHZZ125"]:
            li_ov[a] = f

    f_Data =  TFile(hPath+"/m_Data_"+thissel+".root", "OPEN")
    f_Data = None

    #print  li_allbg
    #li_ov["DATA"] = f_Data

    c1 = TCanvas("c1","small canvas",600,500);

    print "\n\n ******** Make the Yield table ******** \n"

    #u.printYields(li_topbg, li_bg, li_sig1, li_sig2, li_sig3, f_Data, sel, "yields_"+thissel+".html", "VBFZ")
    u.printYields(li_allbg, li_sig1, li_sig2, li_sig3, f_Data, sel, "yields_"+thissel+".html", "HZZ125")
    print "\n **** End of Yield table ****"

    if makeMvaTrees:
        mvaInputsDir = "../mva/mvaInputs_"+hPath+"/"+thissel+"/"
        if not os.path.exists(mvaInputsDir):
            os.makedirs(mvaInputsDir)
        u.makeMvaTrees(mvaInputsDir,li_allbg, li_sig1, li_sig2,sel)
        exit(0)


    if plotMVA:
        sdr = "mvaPresel/"
        u.drawMultiPlot(imgpath+sdr+"mva01_BDT_h125_njet0", "nJets=0","BDT discriminator", "mva_discr_mh0_0", 1, 0.1, 1e5, 0.1,3.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"mva02_BDT_h125_njet1", "nJets=1","BDT discriminator", "mva_discr_mh0_1", 1, 0.1, 1e5, 0.1,3.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"mva03_BDT_h125_njet2", "nJets>1","BDT discriminator", "mva_discr_mh0_2", 1, 0.1, 1e5, 0.1,3.9, li_ov, li_allbg, sel)
        #u.drawMultiPlot(imgpath+sdr+"mva01", "nJets=0","BDT discriminator", "mva_discr_0", 0, 0.1, 80, 0.1,3.9, li_ov, li_allbg, sel)
        #u.drawMultiPlot(imgpath+sdr+"mva02", "nJets=1","BDT discriminator", "mva_discr_1", 0, 0.1, 500, 0.1,3.9, li_ov, li_allbg, sel)
        #u.drawMultiPlot(imgpath+sdr+"mva03", "nJets>1","BDT discriminator", "mva_discr_2", 0, 0.1, 300, 0.1,3.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"mva04", "nJets=0","M(ll)", "di_mass_"+str(FMVA), 0, 0.1, 700, 0.3,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"mva05", "nJets=0","qT(ll)", "di_qt_"+str(FMVA), 0, 0.1, 700, 0.3,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"mva06", "nJets=0","pfMet", "met1_et_"+str(FMVA), 0, 0.1, 700, 0.3,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"mva07", "nJets=0","MT", "mt2_"+str(FMVA), 0, 0.1, 700, 0.3,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"mva08", "","Leading Lepton pt", "l1_pt_"+str(FMVA), 0, 0.1, 700, 0.3,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"mva09", "","Trailing Lepton pt", "l2_pt_"+str(FMVA), 0, 0.1, 700, 0.3,1.9, li_ov, li_allbg, sel)

    if plotSpecial:
        sdr = "Special/"
        u.drawMultiPlot(imgpath+sdr+"di01", "","M(ll)", "di_mass_ext_2", 1, 0.1, 1e10, 0.5,1.49, li_ov, li_allbg, sel)

        '''
        print "* Making pu weights *"
        pileup2011 = TFile("pileup2011.root","open")
        pu_da = pileup2011.Get("pileup").Clone()
        pu_da.Sumw2()
        pu_mc = li_allbg["DYjets"].Get("Andrey/evt_pu_1").Clone()
        sample = li_allbg["DYjets"].Get("Andrey/evt_byCut").GetTitle()
        print "\t\t from sample",sample 
        pu_da.Print()
        pu_mc.Print()
        c3 = TCanvas("c3","small canvas",600,600);
        pu_mc_norm = pu_mc.DrawNormalized("hist")
        pu_da_norm = pu_da.DrawNormalized("same hist")
        
        c3.SaveAs(imgpath+"Special/pileup.png")
        pu_da_norm.Divide(pu_mc_norm)
        h1_PU2011 =  pu_da_norm.Clone()
        h1_PU2011.Draw("hist")
        c3.SaveAs(imgpath+"Special/pileup_weight.png")
        if makePuWeights:
            f = TFile("../data/puReweight.root","recreate")
            f.cd()
            h1_PU2011.Write("h1_PU2011")
            f.Close()
        '''
        
        '''
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
        
         '''
    if plotDiLepton:
        sdr = "dilepron/"
        u.drawMultiPlot(imgpath+sdr+"di01", "","M(ll)", "di_mass_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"di02", "","Leptons in Barrel, |#eta|<1.444, M(ll)", "di_mass_EB_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"di03", "","Leptons in Endcap, |#eta|>1.566, M(ll)", "di_mass_EE_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"di04", "","Leptons in EB/EE, mixed, M(ll)", "di_mass_EX_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"di05T", "","q_{T} (di-lepton p_{T})", "di_qt_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"di06", "","Di-lepton Eta", "di_eta_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"di07", "","dPhi(Di-lep, Met)", "di_dPhiMet_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)

    if plotLepton:
        sdr = "Lepton/"
        u.drawMultiPlot(imgpath+sdr+"l01", "","Leading Lepton pt", "l1_pt_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"l02", "","Trailing Lepton pt", "l2_pt_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"l03", "","Leading Lepton eta", "l1_eta_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"l04", "","Trailing Lepton eta", "l2_eta_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"l05", "","Leading Lepton phi", "l1_phi_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"l06", "","Trailing Lepton phi", "l2_phi_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)

        u.drawMultiPlot(imgpath+sdr+"l07_logAngle", "","log(#pi - Angle(lep1, lep2))", "l0_angleLog_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"l08", "","Angle(lep1, lep2)", "l0_angle_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"l09", "","dPhi(lep1, lep2)", "l0_dPhi_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"l10", "","dEta(lep1, lep2)", "l0_dEta_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"l11", "","dR(lep1, lep2)", "l0_dR_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"l12", "","ratio p_{T2}/p_{T1}", "l0_ptRatio_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        

    if plotJet:
        sdr = "Jet/"
        u.drawMultiPlot(imgpath+sdr+"j01", "","N jets pt>30 | #eta|<4.9", "jet_N_"+str(F0), 1, 0.1, 1e6, 0,2.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"j02", "","N jets pt>30 | #eta|<2.4", "jet_N24_"+str(F0), 1, 0.1, 1e6, 0,2.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"j03", "","#Delta#phi(MET, clos jet), p_{T}>30, |#eta|<4.8", "met1_dPhiClosJet1_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"j04", "No b-veto","N b-jets", "jet_b_N_"+str(FNOB), 1, 0.1, 1e6, 0,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"j05", "","pt of leading jet", "jet_pt1_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"j06", "","eta of leading jet", "jet_eta1_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"j07", "","dR(Lead jet,lep 1)", "jet_dRlep1_"+str(F0), 1, 0.1, 1e6, 0.5,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"j08", "","dR(Lead jet,lep 2)", "jet_dRlep2_"+str(F0), 1, 0.1, 1e6, 0.5,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"j09", "","M(jet1, jet2), (GeV)", "jet_diM_"+str(F0), 1, 0.01, 1e5, 0.1,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"j10", "","dEta(jet1, jet2)", "jet_deltaEta_"+str(F0), 1, 0.01, 1e5, 0.1,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"j11", "","#eta* = #eta_{Z} - 0.5(#eta_{j1} + #eta_{j2})", "jet_zeppZ_"+str(F0), 1, 0.01, 1e5, 0.1,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"j12", "","y* = y_{Z} - 0.5(y_{j1} + y_{j2})", "jet_zeppZy_"+str(F0), 1, 0.01, 1e5, 0.1,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"j13", "","N jets pt>15 | #eta|<4.9", "jet_N15_"+str(F0), 1, 0.01, 1e5, 0.1,1.9, li_ov, li_allbg, sel)

    if plotMet:
        sdr = "Met/"
        u.drawMultiPlot(imgpath+sdr+"m01", "","pfMet, (GeV)", "met1_et_"+str(F0), 1, 0.1, 1e6, 0,2.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"m02", "","pfMet Phi", "met1_phi_"+str(F0), 1, 0.1, 1e6, 0.5,1.49, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"m03", "","MT", "mt_"+str(F0), 1, 0.1, 1e6, 0.5,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"m04", "","Met || q(ll)", "met1_projOnQt_"+str(F0), 1, 0.1, 1e6, 0,2.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"m05", "","Met perp to q(ll)", "met1_perpQt_"+str(F0), 1, 0.1, 1e6, 0,2.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"m06", "","Long Recoil = -(pfMet+Z/G)", "met1_recoil_lg_"+str(F0), 1, 0.1, 1e6, 0.2,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"m07", "","pfMet/q_{T}(ll)", "met1_overQt_"+str(F0), 1, 0.1, 1e6, 0,2.9, li_ov, li_allbg, sel)
        #u.drawMultiPlot(imgpath+sdr+"m06", "","projMET", "met4_et_"+str(F0), 1, 0.1, 1e6, 0,2.9, li_ov, li_allbg, sel)
        #u.drawMultiPlot(imgpath+sdr+"m07", "","ZprojMET", "met5_et_"+str(F0), 1, 0.1, 1e6, 0,2.9, li_ov, li_allbg, sel)
        #u.drawMultiPlot(imgpath+sdr+"m08", "","redMET1", "met6_et_"+str(F0), 1, 0.1, 1e6, 0,2.9, li_ov, li_allbg, sel)
        #u.drawMultiPlot(imgpath+sdr+"m09", "","redMET2", "met7_et_"+str(F0), 1, 0.1, 1e6, 0,2.9, li_ov, li_allbg, sel)

    if plotMisc:
        sdr = "Misc/"
        u.drawMultiPlot(imgpath+sdr+"mis01", "","evts cut by cut", "evt_byCut", 1, 0.1, 1e6, 0,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"mis02", "","nVtx total", "vtx_nPV_tot_"+str(F0), 0, 0.1, 1e4, 0,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"mis03", "","nVtx raw", "vtx_nPV_raw_"+str(F0), 0, 0.1, 1e4, 0,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"mis04", "","nVtx reweighted", "vtx_nPV_weight_"+str(F0), 0, 0.1, 1e4, 0,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"mis05", "","vtx 1 ndof", "vtx_ndof1_"+str(F0), 1, 0.1, 1e6, 0,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"mis06", "","vtx 2 ndof", "vtx_ndof2_"+str(F0), 1, 0.1, 1e6, 0,1.9, li_ov, li_allbg, sel)
        u.drawMultiPlot(imgpath+sdr+"mis07", "","True pile-up", "evt_pu_"+str(F0), 0, 0.1, 1500, 0,1.9, li_ov, li_allbg, sel)


    
    #c2.cd()
    #hh = f_Data.Get('Andrey/di_qt_8')
    #hclone = hh.Clone()
    #hclone.Draw()
    
    #c2.SaveAs(imgpath+"/test.png")
