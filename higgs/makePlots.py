import sys,os

from ROOT import *
import datetime

import top_bkg as t
import utils as u

FMVA = 18
makeMvaTrees = 0
makePuWeights = 0
myParams = u.myParams

def makePlots(sel, dir, ver, F0, plotMuons,plotElectrons,plotSpecial,plotMVA,plotMet,plotJet,plotLepton,plotDiLepton,plotMisc):

    subdir = ver+"_cut"+F0
    hPath = '/eos/uscms/store/user/andreypz/batch_output/8TeV/'+ver
    gROOT.ProcessLine(".L ../data/tdrstyle.C");
    setTDRStyle()
    UP = 1e10
    if F0=="2":
        UP=1e8
    elif F0=="4":
        UP=1e7
    elif F0=="5":
        UP=1e6
    elif F0=="6":
        UP=1e6
        
    period = "2012"
    gStyle.SetPadGridX(1)
    gStyle.SetPadGridY(1);
    gStyle.SetHistLineWidth(2)
    gStyle.SetLegendBorderSize(0);
    TH1.SetDefaultSumw2(kTRUE);
    
    if sel==1:
        thissel = "muon"
        imgpath = dir+"/"+thissel+"/"
    if sel==2:
        thissel = "electron"
        imgpath = dir+"/"+thissel+"/"
    print "Making plots from", hPath, ver
    
    #topbg = ["tW","tbarW"]
    #bg   = ["tW","tbarW","ttbar","WW", "WZ","ZZ","DYjets"]
    #bg   = ["tW","tbarW","ttbar", "WZJetsTo3LNu","ZZJetsTo2L2Nu","DYjets", "DYjets10"]
    bg   = ["tW","tbarW","ttbar","WWJetsTo2L2Nu", "WZJetsTo3LNu","ZZJetsTo2L2Nu","DYjets", "DYjets10"]
    
    sig1 = ["ggHZZ125","ggHZZ200"]#,"ggHZZ300","ggHZZ350","ggHZZ400"]
    sig2 = ["VBFHZZ125","VBFHZZ200"]#,"VBFHZZ300","VBFHZZ350","VBFHZZ400"]
    sig3 = ["ggHWW125","ggHWW200"]#,"ggHWW300","ggHWW350","ggHWW400"]

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
            f = TFile(hPath+"/m_ttbar_"+thissel+"_"+period+".root", "OPEN")
        elif a=='DYjets':
            f = TFile(hPath+"/m_DYjets_"+thissel+"_"+period+".root", "OPEN")
        elif a=='DYjets10':
            f = TFile(hPath+"/m_DYjets10_"+thissel+"_"+period+".root", "OPEN")
        elif a=='vbfZ':
            f = TFile(hPath+"/m_vbfZ_"+thissel+"_"+period+".root", "OPEN")
        else:
            f = TFile(hPath+"/"+thissel+"_"+period+"/hhhh_"+a+"_1.root", "OPEN")
        li_allbg[a] = f
            
    for a in sig1:
        #print a
        f = TFile(hPath+"/"+thissel+"_"+period+"/hhhh_"+a+"_1.root", "OPEN")
        #f.Print()
        li_sig1[a] =f
        if a in ["ggHZZ200"]:
            li_ov[a] = f
    for a in sig2:
        #print a
        f = TFile(hPath+"/"+thissel+"_"+period+"/hhhh_"+a+"_1.root", "OPEN")
        #f.Print()
        li_sig2[a] =f
        if a in ["VBFHZZ200"]:
            li_ov[a] = f
            
    f_Data =  TFile(hPath+"/m_Data_"+thissel+"_"+period+".root", "OPEN")
    #f_Data = None

    #print  li_allbg
    li_ov["DATA"] = f_Data

    c1 = TCanvas("c1","small canvas",600,500);

    print "\n\n ******** Make the Yield table ******** \n"

    #u.printYields(li_topbg, li_bg, li_sig1, li_sig2, li_sig3, f_Data, sel, "yields_"+thissel+"_"+period+".html", "VBFZ")
    u.printYields(li_allbg, li_sig1, li_sig2, li_sig3, f_Data, sel, "yields_"+thissel+"_"+period+".html", "HZZ125", "notscaled")
    print "\n **** End of Yield table ****"


    if makeMvaTrees:
        mvaInputsDir = "../mva/mvaInputs_"+hPath+"/"+thissel+"_"+period+"/"
        if not os.path.exists(mvaInputsDir):
            os.makedirs(mvaInputsDir)
        u.makeMvaTrees(mvaInputsDir,li_allbg, li_sig1, li_sig2,sel)
        exit(0)


    if plotMuons and sel==1:
        path = imgpath+"Muons/"
        u.drawMultiPlot(path+"mu01", "","NumberOfValidMuonHits", "mu_NumberOfValidMuonHits", 1, 0.1,1e9, 0,1.9, li_ov, li_allbg, sel, 0)
        u.drawMultiPlot(path+"mu02", "","NumberOfValidTrackerHits", "mu_NumberOfValidTrackerHits", 1, 0.1,1e9, 0,1.9, li_ov, li_allbg, sel, 0)
        u.drawMultiPlot(path+"mu03", "","NumberOfValidPixelHits", "mu_NumberOfValidPixelHits", 1, 0.1,1e9, 0,1.9, li_ov, li_allbg, sel, 0)
        u.drawMultiPlot(path+"mu04", "","NormalizedChi2", "mu_NormalizedChi2", 1, 0.1,1e9, 0,1.9, li_ov, li_allbg, sel, 0)
        u.drawMultiPlot(path+"mu05", "","dxy", "mu_dxy", 1, 0.1,1e9, 0,1.9, li_ov, li_allbg, sel, 0)
        u.drawMultiPlot(path+"mu06", "","isolation", "mu_iso", 1, 0.1,1e9, 0,1.9, li_ov, li_allbg, sel, 0)
        u.drawMultiPlot(path+"mu07", "","#delta (p_{T})/p_{T}", "mu_ptErrorOverPt", 1, 0.1,1e8, 0,1.9, li_ov, li_allbg, sel, 0)
        u.drawMultiPlot(path+"mu08", "","#delta (p_{T})/p_{T}", "mu_ptErrorOverPt", 1, 0.1,1e8, 0,1.9, li_ov, li_allbg, sel, 0)

        u.drawMultiPlot(path+"mu09", "","triM", "mu_triM", 0, 0.1, 300, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"mu10", "","triM", "mu_triM", 1, 0.1, 1e3, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"mu11", "","triM", "mu_triM_noLowMass", 0, 0.1, 300, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"mu12", "","triM", "mu_triM_noLowMass", 1, 0.1, 1e3, 0.5,1.49, li_ov, li_allbg, sel, 1)

    if plotElectrons and sel==2:
        path = imgpath+"Electrons/"
        u.drawMultiPlot(path+"ele01", "","dR(ele, mu)", "ele_dRmu", 1, 0.1,1e9, 0,1.9, li_ov, li_allbg, sel, 0)
        u.drawMultiPlot(path+"ele02", "","isolation", "ele_iso", 1, 0.1,1e9, 0,1.9, li_ov, li_allbg, sel, 0)
        u.drawMultiPlot(path+"ele03", "","dxy", "ele_dxy", 1, 0.1,1e9, 0,1.9, li_ov, li_allbg, sel, 0)
        u.drawMultiPlot(path+"ele04", "","#delta (p_{T})/p_{T}", "ele_ptErrorOverPt", 1, 0.1,1e9, 0,1.9, li_ov, li_allbg, sel, 0)
        u.drawMultiPlot(path+"ele05", "","e-Mass", "ele_mass", 1, 0.1,1e9, 0,1.9, li_ov, li_allbg, sel, 0)
        u.drawMultiPlot(path+"ele06", "","e-Mass", "ele_mass", 1, 0.1,1e9, 0,1.9, li_ov, li_allbg, sel, 0)
        u.drawMultiPlot(path+"ele07", "","e-Mass", "ele_mass", 1, 0.1,1e9, 0,1.9, li_ov, li_allbg, sel, 0)
        u.drawMultiPlot(path+"ele08", "","e-Mass", "ele_mass", 1, 0.1,1e9, 0,1.9, li_ov, li_allbg, sel, 0)

        u.drawMultiPlot(path+"ele09", "","triM", "ele_triM", 0, 0.1, 300, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"ele10", "","triM", "ele_triM", 1, 0.1, 1e3, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"ele11", "","triM", "ele_triM_noLowMass", 0, 0.1, 300, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"ele12", "","triM", "ele_triM_noLowMass", 1, 0.1, 1e3, 0.5,1.49, li_ov, li_allbg, sel, 1)
       
    if plotMVA:
        path = imgpath+"mvaPresel/"

        u.drawMultiPlot(path+"mva01_BDT_h125_njet0", "nJets=0","BDT discriminator", "mva_discr_mh0_0", 1, 0.1, 1e5, 0.1,3.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"mva02_BDT_h125_njet1", "nJets=1","BDT discriminator", "mva_discr_mh0_1", 1, 0.1, 1e5, 0.1,3.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"mva03_BDT_h125_njet2", "nJets>1","BDT discriminator", "mva_discr_mh0_2", 1, 0.1, 1e5, 0.1,3.9, li_ov, li_allbg, sel, 1)
        #u.drawMultiPlot(path+"mva01", "nJets=0","BDT discriminator", "mva_discr_0", 0, 0.1, 80, 0.1,3.9, li_ov, li_allbg, sel, 1)
        #u.drawMultiPlot(path+"mva02", "nJets=1","BDT discriminator", "mva_discr_1", 0, 0.1, 500, 0.1,3.9, li_ov, li_allbg, sel, 1)
        #u.drawMultiPlot(path+"mva03", "nJets>1","BDT discriminator", "mva_discr_2", 0, 0.1, 300, 0.1,3.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"mva04", "nJets=0","M(ll)", "di_mass_"+str(FMVA), 0, 0.1, 700, 0.3,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"mva05", "nJets=0","qT(ll)", "di_qt_"+str(FMVA), 0, 0.1, 700, 0.3,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"mva06", "nJets=0","pfMet", "met1_et_"+str(FMVA), 0, 0.1, 700, 0.3,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"mva07", "nJets=0","MT", "mt_"+str(FMVA), 0, 0.1, 700, 0.3,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"mva08", "","Leading Lepton pt", "l1_pt_"+str(FMVA), 0, 0.1, 700, 0.3,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"mva09", "","Trailing Lepton pt", "l2_pt_"+str(FMVA), 0, 0.1, 700, 0.3,1.9, li_ov, li_allbg, sel, 1)

    if plotSpecial:

        path = imgpath+"Special/"
        print path

        u.drawMultiPlot(path+"sp01", "","M(ll)", "di_mass_ext_2", 1, 0.1, 1e8, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"sp02", "","M(ll)", "di_mass_ext_2", 1, 0.1, 1e8, 0.5,1.49, li_ov, li_allbg, sel, 1)

        '''
        if sel==1:
            h = f_Data.Get("Muons/mu_triM")
            h.Draw("hist")
            c1.SaveAs(path+"sp05_mu_triM.png")
            h = f_Data.Get("Muons/mu_triM_noLowMass")
            h.Draw("hist")
            c1.SaveAs(path+"sp05_mu_triM_noLowMass.png")
        elif sel==2:

            c1.SaveAs(path+"sp05_ele_triM.png")
        '''
        '''
        lumi = myParams.getLumi(sel) 
        c1.cd()
        h_ttbar_N = li_allbg["ttbar"].Get("Histos/jet_b_N_7")
        h_tW_N    = li_allbg["tW"].Get("Histos/jet_b_N_7")
        h_tbarW_N = li_allbg["tbarW"].Get("Histos/jet_b_N_7")
        u.handleOverflowBinsScaleAndColors(h_ttbar_N, "ttbar", lumi)
        u.handleOverflowBinsScaleAndColors(h_tW_N, "tW", lumi)
        u.handleOverflowBinsScaleAndColors(h_tbarW_N, "tbarW", lumi)

        t.topbg(h_ttbar_N, h_tW_N, h_tbarW_N,  sel, myParams.getS())
        h_ttbar_N.Draw("hist")
        h_tW_N.Add(h_tbarW_N)
        h_tW_N.Draw("hist  same")
    
        c1.SaveAs(path+"sp02.png")

        '''
        
        '''
        presel_ttbar = h_ttbar_N.Integral()
        presel_tW = h_tW_N.Integral()
        for c in range(9,14):
            h_ttbar_c = li_allbg["ttbar"].Get("Histos/jet_b_N_"+str(c))
            h_tW_c    = li_allbg["tW"].Get("Histos/jet_b_N_"+str(c))
            u.handleOverflowBinsScaleAndColors(h_ttbar_c, "ttbar", lumi)
            u.handleOverflowBinsScaleAndColors(h_tW_c, "tW", lumi)
            if h_ttbar_c!=None:
                y = h_ttbar_c.Integral()
                print c, "R = ", y/presel_ttbar
            if h_tW_c!=None:
                y = h_tW_c.Integral()
                print c, "R = ", y/presel_tW
        '''
        '''
        gen_b_Nb30 = {}
        gen_b_Nb30["ttbar"] = li_allbg["ttbar"].Get("gen_b/b_Nb30") 
        gen_b_Nb30["tW"]    = li_allbg["tW"].Get("gen_b/b_Nb30") 
        gen_b_Nb30["tbarW"] = li_allbg["tbarW"].Get("gen_b/b_Nb30") 
        gen_b_fromW = li_allbg["ttbar"].Get("gen_b/b_Nb_fromW") 

        ttbar = gen_b_Nb30["ttbar"].DrawNormalized("hist")
        tW = gen_b_Nb30["tW"].DrawNormalized("same hist")
        ttbar.SetLineColor(kBlue)
        fromW = gen_b_fromW.DrawNormalized("same hist")
        fromW.SetLineColor(kRed)
        ttbar.SetMaximum(0.9)
        c1.SaveAs(path+"sp03.png")
                
        t.CalcAccept(gen_b_Nb30)

        '''

        '''
        # Runs and trigger prescales
        runs = f_Data.Get("DataInfo/run_events_2")
        runs.Draw()
        c1.SaveAs(path+"sp04_runs.png")

        runs = f_Data.Get("DataInfo/run_prescale")
        runs.Draw()
        c1.SaveAs(path+"sp05_prescale.png")

        print "* Making pu weights *"
        pileup2012 = TFile("pileup2012AB.root","open")
        
        pu_da = pileup2012.Get("pileup").Clone()
        pu_da.Sumw2()
        pu_mc  = li_allbg["DYjets"].Get("Histos/evt_puTrue_1").Clone()
        pu_mc2 = li_allbg["DYjets"].Get("Histos/evt_puTrue_2").Clone()
        #pu_mc3 = li_allbg["DYjets"].Get("Histos/evt_pu_1").Clone()
        #pu_da.Print()
        #pu_mc.Print()
        c3 = TCanvas("c3","small canvas",600,600);
        pu_da_norm  = pu_da.DrawNormalized("hist")
        pu_mc_norm  = pu_mc.DrawNormalized("same hist")
        pu_mc_norm2 = pu_mc2.DrawNormalized("same hist")
        #pu_mc_norm3 = pu_mc3.DrawNormalized("same hist")
        pu_mc_norm.SetLineColor(kRed+1)
        #pu_mc_norm3.SetLineColor(kGreen+2)
        pu_da_norm.SetTitle(";True pile-up; normalized")
        leg = TLegend(0.55,0.73,0.85,0.89);
        leg.SetTextSize(0.04);
        leg.AddEntry(pu_da_norm, "data", "f");
        leg.AddEntry(pu_mc_norm, "MC true", "f");
        #leg.AddEntry(pu_mc_norm3, "MC observed", "f");
        leg.AddEntry(pu_mc_norm2, "MC reweighted", "f");
        leg.SetFillColor(kWhite);
        leg.Draw()
        
        c3.SaveAs(path+"pileup.png")
        
        pu_da_norm.Divide(pu_mc_norm)
        h1_PU2012 =  pu_da_norm.Clone()
        h1_PU2012.Draw("hist")
        h1_PU2012.SetTitle("weight; nPu true")
        c3.SaveAs(path+"pileup_weight.png")
        
        if makePuWeights:
            f = TFile("../data/puReweight.root","update")
            f.cd()
            h1_PU2012.Write("h1_PU2012")
            f.Close()


        filt_mc = li_allbg["DYjets"].Get("Histos/evt_byCut_filters").Clone()
        filt_da = f_Data.Get("Histos/evt_byCut_filters").Clone()
        filt_mc_norm = filt_mc.DrawNormalized("hist")
        filt_da_norm = filt_da.DrawNormalized("same hist")
        filt_mc_norm.SetTitle(";cut number;fraction of events")
        filt_mc_norm.SetLineColor(kRed+1)
        c3.SaveAs(path+"filters.png")
        '''
        
        '''
        c1.cd()
        for a in ["ggHZZ400","ggHZZ450","ggHZZ500","ggHZZ550","ggHZZ600"]:
            ind = sig1.index(a)
            #print a, ind
            test1 = li_sig1.At(ind).Get("Andrey/higgs_w_mass_"+F0).Clone()
            u.handleOverflowBinsScaleAndColors(test1,a, 0)
            test2 = li_sig1.At(ind).Get("Andrey/higgs_mass_"+F0).Clone()
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

        path = imgpath+"diLepton/"
        u.drawMultiPlot(path+"di01", "","M(ll)", "di_mass_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"di02", "","Leptons in Barrel, |#eta|<1.444, M(ll)", "di_mass_EB_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"di03", "","Leptons in Endcap, |#eta|>1.566, M(ll)", "di_mass_EE_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"di04", "","Leptons in EB/EE, mixed, M(ll)", "di_mass_EX_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"di05", "","q_{T} (di-lepton p_{T})", "di_qt_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"di06", "","Di-lepton Eta", "di_eta_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"di07", "","dPhi(Di-lep, Met)", "di_dPhiMet_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)

    if plotLepton:
        path = imgpath+"Lepton/"
        u.drawMultiPlot(path+"l01", "","Leading Lepton pt", "l1_pt_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"l02", "","Trailing Lepton pt", "l2_pt_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"l03", "","Leading Lepton eta", "l1_eta_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"l04", "","Trailing Lepton eta", "l2_eta_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"l05", "","Leading Lepton phi", "l1_phi_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"l06", "","Trailing Lepton phi", "l2_phi_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)

        u.drawMultiPlot(path+"l07_logAngle", "","log(#pi - Angle(lep1, lep2))", "l0_angleLog_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"l08", "","Angle(lep1, lep2)", "l0_angle_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"l09", "","dPhi(lep1, lep2)", "l0_dPhi_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"l10", "","dEta(lep1, lep2)", "l0_dEta_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"l11", "","dR(lep1, lep2)", "l0_dR_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"l12", "","ratio p_{T2}/p_{T1}", "l0_ptRatio_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        

    if plotJet:
        path = imgpath+"Jet/"
        u.drawMultiPlot(path+"j01", "","N jets p_{T}>30 | #eta|<4.9", "jet_N_"+F0, 1, 0.1, UP, 0,2.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"j02", "","N jets p_{T}>30 | #eta|<2.4", "jet_N24_"+F0, 1, 0.1, UP, 0,2.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"j03", "","#Delta#phi(MET, closest jet)", "met1_dPhiClosJet1_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"j04", "","M(jet1, jet2), GeV", "jet_diM_"+F0, 1, 0.01, 1e5, 0.1,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"j05", "","dEta(jet1, jet2)", "jet_deltaEta_"+F0, 1, 0.01, 1e5, 0.1,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"j06", "","#eta* = #eta_{Z} - 0.5(#eta_{j1} + #eta_{j2})", "jet_zeppZ_"+F0, 1, 0.01, 1e5, 0.1,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"j07", "","y* = y_{Z} - 0.5(y_{j1} + y_{j2})", "jet_zeppZy_"+F0, 1, 0.01, 1e5, 0.1,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"j08", "","N jets p_{T}>15 | #eta|<4.9", "jet_N15_"+F0, 1, 0.01, 1e5, 0.1,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"j09", "","pt of  leading jet", "jet_pt1_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"j10", "","pt of second jet", "jet_pt2_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"j11", "","eta of leading jet", "jet_eta1_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"j12", "","eta of second jet", "jet_eta2_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"j13", "","dR(Lead jet,lep 1)", "jet_dRlep1_"+F0, 1, 0.1, UP, 0.5,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"j14", "","dR(Lead jet,lep 2)", "jet_dRlep2_"+F0, 1, 0.1, UP, 0.5,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"j15", "","dR(Lead jet, di-Lepton)", "jet_dRDiLept_"+F0, 1, 0.1, UP, 0.5,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"j16", "No b-veto","N b-jets", "jet_b_N_"+F0, 1, 0.1, UP, 0,1.9, li_ov, li_allbg, sel, 1)

    if plotMet:
        path = imgpath+"Met/"
        u.drawMultiPlot(path+"m01", "","pf E_{T}^{miss}, GeV", "met1_et_"+F0, 1, 0.1, UP, 0,2.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"m02", "","pfMet Phi", "met1_phi_"+F0, 1, 0.1, UP, 0.5,1.49, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"m03", "","MT", "mt_"+F0, 1, 0.1, UP, 0.5,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"m04", "","Met || q(ll)", "met1_projOnQt_"+F0, 1, 0.1, UP, 0,2.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"m05", "","Met perp to q(ll)", "met1_perpQt_"+F0, 1, 0.1, UP, 0,2.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"m06", "","Long Recoil = -(pfMet+Z/G)", "met1_recoil_lg_"+F0, 1, 0.1, UP, 0.2,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"m07", "","pfMet/q_{T}(ll)", "met1_overQt_"+F0, 1, 0.1, UP, 0,2.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"m08", "","Sum Et", "met1_SumEt_"+F0, 1, 0.1, UP, 0,2.9, li_ov, li_allbg, sel, 1)
        #u.drawMultiPlot(path+"m06", "","projMET", "met4_et_"+F0, 1, 0.1, UP, 0,2.9, li_ov, li_allbg, sel, 1)
        #u.drawMultiPlot(path+"m07", "","ZprojMET", "met5_et_"+F0, 1, 0.1, UP, 0,2.9, li_ov, li_allbg, sel, 1)
        #u.drawMultiPlot(path+"m08", "","redMET1", "met6_et_"+F0, 1, 0.1, UP, 0,2.9, li_ov, li_allbg, sel, 1)
        #u.drawMultiPlot(path+"m09", "","redMET2", "met7_et_"+F0, 1, 0.1, UP, 0,2.9, li_ov, li_allbg, sel, 1)

    if plotMisc:
        path = imgpath+"Misc/"
        u.drawMultiPlot(path+"mis01", "","evts cut by cut", "evt_byCut", 1, 0.1, UP, 0,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"mis02", "","nVtx total", "vtx_nPV_tot_"+F0, 0, 0.1, 3e5, 0,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"mis03", "","nVtx raw", "vtx_nPV_raw_"+F0, 0, 0.1, 3e5, 0,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"mis04", "","nVtx reweighted", "vtx_nPV_weight_"+F0, 0, 0.1, 3e5, 0,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"mis05", "","vtx 1 ndof", "vtx_ndof1_"+F0, 1, 0.1, UP, 0,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"mis06", "","vtx 2 ndof", "vtx_ndof2_"+F0, 1, 0.1, UP, 0,1.9, li_ov, li_allbg, sel, 1)
        u.drawMultiPlot(path+"mis07", "","Number of photons", "ph_nGamma_"+F0, 1, 0.1, UP, 0,1.9, li_ov, li_allbg, sel, 1)


    
    #c2.cd()
    #hh = f_Data.Get('Andrey/di_qt_8')
    #hclone = hh.Clone()
    #hclone.Draw()
    
    #c2.SaveAs(imgpath+"/test.png")
