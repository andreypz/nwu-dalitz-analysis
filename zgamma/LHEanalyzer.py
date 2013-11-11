#!/usr/bin/env python

import sys,os
import time
from ROOT import *
import utils as u
import makeHTML as ht
gROOT.SetBatch()


subdir = sys.argv[1]

outpath = '/uscms_data/d2/andreypz/html/zgamma/lhe/'
files={}
#files["mad"] = ['/uscms/home/andreypz/work/MadGraph5/PROC_ANO_HEFT_JGF_v2_1/Events/run_01/unweighted_events.root']
files["mcfm"] = ['/uscms_data/d2/andreypz/lhe_mcfm_hzg_dalitz_lord_fixed_unweighted.lhe.root']
#files["mad"] = ['/uscms_data/d2/andreypz/lhe_mad_hzg_dalitz_unweighted.root']
#files["mad"] = ['/uscms_data/d2/andreypz/lhe_mad_hzg_dalitz_unweighted_Mmu.root']
#files["mad"] = ['/uscms_data/d2/andreypz/lhe_mad_LO_HiggsToMuMuGamma.root']
files["mad"] = ['/uscms_data/d2/andreypz/lhe_mad_hzg5.root']
#files["mad"] = ['/uscms/home/andreypz/lhe_higgs_eegamma_dalitz/heeg_m120.root']

print files
#gSystem.Load("/home/andreypz/workspace/MadGraph5/ExRootAnalysis/lib/libExRootAnalysis.so")
gSystem.Load("/uscms/home/andreypz/work/MadGraph5/ExRootAnalysis/lib/libExRootAnalysis.so")
#gSystem.Load("../plugins/HistManager_cc.so")
gROOT.LoadMacro("../plugins/HistManager.cc+");
gROOT.LoadMacro("../plugins/ZGAngles.cc+");

mcfmFile = TFile(outpath+"out_mcfm_"+subdir+".root","RECREATE")
mcfmFile.mkdir("eff")
mcfmFile.cd()
h1 = HistManager(mcfmFile)

madFile = TFile(outpath+"out_mad_"+subdir+".root","RECREATE")
madFile.mkdir("eff")
madFile.cd()
h2 = HistManager(madFile)

ang = ZGAngles()

newDir = outpath 
print newDir
if not os.path.exists(newDir):
    os.makedirs(newDir)

def FillAllHists(files, h):
    # h is a hist manager instance here
    fChain = TChain("LHEF");
    for f in files:
        fChain.Add(f)
        #fChain.Print()

    dcount = 0
    for evt in fChain:

        g1 = TLorentzVector(0)
        g2 = TLorentzVector(0)
        g3 = TLorentzVector(0)
        l1 = TLorentzVector(0)
        l2 = TLorentzVector(0)
        gamma = TLorentzVector(0)
        diLep = TLorentzVector(0)

        hasZ = 0
        hasGlu3=0
        hasGamma=0

        for p in evt.Particle:
            px = p.Px
            py = p.Py
            pz = p.Pz
            E  = p.E
            M  = p.M
        
            if (p.PID == 22 and p.Status==1):
                gamma.SetPxPyPzE(px,py,pz,E)
                hasGamma=1
                
            if (p.PID == 13  or p.PID == 13):
            #if (p.PID == 13 or p.PID == 11):
                l1.SetPxPyPzE(px,py,pz,E)
            if (p.PID == -13 or p.PID ==-13):
            #if (p.PID == -13 or p.PID ==-11):
                l2.SetPxPyPzE(px,py,pz,E)

            if p.PID==23:
                hasZ=1

            if (p.PID==21 and p.Status==-1):
                g1.SetPxPyPzE(px,py,pz,E)
            if (p.PID==21 and p.Status==-1):
                g2.SetPxPyPzE(px,py,pz,E)
            if (p.PID==21 and p.Status==1):
                hasGlu3=1
                g3.SetPxPyPzE(px,py,pz,E)

    

        if l1.Pt()>l2.Pt():
            lPt1 = l1
            lPt2 = l2
        else:
            lPt1 = l2
            lPt2 = l1

    
        diLep = l1+l2
        if diLep.Pt()==0:
            continue

        dcount += 1
        #if dcount >10000:
        #    continue

        #if not hasGamma: continue
    
        tri = diLep + gamma
        
        gammaCM = TLorentzVector(gamma)
        diLepCM = TLorentzVector(diLep)
        b1  = TVector3(-tri.BoostVector())
        gammaCM.Boost(b1)
        diLepCM.Boost(b1)

        c1=Double(0)
        c2=c3=phi=1.1
        ang.SetAngles(l2,l1,gamma)
        c1 = ang.GetCos1()
        c2 = ang.GetCos2()
        c3 = ang.GetCosTheta()
        phi = ang.GetPhi()
        #print dcount, c1, c2, phi, c3
        #if dcount>20: break
        
        h.fill1DHist(c1,  "ang_co1",";gen cos_lp",100,-1,1, 1,"");
        h.fill1DHist(c2,  "ang_co2",";gen cos_lm,", 100,-1,1, 1,"");
        h.fill1DHist(c3,  "ang_co3",";gen cosTheta",100,-1,1, 1,"");
        h.fill1DHist(phi, "ang_phi",";gen phi lp",  100, -TMath.Pi(), TMath.Pi(), 1,"");

                   
        h.fill1DHist(l1.M(),    "l1_mass",  ";l+ mass",    200, -2,2, 1, "")
        h.fill1DHist(l2.M(),    "l2_mass",  ";l- mass",    200, -2,2, 1, "")

        #h.fill1DHist(g1.M(),    "g1_M",  ";g1 M",    200, -2,2, 1, "")
        #h.fill1DHist(g2.M(),    "g2_M",  ";g2 M",    200, -2,2, 1, "")
        
        h.fill1DHist(diLep.M(),     "gen_Mll_0",";gen_Mll",100,0,50, 1,"");
        #h.fill1DHist(diLep.M(),     "gen_Mll_0",";gen_Mll",100,0,50, 1,"eff");
        if lPt1.Pt()>23 and lPt2.Pt()>7 and fabs(lPt1.Eta())<2.4 and  fabs(lPt2.Eta())<2.4 \
               and gamma.Pt()>23 and fabs(gamma.Eta())<2.5:
            h.fill1DHist(diLep.M(),     "gen_Mll_1",";gen_Mll",100,0,50, 1,"");
            h.fill1DHist(diLep.M(),     "gen_Mll_2",";gen_Mll",100,0,50, 1,"");
            
            h.fill1DHist(diLep.M(),     "gen_Mll_3",";gen_Mll",100,0,50, 1,"");

        h.fill1DHist(gamma.M(),"gamma_mass",  ";gamma mass",    200, -2,2, 1, "")
        #h.fill1DHist(g1.Pt(),    "g1_pt",  ";g1 pt",    50, 0,100, 1, "")
        #h.fill1DHist(g2.Pt(),    "g2_pt",  ";g2 pt",    50, 0,100, 1, "")
        if hasGlu3:
            h.fill1DHist(g3.Pt(),   "g3_pt",  ";glu3 pt",   50, 0,100, 1, "")
            h.fill1DHist(g3.Eta(),  "g3_eta", ";glu3 eta",  50, -5,5, 1, "")
            h.fill1DHist(g3.Phi(),  "g3_phi", ";glu3 phi",  50, -4,4, 1, "")
        
        #h.fill1DHist(l1.Pt(),    "l1_pt",  ";l+ pt",    50, 0,100, 1, "")    
        #h.fill1DHist(l1.Eta(),   "l1_eta", ";l+ eta",   50, -3.5,3.5, 1, "")
        #h.fill1DHist(l1.Phi(),   "l1_phi", ";l+ phi",   50, -TMath.Pi(),TMath.Pi(), 1, "")
        #h.fill1DHist(l2.Pt(),    "l2_pt",  ";l- pt",    50, 0,100, 1, "")
        #h.fill1DHist(l2.Eta(),   "l2_eta", ";l- eta",   50, -3.5,3.5, 1, "")
        #h.fill1DHist(l2.Phi(),   "l2_phi", ";l- phi",   50, -TMath.Pi(),TMath.Pi(), 1, "")
        
        h.fill1DHist(diLep.M(),   "diLep_mass",     ";M(ll)", 200, 0,60,  1, "")
        h.fill1DHist(diLep.M(),   "diLep_mass_full",";M(ll)", 200, 0,130, 1, "")
        h.fill1DHist(diLep.M(),   "diLep_mass_low", ";M(ll)", 200, 0,1,   1, "")
        h.fill1DHist(tri.M(),     "h_mass",";M(ll#gamma)",    200, 80,180,1, "")
        h.fill1DHist(tri.M(),     "h_mass_zoom",";M(ll#gamma)",  200, 124,126,  1, "")
        h.fill1DHist(tri.M(),     "h_mass_zoom2",";M(ll#gamma)", 200, 124.9,125.1,  1, "")

        if not hasGlu3:
            h.fill1DHist(tri.Pt(),    "h_pt",";Pt of the Higgs",  200, 0,200,  1, "")
            h.fill1DHist(tri.Pz(),    "h_z", ";Pz of the Higgs",  200, -300,300,  1, "")
        else:
            h.fill1DHist(tri.Pt(),    "h_pt_2","Extra ISR glu;Pt of the Higgs",  200, 0,200,  1, "")
            h.fill1DHist(tri.Pz(),    "h_pz_2","Extra ISR glu;Pz of the Higgs",  200, -300,300,  1, "")
            h.fill1DHist((tri+g3).Pt(),    "h_pt_glu3","Extra ISR glu;Pt of the Higgs+gluon",  200, 0,200,  1, "")
            h.fill1DHist((tri+g3).Pz(),    "h_pz_glu3","Extra ISR glu;Pz of the Higgs+gluon",  200, -300,300,  1, "")
    

        h.fill1DHist(gammaCM.E(), "gamma_Ecom",";E_{#gamma} in CoM",  50, 0,200,  1, "")
        h.fill1DHist(diLepCM.E(), "diLep_Ecom",";Ecom(ll)", 50, 0,100,  1, "")


        h.fill1DHist(diLep.Pt(),    "diLep_pt",  ";diLep_pt",    50, 0,100, 1, "")
        h.fill1DHist(diLep.Eta(),   "diLep_eta", ";diLep_eta",   50, -3.5,3.5, 1, "")
        h.fill1DHist(diLep.Phi(),   "diLep_phi", ";diLep_phi",   50, -TMath.Pi(),TMath.Pi(), 1, "")
        h.fill1DHist(gamma.E(),  "gamma_E",  ";gamma_E",   50, 0,200, 1, "")
        h.fill1DHist(gamma.Pt(), "gamma_pt", ";gamma_pt",  50, 0,100, 1, "")
        h.fill1DHist(gamma.Eta(),"gamma_eta",";gamma_eta", 50, -3.5,3.5, 1, "")
        h.fill1DHist(gamma.Phi(),"gamma_phi",";gamma_phi", 50, -TMath.Pi(),TMath.Pi(), 1, "")
        
        h.fill1DHist(lPt1.Pt(),    "lPt1_pt",  ";Leading lepton pt",    50, 0,100, 1, "")
        h.fill1DHist(lPt1.Eta(),   "lPt1_eta", ";Leading lepton eta",   50, -3.5,3.5, 1, "")
        h.fill1DHist(lPt1.Phi(),   "lPt1_phi", ";Leading lepton phi",   50, -TMath.Pi(),TMath.Pi(), 1, "")
        h.fill1DHist(lPt2.Pt(),    "lPt2_pt",  ";Trailing lepton pt",    50, 0,100, 1, "")
        h.fill1DHist(lPt2.Eta(),   "lPt2_eta", ";Trailing lepton  eta",  50, -3.5,3.5, 1, "")
        h.fill1DHist(lPt2.Phi(),   "lPt2_phi", ";Trailing lepton  phi",  50, -TMath.Pi(),TMath.Pi(), 1, "")
        
        h.fill2DHist(lPt1.Pt(), lPt2.Pt(), "h2D_Pt1_vs_Pt2", ";Leading lepton pt; Trailing lepton pt",    50, 0,100, 50,0,100, 1, "")
        h.fill2DHist(l1.Pt(),   l2.Pt(),   "h2D_l1_vs_l2",   ";l+ pt; l- pt",    50, 0,100, 50,0,100, 1, "")
        h.fill2DHist(diLep.Pt(),  gamma.Pt(),"h2D_diLep_vs_gamma",     ";Pt of ll system; pt of gamma",   50, 0,100, 50,0,100, 1, "")
        h.fill2DHist(gammaCM.E(), gamma.Pt(),"h2D_gamma_Ecom_vs_Pt",   ";E_{#gamma} in CoM; Photon Pt",   50, 0,100, 50,0,100, 1, "")
        h.fill2DHist(gammaCM.E(), tri.M(),   "h2D_gamma_Ecom_vs_triM", ";E_{#gamma} in CoM; M(ll#gamma)", 50, 0,100, 50,0,200, 1, "")
        
        h.fill2DHist(gamma.Pt(), diLep.Eta()-gamma.Eta(),"h2D_gammaPt_vs_deltaEta", ";Pt_{#gamma}; #Delta#eta(ll, #gamma)",    50, 0,100, 50,-5,5, 1, "")
        
        
        h.fill1DHist(diLep.DeltaR(gamma),               "dR_diLep_gamma", ";dR(ll, #gamma)",         50, 0,10, 1, "")
        h.fill1DHist(fabs(diLep.Eta() - gamma.Eta()),   "dEta_diLep_gamma", ";|dEta(ll, #gamma)|",   50, 0,10, 1, "")
        h.fill1DHist((diLep.Vect()+gamma.Vect()).Pt(),  "diff_diLep_gamma_pt", ";four vector sum (diLep+#gamma).Pt()", 50, -20,20, 1, "")

        h.fill1DHist(TVector2.Phi_mpi_pi(diLep.Phi()-gamma.Phi()), "dPhi_diLep_gamma", ";dPhi(ll, #gamma)",            50, -10,10, 1, "")
        
        h.fill1DHist(l1.DeltaR(l2),     "dR_l1_l2",     ";dR(l+, l-)",      50, 0,5, 1, "")
        h.fill1DHist(diLep.DeltaR(l1),  "dR_diLep_l1",  ";dR(diLep, l+)",   50, 0,5, 1, "")
        h.fill1DHist(diLep.DeltaR(l2),  "dR_diLep_l2",  ";dR(diLep, l-)",   50, 0,5, 1, "")

        
    print "Total events = ", dcount
    
if __name__ == "__main__":
    gROOT.ProcessLine(".L ~/tdrstyle.C")
    setTDRStyle()
    TH1.SetDefaultSumw2(kTRUE)
    gStyle.SetOptStat(1)
    
    pathBase = outpath
    path = pathBase+subdir+"/"
    if not os.path.exists(path):
        os.makedirs(path)


    FillAllHists(files["mcfm"], h1)
    FillAllHists(files["mad"],  h2)

    mcfmFile.cd()
    mcfmFile.Write()
    madFile.cd()
    madFile.Write()
    print "Saved files: \n",mcfmFile.GetName(), "\n", madFile.GetName()


    blah = []
    u.drawAllInFile(madFile, "madgra",mcfmFile,"mcfm",None,"","", path, None,"norm")

    plot_types =["mcfm","mad"]
    if subdir not in plot_types:
        plot_types.append(subdir)
    ht.makeHTML("Plots from an lhe file",pathBase, plot_types, blah, "mad")

    mcfmFile.Close()
    madFile.Close()


    print "\n\t\t finita la comedia \n"

