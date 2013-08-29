#!/usr/bin/env python

import sys,os
import time
from ROOT import *
import plot as plot


subdir = sys.argv[1]

outpath = '/uscms_data/d2/andreypz/html/zgamma/lhe/'
files=[]
if subdir=="mad":
    files = ['/home/andreypz/workspace/MadGraph5/PROC_ANO_HEFT_JGF_v2_4/Events/run_01/unweighted_events.root']
else:
    files = ['/uscms_data/d2/andreypz/lhe_mcfm_hzg_dalitz_lord_fixed_unweighted.lhe.root']


print files
#gSystem.Load("/home/andreypz/workspace/MadGraph5/ExRootAnalysis/lib/libExRootAnalysis.so")
gSystem.Load("/uscms/home/andreypz/work/MadGraph5/ExRootAnalysis/lib/libExRootAnalysis.so")

gROOT.LoadMacro("../plugins/HistManager.cc+");

saveFile = TFile(outpath+"out_"+subdir+".root","RECREATE")
saveFile.cd()
h = HistManager(saveFile)

newDir = outpath 
print newDir
if not os.path.exists(newDir):
    os.makedirs(newDir)

fChain = TChain("LHEF");
for f in files:
    fChain.Add(f)
    #fChain.Print()

for evt in fChain:    

    g1 = TLorentzVector(0)
    g2 = TLorentzVector(0)
    g3 = TLorentzVector(0)
    l1 = TLorentzVector(0)
    l2 = TLorentzVector(0)
    gamma = TLorentzVector(0)
    diLep = TLorentzVector(0)

    dcount = 0
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

        if (p.PID == 13 or p.PID == 11):
            l1.SetPxPyPzE(px,py,pz,E)
        if (p.PID == -13 or p.PID ==-11):
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

        dcount += 1


    if l1.Pt()>l2.Pt():
        lPt1 = l1
        lPt2 = l2
    else:
        lPt1 = l2
        lPt2 = l1

    '''
    diLep = l1+l2

    #if not hasGamma: continue

    tri = diLep + gamma

    gammaCM = TLorentzVector(gamma)
    diLepCM = TLorentzVector(diLep)
    b1  = TVector3(-tri.BoostVector())
    gammaCM.Boost(b1)
    diLepCM.Boost(b1)


    h.fill1DHist(l1.M(),    "l1_mass",  ";l+ mass",    200, -2,2, 1, "")
    h.fill1DHist(l2.M(),    "l2_mass",  ";l- mass",    200, -2,2, 1, "")

    #h.fill1DHist(g1.M(),    "g1_M",  ";g1 M",    200, -2,2, 1, "")
    #h.fill1DHist(g2.M(),    "g2_M",  ";g2 M",    200, -2,2, 1, "")



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

    h.fill1DHist(diLep.M(),   "diLep_mass",";M(ll)", 200, 0,60,  1, "")
    h.fill1DHist(tri.M(),     "h_mass",";M(ll#gamma)",  200, 80,200,  1, "")
    if not hasGlu3:
        h.fill1DHist(tri.Pt(),    "h_pt","If NO extra glu;Pt of the Higgs",  200, 0,200,  1, "")
        h.fill1DHist(tri.Pz(),    "h_z", "If NO extra glu;Pz of the Higgs",  200, -300,300,  1, "")
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
    h.fill1DHist(gamma.E(),  "gamma_E", ";gamma_E",    50, 0,200, 1, "")
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
    h.fill2DHist(diLep.Pt(),    gamma.Pt(),"h2D_diLep_vs_gamma", ";Pt of ll system; pt of gamma",    50, 0,100, 50,0,100, 1, "")

    h.fill2DHist(gammaCM.E(), gamma.Pt(),"h2D_gamma_Ecom_vs_Pt", ";E_{#gamma} in CoM; Photon Pt",    50, 0,100, 50,0,100, 1, "")
    h.fill2DHist(gammaCM.E(), tri.M(),"h2D_gamma_Ecom_vs_triM", ";E_{#gamma} in CoM; M(ll#gamma)",    50, 0,100, 50,0,200, 1, "")

    h.fill2DHist(gamma.Pt(), diLep.Eta()-gamma.Eta(),"h2D_gammaPt_vs_deltaEta", ";Pt_{#gamma}; #Delta#eta(ll, #gamma)",    50, 0,100, 50,-5,5, 1, "")


    h.fill1DHist(diLep.DeltaR(gamma),               "dR_diLep_gamma", ";dR(ll, #gamma)",         50, 0,10, 1, "")
    h.fill1DHist(fabs(diLep.Eta() - gamma.Eta()),   "dEta_diLep_gamma", ";|dEta(ll, #gamma)|",   50, 0,10, 1, "")
    h.fill1DHist((diLep.Vect()+gamma.Vect()).Pt(),  "diff_diLep_gamma_pt", ";four vector sum (diLep+#gamma).Pt()",   50, -20,20, 1, "")
    h.fill1DHist(TVector2.Phi_mpi_pi(diLep.Phi()-gamma.Phi()), "dPhi_diLep_gamma", ";dPhi(ll, #gamma)",          50, -10,10, 1, "")

    h.fill1DHist(l1.DeltaR(l2), "dR_l1_l2", ";dR(l+, l-)", 50, 0,5, 1, "")
    h.fill1DHist(diLep.DeltaR(l1),  "dR_diLep_l1",  ";dR(diLep, l+)",   50, 0,5, 1, "")
    h.fill1DHist(diLep.DeltaR(l2),  "dR_diLep_l2",  ";dR(diLep, l-)",   50, 0,5, 1, "")


    '''

def draw(h, path):
    if h.InheritsFrom("TH2"):
        h.Draw("colz")
    else:
        h.Draw()

    c1.SaveAs(path+h.GetName()+".png")
    if h in [diLep_mass]:
        c1.SetLogy()
        c1.SaveAs(path+h.GetName()+"_log.png")
        c1.SetLogy(0)


gROOT.ProcessLine(".L ~/tdrstyle.C")
setTDRStyle()
TH1.SetDefaultSumw2(kTRUE)

pathBase = outpath
path = pathBase+subdir+"/"
if not os.path.exists(path):
    os.makedirs(path)

plot_types =["mcfm","mad"]


saveFile.cd()
saveFile.Write()
print "Saved file: ",saveFile.GetName()


blah = []
plot.drawAllInFile(saveFile, None, path, None)
plot.makeHTML("Plots from an lhe file",pathBase, plot_types, blah, "mad")

saveFile.Close()


print "\n\t\t finita la comedia \n"

