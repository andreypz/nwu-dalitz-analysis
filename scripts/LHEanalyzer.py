#!/usr/bin/env python

import sys,os
import time
from ROOT import *
gROOT.SetBatch()
sys.path.append("../zgamma")
import utils as u
import makeHTML as ht
from optparse import OptionParser

parser = OptionParser(usage="usage: %prog name [options --new]")
parser.add_option("-n","--new", dest="new",action="store_true", default=False, help="Make new hists. Otherwise re-use the existing ones")

(opt, args) = parser.parse_args()


subdir = sys.argv[1]

outpath = '/home/andreypz/workspace/html-lhe/lhe/'
#outpath = '/uscms_data/d2/andreypz/html/zgamma/lhe/'
files={}


#files["one"] = ['../../ano_zeromass_0.2to50.root']
#files["two"] = ['../../heft_zeromass_0.2to50.root']
#files["one"] = ['../../ano_dalitz.root']
#files["two"] = ['../../hc-ufo-dalitz.root']
files["one"] = ['../../hc-ufo-dalitz-ele-8TeV.root']
files["two"] = ['../../cmsgrid_final.lhe.root']
files["three"] = ['../../hc-ufo-dalitz-mmu-b-noT.root']
#files["three"] = ['../../hc-ufo-dalitz-ele-50K.root']
#files["two"] = ['../../ano_dalitz_mod.root']
#files["two"] = ['../../heft_dalitz.root']

#files["one"] = ['/tthome/andrey/LHE_dalitz_Luisa/dygamma.lhe.root']
#files["two"] = ['/tthome/andrey/LHE_dalitz_Luisa/dyjet.lhe.root']
fakeGammaFromJet = 0

MH = 125
LEPID1 = 13
LEPID2 = 11

print files
#gSystem.Load("/home/andreypz/workspace/MadGraph5/ExRootAnalysis/lib/libExRootAnalysis.so")
gSystem.Load("/home/andreypz/workspace/mg5amcnlo/ExRootAnalysis/libExRootAnalysis.so")
#gSystem.Load("/tthome/andrey/workspace/MG5_aMC_v2_1_1/ExRootAnalysis/lib/libExRootAnalysis.so")
#gSystem.Load("../plugins/HistManager_cc.so")
gROOT.LoadMacro("../plugins/HistManager.cc+");
gROOT.LoadMacro("../plugins/ZGAngles.cc+");

openOption = 'OPEN'
if opt.new:
  openOption = 'RECREATE'

oneFile  = TFile(outpath+"out_one_"+subdir+".root",openOption)
twoFile  = TFile(outpath+"out_two_"+subdir+".root",openOption)
testFile = TFile(outpath+"out_test_"+subdir+".root",openOption)

if opt.new:
  oneFile.cd()
  h1 = HistManager(oneFile)
  twoFile.cd()
  h2 = HistManager(twoFile)
  #testFile.mkdir("eff")
  testFile.cd()
  h3 = HistManager(testFile)

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
    # fChain.Print()
    print f

  dcount  = 0
  wpcount = 0
  wmcount = 0
  zcount  = 0

  for evt in fChain:
    g1 = TLorentzVector(0)
    g2 = TLorentzVector(0)
    g3 = TLorentzVector(0)
    l1 = TLorentzVector(0)
    l2 = TLorentzVector(0)
    j1 = TLorentzVector(0)
    j2 = TLorentzVector(0)
    gamma = TLorentzVector(0)
    diLep = TLorentzVector(0)
    trueHiggs = TLorentzVector(0)

    hasZ  = 0
    hasWp = 0
    hasWm = 0
    hasGlu3=0
    hasGamma=0
    hi = 0
    qq = 0
    for p in evt.Particle:
      px = p.Px
      py = p.Py
      pz = p.Pz
      E  = p.E
      M  = p.M

      #print p.PID, p.Status, px,py,pz,E,M

      if (p.PID == 25):
        trueHiggs.SetPxPyPzE(px,py,pz,E)
        hi+=1

      if abs(p.PID) < 6 and p.Status==1:
        qq+=1
        if qq==1:
          j1.SetPxPyPzE(px,py,pz,E)
        if qq==2:
          j2.SetPxPyPzE(px,py,pz,E)

      if (p.PID == 22 and p.Status==1):
        gamma.SetPxPyPzE(px,py,pz,E)
        hasGamma=1

      if fakeGammaFromJet:
        if ((p.PID == 21 and p.Status==1)  or (abs(p.PID) <= 5 and p.Status==1)):
          gamma.SetPxPyPzE(px,py,pz,E)
          hasGamma=1

      if (p.PID == LEPID1 or p.PID==LEPID2):
        l1.SetPxPyPzE(px,py,pz,E)
      if (p.PID == -LEPID1 or p.PID==-LEPID2):
        l2.SetPxPyPzE(px,py,pz,E)

      if p.PID==23:
        hasZ=1
      if p.PID==24:
        hasWp=1
      if p.PID==-24:
        hasWm=1

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
      # print "dLep.Pt = 0"
      continue

      # if qq!=2:
      #     print "Nope, there has to be two of them", qq
      #     sys.exit(0)

    if hasWm:
      wmcount += 1
    if hasWp:
      wpcount += 1
    if hasZ:
      zcount +=1

    # if not hasGamma: continue

    dcount += 1
    # if dcount > 5000: break


    tri = diLep + gamma

    if hi==0:
      # print dcount,"No higgs??? what's up with that??"
      h.fill1DHist(tri.M(),   "LHE_h_mass_noHiggs",     ";m_{ll#gamma}",  200, 80,180, 1, "")
      h.fill1DHist(tri.M(),   "LHE_h_mass_zoom_noHiggs",";m_{ll#gamma}",  200, 124,126,1, "")
      h.fill1DHist(diLep.M(), "LHE_h_mumu_noHiggs",     ";m_{ll}",        200, 80,180, 1, "")
      h.fill1DHist(gamma.Pt(),"LHE_gamma_pt_noHiggs",   ";p_{T}^{#gamma}",    200, 00,180, 1, "")
    else:
      h.fill1DHist(trueHiggs.M(),   "LHE_h_mass_trueHiggs",";m_{H}",   200, 124,126,1, "")
      h.fill1DHist(trueHiggs.Pt(),  "LHE_h_pt_", ";p_{T}^{H} (GeV)",   100, 0,100,  1, "")

    # exit(0)

    gammaCM = TLorentzVector(gamma)
    diLepCM = TLorentzVector(diLep)
    b1  = TVector3(-tri.BoostVector())
    gammaCM.Boost(b1)
    diLepCM.Boost(b1)

    #gamma.Print()

    c1=Double(0)
    c2=c3=phi=1.1
    ang.SetAngles(l2,l1,gamma)
    c1 = ang.GetCos1()
    c2 = ang.GetCos2()
    c3 = ang.GetCosTheta()
    phi = ang.GetPhi()

    #print dcount, c1, c2, phi, c3

    h.fill1DHist(c1,  "LHE_ang_co1",";gen cos_lp",  100,-1,1, 1,"");
    h.fill1DHist(c2,  "LHE_ang_co2",";gen cos_lm",  100,-1,1, 1,"");
    h.fill1DHist(c3,  "LHE_ang_co3",";gen cosTheta",100,-1,1, 1,"");
    h.fill1DHist(phi, "LHE_ang_phi",";gen phi lp",  100, -TMath.Pi(), TMath.Pi(), 1,"");


    h.fill1DHist(l1.M(),    "l1_mass",  ";l+ mass",    200, -2,2, 1, "")
    h.fill1DHist(l2.M(),    "l2_mass",  ";l- mass",    200, -2,2, 1, "")

    # h.fill1DHist(g1.M(),    "g1_M",  ";g1 M",    200, -2,2, 1, "")
    # h.fill1DHist(g2.M(),    "g2_M",  ";g2 M",    200, -2,2, 1, "")

    '''
    Mll_max = 20
    h.fill1DHist(diLep.M(),     "gen_Mll_0",";gen_Mll, GeV",100,0,Mll_max, 1,"");
    if gamma.Pt()>25 and fabs(gamma.Eta())<2.5:
            h.fill1DHist(diLep.M(),     "gen_Mll_1",";gen_Mll, GeV",100,0,Mll_max, 1,"")

            if lPt1.Pt()>23 and lPt2.Pt()>4 and fabs(lPt1.Eta())<2.4 and  fabs(lPt2.Eta())<2.4:
                h.fill1DHist(diLep.M(),     "gen_Mll_2",";gen_Mll, GeV",100,0,Mll_max, 1,"");

                for a in xrange(3,11):
                    if (diLep.Pt()/tri.M() > 0.10+0.05*(a-3) and gamma.Pt()/tri.M() > 0.10 + +0.05*(a-3)):
                        h.fill1DHist(diLep.M(),     "gen_Mll_"+str(a),";gen_Mll, GeV",100,0,Mll_max, 1,"");

                for a in xrange(11,17):
                    if (diLep.Pt()> 25+5*(a-11)):
                        h.fill1DHist(diLep.M(),     "gen_Mll_"+str(a),";gen_Mll, GeV",100,0,Mll_max, 1,"");

                for a in xrange(17,23):
                    if (gamma.Pt()> 25+5*(a-17)):
                        h.fill1DHist(diLep.M(),     "gen_Mll_"+str(a),";gen_Mll, GeV",100,0,Mll_max, 1,"");

                for a in xrange(23,31):
                    if (diLep.Pt()/tri.M() > 0.10+0.05*(a-23)):
                        h.fill1DHist(diLep.M(),     "gen_Mll_"+str(a),";gen_Mll, GeV",100,0,Mll_max, 1,"");

                for a in xrange(31,39):
                    if (gamma.Pt()/tri.M() > 0.10+0.05*(a-31)):
                        h.fill1DHist(diLep.M(),     "gen_Mll_"+str(a),";gen_Mll, GeV",100,0,Mll_max, 1,"");

    '''

    h.fill1DHist(gamma.M(),"gamma_mass",  ";gamma mass",    200, -2,2, 1, "")
        #h.fill1DHist(g1.Pt(),    "g1_pt",  ";g1 pt",    50, 0,100, 1, "")
        #h.fill1DHist(g2.Pt(),    "g2_pt",  ";g2 pt",    50, 0,100, 1, "")

        #h.fill1DHist(l1.Pt(),    "l1_pt",  ";l+ pt",    50, 0,100, 1, "")
        #h.fill1DHist(l1.Eta(),   "l1_eta", ";l+ eta",   50, -3.5,3.5, 1, "")
        #h.fill1DHist(l1.Phi(),   "l1_phi", ";l+ phi",   50, -TMath.Pi(),TMath.Pi(), 1, "")
        #h.fill1DHist(l2.Pt(),    "l2_pt",  ";l- pt",    50, 0,100, 1, "")
        #h.fill1DHist(l2.Eta(),   "l2_eta", ";l- eta",   50, -3.5,3.5, 1, "")
        #h.fill1DHist(l2.Phi(),   "l2_phi", ";l- phi",   50, -TMath.Pi(),TMath.Pi(), 1, "")

    h.fill1DHist(diLep.M(),   "LHE_diLep_mass",     ";m_{ll} (GeV)", 100, 0,60,  1, "")
    h.fill1DHist(diLep.M(),   "LHE_diLep_mass_bins",";m_{ll} (GeV)", 5000,0,55,  1, "")
    h.fill1DHist(diLep.M(),   "LHE_diLep_mass_full",";m_{ll} (GeV)", 100, 0,130, 1, "")
    h.fill1DHist(diLep.M(),   "LHE_diLep_mass_low", ";m_{ll} (GeV)", 100, 0,1,   1, "")
    
    if LEPID1==11 or LEPID2==11:
      h.fill1DHist(diLep.M(),   "LHE_diLep_xxx", ";m_{ll} (GeV)", 200, 0,0.05,   1, "")

    h.fill1DHist(tri.M(),     "LHE_h_mass",    ";m_{ll#gamma} (GeV)",100, 80,180,1, "")
    h.fill1DHist(tri.M(),     "LHE_h_mass_low",";m_{ll#gamma} (GeV)",100,  0,50, 1, "")
    if MH!=0:
      h.fill1DHist(tri.M(),     "LHE_h_mass_zoom", ";m_{ll#gamma} (GeV)", 200, MH-1,MH+1,  1, "")
      h.fill1DHist(tri.M(),     "LHE_h_mass_zoom2",";m_{ll#gamma} (GeV)", 200, MH-0.1,MH+0.1,  1, "")

        ## VBF Plots:
    if qq==2:
      h.fill1DHist((j1+j2).M(),             "LHE_diJet_mass",";M(j1,j2), GeV", 200, 0,600,  1, "")
      h.fill1DHist(fabs(j1.Eta()+j2.Eta()), "LHE_diJet_dEta",";|#Delta#eta(j_{1},j_{2})|",  200,   0,6,  1, "")

      '''
        if not hasGlu3:
            h.fill1DHist(tri.Pt(),    "LHE_h_pt",";Pt of the Higgs",  200, 0,200,  1, "")
            h.fill1DHist(tri.Pz(),    "LHE_h_z", ";Pz of the Higgs",  200, -300,300,  1, "")
        else:
        h.fill1DHist(g3.Pt(),   "LHE_g3_pt",  ";glu3 pt",   50, 0,100, 1, "")
        h.fill1DHist(g3.Eta(),  "LHE_g3_eta", ";glu3 eta",  50, -5,5, 1, "")
        h.fill1DHist(g3.Phi(),  "LHE_g3_phi", ";glu3 phi",  50, -4,4, 1, "")

            h.fill1DHist(tri.Pt(),    "h_pt_2","Extra ISR glu;Pt of the Higgs",  200, 0,200,  1, "")
            h.fill1DHist(tri.Pz(),    "h_pz_2","Extra ISR glu;Pz of the Higgs",  200, -300,300,  1, "")
            h.fill1DHist((tri+g3).Pt(),    "h_pt_glu3","Extra ISR glu;Pt of the Higgs+gluon",  200, 0,200,  1, "")
            h.fill1DHist((tri+g3).Pz(),    "h_pz_glu3","Extra ISR glu;Pz of the Higgs+gluon",  200, -300,300,  1, "")
      '''

    h.fill1DHist(gammaCM.E(), "gamma_Ecom",";E_{#gamma} in CoM, GeV",  50, 0,200,  1, "")
    h.fill1DHist(diLepCM.E(), "diLep_Ecom",";E^{com}_{ll}, GeV", 50, 0,100,  1, "")


    h.fill1DHist(diLep.Pt(), "LHE_diLep_pt",  ";diLep_pt, GeV",    50, 0,100, 1, "")
    h.fill1DHist(diLep.Eta(),"LHE_diLep_eta", ";diLep_eta",     50, -3.5,3.5, 1, "")
    h.fill1DHist(diLep.Phi(),"LHE_diLep_phi", ";diLep_phi",   50, -TMath.Pi(),TMath.Pi(), 1, "")
    h.fill1DHist(gamma.E(),  "LHE_gamma_E",   ";E_{#gamma}, GeV", 50, 0,200, 1, "")
    h.fill1DHist(gamma.Pt(), "LHE_gamma_pt",  ";p_{T}^{#gamma}, GeV",  50, 0,100, 1, "")
    h.fill1DHist(gamma.Eta(),"LHE_gamma_eta", ";#eta_{#gamma}", 50, -3.5,3.5, 1, "")
    h.fill1DHist(gamma.Phi(),"LHE_gamma_phi", ";#phi_{#gamma}", 50, -TMath.Pi(),TMath.Pi(), 1, "")

    h.fill1DHist(lPt1.Pt(),    "LHE_lPt1_pt",  ";Leading lepton p_{T} (GeV)", 50, 0,100, 1, "")
    h.fill1DHist(lPt1.Eta(),   "LHE_lPt1_eta", ";Leading lepton eta",   50, -3.5,3.5, 1, "")
    h.fill1DHist(lPt1.Phi(),   "LHE_lPt1_phi", ";Leading lepton phi",   50, -TMath.Pi(),TMath.Pi(), 1, "")
    h.fill1DHist(lPt2.Pt(),    "LHE_lPt2_pt",  ";Trailing lepton p_{T} (GeV)", 50, 0,100, 1, "")
    h.fill1DHist(lPt2.Eta(),   "LHE_lPt2_eta", ";Trailing lepton  eta",  50, -3.5,3.5, 1, "")
    h.fill1DHist(lPt2.Phi(),   "LHE_lPt2_phi", ";Trailing lepton  phi",  50, -TMath.Pi(),TMath.Pi(), 1, "")

    h.fill2DHist(lPt1.Pt(), lPt2.Pt(), "h2D_Pt1_vs_Pt2", ";Leading lepton pt; Trailing lepton pt",    100, 0,80, 100,0,50, 1, "")
    h.fill2DHist(l1.Pt(),   l2.Pt(),   "h2D_l1_vs_l2",   ";l+ pt; l- pt",    50, 0,100, 50,0,100, 1, "")
    h.fill2DHist(diLep.Pt(),  gamma.Pt(),"h2D_diLep_vs_gamma",     ";Pt of ll system; pt of gamma",   50, 0,100, 50,0,100, 1, "")
    h.fill2DHist(gammaCM.E(), gamma.Pt(),"h2D_gamma_Ecom_vs_Pt",   ";E_{#gamma} in CoM; Photon Pt",   50, 0,100, 50,0,100, 1, "")
    h.fill2DHist(gammaCM.E(), tri.M(),   "h2D_gamma_Ecom_vs_triM", ";E_{#gamma} in CoM; M(ll#gamma)", 50, 0,100, 50,0,200, 1, "")

    h.fill2DHist(gamma.Pt(), diLep.Eta()-gamma.Eta(),"h2D_gammaPt_vs_deltaEta",
                 ";p_{T}_{#gamma}; #Delta#eta(ll, #gamma)",    50, 0,100, 50,-5,5, 1, "")


    h.fill1DHist(diLep.DeltaR(gamma),              "LHE_dR_diLep_gamma", ";dR(ll, #gamma)",         50, 0,10, 1, "")
    h.fill1DHist(fabs(diLep.Eta() - gamma.Eta()),  "LHE_dEta_diLep_gamma", ";|dEta(ll, #gamma)|",   50, 0,10, 1, "")
    h.fill1DHist((diLep.Vect()+gamma.Vect()).Pt(), "LHE_diff_diLep_gamma_pt", ";(diLep+#gamma).Pt()", 50, -20,20, 1, "")
    h.fill1DHist(TVector2.Phi_mpi_pi(diLep.Phi()-gamma.Phi()), "LHE_dPhi_diLep_gamma", ";#Delta#phi(ll, #gamma)", 50, -10,10, 1, "")

    h.fill1DHist(lPt1.DeltaR(gamma),"LHE_dR_lPt1_gamma",";#Delta R(l_{1}, #gamma)",  50, 0,6, 1, "")
    h.fill1DHist(lPt2.DeltaR(gamma),"LHE_dR_lPt2_gamma",";#Delta R(l_{2}, #gamma)",  50, 0,6, 1, "")
    h.fill1DHist(l1.DeltaR(l2),     "LHE_dR_l1_l2",     ";#Delta R(l+, l-)",      50, 0,5, 1, "")
    h.fill1DHist(l1.DeltaR(l2),     "LHE_dR_l1_l2_low", ";#Delta R(l+, l-)",      50, 0,1, 1, "")
    h.fill1DHist(l1.DeltaR(l2),     "LHE_dR_l1_l2_vlow",";#Delta R(l+, l-)",    50, 0,0.3, 1, "")
    h.fill1DHist(diLep.DeltaR(l1),  "LHE_dR_diLep_l1",  ";#Delta R(diLep, l+)",   50, 0,5, 1, "")
    h.fill1DHist(diLep.DeltaR(l2),  "LHE_dR_diLep_l2",  ";#Delta R(diLep, l-)",   50, 0,5, 1, "")


  print "Total events = ", dcount
  print "Has Z:", zcount, "  haswp", wpcount, "  hasWm:", wmcount

if __name__ == "__main__":
  gROOT.ProcessLine(".L ../tdrstyle.C")
  setTDRStyle()
  TH1.SetDefaultSumw2(kTRUE)
    #gStyle.SetOptStat(1)

  pathBase = outpath
  path = pathBase+subdir
  if not os.path.exists(path):
    os.makedirs(path)


  if opt.new:
    FillAllHists(files["one"],  h1)
    FillAllHists(files["two"],  h2)
    FillAllHists(files["three"],  h3)
    oneFile.cd()
    oneFile.Write()
    twoFile.cd()
    twoFile.Write()
    testFile.cd()
    testFile.Write()
    print "Saved files: \n",testFile.GetName(), "\n", oneFile.GetName(), "\n", twoFile.GetName()


  blah = ['LHE for mH = 125 GeV',
          'Comparing ANO and HEFT']
  #blah = ['LHE files comparison for DY+gamma sample',
  #        'Luisa\'s vs mine']


  sigZip = zip(['Ele', 'Grid Mu ','Mu 13TeV'],
                [oneFile, twoFile, testFile])
  #sigZip = zip(['ANO-0','ANO-50'],
  #              [oneFile, twoFile])

  u.drawAllInFile(None, None, None, sigZip, 'HEFT', '', path, None,"norm", isLog=False)


  #myzip = zip(['5<mll<20'],[testFile])
  #u.drawAllInFile(oneFile, "0<mll<2", myzip, twoFile, '2<mll<5', '', path, None,"norm", isLog=True)
  #u.drawAllInFile(oneFile, "ggH to eeg", '', twoFile, 'VH to eeg', '', path, None,"norm", isLog=True)
  #u.drawAllInFile(oneFile, "h #rightarrow ee#gamma", '', twoFile, 'h #rightarrow #mu#mu#gamma', '', path, None,"norm", isLog=False)
  #u.drawAllInFile(oneFile, "3 GeV", [('15 GeV',twoFile)],testFile, '35 GeV', '', path, None,"norm", isLog=True)
  #u.drawAllInFile(oneFile, "DYG", [('DYJ',twoFile)],testFile, 'Sig', '', path, None,"norm", isLog=True)


  u.createDir(pathBase+"/csBR")
  c1 = TCanvas("c1","c1", 600,500)
  h1 = twoFile.Get("LHE_diLep_mass_bins")
  h2 = twoFile.Get("LHE_diLep_mass_bins")
  #h.Print('all')
  h1.Draw()
  h2.Draw()
  c1.SaveAs(pathBase+"/csBR/dilepmass.png")

  '''
  xsbr = 1
  # xsbr = 0.754
  tot  = h.Integral()
  bin0 = h.FindBin(0)
  bin1 = h.FindBin(0.)
  bin2 = h.FindBin(0.5)
  bin3 = h.FindBin(1)
  bin4 = h.FindBin(2)
  bin5 = h.FindBin(4)
  bin6 = h.FindBin(9)
  bin7 = h.FindBin(20)
  bin8 = h.FindBin(50)
  part1 = h.Integral(bin0, bin1)
  part2 = h.Integral(bin1, bin2)
  part3 = h.Integral(bin2, bin3)
  part4 = h.Integral(bin3, bin4)
  part5 = h.Integral(bin4, bin5)
  part6 = h.Integral(bin5, bin6)
  part7 = h.Integral(bin6, bin7)
  part8 = h.Integral(bin7, bin8)
  print 'bin 1 content = ', h.GetBinContent(1),
  print "Integrals: tot = %.3f; part1 = %.3f"%(tot, part1)
  print "In cs: 0-50:  %.3f"%(xsbr)
  print "In cs: x: %.4f"%(xsbr*part1/tot)
  print "In cs: x: %.4f"%(xsbr*part2/tot), bin2
  print "In cs: x: %.4f"%(xsbr*part3/tot), bin3
  print "In cs: x: %.4f"%(xsbr*part4/tot), bin4
  print "In cs: x: %.4f"%(xsbr*part5/tot), bin5
  print "In cs: x: %.4f"%(xsbr*part6/tot), bin6
  print "In cs: x: %.4f"%(xsbr*part7/tot), bin7
  print "In cs: x: %.4f"%(xsbr*part8/tot), bin8
  '''

  '''
  u.createDir(path+"/eff")
  c1.cd()
  hh = []
  ra = []
  for n in xrange(0,39):
      print n
      hh.append(oneFile.Get("gen_Mll_"+str(n)))
      ra.append(hh[-1].Clone())
      ra[-1].Divide(hh[0])

      ra[3].Draw("hist")
    ra[3].SetMinimum(0)
    ra[3].SetMaximum(1)
    ra[3].SetTitle(";Mll gen at LHE; acc")
    ra[3].SetLineColor(29)
    leg = TLegend(0.70,0.65,0.90,0.98);
    leg.AddEntry(ra[3],"pT/Mllg > 0.10", "l")
    for a in xrange(3,11):
      ra[a].Draw('same hist')
      ra[a].SetLineColor(30+a-3)
      leg.AddEntry(ra[a],"pT/mllg > %.2f"%(0.15+0.05*(a-3)), "l")

    leg.SetTextSize(0.03)
    leg.SetFillColor(kWhite)
    leg.Draw()
    c1.SaveAs(path+"/eff/acceptance_Mll_LHE_ptMlg.png")



    ra[11].Draw("hist")
    ra[11].SetMinimum(0)
    ra[11].SetMaximum(1)
    ra[11].SetTitle(";Mll gen at LHE; acc")
    ra[11].SetLineColor(29)
    leg = TLegend(0.70,0.65,0.90,0.98);
    leg.AddEntry(ra[11],"pT(ll) > 25 GeV", "l")
    for a in xrange(12,17):
      ra[a].Draw('same hist')
      ra[a].SetLineColor(30+a-12)
      leg.AddEntry(ra[a],"pT(ll)> %.0f GeV"%(30.+5*(a-12)), "l")

    leg.SetTextSize(0.03)
    leg.SetFillColor(kWhite)
    leg.Draw()
    c1.SaveAs(path+"/eff/acceptance_Mll_LHE_pTll.png")


    ra[17].Draw("hist")
    ra[17].SetMinimum(0)
    ra[17].SetMaximum(1)
    ra[17].SetTitle(";Mll gen at LHE; acc")
    ra[17].SetLineColor(29)
    leg = TLegend(0.70,0.65,0.90,0.98);
    leg.AddEntry(ra[17],"pT(#gamma) > 25 GeV", "l")
    for a in xrange(18,23):
      ra[a].Draw('same hist')
      ra[a].SetLineColor(30+a-18)
      leg.AddEntry(ra[a],"pT(#gamma)> %.0f GeV"%(30.+5*(a-18)), "l")

    leg.SetTextSize(0.03)
    leg.SetFillColor(kWhite)
    leg.Draw()
    c1.SaveAs(path+"/eff/acceptance_Mll_LHE_pTgamma.png")



    ra[23].Draw("hist")
    ra[23].SetMinimum(0)
    ra[23].SetMaximum(1)
    ra[23].SetTitle(";Mll gen at LHE; acc")
    ra[23].SetLineColor(29)
    leg = TLegend(0.70,0.65,0.90,0.98);
    leg.AddEntry(ra[23],"pT(ll)/Mllg > 0.10", "l")
    for a in xrange(24,31):
      ra[a].Draw('same hist')
      ra[a].SetLineColor(30+a-24)
      leg.AddEntry(ra[a],"pT(ll)/Mllg > %.2f"%(0.15 + 0.05*(a-24)), "l")

    leg.SetTextSize(0.03)
    leg.SetFillColor(kWhite)
    leg.Draw()
    c1.SaveAs(path+"/eff/acceptance_Mll_LHE_pTllMllg.png")

    ra[31].Draw("hist")
    ra[31].SetMinimum(0)
    ra[31].SetMaximum(1)
    ra[31].SetTitle(";Mll gen at LHE; acc")
    ra[31].SetLineColor(29)
    leg = TLegend(0.70,0.65,0.90,0.98);
    leg.AddEntry(ra[31],"pT(#gamma)/Mllg > 0.10", "l")
    for a in xrange(32,39):
      ra[a].Draw('same hist')
      ra[a].SetLineColor(30+a-32)
      leg.AddEntry(ra[a],"pT(#gamma)/Mllg > %.2f"%(0.15 + 0.05*(a-31)), "l")

    leg.SetTextSize(0.03)
    leg.SetFillColor(kWhite)
    leg.Draw()
    c1.SaveAs(path+"/eff/acceptance_Mll_LHE_pTgammaMllg.png")
  '''

  plot_types =[]
  list = os.listdir(pathBase)
  for d in list:
    if os.path.isdir(pathBase+"/"+d):
      plot_types.append(d)

  ht.makeHTML("Plots from an lhe file",pathBase, plot_types, blah, subdir)

  testFile.Close()
  oneFile.Close()
  twoFile.Close()

  print "\n\t\t finita la comedia \n"

