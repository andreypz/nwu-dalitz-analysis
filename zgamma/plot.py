#! /usr/bin/env python
from optparse import OptionParser
import sys,os,datetime
from array import *
from ROOT import *
import utils as u
import makeHTML as ht
gROOT.SetBatch()

parser = OptionParser(usage="usage: %prog ver [options -c, -e, -p, -m]")
parser.add_option("-c","--cut",   dest="cut", type="int", default=4,    help="Plots after a certain cut")
parser.add_option("--mass", dest="mass", type="int", default=125,    help="Signal sample (mass)")
parser.add_option("-m","--merge", dest="merge",action="store_true", default=False, help="Do merging?")

parser.add_option("-p", "--period",dest="period", default="2012",  help="Year period; 2011 or 2012")
parser.add_option("--bkg",  dest="bkg",  action="store_true", default=False, help="Make plots from bkg sample")
parser.add_option("--mcfm", dest="mcfm", action="store_true", default=False, help="Use MCFM  as a signal")
parser.add_option("--noeos",dest="noeos",action="store_true", default=False, help="Don't use EOS. pick up the files from nobackup area")

parser.add_option("--mumu",   dest="mumu",    action="store_true", default=False, help="MuMuGamma selection with Double-Mu trigger")
parser.add_option("--mugamma",dest="mugamma", action="store_true", default=False, help="MuMuGamma selection with Mu-Pho trigger")
parser.add_option("--elgamma",dest="elgamma", action="store_true", default=False, help="EEGamma selection")

(options, args) = parser.parse_args()

mass = options.mass

sel = []
if options.mugamma:
    sel.append("mugamma")
if options.mumu:
    sel.append("mumu")
if options.elgamma:
    sel.append("elgamma")


def effPlots2(f1, path):
    c1.cd()
    hh = []
    ra = []
    for n in xrange(0,31):
        print n
        hh.append(f1.Get("eff/gen_Mll_"+str(n)))
        ra.append(hh[-1].Clone())
        ra[-1].Divide(hh[0])

    ra[1].Draw("hist")
    ra[1].SetMinimum(0)
    ra[1].SetMaximum(1)
    ra[1].SetTitle(";Mll; acc * eff")
    bcol = 29
    ra[1].SetLineColor(bcol)
    ra[1].SetLineWidth(2)
    leg = TLegend(0.70,0.65,0.95,0.98);
    leg.AddEntry(ra[1],"pT(ll) > 25 GeV", "l")
    for a in xrange(2,7):
        ra[a].Draw('same hist')
        ra[a].SetLineColor(bcol+2*(a-2))
        ra[a].SetLineWidth(2)
        leg.AddEntry(ra[a],"pT(ll) > %.0f GeV"%(25+5*(a-1)), "l")

    leg.SetTextSize(0.03)
    leg.SetFillColor(kWhite)
    leg.Draw()
    c1.SaveAs(path+"/acc_eff_Mll_ptll.png")

    st = 7
    ed = 12
    ra[st].Draw("hist")
    ra[st].SetMinimum(0)
    ra[st].SetMaximum(1)
    ra[st].SetTitle(";Mll; acc")
    bcol = 33
    ra[st].SetLineColor(bcol)
    ra[st].SetLineWidth(2)
    leg = TLegend(0.70,0.65,0.95,0.98);
    leg.AddEntry(ra[st],"pT(#gamma) > 25 GeV", "l")
    for a in xrange(st+1,ed+1):
        ra[a].Draw('same hist')
        ra[a].SetLineColor(bcol+2*(a-st))
        ra[a].SetLineWidth(2)
        leg.AddEntry(ra[a],"pT(#gamma) > %.0f GeV"%(25+5*(a-st)), "l")

    leg.SetTextSize(0.03)
    leg.SetFillColor(kWhite)
    leg.Draw()
    c1.SaveAs(path+"/acc_eff_Mll_ptgamma.png")

    st = 13
    ed = 18
    ra[st].Draw("hist")
    ra[st].SetMinimum(0)
    ra[st].SetMaximum(1)
    ra[st].SetTitle(";Mll; acc")
    bcol = 38
    ra[st].SetLineColor(bcol)
    ra[st].SetLineWidth(2)
    leg = TLegend(0.70,0.65,0.95,0.98);
    leg.AddEntry(ra[st],"pT(#gamma)/Mllg > 0.20", "l")
    for a in xrange(st+1,ed+1):
        ra[a].Draw('same hist')
        ra[a].SetLineColor(bcol+2*(a-st))
        ra[a].SetLineWidth(2)
        leg.AddEntry(ra[a],"pT(#gamma)/Mllg > %.2f"%(0.20+0.05*(a-st)), "l")

    leg.SetTextSize(0.03)
    leg.SetFillColor(kWhite)
    leg.Draw()
    c1.SaveAs(path+"/acc_eff_Mll_ptgammaMllg.png")

    st = 19
    ed = 24
    ra[st].Draw("hist")
    ra[st].SetMinimum(0)
    ra[st].SetMaximum(1)
    ra[st].SetTitle(";Mll; acc")
    bcol = 41
    ra[st].SetLineColor(bcol)
    ra[st].SetLineWidth(2)
    leg = TLegend(0.70,0.65,0.95,0.98);
    leg.AddEntry(ra[st],"pT(ll)/Mllg > 0.20", "l")
    for a in xrange(st+1,ed+1):
        ra[a].Draw('same hist')
        ra[a].SetLineColor(bcol+2*(a-st))
        ra[a].SetLineWidth(2)
        leg.AddEntry(ra[a],"pT(ll)/Mllg > %.2f"%(0.20+0.05*(a-st)), "l")

    leg.SetTextSize(0.03)
    leg.SetFillColor(kWhite)
    leg.Draw()
    c1.SaveAs(path+"/acc_eff_Mll_ptllMllg.png")


    st = 25
    ed = 30
    ra[st].Draw("hist")
    ra[st].SetMinimum(0)
    ra[st].SetMaximum(1)
    ra[st].SetTitle(";Mll; acc")
    bcol = 43
    ra[st].SetLineColor(bcol)
    ra[st].SetLineWidth(2)
    leg = TLegend(0.70,0.65,0.95,0.98);
    leg.AddEntry(ra[st],"pT(both)/Mllg > 0.20", "l")
    for a in xrange(st+1,ed+1):
        ra[a].Draw('same hist')
        ra[a].SetLineColor(bcol+2*(a-st))
        ra[a].SetLineWidth(2)
        leg.AddEntry(ra[a],"pT(both)/Mllg > %.2f"%(0.20+0.05*(a-st)), "l")

    leg.SetTextSize(0.03)
    leg.SetFillColor(kWhite)
    leg.Draw()
    c1.SaveAs(path+"/acc_eff_Mll_ptbothMllg.png")


def effPlots(f1, path):
    print "Now making efficiency plots"
    #f1.cd(dir)

    for var in ["Mll","dR"]:
    
        h0 = f1.Get("eff/gen_"+var+"_0")
        h1 = f1.Get("eff/gen_"+var+"_acc_gamma")
        h2 = f1.Get("eff/gen_"+var+"_acc_lept")

        h3 = f1.Get("eff/gen_"+var+"_reco_gamma_iso")
        h4 = f1.Get("eff/gen_"+var+"_two_lep_reco")
        

        h0.Draw("hist")
        #h1.Draw("hist same")
        #h2.Draw("hist same")
        #h3.Draw("hist same")
        c1.SaveAs(path+var+".png")
        
        #for acceptance
        r1 = h1.Clone()
        r1.Divide(h0)
        r2 = h2.Clone()
        r2.Divide(h0)
        #for reco eff
        r3 = h3.Clone()
        r3.Divide(h2)
        r4 = h4.Clone()
        r4.Divide(h2)
        #r5 = h5.Clone()
        #r5.Divide(h2)
        #r6 = h6.Clone()
        #r6.Divide(h2)
        #r7 = h7.Clone()
        #r7.Divide(h2)

        
        if var=="Mll":
            xname = ";M(l1,l2)"
        if var=="dR":
            xname = ";dR(l1,l2)"

        r1.Draw("hist")
        r2.Draw("hist same")
        r1.SetMinimum(0)
        r1.SetMaximum(1)
        r1.SetTitle(xname+" gen; acc")
        r1.SetLineColor(kRed+1)
        r2.SetLineColor(kGreen+1)
        leg = TLegend(0.20,0.2,0.90,0.30);
        leg.AddEntry(r1,"photon pt>25, eta<2.5", "l")
        leg.AddEntry(r2,"photon pt>25 and p_{T}(l1)>23, p_{T}(l2)>4", "l")
        leg.SetTextSize(0.04)
        leg.SetFillColor(kWhite)
        leg.Draw()
        c1.SaveAs(path+"acceptance_"+var+".png")
        
        r3.Draw("hist")
        r3.SetMinimum(0)
        r3.SetMaximum(1)
        r3.SetTitle(xname+" gen; reco eff")    
        r4.Draw("hist same")
        #r5.Draw("hist same")
        #r6.Draw("hist same")
        #r7.Draw("hist same")
    
        r3.SetLineColor(kBlack)
        r4.SetLineColor(kOrange+1)
        #r5.SetLineColor(kGreen+1)
        #r6.SetLineColor(kRed+1)
        leg = TLegend(0.35,0.15,0.98,0.30);
        leg.AddEntry(r3,"Reco photon pt>25, Medium WP ID ", "l")
        leg.AddEntry(r4,"Photon and two reco-muons (OUR ID/ISO)", "l")
        #leg.AddEntry(r5,"Photon + 2 electrons NO ID pt(e1,e2) > (23,7)", "l")
        #leg.AddEntry(r6,"Photon + 2 electrons ELE ID pt(e1,e2) > (23,7)", "l")
        #leg.AddEntry(r7,"Photon + ONE electron pt(e) > 30, No ID", "l")
        leg.SetTextSize(0.025)
        leg.SetFillColor(kWhite)
        leg.Draw()
        c1.SaveAs(path+"eff_"+var+".png")
        

if __name__ == "__main__":
    timer = TStopwatch()
    timer.Start()
    
    if len(args) < 1:
        parser.print_usage()
        exit(1)
        
    ver    = sys.argv[1]
    #subdir = sys.argv[3]        
    cut=str(options.cut)
    doMerge = options.merge
    period  = options.period
    doBkg   = options.bkg
    
    gROOT.ProcessLine(".L ~/tdrstyle.C")
    setTDRStyle()
    TH1.SetDefaultSumw2(kTRUE)
    
    
    pathBase = "/uscms_data/d2/andreypz/html/zgamma/dalitz/"+ver+"_cut"+cut
    hPath    = "/eos/uscms/store/user/andreypz/batch_output/zgamma/8TeV/"+ver
    if options.noeos:
        hPath  = "/uscms_data/d2/andreypz/zgamma/"+ver

    u.createDir(pathBase)



    if doMerge:
        os.system("rm "+hPath+"/m_*.root") #removing the old merged files
    yields_data = {}
    yields_bkg  = {}
    yields_sig  = {}
    yields_vbf  = {}
    yields_vh   = {}

    tri_hists = {}
    dataFile  = {}
    bkgFile   = {}


    for thissel in sel:
        if thissel=='none': continue
        u.setSelection(thissel)

        
        if doMerge:
            if thissel =="elgamma":
                os.system("hadd "+hPath+"/m_Data_" +thissel+"_"+period+".root "+hPath+"/"+thissel+"_"+period+"/hhhh_*Run2012D*.root")
            else:
                os.system("hadd "+hPath+"/m_Data_" +thissel+"_"+period+".root "+hPath+"/"+thissel+"_"+period+"/hhhh_*Run20*.root")

            if doBkg:
                os.system("hadd "+hPath+"/m_DY_"   +thissel+"_"+period+".root "+hPath+"/"+thissel+"_"+period+"/hhhh_DYjets*.root")
            #os.system("hadd "+hPath+"/m_ZG_"    +thissel+"_"+period+".root "
            #         +hPath+"/"+thissel+"_"+period+"/hhhh_ZG_*.root ")

        subdir = thissel
        path = pathBase+"/"+subdir
        if doBkg:
            path = pathBase+"/bkg_"+subdir
            u.createDir(path)
        path = pathBase+"/"+subdir
                
        
        #sigFileMCFM = TFile(hPath+"/"+thissel+"_"+period+"/hhhh_ggHMCFM_1.root", "OPEN")
        #sigFileMCFM = TFile(hPath+"/"+thissel+"_"+period+"/hhhh_ggH-mcfm_1.root", "OPEN")
        sigFileMAD  = TFile(hPath+"/"+thissel+"_"+period+"/hhhh_ggH-mad"+str(mass)+"_1.root", "OPEN")
        sigFileVBF  = TFile(hPath+"/"+thissel+"_"+period+"/hhhh_vbf-mad"+str(mass)+"_1.root", "OPEN")
        sigFileVH   = TFile(hPath+"/"+thissel+"_"+period+"/hhhh_vh-mad"+str(mass)+"_1.root", "OPEN")

        if options.mcfm: sigFile = sigFileMCFM
        else:            sigFile = sigFileMAD
        
        dataFile[thissel] = TFile(hPath+"/m_Data_"+thissel+"_"+period+".root","OPEN")
        if doBkg:
            bkgFile[thissel]  = TFile(hPath+"/"+thissel+"_"+period+"/hhhh_DYjets0_1.root", "OPEN")
            #bkgFile[thissel]  = TFile(hPath+"/m_DY_"+thissel+"_"+period+".root","OPEN")


        yields_data[thissel] = u.getYields(dataFile[thissel])
        yields_sig[thissel]  = u.getYields(sigFile, 'ggH-125',True) 
        yields_vbf[thissel]  = u.getYields(sigFileVBF, 'vbfH-125', True) 
        yields_vh[thissel]   = u.getYields(sigFileVH,  'vH-125', True) 
        #yields_sig[thissel]  = u.getYields(sigFile, False)
        #yields_bkg[thissel]  = u.getYields(bkgFile[thissel])

        if doBkg:
            yields_bkg[thissel]  = u.getYields(bkgFile[thissel])

        #if int(cut) >2:
        #tri_hists[thissel]   = dataFile[thissel].Get("tri_mass_cut"+cut).Clone()
        

        u.drawAllInFile(dataFile[thissel], "Data", None, "", sigFile,"Signal",  "",path, cut, "norm")
        #u.drawAllInFile(dataFile[thissel], "data", None, "", sigFile,"signal",  "",path, cut, "norm", doRatio=1)
        #u.drawAllInFile(dataFile[thissel], "Data", None, "", sigFile,"50xSignal",  "",path, cut, "lumi")
        u.drawAllInFile(dataFile[thissel], "data", None, "", sigFile,"50xSignal","EB",pathBase+"/EB", cut, "lumi")
        u.drawAllInFile(dataFile[thissel], "data", None, "", sigFile,"50xSignal","EE",pathBase+"/EE", cut, "lumi")

        if thissel =="mugamma":
            u.drawAllInFile(dataFile[thissel], "data",None, "",sigFile,"signal",  "Muons", pathBase+"/Muons", None,"norm")

            #u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",  "Muons", pathBase+"/Muons/", None,"norm")
            #u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",  "Photon",pathBase+"/Photon/", None,"norm")
        elif thissel =="elgamma":
            if not doBkg:
                bkgFile[thissel]=None


            u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",
                            "Photon",     pathBase+"/Photon",     None,"norm")
            
            u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",
                            "DalitzEle-Before",  pathBase+"/DalitzEle-Before",  None,"norm")
            u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",
                            "DalitzEle-AfterAll",  pathBase+"/DalitzEle-AfterAll",  None,"norm")
            u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",
                            "DalitzEle-Before_tracks",  pathBase+"/DalitzEle-Before_tracks",  None,"norm")
            u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",
                            "DalitzEle-AfterAll_tracks",  pathBase+"/DalitzEle-AfterAll_tracks",  None,"norm")
            
            
            u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",
                            "DalitzEle",  pathBase+"/DalitzEle",  None,"norm")

            u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",
                            "DalitzEleEB",pathBase+"/DalitzEleEB",None,"norm")
            u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",
                            "DalitzEleEE",pathBase+"/DalitzEleEE/",None,"norm")

            u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",
                            "DalitzEleEB_tracks",pathBase+"/DalitzEleEB_tracks",None,"norm")
            u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",
                            "DalitzEleEE_tracks",pathBase+"/DalitzEleEE_tracks",None,"norm")

            
        #dataFile.Close()
    #print yields_data


    #u.drawAllInFile(sigFileMAD, "signal", None, "", None,"","GEN-RECO",pathBase+"/GEN-RECO", None, "lumi")
    
    #u.drawAllInFile(bkgFile[thissel], "DY electrons", sigFile, "Dalitz 2el",  "NewEle-1", pathBase+"/NewEle-1/", None,"norm", isLog=1, )

    
    #sigFileMAD  = TFile(hPath+"/mugamma_"+period+"/hhhh_ggH-mad125_1.root", "OPEN")
    if sel[0]=='none':
        sigFileMAD  = TFile("hhhh_ggH-mad125_1.root", "OPEN")

        u.createDir(pathBase+"/eff")
        effPlots(sigFileMAD, pathBase+"/eff/")
        effPlots2(sigFileMAD, pathBase+"/eff/")



    htmp = sigFileMAD.Get("lPt2_pt__cut10")
    int1 = htmp.Integral()
    int2 = htmp.Integral(0,10)
    print '\n Fraction of event with trailing lepton pt<20: \n', float(int2)/int1

    htmp = sigFileMAD.Get("ll_deltaR__cut10")
    int1 = htmp.Integral()
    int2 = htmp.Integral(0,10)
    print '\n Fraction of event with dR(l1,l2)<0.4: \n', float(int2)/int1

    '''
    c1 = TCanvas("c4","small canvas",600,600);
    c1.cd()
    Nev = sigFileMAD.Get("Counts/evt_byCut").GetBinContent(2)
    cro = u.getCS("ggH-125")
    lumi = u.getLumi("2012")

    scale = float(lumi*cro)/Nev

    hc7 = sigFileMAD.Get("tri_mass80__cut7").Integral(35,41)*scale
    hc8 = sigFileMAD.Get("tri_mass80__cut8").Integral(35,41)*scale
    hc9 = sigFileMAD.Get("tri_mass80__cut9").Integral(35,41)*scale
    print "Yields in bins that supposed to correspond to [122,128] window:\n",hc7, hc8, hc9 

    hc7 = dataFile[thissel].Get("tri_mass80__cut7").Integral(36,40)
    hc8 = dataFile[thissel].Get("tri_mass80__cut8").Integral(36,40)
    hc9 = dataFile[thissel].Get("tri_mass80__cut9").Integral(36,40)
    print "Yields in bins that supposed to correspond to [122,128] window:\n",hc7, hc8, hc9 


    m1 = 110.
    m2 = 170.
    
    for r in ["0"]:
        etaCut = TCut("")
        if r=="EB":
            etaCut = "fabs(ph_eta)<1"
        elif r=="EE":
            etaCut = "fabs(ph_eta)>1"
            
        cut = TCut("m_llg>"+str(m1)+"&&m_llg<"+str(m2))
        cut += etaCut
        treeda = dataFile[thissel].Get("fitTree/fitTree")
        treeda.Draw("m_llg>>hda", cut)
        
        treesi = sigFile.Get("fitTree/fitTree")
        treesi.Draw("m_llg>>hsi", cut)
        
        yda = hda.Integral()
        ysi = hsi.Integral()
        ysi_sc = ysi*scale
        print 'yields inside ', m1,m2, ' r=',r, 'data=', yda, 'signal=',ysi_sc
        print "significance=", ysi_sc/sqrt(ysi_sc+yda)
    
    '''
    '''
    h2da = dataFile[thissel].Get("h2D_dalitzPlot_rotation__cut"+cut).ProjectionX("hda_prx")
    h2si = sigFile.Get("h2D_dalitzPlot_rotation__cut"+cut).ProjectionX("hsi_prx")

    #u.handleOverflowBins(h2da)
    #u.handleOverflowBins(h2si)
        
    h2da.Draw("hist")
    h2si.Draw("same hist")
    h2si.SetLineColor(kRed+1)
    h2si.Scale(10*scale)
    h2da.SetTitle(";projectionX;Events")
    c1.SaveAs("h2da_dalitz.png")
    '''
    
    '''
    nx = h2da.GetNbinsX()
    ny = h2da.GetNbinsY()
    h2da_rot = TH2D("h2da_rot","", nx, )
    h2si_rot = TH2D("h2si_rot","")
    for a in nx:
        for b in ny:
            fda = h2da.GetBinContent(a,b)
            fsi = h2si.GetBinContent(a,b)

            h2da_rot.SetBinContent(a,b, fda)
            h2si_rot.SetBinContent(a,b,fda)
    '''
                        
    plot_types =[]
    list = os.listdir(pathBase)
    for d in list:
        if os.path.isdir(pathBase+"/"+d):
            plot_types.append(d)

    #plot_types
    print yields_sig

    table_all  = u.yieldsTable([yields_data,yields_sig, yields_vbf, yields_vh], sel)

    if doBkg:
        table_bkg = u.yieldsTable(yields_bkg, sel)        

    u.makeTable(table_all,"all", "html")
    u.makeTable(table_all,"all", "twiki")

    if doBkg:
        u.makeTable(table_bkg, "bkg",  "twiki")

    os.system("cat yields_all.html > yields.html")

    comments = ["These plots are made for ...",
                "Blah"]
    
    ht.makeHTML("h &rarr; dalitz decay plots",pathBase, plot_types, comments, "mugamma")

    print "\n\t\t finita la comedia \n"
