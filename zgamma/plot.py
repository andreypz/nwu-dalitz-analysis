#! /usr/bin/env python
from optparse import OptionParser
import sys,os,datetime
from array import *
from ROOT import *
import makeHTML as ht
import utils as u
gROOT.SetBatch()

parser = OptionParser(usage="usage: %prog ver [options -c, -e, -p, -m]")
parser.add_option("-c","--cut",   dest="cut", type="int", default=4,    help="Plots after a certain cut")
parser.add_option("-m","--merge", dest="merge",action="store_true", default=False, help="Do merging?")

parser.add_option("-p", "--period",dest="period", default="2012",  help="Year period; 2011 or 2012")
parser.add_option("--bkg", dest="bkg",  action="store_true", default=False, help="Make plots from bkg sample")
parser.add_option("--mcfm",dest="mcfm", action="store_true", default=False, help="Use MCFM  as a signal")

(options, args) = parser.parse_args()

lumi = u.getLumi(options.period)
cs   = u.getCS(options.mcfm)

#sel = ["electron"]
#sel = ["mugamma"]
#sel = ["mugamma","electron"]
sel = []


def effPlots(f1, path):
    print "Now making efficiency plots"
    #f1.cd(dir)

    for var in ["Mll","dR"]:
    
        h0 = f1.Get("eff/gen_"+var+"_0")
        h1 = f1.Get("eff/gen_"+var+"_acc_gamma")
        h2 = f1.Get("eff/gen_"+var+"_acc_lept")
        h3 = f1.Get("eff/gen_"+var+"_reco_gamma_iso")
        h4 = f1.Get("eff/gen_"+var+"_one_ele_reco_ID")
        h5 = f1.Get("eff/gen_"+var+"_two_ele_reco")
        h6 = f1.Get("eff/gen_"+var+"_two_ele_reco_ID")
        h7 = f1.Get("eff/gen_"+var+"_one_ele_reco")
        

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
        r5 = h5.Clone()
        r5.Divide(h2)
        r6 = h6.Clone()
        r6.Divide(h2)
        r7 = h7.Clone()
        r7.Divide(h2)

        
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
        leg.AddEntry(r1,"photon pt>38, eta<2.5", "l")
        leg.AddEntry(r2,"photon pt>38 and  pt(e1,e2) > (23,7)", "l")
        leg.SetTextSize(0.04)
        leg.SetFillColor(kWhite)
        leg.Draw()
        c1.SaveAs(path+"acceptance_"+var+".png")
        
        r3.Draw("hist")
        r3.SetMinimum(0)
        r3.SetMaximum(1)
        r3.SetTitle(xname+" gen; reco eff")    
        r4.Draw("hist same")
        r5.Draw("hist same")
        r6.Draw("hist same")
        r7.Draw("hist same")
    
        r3.SetLineColor(kBlack)
        r4.SetLineColor(kOrange+1)
        r5.SetLineColor(kGreen+1)
        r6.SetLineColor(kRed+1)
        leg = TLegend(0.35,0.15,0.98,0.30);
        leg.AddEntry(r3,"Reco photon pt>38, tight ID", "l")
        leg.AddEntry(r5,"Photon + 2 electrons NO ID pt(e1,e2) > (23,7)", "l")
        leg.AddEntry(r6,"Photon + 2 electrons ELE ID pt(e1,e2) > (23,7)", "l")
        leg.AddEntry(r7,"Photon + ONE electron pt(e) > 30, No ID", "l")
        leg.AddEntry(r4,"Photon + ONE electron pt(e) > 30, ID", "l")
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
    #hPath  = "/uscms_data/d2/andreypz/zgamma/"+ver

    u.createDir(pathBase)



    if doMerge:
        os.system("rm "+hPath+"/m_*.root") #removing the old merged files
    yields_data = {}
    yields_bkg  = {}
    yields_sig  = {}

    tri_hists = {}
    dataFile  = {}
    bkgFile   = {}


    '''
    #bkgFile = TFile("hhhh_dy.root", "OPEN")
    bkgFile1 = TFile("~/nobackup/v13_DY_electron_2012.root","OPEN")
    sigFile  = TFile("~/nobackup/v13_dalitz_electron_2012.root","OPEN")
    #sigFile = TFile("hhhh_sig.root", "OPEN")
    #effPlots(sigFile, pathBase+"/eff/")


    u.createDir(pathBase+"/eff")
    prof = sigFile.Get("eff/gen_Mll_vs_dR")
    prof.Draw("")
    c1.SaveAs(pathBase+"/eff/gen_Mll_vs_dR.png")

    h1 = sigFile.Get("Electrons/egamma_reco").Clone()
    h2 = bkgFile1.Get("Electrons/egamma_reco").Clone()


    #Nev1 = sigFile.Get("Counts/evt_byCut").GetBinContent(3)
    #Nev2 = bkgFile.Get("Counts/evt_byCut").GetBinContent(3)
    Nev1 = sigFile.Get("Muons/size_mu_cut1").Integral()
    Nev2 = bkgFile1.Get("Muons/size_mu_cut1").Integral()
    print Nev1, Nev2

    h1.Scale(1./Nev1)    
    h2.Scale(1./Nev2)
    h1.SetLineColor(kRed+1)
    h1.Draw("hist")
    h2.Draw("hist same")
    h1.GetXaxis().SetBinLabel(1,"no eele, no gamma");
    h1.GetXaxis().SetBinLabel(2,"electron");
    h1.GetXaxis().SetBinLabel(3,"gamma");
    h1.GetXaxis().SetBinLabel(4,"ele and gamma");
    h1.GetXaxis().SetBinLabel(5,"ele or gamma");
    c1.SaveAs(pathBase+"/eff/reco_egamma.png")

    '''
    


    for thissel in sel:
        if doMerge:
            os.system("hadd "+hPath+"/m_Data_"    +thissel+"_"+period+".root "+hPath+"/"+thissel+"_"+period+"/hhhh_*Run20*.root")
            #os.system("hadd "+hPath+"/m_DY_"    +thissel+"_"+period+".root "+hPath+"/"+thissel+"_"+period+"/hhhh_DYjets*.root")
            #os.system("hadd "+hPath+"/m_ZG_"    +thissel+"_"+period+".root "
            #         +hPath+"/"+thissel+"_"+period+"/hhhh_ZG_*.root ")

        subdir = thissel
        path = pathBase+"/"+subdir
        if doBkg:
            path = pathBase+"/bkg_"+subdir
            u.createDir(path)

        
            sigFileMCFM = TFile(hPath+"/"+thissel+"_"+period+"/hhhh_h-dalitz_1.root", "OPEN")
            sigFileMAD  = TFile(hPath+"/"+thissel+"_"+period+"/hhhh_mad_1.root", "OPEN")
        if options.mcfm: sigFile = sigFileMCFM
        else:            sigFile = sigFileMCFM
        
        dataFile[thissel] = TFile(hPath+"/m_Data_"+thissel+"_"+period+".root","OPEN")
        #bkgFile[thissel]  = TFile(hPath+"/m_DY_"+thissel+"_"+period+".root","OPEN")

        yields_data[thissel] = u.getYields(dataFile[thissel])
        yields_sig[thissel]  = u.getYields(sigFile,True)
        #yields_bkg[thissel]  = u.getYields(bkgFile[thissel])

        #if int(cut) >2:
        #tri_hists[thissel]   = dataFile[thissel].Get("tri_mass_cut"+cut).Clone()
        

        u.drawAllInFile(dataFile[thissel], "data", None, "", sigFile,"10xSignal",  "",path, cut, "lumi")
        u.drawAllInFile(dataFile[thissel], "data", None, "", sigFile,"10xSignal","EB",pathBase+"/EB", cut, "lumi")
        u.drawAllInFile(dataFile[thissel], "data", None, "", sigFile,"10xSignal","EE",pathBase+"/EE", cut, "lumi")
        if thissel =="mugamma":
            u.drawAllInFile(dataFile[thissel], "data",None, "",sigFile,"signal",  "Muons", pathBase+"/Muons", None,"norm")

            #u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",  "Muons", pathBase+"/Muons/", None,"norm")
            #u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",  "Photon",pathBase+"/Photon/", None,"norm")
        elif thissel =="electron":
            u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",
                            "Photon",     pathBase+"/Photon/",     None,"norm")
            u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",
                            "DalitzEle",  pathBase+"/DalitzEle/",  None,"norm")
            u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",
                            "DalitzEleEB",pathBase+"/DalitzEleEB/",None,"norm")
            u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",
                            "DalitzEleEE",pathBase+"/DalitzEleEE/",None,"norm")
            u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",
                            "NewEle-1",   pathBase+"/NewEle-1/",   None,"norm")
            u.drawAllInFile(dataFile[thissel], "data",bkgFile[thissel], "bkg",sigFile,"signal",
                            "NewEle-2",   pathBase+"/NewEle-2/",   None,"norm")

            
        #dataFile.Close()
    #print yields_data


    #u.drawAllInFile(bkgFile[thissel], "DY electrons", sigFile, "Dalitz 2el",  "NewEle-1", pathBase+"/NewEle-1/", None,"norm", isLog=1, )

    sigFileMCFM = TFile(hPath+"/mugamma_"+period+"/hhhh_h-dalitz_1.root", "OPEN")
    sigFileMAD  = TFile(hPath+"/mugamma_"+period+"/hhhh_mad_1.root", "OPEN")

    u.drawAllInFile(sigFileMCFM, "MCFM",None, "",sigFileMAD,"Madgraph",  "Muons", pathBase+"/MC", None,"norm")


    '''
    c1 = TCanvas("c4","small canvas",600,600);
    c1.cd()
    Nev = sigFile.Get("Counts/evt_byCut").GetBinContent(2)
    cro = cs["h"]
    scale = float(1000*lumi*cro)/Nev

    m1 = 121.
    m2 = 129.
    for r in ["","EB","EE"]:
        if int(cut)!=8: break
        pre = r
        if r!="":
            pre=r+"/"
        hda = dataFile[thissel].Get(pre+"tri_mass_"+r+"_cut"+cut).Clone()
        hsi = sigFile.Get(pre+"tri_mass_"+r+"_cut"+cut).Clone()
        yda = hda.Integral()
        ysi = hsi.Integral()/10
        print r, yda,ysi, ysi*scale
        print "significance=", ysi/sqrt(ysi+yda)

        bin1 = hda.FindBin(m1)
        bin2 = hda.FindBin(m2)
        yda_zoom = hda.Integral(bin1,bin2)
        ysi_zoom = hsi.Integral(bin1,bin2)/10
        
        print "in [",m1,", ",m2,"]", "  bins:", bin1, bin2
        print r, "data=", yda_zoom, "sig=", ysi_zoom
        print "  ==> zoomed sign = ", ysi_zoom/sqrt(ysi_zoom + yda_zoom)
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
    table_sig  = u.yieldsTable(yields_sig, sel)
    table_data = u.yieldsTable(yields_data, sel)


    u.makeTable(table_data,"data", "html")
    u.makeTable(table_sig, "sig",  "html")
    u.makeTable(table_data,"data", "twiki")
    u.makeTable(table_sig, "sig",  "twiki")

    os.system("cat yields_data.html yields_sig.html > yields.html")

    comments = ["These plots are made for ...",
                "Blah"]
    
    ht.makeHTML("h &rarr; dalitz decay plots",pathBase, plot_types, comments, "mugamma")

    print "\n\t\t finita la comedia \n"
