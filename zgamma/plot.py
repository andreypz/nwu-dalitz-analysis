#! /usr/bin/env python
from optparse import OptionParser
import sys,os,datetime
from array import *
from ROOT import *
import makeHTML as ht
gROOT.SetBatch()

parser = OptionParser(usage="usage: %prog ver [options -c, -e, -p, -m]")
parser.add_option("-c","--cut",    dest="cut", type="int", default=4,    help="Plots after a certain cut")
parser.add_option("-m", "--merge", dest="merge",action="store_true", default=False, help="Do merging?")
parser.add_option("-e", "--ele",   dest="ele",  action="store_true", default=False, help="Use electron selection")
parser.add_option("-p", "--period",dest="period", default="2012",  help="Year period; 2011 or 2012")
parser.add_option("--bkg", dest="bkg",  action="store_true", default=False, help="Make plots from bkg sample")

(options, args) = parser.parse_args()

import ConfigParser as cp
conf = cp.ConfigParser()
conf.read('config.cfg')
lumi2012 = float(conf.get("lumi","lumi2012A")) + float(conf.get("lumi","lumi2012B"))+\
           float(conf.get("lumi","lumi2012C")) + float(conf.get("lumi","lumi2012D"))
lumi = lumi2012
#sel = ["electron"]
sel = ["mugamma"]
#sel = ["mugamma","muon","electron"]

cuts = []
for key, cut in sorted(conf.items("cuts2")):
    cuts.append(cut)
    #print key, cut

cs = {}
for sample, c in conf.items("cs"):
    print sample,c
    cs[sample] = float(c)
cs["h"]=2*cs["h"]

def handleOverflowBinsScaleAndColors(hist, sample, lumi):
    if hist == None:
        return
    nBins   = hist.GetNbinsX()
    lastBin = hist.GetBinContent(nBins)
    ovflBin = hist.GetBinContent(nBins+1);
    lastBinErr = hist.GetBinError(nBins);
    ovflBinErr = hist.GetBinError(nBins+1);
    firstBin    = hist.GetBinContent(1);
    undflBin    = hist.GetBinContent(0);
    firstBinErr = hist.GetBinError(1);
    undflBinErr = hist.GetBinError(0);
    hist.SetBinContent(nBins, lastBin+ovflBin);
    hist.SetBinError(nBins, sqrt(pow(lastBinErr,2) + pow(ovflBinErr,2)) );
    hist.SetBinContent(1, firstBin+undflBin);
    hist.SetBinError(1, sqrt(pow(firstBinErr,2) + pow(undflBinErr,2)) );
    hist.SetBinContent(0,0);
    hist.SetBinContent(nBins+1,0);

    #if sample!="data":
     #   if lumi!=0: #don't rescale, for the cases we don't want it
      #      hist.Scale(float(lumi*cs(sample))/Nev(sample))
            #hist.SetLineColor(xc[sample][0])
            #hist.SetFillStyle(xc[sample][1])
            #hist.SetFillColor(xc[sample][2])


def draw(h1, path):
    h = h1
    if h.InheritsFrom("TH2"):
        h.Draw("lego")
    else:
        h.Draw("hist")

    c1.SaveAs(path+h.GetName()+".png")

def blindIt(h):
    hbin =  h.FindBin(125)
    for b in xrange(hbin-10,hbin+10):
        h.SetBinContent(b,0)
    
def set_palette(name="palette", ncontours=999):
    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""

    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
    elif name=="signal":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.90, 0.60, 0.40, 0.20]
        green = [0.00, 0.00, 0.00, 0.00, 0.00]
        blue  = [0.00, 0.00, 0.00, 0.00, 0.00]
    # elif name == "whatever":
        # (define more palettes)
    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)

    npoints = len(s)
    TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    gStyle.SetNumberContours(ncontours)
    
def drawAllInFile(f1, name1, f2, name2, dir,path, N, howToScale="none"):
    f1.cd(dir)
    dirList = gDirectory.GetListOfKeys()
    #dirList.Print()

    scale = 1

    if f2!=None and howToScale=="lumi": # only assume signal MC for now
        Nev = f2.Get("Counts/evt_byCut_raw").GetBinContent(1)
        cro = cs["h"]
        scale = float(1000*lumi*cro)/Nev
        print Nev, lumi, cro, scale

    for k1 in dirList:
        if k1.GetName() in ["Counts"]: continue
        if N!=None:
            if not("cut"+N) in k1.GetName(): continue
        h1 = k1.ReadObj()

        h2 = TH1F()
        
        if f2!=None:
            #f2.Print()
            if dir!="":
                h2 = f2.Get(dir+"/"+k1.GetName()) #assumes that the histograms in signal file have the same names  
            else:
                print k1.GetName()
                h2 = f2.Get(k1.GetName()).Clone()
                h2.Scale(float(scale))
                
        print "drawing", h1.GetName()

        if h1.InheritsFrom("TH2"):
            #set_palette("signal")
            h1.Draw("col")
            #set_palette("gray")
            #h1.Draw("col same")
            if "h2D_tri_vs_diLep_mass" in h1.GetName():
                c1.SetLogy()
        else:
            h1.Draw("hist")
            if "tri_mass" in k1.GetName() and name1 not in ["madgra","mcfm"]:
                blindIt(h1)
            if "h_mass" in k1.GetName():
                print "\n *** H-mass RMS:", h1.GetRMS(),h2.GetRMS()
                print "\n *** H-mass Mean:", h1.GetMean(),h2.GetMean()
                
            leg = TLegend(0.70,0.8,0.90,0.90);
            leg.AddEntry(h1,name1, "l")
            leg.SetTextSize(0.04)

            if f2!=None:
                h2.Draw("same hist")
                h2.SetLineColor(kRed+1)
                norm1 = h1.Integral()
                norm2 = h2.Integral()
                if howToScale =="norm":
                    h1.Scale(1./norm1)
                    h2.Scale(1./norm2)
                #print "Integrals = ", norm1, norm2
                if name1 not in ["madgra","mcfm"] and howToScale!="norm":
                    h2.Scale(100)
                leg.AddEntry(h2,name2, "l")
                if h1.GetName() in ["reco_gen_l1_deltaR", "reco_gen_l2_deltaR", "gen_ll_deltaR",
                                "reco_gen_gamma_deltaR", "ll_deltaR"]:
                    c1.SetLogy()
                
            leg.SetFillColor(kWhite)
            leg.Draw()

            c1.SetLogy()

        c1.SaveAs(path+h1.GetName()+".png")
        c1.SetLogy(0)
        

def createDir(dir):
    if not os.path.exists(dir):
        try:
            os.makedirs(dir)
        except OSError:
            if os.path.isdir(dir):
                pass
            else:
                raise
def yieldsTable(yi):
    t = []
    l1 = ["Cut/trigger"]
    l1.extend([a for a in sel])
    print l1
    t.append(l1)
    for line in xrange(len(cuts)):
        l=[]
        l.append(cuts[line])
        for thissel in sel:
            l.append(yi[thissel][line])
            #print l
        t.append(l)
        
    return t
                                                        
def makeTable(table, opt="tex"):
    print "Making sure that the list is alright"
    n_row = len(table)
    n_col = len(table[0])
    for l in table:
        if len(l)!=n_col:
            print "No good, the number of columns is messed up"
            
            
    myTable = ''
    if opt=="tex":
        beginTable = '\\begin{tabular}{|'+n_col*"l|"+'} \n \\hline \n'
        endTable   = '\\end{tabular} \n'

        beginLine  = ''
        endLine    = '\\\\\\hline \n'
        separator  = ' & '


    if opt=="html":
        beginTable = '<table border = "10"    cellpadding="5">'
        endTable = '</table>'

        beginLine  = '\n<tr>\n<td>'
        endLine    = '</td>\n</tr>'
        separator  = '</td><td>'

    if opt=="twiki":
        beginTable = ''
        endTable   = ''

        beginLine  = '| '
        endLine    = ' |\n'
        separator  = ' |  '


    myTable +=beginTable
    for l in range(n_row):
        #print l, table[l]
        myTable+=beginLine
        for c in range(n_col):
            val = table[l][c]
            if not isinstance(val,str):
                myTable+="%.2f" % (table[l][c])
            else:
                myTable+=val
            if c!=n_col-1:
                myTable+=separator
        myTable+=endLine

    myTable +=endTable
    
    ifile = open("yields.html","w")
    ifile.write(myTable)
    ifile.close()
    print myTable

def getYields(f, doLumiScale=False):
    ev = f.Get("Counts/evt_byCut")

    y = []
    scale=1
    if doLumiScale: # only assume signal MC for now
        Nev = ev.GetBinContent(1)
        cro = cs["h"]
        scale = float(1000*lumi*cro)/Nev

    for a in xrange(len(cuts)):
        y.append(scale*ev.GetBinContent(a+1)) # well, that's how the histogram is set up
    return y

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
    period = options.period
    doBkg = options.bkg
    
    gROOT.ProcessLine(".L ~/tdrstyle.C")
    setTDRStyle()
    TH1.SetDefaultSumw2(kTRUE)
    
    
    pathBase = "/uscms_data/d2/andreypz/html/zgamma/dalitz/"+ver+"_cut"+cut
    hPath = "/eos/uscms/store/user/andreypz/batch_output/zgamma/8TeV/"+ver
    hPath = "/uscms_data/d2/andreypz/zgamma/"+ver


    if doMerge:
        os.system("rm "+hPath+"/m_*.root") #removing the old merged files
    yields_data = {}
    yields_bkg = {}
    yields_sig  = {}

    tri_hists = {}
    dataFile = {}
    for thissel in sel:
    #for thissel in ["muon","mugamma","single-mu"]:
        if doMerge:
            os.system("hadd "+hPath+"/m_Data_"    +thissel+"_"+period+".root "+hPath+"/"+thissel+"_"+period+"/hhhh_*Run20*.root")

        subdir = thissel
        path = pathBase+"/"+subdir+"/"
        if doBkg:
            path = pathBase+"/bkg_"+subdir+"/"
        createDir(path)
        createDir(pathBase+"/Muons")

        sigFile = TFile(hPath+"/"+thissel+"_"+period+"/hhhh_h-dalitz_1.root", "OPEN")
        #bkgFile = TFile(hPath+"/"+thissel+"_"+period+"/hhhh_DY-mg5_1.root",   "OPEN")
        dataFile[thissel] = TFile(hPath+"/m_Data_"+thissel+"_"+period+".root","OPEN")

        yields_data[thissel] = getYields(dataFile[thissel])
        yields_sig[thissel]  = getYields(sigFile,True)
        #yields_bkg[thissel]  = getYields(bkgFile)

        if int(cut) >2:
            tri_hists[thissel]   = dataFile[thissel].Get("tri_mass_cut"+cut).Clone()
        

        drawAllInFile(dataFile[thissel], "data", sigFile,"100x h #rightarrow ll#gamma",  "",path, cut, "lumi")
        if thissel =="mugamma":
            drawAllInFile(dataFile[thissel], "data",sigFile,"signal",  "Muons",pathBase+"/Muons/", None,"norm")
        
        #dataFile.Close()
    #print yields_data

    plot_types =[]
    list = os.listdir(pathBase)
    for d in list:
        if os.path.isdir(pathBase+"/"+d):
            plot_types.append(d)

    '''
    if int(cut) >2:
        print tri_hists
        tri_hists["mugamma"].Draw("hist")
        tri_hists["muon"].Draw("same hist")
        #tri_hists["single-mu"].Draw("same hist")
        tri_hists["mugamma"].SetLineColor(kRed+3)
        #tri_hists["single-mu"].SetLineColor(kGreen+2)
        leg = TLegend(0.60,0.73,0.90,0.90);
        leg.AddEntry(tri_hists["mugamma"], "Mu22_Pho22","f")
        leg.AddEntry(tri_hists["muon"],    "Mu17_Mu8",  "f")
        #leg.AddEntry(tri_hists["single-mu"],"IsoMu24",  "f")
        leg.SetTextSize(0.04)
        leg.SetFillColor(kWhite)
        leg.Draw()
    
        c1.SaveAs("tri_plot.png")
    '''


    table_sig  = yieldsTable(yields_sig)
    table_data = yieldsTable(yields_data)

    makeTable(table_data,"twiki")
    makeTable(table_sig,"twiki")

    comments = ["These plots are made for "+thissel+" selection.",
                "Blah"]
    
    ht.makeHTML("h &rarr; dalitz decay plots",pathBase, plot_types, comments, "mugamma")

    print "\n\t\t finita la comedia \n"
