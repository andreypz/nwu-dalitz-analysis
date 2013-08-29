#! /usr/bin/env python
from optparse import OptionParser
import sys,os,datetime
from ROOT import *
gROOT.SetBatch()

parser = OptionParser(usage="usage: %prog ver [options -c, -e, -p, -m]")
parser.add_option("-c","--cut",    dest="cut", type="int", default=3,    help="Plots after a certain cut")
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
cuts = []
sel = ["mugamma","muon","electron"]

for key, cut in sorted(conf.items("cuts_stoyan")):
    cuts.append(cut)
    #print key, cut
def draw(h1, path):
    h = h1
    if h.InheritsFrom("TH2"):
        h.Draw("lego")
    else:
        h.Draw("hist")

    c1.SaveAs(path+h.GetName()+".png")


def drawAllInFile(f, fsig, dir,path, N):
    f.cd(dir)
    dirList = gDirectory.GetListOfKeys()
    dirList.Print()
    
    for k1 in dirList:
        if k1.GetName() in ["Counts"]: continue
        if N!=None:
            if not("cut"+N) in k1.GetName(): continue
        h1 = k1.ReadObj()

        if fsig!=None:
            if dir!="":
                hsig = fsig.Get(dir+"/"+k1.GetName()) #assumes that the histograms in signal file have the same names  
            else:
                hsig = fsig.Get(k1.GetName())
        print "drawing", h1.GetName()

        if h1.InheritsFrom("TH2"):
            h1.Draw("col")
        else:
            h1.Draw("hist")
            leg = TLegend(0.70,0.8,0.90,0.90);
            leg.AddEntry(h1,"data", "l")
            leg.SetTextSize(0.04)
            if fsig!=None and hsig!=None:
                hsig.Draw("same hist")
                hsig.SetLineColor(kRed)
                norm1 = h1.Integral()
                norm2 = hsig.Integral()
                #print "Integrals = ", norm1, norm2
                hsig.Scale(float(norm1)/2/float(norm2))
                leg.AddEntry(hsig,"signal", "l")
            if h1.GetName() in ["reco_gen_l1_deltaR", "reco_gen_l2_deltaR", "gen_ll_deltaR",
                                "reco_gen_gamma_deltaR", "ll_deltaR"]:
                c1.SetLogy()
                
            leg.SetFillColor(kWhite)
            leg.Draw()

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
                                            
def makeHTML(title, htmlDir, plot_types, description, IFRAMEA):

    print "\n\n ******** Now making HTML pages ******** \n"
    menu=""
    #plot_types = ["comb","h_dist1","h_dist2","h_dist3","fit","4mu","4e"]

    fileList = {}
    for x in plot_types:
        newDir = htmlDir+"/"+x
        createDir(newDir)
        #if not os.path.exists(newDir):
        #    print "Creating ", newDir
        #    os.makedirs(newDir)


        fname = x+".html"
        imgfile = open(fname,"w")
        imgfile.write("<html><head><title>"+x+"</title></head><body>\n")

        fileList[x] = os.listdir(newDir)
        #print fileList[x]
        
        count =1
        mod =1
        for pl in  sorted(fileList[x]):
            mod = count % 2
            if mod==1: imgfile.write('<nobr><img src='+x+'/'+pl+' width=45%>')
            if mod==0: imgfile.write('      <img src='+x+'/'+pl+' width=45%></nobr>\n')
            count+=1
        if mod==0: imgfile.write("")
        if mod==1: imgfile.write("</nobr>")

        imgfile.write("</body></html>")
        imgfile.close()
        os.system("mv "+fname+" "+htmlDir)

        menu = menu+"<li><a href=\""+x+".html\" target=\"iframe_a\">"+x+"</a></li>"

    menu += "<li><a href=\"yields.html\" target=\"iframe_a\">yields</a></li>"

    today = datetime.date.today()
    print today

    message = '<h2>Comments</h2>'
    message+='<ul>'
    for d in description:
        message+="<li>"+d+" </li>"
        message+='</ul>'

    tempfile = open("indextemplate.html","r")
    whole_thing = tempfile.read()
    whole_thing = whole_thing.replace("{TITLE}", title)
    whole_thing = whole_thing.replace("{MENU}", menu)
    whole_thing = whole_thing.replace("{DATE}", str(today))
    whole_thing = whole_thing.replace("{ASIDEMESSAGE}", str(message))
    whole_thing = whole_thing.replace("{IFRAMEA}", IFRAMEA)
    tempfile.close()


    ifile = open(htmlDir+"/index.html","w")
    ifile.write(whole_thing)
    ifile.close()

    os.system("cp yields.html "+htmlDir)
    
    print "\n\n *** End of  making HTML pages - all done *** \n"


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
                myTable+="%.0f" % (table[l][c])
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

def getYields(f):
    ev = f.Get("Counts/evt_byCut")

    y = []
    for a in xrange(len(cuts)):
        y.append(ev.GetBinContent(a+1)) # well, that's how the histogram is set up
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


    if doMerge:
        os.system("rm "+hPath+"/m_*.root") #removing the old merged files
    yields_data = {}
    yields_bkg = {}
    yields_sig  = {}

    tri_hists = {}
    h1 = TH1F()
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
        bkgFile = TFile(hPath+"/"+thissel+"_"+period+"/hhhh_DY-mg5_1.root",   "OPEN")
        dataFile[thissel] = TFile(hPath+"/m_Data_"+thissel+"_"+period+".root","OPEN")
        
        yields_data[thissel] = getYields(dataFile[thissel])
        yields_sig[thissel]  = getYields(sigFile)
        #yields_bkg[thissel]  = getYields(bkgFile)

        if int(cut) >2:
            tri_hists[thissel]   = dataFile[thissel].Get("tri_mass_cut"+cut).Clone()
        
        if doBkg:
            drawAllInFile(bkgFile, sigFile, path, cut)
        else:
            drawAllInFile(dataFile[thissel], sigFile, "",path, cut)
            if thissel =="mugamma":
                drawAllInFile(dataFile[thissel], sigFile, "Muons",pathBase+"/Muons/", None)
        
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
    #table = [["Cut/trigger", ]
    table = [["Cut/trigger", "Mu22_Pho22","mumu","phopho"]]

    for line in xrange(len(cuts)):
        l=[]
        l.append(cuts[line])
        for thissel in sel:
            #
            #l.append(yields_sig[thissel][line])
            if doBkg:
                l.append(yields_bkg[thissel][line])
            else:
                #l.append(yields_data[thissel][line])
                l.append(yields_sig[thissel][line])
        table.append(l)

        
    makeTable(table,"twiki")

    comments = ["These plots are made for "+thissel+" selection.",
                "Blah"]
    
    makeHTML("h &rarr; dalitz decay plots",pathBase, plot_types, comments, "mugamma")

    print "\n\t\t finita la comedia \n"
