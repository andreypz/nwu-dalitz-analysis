#! /usr/bin/env python
from optparse import OptionParser
import sys,os,datetime
from ROOT import *
gROOT.SetBatch()

parser = OptionParser(usage="usage: %prog ver [options -c, -e, -p, -m]")
parser.add_option("-c","--cut", type="int", dest="cut", default=3, help="Plots after a certain cut")
parser.add_option("--fsig", type="string", dest="fsig", default=None, help="A file for signal")
parser.add_option("-e", "--ele",   dest="ele", action="store_true", default=False, help="Use electron selection")
parser.add_option("-p", "--period",dest="period", default="2012", help="Year period; 2011 or 2012")
parser.add_option("-m", "--merge", action="store_true", dest="merge", default=False, help="Do merging?")

(options, args) = parser.parse_args()

def draw(h1, path):
    h = h1
    if h.InheritsFrom("TH2"):
        h.Draw("lego")
    else:
        h.Draw("hist")

    c1.SaveAs(path+h.GetName()+".png")


def drawAllInFile(f, fsig, path, N):
    f.cd()
    dirList = gDirectory.GetListOfKeys()
    
    for k1 in dirList:
        if k1.GetName() in ["Counts"]: continue
        if not ("cut"+N) in k1.GetName(): continue
        h1 = k1.ReadObj()

        if fsig!=None:
            hsig = fsig.Get(k1.GetName()).Clone() #assumes that the histograms in signal file have the same names  
        print "drawing", h1.GetName()

        if h1.InheritsFrom("TH2"):
            h1.Draw("col")
        else:
            h1.Draw("hist")
            leg = TLegend(0.70,0.8,0.90,0.90);
            leg.AddEntry(h1,"data", "l")
            leg.SetTextSize(0.04)
            if fsig!=None:
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
                                            
def makeHTML(title, htmlDir, plot_types, description,IFRAMEA):

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

    print "\n\n *** End of  making HTML pages - all done *** \n"



if __name__ == "__main__":
    timer = TStopwatch()
    timer.Start()
    
    if len(args) < 1:
        parser.print_usage()
        exit(1)
        
    ver    = sys.argv[1]
    #subdir = sys.argv[3]        
    cut=str(options.cut)
    fsig = options.fsig
    doMerge = options.merge
    period = options.period
    
    gROOT.ProcessLine(".L ~/tdrstyle.C")
    setTDRStyle()
    TH1.SetDefaultSumw2(kTRUE)
    
    sigFile=None
    if fsig!=None:
        sigFile = TFile(fsig,"OPEN")
      
    
    pathBase = "/uscms_data/d2/andreypz/html/zgamma/dalitz/"+ver+"_cut"+cut
    hPath = "/eos/uscms/store/user/andreypz/batch_output/zgamma/8TeV/"+ver


    if doMerge:
        os.system("rm "+hPath+"/m_*.root") #removing the old merged files
    for thissel in ["muon","electron","mugamma","single-mu"]:
        if doMerge:
            os.system("hadd "+hPath+"/m_Data_"    +thissel+"_"+period+".root "+hPath+"/"+thissel+"_"+period+"/hhhh_*Run20*.root")
            
        subdir = thissel
        path = pathBase+"/"+subdir+"/"
        createDir(path)

        sigFile = TFile(hPath+"/"+thissel+"_"+period+"/hhhh_h-dalitz_1.root")
        inFile = TFile(hPath+"/m_Data_"+thissel+"_"+period+".root","OPEN")
        drawAllInFile(inFile, sigFile, path, cut)
    
        inFile.Close()

    plot_types =[]
    list = os.listdir(pathBase)
    for d in list:
        if os.path.isdir(pathBase+"/"+d):
            plot_types.append(d)
            
    comments = ["These plots are made for "+thissel+" selection."]

    makeHTML("h &rarr; dalitz decay plots",pathBase, plot_types, comments, "muon")
    

    print "\n\t\t finita la comedia \n"
