#! /usr/bin/env python
import sys,os,datetime
from optparse import OptionParser
from ROOT import *
gROOT.SetBatch()

parser = OptionParser(usage="usage: %prog file.root ver subdir [options]")
parser.add_option("--s", dest="s", action="store_true", default=False, help="s")
(options, args) = parser.parse_args()

def draw(h1, path):
    h = h1#inFile.Get(h1)
    if h.InheritsFrom("TH2"):
        h.Draw("lego")
    else:
        h.Draw("hist")

    c1.SaveAs(path+h.GetName()+".png")


def drawAllInFile(f, path):
    f.cd()
    dirList = gDirectory.GetListOfKeys()
    
    for k1 in dirList:
        h1 = k1.ReadObj()
        #print "drawing", h1.GetName()

        if h1.InheritsFrom("TH2"):
            h1.Draw("col")
        else:
            h1.Draw("hist")
            if h1.GetName() in ["reco_gen_l1_deltaR", "reco_gen_l2_deltaR", "gen_ll_deltaR",
                                "reco_gen_gamma_deltaR", "ll_deltaR"]:
                c1.SetLogy()

        c1.SaveAs(path+h1.GetName()+".png")
        c1.SetLogy(0)


def makeHTML(title, htmlDir, plot_types, IFRAMEA):

    print "\n\n ******** Now making HTML pages ******** \n"
    menu=""
    #plot_types = ["comb","h_dist1","h_dist2","h_dist3","fit","4mu","4e"]

    fileList = {}
    for x in plot_types:
        newDir = htmlDir+"/"+x
        if not os.path.exists(newDir):
            print "Creating ", newDir
            os.makedirs(newDir)


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

    message = '<h2>Comments</h2>\
<ul>\
<li>These plots are </li>\
<li></li>\
<li></li>\
</ul>'

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
    
    file   = sys.argv[1]
    ver    = sys.argv[2]
    subdir = sys.argv[3]
    gROOT.ProcessLine(".L ~/tdrstyle.C")
    setTDRStyle()
    TH1.SetDefaultSumw2(kTRUE)
    
    inFile = TFile(file,"OPEN")
    inFile.cd()
    TH1.SetDefaultSumw2(kTRUE)
    
    
    plot_types =["mu","el"]
    
    pathBase = "/uscms_data/d2/andreypz/html/zgamma/dalitz/"
    path = pathBase+ver+"/"+subdir+"/"
    if not os.path.exists(path):
        os.makedirs(path)
        
        
    drawAllInFile(inFile,path)
    
    inFile.Close()
    
    makeHTML("h &rarr; dalitz decay plots",pathBase+ver, plot_types, "mu")
    

    print "\n\t\t finita la comedia \n"
