#!/usr/bin/env python
import sys,os

from ROOT import *
import shutil
import datetime

import config as c

nargs = len(sys.argv)
print sys.argv[0], nargs

sel = 1
doMerge = 0

if (nargs<2):
    print "you have to specify the path where to llok for the files!\nLike v62, or something"
    sys.exit()
if (nargs==3):
    sel = sys.argv[2]

hPath = sys.argv[1]
print hPath, sel


def createDir(dir):
    try:
        os.makedirs(dir)
    except OSError:
        if os.path.isdir(dir):
            pass
        else:
            raise

dirnameOut = "/afs/fnal.gov/files/home/room1/andreypz/public_html/higgs/"+hPath
selection  = ['muon', 'electron']
plot_types = ['diLepton', 'Lepton', 'Jet', 'Met', 'Special', 'Misc']

for x in selection:
    for y in plot_types:
        createDir(dirnameOut+'/'+x+'/'+y)
        
gROOT.LoadMacro("./makePlot.C");
gROOT.LoadMacro("./utils.C");
  
timer = TStopwatch()	
timer.Start()

#TString ssel("none"), gsel("none");
#if (sel==1)  {ssel = "muon";     gsel ="muGamma";}
#    if (sel==2)  {ssel = "electron"; gsel ="eGamma";}

#TString histoPath = Form("%s/%s", hPath.Data(), ssel.Data());
#TString gammaPath = Form("%s/%s", hPath.Data(), gsel.Data());
#cout<<histoPath.Data()<<endl;

 
gROOT.SetBatch()
if (doMerge==0):    
    gROOT.ProcessLine(".x makePlot.C("+str(sel)+", \""+hPath+"\")");
    
    print "\n\nDone!"
    print "CPU Time : ", timer.CpuTime()
    print "RealTime : ", timer.RealTime()
  
else:
    #Same for now, but later do the merging here (t+tW+tt, Data etc)
    gROOT.ProcessLine(".x makePlot.C("+str(sel)+", \""+hPath+"\")");
    
    print "\n\nDone!"
    print "CPU Time : ", timer.CpuTime()
    print "RealTime : ", timer.RealTime()


#os.system("rm ./%s/m_*%i.root "% (hPath.Data(), sel))

print "\n\n ******** Now making HTML pages ******** \n"
menu=""

for x in plot_types:
    fname = x+".html"
    imgfile = open(fname,"w")
    imgfile.write("<html><head><title>"+x+"</title></head><body>\n")
    fileList = {}
    for s in selection:
        fileList[s] = os.listdir(dirnameOut+'/'+s+'/'+x)
    #print fileList
    nmu = len(fileList["muon"])
    nele = len(fileList["electron"])
    if (nmu!=0 and fileList["muon"]==fileList["electron"]):
        for pl in  fileList["muon"]:
            imgfile.write('<nobr><img src=muon/'+x+'/'+pl+' width=47%>')
            imgfile.write('    <img src=electron/'+x+'/'+pl+' width=47%></nobr>\n') 
    elif(nmu!=0 and  nele!=0 and fileList["muon"]!=fileList["electron"]):
        print "Something is wrong: muons and electrons are not symmetric!"
    elif(nmu!=0 and nele==0):
        count =1
        mod =1
        for pl in  fileList["muon"]:
            mod = count % 2
            if mod==1: imgfile.write('<nobr><img src=muon/'+x+'/'+pl+' width=47%>')
            if mod==0: imgfile.write('      <img src=muon/'+x+'/'+pl+' width=47%></nobr>\n')
            count+=1
        if mod==0: imgfile.write("")
        if mod==1: imgfile.write("</nobr>")
    elif(nmu==0):
        print "No plots in", x
        
    imgfile.write("</body></html>")
    imgfile.close()
    os.system("mv "+fname+" "+dirnameOut)
    menu = menu+"<li><a href=\""+x+".html\" target=\"iframe_a\">"+x+"</a></li>"

#print menu

today = datetime.date.today()
print today

tempfile = open("indextemplate.html","r")
whole_thing = tempfile.read()
whole_thing = whole_thing.replace("{MENU}", menu)
whole_thing = whole_thing.replace("{DATE}", str(today))
tempfile.close()

ifile = open("index.html","w")
ifile.write(whole_thing)
ifile.close()

os.system("mv index.html "+dirnameOut)
