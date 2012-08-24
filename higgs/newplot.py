#!/usr/bin/env python
import sys,os

from ROOT import *
import shutil
import datetime

import config as c
import makePlots as mp

nargs = len(sys.argv)
print sys.argv[0], "nargs:", nargs

sel = 1
doMerge = False

if (nargs<2):
    print "you have to specify the path where to llok for the files!\nLike v62, or something"
    sys.exit()
if (nargs>=3):
    sel = int(sys.argv[2])
if (nargs==4):
    doMerge = bool(int(sys.argv[3]))

hPath = sys.argv[1]
print hPath, sel, doMerge

def createDir(dir):
    try:
        os.makedirs(dir)
    except OSError:
        if os.path.isdir(dir):
            pass
        else:
            raise

#baseDir = "/afs/fnal.gov/files/home/room1/andreypz/public_html/test2/"
baseDir = "/uscms_data/d2/andreypz/hzz2l2nu_html/"
dirnameOut = baseDir+hPath
selection  = ['muon', 'electron']
plot_types = ['diLepton', 'Jet', 'Met', 'Misc', "mvaPresel"]
#plot_types = ['diLepton', 'Lepton', 'Jet', 'Met', 'Special', 'Misc', "mvaPresel"]

thissel = selection[sel-1]

for x in selection:
    for y in plot_types:
        createDir(dirnameOut+'/'+x+'/'+y)
        
#gROOT.LoadMacro("./makePlot.C");
  
timer = TStopwatch()	
timer.Start()

gROOT.SetBatch()
print "Merging?", doMerge
if doMerge:
    os.system("rm "+hPath+"/m_*_"+thissel+".root") #removing the old merged files
    os.system("hadd ./"+hPath+"/m_Data_"+thissel+".root  ./"+hPath+"/"+thissel+"/hhhh_Double*.root")
    os.system("hadd ./"+hPath+"/m_ttbar_"+thissel+".root ./"+hPath+"/"+thissel+"/hhhh_ttbar_*.root")
    os.system("hadd ./"+hPath+"/m_DYjets_"+thissel+".root ./"+hPath+"/"+thissel+"/hhhh_DYjets_*.root")
                

mp.makePlots(sel, baseDir, hPath)

print "\n\nDone!"
print "CPU Time : ", timer.CpuTime()
print "RealTime : ", timer.RealTime()  


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
        for pl in  sorted(fileList["muon"]):
            imgfile.write('<nobr><img src=muon/'+x+'/'+pl+' width=47%>')
            imgfile.write('    <img src=electron/'+x+'/'+pl+' width=47%></nobr>\n') 
    elif(nmu!=0 and  nele!=0 and fileList["muon"]!=fileList["electron"]):
        print "Something is wrong: muons and electrons are not symmetric!"
    elif(nmu!=0 and nele==0):
        count =1
        mod =1
        for pl in  sorted(fileList["muon"]):
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

menu += '<li><a href="yields_muon.html" target="iframe_a">Yields mu</a></li>'
menu += '<li><a href="yields_electron.html" target="iframe_a">Yields ele</a></li>'
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

print "\n\n *** End of  making HTML pages - all done *** \n"

os.system("mv index.html "+dirnameOut)
os.system("mv yields* "+dirnameOut)
