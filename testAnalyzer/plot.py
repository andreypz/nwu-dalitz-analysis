#! /usr/bin/env python
from optparse import OptionParser
import sys,os,datetime
from array import *
from ROOT import *
sys.path.append("../zgamma")
import utils as u
from apz import *
import makeHTML as ht
gROOT.SetBatch()

parser = OptionParser(usage="usage: %prog ver [options -c, -e, -p, -m]")
parser.add_option("-c","--cut",   dest="cut", type="int", default=4,    help="Plots after a certain cut")
parser.add_option("--mass", dest="mass", type="int", default=125,    help="Signal sample (mass)")
parser.add_option("-m","--merge", dest="merge",action="store_true", default=False, help="Do merging?")

parser.add_option("-p", "--period",dest="period", default="2012",  help="Year period; 2011 or 2012")
parser.add_option("--bkg",  dest="bkg",  action="store_true", default=False, help="Make plots from bkg sample")

parser.add_option("--zjp",dest="zjp",action="store_true", default=False, help="Study J/Psi region and more (requires apzTree)")
parser.add_option("-s", '--sel', dest="sel", type="string", default='mugamma',
                  help="Selection to be used. Options are: '4mu','2e2mu', 'zee','mugamma', 'egamma'")

(options, args) = parser.parse_args()

mass = options.mass
sel = options.sel

comments = ["These plots are made for z -> J/Psi gamma analysis",
            "MuEG dataset used"]

if __name__ == "__main__":
  timer = TStopwatch()
  timer.Start()

  if len(args) < 1:
    parser.print_usage()
    exit(1)

  ver    = sys.argv[1]
  if 'vv/' in ver: ver = ver[3:].rstrip('/')
  cut=str(options.cut)
  doMerge = options.merge
  period  = options.period
  doBkg   = options.bkg

  gROOT.ProcessLine(".L ../tdrstyle.C")
  setTDRStyle()
  TH1.SetDefaultSumw2(kTRUE)

  pathBase = "/uscms_data/d2/andreypz/html/zgamma/dalitz/"+ver+"_cut"+cut
  hPath    = "/eos/uscms/store/user/andreypz/batch_output/zgamma/8TeV/"+ver

  if '/tthome' in os.getcwd():
    pathBase = "/tthome/andrey/html/zgamma/dalitz/"+ver+"_cut"+cut
    hPath    = "/tthome/andrey/batch_output/zgamma/8TeV/"+ver

  u.createDir(pathBase)

  if doMerge:
    os.system("rm "+hPath+"/m_*.root") #removing the old merged files

  dataFile  = None
  bkgFiles  = []
  bkgNames  = []

  u.setSelection(sel)

  if doMerge:
    os.system("hadd "+hPath+"/m_Data_" +sel+"_"+period+".root "+hPath+"/"+sel+"_"+period+"/hhhh_DoubleMu_Run20*.root")
    # os.system("hadd "+hPath+"/m_Data_" +sel+"_"+period+".root "+hPath+"/"+sel+"_"+period+"/hhhh_MuEG_Run20*.root")
    # os.system("hadd "+hPath+"/m_Data_" +sel+"_"+period+".root "+hPath+"/"+sel+"_"+period+"/hhhh_DoubleElectron_Run20*.root")


    if doBkg:
      os.system("hadd "+hPath+"/m_DYJetsPow20_"+sel+"_"+period+".root "
                +hPath+"/"+sel+"_"+period+"/hhhh_DYJetsPow20-RD1*.root")
      os.system("hadd "+hPath+"/m_DYJets50_"+sel+"_"+period+".root "
                +hPath+"/"+sel+"_"+period+"/hhhh_DYJets50-RD1*.root")
      os.system("hadd "+hPath+"/m_ZG_"+sel+"_"+period+".root "
                +hPath+"/"+sel+"_"+period+"/hhhh_ZGToLLG-RD1*.root ")
      os.system("hadd "+hPath+"/m_QCD_"+sel+"_"+period+".root "
                +hPath+"/"+sel+"_"+period+"/hhhh_QCD_*.root ")

  if doBkg:
    #bkgFiles.append(TFile(hPath+"/m_ZG_"+sel+"_"+period+".root","OPEN"))
    #bkgNames.append('ZG')
    bkgFiles.append(TFile(hPath+"/m_DYJets50_"+sel+"_"+period+".root","OPEN"))
    bkgNames.append('DYJets50')
    bkgFiles.append(TFile(hPath+"/m_QCD_"+sel+"_"+period+".root","OPEN"))
    bkgNames.append('QCD')
    #bkgFiles.append(TFile(hPath+"/"+sel+"_"+period+"/hhhh_ZGDalitz-OLD_1.root","OPEN"))
    #bkgFiles.append(TFile(hPath+"/"+sel+"_"+period+"/hhhh_ZGDalitz_1.root","OPEN"))
    #bkgNames.append('ZGDalitz')
    #bkgFiles.append(TFile(hPath+"/"+sel+"_"+period+"/hhhh_DYJetsDalitz_1.root","OPEN"))
    #bkgNames.append('DYJetsDalitz')

    yields_bkg  = u.getYields(bkgFiles[0],"DYJets50",True)

    bkgZip = zip(bkgNames, bkgFiles)

  else: bkgZip = ''

  print bkgZip

  sigFile  = TFile(hPath+"/"+sel+"_"+period+"/hhhh_ZtoJPsiGamma_1.root",  "OPEN")

  dataFile = TFile(hPath+"/m_Data_"+sel+"_"+period+".root","OPEN")

  yields_data = u.getYields(dataFile)
  yields_zjp  = u.getYields(sigFile, 'ZtoJPsiGamma',  True)

  subdir = sel
  path = pathBase+"/"+subdir

  #sigName= '#splitline{50xSignal}{Z #rightarrow J/Psi #gamma}'

  u.drawAllInFile(dataFile, "Data", bkgZip, None, '', "", path, cut, "lumi")
  #u.drawAllInFile(dataFile, "Data", bkgZip, sigFile, sigName, "", path, cut, "lumi")
  #u.drawAllInFile(dataFile, "Data", bkgZip, sigFile, sigName, "Angles", pathBase+'/Angles', cut, "norm1")
  u.drawAllInFile(dataFile, "Data", bkgZip, None, '', "Muons",  pathBase+'/Muons', None, "norm1")


  plot_types =[]
  list = os.listdir(pathBase)
  for d in list:
    if os.path.isdir(pathBase+"/"+d):
      plot_types.append(d)

    #plot_types
  #print yields_sig
  if options.bkg:
    table_all  = u.yieldsTable([yields_data,yields_bkg], sel)
  else:
    table_all  = u.yieldsTable([yields_data], sel)

  u.makeTable(table_all,"all", "html")
  u.makeTable(table_all,"all", "twiki")
  u.makeTable(table_all,"all", "tex")

  os.system("cat yields_all.html   > yields.html")
  #os.system("cat yields_all.twiki  > yields.html")
  #os.system("cat yields_all.tex    > yields.html")

  defaultPage = sel

  print defaultPage


  ht.makeHTML("h &rarr; dalitz decay plots",pathBase, plot_types, comments, defaultPage)

  print "\n\t\t finita la comedia \n"
