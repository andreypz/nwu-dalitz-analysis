#! /usr/bin/env python
from optparse import OptionParser
import sys,os,datetime
from array import *
from ROOT import *
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
parser.add_option("--mcfm", dest="mcfm", action="store_true", default=False, help="Use MCFM  as a signal")
parser.add_option("--noeos",dest="noeos",action="store_true", default=False, help="Don't use EOS. pick up the files from nobackup area")

parser.add_option("--apz",dest="apz",action="store_true", default=False, help="Discover new particle (requires apzTree produced)")
parser.add_option("--zjp",dest="zjp",action="store_true", default=False, help="Study J/Psi region and more (requires apzTree)")
parser.add_option("--hjp",dest="hjp",action="store_true", default=False, help="Study J/Psi region and more (requires apzTree)")
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
  if options.noeos:
    hPath  = "/uscms_data/d2/andreypz/zgamma/"+ver

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
    if sel =="elgamma":
      os.system("hadd "+hPath+"/m_Data_" +sel+"_"+period+".root "+hPath+"/"+sel+"_"+period+"/hhhh_*Run2012D*.root")
    else:
      # os.system("hadd "+hPath+"/m_Data_" +sel+"_"+period+".root "+hPath+"/"+sel+"_"+period+"/hhhh_DoubleMu_Run20*.root")
      os.system("hadd "+hPath+"/m_Data_" +sel+"_"+period+".root "+hPath+"/"+sel+"_"+period+"/hhhh_MuEG_Run20*.root")
      # os.system("hadd "+hPath+"/m_Data_" +sel+"_"+period+".root "+hPath+"/"+sel+"_"+period+"/hhhh_DoubleElectron_Run20*.root")


    if doBkg:
      os.system("hadd "+hPath+"/m_DYJets50_"   +sel+"_"+period+".root "
                +hPath+"/"+sel+"_"+period+"/hhhh_DYJetsPow20-RD1*.root")
      #os.system("hadd "+hPath+"/m_DYJets50_"   +sel+"_"+period+".root "
      #          +hPath+"/"+sel+"_"+period+"/hhhh_DYJets50*.root")
      os.system("hadd "+hPath+"/m_ZG_"    +sel+"_"+period+".root "
                +hPath+"/"+sel+"_"+period+"/hhhh_ZGToLLG_*.root ")

  if doBkg:
    bkgFiles.append(TFile(hPath+"/m_ZG_"+sel+"_"+period+".root","OPEN"))
    bkgNames.append('ZG')
    #bkgFiles.append(TFile(hPath+"/m_DYJets50_"+sel+"_"+period+".root","OPEN"))
    #bkgNames.append('DYJets50')
    bkgFiles.append(TFile(hPath+"/"+sel+"_"+period+"/hhhh_ZGDalitz_1.root","OPEN"))
    bkgNames.append('ZGDalitz')
    bkgFiles.append(TFile(hPath+"/"+sel+"_"+period+"/hhhh_DYJetsDalitz_1.root","OPEN"))
    bkgNames.append('DYJetsDalitz')

    #yields_bkg  = u.getYields(bkgFiles,"DY",True)

  bkgZip = zip(bkgNames, bkgFiles)

  sigFileMAD  = TFile(hPath+"/"+sel+"_"+period+"/hhhh_ggH-mad"+str(mass)+"_1.root", "OPEN")
  sigFileVBF  = TFile(hPath+"/"+sel+"_"+period+"/hhhh_vbf-mad"+str(mass)+"_1.root", "OPEN")
  sigFileVH   = TFile(hPath+"/"+sel+"_"+period+"/hhhh_vh-mad"+str(mass)+"_1.root",  "OPEN")

  sigFileZjp  = TFile(hPath+"/"+sel+"_"+period+"/hhhh_ZtoJPsiGamma_1.root",  "OPEN")
  sigFileHjp  = TFile(hPath+"/"+sel+"_"+period+"/hhhh_HiggsToJPsiGamma_1.root",  "OPEN")

  if options.mcfm:  sigFile = sigFileMCFM
  elif options.hjp: sigFile = sigFileHjp
  elif options.zjp: sigFile = sigFileZjp
  else:             sigFile = sigFileMAD

  dataFile = TFile(hPath+"/m_Data_"+sel+"_"+period+".root","OPEN")

  if sel=="zee":
    u.drawAllInFile(dataFile, "data",bkgZip,None,"",
                    "Zee",  pathBase+"/Zee",  None,"norm")

  if options.zjp:
    u.setSelection('mugamma')

  yields_data = u.getYields(dataFile)
  yields_zjp  = u.getYields(sigFileZjp, 'ZtoJPsiGamma',  True)
  yields_hjp  = u.getYields(sigFileHjp, 'HtoJPsiGamma',  True)
  yields_ggH  = u.getYields(sigFile,    'ggH-125',  True)
  yields_vbf  = u.getYields(sigFileVBF, 'vbfH-125', True)
  yields_vh   = u.getYields(sigFileVH,  'vH-125',   True)
  yields_sig  = [sum(x) for x in zip(yields_ggH,yields_vbf,yields_vh)]

  print 'ggH yi', yields_ggH
  print 'sig yi', yields_sig

  subdir = sel
  path = pathBase+"/"+subdir
  #if doBkg:
  #  path = pathBase+"/bkg_"+subdir
  #  u.createDir(path)
  #path = pathBase+"/"+subdir


  sigName = '#splitline{100xSignal}{m_{H}=125 GeV}'
  if options.zjp: sigName= '#splitline{50xSignal}{Z #rightarrow J/Psi #gamma}'
  if options.hjp: sigName= '#splitline{100xSignal}{h #rightarrow J/Psi #gamma}'

  if (options.zjp or options.hjp) and cut in ['12']:
    u.drawAllInFile(dataFile, "Data", bkgZip, sigFile, sigName,  "jpsi",path, cut, "lumi")


  if cut not in ['12','14','15']:
    #u.drawAllInFile(dataFile, "Data", bkgZip, None, '', "", path, cut, "lumi")
    u.drawAllInFile(dataFile, "Data", bkgZip, sigFile, sigName, "", path, cut, "lumi")
    u.drawAllInFile(dataFile, "Data", bkgZip, sigFile, sigName, "Angles", pathBase+'/Angles', cut, "norm1")
    u.drawAllInFile(dataFile, "Data", bkgZip, sigFile, sigName, "Muons",  pathBase+'/Muons', None, "norm1")
    #u.drawAllInFile(dataFile, "Data", bkgZip, sigFile, sigName,  "",path, cut, "norm")
    # u.drawAllInFile(dataFile, "data", bkgZip, sigFile,"50xSignal","EB",pathBase+"/EB", cut, "lumi")
    # u.drawAllInFile(dataFile, "data", bkgZip, sigFile,"50xSignal","EE",pathBase+"/EE", cut, "lumi")

  if sel == "mugamma" and not options.zjp:
    if cut not in ['12','14','15','16']:
      # u.drawAllInFile(dataFile, "data",bkgZip,sigFile,"signal",  "Muons", pathBase+"/Muons", None,"norm")
      u.drawAllInFile(dataFile, "data",bkgZip,sigFile,"signal",
                      "Muons", pathBase+"/Muons/", None,"norm")
      u.drawAllInFile(dataFile, "data",bkgZip,sigFile,"#splitline{Signal}{m_{H}=125 GeV}",
                      "Photon",pathBase+"/Photon/", None,"norm")

    elif cut in ['12']:
      u.drawAllInFile(dataFile, "data",bkgZip,sigFile,"#splitline{Signal}{m_{H}=125 GeV}",
                      "jpsi",pathBase+"/jpsi/", cut,"norm2")
    elif cut in ['14','15','16']:
      u.drawAllInFile(dataFile, "data",bkgZip,None,"",
                      "AlphaPiZ",pathBase+"/apz/", cut,"norm2")

  elif sel == "elgamma":
    u.drawAllInFile(dataFile, "data",bkgZip,sigFile,"signal",
                    "Photon",     pathBase+"/Photon",     None,"norm")

    u.drawAllInFile(dataFile, "data",bkgZip,sigFile,"signal",
                    "DalitzEle-Before",  pathBase+"/DalitzEle-Before",  None,"norm")
    u.drawAllInFile(dataFile, "data",bkgZip,sigFile,"signal",
                    "DalitzEle-AfterAll",  pathBase+"/DalitzEle-AfterAll",  None,"norm")
    u.drawAllInFile(dataFile, "data",bkgZip,sigFile,"signal",
                    "DalitzEle-Before_tracks",  pathBase+"/DalitzEle-Before_tracks",  None,"norm")
    u.drawAllInFile(dataFile, "data",bkgZip,sigFile,"signal",
                    "DalitzEle-AfterAll_tracks",  pathBase+"/DalitzEle-AfterAll_tracks",  None,"norm")


    u.drawAllInFile(dataFile, "data",bkgZip,sigFile,"signal",
                    "DalitzEle",  pathBase+"/DalitzEle",  None,"norm")

    u.drawAllInFile(dataFile, "data",bkgZip,sigFile,"signal",
                    "DalitzEleEB",pathBase+"/DalitzEleEB",None,"norm")
    u.drawAllInFile(dataFile, "data",bkgZip,sigFile,"signal",
                    "DalitzEleEE",pathBase+"/DalitzEleEE/",None,"norm")

    u.drawAllInFile(dataFile, "data",bkgZip,sigFile,"signal",
                    "DalitzEleEB_tracks",pathBase+"/DalitzEleEB_tracks",None,"norm")
    u.drawAllInFile(dataFile, "data",bkgZip,sigFile,"signal",
                    "DalitzEleEE_tracks",pathBase+"/DalitzEleEE_tracks",None,"norm")


    # dataFile.Close()
    #print yields_data


    #u.drawAllInFile(sigFileMAD, "signal", bkgZip, None,"","GEN-RECO",pathBase+"/GEN-RECO", None, "lumi")

    #u.drawAllInFile(bkgFiles, "DY electrons", sigFile, "Dalitz 2el",  "NewEle-1", pathBase+"/NewEle-1/", None,"norm", isLog=1, )

    #sigFileMAD  = TFile(hPath+"/mugamma_"+period+"/hhhh_ggH-mad125_1.root", "OPEN")

  ss = options.sel
  if options.apz and doMerge:
    if ss in ['4mu','2e2mu']:
      os.system("hadd -f "+hPath+"/m_Data_apz_DoubleMu_"+ss+"_"+period+".root "+hPath+"/"+ss+"_"+period+"/hhhh_DoubleMu_Run20*.root")

    if ss in ['2e2mu','mugamma']:
      os.system("hadd -f "+hPath+"/m_Data_apz_MuEG_"+ss+"_"+period+".root "+hPath+"/"+ss+"_"+period+"/hhhh_MuEG_Run20*.root")
    if ss =='2e2mu':
      os.system("hadd -f "+hPath+"/m_Data_apz_DoubleElectron_"+ss+"_"+period+".root "+hPath+"/"+ss+"_"+period+"/hhhh_DoubleElectron_Run20*.root")

  c1 = TCanvas("c4","small canvas",600,600);


  if options.apz:
    if ss in ['4mu','2e2mu']:
      if ss=='2e2mu':
        samples = ['DoubleMu','MuEG','DoubleElectron']
      elif ss=='4mu':
        samples = ['DoubleMu']

      for s in samples:
        data = TFile(hPath+"/m_Data_apz_"+s+"_"+ss+"_2012.root","OPEN")

        alphaPiZ2(data, c1, TCut(''),   pathBase+"/"+s+"-alphaPiZ-1/")
        alphaPiZ2(data, c1, TCut('pt12/m4l>0.3 && pt34/m4l>0.3'),                   pathBase+"/"+s+"-alphaPiZ-2/")
        alphaPiZ2(data, c1, TCut('pt12/m4l>0.3 && pt34/m4l>0.3 && m4l>100'),        pathBase+"/"+s+"-alphaPiZ-3/")
        alphaPiZ2(data, c1, TCut('pt12/m4l>0.3 && pt34/m4l>0.3 && m4l>100 && m12>15 && m12<30'),pathBase+"/"+s+"-alphaPiZ-4/")
        alphaPiZ2(data, c1, TCut('m12>15 && m12<30 && m34>15 && m34<30'),pathBase+"/"+s+"-alphaPiZ-5/")

        alphaPiZ2(data, c1, TCut('m12>15 && m12<30 && m34>15 && m34<30 && pt12>20 && pt34>20'),pathBase+"/"+s+"-alphaPiZ-6/")

        alphaPiZ2(data, c1, TCut('m4l>140 && m4l<150'),pathBase+"/"+s+"-alphaPiZ-7/")


        u.drawAllInFile(data, "data",bkgZip,None,"",
                        "apz-plots",pathBase+"/apz-plots/",None,"norm2")


    elif ss=='mugamma':

      data = TFile(hPath+"/m_Data_mugamma_2012.root","OPEN")

      alphaPiZ(data, c1, TCut('ph_pt/m_llg>0.3 && di_pt/m_llg>0.3 && m_llg>100&&m_llg<170'),
               pathBase+"/alphaPiZ-0/")
      alphaPiZ(data, c1, TCut(''),                                                            pathBase+"/alphaPiZ-1/")
      alphaPiZ(data, c1, TCut('ph_pt/m_llg>0.30 && di_pt/m_llg>0.30'),                          pathBase+"/alphaPiZ-2/")
      alphaPiZ(data, c1, TCut('ph_pt/m_llg>0.30 && di_pt/m_llg>0.30 && m_llg>100&&m_llg<150'),  pathBase+"/alphaPiZ-3/")
      alphaPiZ(data, c1, TCut('ph_pt/m_llg>0.35 && di_pt/m_llg>0.35 && m_llg>100&&m_llg<150'),pathBase+"/alphaPiZ-4/")
      alphaPiZ(data, c1, TCut('ph_pt/m_llg>0.35 && di_pt/m_llg>0.35 && m_llg>120&&m_llg<180'),pathBase+"/alphaPiZ-5/")
      alphaPiZ(data, c1, TCut('ph_pt/m_llg>0.35 && di_pt/m_llg>0.35 && m_llg>100&&m_llg<130'),pathBase+"/alphaPiZ-6/")

      alphaPiZ(data, c1, TCut('ph_pt/m_llg>0.30 && di_pt/m_llg>0.30 && m_llg>121&&m_llg<131'),pathBase+"/alphaPiZ-7/")
      alphaPiZ(data, c1, TCut('ph_pt/m_llg>0.30 && di_pt/m_llg>0.30 && m_llg>85&&m_llg<96'),
               pathBase+"/alphaPiZ-8/")



    '''
    if sel[0]=='none':
        sigFileMAD  = TFile("hhhh_ggH-mad125_1.root", "OPEN")

        u.createDir(pathBase+"/eff")
        u.effPlots(sigFileMAD, pathBase+"/eff/")
        u.effPlots2(sigFileMAD, pathBase+"/eff/")



    htmp = sigFileMAD.Get("lPt2_pt__cut10")
    int1 = htmp.Integral()
    int2 = htmp.Integral(0,10)
    print '\n Fraction of event with trailing lepton pt<20: \n', float(int2)/int1

    htmp = sigFileMAD.Get("ll_deltaR__cut10")
    int1 = htmp.Integral()
    int2 = htmp.Integral(0,10)
    print '\n Fraction of event with dR(l1,l2)<0.4: \n', float(int2)/int1
    '''
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

    hc7 = dataFile.Get("tri_mass80__cut7").Integral(36,40)
    hc8 = dataFile.Get("tri_mass80__cut8").Integral(36,40)
    hc9 = dataFile.Get("tri_mass80__cut9").Integral(36,40)
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
        treeda = dataFile.Get("fitTree/fitTree")
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
    h2da = dataFile.Get("h2D_dalitzPlot_rotation__cut"+cut).ProjectionX("hda_prx")
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
  if options.zjp and cut==4:
    data = TFile(hPath+"/m_Data_"+sel+"_2012.root","OPEN")

    ZJPG(data, c1, TCut('pt3/m123>0.3 && pt12/m123>0.3 && m123>110 && m123<170 && dr13>1 && dr23>1'),
         pathBase+"/ZtoJPsiGamma-0/")

    ZJPG(data, c1, TCut(''), pathBase+"/ZtoJPsiGamma-1-full/")
    ZJPG(data, c1, TCut('m12<30'), pathBase+"/ZtoJPsiGamma-2-Mll30/")

    ZJPG(data, c1, TCut('pt3/m123>0.3 && pt12/m123>0.3 && m12<30 && dr1234>1'), pathBase+"/ZtoJPsiGamma-3-ptmllg03/")

    comments.append('ZtoJPsiGamma-0 is Dalitz-like selection (as in HIG-14-003 analysis)')
    comments.append('ZtoJPsiGamma-1 are full events (even without Mll cut)')
    comments.append('ZtoJPsiGamma-2 is with Mll<30 GeV cut')
    comments.append('ZtoJPsiGamma-3: pt3/m123>0.3 && pt12/m123>0.3 && m12<30 && dr1234>1 ')


  plot_types =[]
  list = os.listdir(pathBase)
  for d in list:
    if os.path.isdir(pathBase+"/"+d):
      plot_types.append(d)

    #plot_types
  print yields_sig

  if doBkg:
    #table_all  = u.yieldsTable([yields_data,yields_bkg,yields_sig, yields_ggH,yields_vbf, yields_vh], sel)
    table_all = u.yieldsTable([yields_data, yields_hjp, yields_zjp], sel)
  else:
    if options.zjp or options.hjp:
      table_all = u.yieldsTable([yields_data, yields_hjp, yields_zjp], sel)
    elif not options.apz:
      table_all  = u.yieldsTable([yields_data,yields_sig, yields_ggH,yields_vbf, yields_vh], sel)

  u.makeTable(table_all,"all", "html")
  u.makeTable(table_all,"all", "twiki")
  u.makeTable(table_all,"all", "tex")

  os.system("cat yields_all.html   > yields.html")
  #os.system("cat yields_all.twiki  > yields.html")
  #os.system("cat yields_all.tex    > yields.html")

  defaultPage = sel
  if cut in ['12']: defaultPage = 'jpsi'
  elif cut in ['14','15']: defaultPage = 'apz'
  #elif options.zjp or options.hjp: defaultPage = 'jp-'+sel
  print defaultPage


  ht.makeHTML("h &rarr; dalitz decay plots",pathBase, plot_types, comments, defaultPage)

  print "\n\t\t finita la comedia \n"
