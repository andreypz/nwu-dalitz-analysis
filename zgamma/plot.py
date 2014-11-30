#! /usr/bin/env python
from optparse import OptionParser
import sys,os,datetime
from array import *
from ROOT import *
import utils as u
from apz import *
sys.path.append("../scripts")
import makeHTML as ht

gROOT.SetBatch()

parser = OptionParser(usage="usage: %prog ver [options -c, -e, -p, -m]")
parser.add_option("-c","--cut",   dest="cut", type="int", default=8, help="Plots after a certain cut")
parser.add_option("--mass", dest="mass", type="int", default=125,    help="Signal sample (mass)")
parser.add_option("-m","--merge", dest="merge",action="store_true", default=False, help="Do merging?")

parser.add_option("-p", "--period",dest="period", default="2012",  help="Year period; 2011 or 2012")
parser.add_option("--bkg",  dest="bkg",  action="store_true", default=False, help="Make plots from bkg sample")
parser.add_option("--qcd",  dest="qcd",  action="store_true", default=False, help="Include QCD samples")
parser.add_option("--mcfm", dest="mcfm", action="store_true", default=False, help="Use MCFM  as a signal")

parser.add_option("--evt",dest="evt",action="store_true", default=False, help="Show raw events (not scaled to lumi) for MC samles")
parser.add_option("-e","--extra",dest="extra",action="store_true", default=False, help="Make all extra plots")
parser.add_option("--fit",dest="fit",action="store_true", default=False, help="Do the various fits")
parser.add_option("--apz",dest="apz",action="store_true", default=False, help="Discover new particle (requires apzTree)")
parser.add_option("--zjp",dest="zjp",action="store_true", default=False, help="Study J/Psi region and more (requires apzTree)")
parser.add_option("--hjp",dest="hjp",action="store_true", default=False, help="Study J/Psi region and more (requires apzTree)")
parser.add_option("-s", '--sel', dest="sel", type="string", default='mugamma',
                  help="Selection to be used. Options are: '4mu','2e2mu', 'zee','mugamma', 'elgamma'")

(opt, args) = parser.parse_args()

mass = opt.mass
sel  = opt.sel

comments = ["These plots are made for h -> ll gamma",
#comments = ["These plots are made for z -> J/Psi gamma analysis",
            "MuEG dataset used"]
doLumiScale = 1
if opt.evt: doLumiScale=0

if __name__ == "__main__":
  timer = TStopwatch()
  timer.Start()

  if len(args) < 1:
    parser.print_usage()
    exit(1)

  ver    = sys.argv[1]
  if 'vv/' in ver: ver = ver[3:].rstrip('/')
  cut=str(opt.cut)
  doMerge = opt.merge
  period  = opt.period
  doBkg   = opt.bkg

  #gROOT.ProcessLine(".L ../tdrstyle.C")
  #setTDRStyle()
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
    if sel=="elgamma" or sel=='eegamma':
      os.system("hadd "+hPath+"/m_Data_" +sel+"_"+period+".root "+hPath+"/"+sel+"_"+period+"/hhhh_*Run2012*.root")
    else:
      # os.system("hadd "+hPath+"/m_Data_" +sel+"_"+period+".root "+hPath+"/"+sel+"_"+period+"/hhhh_DoubleMu_Run20*.root")
      os.system("hadd "+hPath+"/m_Data_" +sel+"_"+period+".root "+hPath+"/"+sel+"_"+period+"/hhhh_MuEG_Run2012*.root")
      # os.system("hadd "+hPath+"/m_Data_" +sel+"_"+period+".root "+hPath+"/"+sel+"_"+period+"/hhhh_DoubleElectron_Run20*.root")


    if doBkg:
      #os.system("hadd "+hPath+"/m_DYJetsPow20_"+sel+"_"+period+".root "
      #          +hPath+"/"+sel+"_"+period+"/hhhh_DYJetsPow20-RD1*.root")
      #os.system("hadd "+hPath+"/m_DYJets50_"+sel+"_"+period+".root "
      #          +hPath+"/"+sel+"_"+period+"/hhhh_DYJets50-RD1*.root")
      #os.system("hadd "+hPath+"/m_ZG_"+sel+"_"+period+".root "
      #          +hPath+"/"+sel+"_"+period+"/hhhh_ZGToLLG-RD1*.root ")

      if opt.qcd:
        os.system("hadd "+hPath+"/m_QCD_EM_Pt_20to30_"+sel+"_"+period+".root "
                  +hPath+"/"+sel+"_"+period+"/hhhh_QCD_EM_Pt_20to30_*.root ")
        os.system("hadd "+hPath+"/m_QCD_EM_Pt_30to80_"+sel+"_"+period+".root "
                  +hPath+"/"+sel+"_"+period+"/hhhh_QCD_EM_Pt_30to80_*.root ")
        os.system("hadd "+hPath+"/m_QCD_EM_Pt_80to170_"+sel+"_"+period+".root "
                  +hPath+"/"+sel+"_"+period+"/hhhh_QCD_EM_Pt_80to170_*.root ")
        os.system("hadd "+hPath+"/m_QCD_EM_Pt_170to250_"+sel+"_"+period+".root "
                  +hPath+"/"+sel+"_"+period+"/hhhh_QCD_EM_Pt_170to250_*.root ")
        os.system("hadd "+hPath+"/m_QCD_EM_Pt_250to350_"+sel+"_"+period+".root "
                  +hPath+"/"+sel+"_"+period+"/hhhh_QCD_EM_Pt_350_*.root ")
        os.system("hadd "+hPath+"/m_QCD_EM_Pt_250to350_"+sel+"_"+period+".root "
                  +hPath+"/"+sel+"_"+period+"/hhhh_QCD_EM_Pt_350_*.root ")

        # os.system("hadd "+hPath+"/m_QCD_Mu_Pt_20_"+sel+"_"+period+".root "
        #          +hPath+"/"+sel+"_"+period+"/hhhh_QCD_Mu_Pt_20_*.root ")
        # os.system("hadd "+hPath+"/m_QCD_Mu_Pt_15to30_"+sel+"_"+period+".root "
        #          +hPath+"/"+sel+"_"+period+"/hhhh_QCD_Mu_Pt_5to30_*.root ")
        os.system("hadd "+hPath+"/m_QCD_Mu_Pt_30to50_"+sel+"_"+period+".root "
                  +hPath+"/"+sel+"_"+period+"/hhhh_QCD_Mu_Pt_30to50_*.root ")
        os.system("hadd "+hPath+"/m_QCD_Mu_Pt_50to150_"+sel+"_"+period+".root "
                  +hPath+"/"+sel+"_"+period+"/hhhh_QCD_Mu_Pt_50to150_*.root ")
        os.system("hadd "+hPath+"/m_QCD_Mu_Pt_150_"+sel+"_"+period+".root "
                  +hPath+"/"+sel+"_"+period+"/hhhh_QCD_Mu_Pt_150_*.root ")

  qcdSamples = None
  if doBkg:
    #bkgFiles.append(TFile(hPath+"/m_ZG_"+sel+"_"+period+".root","OPEN"))
    #bkgNames.append('ZG')
    #bkgFiles.append(TFile(hPath+"/m_DYJets50_"+sel+"_"+period+".root","OPEN"))
    #bkgNames.append('DYJets50')
    bkgFiles.append(TFile(hPath+"/"+sel+"_"+period+"/hhhh_DYJetsDalitz_1.root","OPEN"))
    bkgNames.append('DYJetsDalitz')
    bkgFiles.append(TFile(hPath+"/"+sel+"_"+period+"/hhhh_ZGDalitz_1.root","OPEN"))
    bkgNames.append('ZGDalitz')
    if opt.qcd:
      qcdSamples = {'QCD_EM_Pt_20to30':  TFile(hPath+"/m_QCD_EM_Pt_20to30_"  +sel+"_"+period+".root","OPEN"),
                    'QCD_EM_Pt_30to80':  TFile(hPath+"/m_QCD_EM_Pt_30to80_"  +sel+"_"+period+".root","OPEN"),
                    'QCD_EM_Pt_80to170': TFile(hPath+"/m_QCD_EM_Pt_80to170_" +sel+"_"+period+".root","OPEN"),
                    'QCD_EM_Pt_170to250':TFile(hPath+"/m_QCD_EM_Pt_170to250_"+sel+"_"+period+".root","OPEN"),
                    'QCD_EM_Pt_250to350':TFile(hPath+"/m_QCD_EM_Pt_250to350_"+sel+"_"+period+".root","OPEN"),
                    'QCD_EM_Pt_250to350':TFile(hPath+"/m_QCD_EM_Pt_250to350_"+sel+"_"+period+".root","OPEN"),
                    'QCD_EM_Pt_350':     TFile(hPath+"/m_QCD_EM_Pt_350_"     +sel+"_"+period+".root","OPEN"),

                    'QCD_Mu_Pt_20':     TFile(hPath+"/m_QCD_Mu_Pt_20_"      +sel+"_"+period+".root","OPEN"),
                    'QCD_Mu_Pt_15to30': TFile(hPath+"/m_QCD_Mu_Pt_15to30_"  +sel+"_"+period+".root","OPEN"),
                    'QCD_Mu_Pt_30to50': TFile(hPath+"/m_QCD_Mu_Pt_30to50_"  +sel+"_"+period+".root","OPEN"),
                    'QCD_Mu_Pt_50to150':TFile(hPath+"/m_QCD_Mu_Pt_50to150_" +sel+"_"+period+".root","OPEN"),
                    'QCD_Mu_Pt_150':    TFile(hPath+"/m_QCD_Mu_Pt_150_"     +sel+"_"+period+".root","OPEN")
                  }

      bkgNames.append('QCD')
      bkgFiles.append(qcdSamples)

    yields_bkg = []
    for b in xrange(len(bkgNames)):
      yields_bkg.append(u.getYields(bkgFiles[b],bkgNames[b],True))

    bkgZip = zip(bkgNames, bkgFiles)
  else: bkgZip = None

  sigGG  = {}
  sigVBF = {}
  sigVH  = {}
  for m in ['120','125','130','135','140','145','150']:
    sigGG[m]  = TFile(hPath+"/"+sel+"_"+period+"/hhhh_ggH-mad"+m+"_1.root", "OPEN")
    sigVBF[m] = TFile(hPath+"/"+sel+"_"+period+"/hhhh_vbfH-mad"+m+"_1.root", "OPEN")
    sigVH[m]  = TFile(hPath+"/"+sel+"_"+period+"/hhhh_vH-mad"+ m+"_1.root", "OPEN")

  sigFileGG   = sigGG['125']
  sigFileVBF  = sigVBF['125']
  sigFileVH   = sigVH['125']

  sigFileZjp  = TFile(hPath+"/"+sel+"_"+period+"/hhhh_ZtoJPsiGamma_1.root",     "OPEN")
  sigFileHjp  = TFile(hPath+"/"+sel+"_"+period+"/hhhh_HiggsToJPsiGamma_1.root", "OPEN")

  if opt.mcfm:  sigFile = sigFileMCFM
  elif opt.hjp: sigFile = sigFileHjp
  elif opt.zjp: sigFile = sigFileZjp
  else:         sigFile = sigFileGG

  dataFile = TFile(hPath+"/m_Data_"+sel+"_"+period+".root","OPEN")

  if sel=="zee":
    u.drawAllInFile(dataFile, "data",bkgZip,None,"",
                    "Zee",  pathBase+"/Zee",  None,"norm")

  if opt.zjp or opt.hjp:
    u.setSelection('mugamma')

  yields_data = u.getYields(dataFile)
  yields_zjp  = u.getYields(sigFileZjp, 'ZtoJPsiGamma',  doLumiScale)
  yields_hjp  = u.getYields(sigFileHjp, 'HtoJPsiGamma',  doLumiScale)
  yields_ggH  = u.getYields(sigFileGG,  'ggH-125',  doLumiScale)
  yields_vbf  = u.getYields(sigFileVBF, 'vbfH-125', doLumiScale)
  yields_vh   = u.getYields(sigFileVH,  'vH-125',   doLumiScale)
  yields_sig  = [sum(x) for x in zip(yields_ggH,yields_vbf,yields_vh)]

  print 'ggH yi', yields_ggH
  print 'sig yi', yields_sig

  subdir = 'Main'
  path = pathBase+"/"+subdir
  #if doBkg:
  #  path = pathBase+"/bkg_"+subdir
  #  u.createDir(path)
  #path = pathBase+"/"+subdir


  sigName = '#splitline{XX x Signal}{m_{H}=125 GeV}'
  if opt.zjp: sigName= '#splitline{XX x Signal}{Z #rightarrow J/Psi #gamma}'
  if opt.hjp: sigName= 'h #rightarrow J/Psi #gamma'


  ggHZGFile   = TFile(hPath+"/"+sel+"_"+period+"/hhhh_ggHZG-"+str(mass)+"_1.root", "OPEN")
  hackZip = zip(['120','130','135','140','145','150'],
                [sigGG['120'],sigGG['130'],sigGG['135'], sigGG['140'],sigGG['145'],sigGG['150']])
  sigZip = zip(['125','135','145'],
                [sigGG['125'],sigGG['135'],sigGG['145']])
  #sigZip = zip(['125'],
  #             [sigGG['125']])

  #u.drawAllInFile(None, "", hackZip, None, '', "GEN", pathBase+'/GEN-lumi', None, "lumi")

  # u.drawAllInFile(None, "", hackZip, sigGG['125'], '125', "GEN", pathBase+'/GEN-norm', None, "norm", doFits=True)
  # u.drawAllInFile(None, "", '', sigFileGG, "h", "GEN", pathBase+'/GEN-norm', None, "norm")
  # u.drawAllInFile(None, "", '', ggHZGFile, "h", "GEN", pathBase+'/GEN-norm', None, "norm")

  #u.drawAllInFile(None, "", '', sigFileGG, "Dalitz", "eff", pathBase+'/eff', None, "norm")

  if not opt.apz:
    if opt.bkg:
      u.drawAllInFile(dataFile, "Data", bkgZip, None, sigName, "Main", path, cut, "lumi")
      # u.drawAllInFile(dataFile, "Data", bkgZip, sigFile, sigName, "Main", path, cut, "lumi")
    else:
      #u.drawAllInFile(dataFile, "Data", None, None, sigName, "Main", path, cut, "toDataInt", doFits=opt.fit)
      u.drawAllInFile(None, "Data", None, sigZip, sigName, "Main", path, cut, "toDataInt", doFits=opt.fit)
      #u.drawAllInFile(dataFile, "Data", None, sigZip, sigName, "Main", path, cut, "toDataInt", doFits=opt.fit)

    if opt.extra:
      # u.drawAllInFile(None, "", bkgZip, sigFile, sigName, "GEN", pathBase+'/GEN', None, "norm")
      for n in ['Angles','N']:
        u.drawAllInFile(dataFile, "Data",bkgZip,sigFile,"signal", n, pathBase+"/"+n, cut, "lumi")

      if cut in ['8','9']:
        for n in ['Photon','Photon-EGamma']:
          u.drawAllInFile(dataFile, "Data",bkgZip,sigFile,"signal", n, pathBase+"/"+n, None, "norm")

        if sel in ["mugamma","jp-mugamma"]:
          for n in ['Muons','Mu-after']:
            u.drawAllInFile(dataFile, "Data",bkgZip,sigFile,"signal", n, pathBase+"/"+n, None, "norm")

  elif cut in ['14','15','16']:
    u.drawAllInFile(dataFile, "data",bkgZip,None,"",
                    "AlphaPiZ",pathBase+"/apz/", cut,"norm2")

  if sel == "elgamma":
    u.drawAllInFile(dataFile, "Data", bkgZip, sigFile, sigName, "Main-Dale", pathBase+'/Dale/', cut, "toDataInt")
    #u.drawAllInFile(dataFile, "data",bkgZip,sigFile,"signal",
    #                "AnElectron-EGamma",  pathBase+"/AnElectron-EGamma/",  None,"norm")

    if cut in ['9','12','15'] and sel=='elgamma':
      u.drawAllInFile(dataFile, "Data",bkgZip,sigFile,"signal",
                      "AnElectron",  pathBase+"/AnElectron/",  None,"norm")
      dirnameshort = "DalitzEle"
      dirname = dirnameshort+"-cut"+cut
      u.drawAllInFile(dataFile, "Data",bkgZip,sigFile,"signal",
                      dirname,  pathBase+"/"+dirnameshort+"/",  None,"norm")

      for n in ['EGamma','BaseSC-1','BaseSC-2',"track-1","track-2","track-3"]:
        u.drawAllInFile(dataFile, "Data",bkgZip,sigFile,"signal",
                        dirname+"-"+n,  pathBase+"/"+dirnameshort+"-"+n,  None,"norm")

  ss = opt.sel
  if opt.apz and doMerge:
    if ss in ['4mu','2e2mu']:
      os.system("hadd -f "+hPath+"/m_Data_apz_DoubleMu_"+ss+"_"+period+".root "+hPath+"/"+ss+"_"+period+"/hhhh_DoubleMu_Run20*.root")

    if ss in ['2e2mu','mugamma']:
      os.system("hadd -f "+hPath+"/m_Data_apz_MuEG_"+ss+"_"+period+".root "+hPath+"/"+ss+"_"+period+"/hhhh_MuEG_Run20*.root")
    if ss =='2e2mu':
      os.system("hadd -f "+hPath+"/m_Data_apz_DoubleElectron_"+ss+"_"+period+".root "+hPath+"/"+ss+"_"+period+"/hhhh_DoubleElectron_Run20*.root")

  c1 = TCanvas("c4","small canvas",600,600);


  if opt.apz:
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

        alphaPiZ2(data, c1, TCut('m4l>110 && m4l<170'),pathBase+"/"+s+"-alphaPiZ-7/")


        u.drawAllInFile(data, "data",bkgZip,None,"",
                        "apz-plots",pathBase+"/apz-plots/",None,"norm2")


    elif ss=='mugamma':

      data = TFile(hPath+"/m_Data_mugamma_2012.root","OPEN")

      alphaPiZ2(data, c1, TCut('m123>110 && m123<170 && pt12/m123>0.3 && pt3/m123>0.3'),pathBase+"/"+ss+"-APZ-7/")

      alphaPiZ(data, c1, TCut('ph_pt/m_llg>0.3 && di_pt/m_llg>0.3 && m_llg>100&&m_llg<170'),
               pathBase+"/alphaPiZ-0/")
      alphaPiZ(data, c1, TCut(''),                                                            pathBase+"/alphaPiZ-1/")
      alphaPiZ(data, c1, TCut('ph_pt/m_llg>0.30 && di_pt/m_llg>0.30'),                        pathBase+"/alphaPiZ-2/")
      alphaPiZ(data, c1, TCut('ph_pt/m_llg>0.30 && di_pt/m_llg>0.30 && m_llg>100&&m_llg<150'),pathBase+"/alphaPiZ-3/")
      alphaPiZ(data, c1, TCut('ph_pt/m_llg>0.35 && di_pt/m_llg>0.35 && m_llg>100&&m_llg<150'),pathBase+"/alphaPiZ-4/")
      alphaPiZ(data, c1, TCut('ph_pt/m_llg>0.35 && di_pt/m_llg>0.35 && m_llg>120&&m_llg<180'),pathBase+"/alphaPiZ-5/")
      alphaPiZ(data, c1, TCut('ph_pt/m_llg>0.35 && di_pt/m_llg>0.35 && m_llg>100&&m_llg<130'),pathBase+"/alphaPiZ-6/")

      alphaPiZ(data, c1, TCut('ph_pt/m_llg>0.30 && di_pt/m_llg>0.30 && m_llg>121&&m_llg<131'),pathBase+"/alphaPiZ-7/")
      alphaPiZ(data, c1, TCut('ph_pt/m_llg>0.30 && di_pt/m_llg>0.30 && m_llg>85&&m_llg<96'),  pathBase+"/alphaPiZ-8/")

  '''
  htmp = sigFileGG.Get("02_lPt2_pt__cut"+cut)
  int1 = htmp.Integral()
  int2 = htmp.Integral(0,10)
  print '\n Fraction of event with trailing lepton pt<20: \n', float(int2)/int1

  htmp = sigFileGG.Get("04_ll_deltaR__cut"+cut)
  int1 = htmp.Integral()
  int2 = htmp.Integral(0,10)
  print '\n Fraction of event with dR(l1,l2)<0.4: \n', float(int2)/int1


  c1 = TCanvas("c4","small canvas",600,600);
  c1.cd()
  Nev = sigFileGG.Get("Counts/evt_byCut").GetBinContent(2)
  cro  = u.getCS("ggH-125",'mu')
  lumi = u.getLumi("2012")

  scale = float(lumi*cro)/Nev
  '''

  #u.AFB(sigFileZjp, 'signal', 'FB_cosSC-diG_vs_Mllg_cut'+cut, pathBase+'/Angles/')
  #u.AFB(dataFile,   'data',   'FB_cosSC-diG_vs_Mllg_cut'+cut, pathBase+'/Angles/')

  #u.FBAss(sigFileGG, 'signal', 'FB_cosSC-ll_vs_Mll_cut'+cut, pathBase+'/Angles/')
  #u.FBAss(dataFile,   'data',   'FB_cosSC-ll_vs_Mll_cut'+cut, pathBase+'/Angles/')


  '''
  hc7 = sigFileGG.Get("tri_mass80__cut7").Integral(35,41)*scale
  hc8 = sigFileGG.Get("tri_mass80__cut8").Integral(35,41)*scale
  hc9 = sigFileGG.Get("tri_mass80__cut9").Integral(35,41)*scale
  print "Yields in bins that supposed to correspond to [122,128] window:\n",hc7, hc8, hc9

  hc7 = dataFile.Get("tri_mass80__cut7").Integral(36,40)
  hc8 = dataFile.Get("tri_mass80__cut8").Integral(36,40)
  hc9 = dataFile.Get("tri_mass80__cut9").Integral(36,40)
  print "Yields in bins that supposed to correspond to [122,128] window:\n",hc7, hc8, hc9

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

  if opt.zjp and cut=='4':
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




  if sel in ['jp-mugamma']:
    if opt.hjp: u.setCutListFile(hPath+"/"+sel+"_"+period+"/out_cutlist_HiggsToJPsiGamma_1.txt")
    if opt.zjp: u.setCutListFile(hPath+"/"+sel+"_"+period+"/out_cutlist_ZtoJPsiGamma_1.txt")
  else:         u.setCutListFile(hPath+"/"+sel+"_"+period+"/out_cutlist_ggH-mad120_1.txt")

  plot_types =[]
  dirlist = os.listdir(pathBase)
  for d in dirlist:
    if os.path.isdir(pathBase+"/"+d):
      plot_types.append(d)


  print ' Signal yields: '
  print yields_sig

  if doBkg:
    names = ['Data','ZG','ZJet', 'Signal @125']
    table_all  = u.yieldsTable([yields_data,yields_bkg[0], yields_bkg[1],yields_sig], names)
    #table_all  = u.yieldsTable([yields_data,yields_bkg,yields_sig, yields_ggH,yields_vbf, yields_vh], sel)
  else:
    if opt.zjp or opt.hjp:
      names = ['Data','H to J/Psi &gamma;','Z to J/Psi &gamma;', 'ggH-125']
      table_all = u.yieldsTable([yields_data, yields_hjp, yields_zjp, yields_ggH], names)
    #elif not opt.apz:
    else:
      names = ['Data','Sig: total','ggH','vbfH','VH']
      table_all  = u.yieldsTable([yields_data,yields_sig, yields_ggH,yields_vbf, yields_vh], names)


  u.makeTable(table_all,"all", "html", precision='%.2f')
  u.makeTable(table_all,"all", "twiki")
  u.makeTable(table_all,"all", "tex")

  os.system("cat yields_all.html > yields.html")

  defaultPage = 'Main'

  ht.makeHTML("h &rarr; dalitz decay plots",pathBase, plot_types, comments, defaultPage)

  print "\n\t\t finita la comedia \n"
