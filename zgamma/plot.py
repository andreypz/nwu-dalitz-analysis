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
parser.add_option("--sig",  dest="sig",  action="store_true", default=False, help="Signal MC")
parser.add_option("--data", dest="data", action="store_true", default=False, help="Data only")
parser.add_option("--vbf",  dest="vbf",  action="store_true", default=False, help="Use signal samples: ggH, vbf, vH")

parser.add_option("--spec",dest="spec",action="store_true", default=False, help="Make some special plots at the end.")
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

  if '/tthome' in os.getcwd():
    pathBase = "/tthome/andrey/html/zgamma/dalitz/"+ver+"_cut"+cut
    hPath    = "/tthome/andrey/batch_output/zgamma/8TeV/"+ver
  else:
    pathBase = "/uscms_data/d2/andreypz/html/zgamma/dalitz/"+ver+"_cut"+cut
    hPath    = "/eos/uscms/store/user/andreypz/batch_output/zgamma/8TeV/"+ver


  u.createDir(pathBase)

  if doMerge:
    os.system("rm "+hPath+"/m_*.root") #removing the old merged files

  dataFile  = None
  bkgFiles  = []
  bkgNames  = []

  u.setSelection(sel)
  subsel=sel
  if opt.hjp or opt.zjp:
    subsel = 'jp-mugamma'

  if doMerge:
    if sel=="elgamma" or sel=='eegamma':
      os.system("hadd "+hPath+"/m_Data_" +sel+"_"+period+".root "+hPath+"/"+sel+"_"+period+"/hhhh_*Run2012*.root")
    else:
      # os.system("hadd "+hPath+"/m_Data_" +subsel+"_"+period+".root "+hPath+"/"+subsel+"_"+period+"/hhhh_DoubleMu_Run20*.root")
      os.system("hadd "+hPath+"/m_Data_" +subsel+"_"+period+".root "+hPath+"/"+subsel+"_"+period+"/hhhh_MuEG_Run2012*.root")
      # os.system("hadd "+hPath+"/m_Data_" +subsel+"_"+period+".root "+hPath+"/"+subsel+"_"+period+"/hhhh_DoubleElectron_Run20*.root")


    if doBkg:
      #os.system("hadd "+hPath+"/m_DYJetsPow20_"+subsel+"_"+period+".root "
      #          +hPath+"/"+subsel+"_"+period+"/hhhh_DYJetsPow20-RD1*.root")
      #os.system("hadd "+hPath+"/m_DYJets50_"+subsel+"_"+period+".root "
      #          +hPath+"/"+subsel+"_"+period+"/hhhh_DYJets50-RD1*.root")
      #os.system("hadd "+hPath+"/m_ZG_"+subsel+"_"+period+".root "
      #          +hPath+"/"+subsel+"_"+period+"/hhhh_ZGToLLG-RD1*.root ")

      if opt.qcd:
        os.system("hadd "+hPath+"/m_QCD_EM_Pt_20to30_"+subsel+"_"+period+".root "
                  +hPath+"/"+subsel+"_"+period+"/hhhh_QCD_EM_Pt_20to30_*.root ")
        os.system("hadd "+hPath+"/m_QCD_EM_Pt_30to80_"+subsel+"_"+period+".root "
                  +hPath+"/"+subsel+"_"+period+"/hhhh_QCD_EM_Pt_30to80_*.root ")
        os.system("hadd "+hPath+"/m_QCD_EM_Pt_80to170_"+subsel+"_"+period+".root "
                  +hPath+"/"+subsel+"_"+period+"/hhhh_QCD_EM_Pt_80to170_*.root ")
        os.system("hadd "+hPath+"/m_QCD_EM_Pt_170to250_"+subsel+"_"+period+".root "
                  +hPath+"/"+subsel+"_"+period+"/hhhh_QCD_EM_Pt_170to250_*.root ")
        os.system("hadd "+hPath+"/m_QCD_EM_Pt_250to350_"+subsel+"_"+period+".root "
                  +hPath+"/"+subsel+"_"+period+"/hhhh_QCD_EM_Pt_350_*.root ")
        os.system("hadd "+hPath+"/m_QCD_EM_Pt_250to350_"+subsel+"_"+period+".root "
                  +hPath+"/"+subsel+"_"+period+"/hhhh_QCD_EM_Pt_350_*.root ")

        # os.system("hadd "+hPath+"/m_QCD_Mu_Pt_20_"+subsel+"_"+period+".root "
        #          +hPath+"/"+subsel+"_"+period+"/hhhh_QCD_Mu_Pt_20_*.root ")
        # os.system("hadd "+hPath+"/m_QCD_Mu_Pt_15to30_"+subsel+"_"+period+".root "
        #          +hPath+"/"+subsel+"_"+period+"/hhhh_QCD_Mu_Pt_5to30_*.root ")
        os.system("hadd "+hPath+"/m_QCD_Mu_Pt_30to50_"+subsel+"_"+period+".root "
                  +hPath+"/"+subsel+"_"+period+"/hhhh_QCD_Mu_Pt_30to50_*.root ")
        os.system("hadd "+hPath+"/m_QCD_Mu_Pt_50to150_"+subsel+"_"+period+".root "
                  +hPath+"/"+subsel+"_"+period+"/hhhh_QCD_Mu_Pt_50to150_*.root ")
        os.system("hadd "+hPath+"/m_QCD_Mu_Pt_150_"+subsel+"_"+period+".root "
                  +hPath+"/"+subsel+"_"+period+"/hhhh_QCD_Mu_Pt_150_*.root ")

  qcdSamples = None
  if doBkg:
    #bkgFiles.append(TFile(hPath+"/m_ZG_"+subsel+"_"+period+".root","OPEN"))
    #bkgNames.append('ZG')
    #bkgFiles.append(TFile(hPath+"/m_DYJets50_"+subsel+"_"+period+".root","OPEN"))
    #bkgNames.append('DYJets50')
    bkgFiles.append(TFile(hPath+"/"+subsel+"_"+period+"/hhhh_DYJetsDalitz_1.root","OPEN"))
    bkgNames.append('DYJetsDalitz')
    bkgFiles.append(TFile(hPath+"/"+subsel+"_"+period+"/hhhh_ZGDalitz_1.root","OPEN"))
    bkgNames.append('ZGDalitz')
    if opt.qcd:
      qcdSamples = {'QCD_EM_Pt_20to30':  TFile(hPath+"/m_QCD_EM_Pt_20to30_"  +subsel+"_"+period+".root","OPEN"),
                    'QCD_EM_Pt_30to80':  TFile(hPath+"/m_QCD_EM_Pt_30to80_"  +subsel+"_"+period+".root","OPEN"),
                    'QCD_EM_Pt_80to170': TFile(hPath+"/m_QCD_EM_Pt_80to170_" +subsel+"_"+period+".root","OPEN"),
                    'QCD_EM_Pt_170to250':TFile(hPath+"/m_QCD_EM_Pt_170to250_"+subsel+"_"+period+".root","OPEN"),
                    'QCD_EM_Pt_250to350':TFile(hPath+"/m_QCD_EM_Pt_250to350_"+subsel+"_"+period+".root","OPEN"),
                    'QCD_EM_Pt_250to350':TFile(hPath+"/m_QCD_EM_Pt_250to350_"+subsel+"_"+period+".root","OPEN"),
                    'QCD_EM_Pt_350':     TFile(hPath+"/m_QCD_EM_Pt_350_"     +subsel+"_"+period+".root","OPEN"),

                    'QCD_Mu_Pt_20':     TFile(hPath+"/m_QCD_Mu_Pt_20_"      +subsel+"_"+period+".root","OPEN"),
                    'QCD_Mu_Pt_15to30': TFile(hPath+"/m_QCD_Mu_Pt_15to30_"  +subsel+"_"+period+".root","OPEN"),
                    'QCD_Mu_Pt_30to50': TFile(hPath+"/m_QCD_Mu_Pt_30to50_"  +subsel+"_"+period+".root","OPEN"),
                    'QCD_Mu_Pt_50to150':TFile(hPath+"/m_QCD_Mu_Pt_50to150_" +subsel+"_"+period+".root","OPEN"),
                    'QCD_Mu_Pt_150':    TFile(hPath+"/m_QCD_Mu_Pt_150_"     +subsel+"_"+period+".root","OPEN")
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
    sigGG[m]  = TFile(hPath+"/"+subsel+"_"+period+"/hhhh_ggH-mad"+m+"_1.root", "OPEN")
    sigVBF[m] = TFile(hPath+"/"+subsel+"_"+period+"/hhhh_vbfH-mad"+m+"_1.root", "OPEN")
    sigVH[m]  = TFile(hPath+"/"+subsel+"_"+period+"/hhhh_vH-mad"+ m+"_1.root", "OPEN")

  sigFileGG   = sigGG['125']
  sigFileVBF  = sigVBF['125']
  sigFileVH   = sigVH['125']

  sigFileZjp  = TFile(hPath+"/"+subsel+"_"+period+"/hhhh_ZtoJPsiGamma_1.root",     "OPEN")
  sigFileHjp  = TFile(hPath+"/"+subsel+"_"+period+"/hhhh_HiggsToJPsiGamma_1.root", "OPEN")

  if opt.mcfm:  sigFile = sigFileMCFM
  elif opt.hjp: sigFile = sigFileHjp
  elif opt.zjp: sigFile = sigFileZjp
  else:         sigFile = sigFileGG

  dataFile = TFile(hPath+"/m_Data_"+subsel+"_"+period+".root","OPEN")

  if sel=="zee":
    u.drawAllInFile(dataFile, "data",bkgZip,None,"",
                    "Zee",  pathBase+"/Zee",  None,"norm")

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
  if opt.hjp: sigName= 'H #rightarrow J/Psi #gamma'


  ggHZGFile   = TFile(hPath+"/"+sel+"_"+period+"/hhhh_ggHZG-"+str(mass)+"_1.root", "OPEN")
  hackZip = zip(['120','130','135','140','145','150'],
                [sigGG['120'],sigGG['130'],sigGG['135'], sigGG['140'],sigGG['145'],sigGG['150']])
  sigZip = zip(['ggH-125','ggH-135','ggH-145'],
                [sigGG['125'],sigGG['135'],sigGG['145']])

  if opt.vbf:
    sigZipVBF = zip(['ggH-125','vbfH-125','vH-125'],
                    [sigGG['125'],sigVBF['125'],sigVH['125']])
  #sigZip = zip(['125'],
  #             [sigGG['125']])
  if opt.hjp:
    sigZip = zip(['H to J/Psi #gamma'],[sigFileHjp])



  #u.drawAllInFile(None, "", hackZip, None, '', "GEN", pathBase+'/GEN-lumi', None, "lumi")
  # u.drawAllInFile(None, "", hackZip, sigGG['125'], '125', "GEN", pathBase+'/GEN-norm', None, "norm", doFits=True)
  # u.drawAllInFile(None, "", '', sigFileGG, "h", "GEN", pathBase+'/GEN-norm', None, "norm")
  # u.drawAllInFile(None, "", '', ggHZGFile, "h", "GEN", pathBase+'/GEN-norm', None, "norm")

  #u.drawAllInFile(None, "", '', sigFileGG, "Dalitz", "eff", pathBase+'/eff', None, "norm")

  if not opt.apz and not opt.fit:
    if opt.bkg:
      print
      u.drawAllInFile(dataFile, "Data", bkgZip, None, sigName, "Main", path, cut, "lumi")
      # u.drawAllInFile(dataFile, "Data", bkgZip, sigFile, sigName, "Main", path, cut, "lumi")
    else:
      print
      #u.drawAllInFile(dataFile, "Data", None, sigZip, sigName, "Main", path, cut)
      if opt.data:
        u.drawAllInFile(dataFile, "Data", None, None, '', "Main", pathBase+'/Main-Data', cut)
      if opt.sig:
        u.drawAllInFile(None, None, None, sigZip, sigName,"Main", pathBase+'/Main-Sig', cut)

      if opt.vbf:
        u.drawAllInFile(dataFile, "Data", None, sigZipVBF, sigName, "Main", pathBase+'/Main-VBF', cut, 'toData')
    if opt.extra:
      u.drawAllInFile(None, None, None, sigZip, sigName, "GEN", pathBase+'/GEN', None, "norm")
      for n in ['Angles','N']:
        u.drawAllInFile(dataFile, "Data", bkgZip, sigZip,"signal", n, pathBase+"/"+n, cut, "lumi")

      if cut in ['8','9']:
        for n in ['Photon','Photon-EGamma']:
          u.drawAllInFile(dataFile, "Data",bkgZip, sigZip,"signal", n, pathBase+"/"+n, None, "norm")

        if sel in ["mugamma"]:
          for n in ['Muons','Mu-after']:
            u.drawAllInFile(dataFile, "Data",bkgZip,sigZip,"signal", n, pathBase+"/"+n, None, "norm")

  elif cut in ['14','15','16']:
    u.drawAllInFile(dataFile, "data",bkgZip,None,"",
                    "AlphaPiZ",pathBase+"/apz/", cut,"norm2")

  if sel == "elgamma":
    u.drawAllInFile(dataFile, "Data", bkgZip, sigZip, sigName, "Main-Dale", pathBase+'/Dale/', cut)
    if opt.data:
      u.drawAllInFile(dataFile, "Data", None, None, '', "Main-Dale", pathBase+'/Dale-Data', cut)
    if opt.sig:
      u.drawAllInFile(None, None, None, sigZip, sigName,"Main-Dale", pathBase+'/Dale-Sig', cut)

    #u.drawAllInFile(dataFile, "data",bkgZip,sigFile,"signal",
    #                "AnElectron-EGamma",  pathBase+"/AnElectron-EGamma/",  None,"norm")

    if opt.extra:
      u.drawAllInFile(dataFile, "Data",bkgZip,sigZip,"signal",
                      "AnElectron",  pathBase+"/AnElectron/",  None,"norm")
      dirnameshort = "DalitzEle"
      dirname = dirnameshort+"-cut"+cut
      u.drawAllInFile(dataFile, "Data",bkgZip,sigZip,"signal",
                      dirname,  pathBase+"/"+dirnameshort+"/",  None,"norm")

      for n in ['EGamma','BaseSC-1','BaseSC-2',"track-1","track-2","track-3"]:
        u.drawAllInFile(dataFile, "Data",bkgZip,sigZip,"signal",
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


  if opt.hjp and opt.fit:
    u.createDir(pathBase+'/fits/')
    histoName = '01_diLep_mass_jpsi_%s_cut%i'%('Main', opt.cut)
    hDa = dataFile.Get('Main/'+histoName)
    hSi = sigFile.Get('Main/'+histoName)

    #Ntot = hDa.GetEntries()
    #print 'Total Data entries = ', Ntot
    #raw_input('Enter to continue')

    gSystem.Load("libRooFit")
    mll = RooRealVar("mll","mll",2.92,3.28)
    meanSi_mll  = RooRealVar("m",  "massSi", 3.1, 2.92, 3.25)
    meanDa_mll  = RooRealVar("m",  "massDa", 3.1, 2.92, 3.25)
    width_mll = RooConstVar("#Delta","width", 0.01)
    #width_mll = RooRealVar("#Delta","width", 0, 0, 1)
    sigmaSi_mll = RooRealVar("#sigma","sigmaSi", 0.05, 0.0, 0.2)
    sigmaDa_mll = RooRealVar("#sigma","sigmaDa", 0.05, 0.0, 0.2)

    #lamb = RooRealVar("lamb","lamb",-0.5, -20,1)
    #alpha = RooRealVar("#alpha","alpha",5, 0,100)
    #cbN   = RooRealVar("N","cbN",1, 0,100)
    p1 = RooRealVar("p1","p1", 0, -1, 1)
    #p2 = RooRealVar("p2","p2", 1, -5, 5)
    bkg = RooPolynomial("bkg","bkg", mll, RooArgList(p1));
    Nsig = RooRealVar("Nsig","signal fraction",    300, 50.,500);
    Nbkg = RooRealVar("Nbkg","background fraction",150, 10.,300);

    #frac = RooRealVar("frac","frac",0.7, 0.5, 1.0)

    mll.setRange('fitRange',2.9,3.3)
    #mll.setRange('fitRange',3.0,3.3)
    dhDaTmp = RooDataHist("dhDa","dhDa", RooArgList(mll), hDa)
    dhDa = dhDaTmp.reduce(RooFit.CutRange('fitRange'))

    voigtDa = RooVoigtian("voigtDa","voigtDa",mll,meanDa_mll,width_mll,sigmaDa_mll)
    #cbDa = RooCBShape('cbDa','cbDa',mll, mean_mll, sigma_mll, alpha, cbN)
    # expo   = RooExponential("expo", "Bkg Model", mll, lamb)
    SplusB = RooAddPdf("SplusB", "Bkg plus signal peak", RooArgList(voigtDa, bkg), RooArgList(Nsig,Nbkg))

    mllFrame = mll.frame()
    fitDa   = SplusB.fitTo(dhDa, RooFit.Range('fitRange'), RooFit.Save())#, RooFit.SumW2Error(kTRUE))
    #fitDa   = voigtDa.fitTo(dhDa, RooFit.Range('fitRange'), RooFit.Save())
    #fitDa   = cbDa.fitTo(dhDa, RooFit.Range('fitRange'), RooFit.Save())
    #lDa_0 = RooArgSet(mll)
    #parsDa = SplusB.getParameters(RooArgSet(mll))
    #pars = voigtDa.getParameters(RooArgSet(mll))
    #pars = cbDa.getParameters(RooArgSet(mll))

    dhDa.plotOn(mllFrame, RooFit.Name('data'))
    # SplusB.plotOn(mllFrame, RooFit.LineColor(kGreen+3))
    # SplusB.paramOn(mllFrame, RooFit.Label("Data"), RooFit.Layout(0.095,0.45,0.65))
    SplusB.plotOn(mllFrame, RooFit.LineColor(kGreen+3))
    SplusB.paramOn(mllFrame, RooFit.Label("Data"),  RooFit.Layout(0.68, 0.95,0.9),
                   RooFit.Format('N',RooFit.FixedPrecision(4)))

    #voigtDa.plotOn(mllFrame, RooFit.LineColor(kGreen+3))
    #voigtDa.paramOn(mllFrame, RooFit.Label("Data"),  RooFit.Layout(0.65, 1.0,0.95))
    #cbDa.plotOn(mllFrame, RooFit.LineColor(kGreen+3))
    #cbDa.paramOn(mllFrame, RooFit.Label("Data"), RooFit.Layout(0.095,0.45,0.65))

    #lDa_1, l_2 = RooArgList(mll), RooArgList(pars)
    #ff1 = voigtDa.asTF(RooArgList(mll), RooArgList(pars))

    mllFrame.SetTitle(';m_{#mu#mu} (GeV);Events')
    mllFrame.SetMinimum(0.01)
    mllFrame.Draw()
    c1.SaveAs(pathBase+'/fits/JPpeakFit_cut'+cut+'_data.png')

    mllFrame = mll.frame()
    dhSiTmp = RooDataHist("dhSi","dhSi", RooArgList(mll), hSi)
    dhSi = dhSiTmp.reduce(RooFit.CutRange('fitRange'))
    voigtSi = RooVoigtian("voigtSi","voigtSi",mll,meanSi_mll,width_mll,sigmaSi_mll)
    fitSi   = RooFitResult(voigtSi.fitTo(dhSi, RooFit.Range('fitRange'),  RooFit.Save()))
    #cbSi  = RooCBShape('cbSi','cbSi',mll, mean_mll, sigma_mll, alpha, cbN)
    #fitSi = RooFitResult(cbSi.fitTo(dhSi, RooFit.Range('fitRange'),  RooFit.Save()))

    #lSi_0 = RooArgSet(mll)

    #parsSi = voigtSi.getParameters(RooArgSet(mll))
    #pars = cbSi.getParameters(RooArgSet(mll))

    dhSi.plotOn(mllFrame, RooFit.Name('sig'), RooFit.MarkerStyle(23), RooFit.MarkerColor(kBlue), RooFit.LineColor(kBlue))
    voigtSi.plotOn(mllFrame)
    voigtSi.paramOn(mllFrame, RooFit.Label("Signal MC"), RooFit.Layout(0.68, 0.95,0.90),
                    RooFit.Format('N',RooFit.FixedPrecision(4)))
    #cbSi.plotOn(mllFrame)
    #cbSi.paramOn(mllFrame, RooFit.Label("Signal MC"), RooFit.Layout(0.65, 1.0,0.95))
    mllFrame.SetTitle(';m_{#mu#mu} (GeV);Events')
    mllFrame.SetMinimum(0.01)
    mllFrame.Draw()
    c1.SaveAs(pathBase+'/fits/JPpeakFit_cut'+cut+'_sig.png')
    #mllFrame.Print('all')
    #print     mllFrame.FindObject('paramBox')
    #pave = TPaveText(mllFrame.FindObject('voigtSi_paramBox'))
    #pave.SetFillColor(kBlue)

    # mllFrame.SetMaximum(6)
    #leg.Clear()
    #leg.AddEntry(mllFrame.findObject('data'),'Data','e1p')
    #leg.AddEntry(mllFrame.findObject('sig'),'h#rightarrow J/#Psi#gamma#rightarrow#mu#mu#gamma','pl')

    #l3_1, l3_2 = RooArgList(mll), RooArgList(pars)
    #ffSi = voigt3.asTF(RooArgList(mll), RooArgList(pars))

    #ff1.Print('all')
    #ff3.Print('all')
    #NBINS = 100
    #hf1 = TH1F("hf1",";mll;hf1", NBINS, 2.9, 3.3)
    #hf3 = TH1F("hf3",";mll;hf3", NBINS, 2.9, 3.3)

    """
    for b in xrange(NBINS):
      m = 2.9+0.4*(b+0.5)/NBINS
      hf1.SetBinContent(b+1,ff1.Eval(m))
      hf3.SetBinContent(b+1,ff3.Eval(m))
              #print b, ff1.Eval(m), ff3.Eval(m)
      integ1 = hf1.Integral()
      hf1.Scale(1./integ1)
      integ3 = hf3.Integral()
      hf3.Scale(1./integ3)

      hratio = hf1.Clone()
      hratio.Divide(hf3)
      hratio.SetName('jpsi')
      hjpsi = TFile('JPsiMuMu_MCReScale.root','recreate')
      hjpsi.cd()
      hratio.Write()
      hjpsi.Close()

      # hratio.Print('all')
      # dhratio = RooDataHist("dhratio","dhratio", RooArgList(mll), hratio)
      # raw_input('Hit Enter to continue')

      # dhratio.plotOn(mllFrame, RooFit.Name('hratio'), RooFit.MarkerColor(kOrange))

      mllFrame.SetTitle(';m_{#mu#mu} (GeV);Events')
      mllFrame.Draw()
      mllFrame.SetMinimum(0.01)
      # mllFrame.SetMaximum(6)
      leg.Clear()
      leg.AddEntry(mllFrame.findObject('data'),'Data','e1p')
      leg.AddEntry(mllFrame.findObject('sig'),'h#rightarrow J/#Psi#gamma#rightarrow#mu#mu#gamma','pl')

    """


  if opt.spec:

    doLog = 1
    doGen = 1
    sigFactor = 1
    prodSF = 1.13

    for j,h in enumerate(['01_diLep_mass_0to20_b40','01_diLep_mass_0to20_b50','04_ll_gamma_deltaR','02_lPt1_pt','02_lPt2_pt','03_gamma_pt']):
      #   for h in ['01_diLep_mass_0to20_b50','04_ll_gamma_deltaR']:
      #vMu = 'v98-mu-incl-JPsi/'
      #vEl = 'v99-el'
      vMu = 'v03-mll-mu-paper'
      vEl = 'v03-mll-el-paper'
      print
      u.createDir(pathBase+'/Spec')
      sigFileMu = TFile('/tthome/andrey/batch_output/zgamma/8TeV/'+vMu+'/mugamma_2012/hhhh_ggH-mad125_1.root', 'OPEN')
      sigFileEl = TFile('/tthome/andrey/batch_output/zgamma/8TeV/'+vEl+'/elgamma_2012/hhhh_ggH-mad125_1.root', 'OPEN')
      hMu = sigFileMu.Get('Main/'+h+'_Main_cut9')
      hEl = sigFileEl.Get('Main/'+h+'_Main_cut15')


      if doGen and j==0:
        vMuG = 'v03-mll-mu-gen-paper'
        vElG = 'v03-mll-el-gen-paper'
        sigFileMuGen = TFile('/tthome/andrey/batch_output/zgamma/8TeV/'+vMuG+'/mugamma_2012/hhhh_ggH-mad125_1.root', 'OPEN')
        sigFileElGen = TFile('/tthome/andrey/batch_output/zgamma/8TeV/'+vElG+'/elgamma_2012/hhhh_ggH-mad125_1.root', 'OPEN')
        ghMu = sigFileMuGen.Get('GEN/gen_Mll_init_b40')
        ghEl = sigFileElGen.Get('GEN/gen_Mll_init_b40')


      lumi= u.getLumi()
      Nev = u.getTotalEvents(sigFileMu)
      cro = u.getCS('ggH-125', 'mu')
      scale = float(lumi*cro)/Nev
      hMu.Scale(sigFactor*scale)
      if j==0 and doGen: ghMu.Scale(sigFactor*scale*prodSF)
      print Nev, lumi, cro, scale

      Nev = u.getTotalEvents(sigFileEl)
      cro = u.getCS('ggH-125', 'el')
      scale = float(lumi*cro)/Nev
      hEl.Scale(sigFactor*scale)
      if j==0 and doGen: ghEl.Scale(sigFactor*scale*prodSF*(18.5/0.82))
      print Nev, lumi, cro, scale


      if j==0 and doGen:
        print 'mu selection:', hMu.Integral()
        print 'el selection:', hEl.Integral()
        print 'mu gen:', ghMu.Integral()
        print 'el gen:', ghEl.Integral()

      #if j==1:
      #  factor = hMu.Integral()/hEl.Integral()
      #  hEl.Scale(factor)

      hEl.Draw('hist')
      hMu.Draw('hist same')
      if j==0 and doGen:
        ghEl.Draw('hist same')
        ghMu.SetLineColor(kMagenta)
        ghEl.SetLineColor(kBlue-9)
        ghMu.SetLineStyle(7)
        #ghEl.SetLineStyle(4)
        ghEl.SetFillColor(19)
        #ghEl.SetFillColor(kCyan-10)

        #redraw them
        hEl.Draw('hist same')
        hMu.Draw('hist same')
        ghMu.Draw('hist same')
        hMu.SetLineWidth(3)
        ghMu.SetLineWidth(3)
        # hEl.SetFillColor(kRed-10)
        # hMu.SetLineStyle(2)
        hEl.SetFillColor(kOrange-9)

      c1.RedrawAxis()
      hEl.GetYaxis().SetTicks('-')
      hMu.SetLineColor(kGreen+3)
      hEl.SetLineColor(kRed-3)

      if j==0:
        if doLog: hEl.SetMaximum(100)
        else: hEl.SetMaximum(1.2*hEl.Clone().GetMaximum())
        #hEl.SetTitle(';m_{ll} (GeV); Events/0.4 GeV')
        hEl.SetTitle(';m_{\\ell\\ell}\\,\\mbox{(GeV)}; Events/0.5 GeV')
        hEl.GetXaxis().SetNdivisions(405, 0)
      elif j==1:
        hEl.SetMaximum(1.2*hMu.Clone().GetMaximum())
        hEl.SetAxisRange(1,5,"X")
      elif j==0:
        hEl.SetMaximum(1.2*hMu.Clone().GetMaximum())
        hEl.SetTitle(';#Delta R(#gamma,ll); Events')
      else:
        hEl.SetMaximum(1.2*hMu.Clone().GetMaximum())

      CMS_lumi(c1, 2, 11, "Simulation")

      lat = TLatex()
      lat.SetNDC()
      lat.SetTextSize(0.035)

      if j==0 and doGen:
        leg = TLegend(0.55,0.54,0.8,0.75);
        #leg = TLegend(0.55,0.54,0.8,0.88);
      else:
        leg = TLegend(0.61,0.62,0.78,0.85);
      leg.SetFillColor(kWhite)
      leg.SetBorderSize(0)
      leg.SetTextSize(0.035)
      #if j==0:
      #  leg.SetHeader('#splitline{10 x SM Higgs,}{  m_{H} = 125 GeV}')
      #elif j==1:
      #  leg.SetHeader('#splitline{SM Higgs,}{m_{H} = 125 GeV}')
      if sigFactor==1:
        lat.DrawLatex(0.57,0.82, '\\mbox{SM  H} \\rightarrow \\gamma^*\\gamma \\rightarrow \\ell\\ell\\gamma')
        lat.DrawLatex(0.62,0.77, 'm_{H} = 125 GeV')
        #leg.SetHeader('\\mbox{SM H} \\rightarrow \\gamma^*\\gamma \\rightarrow \\ell\\ell\\gamma \\atop \\mathrm{m_{H}} = 125\\, \\mbox{GeV}')
        #leg.SetHeader('#splitline{SM H #rightarrow #gamma*#gamma #rightarrow ll#gamma}{  m_{H} = 125 GeV}')
      else:
        leg.SetHeader('#splitline{'+str(sigFactor)+' x SM H #rightarrow #gamma*#gamma #rightarrow ll#gamma}{  m_{H} = 125 GeV}')
      if j==0 and doGen:
        leg.AddEntry(ghMu,' #mu#mu before selection', 'l')
        leg.AddEntry(ghEl,' ee before selection', 'f')
        leg.AddEntry(hMu,' #mu#mu after selection', 'l')
        leg.AddEntry(hEl,' ee after selection', 'f')
      else:
        leg.AddEntry(hMu,' #mu channel', 'l')
        leg.AddEntry(hEl,' e channel', 'l')

      leg.SetTextFont(42)
      leg.Draw()
      gPad.RedrawAxis()

      for e in ['pdf','png','eps']:
        if doLog: c1.SetLogy()
        c1.SaveAs(pathBase+'/Spec/'+h+'_sig.'+e)

      if j==0:
        c1.SetLogy(0)
        r1 = hMu.Clone()
        r1.Divide(ghMu)
        r1.Draw()
        r1.SetMaximum(0.5)
        r1.SetMinimum(0)
        r1.SetTitle(';m_{ll} (GeV);ratio')
        c1.SaveAs(pathBase+'/Spec/ratio_mu_'+h+'_sig.png')

        r2 = hEl.Clone()
        r2.Divide(ghEl)
        r2.Draw()
        r2.SetMaximum(0.5)
        r2.SetMinimum(0)
        r2.SetTitle(';m_{ll} (GeV);ratio')
        c1.SaveAs(pathBase+'/Spec/ratio_el_'+h+'_sig.png')


      '''
      daFileMu = TFile('/tthome/andrey/batch_output/zgamma/8TeV/'+vMu+'/m_Data_mugamma_2012.root', 'OPEN')
      daFileEl = TFile('/tthome/andrey/batch_output/zgamma/8TeV/'+vEl+'/m_Data_elgamma_2012.root', 'OPEN')
      dhMu = daFileMu.Get('Main/'+h+'_Main_cut9')
      dhEl = daFileEl.Get('Main/'+h+'_Main_cut15')

      dhEl.Draw('e1p')
      dhMu.Draw('same')
      if j==0:
        if doLog: dhEl.SetMaximum(1000)
        else: dhEl.SetMaximum(1.2*dhEl.Clone().GetMaximum())
      else:
        dhEl.SetMaximum(1.2*dhMu.Clone().GetMaximum())
      c1.RedrawAxis()
      dhMu.SetMarkerStyle(21)
      dhEl.SetMarkerStyle(22)
      dhMu.SetMarkerSize(1.2)
      dhEl.SetMarkerSize(1.2)

      dhMu.SetMarkerColor(kBlue+2)
      dhEl.SetMarkerColor(kRed-4)
      dhMu.SetLineColor(kBlue+2)
      dhEl.SetLineColor(kRed-4)
      if j==0:
        dhEl.SetTitle(';m_{ll} (GeV); Events/0.4 GeV')
      elif j==1:
        dhEl.SetAxisRange(1,5,"X")
        dhEl.SetTitle(';#Delta R(#gamma,ll); Events')
        #dhEl.Rebin()
        #dhMu.Rebin()

      leg = TLegend(0.61,0.62,0.78,0.85);
      leg.SetFillColor(kWhite)
      leg.SetBorderSize(0)
      leg.SetTextSize(0.04)
      leg.SetHeader('Data')
      leg.AddEntry(dhMu,' #mu channel', 'lp')
      leg.AddEntry(dhEl,' e channel', 'lp')
      leg.Draw()
      #lat.DrawLatex(0.18,0.95, 'H #rightarrow #gamma*#gamma #rightarrow ll#gamma')
      #CMS_lumi(c1, 2, 11)

      for e in ['pdf','png']:
        if doLog: c1.SetLogy()
        c1.SaveAs(pathBase+'/Spec/'+h+'_data.'+e)

      '''
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




  if sel in ['mugamma']:
    if opt.hjp: u.setCutListFile(hPath+"/"+subsel+"_"+period+"/out_cutlist_HiggsToJPsiGamma_1.txt")
    if opt.zjp: u.setCutListFile(hPath+"/"+subsel+"_"+period+"/out_cutlist_ZtoJPsiGamma_1.txt")
    if opt.spec: u.setCutListFile("localcutlist.txt")
  else:       u.setCutListFile(hPath+"/"+subsel+"_"+period+"/out_cutlist_ggH-mad120_1.txt")


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

  if opt.hjp: precision='%.3f'
  else:       precision='%.2f'

  u.makeTable(table_all,"all", "html", precision)
  u.makeTable(table_all,"all", "html", precision)
  u.makeTable(table_all,"all", "twiki")
  u.makeTable(table_all,"all", "tex")

  os.system("cat yields_all.html > yields.html")

  defaultPage = 'Main'

  ht.makeHTML("h &rarr; dalitz decay plots",pathBase, plot_types, comments, defaultPage)

  print "\n\t\t finita la comedia \n"
