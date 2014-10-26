#! /usr/bin/env python
from optparse import OptionParser
import sys,os,datetime,re
from array import *
from ROOT import *
gROOT.SetBatch()
gROOT.ProcessLine(".L ../tdrstyle.C")
gROOT.LoadMacro("../CMS_lumi.C")
setTDRStyle()
gROOT.ForceStyle()
TH1.SetDefaultSumw2(kTRUE)

import ConfigParser as cp
conf = cp.ConfigParser()
conf.optionxform = str
conf.read('../zgamma/config.cfg')
lumi2012 = float(conf.get("lumi","lumi2012A")) + float(conf.get("lumi","lumi2012B"))+\
    float(conf.get("lumi","lumi2012C")) + float(conf.get("lumi","lumi2012D"))
lumi = lumi2012

roundTo = 50.
# The rooFit stuff is for resolution fits and j/psi mass fits
gSystem.Load("libRooFit")
mll = RooRealVar("mll","mll",2.5,3.7)
mean_mll  = RooRealVar("m",     "mean", 3.0, 2.0, 4.0)
width_mll = RooRealVar("#Delta","width",0.0, 0.0, .1)
sigma_mll = RooRealVar("#sigma","sigma",1.0, 0.0, 2.0)
lamb = RooRealVar("lamb","lamb",-0.5, -20,1)
#res = RooRealVar("res","res",-1,1)
#mean_res  = RooRealVar("mean_res", "mean_res", 3.0, 0.0, 20.0)
#sigma_res = RooRealVar("sigma_res","sigma_res",1.0, 0.0, 10.0)
noOverflowList  = [a.strip() for a in (conf.get("selection","noOverflowList")).split(',')]
print 'No-Overflow list = ', noOverflowList

def mllBins():
  """First number is the mLL cut threshold, second number is the fraction of the cross section"""
  return [(0.200,0), (0.5, 0.145), (1, 0.138),(2, 0.138),(4, 0.134), (9, 0.157), (20, 0.149), (50,0.139)]

def yearToTeV(y):
  if   y=='2011': return '7TeV'
  elif y=='2012': return '8TeV'
  else: return '0TeV'

def setSelection(sel):
  conf.set("selection","sel", sel)
def getSelection():
  return conf.get("selection","sel")

def setCutListFile(fname):
  conf.set("selection","cutlist", fname)

def getCuts():
  fcuts = conf.get("selection","cutlist")
  with open(fcuts) as f:
    lines = f.read().splitlines()

  cuts = ['']*23

  for l in lines:
    n = int(l[0:2])
    c = l[2:].lstrip()
    #print n, c
    cuts[n] = c

  return cuts
# print key, cut

def getLumi(period='2012'):
  if period=="2012":
    return lumi2012
  else:
    print "Sorry only 2012 is considered"
    return 0

def getCS(sample, mySel=None):
  if mySel==None:
    #sel = conf.get("selection", "sel")[0:2]
    cs = float(eval(conf.get(sample, "cs")))
  else:
    sel = mySel
    cs = float(conf.get(sample, "cs-"+sel))

  return cs

def getColors(sample):
  l = 45
  f = 0
  if conf.has_option(sample,'line'):
    l = eval(conf.get(sample, "line"))
  if conf.has_option(sample,'fill'):
    f = eval(conf.get(sample, "fill"))

  return [l,f]

def getTotalEvents(f):
  ev = f.Get("Counts/evt_byCut")
  Nev = ev.GetBinContent(1)
  return Nev


def drange(start, stop, step):
  r = start
  while r <= stop:
    yield r
    r += step


def handleOverflowBins(hist):
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
  if hist.GetName()!="diLep_mass_low":
    hist.SetBinContent(nBins, lastBin+ovflBin);
  hist.SetBinError(nBins, sqrt(pow(lastBinErr,2) + pow(ovflBinErr,2)) );
  hist.SetBinContent(1, firstBin+undflBin);
  hist.SetBinError(1, sqrt(pow(firstBinErr,2) + pow(undflBinErr,2)) );
  hist.SetBinContent(0,0);
  hist.SetBinContent(nBins+1,0);

def draw(h1, path):
  h = h1
  if h.InheritsFrom("TH2"):
    h.Draw("lego")
  else:
    h.Draw("hist")

  c1.SaveAs(path+h.GetName()+".png")

def blindIt(h, nbins=10):
  hbin =  h.FindBin(125)
  for b in xrange(hbin-nbins,hbin+nbins):
    h.SetBinContent(b,0)
  h.SetMinimum(0.001)

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


def createDir(myDir):
  if not os.path.exists(myDir):
    try: os.makedirs(myDir)
    except OSError:
      if os.path.isdir(myDir): pass
      else: raise


def stackQCD(dic, hName, lumi):
  qcd_tot = TH1F()
  isInit = 0
  for key, f in dic.iteritems():
    h1 = f.Get(hName)
    if h1==None: continue
    else: h = h1.Clone()
    #print key, f

    Nev = getTotalEvents(f)
    cro = getCS(key)
    scale = float(lumi*cro)/Nev

    h.Scale(scale)

    #print key, hName, Nev, lumi, cro, scale

    if not isInit:
      qcd_tot = h # first non-epty hist
      isInit  = 1
    else:
      qcd_tot.Add(h)

  return qcd_tot

def makeStack(bZip, histDir, histoName, leg, lumi, howToScale, normToScale=None):
  hs = THStack("temp", "Stacked histo")

  #samplesToAdd = []
  if histDir!="":
    hName = histDir+"/"+histoName
  else:
    hName = histoName

  for n,f in bZip:
    #print n,f
    if n=='QCD':
      #raw_input('QCD is here. Hit enter to continue')
      h1 = stackQCD(f, hName, lumi)
    else:
      h1 = f.Get(hName)

    if h1==None:
      print 'None histogrammm:',hName
      continue
    else:
      h = h1.Clone()

    scale = 1

    h.SetLineColor(int(getColors(n)[0]))
    h.SetLineWidth(2)

    normh = h.Integral()
    #print 'norm =', normh, n
    if howToScale == 'lumi':
      if n!='QCD':  # The QCD hist is already scaled to lumi
        Nev = getTotalEvents(f)
        cro = getCS(n)
        scale = float(lumi*cro)/Nev
        print n, Nev, lumi, cro, scale

        h.Scale(scale)
        # only fill the colors if we are going to stack them (scale to lumi)

      h.SetFillColor( int(getColors(n)[1]))
    elif howToScale == 'norm':
      if normh!=0: h.Scale(1./normh)
    elif howToScale in ['toData','toDataInt'] and normToScale!=None:
      if normh!=0: h.Scale(normToScale/normh)
    else:
      print 'Sorry, there must be a mistake, this norm is not supported: ',howToScale
      sys.exit(0)

    hs.Add(h)
    if conf.has_option(n,"shortName"):
      leg.AddEntry(h, conf.get(n,"shortName") ,"f")
    else:
      leg.AddEntry(h, n, "l")

    #hm.Print('all')

  return hs


def drawAllInFile(f1, name1, bZip, f3, name3, myDir, path, N, howToScale="none", isLog=False, doRatio=False, doFits=False):
  print 'myDir is ', myDir
  if f1!=None and not f1.IsZombie():
    f1.cd(myDir)
  elif bZip!=None:
    print bZip[0]
    bZip[0][1].cd(myDir)
  elif f3!=None:
    f3.cd(myDir)
  else:
    print "Sorry can't draw anything, no files are provided!"
    sys.exit(0)

  dirList = gDirectory.GetListOfKeys()
  # dirList.Print()
  createDir(path)
  split = os.path.split(path.rstrip("/"))
  print "Split the path:", split[0], split[1]

  if doRatio:
    c1 = TCanvas("c2","big canvas",600,700);
  else:
    c1 = TCanvas("c3","small canvas",600,600);
  c1.SetLogy(isLog)
  c1.cd()
  #c1.UseCurrentStyle()

  scale3 = 1
  if f3!=None and howToScale in ['lumi','toData','toDataInt']:
    Nev = getTotalEvents(f3)
    if conf.get("selection","jp")=='1':
      cro = getCS("HtoJPsiGamma")
    else:
      cro = getCS("ggH-"+conf.get("selection","mass"), getSelection())
    scale3 = float(lumi*cro)/Nev
    print Nev, lumi, cro, scale3

  doPdf=0
  histoName = None
  for k in dirList:
    #print k
    if k.GetName() in ["eff"]: continue
    if N!=None:
      if not ("cut"+N) in k.GetName(): continue
      if "DalitzEle"   in k.GetName(): continue
      #if k.ReadObj().isFolder(): continue

    mainHist = k.ReadObj()

    histoName = mainHist.GetName()

    doOverflow = 1
    # Do overflow or do not overflow, that's the question!
    if any(x in histoName for x in noOverflowList):
      doOverflow = 0
    # print 'Do Overflow??', doOverflow

    h1 = TH1F()
    h3 = TH1F()

    hmaxs = []

    if f1!=None:
      if myDir!="":
        h1 = f1.Get(myDir+"/"+histoName)
      else:
        h1 = f1.Get(histoName)
      if h1==None: continue
    else: h1=None

    if f3!=None:
      if myDir!="":
        h3 = f3.Get(myDir+"/"+histoName)
      else:
        h3 = f3.Get(histoName)
        print myDir, histoName
      if h3==None: continue

      h3.Scale(float(scale3))

    print "\t Drawing", histoName

    if mainHist.InheritsFrom("TH2"):
      createDir(split[0]+"/TH2_"+split[1])
      if f1!=None and h1!=None:
        h1.Draw("col")
        prename = split[0]+"/TH2_"+split[1]+"/"+histoName
        if "_dalitzPlot_" in histoName:
          h1.SetNdivisions(505,'X')
          h1.SetNdivisions(505,'Y')

        CMS_lumi(c1, 2, 11)
        c1.SaveAs(prename+"_data.png")


      if f3!=None and h3!=None:
        h3.Draw("col")
        if "_dalitzPlot_" in histoName:
          h3.SetNdivisions(505,'X')
          h3.SetNdivisions(505,'Y')
        CMS_lumi(c1, 2, 11)
        c1.SaveAs(prename+"_sig.png")

      continue

    # set_palette("gray")
    # h1.Draw("col same")
    # if "h2D_tri_vs_diLep_mass" in histoName:
    #    c1.SetLogy()

    else:
      pad1 = TPad("pad1","pad1",0,0.3,1,1);
      pad2 = TPad("pad2","pad2",0,0,1,0.3);

      norm1=0
      if doRatio:
        pad1.SetBottomMargin(0);
        pad1.Draw();
        pad1.cd();
        pad1.SetLogy(isLog)


      if doOverflow:
        handleOverflowBins(h1)

      if f1!=None and h1!=None:
        #h1.Draw("hist")
        h1.Draw("same e1p")
        # h1.SetMinimum()
        h1.SetMarkerStyle(20)
        h1.SetMarkerSize(0.75)
        h1.SetLineColor(kBlack)
        # h1.UseCurrentStyle()
        #h1.SetMarkerStyle(20)
        #h1.SetMarkerSize(0.75)

        # leg.AddEntry(h1,name1, "l")
        norm1 = h1.Integral()
        if howToScale=='norm' and  norm1!=0:
          h1.Scale(1./norm1)
        hmaxs.append(h1.GetMaximum())

        mainHist = h1

        if "phi" in histoName:
          if isLog:
            h1.SetMinimum(1e-6)
          else:
            h1.SetMinimum(0)
        if "ll_deltaR_" in histoName:
          if not 'el' in getSelection():
            h1.SetXTitle("#DeltaR(#mu_{1}, #mu_{2})")
          doPdf=1
        if "ptDaleOverGamma" in histoName:
          doPdf=1
        if "dale_pt" in histoName:
          doPdf=1
        if "gamma_pt_" in histoName:
          h1.SetXTitle("Photon p_{T} (GeV)")
          doPdf=1
        if "lPt1_pt_" in histoName:
          if not 'el' in getSelection():
            h1.SetXTitle("Leading muon p_{T} (GeV)")
          doPdf=1
        if "gamma_deltaR" in histoName:
          h1.SetAxisRange(1,5,"X")
          doPdf=1
        if "deltaR_dale_gamma" in histoName:
          h1.SetAxisRange(1,5,"X")
          doPdf=1

      #if "tri_mass" in k1.GetName() and name1 not in ["madgra","mcfm"]:
      #
      #if "h_mass" in k1.GetName():
      #  print "\n *** H-mass RMS:",  h1.GetRMS(), h3.GetRMS()
      #  print "\n *** H-mass Mean:", h1.GetMean(),h3.GetMean()

      leg = TLegend(0.63,0.72,0.92,0.90)
      if 'LHE' in histoName:
        leg = TLegend(0.67,0.77,0.92,0.90)

      leg.Clear()
      leg.SetTextSize(0.03)
      if bZip!=None and len(bZip)>1: # need more columns if there are backgrounds
        leg = TLegend(0.6,0.7,0.99,0.92)
        leg.SetTextSize(0.02)
        leg.Clear()
        leg.SetNColumns(2)
        leg.SetTextSize(0.04)

      if h1!=None and h1.GetEntries()!=0:
        #print 'DBG lpe option in legentry', h1, h1.GetEntries()
        leg.AddEntry(h1,name1, "lpe")

      if bZip!=None:
        #print ' if bZip!=None'
        #print bZip
        stack = makeStack(bZip, myDir, histoName, leg, lumi, howToScale, norm1)

      if f3!=None and h3!=None:
        if doOverflow:
          handleOverflowBins(h3)
        h3.Draw("sames hist")
        h3.SetLineWidth(2)
        #h3.SetLineColor(kRed-9)
        col = getColors('signal-'+getSelection()[0:2])
        if 'J/Psi' in name3:
          col = getColors('signal-jp')
        h3.SetLineColor(col[0])
        h3.SetFillColor(col[1])

        norm3 = h3.Integral()
        scale = 0

        extract_scale_from_name = re.findall(r'\d+', name3)
        if howToScale=="lumi":
          # print "extract from name =", extract_from_name
          if len(extract_scale_from_name) != 0:
            scale = int(extract_scale_from_name[0])
            h3.Scale(scale)
        elif howToScale=="norm" and norm3!=0:
          h3.Scale(1./norm3)
        elif howToScale=="toData" and norm3!=0:
          h3.Scale(norm1/norm3)
        elif howToScale=="toDataInt" and norm3!=0:
          scale = int(norm1/norm3/roundTo)*roundTo
          if scale==0: scale=roundTo ## Don't want zero scalings...
          h3.Scale(scale)

        #print howToScale, "norm1=%.3f, norm3=%.3f, roundTo=%.3f, scale=%.3f"%(norm1, norm3, roundTo, scale)
        #if 'tri_mass125' in histoName:
        #  h3.Print('all')


        hmaxs.append(h3.GetMaximum())

        if howToScale=="toDataInt":
          leg.AddEntry(h3,name3.replace('XX',str(scale)), "f")
        else:
          leg.AddEntry(h3,name3, "f")
        if h1==None:
          mainHist=h3

        if doRatio:
          c1.cd()
          pad2.SetBottomMargin(0.25);
          pad2.SetTopMargin(0);
          pad2.Draw();
          pad2.cd();

          r = h1.Clone("")
          r.Divide(h3);
          r.GetYaxis().SetTitle("Data/MC");
          r.SetMaximum(2);
          r.SetMinimum(0);
          r.GetYaxis().SetNdivisions(206);
          r.GetYaxis().SetTitleOffset(0.4);
          r.SetTitleSize(0.1,"XYZ");
          r.SetLabelSize(0.1,"XY");

          r.Draw("e1p");

          pad1.cd();

      # mainHist.Draw('hist')
      mainHist.Draw('e1p')
      if bZip!=None:
        if howToScale=='lumi':
          stack.Draw('same hist')
          hmaxs.append(stack.GetMaximum())
        else:
          stack.Draw('same nostack hist')

      if h3!=None:
        h3.Draw("sames hist")
      if h1!=None:
        h1.Draw("same e1p")
        #h1.Draw("same e1p hist")
        # if howToScale=='lumi': h1.Draw("same e1p hist")
        # else: h1.Draw("same hist")

      if len(hmaxs)>0:
        m = max(hmaxs)
        mainHist.SetMaximum(1.2*m)
        if isLog:
          mainHist.SetMaximum(10*m)
      if howToScale == "norm":
        mainHist.GetYaxis().SetTitle("arbitrary units")
      elif howToScale in ["toData",'toDataInt']:
        mainHist.GetYaxis().SetTitle("Events")
        mainHist.SetMaximum(int(1.1*m)+5)
      else:
        mainHist.GetYaxis().SetTitle("Events")

      mainHist.SetMinimum(0)
      # Here we consider particular cases (histograms) that need special care
      # if "tri_mass_longTail" in histoName and howToScale=="lumi":
      # h1.SetMaximum(750)
      if 'tri_mass' in histoName:
        blindIt(h1, 3)
      if '_mDalG_' in histoName:
        blindIt(h1, 4)

      if histoName in ['LHE_diLep_mass_low', 'LHE_diLep_mass',
                       'LHE_dR_l1_l2','LHE_dR_l1_l2_low', 'LHE_dR_l1_l2_vlow']:
        nBins = h1.GetNbinsX()
        int1 = h1.Integral(int(0.25*nBins),nBins)
        int3 = h3.Integral(int(0.25*nBins),nBins)
        print histoName, nBins, int1, int3
        h1.Scale(50./int1)
        h3.Scale(50./int3)
        mainHist.SetMaximum(1.2*h1.GetMaximum())
        mainHist.SetNdivisions(505,'X')
        doPdf=1

      if "diLep_mass_0to20" in histoName:
        mainHist.SetMaximum(1.2*h1.GetMaximum())
        doPdf=1
      if "diLep_mass_jpsi" in histoName:
        mainHist.SetMaximum(1.1*h1.GetMaximum())
        doPdf=1
      if "diLep_mass_low" in histoName:
        mainHist.SetXTitle("m_{ee} (GeV)")
        mainHist.SetMaximum(1.2*h1.GetMaximum())
        doPdf=1

      if histoName in ['gen_co3','gen_phi']:
        mainHist.SetMaximum(0.04)
        mainHist.SetMinimum(0)


      '''
      if "ph_energyCorrection" in histoName:
        gStyle.SetOptStat(1111)
        h1.SetName("Data")
        h3.SetName("ggH-125")
        stats1 = h1.GetListOfFunctions().FindObject("stats");
        stats3 = h3.GetListOfFunctions().FindObject("stats");
        stats3.Print()
        stats1.SetX1NDC(0.7)
        stats1.SetX2NDC(0.95)
        stats1.SetY1NDC(0.9)
        stats1.SetY2NDC(0.7)
        stats3.SetX1NDC(0.7)
        stats3.SetX2NDC(0.95)
        stats3.SetY1NDC(0.7)
        stats3.SetY2NDC(0.5)
        stats3.SetTextColor(kRed+1)
        leg.SetX1(0.19)
        leg.SetX2(0.46)
        leg.SetY1(0.85)
        leg.SetY2(0.7)
      '''
      if doFits:
        if 'diLep_mass_jpsi_' in histoName:
          mllFrame = mll.frame()
          if h1!=None and h3!=None:
            mll.setRange('fitRange',2.9,3.3)

            dh1 = RooDataHist("dh1","dh1", RooArgList(mll), h1)
            voigt1 = RooVoigtian("voigt1","voigt1",mll,mean_mll,width_mll,sigma_mll)
            #expo   = RooExponential("expo", "Bkg Model", mll, lamb)
            #frac = RooRealVar("frac","frac",0.7, 0.5,1)
            #SplusB = RooAddPdf("SplusB", "Bkg plus signal peak", voigt1, expo, frac)

            #fit1   = SplusB.fitTo(dh1, RooFit.Range('fitRange'), RooFit.Save())
            fit1   = voigt1.fitTo(dh1, RooFit.Range('fitRange'), RooFit.Save())
            l1_0 = RooArgSet(mll)
            pars = voigt1.getParameters(l1_0)

            dh1.plotOn(mllFrame, RooFit.Name('data'))
            #SplusB.plotOn(mllFrame, RooFit.LineColor(kGreen+3))
            #SplusB.paramOn(mllFrame, RooFit.Label("Data"), RooFit.Layout(0.095,0.45,0.65))
            voigt1.plotOn(mllFrame, RooFit.LineColor(kGreen+3))
            voigt1.paramOn(mllFrame, RooFit.Label("Data"), RooFit.Layout(0.095,0.45,0.65))

            l1_1, l1_2 = RooArgList(mll), RooArgList(pars)
            ff1 = voigt1.asTF(l1_1, l1_2)

            dh3 = RooDataHist("dh3","dh3", RooArgList(mll), h3)
            voigt3 = RooVoigtian("voigt3","voigt3",mll,mean_mll,width_mll,sigma_mll)
            fit3   = RooFitResult(voigt3.fitTo(dh3, RooFit.Range('fitRange'),  RooFit.Save()))
            l3_0 = RooArgSet(mll)
            pars = voigt3.getParameters(l3_0)
            dh3.plotOn(mllFrame, RooFit.Name('sig'), RooFit.MarkerStyle(23), RooFit.MarkerColor(kBlue), RooFit.LineColor(kBlue))
            voigt3.plotOn(mllFrame)
            voigt3.paramOn(mllFrame, RooFit.Label("Signal MC"), RooFit.Layout(0.65, 1.0,0.65))

            l3_1, l3_2 = RooArgList(mll), RooArgList(pars)
            ff3 = voigt3.asTF(l3_1, l3_2)

            ff1.Print('all')
            ff3.Print('all')
            NBINS = 100
            hf1 = TH1F("hf1",";mll;hf1", NBINS, 2.9, 3.3)
            hf3 = TH1F("hf3",";mll;hf3", NBINS, 2.9, 3.3)
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

            #hratio.Print('all')
            #dhratio = RooDataHist("dhratio","dhratio", RooArgList(mll), hratio)
            #raw_input('Hit Enter to continue')

            #dhratio.plotOn(mllFrame, RooFit.Name('hratio'), RooFit.MarkerColor(kOrange))

            mllFrame.SetTitle(';m_{#mu#mu} (GeV);Events')
            mllFrame.Draw()
            mllFrame.SetMinimum(0.01)
            #mllFrame.SetMaximum(6)
            leg.Clear()
            leg.AddEntry(mllFrame.findObject('data'),'Data','e1p')
            leg.AddEntry(mllFrame.findObject('sig'),'h#rightarrow J/#Psi#gamma#rightarrow#mu#mu#gamma','pl')


        if 'res_' in histoName:
          gStyle.SetOptFit(0)
          func = TF1("func", "gaus", -0.1, 0.1);
          func.SetLineColor(kGreen+2)
          func.SetLineWidth(2)
          r = h3.Fit('func','RS', 'sames')
          h3.Draw('sames')
          lat = TLatex()
          lat.SetNDC()
          lat.DrawLatex(0.2,0.67, '#color[8]{Gaus fit}:')
          lat.DrawLatex(0.2,0.62, '#mu = %.3f'%(r.Parameter(1)))
          lat.DrawLatex(0.2,0.57, '#sigma = %.3f'%(r.Parameter(2)))
          lat.DrawLatex(0.2,0.50, '#chi^{2} = %.0f'%(r.Chi2()))
          # print ' Fit parameters: ', r.Parameter(1), r.Parameter(2)

      if "lPt2_pt_" in histoName:
        if not 'el' in getSelection():
          h1.SetXTitle("Trailing muon p_{T} (GeV)")
        mm = h1.GetMaximum()
        h1.SetMaximum(1.2*mm)
        h1.SetAxisRange(0,60,"X")
        doPdf=1
      '''
      if 'diLep_mass_jpsi_' in histoName:
        gStyle.SetOptStat(1111)
        stats1 = h1.GetListOfFunctions().FindObject("stats");
        stats2 = h3.GetListOfFunctions().FindObject("stats");
        stats2.Print()
        stats1.SetX1NDC(0.7)
        stats1.SetX2NDC(0.95)
        stats1.SetY1NDC(0.9)
        stats1.SetY2NDC(0.7)
        stats2.SetX1NDC(0.7)
        stats2.SetX2NDC(0.95)
        stats2.SetY1NDC(0.7)
        stats2.SetY2NDC(0.5)
        stats2.SetTextColor(kBlue+1)

        leg.SetX1(0.19)
        leg.SetX2(0.46)
        leg.SetY1(0.85)
        leg.SetY2(0.7)
      '''

      if "zee_" in histoName:
        gStyle.SetOptStat(1111)
        h1.SetName("Data")
        h2.SetName("DYjets")
        # h2.Print("all")
        stats1 = h1.GetListOfFunctions().FindObject("stats");
        stats2 = h2.GetListOfFunctions().FindObject("stats");
        stats2.Print()
        stats1.SetX1NDC(0.7)
        stats1.SetX2NDC(0.95)
        stats1.SetY1NDC(0.9)
        stats1.SetY2NDC(0.7)
        stats2.SetX1NDC(0.7)
        stats2.SetX2NDC(0.95)
        stats2.SetY1NDC(0.7)
        stats2.SetY2NDC(0.5)
        stats2.SetTextColor(kBlue+1)
        leg.SetX1(0.19)
        leg.SetX2(0.46)
        leg.SetY1(0.85)
        leg.SetY2(0.7)

      gPad.RedrawAxis()

      lat = TLatex()
      lat.SetNDC()
      lat.SetTextSize(0.035)
      if 'J/Psi' in name3:
        lat.DrawLatex(0.18,0.95, 'H #rightarrow J/#Psi#gamma#rightarrow#mu#mu#gamma')
      else:
        if 'el' in getSelection():
          lat.DrawLatex(0.18,0.95, 'H #rightarrow#gamma*#gamma#rightarrow ee#gamma')
        else:
          lat.DrawLatex(0.18,0.95, 'H #rightarrow#gamma*#gamma#rightarrow#mu#mu#gamma')

      CMS_lumi(c1, 2, 11)

      # print "hmax =", hmaxs

      c1.cd()
      leg.SetFillColor(kWhite)
      leg.Draw()

      c1.SetLogy(int(isLog))
      #print 'Saving PNG:', histoName
      c1.SaveAs(path+"/"+histoName+'.png')
      if doPdf:
        c1.SaveAs(path+"/"+histoName+'.pdf')
        doPdf=0
      gStyle.SetOptStat(0)

    c1.SetLogy(0)

def yieldsTable(yieldList, names, num=True):
  #print sel
  print yieldList
  t = []
  cuts =  getCuts()

  l1 = ["Cut"]
  if num:
    l1.insert(0,"")
  l1.extend(names)

  # print 'First line in yieldsTable:', l1

  t.append(l1)
  for line in xrange(len(cuts)):
    l=[]
    #print line, cuts[line]
    if num:
      l.append(str((line-1)))

    if line==0:l.append('Initial Events')
    else:      l.append(cuts[line-1])

    for yi in yieldList:
      #print line, yi
      l.append(yi[line])
      # print l
    t.append(l)

  # print t
  return t

def makeTable(table, name, opt="tex", precision='%.2f'):
    print "Making sure that the list is alright"
    n_row = len(table)
    if n_row==0: return
    n_col = len(table[0])

    # print table
    for l in table:
      if len(l)!=n_col:
        print "No good, the number of columns is messed up", len(l), n_col


    myTable = ''
    if opt=="tex":
      beginTable = '\\begin{tabular}{|'+n_col*"l|"+'} \n \\hline \n'
      endTable   = '\\end{tabular} \n'

      beginLine  = ''
      endLine    = ' \\\\ \n'
      #endLine    = ' \\\\ \\hline \n'
      separator  = ' & '


    if opt=="html":
      beginTable = '<table border = "10" cellpadding="5">'
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
      # print l, table[l]
      myTable+=beginLine
      for c in range(n_col):
        val = table[l][c]
        if isinstance(val,float):
          if c==2:
            myTable+="%.0f" % (val)
          else:
            myTable+=precision % (val)
        else:
          myTable+=str(val)
        if c!=n_col-1:
          myTable+=separator

      myTable+=endLine

    myTable +=endTable

    ifile = open("yields_"+name+"."+opt,"w")
    ifile.write(myTable)
    ifile.close()

    if opt in ["twiki"]:
      print myTable

def getYields(f, sample='ggH-125', doLumiScale=False):
  print 'Calculating yields for ',sample
  ev   = f.Get("Counts/evt_byCut")
  sel  = conf.get("selection", "sel")[0:2]
  cuts = getCuts()
  if ev==None: return len(cuts)*[0]

  y = []
  scale=1
  if doLumiScale: # only assume signal MC for now
    Nev = getTotalEvents(f)

    if any(substring in sample for substring in ['ggH','vbfH','vH']):
      cro = getCS(sample, sel)
    else:
      cro = getCS(sample)
    scale = float(lumi*cro)/Nev
    print "Lumi scale for sample =", sample, ", sel=",sel, "Nev=",Nev, "cro=",cro, "scale=",scale

  for a in xrange(len(cuts)):
    y.append(scale*ev.GetBinContent(a+1)) # well, that's how the histogram is set up

  return y


def effPlots2(f1, path):
  c1.cd()
  hh = []
  ra = []
  for n in xrange(0,31):
    print n
    hh.append(f1.Get("eff/gen_Mll_"+str(n)))
    ra.append(hh[-1].Clone())
    ra[-1].Divide(hh[0])

  ra[1].Draw("hist")
  ra[1].SetMinimum(0)
  ra[1].SetMaximum(1)
  ra[1].SetTitle(";Mll; acc * eff")
  bcol = 29
  ra[1].SetLineColor(bcol)
  ra[1].SetLineWidth(2)
  leg = TLegend(0.70,0.65,0.95,0.98);
  leg.AddEntry(ra[1],"pT(ll) > 25 GeV", "l")
  for a in xrange(2,7):
    ra[a].Draw('same hist')
    ra[a].SetLineColor(bcol+2*(a-2))
    ra[a].SetLineWidth(2)
    leg.AddEntry(ra[a],"pT(ll) > %.0f GeV"%(25+5*(a-1)), "l")

  leg.SetTextSize(0.03)
  leg.SetFillColor(kWhite)
  leg.Draw()
  c1.SaveAs(path+"/acc_eff_Mll_ptll.png")

  st = 7
  ed = 12
  ra[st].Draw("hist")
  ra[st].SetMinimum(0)
  ra[st].SetMaximum(1)
  ra[st].SetTitle(";Mll; acc")
  bcol = 33
  ra[st].SetLineColor(bcol)
  ra[st].SetLineWidth(2)
  leg = TLegend(0.70,0.65,0.95,0.98);
  leg.AddEntry(ra[st],"pT(#gamma) > 25 GeV", "l")
  for a in xrange(st+1,ed+1):
    ra[a].Draw('same hist')
    ra[a].SetLineColor(bcol+2*(a-st))
    ra[a].SetLineWidth(2)
    leg.AddEntry(ra[a],"pT(#gamma) > %.0f GeV"%(25+5*(a-st)), "l")


  leg.SetTextSize(0.03)
  leg.SetFillColor(kWhite)
  leg.Draw()
  c1.SaveAs(path+"/acc_eff_Mll_ptgamma.png")

  st = 13
  ed = 18
  ra[st].Draw("hist")
  ra[st].SetMinimum(0)
  ra[st].SetMaximum(1)
  ra[st].SetTitle(";Mll; acc")
  bcol = 38
  ra[st].SetLineColor(bcol)
  ra[st].SetLineWidth(2)
  leg = TLegend(0.70,0.65,0.95,0.98);
  leg.AddEntry(ra[st],"pT(#gamma)/Mllg > 0.20", "l")
  for a in xrange(st+1,ed+1):
    ra[a].Draw('same hist')
    ra[a].SetLineColor(bcol+2*(a-st))
    ra[a].SetLineWidth(2)
    leg.AddEntry(ra[a],"pT(#gamma)/Mllg > %.2f"%(0.20+0.05*(a-st)), "l")

  leg.SetTextSize(0.03)
  leg.SetFillColor(kWhite)
  leg.Draw()
  c1.SaveAs(path+"/acc_eff_Mll_ptgammaMllg.png")

  st = 19
  ed = 24
  ra[st].Draw("hist")
  ra[st].SetMinimum(0)
  ra[st].SetMaximum(1)
  ra[st].SetTitle(";Mll; acc")
  bcol = 41
  ra[st].SetLineColor(bcol)
  ra[st].SetLineWidth(2)
  leg = TLegend(0.70,0.65,0.95,0.98);
  leg.AddEntry(ra[st],"pT(ll)/Mllg > 0.20", "l")
  for a in xrange(st+1,ed+1):
    ra[a].Draw('same hist')
    ra[a].SetLineColor(bcol+2*(a-st))
    ra[a].SetLineWidth(2)
    leg.AddEntry(ra[a],"pT(ll)/Mllg > %.2f"%(0.20+0.05*(a-st)), "l")

  leg.SetTextSize(0.03)
  leg.SetFillColor(kWhite)
  leg.Draw()
  c1.SaveAs(path+"/acc_eff_Mll_ptllMllg.png")


  st = 25
  ed = 30
  ra[st].Draw("hist")
  ra[st].SetMinimum(0)
  ra[st].SetMaximum(1)
  ra[st].SetTitle(";Mll; acc")
  bcol = 43
  ra[st].SetLineColor(bcol)
  ra[st].SetLineWidth(2)
  leg = TLegend(0.70,0.65,0.95,0.98);
  leg.AddEntry(ra[st],"pT(both)/Mllg > 0.20", "l")
  for a in xrange(st+1,ed+1):
    ra[a].Draw('same hist')
    ra[a].SetLineColor(bcol+2*(a-st))
    ra[a].SetLineWidth(2)
    leg.AddEntry(ra[a],"pT(both)/Mllg > %.2f"%(0.20+0.05*(a-st)), "l")

  leg.SetTextSize(0.03)
  leg.SetFillColor(kWhite)
  leg.Draw()
  c1.SaveAs(path+"/acc_eff_Mll_ptbothMllg.png")


def effPlots(f1, path):
  print "Now making efficiency plots"
  # f1.cd(myDir)

  for var in ["Mll","dR"]:

    h0 = f1.Get("eff/gen_"+var+"_0")
    h1 = f1.Get("eff/gen_"+var+"_acc_gamma")
    h2 = f1.Get("eff/gen_"+var+"_acc_lept")

    h3 = f1.Get("eff/gen_"+var+"_reco_gamma_iso")
    h4 = f1.Get("eff/gen_"+var+"_two_lep_reco")


    h0.Draw("hist")
        #h1.Draw("hist same")
        #h2.Draw("hist same")
        #h3.Draw("hist same")
    c1.SaveAs(path+var+".png")

        #for acceptance
    r1 = h1.Clone()
    r1.Divide(h0)
    r2 = h2.Clone()
    r2.Divide(h0)
        #for reco eff
    r3 = h3.Clone()
    r3.Divide(h2)
    r4 = h4.Clone()
    r4.Divide(h2)
        #r5 = h5.Clone()
        #r5.Divide(h2)
        #r6 = h6.Clone()
        #r6.Divide(h2)
        #r7 = h7.Clone()
        #r7.Divide(h2)


    if var=="Mll":
      xname = ";M(l1,l2)"
      if var=="dR":
        xname = ";dR(l1,l2)"

      r1.Draw("hist")
      r2.Draw("hist same")
      r1.SetMinimum(0)
      r1.SetMaximum(1)
      r1.SetTitle(xname+" gen; acc")
      r1.SetLineColor(kRed+1)
      r2.SetLineColor(kGreen+1)
      leg = TLegend(0.20,0.2,0.90,0.30);
      leg.AddEntry(r1,"photon pt>25, eta<2.5", "l")
      leg.AddEntry(r2,"photon pt>25 and p_{T}(l1)>23, p_{T}(l2)>4", "l")
      leg.SetTextSize(0.04)
      leg.SetFillColor(kWhite)
      leg.Draw()
      c1.SaveAs(path+"acceptance_"+var+".png")

      r3.Draw("hist")
      r3.SetMinimum(0)
      r3.SetMaximum(1)
      r3.SetTitle(xname+" gen; reco eff")
      r4.Draw("hist same")
        #r5.Draw("hist same")
        #r6.Draw("hist same")
        #r7.Draw("hist same")

      r3.SetLineColor(kBlack)
      r4.SetLineColor(kOrange+1)
        #r5.SetLineColor(kGreen+1)
        #r6.SetLineColor(kRed+1)
      leg = TLegend(0.35,0.15,0.98,0.30);
      leg.AddEntry(r3,"Reco photon pt>25, Medium WP ID ", "l")
      leg.AddEntry(r4,"Photon and two reco-muons (OUR ID/ISO)", "l")
        #leg.AddEntry(r5,"Photon + 2 electrons NO ID pt(e1,e2) > (23,7)", "l")
        #leg.AddEntry(r6,"Photon + 2 electrons ELE ID pt(e1,e2) > (23,7)", "l")
        #leg.AddEntry(r7,"Photon + ONE electron pt(e) > 30, No ID", "l")
      leg.SetTextSize(0.025)
      leg.SetFillColor(kWhite)
      leg.Draw()
      c1.SaveAs(path+"eff_"+var+".png")




def AFB(f, s, name, path):
  c1 = TCanvas("c1","a canvas",600,600);
  cosCS  = f.Get("Angles/"+name)
  nBinsX = cosCS.GetNbinsX()
  nBinsY = cosCS.GetNbinsY()
  #cosCS.Print()
  cosCS.Draw()
  #print "In X: from %.1f to %.1f in %i bins" % (cosCS.GetXaxis().GetBinLowEdge(1), cosCS.GetXaxis().GetBinUpEdge(nBinsX), nBinsX)
  #print "In Y: from %.1f to %.1f in %i bins" % (cosCS.GetYaxis().GetBinLowEdge(1), cosCS.GetYaxis().GetBinUpEdge(nBinsY), nBinsY)
  mMin = cosCS.GetYaxis().GetBinLowEdge(1)
  mMax = cosCS.GetYaxis().GetBinUpEdge(nBinsY)
  mBin = (mMax - mMin)/nBinsY
  pr = []
  step = 10
  nSteps = nBinsY/step
  import numpy as np
  m   = np.zeros(nSteps,  dtype=np.float)
  B   = np.zeros(nSteps,  dtype=np.float)
  F   = np.zeros(nSteps,  dtype=np.float)
  eB  = np.zeros(nSteps,  dtype=np.float)
  eF  = np.zeros(nSteps,  dtype=np.float)
  Afb = np.zeros(nSteps,  dtype=np.float)
  eA  = np.zeros(nSteps,  dtype=np.float)

  histmax = 0
  #pr.append(cosCS.ProfileX('pr_step_'+str(step), 1, 20))
  for i in xrange(nSteps):
    pr.append(cosCS.ProjectionX(s+'pr_step_'+str(i), i*step+1, (i+1)*step))
    histmax = max(histmax,pr[-1].GetMaximum())
    #pr[-1].Rebin()
    m[i] = mMin + (i+0.5)*step*mBin
    e = ROOT.Double()
    #B[i] = pr[-1].IntegralAndError(1,50, e)
    eB[i] = e
    #F[i] = pr[-1].IntegralAndError(51,100,eF[i])
    if (F[i]+B[i])!=0:
      Afb[i] = (F[i]-B[i])/(F[i]+B[i])
      eA[i] = eF[i] + eB[i]
    else:
      Afb[i] = 0
  print 'F, B and  Afb',F,B, Afb
  gr = TGraphErrors(nSteps, m, Afb, eA, eA)
  gr.Draw("APL")
  gr.SetMarkerStyle(23)
  gr.SetMarkerSize(0.66)
  gr.SetTitle(';m (GeV); A_{FB}')
  #gr.Print('all')
  c1.SaveAs(path+'/'+name+'-AFB-vs-M'+s+'.png')


  pr[0].Draw('hist')
  pr[0].SetMaximum(1.2*histmax)
  for i in xrange(1,nBinsY/step):
    pr[i].Draw('same hist')
    pr[i].SetLineColor(22+i)
  c1.SaveAs(path+'/'+name+'-projections'+s+'.png')

if __name__ == "__main__":
  print "This is utils.py script"
