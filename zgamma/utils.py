#! /usr/bin/env python
from optparse import OptionParser
import sys,os,datetime,re
from array import *
from ROOT import *
gROOT.SetBatch()
TH1.SetDefaultSumw2(kTRUE)

import ConfigParser as cp
conf = cp.ConfigParser()
conf.read('../zgamma/config.cfg')
lumi2012 = float(conf.get("lumi","lumi2012A")) + float(conf.get("lumi","lumi2012B"))+\
    float(conf.get("lumi","lumi2012C")) + float(conf.get("lumi","lumi2012D"))
lumi = lumi2012

def setSelection(sel):
  conf.set("selection","sel", sel)
def getSelection():
  conf.get("selection","sel")

def getCuts(cp, name="cuts-mu"):
  cuts = []
  for key, cut in sorted(cp.items(name)):
    cuts.append(cut)
  return cuts
# print key, cut

def getLumi(period):
  if period=="2012":
    return lumi2012
  else:
    print "Sorry only 2012 is considered"
    return 0

def getCS(sample, mySel=None):
  if mySel==None:
    #sel = conf.get("selection", "sel")[0:2]
    cs = float(conf.get(sample, "cs"))
  else:
    sel = mySel
    cs = float(conf.get(sample, "cs-"+sel))

  return cs

def getColors(sample):
  l = conf.get(sample, "line")
  f = conf.get(sample, "fill")
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

def blindIt(h):
  hbin =  h.FindBin(125)
  for b in xrange(hbin-10,hbin+10):
    h.SetBinContent(b,0)

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
    try:
      os.makedirs(myDir)
    except OSError:
      if os.path.isdir(myDir):
        pass
      else:
        raise


def makeStack(bZip, histoName, leg, lumi):
  hs = THStack("temp", "Stacked histo")

  #samplesToAdd = []

  for n,f in bZip:
    print f,n

    h = f.Get(histoName).Clone()
    if h==None:
      print 'None histogrammm:',histoName
      continue

    # if doOverflow:
    # handleOverflowBins(h)

    Nev = getTotalEvents(f)
    cro = getCS(n)
    scale = float(lumi*cro)/Nev
    print n, Nev, lumi, cro, scale

    h.Scale(scale)
    h.SetLineColor( int(getColors(n)[0]))
    h.SetFillColor( int(getColors(n)[1]))
    hs.Add(h)
    leg.AddEntry(h, conf.get(n,"shortName") ,"f")

    #hm.Print('all')

  return hs


def drawAllInFile(f1, name1, bZip, f3, name3, myDir,path, N, howToScale="none", isLog=False, doRatio=False):
    print 'myDir is ', myDir
    f1.cd(myDir)
    dirList = gDirectory.GetListOfKeys()
    # dirList.Print()
    createDir(path)
    scale2 = scale3 = 1
    split = os.path.split(path.rstrip("/"))
    print "Split lit lit", split[0], split[1]

    doOverflow = 0
    if doRatio:
        c1 = TCanvas("c2","big canvas",600,700);
    else:
        c1 = TCanvas("c3","small canvas",600,600);
    c1.SetLogy(isLog)
    c1.cd()

    if f3!=None and howToScale=="lumi": # only assume signal MC for now
        Nev = f3.Get("Counts/evt_byCut_raw").GetBinContent(1)
        print ' Fix ME FIX ME'
        cro = getCS("ZtoJPsiGamma") #
        print 'need more generic way to get cs by sample'
        #cro = getCS("ggH-125")
        scale3 = float(lumi*cro)/Nev
        print Nev, lumi, cro, scale3

    doPdf=0
    histoName = None
    for k1 in dirList:
        if k1.GetName() in ["eff"]: continue
        if N!=None:
            if not("cut"+N) in k1.GetName(): continue
        h1 = k1.ReadObj()
        histoName = h1.GetName()

        h3 = TH1F()

        hmaxs = []

        if f3!=None:
            if myDir!="":
                h3 = f3.Get(myDir+"/"+histoName)
            else:
                h3 = f3.Get(histoName)
            if h3==None: continue

            h3.Scale(float(scale3))


        print "drawing", histoName

        if h1.InheritsFrom("TH2"):
            createDir(split[0]+"/TH2_"+split[1])
            h1.Draw("col")
            prename = split[0]+"/TH2_"+split[1]+"/"+histoName
            c1.SaveAs(prename+"_data.png")
            if f3!=None and h3!=None:
                h3.Draw("col")
                c1.SaveAs(prename+"_sig.png")

            continue

            #set_palette("gray")
            #h1.Draw("col same")
            #if "h2D_tri_vs_diLep_mass" in histoName:
            #    c1.SetLogy()

        else:
            pad1 = TPad("pad1","pad1",0,0.3,1,1);
            pad2 = TPad("pad2","pad2",0,0,1,0.3);

            if doRatio:
                pad1.SetBottomMargin(0);
                pad1.Draw();
                pad1.cd();
                pad1.SetLogy(isLog)


            if doOverflow:
                handleOverflowBins(h1)

            h1.Draw("hist")
            h1.Draw("same e1p")
            #h1.SetMinimum()
            #h1.GetRange
            h1.SetMarkerStyle(20)
            h1.SetMarkerSize(0.75)
            h1.SetLineColor(kBlack)
            #h1.UseCurrentStyle()

            if "phi" in histoName:
              h1.SetMinimum(0)
            if "ll_deltaR_" in histoName:
              # h1.SetAxisRange(0,1,"X")
              h1.SetXTitle("#DeltaR(#mu_{1}, #mu_{2})")
              doPdf=1
            if "gamma_pt__" in histoName:
              h1.SetXTitle("Photon p_{T} (GeV)")
              doPdf=1
            if "lPt1_pt__" in histoName:
              h1.SetXTitle("Leading muon p_{T} (GeV)")
              doPdf=1
            if "lPt2_pt__" in histoName:
              h1.SetXTitle("Trailing muon p_{T} (GeV)")
              doPdf=1


            #if "tri_mass" in k1.GetName() and name1 not in ["madgra","mcfm"]:
            #    blindIt(h1)
            if "h_mass" in k1.GetName():
                print "\n *** H-mass RMS:",  h1.GetRMS(), h2.GetRMS()
                print "\n *** H-mass Mean:", h1.GetMean(),h2.GetMean()

            leg = TLegend(0.63,0.72,0.92,0.90);
            leg.SetTextSize(0.03)

            if len(bZip)!=0:
              leg = TLegend(0.55,0.7,0.99,0.92);
              leg.SetNColumns(2)
              leg.SetTextSize(0.04)
            leg.AddEntry(h1,name1, "lpe")


            #if "phi" in histoName:
            #    h1.SetMinimum(0)
            #    if howToScale=="norm":
            #        h1.SetMaximum(0.035)

            norm1 = h1.Integral()
            if howToScale=='norm' and  norm1!=0:
                h1.Scale(1./norm1)
            hmaxs.append(h1.GetMaximum())

            if len(bZip)!=0:
              print ' if bZip!=None!'
              print bZip

              if myDir!="":
                histoName = myDir+"/"+histoName


              stack = makeStack(bZip, histoName, leg, lumi)
              stack.Draw('same hist')
              hmaxs.append(stack.GetMaximum())

            if f3!=None and h3!=None:
                if doOverflow:
                    handleOverflowBins(h3)
                h3.Draw("sames hist")
                h3.SetLineWidth(2)
                h3.SetLineColor(kRed-9)
                h3.SetFillColor(kYellow-9)
                #h3.SetFillColor(kYellow-9)
                h1.Draw("same hist")
                h1.Draw("same e1p")
                h1.SetMarkerStyle(20)
                h1.SetMarkerSize(0.75)
                extract_from_name = re.findall(r'\d+', name3)
                #print extract_from_name
                if len(extract_from_name) == 1:
                    h3.Scale(int(extract_from_name[0]))
                norm3 = h3.Integral()
                if howToScale=="norm" and  norm3!=0:
                    h3.Scale(1./norm3)
                elif howToScale=="norm2" and  norm3!=0:
                    h3.Scale(norm1/norm3)
                hmaxs.append(h3.GetMaximum())

                leg.AddEntry(h3,name3, "f")

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

            h1.Draw("same hist")
            h1.Draw("same e1p")
            # prelim = TLatex(0.15,0.95, "CMS Preliminary %s #it{L_{int}} = %0.1f fb^{-1}" % (8, lumi))
            prelim = TLatex(0.15,0.95, "CMS Preliminary #sqrt{s} = 8TeV, L = 19.7 fb^{-1}   H#rightarrow#gamma*#gamma#rightarrow#mu#mu#gamma")
            prelim.SetNDC();
            prelim.SetTextSize(0.035);
            prelim.Draw();

            #print "hmax =", hmaxs
            if len(hmaxs)>0:
                m = max(hmaxs)
                h1.SetMaximum(1.1*m)
            if howToScale ==" norm":
                h1.GetYaxis().SetTitle("arbitrary units")
            if howToScale == "norm2":
                h1.GetYaxis().SetTitle("Events")
                h1.SetMaximum(int(1.1*m)+5)

            # Here we consider particular cases (histograms) that need special care
            # if "tri_mass_longTail" in histoName and howToScale=="lumi":
            # h1.SetMaximum(750)

            if "diLep_mass_low_" in histoName:
                h1.Rebin(2)
                h3.Rebin(2)
                h1.SetXTitle("m_{#mu#mu} (GeV)")
                h1.SetMaximum(1.2*h1.GetMaximum())

                doPdf=1

            if "ph_energyCorrection" in histoName:
                gStyle.SetOptStat(1111)
                h1.SetName("Data")
                h3.SetName("ggH-125")
                stats1 = h1.GetListOfFunctions().FindObject("stats");
                stats3 = h3.GetListOfFunctions().FindObject("stats");
                stats1.Print()
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


            if "zee_" in histoName:
                gStyle.SetOptStat(1111)
                h1.SetName("Data")
                h2.SetName("DYjets")
                #h2.Print("all")
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

            c1.cd()
            leg.SetFillColor(kWhite)
            leg.Draw()

            c1.SetLogy(int(isLog))
            #if "diLep_mass_low":
            #    h1.

            c1.SaveAs(path+"/"+histoName+'.png')
            if doPdf:
                c1.SaveAs(path+"/"+histoName+'.pdf')
                doPdf=0
            gStyle.SetOptStat(0)

        c1.SetLogy(0)

def yieldsTable(yieldList, sel, num=True):
  print sel
  t = []
  cuts =  getCuts(conf, "cuts-"+sel[0:2])

  l1 = ["Cut/trigger"]
  if num:
    l1.insert(0,"")
  l1.extend(len(yieldList)*[sel])

  # print 'First line in yieldsTable:', l1

  t.append(l1)
  for line in xrange(len(cuts)):
    l=[]
    if num:
      l.append(str((line-1)))
    l.append(cuts[line])

    for yi in yieldList:
      l.append(yi[line])
      # print l
    t.append(l)

  # print t
  return t

def makeTable(table, name, opt="tex"):
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
      endLine    = ' \\\\ \\hline \n'
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
      # print l, table[l]
      myTable+=beginLine
      for c in range(n_col):
        val = table[l][c]
        if not isinstance(val,str):
          if c==2:
            myTable+="%.0f" % (val)
          else:
            myTable+="%.2f" % (val)
        else:
          myTable+=val
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
    ev = f.Get("Counts/evt_byCut")
    if ev==None: return [0]

    sel = conf.get("selection", "sel")[0:2]
    cuts =  getCuts(conf, "cuts-"+sel)
    y = []
    scale=1
    if doLumiScale: # only assume signal MC for now
        Nev = ev.GetBinContent(1)
        if any(substring in sample for substring in ['ggH','vbfH','vH']):
          cro = getCS(sample, sel)
        else:
          cro = getCS(sample)
        scale = float(lumi*cro)/Nev
        print "Lumi scale for sel=",sel, "Nev=",Nev, "cro=",cro, "scale=",scale

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

if __name__ == "__main__":
    print "This is utils.py script"
