#!/usr/bin/env python
import sys,os

from ROOT import *
import config as c


def  makeMvaTrees(outpath, bg_list, sig1_list, sig2_list, sel):
    print "\n  ** Creating mva trees ***"
    fbg_train = TFile(outpath+"allBg_train.root","RECREATE")
    fbg_test  = TFile(outpath+"allBg_test.root","RECREATE")

    gROOT.ProcessLine(
"struct mva_vars {\
Float_t lep1Pt;\
Float_t lep2Pt;\
Float_t diLepPt;\
Float_t diLepM;\
Float_t diLepEta;\
Float_t met;\
Float_t metOverQt;\
Float_t metProjOnQt;\
Float_t metPerpQt;\
Float_t dPhiMetDiLep;\
Float_t dPhiJetMet;\
Float_t mt;\
Int_t nJets;\
Int_t nJets15;\
Float_t deltaEtaDiJet;\
Float_t massDiJet;\
Float_t zeppDiJetDiLep;\
Float_t evWeight;\
Float_t fullWeight;\
};")

    allVars = 'lep1Pt/F:lep2Pt:diLepPt:diLepM:diLepEta:met:metOverQt:metProjOnQt:metPerpQt:\
dPhiMetDiLep:dPhiJetMet:mt:nJets/I:nJets15/I:deltaEtaDiJet/F:massDiJet:zeppDiJetDiLep:evWeight/F:\
fullWeight/F'
    stuff  = mva_vars()
    fbg_train.cd()
    mvaTree_train = TTree("mvaTree","mvaTree")
    br_train = mvaTree_train.Branch("mvaVars", stuff, allVars)

    fbg_test.cd()
    mvaTree_test = TTree("mvaTree","mvaTree")
    br_test = mvaTree_test.Branch("mvaVars", stuff, allVars)

    for b in bg_list:
        sample = b.Get("Andrey/evt_byCut").GetTitle()
        #if sample!="tW": continue  # for tests
       
        tree = b.Get("mvaTree/mvaTree")
        print sample, "Total events in a tree: ", tree.GetEntries()

        for evt in tree:
            fillStruct(stuff,evt,sel,sample)
            #print "evt", evt
            #stuff2 = stuff
            #print stuff.lep1Pt
            #print stuff2.lep1Pt
            if tree.evNumber % 2 == 0:
                mvaTree_train.Fill()
            else:
                mvaTree_test.Fill()
                #print tree.evNumber
                #print "lepton's pt:",stuff.lep1Pt,stuff.lep2Pt
    print "\t Nevts in train=", mvaTree_train.GetEntries(), "Nevts in test=",mvaTree_test.GetEntries()
        
    fbg_train.cd()      
    mvaTree_train.Write()
    fbg_train.Close()

    fbg_test.cd()      
    mvaTree_test.Write()
    fbg_test.Close()
    
    mhmap = {"125":0,"200":1,"250":2}
    for mh in [125,200,250]:
        fsig_train = TFile(outpath+"hzz"+str(mh)+"_train.root","RECREATE")
        fsig_test  = TFile(outpath+"hzz"+str(mh)+"_test.root","RECREATE")

        fsig_train.cd()
        mvaTree_train = TTree("mvaTree","mvaTree")
        br_train = mvaTree_train.Branch("mvaVars", stuff, allVars)

        fsig_test.cd()
        mvaTree_test = TTree("mvaTree","mvaTree")
        br_test = mvaTree_test.Branch("mvaVars", stuff, allVars)

        sig_list = TList()
        sig_list.Add(sig1_list.At(mhmap[str(mh)]))
        sig_list.Add(sig2_list.At(mhmap[str(mh)]))
        
        for s in sig_list:
            sample = s.Get("Andrey/evt_byCut").GetTitle()
            tree = s.Get("mvaTree/mvaTree")
            
            print sample, "Total events in a tree: ", tree.GetEntries()
            for evt in tree:
                fillStruct(stuff,evt,sel,sample)
                
                #stuff2 = stuff
                #print stuff.lep1Pt
                
                if tree.evNumber % 2 == 0:
                    mvaTree_train.Fill()
                else:
                    mvaTree_test.Fill()
                    #print tree.evNumber
                    #print "lepton's pt:",stuff.lep1Pt,stuff.lep2Pt
        print "\t Nevts in train=", mvaTree_train.GetEntries(), "Nevts in test=",mvaTree_test.GetEntries()
                
        fsig_train.cd()      
        mvaTree_train.Write()
        fsig_train.Close()
        
        fsig_test.cd()      
        mvaTree_test.Write()
        fsig_test.Close()

    print "** End of creating mva trees ***"

def fillStruct(stuff, evt, sel, sample):
    lumi = c.Params().getLumi(sel)
    xc = c.Params().xsec_and_colors()

    
    stuff.lep1Pt   = evt.lep1Pt
    stuff.lep2Pt   = evt.lep2Pt
    stuff.diLepPt  = evt.diLepPt
    stuff.diLepM   = evt.diLepM
    stuff.diLepEta = evt.diLepEta
    stuff.met      = evt.met
    stuff.metOverQt    =  evt.metOverQt
    stuff.metProjOnQt  =  evt.metProjOnQt
    stuff.metPerpQt    =  evt.metPerpQt
    stuff.dPhiMetDiLep = evt.dPhiMetDiLep
    stuff.dPhiJetMet   = evt.dPhiJetMet  
    stuff.mt      = evt.mt
    stuff.nJets   = evt.nJets
    stuff.nJets15 = evt.nJets15
    #print stuff.nJet15
    stuff.deltaEtaDiJet  = evt.deltaEtaDiJet
    stuff.massDiJet      = evt.massDiJet
    stuff.zeppDiJetDiLep = evt.zeppDiJetDiLep
    stuff.evWeight   = evt.weight
    stuff.fullWeight = evt.weight* float(lumi*xc[sample][2])/xc[sample][3]


#def CalculateTransMass(p1,p2):
#    transE    = sqrt(p1.Pt()*p1.Pt() + p2.M()*p2.M()) + sqrt(p2.Pt()*p2.Pt() +p2.M()*p2.M())
#    transPt   = (p1 + p2).Pt()
#    transMass = sqrt(transE*transE - transPt*transPt)
#    return transMass      

"""
def makeWeightBranch(outpatth,bg_list, sig1_list, sig2_list):
    fbg_train = TFile("allBg_train.root","RECREATE")
    #fbg_test =  TFile("allBg_test.root","RECREATE")

    tree = bg_list[0].Get("mvaTree/mvaTree")
    fbg_train.cd()
    tree.Write()
    '''
    gROOT.ProcessLine(
"struct stuff_t {\
Float_t fullWeight;\
};")

    www = stuff_t()

    allList = TList()
    for b in bg_list:
        tree = b.Get("mvaTree/mvaTree")
        newtree = tree.CloneTree(0)
        br = newtree.Branch("fullWeights",www,"fullWeight/F")
        #br.SetFile(fbg_train)

        for evt in tree:
            www.fullWeight = evt.weight * 10
            #print evt.weight
            newtree.Fill()
            
        allList.Add(newtree)
        #newtree.Print()


    allList.Print()
    mvaTree_train = TTree("mvaTree","mvaTree")
    mvaTree_train.MergeTrees(allList)
    #mvaTree_train.Print()

    
    #print "number of merged events:", a 
    '''
    #mvaTree_train = fbg_train.Get("mvaTree")
    
    #for evt in mvaTree_train:
    #    print "weights", evt.weight

    #fbg_train.cd()
    #mvaTree_train.AutoSave()
   
"""        


    

def printYields(topbg_list, bg_list, sig1_list, sig2_list, sig3_list, data, sel, filename, mode="scaled"):
    print "Yield are HTML mode"
    
    cutNames = ["0. Total",
                "1. HLT",
                "2. Two Leptons",
                "3. third lepton veto",
                "4. Z peak",
                "5. Z pt cut",
                "6. dPhi cut",
                "7. Met >70",
                "8. bjet veto",
                "9. H200 cut based",
                "10. H250 cut based",
                "11. H300 cut based",
                "12. H350 cut based",
                "13. H400 cut based",
                "14. H450 cut based",
                "15. H500 cut based",
                "16. H550 cut based",
                "17. H600 cut based",
                "18. H250 MVA",
                "19. 0j bin",
                "20. 1j bin",
                "21. >1j bin",
                "22.",
                "23.",
                "24.",
                "25. mvaTree, nj=0",
                "26. mvaTree, nj=1",
                "27. mvaTree, nj>1",
                ]

    lTot = 28   # Total number of cuts
    lMax = 9 # Lines to print in the main table
    lMax2 = lTot-lMax # additional table
    beginTable = '<table border = "10"    cellpadding="5">'
    endTable = '</table>'

    myTable = ''
    beginLine = '<tr>\n'
    endLine = '</tr>\n'
    beginCell = '<td>'
    endCell = '</td>'
    beginCellH = '<th>'
    endCellH = '</th>'
    beginLineGrey = '<tr style="background-color:LightGrey">'
    beginLineGreen = '<tr style="background-color:LightGreen">'
    beginLineBlue = '<tr style="background-color:LightBlue">'


    lumi = c.Params().getLumi(sel)


    yields_data = []
    yields_data_err = []

    Yields_bg = []
    Yields_sig1 = []
    Yields_sig2 = []
    

    for l in range(0,lTot):

        # Data yields:
        if data != None:
            dd = data.Get("Andrey/met0_et_"+str(l)).Clone()
            derr = Double(0)
            bins = dd.GetNbinsX()
            dy = dd.IntegralAndError(0, bins+1, derr)
            yields_data.append(dy)
            yields_data_err.append(derr)
        else:
            yields_data.append(0)
            yields_data_err.append(0)
            

        Yields_bg_per_cut = {}
        Yields_sig1_per_cut = {}
        Yields_sig2_per_cut = {}

        # Yields of backgrounds:
        #Handle single top backgrounds
        t  = []
        terr =[]
        
        #print "length of topbg list", len(topbg_list)
        for b in topbg_list:
            h = b.Get("Andrey/met0_et_"+str(l)).Clone()
            sample = b.Get("Andrey/evt_byCut").GetTitle()
            #print sample
            if mode=="scaled":
                handleOverflowBinsScaleAndColors(h, sample, lumi)
            err = Double(0)
            bins = h.GetNbinsX()
            y = h.IntegralAndError(0, bins+1, err)
            t.append(y)
            terr.append(err)

        ysum = 0
        toterr = 0
        for a in t:
            ysum += a
        for a in terr:
            toterr += a

        Yields_bg_per_cut["Top"] = ([ysum, toterr])

                                                                    
        # All the other backgrounds
        for b in bg_list:
            h = b.Get("Andrey/met0_et_"+str(l)).Clone()
            sample = b.Get("Andrey/evt_byCut").GetTitle()
            if mode=="scaled":
                handleOverflowBinsScaleAndColors(h, sample, lumi)
            err = Double(0)
            
            bins = h.GetNbinsX()
            y = h.IntegralAndError(0, bins+1, err)
            #print sample, y

            Yields_bg_per_cut[sample] = ([y,err])

        Yields_bg.append(Yields_bg_per_cut)

      
        # Yields of the signals:
        for s in sig1_list:
            h = s.Get("Andrey/met0_et_"+str(l)).Clone()
            sample = s.Get("Andrey/evt_byCut").GetTitle()
            #print l, sample
            if mode=="scaled":
                handleOverflowBinsScaleAndColors(h, sample, lumi)
            err = Double(0)
            
            bins = h.GetNbinsX()
            y = h.IntegralAndError(0, bins+1, err)
            #y = s.Get("Andrey/evt_byCut").GetBinContent(l+1)

            #print "integral = ", y
            Yields_sig1_per_cut[sample] = ([y,err])

        Yields_sig1.append(Yields_sig1_per_cut)

        # Yields of the signals:
        for s in sig2_list:
            h = s.Get("Andrey/met0_et_"+str(l)).Clone()
            sample = s.Get("Andrey/evt_byCut").GetTitle()
            #print l, sample
            if mode=="scaled":
                handleOverflowBinsScaleAndColors(h, sample, lumi)
            err = Double(0)
            
            bins = h.GetNbinsX()
            y = h.IntegralAndError(0, bins+1, err)
            #y = s.Get("Andrey/evt_byCut").GetBinContent(l+1)

            #print "integral = ", y
            Yields_sig2_per_cut[sample] = ([y,err])

        Yields_sig2.append(Yields_sig2_per_cut)

    #print Yields_bg
    myTable += beginTable
    
    myTable += beginLine
    myTable += beginCell
    myTable += endCell
    
    # Make the first line of the table
    for k, v in Yields_bg[0].iteritems():
        #print k, v
        myTable += beginCellH
        myTable += k
        myTable += endCellH
            
    myTable += beginCellH
    myTable += "Total bg"
    myTable += endCellH
    myTable += beginCellH
    myTable += "Data"
    myTable += endCellH
    myTable += beginCellH
    myTable += "ggHZZ 250"
    myTable += endCellH
    myTable += beginCellH
    myTable += "vbfHZZ 250"
    myTable += endCellH
    
    myTable += endLine
            
    for l in range(0, lMax):
        if l%2 == 0:
            myTable += beginLineGrey
        else:
            myTable += beginLine

        myTable += beginCellH
        myTable += cutNames[l]
        myTable += endCellH
        
        total_bg=0
        total_err=0
        for k, v in Yields_bg[l].iteritems():
            #print k, v
            myTable += beginCell
            myTable += '%d'% (v[0])
            myTable += endCell
            
            total_bg += v[0]
            total_err += v[1]
                
        myTable += beginCell
        myTable +='%d'% (total_bg)
        myTable += endCell
        myTable += beginCell
        myTable += '%d'% (yields_data[l])
        myTable += endCell
        myTable += beginCell
        myTable += '%d'% (Yields_sig1[l]["ggHZZ250"][0])
        myTable += endCell
        myTable += beginCell
        myTable += '%d'% (Yields_sig2[l]["VBFHZZ250"][0])
        myTable += endCell

        myTable += endLine


    for l in [25,26,27,10,18,19,20,21]:

        if l in [25,26,27]:
            myTable += beginLineGreen
        elif l in [10,18]:
            myTable += beginLineBlue
        else:
            myTable += beginLine
        
            
        myTable += beginCellH
        myTable += cutNames[l]
        myTable += endCellH
        
        total_bg=0
        total_err=0
        for k, v in Yields_bg[l].iteritems():
            #print k, v
            myTable += beginCell
            myTable += '%0.1f'% (v[0])
            #myTable += '%0.1f &pm; %0.1f'% (v[0],v[1])
            myTable += endCell
        
            total_bg += v[0]
            total_err += v[1]
                
        myTable += beginCell
        myTable +='%0.1f'% (total_bg)
        #myTable +='%0.1f &pm; %0.1f'% (total_bg, total_err)
        myTable += endCell
        myTable += beginCell
        myTable += '%d'% (yields_data[l])
        myTable += endCell
        myTable += beginCell
        myTable += '%0.1f'% (Yields_sig1[l]["ggHZZ250"][0])
        #myTable += '%0.1f &pm; %0.1f'% (Yields_sig1[l]["ggHZZ250"][0], Yields_sig1[l]["ggHZZ250"][1])
        myTable += endCell
        myTable += beginCell
        myTable += '%0.1f'% (Yields_sig2[l]["VBFHZZ250"][0])
        myTable += endCell
        myTable += endLine
        
    myTable += endTable

    hf = open(filename,"w")
    hf.write(myTable)
    hf.close()
    
def handleOverflowBinsScaleAndColors(hist, sample, lumi):
    xc = c.Params().xsec_and_colors()

    nBins   = hist.GetNbinsX()
    lastBin = hist.GetBinContent(nBins)
    ovflBin    = hist.GetBinContent(nBins+1);
    lastBinErr = hist.GetBinError(nBins);
    ovflBinErr = hist.GetBinError(nBins+1);
    firstBin    = hist.GetBinContent(1);
    undflBin    = hist.GetBinContent(0);
    firstBinErr = hist.GetBinError(1);
    undflBinErr = hist.GetBinError(0);
    hist.SetBinContent(nBins, lastBin+ovflBin);
    hist.SetBinError(nBins, sqrt(pow(lastBinErr,2) + pow(ovflBinErr,2)) );
    hist.SetBinContent(1, firstBin+undflBin);
    hist.SetBinError(1, sqrt(pow(firstBinErr,2) + pow(undflBinErr,2)) );
    hist.SetBinContent(0,0);
    hist.SetBinContent(nBins+1,0);

    if sample!="Data":
        if lumi!=0:  #don't rescale, for the cases we don't want it 
            hist.Scale(float(lumi*xc[sample][2])/xc[sample][3])
            #print "rescaling sample", sample, "xsection =", xc[sample][2], " total events: ", xc[sample][3] 
            hist.SetLineColor(xc[sample][0])
        if not "HZZ" in sample: 
            hist.SetFillColor(xc[sample][1])

    hist.SetLineWidth(2)
        
def makeStack(topbglist, otherbglist, histoname, lumi):

    hs = THStack("temp", "Stacked histo")
    #print " * In makeStack python function"

    xc = c.Params().xsec_and_colors()

    htw    = topbglist[0].Get("Andrey/"+histoname).Clone()
    htw.Scale(float(lumi*xc["tW"][2])/xc["tW"][3])
    htbarw = topbglist[1].Get("Andrey/"+histoname).Clone()
    htbarw.Scale(float(lumi*xc["tbarW"][2])/xc["tbarW"][3])
    handleOverflowBinsScaleAndColors(htw, "tW", lumi)
    handleOverflowBinsScaleAndColors(htbarw,"tbarW",lumi)    
    htw.Add(htbarw)
    hs.Add(htw)

    for f in otherbglist:
        #print f.GetName()
        evtHisto = f.Get("Andrey/evt_byCut")
        #evtHisto.Print()
        sample = evtHisto.GetTitle()
        #print "Double check the sample:",sample
        
        hh1 = f.Get("Andrey/"+histoname).Clone()
        handleOverflowBinsScaleAndColors(hh1,sample,lumi)
        hs.Add(hh1)

    return hs

def drawMultiPlot(fname,maintitle, xtitle, h_name, isLog, y1min, y1max, y2min, y2max, ovlist, topbg, bgList, sel):

    
    lumi = c.Params().getLumi(sel)
    
    name = "Andrey/"+h_name
    print "   * Plotting:", h_name
    print "   * Making a stack histo."
    
    hs = makeStack(topbg, bgList, h_name, lumi)
    #hs.Print()
    #hs = THStack()
    
    ff = []
    size = ovlist.GetSize()
    hh = [];
    if (size>4):
        print "To many plots to overlay"
        sys.exit(0)

    print "size of the list", size

    h_data = None
    sig = []
    #ovlist.At(0).Print()
    #ovlist.At(1).Print()
    for n in xrange(size):
        if ovlist.At(n) != None:
            ff.append(ovlist.At(n))
            hh.append(ff[n].Get(name).Clone())
            #print n, "file:", ff[n].GetName(), "   histoname: ", hh[n].GetName()
            
            sample = ff[n].Get("Andrey/evt_byCut").GetTitle()
            #print "Double check the sample:",sample
            
            if sample=="DATA":
                handleOverflowBinsScaleAndColors(hh[n], "Data", lumi)
                h_data = hh[0]
                h_data.SetMarkerStyle(20);
                h_data.SetMarkerSize(0.7)

            else:
                handleOverflowBinsScaleAndColors(hh[n], sample, lumi)
                sig.append(hh[n])

    doRatio = True
    if not isLog or h_data==None:
        doRatio = False

    cc =None
    if doRatio:
        cc = TCanvas("c2","big canvas",600,700);
    else:
        cc = TCanvas("c1","small canvas",600,600);
        cc.SetLogy(isLog)
                
    cc.cd()
    
    #print "drawing"
    if doRatio:
        pad1 = TPad("pad1","pad1",0,0.3,1,1);
        pad1.SetBottomMargin(0);
        pad1.Draw();
        pad1.cd();
        pad1.SetLogy(isLog)
    
    hs.Draw("hist")
    
    if doRatio:
        hs.SetTitle(maintitle+";; Events")
    else:
        hs.SetTitle(maintitle+";"+xtitle+"; Data/MC");
    hs.SetMinimum(y1min);
    hs.SetMaximum(y1max);

    if h_data!=None:
        h_data.Draw("same e1pl");
    for s in sig:
        s.Draw("same hist");
                
    if doRatio:
        cc.cd()
        pad2 = TPad("pad2","pad2",0,0,1,0.3);
        pad2.SetBottomMargin(0.25);
        pad2.SetTopMargin(0);
        pad2.Draw();
        pad2.cd();
  
        h1 = h_data.Clone("")
        stackSum = hs.GetStack().Last().Clone()
    
        h1.Divide(stackSum);
        h1.SetTitle(";"+xtitle+"; Data/MC");
        h1.SetMaximum(y2max);
        h1.SetMinimum(y2min);
        h1.GetYaxis().SetNdivisions(206);
        h1.GetYaxis().SetTitleOffset(0.4);
        h1.SetTitleSize(0.1,"XYZ");
        h1.SetLabelSize(0.1,"XY");

        h1.Draw("e1p");
    
        pad1.cd();

    prelim = TLatex(0.25,0.95, "CMS Preliminary       #it{L_{int}} = %0.1f fb^{-1}" % (lumi/1000.))
    prelim.SetNDC();
    prelim.SetTextSize(0.03); 
    prelim.Draw();


    selection = TLatex()
    if sel==1:
        selection = TLatex(0.70,0.95, "#mu#mu selection");
        selection.SetTextColor(kBlue-3);
    if sel==2:
        selection  = TLatex(0.70,0.95, "#font[32]{ee} selection");
        selection.SetTextColor(kGreen-3);
  
    selection.SetNDC();
    prelim.SetTextSize(0.03); 
    selection.Draw();
        

    leg01 = TLegend(0.53,0.7,0.95,0.90);
    leg01.SetNColumns(2);
    leg01.SetTextSize(0.04)
    leg01.AddEntry(hs.GetStack()[0],"top","f")
    leg01.AddEntry(hs.GetStack()[1],"WW","f")
    leg01.AddEntry(hs.GetStack()[2],"WZ","f")
    leg01.AddEntry(hs.GetStack()[3],"ZZ","f")
    leg01.AddEntry(hs.GetStack()[4],"ttbar","f")
    leg01.AddEntry(hs.GetStack()[5],"Z + jets","f")
    if h_data != None:
        leg01.AddEntry(h_data,"Data","epl")
    for s in sig:
        leg01.AddEntry(s,"H250","l")

    leg01.SetFillColor(kWhite)
    leg01.Draw()

                                            
    cc.SaveAs(fname+".png")
