#!/usr/bin/env python
import sys,os

from ROOT import *
import config as c


myParams = c.Params("8TeV")

myAnaParams =  myParams.analysisParams("HZZ")
cutNames    = myAnaParams[0]
nominalCuts = myAnaParams[1]
extraCuts   = myAnaParams[2]
bgOrder     = myAnaParams[3]

def updateNevents(bg_list, sig1_list, sig2_list, sig3_list):

    if not myParams.isUpdated():
        allsamples = dict(bg_list, **sig1_list)
        allsamples.update(sig2_list)
        allsamples.update(sig3_list)
        #print allsamples
        for a in allsamples:
            #print a, myParams.xsec_and_colors()[a]
            h = allsamples[a].Get("Histos/met0_et_0")
            nn =  h.GetBinContent(1)
            #print "N events ", nn
            
            myParams.UpdateNevents(a, nn)
            #print myParams.xsec_and_colors()[a]
        myParams.SetIsUpdated(True)
        print "  --> End of updating Nevts in config file: updateNevents() function"
        return myParams.xsec_and_colors()
    else:
        return myParams.xsec_and_colors()

    
def makeMvaTrees(outpath, bg_list, sig1_list, sig2_list, sel):
    print "\n  ** Creating mva trees ***"

    updateNevents(bg_list, sig1_list, sig2_list, sig3_list)

    fbg_train = TFile(outpath+"allBg_train.root","RECREATE")
    fbg_test  = TFile(outpath+"allBg_test.root","RECREATE")

    gROOT.ProcessLine(
"struct mva_vars {\
Float_t lep1Pt;\
Float_t lep2Pt;\
Float_t diLepPt;\
Float_t diLepM;\
Float_t diLepEta;\
Float_t lepDeltaPhi;\
Float_t lepDeltaEta;\
Float_t lepDeltaR;\
Float_t lepAngle;\
Float_t lepPtRatio;\
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

    allVars = 'lep1Pt/F:lep2Pt:diLepPt:diLepM:diLepEta:lepDeltaPhi:lepDeltaEta:lepDeltaR:lepAngle:lepPtRatio:\
met:metOverQt:metProjOnQt:metPerpQt:\
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
        sample = b.Get("Histos/evt_byCut").GetTitle()
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
            sample = s.Get("Histos/evt_byCut").GetTitle()
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
    lumi = myParams.getLumi(sel)

    #xc = myParams.xsec_and_colors()

    stuff.lep1Pt   = evt.lep1Pt
    stuff.lep2Pt   = evt.lep2Pt
    stuff.diLepPt  = evt.diLepPt
    stuff.diLepM   = evt.diLepM
    stuff.diLepEta = evt.diLepEta
    stuff.lepDeltaPhi = evt.lepDeltaPhi
    stuff.lepDeltaEta = evt.lepDeltaEta
    stuff.lepDeltaR   = evt.lepDeltaR
    stuff.lepAngle    = evt.lepAngle
    stuff.lepPtRatio  = evt.lepPtRatio
    stuff.met         = evt.met
    stuff.metOverQt    =  evt.metOverQt
    stuff.metProjOnQt  =  evt.metProjOnQt
    stuff.metPerpQt    = evt.metPerpQt
    stuff.dPhiMetDiLep = evt.dPhiMetDiLep
    stuff.dPhiJetMet   = evt.dPhiJetMet  
    stuff.mt      = evt.mt
    stuff.nJets   = evt.nJets
    stuff.nJets15 = evt.nJets15
    stuff.deltaEtaDiJet  = evt.deltaEtaDiJet
    stuff.massDiJet      = evt.massDiJet
    stuff.zeppDiJetDiLep = evt.zeppDiJetDiLep
    stuff.evWeight   = evt.weight
    stuff.fullWeight = evt.weight* float(lumi*myParams.getCS(sample))/myParams.getNev(sample)


def calcYields(bg_list, sig1_list, sig2_list, sig3_list, data, sel, mode="scaled"):

    updateNevents(bg_list, sig1_list, sig2_list, sig3_list)

    lumi = myParams.getLumi(sel)

    lTot = 33   # Total number of cuts


    Yields_data     = []
    Yields_data_err = []

    Yields_bg   = []
    Yields_sig1 = []
    Yields_sig2 = []
    

    for l in range(0,lTot):
        # Data yields:
        if data != None:
            dd = data.Get("Histos/evt_byCut").Clone()
            derr = Double(0)

            if l!=0:
                dy = dd.GetBinContent(l+1)
            else:
                dy = data.Get("Histos/met0_et_0").GetBinContent(1)
                
            Yields_data.append(dy)
            Yields_data_err.append(derr)
        else:
            Yields_data.append(0)
            Yields_data_err.append(0)
            

        Yields_bg_per_cut   = {}
        Yields_sig1_per_cut = {}
        Yields_sig2_per_cut = {}

        # Yields of backgrounds:
        #Handle single top backgrounds
        t  = []
        terr =[]

        top_bg = []
        
        #for b in bg_list:
        for b in bgOrder:

            #print "yields making", b

            if b not in myParams.toMerge().keys():
                h = bg_list[b].Get("Histos/evt_byCut").Clone()
                sample = b
                if mode=="scaled" and l!=0:
                    handleOverflowBinsScaleAndColors(h, b, lumi)
                    #handleOverflowBinsScaleAndColors(h, sample, lumi)
                    y = h.GetBinContent(l+1)
                else:
                    if l!=0:
                        #print sample, mode, l
                        y =  bg_list[b].Get("Histos/evt_byCut_raw").GetBinContent(l+1)
                    else:
                        y = myParams.getNev(b)

                err = Double(0)

                Yields_bg_per_cut[sample] = ([y,err])
            else:
                #print "in to merge dict! ", b
                Yields_bg_per_cut[b] = ([0,0])

                samplesToMerge = myParams.toMerge()[b]
                #print samplesToMerge
                for m in samplesToMerge:
                    hm = bg_list[m].Get("Histos/evt_byCut").Clone()
                    #print m
                    if mode=="scaled" and l!=0:
                        handleOverflowBinsScaleAndColors(hm, m, lumi)
                        #handleOverflowBinsScaleAndColors(h, sample, lumi)
                        y = hm.GetBinContent(l+1)
                        #print l+1,y
                    else:
                        if l!=0:
                            #print m, mode, l
                            y =  bg_list[m].Get("Histos/evt_byCut_raw").GetBinContent(l+1)
                        else:
                            y = myParams.getNev(m)
                    #print "y=",y
                    Yields_bg_per_cut[b][0] += y
                    Yields_bg_per_cut[b][1] = 0


        Yields_bg.append(Yields_bg_per_cut)

        # Yields of the signals:
        for s in sig1_list:
            #print l, s

            h = sig1_list[s].Get("Histos/evt_byCut").Clone()
            sample = s
            if mode=="scaled" and l!=0:
                handleOverflowBinsScaleAndColors(h, sample, lumi)
            y = h.GetBinContent(l+1)
            err = Double(0)

            #print "integral = ", y
            Yields_sig1_per_cut[sample] = ([y,err])

        Yields_sig1.append(Yields_sig1_per_cut)

        # Yields of the signals:
        for s in sig2_list:
            h = sig2_list[s].Get("Histos/evt_byCut").Clone()
            sample = s
            if mode=="scaled" and l!=0:
                handleOverflowBinsScaleAndColors(h, sample, lumi)
            y = h.GetBinContent(l+1)
            err = Double(0)
            
            Yields_sig2_per_cut[sample] = ([y,err])

        Yields_sig2.append(Yields_sig2_per_cut)

    #print Yields_bg
    return [Yields_data, Yields_bg, Yields_sig1, Yields_sig2]

def printYields(bg_list, sig1_list, sig2_list, sig3_list, data, sel, filename, anaType="HZZ250", mode="scaled"):

    y = calcYields(bg_list, sig1_list, sig2_list, sig3_list, data, sel, mode)

    Yields_data = y[0]
    Yields_bg   = y[1]
    Yields_sig1 = y[2]
    Yields_sig2 = y[3]

    print "Yield are HTML mode"
    favMass="0"
    if "VBFZ" not in anaType:
        favMass = anaType[3:]
        anaType ="HZZ"

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

    myTable += beginTable
    
    myTable += beginLine
    myTable += beginCell
    myTable += endCell
    
    # Make the first line of the table

    for b in bgOrder:
        if b not in Yields_bg[0].keys(): continue
        #v = Yields_bg[0][b]
        myTable += beginCellH
        myTable += b[0:5]
        myTable += endCellH
            
    myTable += beginCellH
    myTable += "Total bg"
    myTable += endCellH
    myTable += beginCellH
    myTable += "Data"
    myTable += endCellH

    if "VBFZ" not in anaType:
        myTable += beginCellH
        myTable += "ggHZZ "
        myTable += endCellH
        myTable += beginCellH
        myTable += "vbfHZZ "
        myTable += endCellH
    
    myTable += endLine
    
    for l in nominalCuts:
        if l%2 == 0:
            myTable += beginLineGrey
        else:
            myTable += beginLine

        myTable += beginCellH
        myTable += cutNames[l]
        myTable += endCellH
        
        total_bg=0
        total_err=0

        #print l, Yields_bg[l]
        
        for b in bgOrder:
            if b not in Yields_bg[l].keys(): continue
            v = Yields_bg[l][b]
            #print v
            myTable += beginCell
            myTable += '%d'% (v[0])
            myTable += endCell
            
            total_bg += v[0]
            total_err += v[1]
                
        myTable += beginCell
        myTable +='%d'% (total_bg)
        myTable += endCell
        myTable += beginCell
        myTable += '%d'% (Yields_data[l])
        myTable += endCell
        if "VBFZ" not in anaType:
            myTable += beginCell
            #myTable += ''
            myTable += endCell
            myTable += beginCell
            #myTable += '%d'% (Yields_sig2[l]["VBFHZZ"+favMass][0])
            myTable += endCell
            
        myTable += endLine


    for l in extraCuts:

        if l in [18,19,20]:
            myTable += beginLineGreen
        elif l in [9,10,29,25,21]:
            myTable += beginLineBlue
        else:
            myTable += beginLine
        
            
        myTable += beginCellH
        myTable += cutNames[l]
        myTable += endCellH
        
        total_bg=0
        total_err=0
        for b in bgOrder:
            if b not in Yields_bg[l].keys(): continue
            v = Yields_bg[l][b]
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
        myTable += '%d'% (Yields_data[l])
        myTable += endCell
        if "VBFZ" not in anaType:
            myTable += beginCell
            if l in [21,22,23,24]:
                hMass = "125"
            elif l in [9,25,26,27,28]:
                hMass = "200"
            elif l in [10,29,30,31,32]:
                hMass = "300"
            else: hMass=favMass
            myTable += '%0.2f'% (Yields_sig1[l]["ggHZZ"+hMass][0])
            #myTable += '%0.1f &pm; %0.1f'% (Yields_sig1[l]["ggHZZ250"][0], Yields_sig1[l]["ggHZZ250"][1])
            myTable += endCell
            myTable += beginCell
            myTable += '%0.2f'% (Yields_sig2[l]["VBFHZZ"+hMass][0])
            myTable += endCell
        myTable += endLine
        
    myTable += endTable

    hf = open(filename,"w")
    hf.write(myTable)
    hf.close()
    
def handleOverflowBinsScaleAndColors(hist, sample, lumi):
    if hist == None:
        return 
    xc = myParams.xsec_and_colors()

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
            #if sample=="ttbar":
            #print "rescaling sample", sample, "xsection =", myParams.getCS(sample), " total events: ", myParams.getNev(sample) 
            #    print "lumi=", lumi
                
            hist.Scale(float(lumi*myParams.getCS(sample))/myParams.getNev(sample))
            hist.SetLineColor(xc[sample][0])
        if not "HZZ" in sample: 
            hist.SetFillStyle(xc[sample][1])
            hist.SetFillColor(xc[sample][2])


    hist.SetLineWidth(2)
        
def makeStack(bgList, histoname, lumi):

    #xc = myParams.xsec_and_colors()

    hs = THStack("temp", "Stacked histo")
    #print " * In makeStack python function"

    bgOrder = myParams.analysisParams("HZZ")[3]

    leg01 = TLegend(0.53,0.7,0.95,0.90);
    leg01.SetNColumns(2);
    leg01.SetTextSize(0.04)


    samplesToMerge = myParams.toMerge()

    for b in bgOrder:
        dirname = "Histos/"
        if histoname[0:2] == "mu":
            dirname = "Muons/"
        if histoname[0:2] == "el":
            dirname = "Electrons/"

        if b in samplesToMerge.keys():
            hm = TH1F()
            for i,m in enumerate(samplesToMerge[b]):
                #print "stack for ", b, i,m, samplesToMerge[b]
                if bgList[m].Get(dirname+histoname) != None:
                    h = bgList[m].Get(dirname+histoname).Clone()
                    handleOverflowBinsScaleAndColors(h, m, lumi)
                    if i==0:
                        hm = h 
                    else:
                        hm.Add(h)

            if hm!= None:
                hs.Add(hm)
                leg01.AddEntry(hm,b,"f")

        elif b not in bgList.keys(): continue
        else:
            #print bgList[b].GetName()
            evtHisto = bgList[b].Get("Histos/evt_byCut")
            evtHisto.Print()
            sample = evtHisto.GetTitle()
            #print "Double check the sample:",sample
        
            hh1 = bgList[b].Get(dirname+histoname).Clone()
            handleOverflowBinsScaleAndColors(hh1,b,lumi)
            #handleOverflowBinsScaleAndColors(hh1,sample,lumi)
            hs.Add(hh1)

            leg01.AddEntry(hh1,b[0:2],"f")

    return [hs, leg01]

def drawMultiPlot(fname,maintitle, xtitle, h_name, isLog, y1min, y1max, y2min, y2max, ovList, bgList, sel, doRatio):

    lumi = myParams.getLumi(sel)
    #print lumi
    bgOrder = myParams.analysisParams("HZZ")[3]
    
    dirname = "Histos/"
    if h_name[0:2] == "mu":
        dirname = "Muons/"
    if h_name[0:2] == "el":
        dirname = "Electrons/"
        
    name = dirname+h_name
    print "   * Plotting:", h_name
    print "   * Making a stack histo."
    
    ss = makeStack(bgList, h_name, lumi)
    hs = ss[0]
    leg01 = ss[1]

    #hs.Print()
    #hs = THStack()
    
    ff = []
    size = len(ovList)
    hh = [];
    if (size>4):
        print "To many plots to overlay"
        sys.exit(0)

    #print "size of the overlay list", size

    h_data = None
    sig = []
    n=0
    for o in ovList:
        #print o
        ff.append(ovList[o])
        if ff[n].Get(name) !=None:
            hh.append(ff[n].Get(name).Clone())
            #print n, "file:", ff[n].GetName(), "   histoname: ", hh[n].GetName()
        else:
            hh.append(None)            
        sample = o#ff[n].Get("Histos/evt_byCut").GetTitle()
        #print "Double check the sample:",sample
            
        if sample=="DATA":
            handleOverflowBinsScaleAndColors(hh[n], "Data", lumi)
            h_data = hh[n]
            h_data.SetMarkerStyle(20);
            h_data.SetMarkerSize(0.7)

        else:
            if hh[n]!=None:
                handleOverflowBinsScaleAndColors(hh[n], sample, lumi)
                sig.append(hh[n])
        n+=1
        
    if h_data==None:
        doRatio = False

    cc =None
    if doRatio:
        cc = TCanvas("c2","big canvas",600,700);
    else:
        cc = TCanvas("c3","small canvas",600,600);
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
    if hs.GetHistogram().GetXaxis().GetXmax()>1000:
        hs.GetHistogram().SetNdivisions(505)
    
    if doRatio:
        hs.SetTitle(maintitle+";; Events")
    else:
        hs.SetTitle(maintitle+";"+xtitle+"; Events");
    hs.SetMinimum(y1min);
    hs.SetMaximum(y1max);

    if h_data!=None:
        h_data.Draw("same e1pl");
    for i,s in enumerate(sig):
        #print i,s
        if i==0:
            s.Scale(5)  #for ggH
        elif i==1:
            s.Scale(5)  #for vbf
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

    prelim = TLatex(0.15,0.95, "CMS Preliminary  %s    #it{L_{int}} = %0.1f fb^{-1}" % (myParams.getS(), lumi/1000.))
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
        
    if h_data != None:
        leg01.AddEntry(h_data,"Data","epl")
        
    for i,s in enumerate(sig):
        #print i,s
        if i==0:
            leg01.AddEntry(sig[1],"5x ggH 200","l")
        elif i==1:
            leg01.AddEntry(sig[0],"5x vbfH 200","l")

    leg01.SetFillColor(kWhite)
    leg01.Draw()

                                            
    cc.SaveAs(fname+"_"+h_name+".png")
