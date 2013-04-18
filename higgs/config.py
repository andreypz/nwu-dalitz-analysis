#import sys,os
from ROOT import *

class Params():
    def __init__(self, s="7TeV"):
        self._s=s

        if self._s=="7TeV":
            self._lumi_ele = 2312 + 2739
            self._lumi_ele = 2312 +  (277./312)*2739
            self._lumi_mu  = 2312 + 2739
            self._lumi_mu  = 2312 +  (137./144)*2739
            #self._lumi_mu  = 2312 +  (21./30)*2739
        else:
            self._lumi_ele = 19600
            self._lumi_mu  = 19600
            
            
        self._isUpdated = False
        self._cutNamesHZZ = [
        "0. Total (not scaled)",
        "1. HLT",
        "2. Two Leptons",
        "3. third lepton veto",
        "4. Z peak",
        "5. Met >70",
        "6. dPhi cut",
        "7. Z pt cut",
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
        "18. mvaTree, nj=0",
        "19. mvaTree, nj=1",
        "20. mvaTree, nj>1",
        "21. MVA H125",
        "22. MVA H125 nj=0",
        "23. MVA H125 nj=1",
        "24. MVA H125 nj>1",
        "25. MVA H200",
        "26. MVA H200 nj=0",
        "27. MVA H200 nj=1",
        "28. MVA H200 nj>1",
        "29. MVA H250",
        "30. MVA H250 nj=0",
        "31. MVA H250 nj=1",
        "32. MVA H250 nj>1",
    
        ]
        self._nominalHZZ = [0,1,2,3,4,5,6,7,8]
        #self._bgOrderHZZ = ["Top","ttbar","WW","WZ","ZZ","DYjets"]        
        self._bgOrderHZZ = ["Top","ttbar","DYjets"]        
        #self._bgOrderHZZ = ["Top","ttbar","WWJetsTo2L2Nu","WZJetsTo3LNu","ZZJetsTo2L2Nu","DYjets"]        
        self._extraHZZ = [18,19,20, 9,13]
        #self._extraHZZ = [18,19,20, 9,25,26,27,28, 10,29,30,31,32, 21,22,23,24]
        self._toMerge ={"Top":["tW","tbarW"],
                        "DYjets":["DYjets","DYjets10"]}
                        #"DYjets":["DYjets10"]}
        
        self._cutNamesVBFZ = [
        "0. Total",
        "1. HLT",
        "2. Two Leptons",
        "3. Z peak",
        "4. pt j1>65, j2>40; |eta|<3.6",
        "5. |y*| < 1.2",
        "6. M(j1,j2) > 600",
        "7. Third lepton veto",
        "8. ",
        ]
        self._nominalVBFZ = [3,4,5,6,7]
        self._extraVBFZ = []
        self._bgOrderVBFZ = ["Top","ttbar","WW","ZZ","DYjets"]        
        
        self._xsec_and_colors = {
            "DATA": [-1,-1,-1,-1,-1,-1],
            # schenme: lineColor, fill style, fill color, cs@7, cs@8, Nevts
            #"ZZ":   [kBlue,     3004, kMagenta, 6.46, ,0],  #pythia
            "ZZJetsTo2L2Nu":   [kBlue,    3244, kMagenta, 0.1787, 0.32,  0],  #madgraph
            "WZJetsTo3LNu":    [kYellow-1,1001, kBlue-1,   0.856, 0.75,  0], #or 0.879?
            "WWJetsTo2L2Nu":   [kOrange+1,1001, kGreen+2,  4.783, 6.01,  0],

            "DYjets":[kRed,     3004, kRed+1,   3048.0,   3503.71, 0],
            "DYjets10":[kRed,   3004, kRed+1,   0,   860, 0],
            "Wjets": [kGreen-2, 3004, kCyan+1,  31314.0, 12085.73, 0],
            
            #"ttbar":[kOrange+1,3004, kGreen+2, 16.71, 000, 0], 
            "ttbar": [kOrange+1,1001, kSpring-9, 165, 227.197,  0], #madgraph
            "tW":    [kOrange+9,1001, kOrange+6, 7.87, 11.18, 0],
            "tbarW": [kRed,     1001, kWhite,    7.87, 11.18, 0],

            "vbfZ": [kBlue+3,  3004, -1,   0.248/2, 0.248/2,0],

            "ggHZZ125":[kCyan+1, -1,-1, 0.01651, 0.02396,0],
            "ggHZZ200":[kCyan+1, -1,-1, 0.05425, 0.06876,0],
            "ggHZZ250":[kCyan+1, -1,-1, 0.0398,  0.05760,0],
            "ggHZZ300":[kCyan+1, -1,-1, 0.0301,  0.04471,0],
            "ggHZZ350":[kCyan+1, -1,-1, 0.0286,  0.04223,0],
            "ggHZZ400":[kCyan+1, -1,-1, 0.0221,  0.03177,0],
            "ggHZZ450":[kCyan+1, -1,-1, 0.0143,  0.02095,0],
            "ggHZZ500":[kCyan+1, -1,-1, 0.00896, 0.01352,0],
            "ggHZZ550":[kCyan+1, -1,-1, 0.00568, 0.00878,0],
            "ggHZZ600":[kCyan+1, -1,-1, 0.00359, 0.00448,0],

            "ggHWW125":[kCyan, -1,-1, 0.34723, 1,0],
            "ggHWW200":[kCyan, -1,-1, 0.40840, 1,0],
            "ggHWW250":[kCyan, -1,-1, 0.24378, 1,0],
            "ggHWW300":[kCyan, -1,-1, 0.17598, 1,0],
            "ggHWW350":[kCyan, -1,-1, 0.16368, 1,0],
            "ggHWW400":[kCyan, -1,-1, 0.12418, 1,0],
            "ggHWW450":[kCyan, -1,-1, 0.07880, 1,0],
            "ggHWW500":[kCyan, -1,-1, 0.04868, 1,0],
            "ggHWW550":[kCyan, -1,-1, 0.03053, 1,0],
            "ggHWW600":[kCyan, -1,-1, 0.01914, 1,0],

            "VBFHZZ125":[kMagenta+2, -1,-1, 1.310e-3, 1.715e-3, 0],
            "VBFHZZ200":[kMagenta+2, -1,-1, 6.562e-3, 8.728e-3, 0],
            "VBFHZZ250":[kMagenta+2, -1,-1, 5.163e-3, 6.960e-3, 0],
            "VBFHZZ300":[kMagenta+2, -1,-1, 3.734e-3, 5.466e-3, 0],
            "VBFHZZ350":[kMagenta+2, -1,-1, 2.644e-3, 3.968e-3, 0],
            "VBFHZZ400":[kMagenta+2, -1,-1, 1.760e-3, 2.763e-3, 0],
            "VBFHZZ450":[kMagenta+2, -1,-1, 1.299e-3, 1.832e-3, 0],
            "VBFHZZ500":[kMagenta+2, -1,-1, 1.001e-3, 1.646e-3, 0],
            "VBFHZZ550":[kMagenta+2, -1,-1, 0.876e-3, 1.119e-3, 0],
            "VBFHZZ600":[kMagenta+2, -1,-1, 0.634e-3, 1.064e-3, 0],
            }
        
    def xsec_and_colors(self):
        return  self._xsec_and_colors

    def SetIsUpdated(self, bit=True):
        self._isUpdated = bit

    def isUpdated(self):
        return self._isUpdated

    def UpdateNevents(self, sample, n=1):
        if  self._xsec_and_colors[sample][5] != n:
            #if sample in ["ttbar","tW","tbarW", "DATA"]:
            #    print " --->  Changing the initial number of events in config file for sample", sample,"to", n
            self._xsec_and_colors[sample][5] = n 
        if sample not in self._xsec_and_colors.keys():
            print "Warning!! the sample name was no in the original dictionary"
            
    def getNev(self, sample):
        return self._xsec_and_colors[sample][5]
    
    def getCS(self, sample):
        if self._s == "7TeV": 
            return self._xsec_and_colors[sample][3]
        else: 
            return self._xsec_and_colors[sample][4]

    def getS(self):
        return self._s

    def toMerge(self):
        return self._toMerge

    def getLumi(self, sel=1):
        if sel==1:
            return self._lumi_mu
        elif sel==2:
            return self._lumi_ele
        else:
            print "No lumi for selection", sel
            return 0

    def analysisParams(self,anatype="HZZ"):
        if anatype=="HZZ":
            return [self._cutNamesHZZ,  self._nominalHZZ, self._extraHZZ, self._bgOrderHZZ]
        elif anatype=="VBFZ":
            return [self._cutNamesVBFZ,  self._nominalVBFZ, self._extraVBFZ,self._bgOrderVBFZ]
        else:
            return []
