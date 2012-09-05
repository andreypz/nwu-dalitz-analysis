#import sys,os
from ROOT import *

class Params():
    def __init__(self):

        self._lumi_ele = 2312 + 2739
        self._lumi_mu  = 2312 + 2739
        #self._lumi_mu  = 30000

        self._cutNamesHZZ = [
        "0. Total",
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
        "18. mvaTree, nj=0",
        "19. mvaTree, nj=1",
        "20. mvaTree, nj>1",
        "21. MVA H",
        "22. 0j bin",
        "23. 1j bin",
        "24. >1j bin",
        "25. MVA H250",
        "26. 0j",
        "27. 1j",
        "28. >1j",
        ]
        self._nominalHZZ = [1,2,3,4,5,6,7,8]
        self._extraHZZ = [10,25,26,27,28,18,19,20]

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
        
        self._xsec_and_colors = {
            "DATA": [-1,-1,-1,-1],
            #"ZZ":   [kBlue,     kMagenta, 6.46, 4191.045e3],  #pythia
            "ZZ":   [kBlue,     kMagenta, 0.1787, 1103.468e3],  #madgraph
            "WZ":   [kYellow-1,    kBlue-1,  0.856, 1221.134e3], #or 0.879?
            "WW":   [kOrange+1, kGreen+2,  4.783, 1197.558e3],

            "DYjets":[kRed,  kWhite,  3048.0, 36.278e6],
            "Wjets": [kGreen-2,  kCyan+1,  31314.0, 81.251e6],
            
            #"ttbar":[kOrange+1, kGreen+2, 16.71, 10.34e6],
            "ttbar": [kOrange+1, kSpring-9, 165, 3701.947e3], #madgraph
            "tW":    [kOrange+9, kOrange+6, 7.87, 814.390e3],
            "tbarW": [kRed,      kWhite,   7.87, 809.984e3],

            "vbfZ": [kBlue+3,  -1,   0.248/2, 1732.824e3],

            "ggHZZ125":[kCyan+1, -1, 0.01651, 292.174e3],
            "ggHZZ200":[kCyan+1, -1, 0.05425, 96.215e3],
            "ggHZZ250":[kCyan+1, -1, 0.0398,  99.993e3],
            "ggHZZ300":[kCyan+1, -1, 0.0301,  99.990e3],
            "ggHZZ350":[kCyan+1, -1, 0.0286,  99.998e3],
            "ggHZZ400":[kCyan+1, -1, 0.0221,  93.854e3],
            "ggHZZ450":[kCyan+1, -1, 0.0143,  99.989e3],
            "ggHZZ500":[kCyan+1, -1, 0.00896, 99.980e3],
            "ggHZZ550":[kCyan+1, -1, 0.00568, 97.824e3],
            "ggHZZ600":[kCyan+1, -1, 0.00359, 99.972e3],

            "ggHWW125":[kCyan, -1, 0.34723, 0e3],
            "ggHWW200":[kCyan, -1, 0.40840, 109.993e3],
            "ggHWW250":[kCyan, -1, 0.24378, 106.038e3],
            "ggHWW300":[kCyan, -1, 0.17598, 109.989e3],
            "ggHWW350":[kCyan, -1, 0.16368, 109.991e3],
            "ggHWW400":[kCyan, -1, 0.12418, 109.993e3],
            "ggHWW450":[kCyan, -1, 0.07880, 109.977e3],
            "ggHWW500":[kCyan, -1, 0.04868, 109.973e3],
            "ggHWW550":[kCyan, -1, 0.03053, 109.973e3],
            "ggHWW600":[kCyan, -1, 0.01914, 109.972e3],

            "VBFHZZ125":[kMagenta+2, -1, 1.310e-3, 49.347e3],
            "VBFHZZ200":[kCyan+1, -1, 6.562e-3, 49.943e3],
            "VBFHZZ250":[kCyan+1, -1, 5.163e-3, 49.957e3],
            "VBFHZZ300":[kCyan+1, -1, 3.734e-3, 49.935e3],
            "VBFHZZ350":[kCyan+1, -1, 2.644e-3, 46.002e3],
            "VBFHZZ400":[kCyan+1, -1, 1.760e-3, 49.851e3],
            "VBFHZZ450":[kCyan+1, -1, 1.299e-3, 46.065e3],
            "VBFHZZ500":[kCyan+1, -1, 1.001e-3, 46.049e3],
            "VBFHZZ550":[kCyan+1, -1, 0.876e-3, 49.729e3],
            "VBFHZZ600":[kCyan+1, -1, 0.634e-3, 49.698e3],
            }
        
    def xsec_and_colors(self):
        return  self._xsec_and_colors

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
            return [self._cutNamesHZZ,  self._nominalHZZ, self._extraHZZ]
        elif anatype=="VBFZ":
            return [self._cutNamesVBFZ,  self._nominalVBFZ, self._extraVBFZ]
        else:
            return []
