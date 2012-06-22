#import sys,os
from ROOT import *

class SampleParams():
    def __init__(self):
        
        self._xsec_and_colors = {
            "ZZ":[kBlue, kMagenta, 0.1787, 1103468],
            "WZ":[kBlue, kMagenta, 0.1787, 1103468],
            "ggHZZ200":[kPink, -1, 0.05406, 96.215e3],
            "ggHZZ250":[kPink, -1, 0.0398,  99.993e3],
            "ggHZZ300":[kPink, -1, 0.0301,  99.990e3],
            "ggHZZ350":[kPink, -1, 0.0286,  99.998e3],
            "ggHZZ400":[kPink, -1, 0.0221,  93.854e3],
            "ggHZZ450":[kPink, -1, 0.0143,  99.989e3],
            "ggHZZ500":[kPink, -1, 0.00896, 99.980e3],
            "ggHZZ550":[kPink, -1, 0.00568, 97.824e3],
            "ggHZZ600":[kPink, -1, 0.00359, 99.972e3],
            "ttbar":[kOrange+1, kGreen+2, 16.71, 10.34e6],
            }
        
    def xsec_and_colors(self):
        return  self._xsec_and_colors

