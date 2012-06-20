#import sys,os
from ROOT import *

class SampleParams():
    def __init__(self):
        
        self._xsec_and_colors = {
            "ZZ":[kBlue, kMagenta, 0.1787, 1103468],
            "WZ":[kBlue, kMagenta, 0.1787, 1103468],
            "ggHZZ200":[kPink, -1, 0.05406, 96.20e3],
            "ggHZZ600":[kPink, -1, 0.05406, 96.20e3],
            }
        
    def xsec_and_colors(self):
        return  self._xsec_and_colors

