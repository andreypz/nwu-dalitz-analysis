#!/usr/bin/env python
import sys,os

from ROOT import *

effic = {
    "ee":  [0.8469, 0.0369],
    "mumu":[0.8541, 0.0355],
    } 
accep = {
    "ee": {
        "A2tt": [0.001363, 0.000073],
        "A1tt": [0.001254, 0.000011],
        "A0tt": [0.000255, 0.000005],
        "A2t":  [0.000054, 0.000006],
        "A1t":  [0.000173, 0.000010],
        "A0t":  [0.000061, 0.000006],
        },
    "mumu": {
        "A2tt": [0.001956, 0.000087],
        "A1tt": [0.001798, 0.000013],
        "A0tt": [0.000385, 0.000006],
        "A2t":  [0.000057, 0.000006],
        "A1t":  [0.000215, 0.000008],
        "A0t":  [0.000059, 0.000004],
        }
    } 


def topbg(h_ttbar, h_tW, h_tbarW, sel):
    # Input hystograms are the b-jet multiplicity hists
    N0 = N1 = N2 = 0
    if sel ==1: selstr = "mumu"
    else: selstr = "ee"
    if h_ttbar!=None:
        N0 = h_ttbar.GetBinContent(1)
        N1 = h_ttbar.GetBinContent(2)
        N2 = h_ttbar.GetBinContent(3)
    
    print "ttbar Numbers:\n\
    N2 = %.2f\n\
    N1 = %.2f\n\
    N0 = %.2f\n\
    "  % (N2,N1,N0)

    A2tt = accep[selstr]["A2tt"][0]
    A1tt = accep[selstr]["A1tt"][0]
    A0tt = accep[selstr]["A0tt"][0]
 
    eps = N2*(2*A2tt+A1tt)/(N1+2*N2)/A2tt
    Ntt = N2/(eps**2*A2tt)

    N0_calc = Ntt*(A2tt*(1-eps)**2 + A1tt*(1-eps) + A0tt) 
    
    print "Calculated:\n\
    Ntt = %.2f\n\
    eff = %.2f\n\
    N0 = %.2f\n\
    " % (Ntt, eps, N0_calc)
