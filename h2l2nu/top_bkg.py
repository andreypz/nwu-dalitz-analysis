#!/usr/bin/env python
import sys,os

#from ROOT import *
import utils as u

effic = {
    "ee":  [0.8469, 0.0369],
    "mumu":[0.8541, 0.0355],
    } 
accep = {
    "mumu8TeV": {
        "A0tt": [ 2.27e-5, 0.000087],
        "A1tt": [11.27e-5, 0.000013],
        "A2tt": [11.69e-5, 0.000006],
        "A3tt": [ 0.22e-5, 0.000006],
        "A4tt": [ 0.04e-5, 0.000006],

        "A0t":  [ 7.19e-5, 0.000008],
        "A1t":  [15.87e-5, 0.000004],
        "A2t":  [ 3.95e-5, 0.000008],
        "A3t":  [ 0.10e-5, 0.000004],

        },
        "ee8TeV": {
        "A0tt": [ 1.87e-5, 0.000087],
        "A1tt": [ 8.19e-5, 0.000013],
        "A2tt": [ 8.52e-5, 0.000006],
        "A3tt": [ 0.13e-5, 0.000006],
        "A4tt": [ 0.01e-5, 0.000006],

        "A0t":  [ 4.95e-5, 0.000008],
        "A1t":  [12.53e-5, 0.000004],
        "A2t":  [ 3.85e-5, 0.000008],
        "A3t":  [ 0.00e-5, 0.000004],
        },

    "mumu7TeV": {
        "A0tt": [ 4.36e-5, 0.000087],
        "A1tt": [22.81e-5, 0.000013],
        "A2tt": [25.76e-5, 0.000006],
        "A3tt": [ 0.33e-5, 0.000006],
        "A4tt": [ 0.04e-5, 0.000006],

        "A0t":  [0.000215, 0.000008],
        "A1t":  [0.000059, 0.000004],
        "A2t":  [0.000215, 0.000008],
        "A3t":  [0.000059, 0.000004],
        },
    "ee7TeV": {
        "A2tt": [0.001363, 0.000073],
        "A1tt": [0.001254, 0.000011],
        "A0tt": [0.000255, 0.000005],
        "A2t":  [0.000054, 0.000006],
        "A1t":  [0.000173, 0.000010],
        "A0t":  [0.000061, 0.000006],
        },

        #Radek's:
        #"A2tt": [0.001956, 0.000087],
        #"A1tt": [0.001798, 0.000013],
        #"A0tt": [0.000385, 0.000006],
        #"A2t":  [0.000057, 0.000006],
        #"A1t":  [0.000215, 0.000008],
        #"A0t":  [0.000059, 0.000004],
        
    } 

myParams = u.myParams

def topbg(h_ttbar, h_tW, h_tbarW, sel, en):
    # Input hystograms are the b-jet multiplicity hists
    N0 = N1 = N2 = 0
    if sel ==1: selstr = "mumu"+en
    else: selstr = "ee"+en
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

    A2t = accep[selstr]["A2t"][0]
    A1t = accep[selstr]["A1t"][0]
    A0t = accep[selstr]["A0t"][0]
 
    eps = N2*(2*A2tt+A1tt)/(N1+2*N2)/A2tt
    Ntt = N2/(eps**2*A2tt)

    N0_calc = Ntt*(A2tt*(1-eps)**2 + A1tt*(1-eps) + A0tt) 
    
    print "Calculated:\n\
    Ntt = %.2f\n\
    eff = %.2f\n\
    N0 = %.2f\n\
    " % (Ntt, eps, N0_calc)

    sigma_tt = myParams.getCS("ttbar")
    sigma_t  = 2*myParams.getCS("tW") # multiply by two because of tbarW 

    if h_tW!=None and h_tbarW!=None:
        N0 += (h_tW.GetBinContent(1) + h_tbarW.GetBinContent(1)) 
        N1 += (h_tW.GetBinContent(2) + h_tbarW.GetBinContent(2))
        N2 += (h_tW.GetBinContent(3) + h_tbarW.GetBinContent(3))
    
    print " tt + top Numbers:\n\
    N2 = %.2f\n\
    N1 = %.2f\n\
    N0 = %.2f\n\
    "  % (N2,N1,N0)

    print "*** Including single top. cs(tt) = ", sigma_tt, " cs(t) =", sigma_t


    eps = N2*(2*A2tt+A1tt + (sigma_t/sigma_tt)*A1t)/(N1+2*N2)/A2tt
    Ntt = N2/(eps**2*A2tt)

    N0_calc = Ntt*(A2tt*(1-eps)**2 + A1tt*(1-eps) + A0tt + (sigma_t/sigma_tt)*(A1t*(1-eps)+ A0t)) 

    print "Calculated with signle top:\n\
    Ntt = %.2f\n\
    eff = %.2f\n\
    N0 = %.2f\n\
    " % (Ntt, eps, N0_calc)

    eps = N2*(2*A2tt+A1tt + (sigma_t/sigma_tt)*(2*A2t+A1t))/(N1+2*N2)/(A2tt + (sigma_t/sigma_tt)*A2t)
    Ntt = N2/(A2tt + (sigma_t/sigma_tt)*A2t)/eps**2
    N0_calc = Ntt*(A2tt*(1-eps)**2 + A1tt*(1-eps) + A0tt + (sigma_t/sigma_tt)*(A2t*(1-eps)**2 + A1t*(1-eps) + A0t)) 

    print "With a contribution from A2t:\n\
    Ntt = %.2f\n\
    eff = %.2f\n\
    N0 = %.2f\n\
    " % (Ntt, eps, N0_calc)

def CalcAccept(gen_b_Nb30):    
    nev1 = myParams.getNev("ttbar")
    
    sigma_tt = myParams.getCS("ttbar")
    sigma_t  = myParams.getCS("tW")

    gen_b_Nb30["ttbar"].Scale(1./nev1)
    A0tt = (gen_b_Nb30["ttbar"].GetBinContent(1))
    A1tt = (gen_b_Nb30["ttbar"].GetBinContent(2))
    A2tt = (gen_b_Nb30["ttbar"].GetBinContent(3))
    A3tt = (gen_b_Nb30["ttbar"].GetBinContent(4))
    A4tt = (gen_b_Nb30["ttbar"].GetBinContent(5))
    A5tt = (gen_b_Nb30["ttbar"].GetBinContent(6))
    print "nev ttbar =", nev1
    print "Acceptances for ttbar, using 1fb lumi:"
    print "A0tt = %.7f\n A1tt = %.7f \n A2tt = %.7f \n A3tt = %.7f\n  A4tt = %.7f \nA5tt = %.7f" % (A0tt, A1tt, A2tt, A3tt, A4tt, A5tt)
    
    nev2 = myParams.getNev("tW")
    nev3 = myParams.getNev("tbarW")

    gen_b_Nb30["tW"].Scale(1./nev2)
    gen_b_Nb30["tbarW"].Scale(1./nev3)
    
    A0t = 0.5*(gen_b_Nb30["tW"].GetBinContent(1) + gen_b_Nb30["tbarW"].GetBinContent(1))
    A1t = 0.5*(gen_b_Nb30["tW"].GetBinContent(2) + gen_b_Nb30["tbarW"].GetBinContent(2))
    A2t = 0.5*(gen_b_Nb30["tW"].GetBinContent(3) + gen_b_Nb30["tbarW"].GetBinContent(3))
    A3t = 0.5*(gen_b_Nb30["tW"].GetBinContent(4) + gen_b_Nb30["tbarW"].GetBinContent(4))
    A4t = 0.5*(gen_b_Nb30["tW"].GetBinContent(5) + gen_b_Nb30["tbarW"].GetBinContent(5))
    A5t = 0.5*(gen_b_Nb30["tW"].GetBinContent(6) + gen_b_Nb30["tbarW"].GetBinContent(6))
    print "At 1 fb, nev tW", nev2, "nev tbarW", nev3
    print "\n tW\n A0t = %.7f\n A1t = %.7f \n A2t = %.7f \n A3t = %.7f\n  A4t = %.7f \nA5t = %.7f" % (A0t, A1t, A2t, A3t, A4t, A5t)
    print "nev total", nev1+nev2+nev3
