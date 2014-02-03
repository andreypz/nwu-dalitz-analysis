#!/usr/bin/env python
import sys
from ROOT import *
gROOT.SetBatch()
from systematics import *
from rooFitBuilder import *
sys.path.append("../zgamma")
import utils as u
#from collections import defaultdict
import ConfigParser as cp
cf = cp.ConfigParser()
cf.read('config.cfg')
s = cf.get("fits","ver")

#################################################
# We're finally ready to make the datacards     #
# Pull the info out, make the cards, write them #
# Run with: combine -M Asymptotic datacard.txt  #
#################################################

makeSys = 1

def makeCards(subdir):

  yearList   = [a.strip() for a in (cf.get("fits","yearList")).split(',')]
  leptonList = [a.strip() for a in (cf.get("fits","leptonList")).split(',')]
  catList    = [a.strip() for a in (cf.get("fits","catList")).split(',')]
  massList   = [a.strip() for a in (cf.get("fits","massList")).split(',')]
    
  for year in yearList:
    for lepton in leptonList:
      for cat in catList:
        bgFileName = subdir+'/testCardBackground_Dalitz.root'
        bgFile = TFile(bgFileName)
        bgWs = bgFile.Get('ws_card')
       
        sigNameList = ['gg','vbf']
        channel = '_'.join([lepton,year,'cat'+cat])
        #bkgParams = ['sigma','mean','tau','norm']
        bkgParams = ['p0','p1','p2','p3','p4','sigma','step','mean','norm']

        for mass in massList:
          
          sigFileName = subdir+'/'+'_'.join(['SignalOutput',lepton,year,'cat'+cat,mass])+'.root'
          sigFile = TFile(sigFileName)
          sigWs = sigFile.Get('ws_card')
          prefixSigList = ['sig_'+sig for sig in sigNameList]

          u.createDir(subdir+'/output_cards/')
          card = open(subdir+'/output_cards/'+'_'.join(['hzg',lepton,year,'cat'+cat,'M'+mass,'Dalitz'])+'.txt','w')
          
          card.write('#This is a card produced by a cardMaker.py script using the information from a root file, which containes PDFs for Data model ans signal shapes. Normalization (yields) in signal are taken from those PDFs\n')
          
          card.write('imax *\n')
          card.write('jmax *\n')
          card.write('kmax *\n')
          card.write('---------------\n')
          card.write('shapes {0:<8} * {1:<20} ws_card:$PROCESS_$CHANNEL\n'.format('*',bgFileName))
          card.write('shapes {0:<8} * {1:<20} ws_card:bkg_$CHANNEL\n'.format('bkg',bgFileName))
          for sig in prefixSigList:
            card.write('shapes {0:<8} * {1:<20} ws_card:{2}_$CHANNEL\n'.format(sig,sigFileName,sig))
          card.write('---------------\n')
          bgYield = bgWs.var('data_yield_'+channel).getVal()
          print cat, "bg yield:", bgYield
          
          card.write('{0:<12} {1}\n'.format('bin',channel))
          card.write('{0:<12} {1}\n'.format('observation',int(bgYield)))
          card.write('------------------------------\n')

          print "WTF is [::-1] ??? ->"
          print channel, prefixSigList, prefixSigList[::-1]
          
          card.write('{0:<25} {1:^15} {2:^15} {3:^15}\n'.format(*(['bin']+[channel]*3)))
          card.write('{0:<25} {1:^15} {2:^15} {3:^15}\n'.format(*(['process']+prefixSigList[::-1]+['bkg'])))
          card.write('{0:<25} {1:^15} {2:^15} {3:^15}\n'.format(*(['process', -1,0,1])))
          
          card.write('--------------------------------------------------------------\n')
          sigYields = []
          for sig in prefixSigList[::-1]:
            sigYields.append(sigWs.var(sig+'_yield_'+channel).getVal())

          print sigYields
          card.write('{0:<25} {1:^15.5} {2:^15} {3:^15}\n'.format(*(['rate']+sigYields+[1])))

          card.write('-------------------------------------------------------------\n')
          card.write(' \n')
          
          if makeSys:
            
            card.write('lumi         lnN     1.026        1.026   1.026\n')
            card.write('brLLG        lnN     1.10         1.10  -     \n')
            card.write('muonID       lnN     1.11         1.11  -     \n')
            card.write('pdf_vbf       lnN   '+pdf_vbf[year][mass]+'  -   -  \n')
            card.write('QCDscale_vbfH lnN   '+qcd_vbf[year][mass]+'  -   -  \n')
            card.write('pdf_gg       lnN     -   '+pdf_gg[year][mass]+'  -     \n')
            card.write('QCDscale_ggH lnN     -   '+qcd_gg[year][mass]+'  -     \n')

            # card.write('unc_Bkg1  gmN 72      -           1.00  \n')
            
            
            for sig in prefixSigList:
              card.write('{0:<40} {1:<10} {2:^10} {3:^10}\n'.format(sig+'_mShift_'+channel,'param', 1, 0.01))
              card.write('{0:<40} {1:<10} {2:^10} {3:^10}\n'.format(sig+'_sigmaShift_'+channel,'param', 1, 0.05))
              
            for param in bkgParams[:-1]:
              card.write('{0:<45} {1:<15}\n'.format('bkg_'+param+'_'+channel,'flatParam'))
            card.write('{0:<45} {1:<15}\n'.format('bkg_'+channel+'_'+bkgParams[-1],'flatParam'))


          card.close()




if __name__=="__main__":
  
  print len(sys.argv), sys.argv
  
  makeCards(s)
