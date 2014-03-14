#!/usr/bin/env python
import sys
from ROOT import *
gROOT.SetBatch()
from systematics import *
sys.path.append("../zgamma")
import utils as u
from optparse import OptionParser
parser = OptionParser(usage="usage: %prog  [options --brSyst 1.20]")
parser.add_option("--brSyst", dest="brSyst", default='1.10', help="error on the branching ratio")
parser.add_option("-b", dest="bach", action="store_true", default=True, help="batch")
parser.add_option("--ext",dest="ext", action="store_true", default=False, help="Extended pdf for bkg")
parser.add_option("--br",dest="br", action="store_true", default=False, help="Do the limit on BR instead of the mu")
(options, args) = parser.parse_args()

import ConfigParser as cp
cf = cp.ConfigParser()
cf.read('config.cfg')
s = cf.get("path","ver")

# ################################################
# We're finally ready to make the datacards     #
# Pull the info out, make the cards, write them #
# Run with: combine -M Asymptotic datacard.txt  #
# ################################################
brTag = ''

Ext = options.ext

brLLG = '1.10'
# Unsertainities to be put in the datacard:
lumi    = '1.026'
muID    = '1.110'
muISO   = '1.003'
muTRIG  = '1.040'
phoID   = '1.006'
phoTRIG = '1.020'
PU      = '1.008'

proc = {'gg':'ggH', 'vbf':'qqH','v':'WH'}
massList   = ['%.1f'%(a) for a in u.drange(120,150,.5)]

def makeCards(subdir):
  yearList   = [a.strip() for a in (cf.get("fits","yearList")).split(',')]
  leptonList = [a.strip() for a in (cf.get("fits","leptonList")).split(',')]
  catList    = [a.strip() for a in (cf.get("fits","catList")).split(',')]
  sigNameList= [a.strip() for a in (cf.get("fits","signameList")).split(',')]

  for year in yearList:
    for lepton in leptonList:
      for cat in catList:
        bgFileName = subdir+'/testCardBackground_Dalitz.root'
        bgFile = TFile(bgFileName)
        bgWs   = bgFile.Get('ws_card')
        bgWs.Print()

        channel = '_'.join([lepton,year,'cat'+cat])
        #bkgParams = ['sigma','mean','tau','norm']
        bkgParams = ['p1','p2','p3','p4','norm']

        for mass in massList:

          sigFileName = subdir+'/'+'_'.join(['SignalOutput',lepton,year,'cat'+cat,mass])+'.root'
          sigFile = TFile(sigFileName)
          sigWs = sigFile.Get('ws_card')
          procList  = [proc[a] for a in sigNameList]
          print procList
          if options.br and float(mass)%5!=0:
            print "Sorry can't do this with those mass points yet, ",mass
            sys.exit(0)
          if options.br:
            croDict = {a:float(u.conf.get(a+"H-"+mass[0:3], "cs-mu")) for a in sigNameList}
            print croDict

          u.createDir(subdir+'/output_cards/')
          card = open(subdir+'/output_cards/'+'_'.join(['hzg',lepton,year,'cat'+cat,'M'+mass,'Dalitz',brTag])+'.txt','w')

          card.write('# This is a card produced by a cardMaker.py script using the information from a workspace in a root file \n# That file containes PDFs for Data model and the signals\' shapes. \n# Normalization (yields) for the signal are taken from that workspace.\n')

          card.write('imax *\n')
          card.write('jmax *\n')
          card.write('kmax *\n')
          card.write('---------------\n')
          card.write('shapes {0:<8} * {1:<20} ws_card:$PROCESS_$CHANNEL\n'.format('*',bgFileName))
          card.write('shapes {0:<8} * {1:<20} ws_card:bkg_$CHANNEL\n'.format('bkg',bgFileName))
          for sig in sigNameList:
            card.write('shapes {0:<8} * {1:<20} ws_card:{2}_$CHANNEL\n'.format(proc[sig],sigFileName,'sig_'+sig))

          card.write('---------------\n')
          bgYield = bgWs.var('data_yield_'+channel).getVal()
          print cat, "bg yield:", bgYield
          bkgRate = 1
          #if Ext: bkgRate = '1'
          #else: bkgRate = str(bgYield)
          card.write('{0:<12} {1}\n'.format('bin',channel))
          card.write('{0:<12} {1}\n'.format('observation',int(bgYield)))
          card.write('------------------------------\n')

          #print "WTF is [::-1] ??? ->"
          #print channel, prefixSigList, prefixSigList[::-1]

          card.write('{0:<25} {1:^15} {2:^15} {3:^15} {4:^15}\n'.format(*(['bin']+[channel]*4)))
          card.write('{0:<25} {1:^15} {2:^15} {3:^15} {4:^15}\n'.format(*(['process']+procList[::-1]+['bkg'])))
          card.write('{0:<25} {1:^15} {2:^15} {3:^15} {4:^15}\n'.format(*(['process', -2,-1,0,1])))

          card.write('--------------------------------------------------------------\n')
          sigYields = []
          for s in sigNameList[::-1]:
            cs = 1
            if options.br:
              cs = croDict[s]
            print s, cs
            sigYields.append(sigWs.var('sig_'+s+'_yield_'+channel).getVal() /cs)

          print 'mh=',mass, sigYields
          sigRate = sigYields
          card.write('{0:<25} {1:^15.4} {2:^15.4} {3:^15.5} {4:^15}\n'.format(*(['rate']+sigYields+[bkgRate])))

          card.write('-------------------------------------------------------------\n')
          card.write(' \n')

          card.write('{0} {1}  {2} {3} {4} -  \n'.format(*(['lumi_8TeV',      'lnN']+3*[lumi])))
          card.write('{0} {1}  {2} {3} {4} -  \n'.format(*(['CMS_hllg_brLLG', 'lnN']+3*[brLLG])))
          card.write('{0} {1}  {2} {3} {4} -  \n'.format(*(['CMS_eff_m_ID',   'lnN']+3*[muID])))
          card.write('{0} {1}  {2} {3} {4} -  \n'.format(*(['CMS_eff_m_ISO',  'lnN']+3*[muISO])))
          card.write('{0} {1}  {2} {3} {4} -  \n'.format(*(['CMS_eff_m_TRIG', 'lnN']+3*[muTRIG])))
          card.write('{0} {1}  {2} {3} {4} -  \n'.format(*(['CMS_eff_g_ID',   'lnN']+3*[phoID])))
          card.write('{0} {1}  {2} {3} {4} -  \n'.format(*(['CMS_eff_g_TRIG', 'lnN']+3*[phoTRIG])))
          card.write('{0} {1}  {2} {3} {4} -  \n'.format(*(['CMS_hllg_PU',    'lnN']+3*[PU])))
          mmm = mass
          if float(mass)>140:
            mmm = mass[0:4]+'0'

          if not options.br:
            card.write('pdf_WH        lnN     '+pdf_wh[year][mmm]+'  -  -   -  \n')
            card.write('QCDscale_WH   lnN     '+qcd_wh[year][mmm]+'  -  -   -  \n')
            card.write('pdf_qqH       lnN     -   '+pdf_vbf[year][mmm]+' -  -  \n')
            card.write('QCDscale_qqH  lnN     -   '+qcd_vbf[year][mmm]+' -  -  \n')
            card.write('pdf_ggH       lnN     -   -  '+pdf_gg[year][mmm]+'  -  \n')
            card.write('QCDscale_ggH  lnN     -   -  '+qcd_gg[year][mmm]+'  -  \n')

          for sig in sigNameList:
            card.write('{0:<40} {1:<10} {2:^10} {3:^10}\n'.format('sig_'+sig+'_mShift_'    +channel,'param', 1, 0.005))
            card.write('{0:<40} {1:<10} {2:^10} {3:^10}\n'.format('sig_'+sig+'_sigmaShift_'+channel,'param', 1, 0.05))

          for param in bkgParams[:-1]:
            card.write('{0:<45} {1:<15}\n'.format('bkg_'+param+'_'+channel,'flatParam'))
          card.write('{0:<45} {1:<15}\n'.format('bkg_'+channel+'_'+bkgParams[-1],'flatParam'))


          card.close()




if __name__=="__main__":

  print len(sys.argv), sys.argv



  mmm = massList[0]
  for m in massList[1:]:
    mmm += ','+m
  cf.set("fits","masslist-more",mmm)
  with open(r'config.cfg', 'wb') as configfile:
    cf.write(configfile)

  makeCards(s)
