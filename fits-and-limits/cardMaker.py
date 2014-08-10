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
parser.add_option("--br",dest="br", action="store_true", default=False, help="Do the limit on BR*cs instead of the mu")
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

proc = {'gg':'ggH', 'vbf':'qqH','v':'WH','hjp':'hjp'}
#massList   = ['%.1f'%(a) for a in u.drange(120,150,5)]
#massList   = ['%.1f'%(a) for a in u.drange(120,150,.5)]
massList = ['125.0']
csBR = {}

#for i,m in enumerate(massList):
#  # for now this numbers are retrieved from systematics.py file
#  ccc = cs_tot[i]*BR_mu[i]*1000
#  print i,m,"CS, BR and cs*Br = ", cs_tot[i],BR_mu[i], ccc
#  ## this is wrong really, cant use total cross section..
#  ## don't use this
#  csBR[m] = ccc

if options.br:
  f = TFile('../data/Dalitz_BR50.root','READ')
  g = f.Get('csbr_mu')
  fit = g.GetFunction('pol4')
  #g.Print('all')

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
        bkgParams = ['p1','p2','p3','norm']
        #bkgParams = ['p1','p2','p3','p4','norm']

        for mass in massList:

          sigFileName = subdir+'/'+'_'.join(['SignalOutput',lepton,year,'cat'+cat,mass])+'.root'
          sigFile = TFile(sigFileName)
          sigWs = sigFile.Get('ws_card')
          procList  = [proc[a] for a in sigNameList]
          if options.br:
            procList=['hjp']
            #procList=['ggH']

          nProc = len(procList)
          # print procList
          #if options.br and float(mass)%5!=0:
          #  print "Sorry can't do this with those mass points yet, ",mass
          #  sys.exit(0)

          #  croDict = {a:float(u.conf.get(a+"H-"+mass[0:3], "cs-mu")) for a in sigNameList}
          #  print croDict

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
          nSubLine1 = '{0:<25}'
          nSubLine2 = '{0:<30}'
          for i in xrange(nProc):
            nSubLine1 +=' {'+str(i+1)+':^15}'
            nSubLine2 +=' {'+str(i+1)+':^5}'
            #print i, nSubLine

          nSubLine1 +=' {'+str(i+2)+':^15}\n'
          nSubLine2 +=' {'+str(i+2)+':5} - \n'

          print nSubLine1

          card.write(nSubLine1.format(*(['bin']+[channel]*(nProc+1))))
          card.write(nSubLine1.format(*(['process']+procList[::-1]+['bkg'])))
          card.write(nSubLine1.format(*(['process']+range(-(nProc-1), 2, 1))))

          card.write('--------------------------------------------------------------\n')
          sigYields = []
          for s in sigNameList[::-1]:
            cs = 1
            if options.br:
              #cs = fit(float(mass))
              cs = 0.5
              # print 'from the fit', fit(float(mass))
              # dont use this:: cs = csBR[mass]
              #print 'from csbr:',cs
            print mass, s, cs
            #sigYields.append(sigWs.var('sig_'+s+'_yield_'+channel).getVal()*cs)
            sigYields.append(sigWs.var('sig_'+s+'_yield_'+channel).getVal() /cs)

          print 'mh=',mass, sigYields
          sigRate = sigYields
          card.write(nSubLine1.format(*(['rate']+sigRate+[bkgRate])))

          card.write('-------------------------------------------------------------\n')
          card.write(' \n')

          card.write(nSubLine2.format(*(['lumi_8TeV',      'lnN']+nProc*[lumi])))
          mmm = mass

          if float(mass)>140:
            # after mH=140 the syst only available with 1GeV intervals
            mmm = mass[0:4]+'0'

          if not options.br:
            card.write(nSubLine2.format(*(['CMS_hllg_brLLG', 'lnN']+nProc*[brLLG])))
            card.write('pdf_WH        lnN     '+pdf_wh[year][mmm]+'  -  -   -  \n')
            card.write('QCDscale_WH   lnN     '+qcd_wh[year][mmm]+'  -  -   -  \n')
            card.write('pdf_qqH       lnN     -   '+pdf_vbf[year][mmm]+' -  -  \n')
            card.write('QCDscale_qqH  lnN     -   '+qcd_vbf[year][mmm]+' -  -  \n')
            card.write('pdf_ggH       lnN     -   -  '+pdf_gg[year][mmm]+'  -  \n')
            card.write('QCDscale_ggH  lnN     -   -  '+qcd_gg[year][mmm]+'  -  \n')

          card.write(nSubLine2.format(*(['CMS_eff_m_ID',   'lnN']+nProc*[muID])))
          card.write(nSubLine2.format(*(['CMS_eff_m_ISO',  'lnN']+nProc*[muISO])))
          card.write(nSubLine2.format(*(['CMS_eff_m_TRIG', 'lnN']+nProc*[muTRIG])))
          card.write(nSubLine2.format(*(['CMS_eff_g_ID',   'lnN']+nProc*[phoID])))
          card.write(nSubLine2.format(*(['CMS_eff_g_TRIG', 'lnN']+nProc*[phoTRIG])))
          card.write(nSubLine2.format(*(['CMS_hllg_PU',    'lnN']+nProc*[PU])))

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
