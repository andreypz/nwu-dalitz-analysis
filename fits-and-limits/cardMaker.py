#!/usr/bin/env python
import sys
from ROOT import *
gROOT.SetBatch()
sys.path.append("../zgamma")
sys.path.append("../scripts")
from xsBrReader import *
import utils as u
from optparse import OptionParser
parser = OptionParser(usage="usage: %prog  [options --brSyst 1.20]")
parser.add_option("--brSyst", dest="brSyst", default='1.10', help="error on the branching ratio")
parser.add_option("-b", dest="bach", action="store_true", default=True, help="batch")
parser.add_option("--ext",dest="ext", action="store_true", default=False, help="Extended pdf for bkg")
parser.add_option("--br",dest="br", action="store_true", default=False, help="Do the limit on BR*cs instead of the mu")
(opt, args) = parser.parse_args()

import ConfigParser as cp
cf = cp.ConfigParser()
cf.read('config.cfg')
s = cf.get("path","ver")

# ################################################
# We're finally ready to make the datacards     #
# Pull the info out, make the cards, write them #
# Run with: combine -M Asymptotic datacard.txt  #
# ################################################
lumi2012  = 19.703

massList   = ['%.1f'%(a) for a in u.drange(120,150,1.)]
sigNameList= [a.strip() for a in (cf.get("fits","signameList")).split(',')]
hjp = 0
if 'hjp' in sigNameList:
  hjp = 1
  massList = ['125.0']

brTag = ''

Ext = opt.ext

brLLG = '1.10'
# Unsertainities to be put in the datacard:
lumi     = '1.026'
elID     = '1.035'
muID     = '1.110'
muISO    = '1.003'
muTRIG   = '1.040'
phoID    = '1.006'
phoTRIG_mu = '1.020'
phoTRIG_el = '1.020'
PU       = '1.008'
meanUnc  = {'mu':'0.001', 'el':'0.005'}
sigmaUnc = {'mu':'0.100', 'el':'0.100'}

yearToTeV ={'2011':'7TeV', '2012':'8TeV'}
proc = {'gg':'ggH', 'vbf':'qqH','v':'VH','hjp':'hjp'}
proc_conf = {'gg':'ggH', 'vbf':'vbfH','v':'vH'}
#proc = {'gg':'ggH', 'vbf':'qqH','v':'WH','hjp':'hjp'}
#prYR = {'gg':'ggF', 'vbf':'VBF','v':'WH'}
mllBins = u.mllBins()

AEFF = u.AutoVivification()
#csBR = {}
#for i,m in enumerate(massList):
#  # for now this numbers are retrieved from systematics.py file
#  ccc = cs_tot[i]*BR_mu[i]*1000
#  print i,m,"CS, BR and cs*Br = ", cs_tot[i],BR_mu[i], ccc
#  ## this is wrong really, cant use total cross section..
#  ## don't use this
#  csBR[m] = ccc

introMessage = ('# This is a card produced by a cardMaker.py script using the information from a workspace in a root file\n'\
                  '# That file contains the PDFs for Data model and the signals\' shapes. \n#'\
                  '# Normalization (yields) for the signal are taken from that workspace.\n')


if opt.br:
  f = TFile('../data/Dalitz_BR20.root','READ')
  gMu = f.Get('csbr_mu')
  fitMu = gMu.GetFunction('pol4')

  gEl = f.Get('csbr_el')
  fitEl = gEl.GetFunction('pol4')
  # g.Print('all')

def makeCards(subdir):
  yearList   = [a.strip() for a in (cf.get("fits","yearList")).split(',')]
  leptonList = [a.strip() for a in (cf.get("fits","leptonList")).split(',')]
  catList    = [a.strip() for a in (cf.get("fits","catList")).split(',')]
  if hjp: leptonList =['mu']

  for year in yearList:
    for lep in leptonList:
      for cat in catList:
        if lep=='el' and cat!='EB': continue
        bgFileName = subdir+'/testCardBackground_Dalitz.root'
        bgFile = TFile(bgFileName)
        bgWs   = bgFile.Get('ws_card')
        bgWs.Print()

        channel = '_'.join([lep,year,'cat'+cat])
        print channel
        #bkgParams = ['sigma','mean','tau','norm']
        if hjp:
          bkgParams = ['p1','p2','norm']
        else:
          bkgParams = ['p1','p2','p3','p4','norm']

        for mass in massList:
          if cat in ['m1','m2','m3','m4','m5','m6','m7'] and (mass!='125.0' or not opt.br): continue

          sigFileName = subdir+'/'+'_'.join(['SignalOutput',lep,year,'cat'+cat,mass])+'.root'
          sigFile = TFile(sigFileName)
          sigWs = sigFile.Get('ws_card')
          #sigWs.Print()

          procList  = []
          for s in sigNameList:
            if lep=='el' and s=='v': continue
            elif cat in ['m1','m2','m3','m4','m5','m6','m7'] and s!='gg': continue
            elif s=='v':
              procList.extend(['WH','ZH'])
            else:
              procList.append(proc[s])

          if opt.br:
            if hjp:
              procList=['hjp']
            else:
              procList=['ggH']

          print 'List of procs:', procList

          nProc = len(procList)
          # print procList
          #if opt.br and float(mass)%5!=0:
          #  print "Sorry can't do this with those mass points yet, ",mass
          #  sys.exit(0)

          #  croDict = {a:float(u.conf.get(a+"H-"+mass[0:3], "cs-mu")) for a in sigNameList}
          #  print croDict
          if opt.br: outCardDir = '/output_cards_xsbr/'
          else: outCardDir = '/output_cards/'
          u.createDir(subdir+outCardDir)
          card = open(subdir+outCardDir+'_'.join(['hzg',lep,year,'cat'+cat,'M'+mass,'Dalitz',brTag])+'.txt','w')

          card.write(introMessage)
          card.write('imax *\n')
          card.write('jmax *\n')
          card.write('kmax *\n')
          card.write('---------------\n')
          card.write('shapes {0:<8} * {1:<20} ws_card:$PROCESS_$CHANNEL\n'.format('*',bgFileName))
          card.write('shapes {0:<8} * {1:<20} ws_card:bkg_$CHANNEL\n'.format('bkg',bgFileName))
          for sig in sigNameList:
            if sig=='v' and lep=='el': continue
            if opt.br and sig!='gg' and sig!='hjp': continue
            elif sig=='v':
              card.write('shapes {0:<8} * {1:<20} ws_card:{2}_$CHANNEL\n'.format('ZH',sigFileName,'sig_'+sig))
              card.write('shapes {0:<8} * {1:<20} ws_card:{2}_$CHANNEL\n'.format('WH',sigFileName,'sig_'+sig))
            else:
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

          #print channel, prefixSigList, prefixSigList[::-1]
          nSubLine1 = '{0:<25}'
          nSubLine2 = '{0:<30}'
          for i in xrange(nProc):
            nSubLine1 +=' {'+str(i+1)+':^15}'
            nSubLine2 +=' {'+str(i+1)+':^5}'

          # print i, nSubLine1

          nSubLine1 +=' {'+str(i+2)+':^15} \n'
          nSubLine2 +=' {'+str(i+2)+':5} - \n'

          # print 'nSubLine1=', nSubLine1

          card.write(nSubLine1.format(*(['bin']+[channel]*(nProc+1))))
          card.write(nSubLine1.format(*(['process']+procList[::-1]+['bkg'])))
          card.write(nSubLine1.format(*(['process']+range(-(nProc-1), 2, 1))))

          card.write('--------------------------------------------------------------\n')
          sigYields = []
          for s in sigNameList[::-1]:
            if lep=='el' and s=='v': continue
            if opt.br and s!='gg' and s!='hjp': continue
            cs = 1

            procYield = sigWs.var('sig_'+s+'_yield_'+channel).getVal()

            if opt.br:
              if hjp:
                cs = u.getCS("HtoJPsiGamma")
              else:
                if lep=='mu':   cs = fitMu(float(mass))
                elif lep=='el': cs = fitEl(float(mass))
                else: print 'This fit does not exist'

              if cat in ['m1','m2','m3','m4','m5','m6','m7']:
                cs = cs*mllBins[int(cat[1])][1]

              # When doing limit on cs*BR, need to divide by cs:
              procYield = sigWs.var('sig_'+s+'_yield_'+channel).getVal()/cs

              print '\t Category:', cat, '  signal:', s, 'm=', mass
              print '\t ** From the fit: %.3f' % cs
              print mass, s, lep, 'cs for scale=', cs


            print mass, s, lep
            if opt.br:
              print '\t \t raw sig= ', procYield*cs, 'Acc*Eff = %.3f' % (procYield/lumi2012)
            else:
              # cs = u.getCS('%s-%i' %(proc_conf[s], float(mass)), lep)
              print '\t \t raw sig= ', procYield, 'Acc*Eff = %.3f' % (procYield/lumi2012/cs)


            if s=='v':
              zh_frac = float(xsDict[YR][TeV]['ZH'][mass])/(float(xsDict[YR][TeV]['ZH'][mass])+float(xsDict[YR][TeV]['WH'][mass]))
              wh_frac = 1-zh_frac
              print s, 'YR=',YR, ' ZH fraction =', zh_frac, '  WH fraction =',wh_frac
              sigYields.append('%.4f'%(procYield*zh_frac))
              sigYields.append('%.4f'%(procYield*wh_frac))
            else:
              sigYields.append('%.4f'%procYield)

          print 'M =', mass, 'YIELDS=', sigYields

          sigRate = sigYields
          card.write(nSubLine1.format(*(['rate']+sigRate+[bkgRate])))

          card.write('-------------------------------------------------------------\n')
          card.write(' \n')

          card.write(nSubLine2.format(*(['lumi_8TeV', 'lnN']+nProc*[lumi])))
          mmm = mass

          if float(mass)>140:
            # after mH=140 the syst only available with 1GeV intervals
            mmm = mass[0:4]+'0'

          if not opt.br:
            card.write(nSubLine2.format(*(['CMS_hllg_brLLG', 'lnN']+nProc*[brLLG])))

            if lep=='el':
              card.write('pdf_qqbar     lnN       '+xsPDFErrDict['YR3'][yearToTeV[year]]['VBF'][mmm]+' -  -  \n')
              card.write('QCDscale_qqH  lnN       '+xsScaleErrDict['YR3'][yearToTeV[year]]['VBF'][mmm]+' -  -  \n')
              card.write('pdf_gg        lnN       -  '+xsPDFErrDict['YR3'][yearToTeV[year]]['ggF'][mmm]+'  -  \n')
              card.write('QCDscale_ggH  lnN       -  '+xsScaleErrDict['YR3'][yearToTeV[year]]['ggF'][mmm]+'  -  \n')
            else:
              card.write('QCDscale_ZH   lnN     '+xsScaleErrDict['YR3'][yearToTeV[year]]['ZH'][mmm]+'  - -  -   -  \n')
              card.write('QCDscale_WH   lnN     - '+xsScaleErrDict['YR3'][yearToTeV[year]]['WH'][mmm]+'  -  -   -  \n')
              card.write('pdf_qqbar     lnN     '+
                         xsPDFErrDict['YR3'][yearToTeV[year]]['ZH'][mmm]+' '+
                         xsPDFErrDict['YR3'][yearToTeV[year]]['WH'][mmm]+' '+
                         xsPDFErrDict['YR3'][yearToTeV[year]]['VBF'][mmm]+' -  -  \n')
              card.write('QCDscale_qqH  lnN      - - '+xsScaleErrDict['YR3'][yearToTeV[year]]['VBF'][mmm]+' -  -  \n')
              card.write('pdf_gg        lnN      - - -  '+xsPDFErrDict['YR3'][yearToTeV[year]]['ggF'][mmm]+'  -  \n')
              card.write('QCDscale_ggH  lnN      - - -  '+xsScaleErrDict['YR3'][yearToTeV[year]]['ggF'][mmm]+'  -  \n')

          if lep=='mu':
            card.write(nSubLine2.format(*(['CMS_eff_m',   'lnN']+nProc*[muID])))
            #card.write(nSubLine2.format(*(['CMS_eff_m_ISO',  'lnN']+nProc*[muISO])))
            card.write(nSubLine2.format(*(['CMS_eff_m_MuEgTRIG', 'lnN']+nProc*[muTRIG])))
            card.write(nSubLine2.format(*(['CMS_eff_g',   'lnN']+nProc*[phoID])))
            card.write(nSubLine2.format(*(['CMS_eff_g_MuEgTRIG', 'lnN']+nProc*[phoTRIG_mu])))

          else:
            card.write(nSubLine2.format(*(['CMS_eff_e',   'lnN']+nProc*[elID])))
            card.write(nSubLine2.format(*(['CMS_eff_g',   'lnN']+nProc*[phoID])))
            card.write(nSubLine2.format(*(['CMS_eff_g_GGTRIG', 'lnN']+nProc*[phoTRIG_el])))

          card.write(nSubLine2.format(*(['CMS_hllg_PU',    'lnN']+nProc*[PU])))

          for sig in sigNameList:
            if lep=='el' and sig=='v': continue
            card.write('{0:<40} {1:<10} {2:^10} {3:^10}\n'.format('sig_'+sig+'_mShift_'    +channel,'param', 1, meanUnc[lep]))
            card.write('{0:<40} {1:<10} {2:^10} {3:^10}\n'.format('sig_'+sig+'_sigmaShift_'+channel,'param', 1, sigmaUnc[lep]))

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
