#!/usr/bin/env python
import os,sys
sys.path.append("../zgamma")
import utils as u
from optparse import OptionParser
parser = OptionParser(usage="usage: %prog ver [options --noSys]")
parser.add_option("--noSyst", dest="noSyst",action="store_true", default=False,
                  help="Don't include systematics in the limit calculation")
parser.add_option("--comb", dest="comb",action="store_true", default=False,
                  help="Combine with electrons")
parser.add_option("--prof", dest="prof",action="store_true", default=False,
                  help="Use -M ProfileLikelihood instead of Asymptotic")

(options, args) = parser.parse_args()

import ConfigParser as cp
cf = cp.ConfigParser()
cf.read('config.cfg')

method = "Asymptotic"
if options.prof:
  method = "ProfileLikelihood"

if __name__ == "__main__":

  #massList   = ['%.1f'%(a) for a in u.drange(120,150,1)]
  massList   = [a.strip() for a in (cf.get("fits","massList-more")).split(',')]
  catList    = [a.strip() for a in (cf.get("fits","catList")).split(',')]

  catList = ['EB','EE','mll50']

  s = cf.get("path","ver")
  myDir = s

  if options.noSyst:
    myDir = s.replace("/","")+"-noSyst"
    print s, myDir
    os.system("cp -r "+s+"  "+myDir)
    cf.set("path","ver", myDir)
    with open(r'config.cfg', 'wb') as configfile:
      cf.write(configfile)

  for m in massList:
    cardNames = ''
    for cat in catList:
      print "\n**\t Making limits for M=",m,"  cat = ",cat, "\t**\n"
      # card = "output_cards/hzg_el_2012_cat0_M"+m+"_Dalitz.txt"
      card_mu  = myDir+"/output_cards/hzg_mu_2012_cat"+cat+"_M"+m+"_Dalitz_.txt"
      card_el  = myDir+"/Dalitz_electron2/datacard_hzg_eeg_cat12_8TeV_"+m[:3]+".txt"

      if cat in ['EB','EE','mll50']:
        cardNames = cardNames+' '+card_mu
      if cat in ['el']:
        cardNames = cardNames+' '+card_el
      # card_ele = "./cards_ele2/realistic-counting-experiment_"+m[0:3]+".txt"
      card = ""

      if options.comb:
        print "not yet"
      elif cat in ['el']:
        card = card_el
      else:
        card = card_mu

      print 'Using CARD = \n',card

      if options.noSyst:
        os.system('combine -M '+method+' '+card + ' -m '+m+' -S 0')
      else:
        print 'not doing it now'
        os.system('combine -M '+method+' '+card + ' -m '+m)

      if m[-1]=='5':
        fname = "higgsCombineTest."+method+".mH"+m
        os.system("mv "+fname+".root "+myDir+'/'+fname+'_cat'+cat+'.root')
      else:
        fname = "higgsCombineTest."+method+".mH"+m[0:3]
        os.system("mv "+fname+".root "+myDir+'/'+fname+'.0_cat'+cat+'.root')

    print '\t Combining the cards: \n', cardNames
    comboName = myDir+'/output_cards/combo_mH'+m+'.txt'
    os.system('combineCards.py '+cardNames+' > '+comboName)
    os.system("sed -i 's:"+myDir+"/output_cards/"+s[:3]+":"+s[:3]+":g' "+comboName)

    os.system('combine -M '+method+' '+comboName + ' -m '+m)

    if m[-1]=='5':
      fname = "higgsCombineTest."+method+".mH"+m
      os.system("mv "+fname+".root "+myDir+'/'+fname+'_catCombo.root')
    else:
      fname = "higgsCombineTest."+method+".mH"+m[0:3]
      os.system("mv "+fname+".root "+myDir+'/'+fname+'.0_catCombo.root')


    print myDir, s[:3]

  for cat in catList:
    os.system("./limitPlotter.py "+myDir+" --cat="+cat)
  os.system("./limitPlotter.py "+myDir+" --cat=Combo")
    # os.system('combine -M ProfileLikelihood '+sys.argv[1]+ ' -m 125 -t 100')
