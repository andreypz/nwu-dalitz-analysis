#!/usr/bin/env python
import os,sys
from optparse import OptionParser
parser = OptionParser(usage="usage: %prog ver [options --noSys]")
parser.add_option("--noSyst", dest="noSyst",action="store_true", default=False, help="Don't include systematics in the limit calculation")
parser.add_option("--comb", dest="comb",action="store_true", default=False, help="Combine with electrons")
(options, args) = parser.parse_args()

import ConfigParser as cp
cf = cp.ConfigParser()
cf.read('config.cfg')

def produceLimits(inputFolder = '', outPutFolder = 'limitOutputs/', mass = '125.0'):
  print None

if __name__ == "__main__":
  
  massList   = [a.strip() for a in (cf.get("fits","massList-more")).split(',')]
  s = cf.get("fits","ver")
  dir = s
  if options.noSyst:
    dir = s.replace("/","")+"-noSyst"    
    print s, dir
    os.system("cp -r "+s+"  "+dir)
    
  for m in massList:
    print "\n**\t Making limits for M=",m,"\t**\n"
    #card = "output_cards/hzg_el_2012_cat0_M"+m+"_Dalitz.txt"
    card_mu  = dir+"/output_cards/hzg_mu_2012_cat0_M"+m+"_Dalitz.txt"
    card_ele = "./cards_ele2/realistic-counting-experiment_"+m[0:3]+".txt"
    card = ""
    if options.comb:
      print "not yet"
    else:
      card = card_mu

    if options.noSyst:
      os.system('combine -M Asymptotic '+card + ' -m '+m+' -S 0')
    else:
      os.system('combine -M Asymptotic '+card + ' -m '+m)
      
  os.system("mv higgsCombineTest.Asymptotic* "+dir)
  os.system("./limitPlotter.py "+dir)
  #os.system('combine -M ProfileLikelihood '+sys.argv[1]+ ' -m 125 -t 100')
