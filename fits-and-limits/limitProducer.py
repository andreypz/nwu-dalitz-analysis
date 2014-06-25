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

def produceLimits(inputFolder = '', outPutFolder = 'limitOutputs/', mass = '125.0'):
  print None

if __name__ == "__main__":

  #massList   = ['%.1f'%(a) for a in u.drange(120,150,5)]
  massList   = [a.strip() for a in (cf.get("fits","massList-more")).split(',')]

  s = cf.get("path","ver")
  dir = s
  if options.noSyst:
    dir = s.replace("/","")+"-noSyst"
    print s, dir
    os.system("cp -r "+s+"  "+dir)
    cf.set("path","ver", dir)
    with open(r'config.cfg', 'wb') as configfile:
      cf.write(configfile)

  for m in massList:
    print "\n**\t Making limits for M=",m,"\t**\n"
    #card = "output_cards/hzg_el_2012_cat0_M"+m+"_Dalitz.txt"
    card_mu  = dir+"/output_cards/hzg_mu_2012_cat0_M"+m+"_Dalitz_.txt"

    #card_ele = "./cards_ele2/realistic-counting-experiment_"+m[0:3]+".txt"
    card = ""

    if options.comb:
      print "not yet"
    else:
      card = card_mu

    print 'Using CARD = \n',card


    method = "Asymptotic"
    if options.prof:
      method = "ProfileLikelihood"
    if options.noSyst:
      os.system('combine -M '+method+' '+card + ' -m '+m+' -S 0')
    else:
      os.system('combine -M '+method+' '+card + ' -m '+m)

  os.system("mv higgsCombineTest."+method+"* "+dir)
  print dir
  os.system("./limitPlotter.py "+dir)
  #os.system('combine -M ProfileLikelihood '+sys.argv[1]+ ' -m 125 -t 100')
