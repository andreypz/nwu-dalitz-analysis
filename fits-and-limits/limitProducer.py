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
parser.add_option("--plot", dest="plot",action="store_true", default=False,
                  help="Run the plotter script")
parser.add_option("--dont", dest="dont",action="store_true", default=False,
                  help="Don't run it, hust combine the cards")

(opt, args) = parser.parse_args()

import ConfigParser as cp
cf = cp.ConfigParser()
cf.read('config.cfg')

method = "Asymptotic"
if opt.prof:
  method = "ProfileLikelihood"

if __name__ == "__main__":
  #massList = ['125.0']
  #massList   = ['%.1f'%(a) for a in u.drange(120,150,1)]
  massList   = [a.strip() for a in (cf.get("fits","massList-more")).split(',')]
  catList    = [a.strip() for a in (cf.get("fits","catList")).split(',')]
  leptonList = [a.strip() for a in (cf.get("fits","leptonList")).split(',')]

  s = cf.get("path","ver")
  #s = "MingDataCards/Dalitz_electron/"
  myDir = s

  if opt.noSyst:
    myDir = s.replace("/","")+"-noSyst"
    print s, myDir
    os.system("cp -r "+s+"  "+myDir)
    cf.set("path","ver", myDir)
    with open(r'config.cfg', 'wb') as configfile:
      cf.write(configfile)

  for m in massList:
    cardNames = ''
    for lep in leptonList:
      for cat in catList:
        print "\n**\t Making limits for M=",m,"  cat = ",cat, "\t**\n"

        #card  = myDir+"/datacard_hzg_eeg_cat12_8TeV_"+m[:3]+".txt"
        #card  = myDir+"/mu_plus_el_M"+m+".txt"
        card  = myDir+"/output_cards/hzg_"+lep+"_2012_cat"+cat+"_M"+m+"_Dalitz_.txt"

        if cat in ['m1','m2','m3','m4','m5','m6']:
          cardNames = cardNames+' '+card
        if cat in ['EB','EE','mll50']:
          cardNames = cardNames+' '+card


        print 'Using CARD = \n',card

        if opt.noSyst:
          os.system('combine -M '+method+' '+card + ' -m '+m+' -S 0')
        else:
          print
          if not opt.dont:
            os.system('combine -M '+method+' '+card + ' -m '+m)

        if m[-1]=='5':
          fname = "higgsCombineTest."+method+".mH"+m
          os.system("mv "+fname+".root "+myDir+'/'+fname+'_cat_'+cat+'_'+lep+'.root')
        else:
          fname = "higgsCombineTest."+method+".mH"+m[0:3]
          os.system("mv "+fname+".root "+myDir+'/'+fname+'.0_cat_'+cat+'_'+lep+'.root')

      if opt.comb:
        print '\t Combining the cards: \n', cardNames
        comboName = myDir+'/output_cards/combo_mH'+m+'.txt'
        os.system('combineCards.py '+cardNames+' > '+comboName)
        os.system("sed -i 's:"+myDir+"/output_cards/"+s[:3]+":"+s[:3]+":g' "+comboName)

        if opt.noSyst:
          os.system('combine -M '+method+' '+comboName + ' -m '+m+' -S 0')
        else:
          print
          if not opt.dont:
            os.system('combine -M '+method+' '+comboName + ' -m '+m)

        if m[-1]=='5':
          fname = "higgsCombineTest."+method+".mH"+m
          os.system("mv "+fname+".root "+myDir+'/'+fname+'_cat_Combo.root')
        else:
          fname = "higgsCombineTest."+method+".mH"+m[0:3]
          os.system("mv "+fname+".root "+myDir+'/'+fname+'.0_cat_Combo.root')


      print myDir, s[:3]

  if opt.plot:
    for cat in catList:
      os.system("./limitPlotter.py "+myDir+" --cat="+cat+' --lep='+lep)
    if opt.comb:
      os.system("./limitPlotter.py "+myDir+" --cat=comb")
      # os.system('combine -M ProfileLikelihood '+sys.argv[1]+ ' -m 125 -t 100')
