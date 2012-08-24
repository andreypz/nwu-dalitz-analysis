#!/usr/bin/env python

import sys,os

'''
print "do bin", 0
os.system('root -l -b -q \'TrainMva.C("BDT", "hzz125", "allBg", 0, 0,0,0,0, "weights")\' > hzz125_allBg_0j.log')
print "do bin", 1
os.system('root -l -b -q \'TrainMva.C("BDT", "hzz125", "allBg", 1, 0,0,0,0, "weights")\' > hzz125_allBg_1j.log')
print "do bin", 2
os.system('root -l -b -q \'TrainMva.C("BDT", "hzz125", "allBg", 2, 0,0,0,0, "weights")\' > hzz125_allBg_2j.log')
'''
os.system('root -l -b -q \'TrainMva.C("BDT", "hzz250", "allBg", 0, 0,0,0,0, "weights")\' > hzz250_allBg_0j.log')
os.system('root -l -b -q \'TrainMva.C("BDT", "hzz250", "allBg", 1, 0,0,0,0, "weights")\' > hzz250_allBg_1j.log')
os.system('root -l -b -q \'TrainMva.C("BDT", "hzz250", "allBg", 2, 0,0,0,0, "weights")\' > hzz250_allBg_2j.log')
