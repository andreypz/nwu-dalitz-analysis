#!/usr/bin/env python

import sys,os

for m in [125,200]:
    print "Doing mass", m
    for j in [0,1,2]:
        print "doing jet bin",j
        os.system('root -l -b -q \'TrainMva.C("BDT", "hzz'+str(m)+'", "allBg", '+str(j)+', 0,0,0,0, "weights")\' > hzz'+str(m)+'_allBg_'+str(j)+'j.log')
        print"\n"
