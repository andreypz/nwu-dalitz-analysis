#!/usr/bin/env python

import sys,os

os.system('root -l -b -q \'TrainMva.C("BDT", "hzz250", "allBg", 0, 0,0,0,0, "weights")\' > hzz250_allBg_0j.log')
#os.system('root -l -b -q \'TrainMva.C("BDT", "hzz250_all", "allBgPhoton", 1, 0,0,0,0, "weights")\' > hzz250_allBgPhoton_1j.log')
#os.system('root -l -b -q \'TrainMva.C("BDT", "hzz250_all", "allBgPhoton", 2, 0,0,0,0, "weights")\' > hzz250_allBgPhoton_2j.log')
