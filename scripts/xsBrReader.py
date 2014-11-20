#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv

#from array import *
class AutoVivification(dict):
  """Implementation of perl's autovivification feature."""
  def __getitem__(self, item):
    try:
      return dict.__getitem__(self, item)
    except KeyError:
      value = self[item] = type(self)()
      return value


xsDict = AutoVivification()
xsScaleErrDict = AutoVivification()
xsPDFErrDict = AutoVivification()

YR = 'YR3'
TeV = '8TeV'
xsPrecision = '%.4f'
spamReader = csv.reader(open('../data/Higgs_XSBR_'+YR+'_SM_'+TeV+'.csv', 'rb'), delimiter=',')

for i, row in enumerate(spamReader):
  #if i<10: print i, row[0],row[1],row[8],row[15],row[24],row[32]
  #if i<10: print i, row[0],row[1],row[2],row[3],row[4],row[5]
  #if i<10: print i, row[0],row[8],row[9],row[10],row[11],row[12]
  #if i<10: print i, row[0],row[15],row[16],row[17],row[18],row[19]
  #if i<10: print i, row[0],row[24],row[25],row[26],row[27],row[28]
  if i < 5: continue
  if i > 289: break
  mass = '%.1f'% float(row[0])
  xsDict[YR][TeV]['ggF'][mass] = xsPrecision % float(row[1])
  xsScaleErrDict[YR][TeV]['ggF'][mass] = '%.3f/%.3f' % (1-0.01*abs(float(row[3])),1+0.01*abs(float(row[2])))
  xsPDFErrDict[YR][TeV]['ggF'][mass]   = '%.3f/%.3f' % (1-0.01*abs(float(row[5])),1+0.01*abs(float(row[4])))

  xsDict[YR][TeV]['VBF'][mass] = xsPrecision % float(row[8])
  xsScaleErrDict[YR][TeV]['VBF'][mass] = '%.3f/%.3f' % (1-0.01*abs(float(row[10])),1+0.01*abs(float(row[9])))
  xsPDFErrDict[YR][TeV]['VBF'][mass]   = '%.3f/%.3f' % (1-0.01*abs(float(row[12])),1+0.01*abs(float(row[11])))

  xsDict[YR][TeV]['WH'][mass]  = xsPrecision % float(row[15])
  xsScaleErrDict[YR][TeV]['WH'][mass] = '%.3f/%.3f' % (1-0.01*abs(float(row[17])),1+0.01*abs(float(row[16])))
  xsPDFErrDict[YR][TeV]['WH'][mass]   = '%.3f/%.3f' % (1-0.01*abs(float(row[19])),1+0.01*abs(float(row[18])))

  xsDict[YR][TeV]['ZH'][mass]  = xsPrecision % float(row[24])
  xsScaleErrDict[YR][TeV]['ZH'][mass] = '%.3f/%.3f' % (1-0.01*abs(float(row[26])),1+0.01*abs(float(row[25])))
  xsPDFErrDict[YR][TeV]['ZH'][mass]   = '%.3f/%.3f' % (1-0.01*abs(float(row[28])),1+0.01*abs(float(row[27])))

  xsDict[YR][TeV]['ttH'][mass] = xsPrecision % float(row[32])

print 'The end of CSV file reader'

sig = 'WH'
m = '125.0'
print 'xs =    ', xsDict[YR][TeV][sig][m]
print 'scale = ', xsScaleErrDict[YR][TeV][sig][m]
print 'pdf =   ', xsPDFErrDict[YR][TeV][sig][m]
