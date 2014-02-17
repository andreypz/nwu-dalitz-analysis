#!/usr/bin/env python
from optparse import OptionParser
import sys,os,re
import ConfigParser as cp
cf = cp.ConfigParser()
cf.read('config.cfg')
s = cf.get("fits","ver")
verbose = 0

from ROOT import *
from toyStructs import makeToyStucts
gSystem.SetIncludePath( "-I$ROOFITSYS/include/" );
#gROOT.ProcessLine('.x RooStepBernstein.cxx+')
#gROOT.ProcessLine('.x RooGaussStepBernstein.cxx+')

gROOT.SetBatch()

parser = OptionParser(usage="usage: %prog [options -v version]")
parser.add_option("-v", "--ver",dest="ver", default="none", help="version - a subdir name where the input files are stored.")
parser.add_option("-b", dest="bullshit", action="store_true",default=False, help="Ys it is true")
parser.add_option("-j", dest="job",    default=0,     help="Job number")
parser.add_option("-t", dest="trials", default=5,     help="Number of trials")
parser.add_option("-m", dest="mass",   default='125', help="Number of trials")

(options, args) = parser.parse_args()

colors = {}
for f,col in cf.items("colors"): colors[f] = int(col)

fileBase = './'
plotBase = './'
if options.ver=='none':
  fileBase = s
  plotBase = '/uscms_data/d2/andreypz/html/zgamma/dalitz/fits-'+s
else:
  fileBase = options.ver
  plotBase = options.ver


def doBiasStudy(year = '2012', lepton = 'mu', cat = '0', genFunc = 'Bern3', mass = '125', trials = 5, job = 0, plotEvery = 50):
  #get all the starting objects
  c = TCanvas("c","c",0,0,500,400)
  c.cd()
  testFuncs = ['Bern2','Bern3','Bern4','Bern5']

  if not os.path.isdir(plotBase):
    os.makedirs(plotBase)
    print 'bias.py: making directory',plotBase

  if verbose:
    print 'DBG bias ', fileBase
    print '  RAW INPUT HERE:'
    raw_input()
  
  rooWsFile = TFile(fileBase.replace('/','')+'-testRooFitOut_Dalitz.root','r')
  myWs = rooWsFile.Get('ws')
  sigRangeName = '_'.join(['range',lepton,year,'cat'+cat,'M'+mass])

  myWs.Print()
  # get the x-axis
  mzg = myWs.var('CMS_hzg_mass')
  if verbose:
    mzg.Print()
    print sigRangeName, mzg.getMin(sigRangeName), mzg.getMax(sigRangeName)
    #print " RAW  INPUT NEEDED HERE"
    #raw_input()
    
  # get the data
  dataName = '_'.join(['data',lepton,year,'cat'+cat])
  data = myWs.data(dataName)
  realDataYield = data.sumEntries()
  bkgInSigWin= data.sumEntries('1',sigRangeName)
  if verbose: data.Print()
  if verbose: print 'total data:', realDataYield, 'total data2:', data.numEntries(),'total data in sig window:', bkgInSigWin


  # get the gen pdf
  genFitName = '_'.join([genFunc,year,lepton,'cat'+cat])
  genFit = myWs.pdf(genFitName)
  if verbose:
    print genFitName
    genFit.Print()

  # get the signal
  sigName = '_'.join(['pdf','sig',lepton,year,'cat'+cat,'M'+mass])
  sig = myWs.pdf(sigName)
  if verbose:
    print sigName
    sig.Print()

  # get the test functions, turn them into extended pdfs with signal models.  Also include the gen function for closure tests
  testPdfs     = []
  testBkgNorms = []
  testPdfs_ext = []
  testSigNorms = []
  testSig_ext  = []
  testModels   = []

  testFuncs.append('gen')

  for n,func in enumerate(testFuncs):
    print n, 'Doing func =', func
    fitName = '_'.join([func,year,lepton,'cat'+cat])
    if func=='gen':  # the last one appended to testFuncs
      testPdfs.append(genFit)
      print 'Found it !', func
    else:
      testPdfs.append(myWs.pdf(fitName))
    if verbose:
      print fitName
      testPdfs[-1].Print()
      
    testBkgNorms.append(RooRealVar(   'norm'+fitName,'norm'   +fitName,bkgInSigWin,0,5*bkgInSigWin))
    testPdfs_ext.append(RooExtendPdf(  'ext'+fitName,'ext'    +fitName,testPdfs[-1],testBkgNorms[-1], sigRangeName))
    testSigNorms.append(RooRealVar('normSig'+fitName,'normSig'+fitName,0,-80,80))
    testSig_ext.append(RooExtendPdf('extSig'+fitName,'ext'    +fitName,sig,testSigNorms[-1]))

    #testFrame = mzg.frame()
    #testPdfs_ext[-1].plotOn(testFrame)
    #testSig_ext[-1].plotOn(testFrame)
    ##testFrame.Draw()
    #c.SaveAs(s+'/pdfs.png')
    
    if func=='gen':
      testModels.append(testPdfs_ext[-1])
    else:
      testModels.append(RooAddPdf('model'+fitName,'model'+fitName,RooArgList(testSig_ext[-1],testPdfs_ext[-1])))

    if verbose:
      print 'model'+fitName
      testModels[-1].Print()

  testModelsDict   = dict(zip(testFuncs,testModels))
  testBkgNormsDict = dict(zip(testFuncs,testBkgNorms))
  testSigNormsDict = dict(zip(testFuncs,testSigNorms))

  # prep the outputs
  outName = fileBase+'/'+'_'.join(['biasToys',year,lepton,'cat'+cat,genFunc,mass,'job'+str(job)])+'.root'
  outFile = TFile(outName, 'RECREATE')
  tree = TTree('toys','toys')
  makeToyStucts()


  #set up branches
  from ROOT import TOYDATA
  toyDataStruct = TOYDATA()
  tree.Branch('toyData', toyDataStruct, 'totalData/I:sigWindowData')
  truthStruct = TRUTH()
  tree.Branch('truth',truthStruct, 'yieldBkg/D:yieldBkgErr:yieldSig:yieldSigErr')


  from ROOT import GEN
  for n,func in enumerate(testFuncs):
    print n,func

    if verbose:
      print "RAW INPUT"
      raw_input()


    if func=='gen':
      genStruct = GEN()
      tree.Branch('gen',genStruct, 'yieldBkg/D:yieldBkgErr:yieldSig:yieldSigErr:edm:minNll:statusAll/I:statusMIGRAD:statusHESSE:covQual:numInvalidNLL')
    elif func=='Bern2':
      from ROOT import BERN2
      Bern2Struct = BERN2()
      tree.Branch('Bern2',Bern2Struct,'yieldBkg/D:yieldBkgErr:yieldSig:yieldSigErr:paramP1:paramP1Err:paramP2:paramP2Err:edm:minNll:statusAll/I:statusMIGRAD:statusHESSE:covQual:numInvalidNLL')
    elif func=='Bern3':
      from ROOT import BERN3
      Bern3Struct = BERN3()
      tree.Branch('Bern3',Bern3Struct,'yieldBkg/D:yieldBkgErr:yieldSig:yieldSigErr:paramP1:paramP1Err:paramP2:paramP2Err:paramP3:paramP3Err:edm:minNll:statusAll/I:statusMIGRAD:statusHESSE:covQual:numInvalidNLL')
    elif func=='Bern4':
      from ROOT import BERN4
      Bern4Struct = BERN4()
      tree.Branch('Bern4',Bern4Struct,'yieldBkg/D:yieldBkgErr:yieldSig:yieldSigErr:paramP1:paramP1Err:paramP2:paramP2Err:paramP3:paramP3Err:paramP4:paramP4Err:edm:minNll:statusAll/I:statusMIGRAD:statusHESSE:covQual:numInvalidNLL')
    elif func=='Bern5':
      from ROOT import BERN5
      Bern5Struct = BERN5()
      tree.Branch('Bern5',Bern5Struct,'yieldBkg/D:yieldBkgErr:yieldSig:yieldSigErr:paramP1:paramP1Err:paramP2:paramP2Err:paramP3:paramP3Err:paramP4:paramP4Err:paramP5:paramP5Err:edm:minNll:statusAll/I:statusMIGRAD:statusHESSE:covQual:numInvalidNLL')
    elif func=='Bern6':
      from ROOT import BERN6
      Bern6Struct = BERN6()
      tree.Branch('Bern6',Bern6Struct,'yieldBkg/D:yieldBkgErr:yieldSig:yieldSigErr:paramP1:paramP1Err:paramP2:paramP2Err:paramP3:paramP3Err:paramP4:paramP4Err:paramP5:paramP5Err:paramP6:paramP6Err:edm:minNll:statusAll/I:statusMIGRAD:statusHESSE:covQual:numInvalidNLL')

    elif func=='GaussBern3':
      from ROOT import GAUSSBERN3
      GaussBern3Struct = GAUSSBERN3()
      tree.Branch('GaussBern3',GaussBern3Struct,'yieldBkg/D:yieldBkgErr:yieldSig:yieldSigErr:paramSigma:paramSigmaErr:paramStep:paramStepErr:paramP1:paramP1Err:paramP2:paramP2Err:paramP3:paramP3Err:edm:minNll:statusAll/I:statusMIGRAD:statusHESSE:covQual:numInvalidNLL')
    elif func=='GaussBern4':
      from ROOT import GAUSSBERN4
      GaussBern4Struct = GAUSSBERN4()
      tree.Branch('GaussBern4',GaussBern4Struct,'yieldBkg/D:yieldBkgErr:yieldSig:yieldSigErr:paramSigma:paramSigmaErr:paramStep:paramStepErr:paramP1:paramP1Err:paramP2:paramP2Err:paramP3:paramP3Err:paramP4:paramP4Err:edm:minNll:statusAll/I:statusMIGRAD:statusHESSE:covQual:numInvalidNLL')
    elif func=='GaussBern5':
      from ROOT import GAUSSBERN5
      GaussBern5Struct = GAUSSBERN5()
      tree.Branch('GaussBern5',GaussBern5Struct,'yieldBkg/D:yieldBkgErr:yieldSig:yieldSigErr:paramSigma:paramSigmaErr:paramStep:paramStepErr:paramP1:paramP1Err:paramP2:paramP2Err:paramP3:paramP3Err:paramP4:paramP4Err:paramP5:paramP5Err:edm:minNll:statusAll/I:statusMIGRAD:statusHESSE:covQual:numInvalidNLL')
    elif func=='GaussBern6':
      # for now just use the same as Bern5, fix it later
      from ROOT import GAUSSBERN6
      GaussBern6Struct = GAUSSBERN6()
      tree.Branch('GaussBern6',GaussBern6Struct,'yieldBkg/D:yieldBkgErr:yieldSig:yieldSigErr:paramSigma:paramSigmaErr:paramStep:paramStepErr:paramP1:paramP1Err:paramP2:paramP2Err:paramP3:paramP3Err:paramP4:paramP4Err:paramP5:paramP5Err:paramP6:paramP6Err:edm:minNll:statusAll/I:statusMIGRAD:statusHESSE:covQual:numInvalidNLL')
    else:
      print "No struct exists for ", func

  structDict = dict(zip(testFuncs,[Bern2Struct,Bern3Struct,Bern4Struct,Bern5Struct, genStruct]))
  #structDict = dict(zip(testFuncs,[GaussBern4Struct,GaussBern5Struct,GaussBern6Struct]))
  print structDict
  #print "  RAW INPUT HERE "
  #raw_input()
  
  r = TRandom3(1024+31*job)
  RooRandom.randomGenerator().SetSeed(1024+31*job)


  ############################
  # time to throw some toys! #
  ############################

  for i in range(0,trials):
    print 'doing trial:',i

    genBkgYield = r.Poisson(data.numEntries())
    toyData = genFit.generate(RooArgSet(mzg),genBkgYield)
    #toyData = genFit.generate(RooArgSet(mzg), RooFit.Extended(kTRUE))
    bkg_est = toyData.sumEntries('1',sigRangeName)
    if verbose: print 'bkg_est',bkg_est
    #print "RAW INPUT"
    #raw_input()
    truthStruct.yieldBkg = bkg_est
    truthStruct.yieldBkgErr = sqrt(bkg_est)
    truthStruct.yieldSig = 0
    truthStruct.yieldSigErr = 0
    
    for func in testFuncs:
      testSigNormsDict[func].setVal(0)
      testBkgNormsDict[func].setVal(bkg_est)

      nll = testModelsDict[func].createNLL(toyData,RooFit.Extended())
      m = RooMinuit(nll)
      m.migrad()
      resMigrad = m.save()
      m.hesse()
      resHesse = m.save()

      res = testModelsDict[func].fitTo(toyData,RooFit.Save(),RooFit.PrintLevel(-1))

      statusAll    = res.status()
      statusMIGRAD = resMigrad.status()
      statusHESSE  = resHesse.status()
      numInvalidNLL = res.numInvalidNLL()
      edm     = res.edm()
      minNll  = res.minNll()
      covQual = res.covQual()

      if verbose:
        print 'statusAll', statusAll
        print 'statusMIGRAD', statusMIGRAD
        print 'statusHESSE', statusHESSE
        print 'numInvalidNLL', numInvalidNLL
        print 'edm', edm
        print 'minNll', minNll
        print 'covQual', covQual
        print 'yieldBkg', testBkgNormsDict[func].getVal()
        print 'yieldSig', testSigNormsDict[func].getVal()
        testModelsDict[func].getParameters(toyData).Print('v')

      suffix = '_'.join([func,year,lepton,'cat'+cat])

      structDict[func].yieldBkg    = testBkgNormsDict[func].getVal()
      structDict[func].yieldBkgErr = testBkgNormsDict[func].getError()
      structDict[func].yieldSig    = testSigNormsDict[func].getVal()
      structDict[func].yieldSigErr = testSigNormsDict[func].getError()

      if func in ['Bern2','Bern3','Bern4','Bern5','Bern6','GaussBern3','GaussBern4','GaussBern5','GaussBern6']:
        structDict[func].paramP1       = testModelsDict[func].getParameters(toyData)['p1'+suffix].getVal()
        structDict[func].paramP1Err    = testModelsDict[func].getParameters(toyData)['p1'+suffix].getError()
        structDict[func].paramP2       = testModelsDict[func].getParameters(toyData)['p2'+suffix].getVal()
        structDict[func].paramP2Err    = testModelsDict[func].getParameters(toyData)['p2'+suffix].getError()
      if func in ['Bern3','Bern4','Bern5','Bern6','GaussBern3','GaussBern4','GaussBern5','GaussBern6']:
        structDict[func].paramP3       = testModelsDict[func].getParameters(toyData)['p3'+suffix].getVal()
        structDict[func].paramP3Err    = testModelsDict[func].getParameters(toyData)['p3'+suffix].getError()
      if func in ['Bern4','Bern5','Bern6','GaussBern4','GaussBern5','GaussBern6']:
        structDict[func].paramP4    = testModelsDict[func].getParameters(toyData)['p4'+suffix].getVal()
        structDict[func].paramP4Err = testModelsDict[func].getParameters(toyData)['p4'+suffix].getError()
      if func in ['Bern5','Bern6','GaussBern5','GaussBern6']:
        structDict[func].paramP5    = testModelsDict[func].getParameters(toyData)['p5'+suffix].getVal()
        structDict[func].paramP5Err = testModelsDict[func].getParameters(toyData)['p5'+suffix].getError()
      if func in ['Bern6','GaussBern6']:
        #print "do nothing for now"
        structDict[func].paramP6    = testModelsDict[func].getParameters(toyData)['p6'+suffix].getVal()
        structDict[func].paramP6Err = testModelsDict[func].getParameters(toyData)['p6'+suffix].getError()
      if func in ['GaussBern3','GaussBern4','GaussBern5','GaussBern6']:
          structDict[func].paramSigma    = testModelsDict[func].getParameters(toyData)['sigma'+suffix].getVal()
          structDict[func].paramSigmaErr = testModelsDict[func].getParameters(toyData)['sigma'+suffix].getError()
          structDict[func].paramStep     = testModelsDict[func].getParameters(toyData)['step'+suffix].getVal()
          structDict[func].paramStepErr  = testModelsDict[func].getParameters(toyData)['step'+suffix].getError()


      structDict[func].statusAll     = statusAll
      structDict[func].statusMIGRAD  = statusMIGRAD
      structDict[func].statusHESSE   = statusHESSE
      structDict[func].numInvalidNLL = numInvalidNLL
      structDict[func].edm     = edm
      structDict[func].minNll  = minNll
      structDict[func].covQual = covQual



    toyDataStruct.totalData = toyData.numEntries()
    toyDataStruct.sigWindowData = bkg_est

    if i%plotEvery==0:
      testFrame = mzg.frame()
      toyData.plotOn(testFrame)
      testFrame.SetTitle(";m_{H} (GeV);Events/1 GeV")
      genFit.plotOn(testFrame, RooFit.Name(genFunc))
                                    
      legendTmp = TLegend(0.6,0.65,0.9,0.9)
      legendTmp.SetFillColor(0)
      #legendTmp.SetFillStyle(0)
      legendTmp.SetTextSize(0.03)

      legendTmp.AddEntry(testFrame.findObject(genFunc), genFunc+': used for GEN toys' ,'l')
      for func in testFuncs:
        testModelsDict[func].plotOn(testFrame, RooFit.LineColor(colors[func.lower()]), RooFit.Range('fullRange'), RooFit.Name(func))
        legendTmp.AddEntry(testFrame.findObject(func), func ,'l')  # trick to get the colors on the legend
        
      testFrame.Draw()
      legendTmp.Draw()

      c.SaveAs(plotBase+'/'+'_'.join(['toyFits',lepton,year,'cat'+cat,genFunc,'M'+mass,'job'+str(job),'trial'+str(i)])+'.png')

    tree.Fill()

  outFile.cd()
  tree.Write()
  outFile.Close()

  print 'so many toys!'
  
  
if __name__=="__main__":
  print len(sys.argv), sys.argv
  
  job = int(options.job)
  trials = int(options.trials)
  mass = str(options.mass)
  
  for f in ['Exp','Pow','Bern3']:
    print 'Starting'
    doBiasStudy(trials = trials, genFunc=f, job=job, mass=mass)
