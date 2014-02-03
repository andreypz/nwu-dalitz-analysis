#!/usr/bin/env python
import sys
from ROOT import *
gROOT.SetBatch()
from toyStructs import makeToyStucts
sys.path.append("../zgamma")
import utils as u
import ConfigParser as cp
cf = cp.ConfigParser()
cf.read('config.cfg')
s = cf.get("fits","ver")

gSystem.SetIncludePath( "-I$ROOFITSYS/include/" );
gROOT.ProcessLine('.x RooStepBernstein.cxx+')
gROOT.ProcessLine('.x RooGaussStepBernstein.cxx+')

verbose = 0
testFuncs = ['GaussBern4', 'GaussBern5', 'GaussBern6']
colors = [kRed, kOrange, kGreen]
plotBase = '/uscms_data/d2/andreypz/html/zgamma/dalitz/fits-'+s+'/toyFits/'
u.createDir(plotBase)
u.createDir(s+'/biasToys/')
  
def doBiasStudy(year = '2012', lepton = 'mu', cat = '0', genFunc = 'GaussExp', mass = '125', trials = 100, job = 0, plotEvery = 5):
  #get all the starting objects
  c = TCanvas("c","c",0,0,500,400)
  c.cd()

  rooWsFile = TFile(s+'/testRooFitOut_Dalitz.root','r')
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

  #testFuncs.append(genFunc)
  print  testFuncs

  for func in testFuncs:
    fitName = '_'.join([func,year,lepton,'cat'+cat])
    if func==genFunc:
      testPdfs.append(genFit)
    else:
      testPdfs.append(myWs.pdf(fitName))
    if verbose:
      print fitName
      testPdfs[-1].Print()
      
    testBkgNorms.append(RooRealVar(   'norm'+fitName,'norm'   +fitName,bkgInSigWin,0,3*bkgInSigWin))
    testPdfs_ext.append(RooExtendPdf(  'ext'+fitName,'ext'    +fitName,testPdfs[-1],testBkgNorms[-1],sigRangeName))
    testSigNorms.append(RooRealVar('normSig'+fitName,'normSig'+fitName,0,-30,30))
    testSig_ext.append(RooExtendPdf('extSig'+fitName,'ext'    +fitName,sig,testSigNorms[-1]))

    #testFrame = mzg.frame()
    #testPdfs_ext[-1].plotOn(testFrame)
    #testSig_ext[-1].plotOn(testFrame)
    ##testFrame.Draw()
    #c.SaveAs(s+'/pdfs.png')
    
    if func==genFunc:
      testModels.append(testPdfs_ext[-1])
    else:
      testModels.append(RooAddPdf('model'+fitName,'model'+fitName,RooArgList(testSig_ext[-1],testPdfs_ext[-1])))

    if verbose:
      print 'model'+fitName
      testModels[-1].Print()

  testModelsDict   = dict(zip(testFuncs,testModels))
  testBkgNormsDict = dict(zip(testFuncs,testBkgNorms))
  testSigNormsDict = dict(zip(testFuncs,testSigNorms))
  ColorDict = dict(zip(testFuncs,colors))


  # prep the outputs
  outName = s+'/biasToys/'+'_'.join(['biasToys',year,lepton,'cat'+cat,genFunc,mass,'job'+str(job)])+'.root'
  outFile = TFile(outName, 'RECREATE')
  tree = TTree('toys','toys')
  makeToyStucts()


  #set up branches
  from ROOT import TOYDATA
  toyDataStruct = TOYDATA()
  tree.Branch('toyData', toyDataStruct, 'totalData/I:sigWindowData')
  from ROOT import GEN
  genStruct = GEN()
  tree.Branch('gen',genStruct, 'yieldBkg/D:yieldBkgErr:yieldSig:yieldSigErr:edm:minNll:statusAll/I:statusMIGRAD:statusHESSE:covQual:numInvalidNLL')

  truthStruct = TRUTH()
  tree.Branch('truth',truthStruct, 'yieldBkg/D:yieldBkgErr:yieldSig:yieldSigErr')

  for func in testFuncs:
    print func
    
    if func=='GaussBern3':
      from ROOT import GAUSSBERN3
      GaussBern3Struct = GAUSSBERN3()
      tree.Branch('GaussBern3',GaussBern3Struct,'yieldBkg/D:yieldBkgErr:yieldSig:yieldSigErr:paramSigma:paramSigmaErr:paramStep:paramStepErr:paramP1:paramP1Err:paramP2:paramP2Err:paramP3:paramP3Err:edm:minNll:statusAll/I:statusMIGRAD:statusHESSE:covQual:numInvalidNLL')
    if func=='GaussBern4':
      from ROOT import GAUSSBERN4
      GaussBern4Struct = GAUSSBERN4()
      tree.Branch('GaussBern4',GaussBern4Struct,'yieldBkg/D:yieldBkgErr:yieldSig:yieldSigErr:paramSigma:paramSigmaErr:paramStep:paramStepErr:paramP1:paramP1Err:paramP2:paramP2Err:paramP3:paramP3Err:paramP4:paramP4Err:edm:minNll:statusAll/I:statusMIGRAD:statusHESSE:covQual:numInvalidNLL')
    if func=='GaussBern5':
      from ROOT import GAUSSBERN5
      GaussBern5Struct = GAUSSBERN5()
      tree.Branch('GaussBern5',GaussBern5Struct,'yieldBkg/D:yieldBkgErr:yieldSig:yieldSigErr:paramSigma:paramSigmaErr:paramStep:paramStepErr:paramP1:paramP1Err:paramP2:paramP2Err:paramP3:paramP3Err:paramP4:paramP4Err:paramP5:paramP5Err:edm:minNll:statusAll/I:statusMIGRAD:statusHESSE:covQual:numInvalidNLL')

      
    if func=='GaussBern6':
      # for now just use the same as Bern5, fix it later
      from ROOT import GAUSSBERN5
      GaussBern6Struct = GAUSSBERN5()
      tree.Branch('GaussBern6',GaussBern6Struct,'yieldBkg/D:yieldBkgErr:yieldSig:yieldSigErr:paramSigma:paramSigmaErr:paramStep:paramStepErr:paramP1:paramP1Err:paramP2:paramP2Err:paramP3:paramP3Err:paramP4:paramP4Err:paramP5:paramP5Err:edm:minNll:statusAll/I:statusMIGRAD:statusHESSE:covQual:numInvalidNLL')
      

  structDict = dict(zip(testFuncs,[GaussBern4Struct,GaussBern5Struct,GaussBern6Struct]))
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

        ##print "RAW INPUT"
        #raw_input()

      suffix = '_'.join([func,year,lepton,'cat'+cat])

      structDict[func].yieldBkg    = testBkgNormsDict[func].getVal()
      structDict[func].yieldBkgErr = testBkgNormsDict[func].getError()
      structDict[func].yieldSig    = testSigNormsDict[func].getVal()
      structDict[func].yieldSigErr = testSigNormsDict[func].getError()
      if func in ['GaussBern3','GaussBern4','GaussBern5','GaussBern6']:
        structDict[func].paramSigma    = testModelsDict[func].getParameters(toyData)['sigma'+suffix].getVal()
        structDict[func].paramSigmaErr = testModelsDict[func].getParameters(toyData)['sigma'+suffix].getError()
        structDict[func].paramStep     = testModelsDict[func].getParameters(toyData)['step'+suffix].getVal()
        structDict[func].paramStepErr  = testModelsDict[func].getParameters(toyData)['step'+suffix].getError()
        structDict[func].paramP1       = testModelsDict[func].getParameters(toyData)['p1'+suffix].getVal()
        structDict[func].paramP1Err    = testModelsDict[func].getParameters(toyData)['p1'+suffix].getError()
        structDict[func].paramP2       = testModelsDict[func].getParameters(toyData)['p2'+suffix].getVal()
        structDict[func].paramP2Err    = testModelsDict[func].getParameters(toyData)['p2'+suffix].getError()
        structDict[func].paramP3       = testModelsDict[func].getParameters(toyData)['p3'+suffix].getVal()
        structDict[func].paramP3Err    = testModelsDict[func].getParameters(toyData)['p3'+suffix].getError()
      if func in ['GaussBern4','GaussBern5','GaussBern6']:
        structDict[func].paramP4    = testModelsDict[func].getParameters(toyData)['p4'+suffix].getVal()
        structDict[func].paramP4Err = testModelsDict[func].getParameters(toyData)['p4'+suffix].getError()
      if func in ['GaussBern5','GaussBern6']:
        structDict[func].paramP5    = testModelsDict[func].getParameters(toyData)['p5'+suffix].getVal()
        structDict[func].paramP5Err = testModelsDict[func].getParameters(toyData)['p5'+suffix].getError()
      if func in ['GaussBern6']:
        print "Do nothing else, add later"

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
        testModelsDict[func].plotOn(testFrame, RooFit.LineColor(ColorDict[func]), RooFit.Range('fullRange'), RooFit.Name(func))
        legendTmp.AddEntry(testFrame.findObject(func), func ,'l')  # trick to get the colors on the legend
        
      testFrame.Draw()
      legendTmp.Draw()

      c.SaveAs(plotBase+'_'.join(['toyFits',lepton,year,'cat'+cat,genFunc,'M'+mass,'job'+str(job),'trial'+str(i)])+'.png')

    tree.Fill()

  outFile.cd()
  tree.Write()
  outFile.Close()

  print 'so many toys!'

  
if __name__=="__main__":
  print len(sys.argv)
  print sys.argv
  if len(sys.argv) != 3:
    sys.exit()
    
  n = sys.argv[1]
              
  doBiasStudy(job=int(n))


