// Modified from TMVAClassification.C 
// last  version with variables specific for HZZ
//                      Dec 05, 2011 -  Anton A.

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

// temporary comment out for batch jobs
// to launch it at the end of the job, put back the call at the end of this file
//#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif


TString inputFilesDir = "mvaInputs_v69/muon/";
TString outputWeightsDir = "";



Bool_t USE_SEPARATE_TRAIN_TEST_FILES = true;
// if selecting training and testing events from the same file
// one has to enter specify the number of events

void TrainMva(TString myMethodList = "",
        TString _signalName,
        TString _bgName,
        Int_t _jetMulti,
        Int_t _numSignalTrain = 0,
        Int_t _numBgTrain = 0,
        Int_t _numSignalTest = 0,
        Int_t _numBgTest = 0,
        TString _weightsSubDir = "") {


    std::cout << "==> Starting TrainMva " << std::endl;
    // subdirectory where output weights for classifiers will go
    outputWeightsDir += _weightsSubDir;


    // Common weights for the entire sample
    // use with event weights that do not account for x-section,
    // or with fullWeight to avoid too small weights due to small sample scale factors (if such issues occur)
    Double_t signalWeight = 1.0;
    Double_t backgroundWeight = 1.0;

    TString sigFileName_train;
    TString sigFileName_test;

    TString bgFileName_train;
    TString bgFileName_test;

    // when using common file for training and testing
    TString sigFileName;
    TString bgFileName;

    // signal names are hzzXXX
    // bgName: zz, ttbar, dymm,.... see below
    // this is for individual training. 
    // When bg's are combined we use a variation of  allBg.
    // For the case when the gamma+jet sample was used for DY
    // the name is allBgPhoton.
    //  Define bg's inputs as needed

    // CHECK THE FILENAMES BELOW AGAIN!!!

    if (_bgName == "zz") {
        bgFileName_train = inputFilesDir + "zz_train.root";
        bgFileName_test = inputFilesDir + "zz_test.root";
        bgFileName = inputFilesDir + "zz.root"; // when it is common.
    } else

        if (_bgName == "ttbar") {
        bgFileName_train = inputFilesDir + "ttbar_train.root";
        bgFileName_test = inputFilesDir + "ttbar_test.root";
        bgFileName = inputFilesDir + "ttbar.root"; // when it is common. 
    } else

        if (_bgName == "dymm") {
        bgFileName_train = inputFilesDir + "dymm_train.root";
        bgFileName_test = inputFilesDir + "dymm_test.root";
        bgFileName = inputFilesDir + "dymm.root"; // when it is common.
    } else

        if (_bgName == "dyee") {
        bgFileName_train = inputFilesDir + "dyee_train.root";
        bgFileName_test = inputFilesDir + "dyee_test.root";
        bgFileName = inputFilesDir + "dyee.root"; // when it is common.
    } else

        if (_bgName == "dytt") {
        bgFileName_train = inputFilesDir + "dytt_train.root";
        bgFileName_test = inputFilesDir + "dytt_test.root";
        bgFileName = inputFilesDir + "dytt.root"; // when it is common.
    } else

        if (_bgName == "allBg") {
        bgFileName_train = inputFilesDir + "allBg_train.root";
        bgFileName_test = inputFilesDir + "allBg_test.root";
        bgFileName = inputFilesDir + "allBg.root"; // when it is common.
    } else

        if (_bgName == "allBg2") {
        bgFileName_train = inputFilesDir + "allBg2_train.root";
        bgFileName_test = inputFilesDir + "allBg2_test.root";
        bgFileName = inputFilesDir + "allBg2.root"; // when it is common.
    } else
        if (_bgName == "allBgPhoton") {
        bgFileName_train = inputFilesDir + "allBgPhoton_train.root";
        bgFileName_test = inputFilesDir + "allBgPhoton_test.root";
        bgFileName = inputFilesDir + "allBgPhoton.root"; // when it is common.
    } else {
        std::cout << "\n\nUnknown background \"" << _bgName << "\". Check Input!\n\n";
        return;
    }


    if (_signalName == "hzz125") {
        sigFileName_train = inputFilesDir + "hzz125_train.root";
        sigFileName_test = inputFilesDir + "hzz125_test.root";
        sigFileName = inputFilesDir + "hzz125.root";
    } else
    if (_signalName == "hzz200") {
        sigFileName_train = inputFilesDir + "hzz200_train.root";
        sigFileName_test = inputFilesDir + "hzz200_test.root";
        sigFileName = inputFilesDir + "hzz200.root";
    } else

        if (_signalName == "hzz250") {
        sigFileName_train = inputFilesDir + "hzz250_train.root"; /// TMP (AA) select gg explicitly for testing
        sigFileName_test = inputFilesDir + "hzz250_test.root";
        sigFileName = inputFilesDir + "hzz250.root";
    } else

        if (_signalName == "hzz300") {
        sigFileName_train = inputFilesDir + "hzz300_train.root";
        sigFileName_test = inputFilesDir + "hzz300_test.root";
        sigFileName = inputFilesDir + "hzz300.root";
    } else

        if (_signalName == "hzz350") {
        sigFileName_train = inputFilesDir + "hzz350_train.root";
        sigFileName_test = inputFilesDir + "hzz350_test.root";
        sigFileName = inputFilesDir + "hzz350.root";
    } else

        if (_signalName == "hzz400") {
        sigFileName_train = inputFilesDir + "hzz400_train.root";
        sigFileName_test = inputFilesDir + "hzz400_test.root";
        sigFileName = inputFilesDir + "hzz400.root";
        ////        signalWeight = 0.000284 * 1.1;
    } else

        if (_signalName == "hzz450") {
        sigFileName_train = inputFilesDir + "hzz450_train.root";
        sigFileName_test = inputFilesDir + "hzz450_test.root";
        sigFileName = inputFilesDir + "hzz450.root";
    } else

        if (_signalName == "hzz500") {
        sigFileName_train = inputFilesDir + "hzz500_train.root";
        sigFileName_test = inputFilesDir + "hzz500_test.root";
        sigFileName = inputFilesDir + "hzz500.root";
    } else

        if (_signalName == "hzz550") {
        sigFileName_train = inputFilesDir + "hzz550_train.root";
        sigFileName_test = inputFilesDir + "hzz550_test.root";
        sigFileName = inputFilesDir + "hzz550.root";
    } else

        if (_signalName == "hzz600") {
        sigFileName_train = inputFilesDir + "hzz600_train.root";
        sigFileName_test = inputFilesDir + "hzz600_test.root";
        sigFileName = inputFilesDir + "hzz600.root";

    }// this for the case when the gg+VBF is used in the training

	else


        if (_signalName == "hzz125_all") {
        sigFileName_train = inputFilesDir + "hzz125_all_train.root";
        sigFileName_test = inputFilesDir + "hzz125_all_test.root";
        sigFileName = inputFilesDir + "hzz125_all.root";
    } else

        if (_signalName == "hzz200_all") {
        sigFileName_train = inputFilesDir + "hzz200_all_train.root";
        sigFileName_test = inputFilesDir + "hzz200_all_test.root";
        sigFileName = inputFilesDir + "hzz200_all.root";
    } else

        if (_signalName == "hzz250_all")) {
            sigFileName_train = inputFilesDir + "hzz250_all_train.root";
            sigFileName_test = inputFilesDir + "hzz250_all_test.root";
            sigFileName = inputFilesDir + "hzz250_all.root";
        } else

        if (_signalName == "hzz300_all")) {
            sigFileName_train = inputFilesDir + "hzz300_all_train.root";
            sigFileName_test = inputFilesDir + "hzz300_all_test.root";
            sigFileName = inputFilesDir + "hzz300_all.root";
        } else

        if (_signalName == "hzz350_all")) {
            sigFileName_train = inputFilesDir + "hzz350_all_train.root";
            sigFileName_test = inputFilesDir + "hzz350_all_test.root";
            sigFileName = inputFilesDir + "hzz350_all.root";
        } else

        if (_signalName == "hzz400_all") {
        sigFileName_train = inputFilesDir + "hzz400_all_train.root";
        sigFileName_test = inputFilesDir + "hzz400_all_test.root";
        sigFileName = inputFilesDir + "hzz400_all.root";
    } else

        if (_signalName == "hzz450_all") {
        sigFileName_train = inputFilesDir + "hzz450_all_train.root";
        sigFileName_test = inputFilesDir + "hzz450_all_test.root";
        sigFileName = inputFilesDir + "hzz450_all.root";
    } else

        if (_signalName == "hzz500_all") {
        sigFileName_train = inputFilesDir + "hzz500_all_train.root";
        sigFileName_test = inputFilesDir + "hzz500_all_test.root";
        sigFileName = inputFilesDir + "hzz500_all.root";
    } else

        if (_signalName == "hzz550_all") {
        sigFileName_train = inputFilesDir + "hzz550_all_train.root";
        sigFileName_test = inputFilesDir + "hzz550_all_test.root";
        sigFileName = inputFilesDir + "hzz550_all.root";
    } else

        if (_signalName == "hzz600_all") {
        sigFileName_train = inputFilesDir + "hzz600_all_train.root";
        sigFileName_test = inputFilesDir + "hzz600_all_test.root";
        sigFileName = inputFilesDir + "hzz600_all.root";

    } else {
        std::cout << "\n\nUnknown signal sample \"" << _signalName << "\". Check Input!\n\n";
        return;
    }


    // CUTS applied on top of the pre-selections in the input ntuples
    // If there there is different MVA selector for different jet bins the 
    // cuts should reflect the selections.

    TString jetMultiName = "0j"; //used for inclusion in file names   
    TCut addCutsSig = "nJets==0";
    TCut addCutsBg = "nJets==0";

    if (_jetMulti == 1) {
        jetMultiName = "1j";
        addCutsSig = "nJets==1";
        addCutsBg = "nJets==1";
    } else if (_jetMulti > 1) {
        jetMultiName = "gt1j";
        addCutsSig = "nJets>1";
        addCutsBg = "nJets>1";
    }

    //   used for the filenames of the MVA weights xml files, etc.
    TString sampleNames = TString::Format("%s_%s_%s_", _bgName.Data(), _signalName.Data(), jetMultiName.Data());

    // contains the performance histograms from the training
    // and the input variables
    TString outFileName = sampleNames + "TMVA.root";

    //-----------------------------------------------------------

    TMVA::Tools::Instance();

    // Default MVA methods to be trained + tested
    std::map<std::string, int> Use;

    // See available methods below. 
    // If the list of methods passed to TrainMva is empty (=""), the switches below determine what methods are used
    // else they are ignored.
    
    // --- Cut optimisation
    Use["Cuts"] = 0; 
    Use["CutsD"] = 0;
    Use["CutsPCA"] = 0;
    Use["CutsGA"] = 0;
    Use["CutsSA"] = 0;
    // 
    // --- 1-dimensional likelihood ("naive Bayes estimator")
    Use["Likelihood"] = 0; 
    Use["LikelihoodD"] = 0; 
    Use["LikelihoodPCA"] = 0;
    Use["LikelihoodKDE"] = 0;
    Use["LikelihoodMIX"] = 0;
    //
    // --- Mutidimensional likelihood and Nearest-Neighbour methods
    Use["PDERS"] = 0; 
    Use["PDERSD"] = 0;
    Use["PDERSPCA"] = 0;
    Use["PDEFoam"] = 0; 
    Use["PDEFoamBoost"] = 0; 
    Use["KNN"] = 0; 
    //
    // --- Linear Discriminant Analysis
    Use["LD"] = 0;
    Use["Fisher"] = 0;
    Use["FisherG"] = 0;
    Use["BoostedFisher"] = 0;
    Use["HMatrix"] = 0;
    //
    // --- Function Discriminant analysis
    Use["FDA_GA"] = 0; 
    Use["FDA_SA"] = 0;
    Use["FDA_MC"] = 0;
    Use["FDA_MT"] = 0;
    Use["FDA_GAMT"] = 0;
    Use["FDA_MCMT"] = 0;
    //
    // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
    Use["MLP"] = 0; // Recommended ANN
    Use["MLPBFGS"] = 0; // Recommended ANN with optional training method
    Use["MLPBNN"] = 0; // Recommended ANN with BFGS training method and bayesian regulator
    Use["CFMlpANN"] = 0; // Depreciated ANN from ALEPH
    Use["TMlpANN"] = 0; // ROOT's own ANN
    //
    // --- Support Vector Machine 
    Use["SVM"] = 0; 
    // 
    // --- Boosted Decision Trees
    Use["BDT"] = 1; // uses Adaptive Boost
    Use["BDTG"] = 0; // uses Gradient Boost
    Use["BDTB"] = 0; // uses Bagging
    Use["BDTD"] = 0; // decorrelation + Adaptive Boost
    Use["BDTF"] = 0; // allow usage of fisher discriminant for node splitting 
    // 
    // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
    Use["RuleFit"] = 0; // problem with DY (AA)
    // ---------------------------------------------------------------

   
    std::cout << std::endl;
    std::cout << "==> Start TrainMva" << std::endl;

    // Select methods (don't look at this code - not of interest)
    if (myMethodList != "") {
        for (std::map<std::string, int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

        std::vector<TString> mlist = TMVA::gTools().SplitString(myMethodList, ',');
        for (UInt_t i = 0; i < mlist.size(); i++) {
            std::string regMethod(mlist[i]);

            if (Use.find(regMethod) == Use.end()) {
                std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
                for (std::map<std::string, int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
                std::cout << std::endl;
                return;
            }
            Use[regMethod] = 1;
        }
    }

    // --------------------------------------------------------------------------------------------------

    // --- Here the preparation phase begins

    // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
    TFile* outputFile = TFile::Open(outFileName, "RECREATE");
    TString classificationBaseName = TString("discr_") + sampleNames + TString("_");

    TMVA::Factory *factory = new TMVA::Factory(classificationBaseName.Data(), outputFile,
            "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

    (TMVA::gConfig().GetIONames()).fWeightFileDir = outputWeightsDir;

    // Define the input variables that shall be used for the MVA training

    factory->AddVariable("lep1Pt", "p_{T}(l_{1})", "GeV", 'F');
    factory->AddVariable("lep2Pt", "p_{T}(l_{2})", "GeV", 'F');

    factory->AddVariable("diLepPt", "p_{T}(ll)", "GeV", 'F');
    factory->AddVariable("diLepEta", "|#eta(ll)|", "", 'F');
    factory->AddVariable("diLepM", "m_{ll}", "GeV", 'F'); // do not use when using individual bg's that have real Z->ll
    factory->AddVariable("met", "MET", "GeV", 'F');
    factory->AddVariable("metOverQt", "MET/p_{T}(ll)", "", 'F');
    factory->AddVariable("metProjOnQt", "MET || p_{T}(ll)", "GeV", 'F');
    factory->AddVariable("metPerpQt", "MET perp  p_{T}(ll)", "GeV", 'F');
    factory->AddVariable("dPhiMetDiLep", "#Delta#phi(p_{ll},MET)", "rad", 'F');
    factory->AddVariable("mt", "m_{T}(HZZ hypothesis)", "GeV", 'F');

    factory->AddVariable("nJets15", "Number of jets with p_{T}>15", "", 'I');


//if (_jetMulti > 0) 
//       factory->AddVariable("dPhiJetMet", "#Delta#phi(Closest jet,MET)", "rad", 'F');
   if (_jetMulti > 1) 
     {
       factory->AddVariable("deltaEtaDiJet", "#Delta#eta(jet1,jet2)", "rad", 'F');
       factory->AddVariable("massDiJet", "M(jet1,jet2)", "GeV", 'F');
       factory->AddVariable("zeppDiJetDiLep", "Zeppelfeld variable", "rad", 'F');
     }

 
    //    if (_bgName.EqualTo("TTBAR", TString::kIgnoreCase))
//if (_bgName == "ttbar" || _bgName.Contains("allBg", TString::kIgnoreCase))
//     factory->AddVariable("diLepM", "m_{ll}", "GeV", 'F'); // do not use when using individual bg's that have real Z->ll

    // angle between MET and the jets
    //    if (_jetMulti > 0) //  && _bgName.EqualTo("DYMM", TString::kIgnoreCase)) 
    //        factory->AddVariable("dPhiJetMet", "#Delta#phi(jet,MET)", "rad", 'F');


    if (!USE_SEPARATE_TRAIN_TEST_FILES) {
        TFile *sigFile = TFile::Open(sigFileName);
        TFile *bgFile = TFile::Open(bgFileName);
        std::cout << "--- TrainMva       : Using input files: " << sigFile->GetName() << " and  " << bgFile->GetName() << std::endl;
        // --- Register the training and test trees
        TTree *signal = (TTree*) sigFile->Get("mvaTree");
        TTree *background = (TTree*) bgFile->Get("mvaTree");

        // You can add an arbitrary number of signal or background trees
        factory->AddSignalTree(signal, signalWeight);
        factory->AddBackgroundTree(background, backgroundWeight);
    } else {
        TFile *sigFile_train = TFile::Open(sigFileName_train);
        TFile *bgFile_train = TFile::Open(bgFileName_train);
        TFile *sigFile_test = TFile::Open(sigFileName_test);
        TFile *bgFile_test = TFile::Open(bgFileName_test);

        std::cout << "--- TrainMva       : Using input files: "
                << sigFile_train->GetName() << " and  " << bgFile_train->GetName() << "\n"
                << sigFile_test->GetName() << " and  " << bgFile_test->GetName() << std::endl;

        TTree *signal_train = (TTree*) sigFile_train->Get("mvaTree");
        TTree *background_train = (TTree*) bgFile_train->Get("mvaTree");
        TTree *signal_test = (TTree*) sigFile_test->Get("mvaTree");
        TTree *background_test = (TTree*) bgFile_test->Get("mvaTree");

        factory->AddSignalTree(signal_train, signalWeight * 2, "Training");
        factory->AddBackgroundTree(background_train, backgroundWeight * 2, "Training");
        factory->AddSignalTree(signal_test, signalWeight * 2, "Test");
        factory->AddBackgroundTree(background_test, backgroundWeight * 2, "Test");
    }

    // Weights stored in the input ntuples to be used for event weighting.
    // For mixed bg samples these will contain scaling factors that ensure the 
    // different components are mixed according to the xsections.
    // Vertex multiplicity weighting (where applicable) should also be included.
    
    factory->SetSignalWeightExpression("fullWeight");
    factory->SetBackgroundWeightExpression("fullWeight");

    if (!USE_SEPARATE_TRAIN_TEST_FILES) {
        factory->PrepareTrainingAndTestTree(addCutsSig, addCutsBg,
                "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V"); // I used this one when reading train/test from same file
    } else {
        // I prefer using separate files for training/testing, so this is just for completeness.
        
        // note that here the signal and the bg cuts are the same, so we can use this form
        factory->PrepareTrainingAndTestTree(addCutsSig,
                _numSignalTrain, _numBgTrain, _numSignalTest, _numBgTest,
                "SplitMode=Random:NormMode=NumEvents:!V");
    }


    std::cout << "==> Booking the methods " << std::endl;
    // ---- Book MVA methods
    //
    // Please lookup the various method configuration options in the corresponding cxx files, eg:
    // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
    // it is possible to preset ranges in the option string in which the cut optimisation should be done:
    // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

    // Cut optimisation
    if (Use["Cuts"])
        factory->BookMethod(TMVA::Types::kCuts, "Cuts",
            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart");

    if (Use["CutsD"])
        factory->BookMethod(TMVA::Types::kCuts, "CutsD",
            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate");

    if (Use["CutsPCA"])
        factory->BookMethod(TMVA::Types::kCuts, "CutsPCA",
            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA");

    if (Use["CutsGA"])
        factory->BookMethod(TMVA::Types::kCuts, "CutsGA",
            "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95");

    if (Use["CutsSA"])
        factory->BookMethod(TMVA::Types::kCuts, "CutsSA",
            "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale");

    // Likelihood ("naive Bayes estimator")
    if (Use["Likelihood"])
        factory->BookMethod(TMVA::Types::kLikelihood, "Likelihood",
            "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50");

    // Decorrelated likelihood
    if (Use["LikelihoodD"])
        factory->BookMethod(TMVA::Types::kLikelihood, "LikelihoodD",
            "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate");

    // PCA-transformed likelihood
    if (Use["LikelihoodPCA"])
        factory->BookMethod(TMVA::Types::kLikelihood, "LikelihoodPCA",
            "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA");

    // Use a kernel density estimator to approximate the PDFs
    if (Use["LikelihoodKDE"])
        factory->BookMethod(TMVA::Types::kLikelihood, "LikelihoodKDE",
            "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50");

    // Use a variable-dependent mix of splines and kernel density estimator
    if (Use["LikelihoodMIX"])
        factory->BookMethod(TMVA::Types::kLikelihood, "LikelihoodMIX",
            "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50");

    // Test the multi-dimensional probability density estimator
    // here are the options strings for the MinMax and RMS methods, respectively:
    //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
    //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
    if (Use["PDERS"])
        factory->BookMethod(TMVA::Types::kPDERS, "PDERS",
            "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600");

    if (Use["PDERSD"])
        factory->BookMethod(TMVA::Types::kPDERS, "PDERSD",
            "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate");

    if (Use["PDERSPCA"])
        factory->BookMethod(TMVA::Types::kPDERS, "PDERSPCA",
            "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA");

    // Multi-dimensional likelihood estimator using self-adapting phase-space binning
    if (Use["PDEFoam"])
        factory->BookMethod(TMVA::Types::kPDEFoam, "PDEFoam",
            "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T");

    if (Use["PDEFoamBoost"])
        factory->BookMethod(TMVA::Types::kPDEFoam, "PDEFoamBoost",
            "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T");

    // K-Nearest Neighbour classifier (KNN)
    if (Use["KNN"])
        factory->BookMethod(TMVA::Types::kKNN, "KNN",
            "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim");

    // H-Matrix (chi2-squared) method
    if (Use["HMatrix"])
        factory->BookMethod(TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None");

    // Linear discriminant (same as Fisher discriminant)
    if (Use["LD"])
        factory->BookMethod(TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");

    // Fisher discriminant (same as LD)
    if (Use["Fisher"])
        factory->BookMethod(TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");

    // Fisher with Gauss-transformed input variables
    if (Use["FisherG"])
        factory->BookMethod(TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss");

    // Composite classifier: ensemble (tree) of boosted Fisher classifiers
    if (Use["BoostedFisher"])
        factory->BookMethod(TMVA::Types::kFisher, "BoostedFisher",
            "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring");

    // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
    if (Use["FDA_MC"])
        factory->BookMethod(TMVA::Types::kFDA, "FDA_MC",
            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1");

    if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
        factory->BookMethod(TMVA::Types::kFDA, "FDA_GA",
            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1");

    if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
        factory->BookMethod(TMVA::Types::kFDA, "FDA_SA",
            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale");

    if (Use["FDA_MT"])
        factory->BookMethod(TMVA::Types::kFDA, "FDA_MT",
            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch");

    if (Use["FDA_GAMT"])
        factory->BookMethod(TMVA::Types::kFDA, "FDA_GAMT",
            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim");

    if (Use["FDA_MCMT"])
        factory->BookMethod(TMVA::Types::kFDA, "FDA_MCMT",
            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20");

    // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
    if (Use["MLP"])
      factory->BookMethod(TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=tanh:VarTransform=N:NCycles=500:HiddenLayers=N:TestRate=10");

    if (Use["MLPBFGS"])
        factory->BookMethod(TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator");

    if (Use["MLPBNN"])
        factory->BookMethod(TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator"); // BFGS training with bayesian regulators
    // reduced epochs from 600 (AA)


    // CF(Clermont-Ferrand)ANN
    if (Use["CFMlpANN"])
        factory->BookMethod(TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"); // n_cycles:#nodes:#nodes:...  

    // Tmlp(Root)ANN
    if (Use["TMlpANN"])
        factory->BookMethod(TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"); // n_cycles:#nodes:#nodes:...

    // Support Vector Machine
    if (Use["SVM"])
        factory->BookMethod(TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm");

    // Boosted Decision Trees - DEFAULT
    //    if (Use["BDTG"]) // Gradient Boost -> Default parameters
    //        factory->BookMethod(TMVA::Types::kBDT, "BDTG",
    //            "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5");


    // here is a version with modified parameters
    if (Use["BDTG"]) // Gradient Boost 
        factory->BookMethod(TMVA::Types::kBDT, "BDTG",
            "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:NNodesMax=8:nEventsMin=12");


    if (Use["BDT"]) // Adaptive Boost
        factory->BookMethod(TMVA::Types::kBDT, "BDT",
            "!H:!V:NTrees=400:nEventsMin=40:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=1:SeparationType=GiniIndex:nCuts=20");


    if (Use["BDTB"]) // Bagging
        factory->BookMethod(TMVA::Types::kBDT, "BDTB",
            "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");

    if (Use["BDTD"]) // Decorrelation + Adaptive Boost
        factory->BookMethod(TMVA::Types::kBDT, "BDTD",
            "!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate");

    if (Use["BDTF"]) // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
        factory->BookMethod(TMVA::Types::kBDT, "BDTMitFisher",
            "!H:!V:NTrees=50:nEventsMin=150:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");

    // RuleFit -- TMVA implementation of Friedman's method
    if (Use["RuleFit"])
        factory->BookMethod(TMVA::Types::kRuleFit, "RuleFit",
            "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02");

    // For an example of the category classifier usage, see: TMVAClassificationCategory

    // --------------------------------------------------------------------------------------------------

    // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

    // factory->OptimizeAllMethods("SigEffAt001","Scan");
    // factory->OptimizeAllMethods("ROCIntegral","GA");

    // --------------------------------------------------------------------------------------------------

    // ---- Now you can tell the factory to train, test, and evaluate the MVAs

    // Train MVAs using the set of training events
    factory->TrainAllMethods();

    // ---- Evaluate all MVAs using the set of test events
    factory->TestAllMethods();

    // ----- Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();

    // --------------------------------------------------------------

    // Save the output
    outputFile->Close();

    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TrainMva is done!" << std::endl;

    delete factory;

    // Launch the GUI for the root macro
    // make it lightweight for batch jobs and skip loading this script -> for interactive include
    // TMVAGui.C which is currently commented out.
    // if (!gROOT->IsBatch()) TMVAGui(outFileName);

}

