The log of changes:
==================

2013-Nov-13
  * Muon soft + (tight vertex) ID is synchronized with Stoyan (almost)
  * Low-pt muon category introduces
  * A tree for fitting input is optimized
  * An option to use MCFM or MAdgraph signal sample
  
  * Also, improved fits-and-limits: 
    * loop over tree python way
    * introduced 'Low' category
    * GaussBern3 used instead of GaussBern6

2013-Oct-30
  * Dalitz plots added
  * EB/EE categories introduced
  * Plotting scripts improved: significance and yields added

2013-Oct-18
  * Mva Tree is saved for further electron training
  * Options in TSelector class are used (instead of template files)

2013-Oct-17
  * cmsenv.sh script is added: no need to install CMSSW anymore
  * First selection of merged dalitz electrons is implemented
  * Photon ID and iso synchronized with official recommendations
  
2013-Oct-12
  * Running on v8 ntuples, tag v4.3
  * Interface and linkdefs implemented just like in ntuple producer
  * Analysis on DY and dalitz electrons

2013-Oct-03
  * DiscoverGeneology()
  * Figured out mothers for leptons and photons
  * Implementing synch cuts with Stoyan (triangle, muon ID, etc)

2013-Oct-02
  * Updated to run over v8 ntuples:
    * Electron, Photon, Met classes changed (Mva Id, regression etc,)
    * More Muon variables added (soft muons)
    
2013-Sep-11
  * Added ZGAngles plugin
    * Need to be validated against Kristian's code (from Brian)
    * Different approached don't give the same results in angles, need to understand why.

2013-Sep-09
  * Batch master improvements:
    * Output root file is now copied into Stage directory along with reports
  * Histograms for acceptance vs Mll added
   
2013-Sep-05
  * Improvements to zgamma: 
    * more plots
    * selection for muons updated according to POG (new TCMuon object added)
     
2013-Aug-28
  * Adding LHEanalyzer (in zgamma) form making plots of LHE root files produced with ExRootLHEFConverter. Need a shared library to run it.
  * Switched to Stoyan's cut ordering (temporarily)  	
 
2013-Aug-14
  * Adding a config file as well as the parser

2013-Aug-13
  * Make Id/Iso plots for muons.
  * Photon eta < 2.5 cut (apparently the go up to 3.0 by default. Why?)
  * Run with different trigger selections implemented: 
        double-mu, single-mu, muon-photon, double-photon (for electron selection).
  * Data samples added to run on batch.
  * Plotting script updated: a mix of h2l2nu and a simple plotting from a file.

2013-Aug-06
  * zgamma analyzer developments:
    * Trigger efficiency study for h->llgamma dalitz decay
    * Selections for mu, el, and mugamma introduced
    * Filling histograms in a similar way as in h2l2nu analysis (by cut)
  * batch master revived and made to work with zgamma sample.

2013-Jul-04
  * Ported to Git!
  * Updated the TC classes with the latest version
  * Introduced Zgamma analyzer

2013-May-10
  * Apparently, forgot to commit a TCPhysObject class, doing it now.
  * Tag with V02-16    

2013-Apr-18
  * Tagged with V02-15
    * See previous logs for changes
    
2013-Apr-17
  * Nate's tri-lepton plot is produced
  * Trigger Selector is modified to return: isFound, isPassed and a prescale
  * Modified the Batch master to create a tar.gz before submission and to write output to EOS
  
2012-Dec-06
  * Introduced option parser to run.py script 
  * Synchronization with Nate is more or less successful
    * Ulong64 bug in trigger selector is fixed 
    * Conversion veto for electrons

2012-Nov-13
  * V02-14 tag
    * Pile-up for 2012 samples is included
    * More plots in Muons and electrons
    * Cut on dR(ele, mu) to reject fake electrons
    * Top background estimation scripts are added

2012-Oct-31
  * Made it to run with 2012 and 2011 ntuples using same code (although, the ntuples are different)
    * Lepton selection is modified
  * b quarks study and top-background estimation scripts are updated

2012-Oct-07
   * 8 TeV cross section in config file
    
2012-Oct-06
  * V02-13 tag for the latest version of 7TeV analyzer
  * Moving to struct for Mu/Ele isolation and Id.

2012-Sep-29
  * V02-12 tag
    * HistManager added
    * plugins/lib for lineshape library
    * synchronization for 8 TeV
    
2012-Sep-23
  * Adding HistManager class from Andy. Implemented it in higgs analyzer
  * Using N events from Notify to scale the samples to lnmi
  * Removing soucefiles - no use anymore

2012-Sep-21
 * Adding event counts in Notify()
 * Cleaning sourcefiles
 * PU reweighting with fine binning for 2011

2012-Sep-13
 * V02-09 tag
   * k-factors for mH 125 and 200
   * commit before moving to a new ntuple format
   
2012-Sep-06
 * V02-08 tag:
   * MVA for low masses is done (weights are in data). 
   * Added variables for mva trees. 
   * Plotting script improved and more...

2012-Aug-23
 * Recent tags: 
   * V02-07 - Updated MVA cuts, tree maker and training  
   * V02-05 - Workable MVA code 
   * V02-03 - Added rochcor function for correcting muons, nJets histogram for jets in eta<2.5; 
         batch maser is modified to copy the files to a local directory before sending them over to batch 

2012-Jul-18
  * V02-02 tag
  * Converted scripts and plotting macros to python
  * Updates in analyzer in order to synchronize with Common analysis:
    * b-tagging
    * jet cuts
    * lepton iso and third lepton veto
    * Met and MT cuts re-optimized
    * etc
    
2012-Jun-14
  * V02-01 tag
  * Lineshape reweighting plugins
  * Updating to the recent TC object classes

2011-Dec-17
  * New version of Weight Utils and re-weight files.
      * Fixed k-nlo for ggH samples 
  * Third muon veto updated (closer to Nate/ Pas-like)
  
2011-Dec-14
  * simplified drawing macro (no need to make all stacks - they are produced in drawMulti function)
  
2011-Dec-13
  V01-04 tag
  * Z library class is working

2011-Dec-09
  V01-03 tag
  * Switching to B-master for condor jobs submissions
  * Created plugins directory where all new classes will go
  (moved there weight utils and trigger selector)
  * Starting ZedEvensLibrary class

2011-Dec-05
  V01-02 tag
  * gluglu Higgs reweighting is implemented, kinTree ntuples are reproduced for higgs
  
2011-Dec-04
  V01-01 tag
  * TriggerSelector added (Nate's)
  * Z+jets Library update: more eta-bins
  * Overflow and underflow bins are now added to the last/first bins of all hists
  
2011-Dec-03
  * Weight utils updated: GetPhotonMass moved to Weight Utils
  * Changed the colors for the histograms   

2011-Dec-01 
  * Weighting utils class is implemented (from Nate) 
  * dPhi is changed for 0-jet bin (adding soft jets)
  * WW,WZ,ZZ samples switching to Madgraph (xSec and source files updated)
  * Trigger prescale bug fixed
  
2011-Nov-24
  The following updates are added to the code:
  * Some more histograms added: angle correlations of leptons, dd-lepton, met
  * Z-library implementation
  * Photon+gets reweighting for Z+jets background estimation
  * Ntuples for MVA, in "Anton" directory of a root file
  * The H-mass dependent cuts are changed to new PAS-like cuts
  
2011-Oct-17
 *  VBF signal included

2011-Oct-16
  * V00-08 - for PAS-like yields
  * PAS-like selection restored

2011-Sep-21
*  Created MetDefinitions.h where to put all the Met functions (to share with others)
*  Put a several macros from makeplot to utils - this will be a general place to put such things
*  Deleted merge.C (included in utils.C)
  
2011-Sep-16
*  V00-07 	Fixes, matching Nate's yields; third lepton electron updated
*  V00-06 	typos, broken 

2011-Aug-26
*  tagged with V00-06
*  Moved to higgs7 directory
*  Updated cs and Nev (few things were missed)
*  Mote ttbar samples to run
*  Renamed a bunch of files

  
2011-Aug-25
*  V00-05 tagged
*  Removed old analyzer and scripts
*  fixed trigger selection; updated the plotting macros


2011-Aug-23
*  Removed Ntuplizer code
*  Added run.csh that was missing

V00-04 - synchronized with Nate's code

V00-03 - latest version based on old code.
