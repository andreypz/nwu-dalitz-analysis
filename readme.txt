Here is Higgs to ZZ to 2l2nu Analysis.
More details on this code with the instructions how to run it can be found here:
https://twiki.cern.ch/twiki/bin/view/CMS/UserCodeNWUHiggs


The log of changes:

2012-Jul-18
  * V02-02 tag
  * Converted scripts and plotting macros to python
  * Updates in analyzer in order to synchronize with Common analysis:
    * b-tagging
    * jet cuts
    * lepton iso and third lepton veto
    * Met and MT cuts reoptimised
    * etc
    
2012-Jun-14
  * V02-01 tag
  * Lineshape reweighting plugins
  * Updating to the recent TC object classes

2011-Dec-17
  * New version of Weight Utils and reweight files.
      * Fixed k-nlo for ggH samples 
  * Third muon veto updated (closer to Nate/ Pas-like)
  
2011-Dec-14
  * simplified drawing macro (no need to make all stacks - they are produced in drawMulti function)
  
2011-Dec-13
  V01-04 tag
  * Z library class is working

2011-Dec-09
  V01-03 tag
  * Switching to B-master for condor jobs submitions
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
  * Some more histograms added: angle correlationsleptons of leptons, di-lepton, met
  * Z-library implementation
  * Photon+gets reweighting for Z+jets background estimation
  * Ntuples for MVA, in "Anton" directory of a root file
  * The H-mass dependent cuts are changed to new PAS-like cuts
  
2011-Oct-17
  VBF signal included

2011-Oct-16
  V00-08 - for PAS-like yields
  PAS-like selection restored

2011-Sep-21
  Created MetDefinitions.h where to put all the Met functions (to share with others)
  Put a several macros from makeplot to utils - this will be a general place to put such things
  Deleted merge.C (included in utils.C)
  
2011-Sep-16
  V00-07 	Fixes, matching Nate's yields; third lepton electron updated
  V00-06 	typos, broken 

2011-Aug-26
  tagged with V00-06
  Moved to higgs7 directory
  Updated cs and Nev (few things were missed)
  Mote ttbar samples to run
  Renamed a banch of files

  
2011-Aug-25
  V00-05 tagged
  Removed old analyzer and scripts
  fixed trigger selection; updated the plotting macros


2011-Aug-23
  Removed Ntuplizer code
  Added run.csh that was missing

V00-04 - synchronized with Nate's code

V00-03 - latest version based on old code.

