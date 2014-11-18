## Fits,  Limit calculation and  Bias study
These codes are meant for Higgs to l+l-gamma Dalitz mode analysis,  but it could be used for any other similar tasks: fit the final mass distribution and search for peaks.
Most of it was orriginally taken from Brover Clevelend [BiasAndLimits code][brover].
Debugged and improved.

### How to run it. 
#### Fits and limits
 0. ```config.cfg``` sepcifies various parameters to be used: signal mass points and production channels, paths to output directories, colors of the fits, etc.
 1. ```apzFitProducer.py  path/to/rootfile``` to perform initial fits. Where _path/to/rootfile_ is specific to my setup (look into te code), and the root files in there have to contain a tree called 'fitTree/fitTree'. That tree has to contain the following variables: **m_llg**, **m_ll**, **di_pt**, **ph_pt**, **ph_eta** and **weight**. The **m_llg** is the most important one: three-body invariant mass that we are going to fit. 
 It performs various different fits to Data (Exponential, Bernsteins etc) and also fits all signal MC to a  Gaussian plus a Crystall Ball function with the same mean.
 2. ```bgCardPrep.py``` makes a final fit of the Data to a choosen background model (e.g. Bern5 polynomials). It also renames the variables and created the RooFit Workspace for the use by Combination tool. 
 3. ```signalCBFits.py```  takes the fits to the signal perforemed at step 1 and performs an interpolation to other masses in beteen. Renames the variables and creates the RooFit workspaces to be used by Combination tool. 
 4. ```cardMaker.py``` makes the datacards.
 5. ```limitProducer.py``` produces the limit by runnning the Combination tool. *The proper CMSSW version, recommended by Combination guys, has to be used at this stage!*. 
 6. Then ``limitPlotter.py``` makes the plots for the limits.

#### Bias study
 1. ```biasStudy_toyMaker.py``` makes the toys from a given PDF and performes the fits with various test function.
 2. ```plotPulls.py``` makes the pull plots from all the fits
 
[brover]: https://github.com/brovercleveland/BiasAndLimits
