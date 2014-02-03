## Fits,  Limit calculation and  Bias study
These codes are meant for Higgs to l+l-gamma Dalitz mode analyzis.
Most of it was orriginally taken from Brover Clevelend [BiasAndLimits code][brover].
Debugged and improved.

### How to run it. 
#### Fits and limits
 1. ```initialFitProducer.py``` to perform initial fits.
 2. ```bgCardPrep.py``` and ```signalCBFits.py``` for the final fits and the info for the datacards
 3. ```cardMaker.py``` makes the datacards
 4. ```limitProducer.py``` produces the limit. Then ``limitPlotter.py``` makes the plots out of the result.

#### Bias study
 1. ```biasStudy_toyMaker.py``` makes the toys from a given PDF and performes the fits with various test function.
 2. ```plotPulls.py``` make the pull plots framm all the toys
 
[brover]: https://github.com/brovercleveland/BiasAndLimits
