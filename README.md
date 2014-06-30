My analysis at CMS
------------------

Here is my analyses that run on the nuTuples produced with the dedicated [ntupleProducer][1], by NWU group.


Multiple analyses are storred inside:
 * **higgs &rarr; gamma* gamma &rarr; ll gammma** as well as other chanels with the same final state. E.g. h -> J/Psi + gamma. etc
 * **higgs -> 4 leptons**, some basic analyzer available
 * **higgs -> ZZ -> 2l2nu**, an old analysis kept for history

### zgamma
Inside this directory there are multiple analyzerslook into a similar final states: l+l-gamma.
Those are:
  * *zgamma.C* for muon channel,
  * *egamma.C* for electron channel,
  * *jpsiGamma.C* looking into J/Psi+gamma final state, where J/Psi goes into two muons.
  * *apz.py* - searching for a new particle at m_a = 18 GeV or so
  * *fourLeptons.C* - a simple 4l analyzer also searching for an apz particle
In order to run all of those I use *run.py* script. Here is how it goes for the local job to run jpsiGamma analyzer:
```run.py myTag file -s jp-mugamma -t mugamma```



### fits-and-limits
Scripts to set the limits on h->llgamma channel (may be extended to other situations). See the readme file isede it.

### batch-condor
BatchMaster scripts to submit jons. Works both at NU-tier3 and FNAL-lpc


### h2l2nu
This is  more or less complete analysis on 2l2nu selection.

### mva
MVA tools for h2l2nu analysis. Old, but maybe useful for a new analyses.

### interface ans src
That is where I keep all the TC object classes and definitions
### plugind
A bunch of very useful plugins are stored here:
  * *HistManager.C* - the best contribution to HEP by A.K.19, makes the histograms in one line code.
  * *ZGAngles.C* - does the angles calculation in the llgamma final states
  * See its readme file for more details

### misc
In order to set up root I use *cmsenv.sh* script:
```source cmsenv.sh```

* *zgamma/LHEanalyzer.py* - does the analysis on LHE files converted into root trees using ExRootAnalysis tool

[1]: https://github.com/NWUHEP/ntupleProducer
