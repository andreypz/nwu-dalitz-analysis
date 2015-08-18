My analysis at CMS
------------------

Here is my analyses that run on the nuTuples produced with the dedicated [ntupleProducer][1], by NWU group.


Multiple analyses are stored inside:
 * **higgs &rarr; gamma* gamma &rarr; ll gamma** as well as other channels with the same final state. E.g. h &rarr; J/Psi + gamma. etc
 * **higgs &rarr; 4 leptons**, a basic analyzer available
 * **higgs &rarr; ZZ &rarr; 2l2nu**, an old analysis kept for history

### zgamma
Inside this directory there are multiple analyzers looking into a similar final states: l+l-gamma.
Those are:
  * *zgamma.C* for muon channel,
  * *egamma.C* for electron channel,
  * *egamma2.C* for electron channel in category 2 (2 separate electrons),
  * *jpsiGamma.C* - looking into J/Psi+gamma final state, where J/Psi goes into two muons,
  * *apz.py* - searching for a new particle at m_a = 18 GeV or so

In order to run all of those I use *run.py* script. Here is how it goes for the local job to run the jpsiGamma analyzer:
```run.py myTag file-list.txt -s jp-mugamma -t mugamma```,
where *myTag* is appended to the output file name. It is  also passed to the analyzer itself, so that one could
distinguish samples and apply different selections etc. *file-list.txt* is  a text file thst contains a list of input root files
(full path names).

### testAnalyzer
  * *natenalyzer.C* - same sign selection in di-muon channel
  * *fourLeptons.C* - a simple 4l analyzer also searching for an apz particle

### fits-and-limits
Scripts to set the limits on *h&rarr;llgamma* channel (may be extended to other situations). See the readme file inside it.

### batch-condor
BatchMaster scripts to submit jobs to condor. Works both at NU-tier3 and FNAL-lpc


### h2l2nu
This is  more or less complete analysis on 2l2nu selection.

### mva
MVA tools for h2l2nu analysis. Old, but maybe useful for new analyses.

### interface and src
That is where I keep all the TC object classes and definitions

### plugins
A bunch of very useful plugins are stored here:
  * *HistManager.C* - the most significant contribution to HEP by AK-19. This plugin makes the histograms in-one-line code.
  * *ZGAngles.C* - this does the angles calculation in the llgamma final states
  * See its readme file for more details

### misc
In order to set up root I use *cmsenv.sh* script:
```source cmsenv.sh```

* *zgamma/LHEanalyzer.py* - does the analysis on LHE files converted into root trees using ExRootAnalysis tool

[1]: https://github.com/NWUHEP/ntupleProducer
