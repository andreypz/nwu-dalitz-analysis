HistManager
-----------
It's for filling/making all the histograms in one line.


TriggerSelector
------------------
Selects HLT trigger path, stored in ntuples.

ZedEventsLibrary
-------------------
A tool to create a library of di-lepton pairs from data in order to use for Z+jets background estimate in h2l2nu analysis.


WeightUtils
-------------
Weigts for pile-up and lepton efficiency. 

rochcor.cc
----------
Rochester muon corrections. It's old need to be updated for future use.


pwhg_cpHTO_reweight.f
--------------------
Code for creating higgs lineShape library.
The library itself is located in /lib

In order to use the weighting in root, one needs to create a shared library:
'''
gfortran -shared -O2 *.f -o libShapeLine.so -fPIC
'''

Then, load this library to root by gSystem->Load("libShapeLine.so");
''' C 
Define an external function in .h file: 
 extern"C" {
   void pwhg_cphto_reweight_(double *mh, double *gh, double *mt, int *BWflag, double *m, double *w);
           }
'''
