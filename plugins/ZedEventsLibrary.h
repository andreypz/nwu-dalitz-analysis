/*
    Utilities for making a library of Z events and retrieving them.
*/

#ifndef _ZedEventsLibrary_H
#define _ZedEventsLibrary_H

#include <string>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "TROOT.h"
#include "TTree.h"
#include "TObject.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TLorentzVector.h"

using namespace std;

class ZedEventsLibrary: public TObject {
    public:
        ZedEventsLibrary();
        virtual ~ZedEventsLibrary();
        ZedEventsLibrary(string selection, bool makeFromData);

        void  CreateLibrary();        
	pair <int, int> GetEtaQtBin(Float_t eta, Float_t qt);
        pair<TLorentzVector, TLorentzVector> GetLeptonsFromLibrary(TLorentzVector origZP4);

        ClassDef(ZedEventsLibrary, 0);

    private:
        //input parameters
        string _selection;
        bool   _makeFromData;

	int    _qt_bins[10];
	float  _eta_bins[7];

        TFile * _inFile;

	TRandom3 *_myRandom;

	//TLorentzVector *_lep1, *_lep2;
	std::vector< pair<TLorentzVector, TLorentzVector> >  _diMap[6][10];

};

#endif





#if !defined(__CINT__)
ClassImp(ZedEventsLibrary);
#endif
