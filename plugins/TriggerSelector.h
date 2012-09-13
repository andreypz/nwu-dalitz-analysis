/*
    Utilities for retrieving weights for PU,etc.
*/

#ifndef _TriggerSelector_H
#define _TriggerSelector_H

#include <string>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "TROOT.h"
#include "TObject.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

using namespace std;

class TriggerSelector: public TObject {
    public:
        TriggerSelector();
        virtual ~TriggerSelector();
        TriggerSelector(string selection, string dataPeriod, vector<int> triggers);
        bool SelectTriggers(unsigned int triggerStatus, UInt_t hltPrescale[]);
        bool PhotonTriggerBins(float photonPt, bool isoTriggers) const;
        int GetEventPrescale() const;
        int GetPassTrigger() const;

        ClassDef(TriggerSelector, 0);

    private:
        //input parameters
        string _dataPeriod;
        string _selection;
        vector<int> _triggers;

        //trigger info
        bool   _eventPass;
        int    _eventPrescale;
        int    _passTrigger;
        vector<int> _passTriggers;
};

#endif

#if !defined(__CINT__)
ClassImp(TriggerSelector);
#endif
