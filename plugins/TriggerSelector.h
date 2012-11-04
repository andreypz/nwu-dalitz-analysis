/*
    Utilities for retrieving weights for PU,etc.
*/

#ifndef _TriggerSelector_H
#define _TriggerSelector_H

// c++ libraries
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

// ROOT libraries
#include "TROOT.h"
#include "TObject.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

using namespace std;

typedef vector<string> vstring;

class TriggerSelector: public TObject {
    public:
        TriggerSelector();
        virtual ~TriggerSelector();
        TriggerSelector(string, string, vstring);

        bool    SelectTrigger(string,unsigned, UInt_t*);
        bool    CheckOverlap();
        bool    CheckPrescales(unsigned, UInt_t*);

        int     GetEventPrescale() const;

        void    SetPassNames(unsigned);
        void    SetDataBit(bool);
        void    SetSelectedBits();
        void    TriggerDefaults();
        void    AddTriggers(vstring);

        ClassDef(TriggerSelector, 0);

    private:
        // input parameters
        vstring         _triggerNames;
        vstring         _triggers;
        string          _dataPeriod;
        string          _type;
        bool            _isRealData;

        // trigger info
        vstring         _passNames;
        bool            _eventPass;
        int             _eventPrescale;
        int             _passTriggers;
};

#endif

#if !defined(__CINT__)
ClassImp(TriggerSelector);
#endif
