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
#include "TObject.h"
#include "TFile.h"

using namespace std;

class ZedEventsLibrary: public TObject {
    public:
        ZedEventsLibrary();
        virtual ~ZedEventsLibrary();
        ZedEventsLibrary(string selection, bool makeFromData);

        bool GetBin() const;

        ClassDef(ZedEventsLibrary, 0);

    private:
        //input parameters
        string _selection;
        bool   _makeFromData;
};

#endif

#if !defined(__CINT__)
ClassImp(ZedEventsLibrary);
#endif
