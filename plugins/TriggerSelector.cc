#include "TriggerSelector.h"

TriggerSelector::TriggerSelector() {}

TriggerSelector::~TriggerSelector() {}

TriggerSelector::TriggerSelector(string selection, string dataPeriod, vector<int> triggers)
{
    _selection  = selection;
    _dataPeriod = dataPeriod;
    _triggers   = triggers;
}

bool TriggerSelector::SelectTriggers(unsigned int triggerStatus, int hltPrescale[])
{
    _eventPrescale = -1;
    _passTrigger   = 0;
    _eventPass     = false; 

    if (_triggers[0] != 0) { 
        int prescale = 1e9;
        for (int i = 0; i < _triggers.size(); ++i) {
            unsigned int iHLT = 0x01 << (_triggers[i] - 1);  

            // Check that analysis triggers are unprescaled, save lowest prescale for photon triggers
            if (_selection == "muEG" || _selection == "muon" || _selection == "electron") 
                if (hltPrescale[_triggers[i]-1] != 1) continue;

            if ((triggerStatus & iHLT) == iHLT) {
                if (prescale > hltPrescale[_triggers[i]-1]) { 
                    _passTrigger = _triggers[i];
                    prescale = hltPrescale[_triggers[i]-1];
                }
                _passTriggers.push_back(_triggers[i]);
                _eventPass = true; 
                if (prescale == 1) break;
            }
        }
        _eventPrescale = prescale;
    } else {_eventPass = true;}

    return _eventPass;
}

bool TriggerSelector::PhotonTriggerBins(float photonPt, bool isoTrigger) const
{
    bool photonPass = false;
    //float ptBins[] = { 25., 35., 45., 55., 80., 95., 140.};

    if ( _dataPeriod == "2011B"
            && _passTrigger == _triggers[_triggers.size()-1]
            && photonPt > 140.
       ) photonPass = true;
    else if ((_passTrigger == _triggers[8] || (!isoTrigger && _passTrigger == _triggers[9]))
            && (photonPt > 95. && (photonPt < 140. && _dataPeriod == "2011B"))
            ) photonPass = true;
    else if ((_passTrigger == _triggers[6] || (!isoTrigger && _passTrigger == _triggers[7]))
            && (photonPt > 80. && photonPt < 95.)
            ) photonPass = true;
    else if ((_passTrigger == _triggers[4] || (!isoTrigger && _passTrigger == _triggers[5]))
            && (photonPt > 55. && photonPt < 80.)
            ) photonPass = true;
    else if ((_passTrigger == _triggers[2] || (!isoTrigger && _passTrigger == _triggers[3]))
            && (photonPt > 35. && photonPt < 55.)
            ) photonPass = true;
    else if (_passTrigger == _triggers[0]
            && (photonPt > 25. && photonPt < 35.)
            ) photonPass = true;

    return photonPass;
}

int TriggerSelector::GetEventPrescale() const
{
    return _eventPrescale;
}

int TriggerSelector::GetPassTrigger() const
{
    return _passTrigger;
}
