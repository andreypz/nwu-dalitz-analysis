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
                if (hltPrescale[_triggers[i]-1] == 1) {
                    _passTrigger = _triggers[i];
                    prescale = 1;
                } else if (prescale > hltPrescale[_triggers[i]-1]) { 
                    _passTrigger = _triggers[i];
                    prescale = hltPrescale[_triggers[i]-1];
                }
                _passTriggers.push_back(_triggers[i]);
                _eventPass = true; 
            }
        }
        _eventPrescale = prescale;
    } else {_eventPass = true;}

    return _eventPass;
}

int TriggerSelector::GetEventPrescale() const
{
    return _eventPrescale;
}

int TriggerSelector::GetPassTrigger() const
{
    return _passTrigger;
}
