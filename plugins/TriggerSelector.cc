#include "TriggerSelector.h"

TriggerSelector::TriggerSelector() {}

TriggerSelector::~TriggerSelector() {}

TriggerSelector::TriggerSelector(string type, string dataPeriod, vstring triggerNames)
{
    _type           = type;
    _dataPeriod     = dataPeriod;
    _triggerNames   = triggerNames;
    _isRealData     = false;

    TriggerDefaults();
    SetSelectedBits();
}

int TriggerSelector::GetEventPrescale() const
{
    return _eventPrescale;
}

void TriggerSelector::SetDataBit(bool b) 
{
    _isRealData = b;
}

void TriggerSelector::AddTriggers(vstring trigs)
{
    _triggers.insert(_triggers.begin(), trigs.begin(), trigs.end());
}

void TriggerSelector::SetPassNames(unsigned trigs)
{
    unsigned count = 0;
    vstring test;

    while (count < sizeof(trigs)*8) {
        if (trigs & (0x1 << count)) {
            string passName = _triggerNames[count];
            _passNames.push_back(passName);
        }
        ++count;
    }
}

bool TriggerSelector::CheckPrescales(unsigned trigs, UInt_t* hltPrescales)
{
    unsigned count = 0;
    bool isPrescaled = true;

    while (count < sizeof(trigs)*8) {
        if (trigs & (0x1 << count)) {
            if (hltPrescales[count] == 1) {
                isPrescaled = false;
                break;
            }
        }
        ++count;
    }
    return isPrescaled;
}

void TriggerSelector::TriggerDefaults()
{
    /*
       Sets triggers for analysis that mixes dilepton
       stream data.  During selection, priority is given
       to sf streams with dimuons taking highest precedent. */

    // double muon triggers
  //_triggers.push_back("HLT_Mu13_Mu8_v");
    /*
    _triggers.push_back("HLT_Mu17_Mu8_v");
    _triggers.push_back("HLT_Mu17_TkMu8_v");
    _triggers.push_back("HLT_Mu22_TkMu8_v");    
    _triggers.push_back("HLT_Mu22_TkMu22_v");
    _triggers.push_back("HLT_DoubleMu7_v");
    _triggers.push_back("HLT_DoubleMu6_v");

    // double electron triggers
    _triggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v"); 
    _triggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");

    // muEG triggers
    _triggers.push_back("HLT_Mu17_Ele8_CaloIdL_v");
    _triggers.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v");
    _triggers.push_back("HLT_Mu8_Ele17_CaloIdL_v");
    _triggers.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v");

    */
}

void TriggerSelector::SetSelectedBits()
{
  /* Convert list of trigger names to a corresponding
     bit array by comparing to stored trigger names. */
  
  _passTriggers = 0x0;
  
  for (vstring::const_iterator iTrig = _triggers.begin(); iTrig != _triggers.end(); ++iTrig) {
    unsigned count = 0;
    for (vstring::const_iterator iName = _triggerNames.begin(); iName != _triggerNames.end(); ++iName) {
      if (*iTrig == *iName) {
        //cout<<"TRG dbg: iName = "<<*iName<<"   matches to iTrig = "<<*iTrig<<endl;
        _passTriggers |= 0x1 << count;
        break;
      }
      ++count;
    }
  }
  //cout<<"pass triggers  "<<_passTriggers<<endl;
}

bool TriggerSelector::CheckOverlap() {

    /* Depending on data stream, this method will
       return true if the event should be vetoed, i.e.,
       if a double muon trigger fires for an event from 
       the MuEG stream. */

    bool overlapVeto = false;

    if ( 
            binary_search(_passNames.begin(), _passNames.end(), "HLT_Mu13_Mu8_v")
            || binary_search(_passNames.begin(), _passNames.end(), "HLT_Mu17_Mu8_v")
            || binary_search(_passNames.begin(), _passNames.end(), "HLT_DoubleMu3_v")
            || binary_search(_passNames.begin(), _passNames.end(), "HLT_DoubleMu6_v")
            || binary_search(_passNames.begin(), _passNames.end(), "HLT_DoubleMu7_v")
       ) 
        if (_type == "muEG" || _type == "electron") overlapVeto = true;

    if (
            (binary_search(_passNames.begin(), _passNames.end(), "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v")
             || binary_search(_passNames.begin(), _passNames.end(), "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v"))
       ) 
        if (_type == "muEG") overlapVeto = true;

    return overlapVeto;
}

bool TriggerSelector::SelectTrigger(string name, unsigned triggerStatus, UInt_t* hltPrescales)
{

  bool isFound  = false;
  bool passed = false;
  bool isPrescaled = true;

  int count =0;
  for (vstring::const_iterator iName = _triggerNames.begin(); iName != _triggerNames.end(); ++iName) {
    if (*iName == name) {
      isFound = true;
      if (triggerStatus & (0x1 << count)) 
        passed = true;
      if (hltPrescales[count] == 1) 
        isPrescaled = false;

      //cout<<"TRG dbg: iName = "<<name<<" matches to the name  "<<name<<endl;
      //cout<<"which is   number "<<count<<"  in the array of trigger names"<<endl;
      //cout<<"Is it passed? = "<<passed<<"    is it prescaled? = "<<isPrescaled<<endl;
      
      break;
    }
    ++count;
  }

  if (passed && isPrescaled)
    cout<<count<<"  Caution **: the trigger "<<name<<"   passed but it is Prescaled!  with a prescale"<<hltPrescales[count]<<endl;
  if(!isFound)
    cout<<"TRG **** Warning ***\n The trigger name you specified is not in the list of trigger names"<<endl;


  return passed;
}
