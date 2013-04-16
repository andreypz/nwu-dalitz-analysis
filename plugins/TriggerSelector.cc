#include "TriggerSelector.h"

TriggerSelector::TriggerSelector() {}

TriggerSelector::~TriggerSelector() {}

TriggerSelector::TriggerSelector(string type, string dataPeriod, vstring triggerNames)
{
    _type           = type;
    _dataPeriod     = dataPeriod;
    _triggerNames   = triggerNames;
    _isRealData     = false;

}

void TriggerSelector::SelectTrigger(string name, unsigned triggerStatus, UInt_t* hltPrescales, Bool_t& isFound, Bool_t& isPassed, Int_t& prescale )
{

  prescale = 90; 
  isFound  = false;
  isPassed = false;
  //bool isPrescaled = true;

  int count = 0;
  for (vstring::const_iterator iName = _triggerNames.begin(); iName != _triggerNames.end(); ++iName) {
    if (*iName == name) {
      isFound = true;
      if (triggerStatus & (ULong64_t(0x1) << count)) 
        isPassed = true;
      
      //cout<<"TRG dbg: iName = "<<name<<" matches to the name  "<<name<<endl;
      //cout<<"which is   number "<<count<<"  in the array of trigger names"<<endl;
      //cout<<"Is it passed? = "<<passed<<"    is it prescaled? prescsle = "<<hltPrescales[count]<<endl;
      
      break;
    }
    ++count;
  }

  if (isPassed)
    prescale = hltPrescales[count]; 

  //if (passed && prescale!=1)
  //cout<<count<<"  Caution **: the trigger "<<name<<"   passed but it is Prescaled!  with a prescale = "<<hltPrescales[count]<<"  "<<prescale<<endl;

  //if(!isFound)
  // cout<<"TRG **** Warning ***\n The trigger name you specified is not in the list of trigger names"<<endl;


  //return prescale;
}
