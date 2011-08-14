#ifndef _NUTRIGGER_H
#define	_NUTRIGGER_H

#include "TObject.h"
#include "TLorentzVector.h"

class TCTrigger : public TObject {
private:
  int valid_;
  int fired_;
  int prescale_;

public:
    TCTrigger();
    virtual ~TCTrigger();

    int isValid() {return valid_;}
    int hasFired() {return fired_;}
    int prescale() {return prescale_;}

    void setValid(int x) {valid_ = x;}
    void setFired(int x) {fired_ = x;}
    void setPrescale(int x) {prescale_ = x;}

    ClassDef(TCTrigger, 1);

};

#endif	/* _NUTRIGGER_H */


