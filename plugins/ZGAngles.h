#ifndef _ZGAngles_H
#define _ZGAngles_H


#include "TLorentzVector.h"

// c++ libraries
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

class ZGAngles {
 public:
  ZGAngles();
  virtual ~ZGAngles();
  void SetAngles(TLorentzVector l1, TLorentzVector l2, TLorentzVector g);
  void GetAngles(float& cos1, float& cos2, float& phi, float& cos3);
  void GetAngles(TLorentzVector l1, TLorentzVector l2, TLorentzVector g, float& cos1, float& cos2, float& phi, float& cos3);
  
 private:
  float _costheta_lplus, _costheta_lminus, _phi, _cosTheta;
};


#endif
