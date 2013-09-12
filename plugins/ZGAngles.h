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
  void GetAngles(double& cos1, double& cos2, double& phi, double& cos3);
  void GetAngles(TLorentzVector l1, TLorentzVector l2, TLorentzVector g, double& cos1, double& cos2, double& phi, double& cos3);
  double GetCos1();
  double GetCos2();
  double GetCosTheta();
  double GetPhi();

 private:
  double _costheta_lplus, _costheta_lminus, _phi, _cosTheta;
};


#endif
