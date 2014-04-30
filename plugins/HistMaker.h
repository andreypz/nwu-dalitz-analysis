#ifndef _HistMaker_H
#define _HistMaker_H

#include <cstring>
#include "HistManager.h"

#include "../interface/TCPhoton.h"
#include "../interface/TCEGamma.h"
#include "../interface/TCElectron.h"
#include "../interface/TCPhoton.h"
#include "../interface/TCMuon.h"


class HistMaker {
 public:
  HistMaker(HistManager * h);
  virtual ~HistMaker();

   virtual void FillHistosFull(Int_t n, Double_t w, string s="");
   virtual void MakeMuonPlots(TCMuon mu, TVector3 *pv);
   virtual void MakePhotonPlots(TCPhoton ph);
   virtual void MakeZeePlots(TCPhoton , TCPhoton );
   virtual void MakePhotonEnergyCorrPlots(TCPhoton , Float_t , Float_t );
   virtual void SetLeptons(TCPhysObject l1, TCPhysObject l2);
   virtual void SetGamma(TCPhysObject g);

 private:
  HistManager * hists;

  TCPhysObject _lPt1, _lPt2, _gamma;

};


#endif
