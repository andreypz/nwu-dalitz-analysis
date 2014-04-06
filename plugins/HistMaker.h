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

   virtual void FillHistosFull(Int_t n, Double_t w,
			       TCPhysObject , TCPhysObject , TCPhysObject , TCPhysObject , TCPhysObject,
			       string s="");
   virtual void MakeMuonPlots(TCMuon mu, TVector3 *pv);
   virtual void MakePhotonPlots(TCPhoton ph);
   virtual void MakeZeePlots(TCPhoton , TCPhoton );
   virtual void MakePhotonEnergyCorrPlots(TCPhoton , Float_t , Float_t );

 private:
  HistManager * hists;

};


#endif
