// -*- C++ -*-
//
// Package:    Dummy
// Class:      Dummy
// 
/**\class Dummy Dummy.cc UserCode/Dummy/src/Dummy.cc

 Description: Dummy analyzer. For CRAB testing.

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Andrey Pozdnyakov
//         Created:  Wed Jul  8 11:08:04 CEST 2009
// $Id: Dummy.cc,v 1.4 2011/05/08 16:47:53 andrey Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"

//#include <fstream>
//#include <map>

using namespace edm;
using namespace std;


class Dummy : public edm::EDAnalyzer {
   public:
      explicit Dummy(const edm::ParameterSet&);
      ~Dummy();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------


  //string outputFileName_;
  TH1F *EV;
  //  TTree *gTree;
  float cs;
  edm::Service<TFileService> fs;

};





Dummy::Dummy(const edm::ParameterSet& iConfig)

{

//  outputFileName_=iConfig.getParameter<std::string>("outputFileName");
}


Dummy::~Dummy()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
Dummy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Bool_t isRealData = iEvent.isRealData();

   double sigma =10.;

   if (!isRealData)
     {
       //That stuff doesn't wark at 310pre10. Uncomment later.
       edm::Handle<GenRunInfoProduct> gi;
       //iEvent.getRun().getByType(gi);
       iEvent.getRun().getByLabel("generator", gi );
       sigma = gi->externalXSecLO();
       
       //  double external_cross_section = gi->external_cross_section();
       //  double filter_eff = gi->filter_efficiency();
       //  cout  << "Cross Section:  "<< sigma <<"  fil_eff: "<<filter_eff<<"   ext_cs: "<<external_cross_section<<endl;
       //double sigma = gi->cross_section();
       
     }
  cs= sigma;
  //cout<<"Hello World  "<<sigma<<endl;
  EV -> Fill(log10(sigma));
  //  gTree -> Fill();
  

}


// ------------ method called once each job just before starting event loop  ------------
void 
Dummy::beginJob()
{
 // Open the histogram file and book some associated histograms
  //m_file = new TFile(outputFileName_.c_str(),"RECREATE");

  //  gTree = fs->make<TTree>("gTree","gTree");
  // gTree->Branch("cs", &cs, "cs/F");
  EV = fs->make<TH1F>("EV", "cross section, log scale", 5000, -10., 2.);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
Dummy::endJob() {
  //Write out the histogram file.
  //EV -> Write();
  //m_file->Write();
  //m_file->Close();

}

//define this as a plug-in
DEFINE_FWK_MODULE(Dummy);
