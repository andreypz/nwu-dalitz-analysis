#include "TFile.h"
#include "TString.h" 
#include "TH1F.h"
#include "TCanvas.h"
#include <iostream>

Float_t PT_LOW = 30;
Float_t PT_HIGH = 400;

Int_t REBIN_GAMMA = 10;

// files with histograms from mumu,ee needed for scaling
//TString eeFileName = "forAnton_Data_ee.root";
//TString mumuFileName = "forAnton_Data_mumu.root";
//TString gammaFileName = "photonResult_070811_1.root";
//TString outFileName = "photonWeight_070811_1.root";

//TString eeFileName = "forAnton_Data_sbtr_ee.root";
//TString mumuFileName = "forAnton_Data_sbtr_mumu.root";

TString eeFileName = "forAnton_Data_sbtr_2.root";
TString mumuFileName = "forAnton_Data_sbtr_1.root";

TString gammaFileName = "PhotonRun_110714-1.root";

TString outFileName = "photonWeight_110714-1.root";




TString eeHistoName = "di_qt_4";
TString mumuHistoName = "di_qt_4";
TString gammaHistoName = "photon_pt"; // for gettin the weights
//TString gammaHistoName = "photon1_weighted_pt";

TString llMassHistoName = "di_mass_4";

using namespace std;

void scale() {

TFile* eeFile = new TFile(eeFileName);
TFile* mumuFile = new TFile(mumuFileName);
TFile* gammaFile = new TFile(gammaFileName);

TH1F* h1_ee = (TH1F*) eeFile->Get(eeHistoName)->Clone();
TH1F* h1_mumu = (TH1F*) mumuFile->Get(mumuHistoName)->Clone();

// combined ee+mumu for better statistics
// may be useful if the ee, mumu do not show difference in shape
TH1F* h1_ll = (TH1F*) eeFile->Get(eeHistoName)->Clone(); 
h1_ll->Add(h1_mumu);

// Gamma shape for normalization of the pt spectrum
// allow for the possibility to have separate normalization shapes
// for use in ee, mumu. For now they are common and the same histogram is used.
TH1F* h1_gamma1 = (TH1F*) gammaFile->Get(gammaHistoName)->Clone();
TH1F* h1_gamma2 = (TH1F*) gammaFile->Get(gammaHistoName)->Clone();


// get the ee, mumu mass distributions so that they can be copied to the
// weight file. These are used for "photon mass" assignment (to mimic Z) 
// if we do not use an analytic function. 

TH1F* h1_photonMassFromEe = (TH1F*) eeFile->Get(llMassHistoName)->Clone();
TH1F* h1_photonMassFromMumu = (TH1F*) mumuFile->Get(llMassHistoName)->Clone();
h1_photonMassFromEe->SetName("h1_photonMassFromEe");
h1_photonMassFromMumu->SetName("h1_photonMassFromMumu");

h1_photonMassFromEe->SetTitle("h1_photonMassFromEe");
h1_photonMassFromMumu->SetTitle("h1_photonMassFromMumu");

// Just in case...
h1_ee->Sumw2();
h1_mumu->Sumw2();
h1_gamma1->Sumw2();
h1_gamma2->Sumw2();

TH1F* h1_eeScale = (TH1F*) h1_ee->Clone();
TH1F* h1_mumuScale = (TH1F*) h1_mumu->Clone();
// if the ee, mumu are combined to get shape
// Additional ovarall normalization is needed if this is used for 
// ee, mumu based on number of events.
// Thi is not used for now
TH1F* h1_llScale = (TH1F*) h1_ll->Clone();

h1_eeScale->Divide(h1_gamma1);
h1_mumuScale->Divide(h1_gamma1);
h1_llScale->Divide(h1_gamma1);

h1_eeScale->SetName("h1_eeScale");
h1_eeScale->SetTitle("weight for ee sample; p_{T} (GeV); ratio photon/DiLepton");
h1_mumuScale->SetName("h1_mumuScale");
h1_mumuScale->SetTitle("weight for #mu#mu sample; p_{T} (GeV); ratio photon/DiLepton");
h1_llScale->SetName("h1_llScale");
h1_llScale->SetTitle("Combined ee+#mu#mu; p_{T} (GeV); ratio photon/DiLepton");


// Sanity check
TCanvas* eeCan = new TCanvas();
eeCan->cd();
h1_eeScale->Draw("E");

TCanvas* mumuCan = new TCanvas();
mumuCan->cd();
h1_mumuScale->Draw("E");

TCanvas* llCan = new TCanvas();
llCan->cd();
h1_llScale->Draw("E");


// compare the distributions for ee,mumu
TH1F* h1_eeMummuRatio = (TH1F*) h1_ee->Clone();
h1_eeMummuRatio->Divide(h1_mumu);
h1_eeMummuRatio->SetTitle("Ratio of ee/mumu distributions; p_{T} (GeV); ratio");

TCanvas* canEeMumuRatio = new TCanvas();
canEeMumuRatio->cd();
h1_eeMummuRatio->Draw();




// This file contains the weights and all other information
// needed to plot distributions with weighted events
// Amke sure you include the proper name in the analyzer
TFile* outFile = new TFile(outFileName.Data(), "RECREATE");
outFile->cd();
h1_eeScale->Write();
h1_mumuScale->Write();
h1_llScale->Write();
// these are the copies from Andrey's files
h1_photonMassFromEe->Write();
h1_photonMassFromMumu->Write();

outFile->Close();

}
