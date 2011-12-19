// $Id: ZedEventsLibrary.cc,v 1.4 2011/12/14 17:30:41 andrey Exp $

#include "ZedEventsLibrary.h"

ZedEventsLibrary::ZedEventsLibrary() {}

ZedEventsLibrary::~ZedEventsLibrary() {}

ZedEventsLibrary::ZedEventsLibrary(string selection, bool makeFromData)
{
  _selection    = selection;
  //  _dataPeriod   = dataPeriod;
  _makeFromData = makeFromData;
  
  if(_makeFromData){
    if(_selection=="muGamma" || _selection=="muon")
      _inFile = new TFile("/uscms_data/d2/andreypz/ZedLib/m_Data_1.root", "OPEN");
    else if(_selection=="eGamma" || _selection=="electron")
      _inFile = new TFile("/uscms_data/d2/andreypz/ZedLib/m_Data_2.root", "OPEN");
    else cout<<"this is not supported: "<<_selection<<endl;
  }
  else{
    if(_selection=="muGamma" || _selection=="muon")
      _inFile = new TFile("/uscms_data/d2/andreypz/ZedLib/m_Zjets_1.root", "OPEN");
    else if(_selection=="eGamma" || _selection=="electron")
      _inFile = new TFile("/uscms_data/d2/andreypz/ZedLib/m_Zjets_2.root", "OPEN");
    else cout<<"this is not supported: "<<_selection<<endl;
  }
  

  _myRandom = new TRandom3();                                                                                                                             

  _qt_bins[0] =  0;
  _qt_bins[1] =  25;
  _qt_bins[2] =  26;
  _qt_bins[3] =  28;
  _qt_bins[4] =  32;
  _qt_bins[5] =  40;
  _qt_bins[6] =  56;
  _qt_bins[7] =  88;
  _qt_bins[8] =  152;
  _qt_bins[9] =  280;

  _eta_bins[0] =  -15.0;
  _eta_bins[1] =  -1.5;
  _eta_bins[2] =  -1.0;
  _eta_bins[3] =  0.0;
  _eta_bins[4] =  1.0;
  _eta_bins[5] =  1.5;
  _eta_bins[6] =  15.0;
  
}


pair <int, int> ZedEventsLibrary::GetEtaQtBin(Float_t myEta, Float_t myQT){
  
  //  _qt_bins[10] =  {0,25,26,28,32,40,56,88,152,280};
  //_eta_bins[7] =  {-15, -1.5, -1.0, 0, 1.0, 1.5, 15};
  
  pair <int, int> myBins = make_pair(0,0); 
  // Need to know the size of qt_bins array
  int len = (int)sizeof(_eta_bins)/sizeof(float); 
  Int_t ebin = -1;
  for(Int_t e = len-2; e>=0; e--){
    if(myEta>=_eta_bins[e]){
      ebin = e;
      //bin is assigned - exit for loop
      break;
    }
  }
  
  //cout<<"Assigning ebin = "<<ebin<<"   in eta = ["<< eta_bins[ebin]<<", "<<eta_bins[ebin+1]<<"]"<<endl;
  //loop over qt_bins and see where our events belongs  
  
  len=sizeof(_qt_bins)/sizeof(int); 
  Int_t qbin = -1;
  for(Int_t q = len-1; q>=0; q--){
    if(myQT>=_qt_bins[q]){
      qbin = q;
      //bin is assigned - exit for loop
      break;
    }
  }
  //cout<<"Assigning qbin = "<<qbin<<"   in qt = ["<< qt_bins[qbin]<<", ]"<<endl;

  if(ebin==-1 || qbin==-1) cout<<"Warning: ebin or qbin is not assigned!  eta="<<myEta<<"   qt="<<myQT<<endl;
  myBins = make_pair(ebin,qbin); 
  return myBins;

}


void ZedEventsLibrary::CreateLibrary(){

  cout<<"\n    Doing the library thing now"<<endl;
  cout<<"input file:  "<<endl;
  _inFile->Print();

  TTree *libTree = (TTree*)_inFile->Get("Anton/kinTree");

  TLorentzVector  *lep1=0, *lep2=0;
  //libTree->SetBranchStatus("*",0);
  //libTree->SetBranchStatus("lep1", 1); //doesn't work if *,0 set
  //libTree->SetBranchStatus("lep2", 1);
  libTree->SetBranchAddress("lep1", &lep1);
  libTree->SetBranchAddress("lep2", &lep2);
  //libTree->Print();
      
  Int_t nentries = libTree->GetEntries();

  for(Int_t j = 0; j<nentries; j++){
    libTree-> GetEntry(j);

    Float_t myQT  = (*lep1 + *lep2).Pt() ;
    Float_t myEta = (*lep1 + *lep2).Eta() ;
    //cout<<j<<"  myQT= "<<myQT<<endl;
    pair <int, int> myBin = GetEtaQtBin(myEta,myQT);

    pair <TLorentzVector, TLorentzVector> diLeptonPair;
    diLeptonPair = make_pair (*lep1, *lep2);
    //cout<<"Orig 1: "<<lep1->Pt()<<" 2: "<<lep2->Pt()<<endl;
    //cout<<"Pair 1: "<< diLeptonPair.first.Pt()<<" 2: "<< diLeptonPair.second.Pt()<<endl;

    _diMap[myBin.first][myBin.second].push_back(diLeptonPair);
  }

  //Make a hist with events in the library
  //for(Int_t q=0; q<10; q++){
  //  for(Int_t e=0; e<4; e++)
  //    evt_libQt -> SetBinContent(q+1, e+1,diMap[e][q].size());
  // }
  cout<<"\n  Total events in the library: "<<nentries<<endl;
  //cout<<"in bins:  "<<endl;
  _inFile->Close();
  cout<<"\nDone with the library. Moving on. \n"<<endl;
}

pair<TLorentzVector, TLorentzVector> ZedEventsLibrary::GetLeptonsFromLibrary(TLorentzVector origZP4){

  //Find out which qT bin we are
  pair <int, int> myBin = GetEtaQtBin(origZP4.Eta(), origZP4.Pt());
  Int_t etaBin = myBin.first;
  Int_t qtBin  = myBin.second;

  //Look how many events are in this bin:
  Int_t mySize = _diMap[etaBin][qtBin].size();                                                                                                            
  //cout<<"There are --"<<mySize<<"--  events in bin[eta,pt] = ["<<etaBin<<", "<<qtBin<<"]"<<endl;                                                       
  if (mySize==0) {cout<<"No events found in the library for that bin - ["<<etaBin<<", "<<qtBin<<"]  -! return kTRUE"<<endl;}               

  TLorentzVector  lep1, lep2, diLepton;
 
  //cout<<l<<"\n  ---- Original Gamma/Z ----"<<endl;                                                                                                 
  //cout<<"\t qt= "<<origZP4.Pt()<<"\t eta= "<<origZP4.Eta()<<"\t phi= "<<origZP4.Phi()<<"\t M= "<<origZP4.M()<<endl;
  
  //Try 20 times (need both leptons within eta<2.4 and pt>20 after the boosting):
  for(Int_t l=0; l<20; l++){                                                                                                                             
    //Throw a random number 
    Int_t rndm = _myRandom->Integer(mySize);                                                                                                           
    //cout<<l<<"  Random number out of --"<<mySize<<"--   is going to be: *"<<rndm<<"*"<<endl;
    pair<TLorentzVector, TLorentzVector> myPair = _diMap[etaBin][qtBin].at(rndm);                                                                         

    lep1 = myPair.first;
    lep2 = myPair.second;
    diLepton = lep1 + lep2;

    //First, boost them to rest-frame:
    TVector3 b1 = diLepton.BoostVector();                                                                                                                
    lep1.Boost(-b1);                                                                                                                        
    lep2.Boost(-b1);                                                                                                                                  
    //cout<<"\nAfter boost to rest-frame:\n"<<endl;                                                                                                      

    //Make a vector with Z/Gamma direction but the mass of  diLepton
    TLorentzVector origZP4_Mass(0,0,0,0);
    origZP4_Mass.SetPtEtaPhiM(origZP4.Pt(), origZP4.Eta(), origZP4.Phi(), diLepton.M());

    //Now, boost the leptons in Gamma/Z direction:
    TVector3 b2 = origZP4_Mass.BoostVector();

    lep1.Boost(b2);                                                                                                                                   
    lep2.Boost(b2);
    //Check if it satisfy our main selection cuts                                                                                  
    if(lep1.Pt()>20 && lep2.Pt()>20 && fabs(lep1.Eta())<2.4 && fabs(lep2.Eta())<2.4)                                                         
      break; //found a good pair of leptons
    //if not found in 20 attempts - it will use the latest tried
  }

  //diLepton = lep1 + lep2;
  //cout<<"\n After final boost:"<<endl;                              
  //cout<<"\t qt= "<<diLepton.Pt()<<"\t eta= "<<diLepton.Eta()<<"\t phi= "<<diLepton.Phi()<<"\t M= "<<diLepton.M()<<endl; 
  
  //Swap the order of leptons if pt's changed                                                                                                            
  if (lep2.Pt() > lep1.Pt()) {                                                                                                                       
    //cout<<"swaping leptons pt: "<<endl;
    TLorentzVector tempLep = lep2;                                                                                                                    
    lep2 = lep1;                                                                                                                                   
    lep1 = tempLep;                                                                                                                                   
  }                
  return make_pair(lep1, lep2);
}
