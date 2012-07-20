#include <iostream>
#include <TChain.h>
#include <TClonesArray.h>
#include <TString.h>
#include <map>

#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TRandom3.h>


class rochcor {
 public:
  rochcor();
  ~rochcor();
  
  void momcor_mc(TLorentzVector& , TLorentzVector& , int, float, int);
  void momcor_data(TLorentzVector& , TLorentzVector& , int, float, int);
  
  void momcor_mc_muscle(TLorentzVector& , TLorentzVector& , int, float);
  void momcor_data_muscle(TLorentzVector& , TLorentzVector& , int, float);
  
  void momcor_mc_sidra(TLorentzVector& , TLorentzVector& , int, float);
  void momcor_data_sidra(TLorentzVector& , TLorentzVector& , int, float);
  
  void musclefit_data(TLorentzVector& , TLorentzVector&);
  
  float zptcor(float);
  int etabin(float);
  int phibin(float);
  
 private:
  
  TRandom3 eran1;
  TRandom3 eran2;

  TRandom3 sran1;
  TRandom3 sran2;
  
  
  //  static float netabin[9] = {-2.4,-2.1,-1.4,-0.7,0.0,0.7,1.4,2.1,2.4};
  static const double pi = 3.14159265358979323846;
  static const float netabin[9];

  static const float genm_smr = 9.09956e+01; //gen mass peak with eta dependent gaussian smearing => better match in Z mass profile vs. eta/phi
  static const float genm = 91.06; //gen mass peak without smearing => Z mass profile vs. eta/phi in CMS note

  static const float recmA = 9.09971e+01; //rec mass peak in MC (2011A)
  static const float drecmA = 9.09103e+01; //rec mass peak in data (2011A)
  static const float recmB = 9.10028e+01; //rec mass peak in MC (2011B)
  static const float drecmB = 9.09105e+01; //rec mass peak in data (2011B)
  
  //iteration2 after FSR : after Z Pt correction
  static const float deltaA = -3.04094e-06;
  static const float deltaAer = 7.68168e-07;
  
  static const float sfA = 43.4069;
  static const float sfAer = 1.50536;

  static const float deltaB = -3.04094e-06;
  static const float deltaBer = 7.68168e-07;
  
  static const float sfB = 43.4069;
  static const float sfBer = 1.50536;

  static const float delta_muscle = -2.89579e-06;
  static const float deltaer_muscle = 7.86676e-07;
  
  static const float sf_muscle = 44.3096;
  static const float sfer_muscle = 1.66431;
  
  static const float delta_sidra = -3.46001e-06;
  static const float deltaer_sidra = 7.7276e-07;
  
  static const float sf_sidra = 40.2147;
  static const float sfer_sidra = 1.34094;
 
  static const float apar = 1.0; //+- 0.002
  static const float bpar = -5.03313e-06; //+- 1.57968e-06
  static const float cpar = -4.41463e-05; //+- 1.92775e-06
  static const float d0par = -0.000148871; //+- 3.16301e-06
  static const float e0par = 1.59501; //+- 0.0249021
  static const float d1par = 7.95495e-05; //+- 1.12386e-05
  static const float e1par = -0.364823; //+- 0.17896
  static const float d2par = 0.000152032; //+- 5.68386e-06
  static const float e2par = 0.410195; //+- 0.0431732
 
  //---------------------------------------------------------------------------------------------
  
  static const float dcor_bfA[8][8];  
  static const float dcor_maA[8][8];
  static const float mcor_bfA[8][8];
  static const float mcor_maA[8][8];
  static const float dcor_bfAer[8][8];  
  static const float dcor_maAer[8][8];
  static const float mcor_bfAer[8][8];
  static const float mcor_maAer[8][8];

  static const float dcor_bfB[8][8];  
  static const float dcor_maB[8][8];
  static const float mcor_bfB[8][8];
  static const float mcor_maB[8][8];
  static const float dcor_bfBer[8][8];  
  static const float dcor_maBer[8][8];
  static const float mcor_bfBer[8][8];
  static const float mcor_maBer[8][8];

  static const float dcor_bf_muscle[8][8];  
  static const float dcor_ma_muscle[8][8];
  static const float dcor_bfer_muscle[8][8];  
  static const float dcor_maer_muscle[8][8];

  static const float dcor_bf_sidra[8][8];  
  static const float dcor_ma_sidra[8][8];
  static const float dcor_bfer_sidra[8][8];  
  static const float dcor_maer_sidra[8][8];

  //=======================================================================================================
  
  static const float dmavgA[8][8];  
  static const float dpavgA[8][8];  
  static const float mmavgA[8][8];  
  static const float mpavgA[8][8];

  static const float dmavgB[8][8];  
  static const float dpavgB[8][8];  
  static const float mmavgB[8][8];  
  static const float mpavgB[8][8];
  
  static const float dmavg_muscle[8][8];  
  static const float dpavg_muscle[8][8];  
  static const float dmavg_sidra[8][8];  
  static const float dpavg_sidra[8][8];  
  
  //===============================================================================================
  //parameters for Z pt correction
  static const int nptbins=84;
  static const float ptlow[85];    
  
  static const float zptscl[84];
  static const float zptscler[84];

};
  
