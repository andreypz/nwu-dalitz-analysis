//#include "TVector2.h"
//#include "TMath.h"


double projectedMET(Float_t f_MET, Float_t f_MET_phi, TLorentzVector f_Lepton1, TLorentzVector f_Lepton2)
{
  Float_t deltaPhiMin =   TMath::Min( fabs(TVector2::Phi_mpi_pi(f_Lepton1.Phi() - f_MET_phi)), fabs(TVector2::Phi_mpi_pi(f_Lepton2.Phi() - f_MET_phi)));
  //cout<<"dbg:  "<<deltaPhiMin<<endl;

  Float_t f_projMET = 0;
  if (deltaPhiMin > TMath::Pi()/2) f_projMET = f_MET;
  else f_projMET = f_MET*sin(deltaPhiMin);
  return f_projMET;
}

double ZprojectedMET(Float_t f_MET, Float_t f_MET_phi, TLorentzVector f_Lepton1, TLorentzVector f_Lepton2)
{
  TLorentzVector di = f_Lepton1 + f_Lepton2;
  Float_t deltaPhiMin =   fabs(TVector2::Phi_mpi_pi(di.Phi() - f_MET_phi));

  Float_t f_projMET = 0;
  if (deltaPhiMin > TMath::Pi()/2) f_projMET = f_MET;
  else f_projMET = f_MET*sin(deltaPhiMin);
  return f_projMET;
}

double PUCorrectedMET(float met, int nPV, string metType, Bool_t isData)
{
  Float_t sigma0=0, shift=0, sigmaT=0, shift0=0, puCorMet=0;
  if (metType == "pfMet") {
    if (isData) {
      sigma0 = 5.868;
      shift0 = 8.231;
      shift  = 0.5761*(nPV - 1);
      sigmaT = sqrt(pow(sigma0, 2) + (nPV - 1)*pow(2.112, 2));
    } else {
      sigma0 = 4.569;
      shift0 = 6.282;
      shift  = 0.8064*(nPV - 1);
      sigmaT = sqrt(pow(sigma0, 2) + (nPV - 1)*pow(2.387, 2));
    }
  }
  puCorMet = (met - (shift0 + shift))*sigma0/sigmaT + shift0;
  return puCorMet;
}

TLorentzVector GetReducedMET(TLorentzVector sumJet, TLorentzVector lep1, TLorentzVector lep2, TLorentzVector metP4, int version)
{
  TLorentzVector Q  = lep1 + lep2;
  float bisectorPhi = min(lep1.Phi(), lep2.Phi()) + lep1.DeltaPhi(lep2)/2;

  TVector2 projectDilepton(Q.Px(), Q.Py());
  TVector2 projectSumJet(sumJet.Px(), sumJet.Py());
  TVector2 projectMET(metP4.Pt()*cos(metP4.Phi()), metP4.Pt()*sin(metP4.Phi()));

  //TVector2 delta;                                                                                                                                        
  //TLorentzVector Thrust = lep1 - lep2;                                                                                                                   
  //if (fabs(lep1.DeltaPhi(lep2)) > TMath::Pi()/2) {                                                                                                       
  //  delta.Set(0., projDilepton.Rotate(-Thrust.Phi()).Py() + projDilepton.Py());                                                                          
  //} else {                                                                                                                                               
  //  delta.Set(0, projDilepton.Rotate(-Q.Phi()).Py() + projDilepton.Py());                                                                                
  //}                                                                                                                                                      

  if (version == 1) {
    projectDilepton = projectDilepton.Rotate(-bisectorPhi);
    projectSumJet   = projectSumJet.Rotate(-bisectorPhi);
    projectMET      = projectMET.Rotate(-bisectorPhi);
  } else if (version == 2) {
    projectDilepton = projectDilepton.Rotate(-Q.Phi());
    projectSumJet   = projectSumJet.Rotate(-Q.Phi());
    projectMET      = projectMET.Rotate(-Q.Phi());
  }
  TVector2 unclustered = -1*projectMET;                        //projectDilepton - 1.*(projectMET + projectDilepton) + delta;                              
  TVector2 clustered   = projectDilepton + 1.*projectSumJet;   // + delta;                                                                                 

  TVector2 reducedMET = TVector2((fabs(unclustered.Px()) < fabs(clustered.Px()) ? unclustered.Px() : clustered.Px()),
				 (fabs(unclustered.Py()) < fabs(clustered.Py()) ? unclustered.Py() : clustered.Py()));

  return TLorentzVector(reducedMET.Px(), reducedMET.Py(), 0, reducedMET.Mod());
}


