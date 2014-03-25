#ifndef _histo_H
#define _histo_H

#include "TLorentzVector.h"
#include "TVector3.h"
#include "../interface/TCPhoton.h"
#include "../interface/TCElectron.h"
#include "../interface/TCMuon.h"
#include "HistManager.h"

// This macro defines various functions to fill the histogramms
// To be used in multiple analyzers.

void MakeElectronPlots(HistManager *hists, TCElectron el, string dir, Double_t rho25Factor)
{
  //cout<<"Making electron plots in dir="<<dir<<endl;


  hists->fill1DHist(el.InverseEnergyMomentumDiff(), dir+"_el_fabsEPDiff",";|1/E - 1/p|", 100, 0, 0.1, 1, dir);

  hists->fill1DHist(el.R9(),    dir+"_el_R9",   ";R9",     100, 0, 1, 1,dir);
  hists->fill1DHist(el.MvaID(), dir+"_el_mvaID",";mva ID", 100, -1,1, 1,dir);
  hists->fill1DHist(el.FBrem(), dir+"_el_fbrem","; fbrem", 100, 0, 1, 1,dir);
  hists->fill1DHist(el.HadOverEm(),  dir+"_el_HadOverEm", ";HadOverEm",  100, 0, 0.1,     1,dir);
  hists->fill1DHist(el.SCDeltaPhi(), dir+"_el_SCDeltaPhi",";SCDeltaPhi", 100, -0.2, 0.2,  1,dir);
  hists->fill1DHist(el.SCDeltaEta(), dir+"_el_SCDeltaEta",";SCDeltaEta", 100, -0.02, 0.02,1,dir);

  hists->fill1DHist(el.PtError()/el.Pt(),     dir+"_el_ptErrorOverPt",";ptErrorOverPt", 100, 0, 1,   1,dir);
  hists->fill1DHist(el.SigmaIEtaIEta(),       dir+"_el_SigmaIEtaIEta",";SigmaIEtaIEta", 100, 0, 0.1, 1,dir);
  hists->fill1DHist(el.SigmaIPhiIPhi(),       dir+"_el_SigmaIPhiIPhi",";SigmaIPhiIPhi", 100, 0, 0.1, 1,dir);
  hists->fill1DHist(el.NormalizedChi2(),      dir+"_el_gsfChi2",      ";gsfChi2",       100, 0, 3,   1,dir);
  hists->fill1DHist(el.NormalizedChi2Kf(),    dir+"_el_kfChi2",       ";kfChi2",        100, 0, 3,   1,dir);
  hists->fill1DHist(el.DeltaEtaSeedCluster(), dir+"_el_dEtaAtCalo",   ";dEtaAtCalo",    100, -0.05, 0.05, 1,dir);
  hists->fill1DHist(el.DeltaPhiSeedCluster(), dir+"_el_dPhiAtCalo",   ";dPhiAtCalo",    100, -0.05, 0.05, 1,dir);
  hists->fill1DHist(el.PreShowerOverRaw(),    dir+"_el_preShowerOverRaw",";preShowerOverRaw", 100, 0, 0.3,1,dir);
  hists->fill1DHist(el.SCEtaWidth(),          dir+"_el_SCEtaWidth",   ";SCEtaWidth",    100, 0, 0.04, 1,dir);
  hists->fill1DHist(el.SCPhiWidth(),          dir+"_el_SCPhiWidth",   ";SCPhiWidth",    100, 0, 0.1,  1,dir);
  hists->fill1DHist(1 - el.E1x5()/el.E5x5(),  dir+"_el_ome1x5oe5x5",  ";1 - e1x5/e5x5", 100, 0, 1,    1,dir);
  hists->fill1DHist(el.EoP(),                 dir+"_el_EoP",    ";EoP",     100, 0, 6,   1,dir);
  hists->fill1DHist(el.EoPout(),              dir+"_el_EoPout", ";EoPout",  100, 0, 6,   1,dir);
  hists->fill1DHist(el.IP3d(),                dir+"_el_ip3d",   ";ip3d",    100, 0, 0.04,1,dir);
  hists->fill1DHist(el.IP3dSig(),             dir+"_el_ip3dSig",";ip3dSig", 100, 0, 3,   1,dir);
  hists->fill1DHist(el.NumberOfValidHits(),   dir+"_el_NumberOfValidHits", ";NumberOfValidHits",   25, 0,25,   1,dir);
  hists->fill1DHist(el.TrackerLayersWithMeasurement(), dir+"_el_TrackerLayersWithMeasurement",
		    ";TrackerLayersWithMeasurement", 20, 0,20,1,dir);

  hists->fill1DHist(el.PassConversionVeto(), dir+"_el_passConv",";pass conv veto", 3, 0, 3, 1, dir);
  hists->fill1DHist(el.ConversionMissHits(), dir+"_el_convMissHits",";conv miss hits", 10, 0, 10, 1, dir);

  hists->fill1DHist(el.ConversionDist(),   dir+"_el_convDist",  ";conv Dist",   200, 0, 50, 1, dir);
  hists->fill1DHist(el.ConversionDcot(),   dir+"_el_convDcot",  ";conv Dcot",   200, -5, 5, 1, dir);
  hists->fill1DHist(el.ConversionRadius(), dir+"_el_convRadius",";conv Radius", 200, 0, 50, 1, dir);

  hists->fill1DHist(el.GetTracks().size(),   dir+"_el_nTracks", ";nTracks",   10, 0,10,   1,dir+"_tracks");

  if (el.GetTracks().size()>=2)

    {
      Float_t mll = (el.GetTracks()[0]+el.GetTracks()[1]).M();
      hists->fill1DHist(mll, dir+"_el_mll_full", ";mll",  200, 0,120, 1,dir);
      hists->fill1DHist(mll, dir+"_el_mll_low",  ";mll",  200, 0,20,  1,dir);
      hists->fill1DHist(mll, dir+"_el_mll_jpsi", ";mll",  200, 1, 6,  1,dir);

    }
  for (Int_t i=0; i<el.GetTracks().size() && i<3; i++){

    string subdir = dir+"_tracks";
    hists->fill1DHist(el.GetTracks()[i].Pt(), subdir+Form("_el_track_%i",i)+"_Pt","Pt", 200, 0,100, 1,subdir);
    hists->fill1DHist(el.GetTracks()[i].PtError()/el.GetTracks()[i].Pt(), subdir+Form("_el_track_%i",i)+"_PtError","Pt Error /Pt", 200, 0,100, 1,subdir);
    hists->fill1DHist(el.GetTracks()[i].NormalizedChi2(), subdir+Form("_el_track_%i",i)+"_NormalizedChi2","NormalizedChi2", 100, 0,3, 1,subdir);

    TCTrack::ConversionInfo inf = el.GetTracks()[i].GetConversionInfo();
    hists->fill1DHist(inf.isValid, subdir+Form("_el_conv_%i",i)+"_isValid","isValid", 3, 0,3, 1,subdir);
    if (inf.isValid){
      hists->fill1DHist(inf.nHitsMax, subdir+Form("_el_conv_%i",i)+"_nHitsMax","nHitsMax", 20, 0,20, 1,subdir);
      hists->fill1DHist(inf.vtxProb,  subdir+Form("_el_conv_%i",i)+"_vtxProb", "vtxProb",  200,0, 1, 1,subdir);
      hists->fill1DHist(inf.lxyPV,    subdir+Form("_el_conv_%i",i)+"_lxyPV",   "lxyPV",    200,0,50, 1,subdir);
      hists->fill1DHist(inf.lxyBS,    subdir+Form("_el_conv_%i",i)+"_lxyBS",   "lxyBS",    200,0,50, 1,subdir);
    }

  }

  Double_t eleISO = (el.PfIsoCharged() + TMath::Max(0.0, (Double_t)(el.PfIsoPhoton() + el.PfIsoNeutral() - rho25Factor*el.EffArea())))/el.Pt();

  hists->fill1DHist(el.PfIsoCharged(), dir+"_el_iso_Charged",";pfChIso_R04",  100, 0, 30, 1, dir);
  hists->fill1DHist(el.PfIsoPhoton(),  dir+"_el_iso_Photon", ";pfPhoIso_R04", 100, 0, 20, 1, dir);
  hists->fill1DHist(el.PfIsoNeutral(), dir+"_el_iso_Neutral",";pfNeuIso_R04", 100, 0, 50, 1, dir);
  hists->fill1DHist(eleISO, dir+"_el_iso",";iso", 100, 0, 1, 1, dir);
  //hists->fill1DHist(el.IdMap(""), dir+"_el_",";", 3, 0, 3, 1, dir);
  //hists->fill1DHist(el., dir+"_el_",";", 3, 0, 3, 1, dir);

}

void FillHistosFull(HistManager *hists, Int_t num, Double_t weight,
				 TCPhysObject l1, TCPhysObject l2,
				 TCPhysObject lPt1, TCPhysObject lPt2, TCPhysObject gamma, string dir)
{
  Float_t mllMax = 25;
  char * d = new char [dir.length()+1];
  std::strcpy (d, dir.c_str());

  TLorentzVector diLep = l1+l2;
  TLorentzVector tri = diLep + gamma;

  TLorentzVector lPt1Gamma = lPt1+gamma;
  TLorentzVector lPt2Gamma = lPt2+gamma;

  TLorentzVector gammaCM = gamma;
  TLorentzVector diLepCM = diLep;
  TVector3    b1  = -tri.BoostVector();
  gammaCM.Boost(b1);
  diLepCM.Boost(b1);


  hists->fill1DHist(weight,     Form("weight_%s_cut%i",         d, num),";weight",  200, 0,3, 1, dir);

  hists->fill1DHist(tri.M(),     Form("tri_mass_%s_cut%i",         d, num),";M(ll#gamma)",  100,  0,200, weight, dir);
  hists->fill1DHist(tri.M(),     Form("tri_mass80_%s_cut%i",       d, num),";M(ll#gamma)",  100, 80,200, weight, dir);
  hists->fill1DHist(tri.M(),     Form("tri_mass_longTail_%s_cut%i",d, num),";M(ll#gamma)",  100,  0,400, weight, dir);
  hists->fill1DHist(tri.M(),     Form("tri_mass_Z_%s_cut%i",       d, num),";M(ll#gamma)",  100, 70,110, weight, dir);
  hists->fill1DHist(tri.M(),     Form("tri_mass_FIT_%s_cut%i",     d, num),";M(ll#gamma)",  100,110,170, weight, dir);
  hists->fill1DHist(tri.M(),     Form("tri_mass_mH125_%s_cut%i",   d, num),";M(ll#gamma)",  100,115,135, weight, dir);

  hists->fill1DHist(lPt1Gamma.M(), Form("lPt1Gamma_%s_cut%i",d, num),";M(l1#gamma)",  100, 0,400,  weight, dir);
  hists->fill1DHist(lPt2Gamma.M(), Form("lPt2Gamma_%s_cut%i",d, num),";M(l2#gamma)",  100, 0,400,  weight, dir);

  hists->fill2DHist(lPt1Gamma.M2(), lPt2Gamma.M2(),Form("h2D_dalitzPlot_%s_cut%i",  d, num),
		    ";M(l1#gamma)^{2};M(l2#gamma)^{2}", 100,0,20000, 100, 0,10000,  weight, dir);

  double rotation = atan(-1.);
  double sin_rotation = sin(rotation);
  double cos_rotation = cos(rotation);
  double x_prime = (lPt1Gamma.M2()*cos_rotation - lPt2Gamma.M2()*sin_rotation) + (1.-cos_rotation)*pow(125.,2);
  double y_prime = lPt1Gamma.M2()*sin_rotation + lPt2Gamma.M2()*cos_rotation;

  hists->fill2DHist(x_prime, y_prime, Form("h2D_dalitzPlot_rotation_%s_cut%i",  d, num),
		    ";M(l1#gamma)^{2};M(l2#gamma)^{2}", 100,0,22000, 100, -14000,5000,  weight, dir);

  hists->fill2DHist(tri.M(), diLep.M(),Form("h2D_tri_vs_diLep_mass0_%s_cut%i", d, num),
		    ";m_{#mu#mu#gamma} (GeV);m_{#mu#mu} (GeV)", 100,0,200, 100, 0,mllMax,  weight, dir);
  hists->fill2DHist(tri.M(), diLep.M(),Form("h2D_tri_vs_diLep_mass1_%s_cut%i", d, num),
		    ";m_{#mu#mu#gamma} (GeV);m_{#mu#mu} (GeV)", 100,0,200, 100, 0,1.5, weight, dir);
  hists->fill2DHist(tri.M(), diLep.M(),Form("h2D_tri_vs_diLep_mass2_%s_cut%i", d, num),
		    ";m_{#mu#mu#gamma} (GeV);m_{#mu#mu} (GeV)", 100,0,200, 100, 1.5,8, weight, dir);
  hists->fill2DHist(tri.M(), diLep.M(),Form("h2D_tri_vs_diLep_mass3_%s_cut%i", d, num),
		    ";m_{#mu#mu#gamma} (GeV);m_{#mu#mu} (GeV)", 100,0,200, 20, 8,mllMax,  weight, dir);
  hists->fill2DHist(tri.M(), diLep.M(),Form("h2D_tri_vs_diLep_Jpsi_%s_cut%i", d, num),
		    ";m_{#mu#mu#gamma} (GeV);m_{#mu#mu} (GeV)", 100,0,200, 100, 2.7,3.7,  weight, dir);
  hists->fill2DHist(tri.M(), diLep.M(),Form("h2D_tri_vs_diLep_apz_%s_cut%i", d, num),
		    ";m_{#mu#mu#gamma} (GeV);m_{#mu#mu} (GeV)", 100,0,200, 20, 15,mllMax,  weight, dir);

  hists->fill1DHist(diLep.M(),   Form("diLep_mass_low1_%s_cut%i",  d, num),";m_{#mu#mu} (GeV)", 20, 0,1.5,     weight, dir);
  hists->fill1DHist(diLep.M(),   Form("diLep_mass_low2_%s_cut%i",  d, num),";m_{#mu#mu} (GeV)", 20, 1.5,8,     weight, dir);
  hists->fill1DHist(diLep.M(),   Form("diLep_mass_low3_%s_cut%i",  d, num),";m_{#mu#mu} (GeV)", 20, 8,mllMax,  weight, dir);
  hists->fill1DHist(diLep.M(),   Form("diLep_mass_jpsi_%s_cut%i",  d, num),";m_{#mu#mu} (GeV)",100, 2.7,3.5,   weight, dir);
  hists->fill1DHist(diLep.M(),   Form("diLep_mass_apz_%s_cut%i",   d, num),";m_{#mu#mu} (GeV)", 20, 15,mllMax, weight, dir);
  hists->fill1DHist(diLep.M(),   Form("diLep_mass_low_%s_cut%i",   d, num),";m_{#mu#mu} (GeV)", 50, 0,mllMax,  weight, dir);
  hists->fill1DHist(diLep.M(),   Form("diLep_mass_high_%s_cut%i",  d, num),";m_{#mu#mu} (GeV)",100, 0,120,     weight, dir);

  //hists->fill1DHist(diLep.M(),   Form("diLep_mass_jpsi_%s_cut%i",d, dir, num),";m_{#mu#mu} (GeV)", 100, 2,5,   weight, "jpsi");
  hists->fill1DHist(diLep.Pt(),  Form("diLep_pt_%s_cut%i",   d, num),";di-Lepton p_{T}", 50, 0,120, weight, dir);
  hists->fill1DHist(diLep.Eta(), Form("diLep_eta_%s_cut%i",  d, num),";di-Lepton eta",   50, -3.5,3.5, weight, dir);
  hists->fill1DHist(diLep.Phi(), Form("diLep_phi_%s_cut%i",  d, num),";di-Lepton phi",   50, -TMath::Pi(),TMath::Pi(), weight, dir);
  hists->fill1DHist(diLepCM.E(), Form("diLep_Ecom_%s_cut%i", d, num),";E(ll) in CoM",    50, 0,200,  weight, dir);

  hists->fill1DHist(diLep.Pt()/tri.M(),  Form("diLep_ptOverMllg_%s_cut%i",   d, num),
		    ";di-Lepton p_{T}/m_{ll#gamma}",100, 0,1, weight, dir);
  hists->fill1DHist(gamma.Pt()/tri.M(), Form("gamma_ptOverMllg_%s_cut%i",   d, num),
		    ";Photon p_{T}/m_{ll#gamma}",  100, 0,1, weight, dir);
  hists->fill1DHist(gammaCM.E(),Form("gamma_Ecom_%s_cut%i", d, num),";E_{#gamma} in CoM",  50, 0,120,  weight, dir);
  hists->fill1DHist(gamma.E(),  Form("gamma_E_%s_cut%i",    d, num),";Photon Energy", 50, 0,200, weight, dir);
  hists->fill1DHist(gamma.Pt(), Form("gamma_pt_%s_cut%i",   d, num),";Photon p_{T} (GeV)",  50, 0,120, weight, dir);
  hists->fill1DHist(gamma.Eta(),Form("gamma_eta_%s_cut%i",  d, num),";Photon eta", 50, -3.5,3.5, weight, dir);
  hists->fill1DHist(gamma.Phi(),Form("gamma_phi_%s_cut%i",  d, num),";Photon phi", 50, -TMath::Pi(),TMath::Pi(), weight, dir);


  hists->fill1DHist(lPt1.Pt(),  Form("lPt1_pt_%s_cut%i", d, num), ";Leading lepton p_{T} (GeV)", 50, 0,100,    weight, dir);
  hists->fill1DHist(lPt1.Eta(), Form("lPt1_eta_%s_cut%i",d, num), ";Leading lepton eta",   50, -3.5,3.5, weight, dir);
  hists->fill1DHist(lPt1.Phi(), Form("lPt1_phi_%s_cut%i",d, num), ";Leading lepton phi",   50, -TMath::Pi(),TMath::Pi(), weight, dir);
  hists->fill1DHist(lPt2.Pt(),  Form("lPt2_pt_%s_cut%i", d, num), ";Trailing lepton p_{T} (GeV)",50, 0,100,    weight, dir);
  hists->fill1DHist(lPt2.Eta(), Form("lPt2_eta_%s_cut%i",d, num), ";Trailing lepton  eta", 50, -3.5,3.5, weight, dir);
  hists->fill1DHist(lPt2.Phi(), Form("lPt2_phi_%s_cut%i",d, num), ";Trailing lepton  phi", 50, -TMath::Pi(),TMath::Pi(), weight, dir);

  //hists->fill2DHist(l1.Pt(),   l2.Pt(),   Form("h2D_l1_vs_l2_%s_cut%i",  d, num),
  //";l+ pt; l- pt",    50, 0,100, 50,0,100, weight, dir);

  hists->fill2DHist(lPt1.Pt(), lPt2.Pt(), Form("h2D_Pt2_vs_Pt1_%s_cut%i", d, num),
		    ";Leading lepton p_{T};Trailing lepton p_{T}",  50, 0,100, 50,0,100, weight, dir);
  hists->fill2DHist(diLep.Pt(),gamma.Pt(),Form("h2D_diLep_vs_gamma_%s_cut%i",d, num),
		    ";q_{T}^{ll} (GeV); p_{T}^{#gamma} (GeV)",      50, 0,100, 50,0,100, weight, dir);

  hists->fill2DHist(diLep.Pt()/tri.M(),gamma.Pt()/tri.M(),Form("h2D_diLep_vs_gamma_ptOverMllg_%s_cut%i",d, num),
		    ";q_{T}^{ll}/m_{ll#gamma}; p_{T}^{#gamma}/m_{ll#gamma}",  100, 0,1, 100,0,1, weight, dir);

  hists->fill2DHist(gammaCM.E(), gamma.Pt(),Form("h2D_gamma_Ecom_vs_Pt_%s_cut%i",  d, num),
		    ";E_{#gamma} in CoM; p_{T}(#gamma)", 50, 0,100, 50,0,100, weight, dir);
  hists->fill2DHist(gammaCM.E(), tri.M(),   Form("h2D_gamma_Ecom_vs_triM_%s_cut%i",d, num),
		    ";E_{#gamma} in CoM; M(ll#gamma)",   50, 0,100, 50,0,200, weight, dir);

  hists->fill2DHist(diLep.Eta()-gamma.Eta(), gamma.Pt(), Form("h2D_deltaEta_vs_gammaPt_deltaEta_%s_cut%i",d, num),
		    ";#Delta#eta(ll, #gamma);p_{T} of {#gamma}",    100, -5,5, 100,0,130, weight, dir);

  hists->fill1DHist(lPt1.DeltaR(lPt2), Form("ll_deltaR_%s_cut%i",        d, num),";#Delta R(l_{1}, l_{2})",100,0,2, weight,dir);
  hists->fill1DHist(lPt1.DeltaR(gamma),Form("lPt1_gamma_deltaR_%s_cut%i",d, num),";#Delta R(#gamma,l_{1})",100,0,5, weight,dir);
  hists->fill1DHist(lPt2.DeltaR(gamma),Form("lPt2_gamma_deltaR_%s_cut%i",d, num),";#Delta R(#gamma,l_{2})",100,0,5, weight,dir);

  //ZGanlgles:
  /*
  if(selection!="el"){

    double co1,co2,phi,co3;
    ang->GetAngles(l1,l2,gamma,co1,co2,phi,co3);
    //cout<<eventNumber<<"RECO  Angles: c1= "<<co1<<"  c2="<<co2<<"   phi="<<phi<<"   coco="<<co3<<endl;

    hists->fill1DHist(co1, Form("co1_cut%i", num), ";cos_lp",  100,-1,1, 1,"Angles");
    hists->fill1DHist(co2, Form("co2_cut%i", num), ";cos_lm,", 100,-1,1, 1,"Angles");
    hists->fill1DHist(co3, Form("co3_cut%i", num), ";cosTheta",100,-1,1, 1,"Angles");
    hists->fill1DHist(phi, Form("phi_cut%i", num), ";phi lp",  100, -TMath::Pi(), TMath::Pi(), 1,"Angles");
  }
  */
  //hists->fill1DHist(nVtxTotal, Form("vtx_nPV_tot_%s_cut%i",   d, num), "vtx_nPV_tot",    40, 0,40, weight, dir);
  //hists->fill1DHist(nVtx,      Form("vtx_nPV_raw_%s_cut%i",   d, num), "Raw;nPV",    40, 0,40,      1, dir);
  //hists->fill1DHist(nVtx,      Form("vtx_nPV_weight_%s_cut%i",d, num), "Re-weighted;nPV", 40, 0,40, weight, dir);

  //hists->fill1DHist(nDofVtx1, Form("vtx_ndof1_%s_cut%i",  d, num), ";First vertex ndof", 50, 0,200, weight, dir);
  //hists->fill1DHist(nDofVtx2, Form("vtx_ndof2_%s_cut%i",  d, num), ";Second vertex ndof", 50, 0,200, weight, dir);

}

void MakeMuonPlots(HistManager *hists, TCMuon mu, TVector3 *pv)
{
  hists->fill1DHist(mu.PixelLayersWithMeasurement(),"mu_PixelLayersWithMeasurement",
		    ";PixelLayersWithMeasurement",   15, 0,15, 1, "Muons");
  hists->fill1DHist(mu.TrackLayersWithMeasurement(),"mu_TrackLayersWithMeasurement",
		    ";TrackerLayersWithMeasurement", 30, 0,30, 1, "Muons");
  hists->fill1DHist(mu.NumberOfMatchedStations(),   "mu_NumberOfMatchedStations",
		    ";NumberOfMatchedStations",  10, 0,10, 1, "Muons");
  hists->fill1DHist(mu.NumberOfValidMuonHits(),     "mu_NumberOfValidMuonHits",
		    ";NumberOfValidMuonHits",    60, 0,60, 1, "Muons");
  hists->fill1DHist(mu.NumberOfValidTrackerHits(),  "mu_NumberOfValidTrackerHits",
		    ";NumberOfValidTrackerHits", 40, 0,40, 1, "Muons");
  hists->fill1DHist(mu.NumberOfValidPixelHits(),    "mu_NumberOfValidPixelHits",
		    ";NumberOfValidPixelHits",   15, 0,15, 1, "Muons");
  hists->fill1DHist(mu.NormalizedChi2(),            "mu_NormalizedChi2",
		    ";NormalizedChi2",         100, 0,12,  1, "Muons");
  hists->fill1DHist(mu.NormalizedChi2_tracker(),    "mu_NormalizedChi2_tracker",
		    ";NormalizedChi2_tracker", 100, 0,6,   1, "Muons");
  hists->fill1DHist(mu.Dxy(pv), "mu_dxy", ";dxy", 50, -0.1,0.1, 1, "Muons");
  hists->fill1DHist(mu.Dz(pv),  "mu_dxz", ";dz",  50, -0.3,0.3, 1, "Muons");
  hists->fill1DHist(mu.PtError()/mu.Pt(), "mu_ptErrorOverPt", ";ptErrorOverPt", 50, 0,0.2, 1, "Muons");

  //Float_t muIso = 0;//CalculateMuonIso(&mu);
  //hists->fill1DHist(muIso, "mu_iso", ";rel isolation, pu corr", 50, 0,0.7, 1, "Muons");
  //hists->fill2DHist(muIso, global_Mll, "mu_iso_vs_Mll", ";pfIsolationR04;M(l_{1},l_{2})", 50, 0,0.7, 50, 0,20, 1, "Muons");
  //hists->fill1DHist(iso2, "mu_iso_official", ";pfIsolationR04", 50, 0,0.7, 1, "Muons");


}



void MakePhotonPlots(HistManager *hists, TCPhoton ph)
{

  hists->fill1DHist(ph.R9(),            "ph_R9",           ";R9",        100,    0, 1,    1,"Photon");
  hists->fill1DHist(ph.E2x5()/ph.R9(),  "ph_E2OverE9",     ";E2OverE9",  100,    0, 0.5,  1,"Photon");
  hists->fill1DHist(ph.TrackVeto(),     "ph_trackVeto",    ";track veto",  3,    0, 3,    1,"Photon");
  hists->fill1DHist(ph.HadOverEm(),     "ph_HadOverEm",    ";HadOverEm", 100,    0, 0.05, 1,"Photon");
  //hists->fill1DHist(ph.NormChi2(),    "ph_NormChi2",     ";NormChi2",  100,    0, 0.01, 1,"Photon");
  hists->fill1DHist(ph.SCDeltaPhi(),    "ph_SCDeltaPhi",   ";SCDeltaPhi", 100, -0.2, 0.2, 1,"Photon");
  hists->fill1DHist(ph.SCDeltaEta(),    "ph_SCDeltaEta",   ";SCDeltaEta", 100,-0.02, 0.02,1,"Photon");
  hists->fill1DHist(ph.GetNCrystals(),  "ph_nCrystals",    ";nCrystals",     160, 0, 160, 1,"Photon");
  hists->fill1DHist(ph.SigmaIEtaIEta(), "ph_SigmaIEtaIEta",";SigmaIEtaIEta", 100, 0, 0.1, 1,"Photon");
  hists->fill1DHist(ph.SigmaIPhiIPhi(), "ph_SigmaIPhiIPhi",";SigmaIPhiIPhi", 100, 0, 0.1, 1,"Photon");
  hists->fill1DHist(ph.ConversionVeto(),"ph_ConversionVeto",";ConversionVeto", 3, 0, 3,   1,"Photon");
  //hists->fill1DHist(ph., "ph_",";", 3, 0, 3, 1, "Photon");
  //hists->fill1DHist(ph., "ph_",";", 3, 0, 3, 1, "Photon");

}

#endif
