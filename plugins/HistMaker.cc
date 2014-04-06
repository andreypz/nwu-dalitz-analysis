#include "HistMaker.h"

HistMaker::HistMaker(HistManager *h)
{
  hists = h;
}

HistMaker::~HistMaker(){}


void HistMaker::MakeMuonPlots(TCMuon mu, TVector3 *pv)
{
  hists->fill1DHist(mu.PixelLayersWithMeasurement(),"mu_PixelLayersWithMeasurement",";PixelLayersWithMeasurement",   15, 0,15, 1, "Muons");
  hists->fill1DHist(mu.TrackLayersWithMeasurement(),"mu_TrackLayersWithMeasurement",";TrackerLayersWithMeasurement", 30, 0,30, 1, "Muons");
  hists->fill1DHist(mu.NumberOfMatchedStations(),   "mu_NumberOfMatchedStations",   ";NumberOfMatchedStations",  10, 0,10, 1, "Muons");
  hists->fill1DHist(mu.NumberOfValidMuonHits(),     "mu_NumberOfValidMuonHits",     ";NumberOfValidMuonHits",    60, 0,60, 1, "Muons");
  hists->fill1DHist(mu.NumberOfValidTrackerHits(),  "mu_NumberOfValidTrackerHits",  ";NumberOfValidTrackerHits", 40, 0,40, 1, "Muons");
  hists->fill1DHist(mu.NumberOfValidPixelHits(),    "mu_NumberOfValidPixelHits",    ";NumberOfValidPixelHits",   15, 0,15, 1, "Muons");
  hists->fill1DHist(mu.NormalizedChi2(),            "mu_NormalizedChi2",            ";NormalizedChi2",         100, 0,12,  1, "Muons");
  hists->fill1DHist(mu.NormalizedChi2_tracker(),    "mu_NormalizedChi2_tracker",    ";NormalizedChi2_tracker", 100, 0,6,   1, "Muons");
  hists->fill1DHist(mu.Dxy(pv), "mu_dxy", ";dxy", 50, -0.1,0.1, 1, "Muons");
  hists->fill1DHist(mu.Dz(pv),  "mu_dxz", ";dz",  50, -0.3,0.3, 1, "Muons");
  hists->fill1DHist(mu.PtError()/mu.Pt(), "mu_ptErrorOverPt", ";ptErrorOverPt", 50, 0,0.2, 1, "Muons");

  //Float_t muIso = ObjID->CalculateMuonIso(&mu);
  //hists->fill1DHist(muIso, "mu_iso", ";rel isolation, pu corr", 50, 0,0.7, 1, "Muons");
  //hists->fill2DHist(muIso, global_Mll, "mu_iso_vs_Mll", ";pfIsolationR04;M(l_{1},l_{2})", 50, 0,0.7, 50, 0,20, 1, "Muons");
  //Float_t iso2 = (mu.IsoMap("pfChargedHadronPt_R04") + mu.IsoMap("pfNeutralHadronEt_R04") + mu.IsoMap("pfPhotonEt_R04"))/mu.Pt();
  //hists->fill1DHist(iso2, "mu_iso_official", ";pfIsolationR04", 50, 0,0.7, 1, "Muons");


}

void HistMaker::MakePhotonPlots(TCPhoton ph)
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

void HistMaker::FillHistosFull(Int_t num, Double_t weight,
			    TCPhysObject l1, TCPhysObject l2, TCPhysObject lPt1, TCPhysObject lPt2, TCPhysObject gamma, string dir)
{

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

  hists->fill1DHist(tri.M(),     Form("tri_mass_%s_cut%i",         d, num),";M(ll#gamma)",  100, 0,200,  weight, dir);
  hists->fill1DHist(tri.M(),     Form("tri_mass80_%s_cut%i",       d, num),";M(ll#gamma)",  100,80,200,  weight, dir);
  hists->fill1DHist(tri.M(),     Form("tri_mass_longTail_%s_cut%i",d, num),";M(ll#gamma)",  100, 0,400,  weight, dir);

  hists->fill1DHist(lPt1Gamma.M(), Form("lPt1Gamma_%s_cut%i",d, num),";M(l1#gamma)",  100, 0,400,  weight, dir);
  hists->fill1DHist(lPt2Gamma.M(), Form("lPt2Gamma_%s_cut%i",d, num),";M(l2#gamma)",  100, 0,400,  weight, dir);

  hists->fill2DHist(lPt1Gamma.M2(), lPt2Gamma.M2(),
		    Form("h2D_dalitzPlot_%s_cut%i",  d, num),";M(l1#gamma)^{2};M(l2#gamma)^{2}", 100,0,20000, 100, 0,10000,  weight, dir);

  double rotation = atan(-1.);
  double sin_rotation = sin(rotation);
  double cos_rotation = cos(rotation);
  double x_prime = (lPt1Gamma.M2()*cos_rotation - lPt2Gamma.M2()*sin_rotation) + (1.-cos_rotation)*pow(125.,2);
  double y_prime = lPt1Gamma.M2()*sin_rotation + lPt2Gamma.M2()*cos_rotation;

  hists->fill2DHist(x_prime, y_prime,
		    Form("h2D_dalitzPlot_rotation_%s_cut%i",  d, num),";M(l1#gamma)^{2};M(l2#gamma)^{2}", 100,0,22000, 100, -14000,5000,  weight, dir);

  hists->fill2DHist(tri.M(), diLep.M(),
		    Form("h2D_tri_vs_diLep_mass0_%s_cut%i",  d, num),";M(llgamma);M(ll)", 100,0,200, 100, 0,20,  weight, dir);
  hists->fill2DHist(tri.M(), diLep.M(),
		    Form("h2D_tri_vs_diLep_mass1_%s_cut%i",  d, num),";M(llgamma);M(ll)", 100,0,200, 100, 0,1.5, weight, dir);
  hists->fill2DHist(tri.M(), diLep.M(),
		    Form("h2D_tri_vs_diLep_mass2_%s_cut%i",  d, num),";M(llgamma);M(ll)", 100,0,200, 100, 1.5,8, weight, dir);
  hists->fill2DHist(tri.M(), diLep.M(),
		    Form("h2D_tri_vs_diLep_mass3_%s_cut%i",  d, num),";M(llgamma);M(ll)", 100,0,200, 100, 8,20,  weight, dir);

  hists->fill1DHist(diLep.M(),   Form("diLep_mass_low1_%s_cut%i",  d, num),";M(ll)", 100, 0,1.5,  weight, dir);
  hists->fill1DHist(diLep.M(),   Form("diLep_mass_low2_%s_cut%i",  d, num),";M(ll)", 100, 1.5,8,  weight, dir);
  hists->fill1DHist(diLep.M(),   Form("diLep_mass_low3_%s_cut%i",  d, num),";M(ll)", 100, 8,20,   weight, dir);
  hists->fill1DHist(diLep.M(),   Form("diLep_mass_jpsi_%s_cut%i",  d, num),";M(ll)", 100, 2.7,3.5,weight, dir);

  hists->fill1DHist(diLep.M(),   Form("diLep_mass_low_%s_cut%i",  d, num),";M_{#mu#mu} (GeV)", 100, 0,20,  weight, dir);
  hists->fill1DHist(diLep.M(),   Form("diLep_mass_high_%s_cut%i", d, num),";M_{#mu#mu} (GeV)", 100, 0,120, weight, dir);
  //hists->fill1DHist(diLep.M(),   Form("diLep_mass_jpsi_%s_cut%i",d, dir, num),";M(ll)", 100, 2,5,   weight, "jpsi");
  hists->fill1DHist(diLep.Pt(),  Form("diLep_pt_%s_cut%i",   d, num),";di-Lepton p_{T}", 50, 0,120, weight, dir);
  hists->fill1DHist(diLep.Eta(), Form("diLep_eta_%s_cut%i",  d, num),";di-Lepton eta",   50, -3.5,3.5, weight, dir);
  hists->fill1DHist(diLep.Phi(), Form("diLep_phi_%s_cut%i",  d, num),";di-Lepton phi",   50, -TMath::Pi(),TMath::Pi(), weight, dir);
  hists->fill1DHist(diLepCM.E(), Form("diLep_Ecom_%s_cut%i", d, num),";E(ll) in CoM",    50, 0,200,  weight, dir);

  hists->fill1DHist(diLep.Pt()/tri.M(),  Form("diLep_ptOverMllg_%s_cut%i",   d, num),
		    ";di-Lepton p_{T}/M_{ll#gamma}",100, 0,1, weight, dir);
  hists->fill1DHist(gamma.Pt()/tri.M(), Form("gamma_ptOverMllg_%s_cut%i",   d, num),
		    ";Photon p_{T}/M_{ll#gamma}",  100, 0,1, weight, dir);
  hists->fill1DHist(gammaCM.E(),Form("gamma_Ecom_%s_cut%i", d, num),";E_{#gamma} in CoM",  50, 0,120,  weight, dir);
  hists->fill1DHist(gamma.E(),  Form("gamma_E_%s_cut%i",    d, num),";Photon Energy", 50, 0,200, weight, dir);
  hists->fill1DHist(gamma.Pt(), Form("gamma_pt_%s_cut%i",   d, num),";Photon p{T} (GeV)",  50, 0,120, weight, dir);
  hists->fill1DHist(gamma.Eta(),Form("gamma_eta_%s_cut%i",  d, num),";Photon eta", 50, -3.5,3.5, weight, dir);
  hists->fill1DHist(gamma.Phi(),Form("gamma_phi_%s_cut%i",  d, num),";Photon phi", 50, -TMath::Pi(),TMath::Pi(), weight, dir);
  /*
  hists->fill1DHist(l1.Pt(),    Form("l1_pt_%s_cut%i",  d, num), ";l+ pt",    50, 0,100, weight, dir);
  hists->fill1DHist(l1.Eta(),   Form("l1_eta_%s_cut%i", d, num), ";l+ eta",   50, -3.5,3.5, weight, dir);
  hists->fill1DHist(l1.Phi(),   Form("l1_phi_%s_cut%i", d, num), ";l+ phi",   50, -TMath::Pi(),TMath::Pi(), weight, dir);
  hists->fill1DHist(l2.Pt(),    Form("l2_pt_%s_cut%i",  d, num), ";l- pt",    50, 0,100, weight, dir);
  hists->fill1DHist(l2.Eta(),   Form("l2_eta_%s_cut%i", d, num), ";l- eta",   50, -3.5,3.5, weight, dir);
  hists->fill1DHist(l2.Phi(),   Form("l2_phi_%s_cut%i", d, num), ";l- phi",   50, -TMath::Pi(),TMath::Pi(), weight, dir);
  */
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
		    ";q_{T}^{ll}/M_{ll#gamma}; p_{T}^{#gamma}/M_{ll#gamma}",  100, 0,1, 100,0,1, weight, dir);

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


void HistMaker::MakePhotonEnergyCorrPlots(TCPhoton pho, Float_t corrPhoEnReco, Float_t corrPhoEnSC)
{
  hists->fill1DHist(pho.M(), "ph_recoMass",";M_{#gamma}", 100,-0.05,0.05, 1,"Photon");

  hists->fill1DHist((pho.E() - corrPhoEnReco)/pho.E(),"ph_energyCorrectionReco",
		    ";(E_{#gamma}^{reco} - E_{#gamma}^{reco corr})/E_{#gamma}^{reco}", 100,-0.1,0.1, 1,"Photon");

  hists->fill1DHist((pho.SCEnergy() - corrPhoEnSC)/pho.SCEnergy(),"ph_energyCorrectionSC",
		    ";(E_{#gamma}^{SC} - E_{#gamma}^{SC corr})/E_{#gamma}^{SC}", 100,-0.1,0.1, 1,"Photon");

  hists->fill1DHist((pho.E() - pho.SCEnergy())/pho.E(),"ph_energyRecoVsSC",
		    ";(E_{#gamma}^{reco} - E_{#gamma}^{SC})/E_{#gamma}^{reco}", 100,-0.1,0.1, 1,"Photon");

  if (pho.Pt()>40 && fabs(pho.SCEta())<1.444){
    hists->fill1DHist((pho.E() - pho.SCEnergy())/pho.E(),"ph_energyRecoVsSC_PtEtaCut",
		      ";(E_{#gamma}^{reco} - E_{#gamma}^{SC})/E_{#gamma}^{reco}", 100,-0.1,0.1, 1,"Photon");

    hists->fill1DHist((pho.E() - corrPhoEnReco)/pho.E(),"ph_energyCorrectionReco_PtEtaCut",
		      ";(E_{#gamma}^{reco} - E_{#gamma}^{reco corr})/E_{#gamma}^{reco}", 100,-0.1,0.1, 1,"Photon");

    hists->fill1DHist((pho.SCEnergy() - corrPhoEnSC)/pho.SCEnergy(),"ph_energyCorrectionSC_PtEtaCut",
		      ";(E_{#gamma}^{SC} - E_{#gamma}^{SC corr})/E_{#gamma}^{SC}", 100,-0.1,0.1, 1,"Photon");

    if (pho.R9() > 0.94)
      {
	hists->fill1DHist((pho.E() - corrPhoEnReco)/pho.E(),"ph_energyCorrectionReco_PtEtaCut_hightR9",
			  ";(E_{#gamma}^{reco} - E_{#gamma}^{reco corr})/E_{#gamma}^{reco}", 100,-0.1,0.1, 1,"Photon");

	hists->fill1DHist((pho.SCEnergy()-corrPhoEnSC)/pho.SCEnergy(),"ph_energyCorrectionSC_PtEtaCut_highR9",
			  ";(E_{#gamma}^{SC} - E_{#gamma}^{SC corr})/E_{#gamma}^{SC}", 100,-0.1,0.1, 1,"Photon");
      }
    else
      {
	hists->fill1DHist((pho.E() - corrPhoEnReco)/pho.E(),"ph_energyCorrectionReco_PtEtaCut_lowR9",
			  ";(E_{#gamma}^{reco} - E_{#gamma}^{reco corr})/E_{#gamma}^{reco}", 100,-0.1,0.1, 1,"Photon");

	hists->fill1DHist((pho.SCEnergy()-corrPhoEnSC)/pho.SCEnergy(),"ph_energyCorrectionSC_PtEtaCut_lowR9",
			  ";(E_{#gamma}^{SC} - E_{#gamma}^{SC corr})/E_{#gamma}^{SC}", 100,-0.1,0.1, 1,"Photon");

      }
  }
}

void HistMaker::MakeZeePlots(TCPhoton p1, TCPhoton p2)
{
  hists->fill2DHist(p1.R9(), p2.R9(), "zee_p1R9_p2R9",";#gamma_{1} R9; #gamma_{2} R9", 100, 0,1, 100,0,1,  1, "Zee");
  Float_t mZ = (p1+p2).M();
  hists->fill1DHist(mZ, "zee_M",";", 100, 70, 110, 1, "Zee");
  if (p1.R9() > 0.94 && p2.R9() > 0.94)
    hists->fill1DHist(mZ, "zee_M_highR9",";", 100, 70, 110, 1, "Zee");
  else
    hists->fill1DHist(mZ, "zee_M_lowR9",";", 100, 70, 110, 1, "Zee");

  //hists->fill1DHist(mZ, "zee_M",";", 200, 60, 130, 1, "Zee");
  //hists->fill1DHist(mZ, "zee_M",";", 200, 60, 130, 1, "Zee");

}
