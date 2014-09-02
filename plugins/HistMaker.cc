#include "HistMaker.h"

HistMaker::HistMaker(HistManager *h):
  _isVtxSet(false),
  _isLepSet(false),
  _isGammaSet(false),
  _isGamma2Set(false),
  _rhoFactor(-1),
  _nVtx(0)
{
  hists = h;
  angles = new ZGAngles();
}

HistMaker::~HistMaker(){}

void HistMaker::Reset(float rho, UInt_t n){
  _isLepSet = 0; _isGammaSet=0; _isGamma2Set=0;
  _rhoFactor = rho; _nVtx=n;}

void HistMaker::SetVtx(TCPrimaryVtx v){
  _pv=v; _isVtxSet = 1;}
void HistMaker::SetLeptons(TCPhysObject l1, TCPhysObject l2){
  _lPt1 = l1;  _lPt2 = l2; _isLepSet = 1;}
void HistMaker::SetGamma(TCPhysObject g){
  _gamma = g; _isGammaSet=1;}
void HistMaker::SetGamma2(TCPhysObject g){
  _gamma2 = g; _isGamma2Set=1;}

void HistMaker::MakeMuonPlots(const TCMuon& mu, string dir)
{

  hists->fill1DHist(mu.PixelLayersWithMeasurement(),dir+"_mu_PixelLayersWithMeasurement",";PixelLayersWithMeasurement",   15, 0,15, 1, dir);
  hists->fill1DHist(mu.TrackLayersWithMeasurement(),dir+"_mu_TrackLayersWithMeasurement",";TrackerLayersWithMeasurement", 30, 0,30, 1, dir);
  hists->fill1DHist(mu.NumberOfMatchedStations(), dir+"_mu_NumberOfMatchedStations", ";NumberOfMatchedStations", 10, 0,10, 1, dir);
  hists->fill1DHist(mu.NumberOfValidMuonHits(),   dir+"_mu_NumberOfValidMuonHits",   ";NumberOfValidMuonHits",   60, 0,60, 1, dir);
  hists->fill1DHist(mu.NumberOfValidTrackerHits(),dir+"_mu_NumberOfValidTrackerHits",";NumberOfValidTrackerHits",40, 0,40, 1, dir);
  hists->fill1DHist(mu.NumberOfValidPixelHits(),  dir+"_mu_NumberOfValidPixelHits",  ";NumberOfValidPixelHits",  15, 0,15, 1, dir);
  hists->fill1DHist(mu.NormalizedChi2_tracker(),  dir+"_mu_NormalizedChi2_tracker",  ";NormalizedChi2_tracker", 100, 0,4,  1, dir);
  hists->fill1DHist(mu.NormalizedChi2(),          dir+"_mu_NormalizedChi2",          ";NormalizedChi2",         100, 0,4,  1, dir);
  hists->fill1DHist(mu.Dxy(&_pv), dir+"_mu_dxy", ";dxy", 50, -0.02,0.02, 1, dir);
  hists->fill1DHist(mu.Dz(&_pv),  dir+"_mu_dxz", ";dz",  50, -0.1,0.1, 1, dir);
  hists->fill1DHist(mu.PtError()/mu.Pt(), dir+"_mu_ptErrorOverPt", ";ptErrorOverPt", 50, 0,0.1, 1, dir);

  //Float_t muIso = ObjID->CalculateMuonIso(&mu);
  //hists->fill1DHist(muIso, dir+"_mu_iso", ";rel isolation, pu corr", 50, 0,0.7, 1, dir);
  //hists->fill2DHist(muIso, global_Mll, dir+"_mu_iso_vs_Mll", ";pfIsolationR04;M(l_{1},l_{2})", 50, 0,0.7, 50, 0,20, 1, dir);
  //Float_t iso2 = (mu.IsoMap("pfChargedHadronPt_R04") + mu.IsoMap("pfNeutralHadronEt_R04") + mu.IsoMap("pfPhotonEt_R04"))/mu.Pt();
  //hists->fill1DHist(iso2, dir+"_mu_iso_official", ";pfIsolationR04", 50, 0,0.7, 1, dir);

}

void HistMaker::MakeEGammaCommonPlots(const TCEGamma& egm, TString n)
{
  string dir = n.Data();

  hists->fill1DHist(egm.SCEta(),           dir+"-egm_SCEta",        ";SCEta",        100, -2.5, 2.5,  1,dir);
  hists->fill1DHist(egm.R9(),              dir+"-egm_R9",           ";R9",           100,    0.2, 1,  1,dir);
  hists->fill1DHist(egm.E1x3()/egm.E5x5(), dir+"-egm_E1x3OverE5x5", ";E1x3OverE5x5", 100,    0.1, 1,  1,dir);
  hists->fill1DHist(egm.E2x2()/egm.E5x5(), dir+"-egm_E2x2OverE5x5", ";E2x2OverE5x5", 100,    0.1, 1,  1,dir);
  hists->fill1DHist(egm.E1x5()/egm.E5x5(), dir+"-egm_E1x5OverE5x5", ";E1x5OverE5x5", 100,    0.1, 1,  1,dir);
  hists->fill1DHist(egm.E2x5()/egm.E5x5(), dir+"-egm_E2x5OverE5x5", ";E2x5OverE5x5", 100,    0.5, 1,  1,dir);
  //hists->fill1DHist(egm.E2x5Max()/egm.R9(),dir+"-egm_E2x5MaxOverR9",";E2x5MaxOverR9",100,    0, 1,  1,dir);

  hists->fill1DHist(egm.SigmaIEtaIPhi(),   dir+"-egm_SigmaIEtaIPhi",";SigmaIEtaIPhi",100,-2e-4,2e-4,  1,dir);
  hists->fill1DHist(egm.SigmaIEtaIEta(),   dir+"-egm_SigmaIEtaIEta",";SigmaIEtaIEta",100, 0, 0.016,   1,dir);
  hists->fill1DHist(egm.SigmaIPhiIPhi(),   dir+"-egm_SigmaIPhiIPhi",";SigmaIPhiIPhi",100, 0, 0.04,    1,dir);
  hists->fill1DHist(egm.SCEnergy(),        dir+"-egm_SCEnergy",     ";SCEnergy",     100, 10, 120,    1,dir);
  //  hists->fill1DHist(egm.GetNCrystals(),    dir+"-egm_nCrystals",    ";nCrystals",    160, 0, 160,     1,dir);

  if (!n.Contains("BaseSC")) {
    hists->fill1DHist(egm.HadOverEm(),       dir+"-egm_HadOverEm",    ";HadOverEm",    100,0.001, 0.05, 1,dir);
    hists->fill1DHist(egm.SCEtaWidth(),      dir+"-egm_SCEtaWidth",   ";SCEtaWidth",   100, 0,  0.03,   1,dir);
    hists->fill1DHist(egm.SCPhiWidth(),      dir+"-egm_SCPhiWidth",   ";SCPhiWidth",   100, 0,  0.12,   1,dir);
    hists->fill1DHist(egm.SCRawEnergy(),     dir+"-egm_SCRawEnergy",  ";SCRawEnergy",  100, 10, 120,    1,dir);
    if (fabs(egm.SCEta()) > 1.56)
      hists->fill1DHist(egm.PreShowerOverRaw(),dir+"-egm_PreShowerOverRaw",";PreShowerOverRaw",100,0,0.2, 1,dir);
  }

  if (n.Contains("Ele") && !n.Contains("BaseSC")){ //Not filled for photons
    hists->fill1DHist(egm.SCDeltaPhi(),   dir+"-egm_SCDeltaPhi",   ";SCDeltaPhiAtVtx", 100, -0.1, 0.1, 1, dir);
    hists->fill1DHist(egm.SCDeltaEta(),   dir+"-egm_SCDeltaEta",   ";SCDeltaEtaAtVtx", 100,-0.01,0.01, 1, dir);
  }

  //if (n.Contains("BaseSC")){
  // egm.Dump();
  //}

}


void HistMaker::MakePhotonPlots(const TCPhoton& ph, string dir)
{
  //cout<<"DBG make photon plots  "<<dir<<endl;

  MakeEGammaCommonPlots(ph, dir+"-EGamma");
  hists->fill1DHist(ph.TrackVeto(),      dir+"ph_trackVeto",     ";track veto",      3, 0, 3,   1,dir);
  hists->fill1DHist(ph.ConversionVeto(), dir+"ph_ConversionVeto",";ConversionVeto",  3, 0, 3,   1,dir);
  hists->fill1DHist(ph.PfIsoPhoton(),    dir+"ph_PfIsoPhoton",   ";PfIsoPhoton",   100, 0.01,5, 1,dir);
  hists->fill1DHist(ph.PfIsoCharged(),   dir+"ph_PfIsoCharged",  ";PfIsoCharged",  100, 0.01,2, 1,dir);
  if (fabs(ph.SCEta())>1.6){ //These are only for EE
    hists->fill1DHist(ph.ESEffSigmaRR()[0],dir+"ph_ESEffSigmaRR_x",";ESEffSigmaRR_x",100, -5e-8,5e-8, 1,dir);
    hists->fill1DHist(ph.ESEffSigmaRR()[1],dir+"ph_ESEffSigmaRR_y",";ESEffSigmaRR_y",100, -5e-8,5e-8, 1,dir);
  }
  hists->fill1DHist(ph.IdMap("mvaScore"),    dir+"ph_mvaScore",      ";MVA score",     100, -1,1,   1,dir);
  hists->fill1DHist(ph.CiCPF4chgpfIso02()[0],dir+"ph_CiCPF4chgpfIso02", ";ph_CiCPF4chgpfIso02", 100, 0,5, 1,dir);

  if (ph.R9() > 0.9){
    hists->fill1DHist(ph.IdMap("HadIso_R03")-0.005*ph.Et(), dir+"ph_HadIso_R03_R9high", ";ph_HadIso_R03_R9high", 100, 0, 60, 1,dir);
    hists->fill1DHist(ph.IdMap("TrkIso_R03")-0.002*ph.Et(), dir+"ph_TrkIso_R03_R9high", ";ph_HadIso_R03_R9high", 100, 0, 60, 1,dir);
  }
  else{
    hists->fill1DHist(ph.IdMap("HadIso_R03")-0.005*ph.Et(), dir+"ph_HadIso_R03_R9low", ";ph_HadIso_R03_R9low", 100, 0, 6, 1,dir);
    hists->fill1DHist(ph.IdMap("TrkIso_R03")-0.002*ph.Et(), dir+"ph_TrkIso_R03_R9low", ";ph_HadIso_R03_R9low", 100, 0, 6, 1,dir);
  }

  //hists->fill1DHist(ph.(), dir+"ph_", ";", 100, 0, 100,1, dir);
}


void HistMaker::MakeElectronPlots(const TCElectron& el, string dir)
{

  MakeEGammaCommonPlots(el, dir+"-EGamma");

  hists->fill1DHist(el.IdMap("mvaScore"),dir+"el_mvaScore",";Dalitz MVA score",100, -1,1, 1,dir);
  hists->fill1DHist(el.MvaID(), dir+"_el_mvaID",";HZZ mva ID", 100, -1,1, 1,dir);
  hists->fill1DHist(el.FBrem(), dir+"_el_fbrem",";fbrem",      100, 0, 1, 1,dir);

  hists->fill1DHist(el.InverseEnergyMomentumDiff(), dir+"_el_fabsEPDiff",";|1/E - 1/p|", 100, 0, 0.1, 1, dir);
  hists->fill1DHist(el.PtError()/el.Pt(),     dir+"_el_ptErrorOverPt",";ptErrorOverPt", 100, 0, 1, 1,dir);
  hists->fill1DHist(el.NormalizedChi2(),      dir+"_el_gsfChi2", ";gsfChi2", 100, 0, 3, 1,dir);
  hists->fill1DHist(el.NormalizedChi2Kf(),    dir+"_el_kfChi2",  ";kfChi2",  100, 0, 3, 1,dir);
  hists->fill1DHist(el.DeltaEtaSeedCluster(), dir+"_el_deltaEtaSeedAtCala", ";DeltaEtaSeedAtCalo", 100, -0.05, 0.05, 1,dir);
  hists->fill1DHist(el.DeltaPhiSeedCluster(), dir+"_el_deltaPhiSeedAtCalo", ";DeltaPhiSeedAtCalo", 100, -0.05, 0.05, 1,dir);
  hists->fill1DHist(1 - el.E1x5()/el.E5x5(),  dir+"_el_ome1x5oe5x5",";1 - e1x5/e5x5",100, 0, 1,      1,dir);
  hists->fill1DHist(el.EoP(),    dir+"_el_EoP",    ";EoP",    100, 0, 6,   1,dir);
  hists->fill1DHist(el.EoPout(), dir+"_el_EoPout", ";EoPout", 100, 0, 6,   1,dir);
  hists->fill1DHist(el.IP3d(),   dir+"_el_ip3d",   ";ip3d",   100, 0, 0.04,1,dir);
  hists->fill1DHist(el.IP3dSig(),dir+"_el_ip3dSig",";ip3dSig",100, 0, 3,   1,dir);
  hists->fill1DHist(el.NumberOfValidHits(), dir+"_el_NumberOfValidHits", ";NumberOfValidHits", 25, 0,25, 1,dir);
  hists->fill1DHist(el.TrackerLayersWithMeasurement(), dir+"_el_TrackerLayersWithMeasurement",
		    ";TrackerLayersWithMeasurement", 20, 0,20,1,dir);

  hists->fill1DHist(el.PassConversionVeto(), dir+"_el_conv_passConv",";pass conv veto",  3, 0,  3, 1, dir);
  hists->fill1DHist(el.ConversionMissHits(), dir+"_el_conv_MissHits",";conv miss hits", 10, 0, 10, 1, dir);

  hists->fill1DHist(el.ConversionDist(),   dir+"_el_convDist",  ";conv Dist",   100, 0, 1,  1, dir);
  hists->fill1DHist(el.ConversionDcot(),   dir+"_el_convDcot",  ";conv Dcot",   100, -1, 1, 1, dir);
  hists->fill1DHist(el.ConversionRadius(), dir+"_el_convRadius",";conv Radius", 100, 0, 50, 1, dir);

  hists->fill1DHist(el.GetTracks().size(), dir+"_el_nTracks", ";nTracks", 10, 0,10, 1,dir);

  if (el.GetTracks().size()>=2){
    Float_t mll = (el.GetTracks()[0]+el.GetTracks()[1]).M();
    Float_t dR  = el.GetTracks()[0].DeltaR(el.GetTracks()[1]);
    hists->fill1DHist(dR,  dir+"_el_dR_all", ";dR(gsf1, gsf2)", 100,   0,4, 1,dir);
    hists->fill1DHist(dR,  dir+"_el_dR_low", ";dR(gsf1, gsf2)", 100,   0,1, 1,dir);
    hists->fill1DHist(mll, dir+"_el_mll_full", ";m_{ll} (GeV)", 100, 0,120, 1,dir);
    hists->fill1DHist(mll, dir+"_el_mll_low",  ";m_{ll} (GeV)", 100,  0,20, 1,dir);
    hists->fill1DHist(mll, dir+"_el_mll_jpsi", ";m_{ll} (GeV)", 100,  1, 5, 1,dir);
  }

  for (UInt_t i=0; i<el.GetTracks().size() && i<3; i++){
    string subdir = dir+Form("-track-%i",i+1);
    hists->fill1DHist(el.GetTracks()[i].Pt(), subdir+Form("_el_01_track_%i",i+1)+"_Pt",
		      Form(";gsf track %i Pt (GeV)",i+1), 200, 0,100, 1,subdir);
    hists->fill1DHist(el.GetTracks()[i].PtError()/el.GetTracks()[i].Pt(), subdir+Form("_el_01_track_%i",i+1)+"_PtErrorOverPt",
		      Form(";Pt Error /Pt for track %i",i+1), 200, 0,1.5, 1,subdir);
    hists->fill1DHist(el.GetTracks()[i].NormalizedChi2(), subdir+Form("_el_track_%i",i+1)+"_01_NormalizedChi2",
		      Form(";NormalizedChi2 for track %i",i+1), 100, 0,3, 1,subdir);

    TCTrack::ConversionInfo inf = el.GetTracks()[i].GetConversionInfo();
    hists->fill1DHist(inf.isValid, subdir+Form("_el_conv_%i",i+1)+"_isValid",Form(";isValid trk %i", i+i),3, 0,3, 1,subdir);
    if (inf.isValid){
      hists->fill1DHist(inf.nHitsMax,subdir+Form("_el_conv_%i",i+1)+"_nHitsMax",Form(";nHitsMax trk %i",i+i),  20, 0,20, 1,subdir);
      hists->fill1DHist(inf.vtxProb, subdir+Form("_el_conv_%i",i+1)+"_vtxProb", Form(";vtxProb trk %i", i+i), 200,0,0.01,1,subdir);
      //hists->fill1DHist(inf.lxyPV,   subdir+Form("_el_conv_%i",i+1)+"_lxyPV",   Form(";lxyPV trk %i",   i+i), 200, 0,50, 1,subdir);
      hists->fill1DHist(inf.lxyBS,   subdir+Form("_el_conv_%i",i+1)+"_lxyBS",   Form(";lxyBS trk %i",   i+i), 200, 0,50, 1,subdir);
    }
  }

  hists->fill1DHist(el.BaseSC().size(), dir+"_el_nBaseSC", ";nBaseSC", 10, 0,10, 1,dir);

  for (UInt_t i=0; i<el.BaseSC().size() && i<2; i++)
    MakeEGammaCommonPlots(el.BaseSC()[i], dir+Form("-BaseSC-%i",i+1));
}

void HistMaker::FillHistosFull(Int_t num, Double_t weight, string dir)
{
  //cout<<num<<" Making plots "<<dir<<endl;
  char * d = new char [dir.length()+1];
  std::strcpy (d, dir.c_str());

  TLorentzVector diLep = _lPt1 + _lPt2;
  TLorentzVector tri   = diLep + _gamma;
  TLorentzVector four  = diLep + _gamma + _gamma2;

  TLorentzVector lPt1Gamma = _lPt1+_gamma;
  TLorentzVector lPt2Gamma = _lPt2+_gamma;

  TLorentzVector gammaCM = _gamma;
  TLorentzVector diLepCM = diLep;
  TVector3 b1  = -tri.BoostVector();
  if (_isGammaSet)
    gammaCM.Boost(b1);
  if (_isLepSet)
    diLepCM.Boost(b1);

  hists->fill1DHist(weight,  Form("weight_%s_cut%i",         d, num),";weight",  200, 0,3, 1, dir);

  if (_isGammaSet && _isGamma2Set)
    hists->fill1DHist(four.M(), Form("four_massZ_%s_cut%i",d, num),";m_{ll#gamma#gamma}",  100,70,120,  weight, dir);

  hists->fill1DHist(tri.M(), Form("00_tri_mass_%s_cut%i",         d, num),";m_{ll#gamma}",  100, 0,200,  weight, dir);
  hists->fill1DHist(tri.M(), Form("00_tri_mass_3GeVBin_%s_cut%i", d, num),";m_{ll#gamma}",  20,110,170,  weight, dir);

  hists->fill1DHist(tri.M(), Form("00_tri_mass125_%s_cut%i",      d, num),";m_{ll#gamma}",  30,110,140,  weight, dir);
  hists->fill1DHist(tri.M(), Form("00_tri_mass80_%s_cut%i",       d, num),";m_{ll#gamma}",  100,80,200,  weight, dir);
  hists->fill1DHist(tri.M(), Form("00_tri_massZ_%s_cut%i",        d, num),";m_{ll#gamma}",  100,60,120,  weight, dir);
  hists->fill1DHist(tri.M(), Form("00_tri_mass_longTail_%s_cut%i",d, num),";m_{ll#gamma}",  100, 0,400,  weight, dir);

  hists->fill1DHist(lPt1Gamma.M(), Form("lPt1Gamma_%s_cut%i",d, num),";M(l1#gamma)",  100, 0,400,  weight, dir);
  hists->fill1DHist(lPt2Gamma.M(), Form("lPt2Gamma_%s_cut%i",d, num),";M(l2#gamma)",  100, 0,400,  weight, dir);

  hists->fill2DHist(diLep.M2(), lPt2Gamma.M2(), Form("h2D_dalitzPlot_ll_vs_gl_%s_cut%i",  d, num),
		    ";m^{2}_{ll} (GeV^{2});m^{2}_{l_{2}#gamma} (GeV^{2})", 100,0,20000, 100, 0,10000,  weight, dir);

  hists->fill2DHist(lPt1Gamma.M2(), lPt2Gamma.M2(), Form("h2D_dalitzPlot_l1g_vs_l2g_%s_cut%i",  d, num),
		    ";m^{2}_{l_{1}#gamma} (GeV^{2});m^{2}_{l_{2}#gamma} (GeV^{2})", 100,0,20000, 100, 0,10000,  weight, dir);

  double rotation = atan(-1.);
  double sin_rotation = sin(rotation);
  double cos_rotation = cos(rotation);
  double x_prime = (lPt1Gamma.M2()*cos_rotation - lPt2Gamma.M2()*sin_rotation) + (1.-cos_rotation)*pow(125.,2);
  double y_prime =  lPt1Gamma.M2()*sin_rotation + lPt2Gamma.M2()*cos_rotation;

  if (_isGammaSet){
    hists->fill2DHist(x_prime, y_prime, Form("h2D_dalitzPlot_l1g_vs_l2g_rotation_%s_cut%i",  d, num),
		      ";M(l1#gamma)^{2};M(l2#gamma)^{2}", 100,0,22000, 100, -14000,5000,  weight, dir);

    const UInt_t nBins1 = 60;
    hists->fill2DHist(tri.M(), diLep.M(), Form("h2D_tri_vs_diLep_mass0_%s_cut%i",  d, num),
		      ";m_{ll#gamma} (GeV);m_{ll} (GeV)", nBins1,0,200, nBins1, 0,20,  weight, dir);
    hists->fill2DHist(tri.M(), diLep.M(), Form("h2D_tri_vs_diLep_mass1_%s_cut%i",  d, num),
		      ";m_{ll#gamma} (GeV);m_{ll} (GeV)", nBins1,0,200, nBins1, 0,1.5, weight, dir);
    hists->fill2DHist(tri.M(), diLep.M(), Form("h2D_tri_vs_diLep_mass2_%s_cut%i",  d, num),
		      ";m_{ll#gamma} (GeV);m_{ll} (GeV)", nBins1,0,200, nBins1, 1.5,8, weight, dir);
    hists->fill2DHist(tri.M(), diLep.M(), Form("h2D_tri_vs_diLep_mass3_%s_cut%i",  d, num),
		      ";m_{ll#gamma} (GeV);m_{ll} (GeV)", nBins1,0,200, nBins1, 8,20,  weight, dir);

    hists->fill1DHist(diLep.Pt()/tri.M(),  Form("diLep_ptOverMllg_%s_cut%i",   d, num),
		      ";di-Lepton p_{T}/m_{ll#gamma}",100, 0,1, weight, dir);
  }

  hists->fill1DHist(diLep.M(), Form("01_diLep_mass_20_%s_cut%i",   d, num),";m_{#mu#mu} (GeV)", 100, 0,20,  weight, dir);
  hists->fill1DHist(diLep.M(), Form("01_diLep_mass_50_%s_cut%i",   d, num),";m_{#mu#mu} (GeV)", 100, 0,50,  weight, dir);
  hists->fill1DHist(diLep.M(), Form("01_diLep_mass_full_%s_cut%i", d, num),";m_{#mu#mu} (GeV)", 100, 0,120, weight, dir);
  hists->fill1DHist(diLep.M(), Form("01_diLep_mass_low1_%s_cut%i", d, num),";m_{ll} (GeV)",  50,   0,1.5,   weight, dir);
  hists->fill1DHist(diLep.M(), Form("01_diLep_mass_low2_%s_cut%i", d, num),";m_{ll} (GeV)", 100, 1.5,  8,   weight, dir);
  hists->fill1DHist(diLep.M(), Form("01_diLep_mass_low3_%s_cut%i", d, num),";m_{ll} (GeV)",  30,   8, 20,   weight, dir);
  hists->fill1DHist(diLep.M(), Form("01_diLep_mass_bump1_%s_cut%i",d, num),";m_{ll} (GeV)",  40,  15, 35,   weight, dir);
  hists->fill1DHist(diLep.M(), Form("01_diLep_mass_bump2_%s_cut%i",d, num),";m_{ll} (GeV)",  30,  35, 50,   weight, dir);
  hists->fill1DHist(diLep.M(), Form("01_diLep_mass_jpsi_%s_cut%i", d, num),";m_{ll} (GeV)",  50, 2.7,3.5,   weight, dir);

  hists->fill1DHist(TMath::Log10(diLep.M()), Form("01_diLep_mass_log1_%s_cut%i", d, num),";log10(m_{ll})", 100,   -1,7, weight, dir);
  hists->fill1DHist(TMath::Log10(diLep.M()), Form("01_diLep_mass_log2_%s_cut%i", d, num),";log10(m_{ll})", 100,   -1,3, weight, dir);
  hists->fill1DHist(TMath::Log10(diLep.M()), Form("01_diLep_mass_log3_%s_cut%i", d, num),";log10(m_{ll})", 100,-0.7,1.3,weight, dir);

  hists->fill1DHist(diLep.Pt(),  Form("diLep_pt_%s_cut%i",   d, num),";di-Lepton p_{T}", 50, 0,120, weight, dir);
  hists->fill1DHist(diLep.Eta(), Form("diLep_eta_%s_cut%i",  d, num),";di-Lepton eta",   50, -3.5,3.5, weight, dir);
  hists->fill1DHist(diLep.Phi(), Form("diLep_phi_%s_cut%i",  d, num),";di-Lepton phi",   50, -TMath::Pi(),TMath::Pi(), weight, dir);
  hists->fill1DHist(diLepCM.E(), Form("diLep_Ecom_%s_cut%i", d, num),";E(ll) in CoM",    50, 0,200,  weight, dir);


  hists->fill1DHist(_lPt1.Pt(),  Form("02_lPt1_pt_%s_cut%i", d, num), ";Leading lepton p_{T} (GeV)", 50, 0,100,    weight, dir);
  hists->fill1DHist(_lPt1.Eta(), Form("02_lPt1_eta_%s_cut%i",d, num), ";Leading lepton eta",   50, -3.5,3.5, weight, dir);
  hists->fill1DHist(_lPt1.Phi(), Form("02_lPt1_phi_%s_cut%i",d, num), ";Leading lepton phi",   50, -TMath::Pi(),TMath::Pi(), weight, dir);
  hists->fill1DHist(_lPt2.Pt(),  Form("02_lPt2_pt_%s_cut%i", d, num), ";Trailing lepton p_{T} (GeV)",50, 0,100,    weight, dir);
  hists->fill1DHist(_lPt2.Eta(), Form("02_lPt2_eta_%s_cut%i",d, num), ";Trailing lepton eta",  50, -3.5,3.5, weight, dir);
  hists->fill1DHist(_lPt2.Phi(), Form("02_lPt2_phi_%s_cut%i",d, num), ";Trailing lepton phi",  50, -TMath::Pi(),TMath::Pi(), weight, dir);


  //hists->fill2DHist(l1.Pt(),   l2.Pt(),   Form("h2D_l1_vs_l2_%s_cut%i",  d, num),
  //";l+ pt; l- pt",    50, 0,100, 50,0,100, weight, dir);

  hists->fill2DHist(_lPt1.Pt(), _lPt2.Pt(), Form("h2D_Pt2_vs_Pt1_%s_cut%i", d, num),
		    ";Leading lepton p_{T};Trailing lepton p_{T}",  50, 0,100, 50,0,100, weight, dir);
  if(_isGammaSet)
    {
      hists->fill1DHist(_gamma.Pt()/tri.M(), Form("03_gamma_ptOverMllg_%s_cut%i",   d, num),
			";Photon p_{T}/M_{ll#gamma}",  100, 0,1, weight, dir);
      hists->fill1DHist(gammaCM.E(), Form("03_gamma_Ecom_%s_cut%i",d, num),";Photon Energy in CoM",50, 0,120, weight, dir);
      hists->fill1DHist(_gamma.E(),  Form("03_gamma_E_%s_cut%i",   d, num),";Photon Energy",       50, 0,200, weight, dir);
      hists->fill1DHist(_gamma.Pt(), Form("03_gamma_pt_%s_cut%i",  d, num),";Photon p_{T} (GeV)",  50, 0,120, weight, dir);
      hists->fill1DHist(_gamma.Eta(),Form("03_gamma_eta_%s_cut%i", d, num),";Photon eta", 50, -3.5,3.5,       weight, dir);
      hists->fill1DHist(_gamma.Phi(),Form("03_gamma_phi_%s_cut%i", d, num),";Photon phi", 50, -TMath::Pi(),TMath::Pi(), weight, dir);

      hists->fill2DHist(diLep.Pt(),_gamma.Pt(),Form("h2D_diLep_vs_gamma_%s_cut%i",d, num),
			";q_{T}^{ll} (GeV); p_{T}^{#gamma} (GeV)", 50, 0,100, 50,0,100, weight, dir);
      hists->fill2DHist(diLep.Pt()/tri.M(),_gamma.Pt()/tri.M(),Form("h2D_diLep_vs_gamma_ptOverMllg_%s_cut%i",d, num),
			";q_{T}^{ll}/M_{ll#gamma}; p_{T}^{#gamma}/M_{ll#gamma}",  100, 0,1, 100,0,1, weight, dir);
      hists->fill2DHist(gammaCM.E(), _gamma.Pt(),Form("h2D_gamma_Ecom_vs_Pt_%s_cut%i",  d, num),
			";E_{#gamma} in CoM; p_{T}(#gamma)", 50, 0,100, 50,0,100, weight, dir);
      hists->fill2DHist(gammaCM.E(), tri.M(),   Form("h2D_gamma_Ecom_vs_triM_%s_cut%i",d, num),
			";E_{#gamma} in CoM; m_{ll#gamma}",   50, 0,100, 50,0,200, weight, dir);
      hists->fill2DHist(diLep.Eta()-_gamma.Eta(), _gamma.Pt(), Form("h2D_deltaEta_vs_gammaPt_deltaEta_%s_cut%i",d, num),
			";#Delta#eta(ll, #gamma);p_{T} of #gamma", 100, -5,5, 100,0,130, weight, dir);
    }

  hists->fill1DHist(_lPt1.DeltaR(_lPt2), Form("04_ll_deltaR_full_%s_cut%i",   d, num),";#Delta R(l_{1}, l_{2})",50,0,4,  weight,dir);
  hists->fill1DHist(_lPt1.DeltaR(_lPt2), Form("04_ll_deltaR_%s_cut%i",        d, num),";#Delta R(l_{1}, l_{2})",50,0,1,  weight,dir);
  if (_isGammaSet){
    hists->fill1DHist(_lPt1.DeltaR(_gamma),Form("04_lPt1_gamma_deltaR_%s_cut%i",d, num),";#Delta R(#gamma,l_{1})",100,0,5, weight,dir);
    hists->fill1DHist(_lPt2.DeltaR(_gamma),Form("04_lPt2_gamma_deltaR_%s_cut%i",d, num),";#Delta R(#gamma,l_{2})",100,0,5, weight,dir);
  }

  hists->fill1DHist(_rhoFactor, Form("rhoFactor_%s_cut%i",     d, num), ";#rho-factor",     100, 0,40, weight, dir);
  hists->fill1DHist(_nVtx,      Form("vtx_nPV_weight_%s_cut%i",d, num), ";nPV, re-weighted", 40, 0,40, weight, dir);
  hists->fill1DHist(_nVtx,      Form("vtx_nPV_raw_%s_cut%i",   d, num), ";nPV, Raw",         40, 0,40,      1, dir);
  //hists->fill1DHist(nVtxTotal, Form("vtx_nPV_tot_%s_cut%i",   d, num), "vtx_nPV_tot",    40, 0,40, weight, dir);

  hists->fill1DHist(_pv.NDof(), Form("vtx_ndof1_%s_cut%i",  d, num), ";First vertex nDof", 50, 0,200, weight, dir);
  //hists->fill1DHist(nDofVtx2, Form("vtx_ndof2_%s_cut%i",  d, num), ";Second vertex ndof", 50, 0,200, weight, dir);



  if (_isLepSet){
    TLorentzVector l1,l2;
    // l1 has to be positive and l2 is negative! charge!

    if (_lPt1.Charge()==1 && _lPt2.Charge()==-1)
      {l1 = _lPt1; l2 = _lPt2;}
    else if (_lPt1.Charge()==-1 && _lPt2.Charge()==1)
      {l1 = _lPt2; l2 = _lPt1;}
    else{
      return;
      //co1=co2=phi=co3=-5;
      //cout<<"Something is wrong: leptons have the same charge!"<<endl;
    }
    //cout<<eventNumber<<"RECO  Angles: c1= "<<co1<<"  c2="<<co2<<"   phi="<<phi<<"   coco="<<co3<<endl;


    //ZGanlgles:
    double co1,co2,phi,co3;
    if(_isGammaSet){
      angles->GetAngles(l1, l2, _gamma,co1,co2,phi,co3);

      hists->fill1DHist(co1, Form("co1_cut%i", num), ";cos(#vec{l-},#vec{Z}) in CM", 100,-1,1, weight,"Angles");
      hists->fill1DHist(co2, Form("co2_cut%i", num), ";cos(#vec{l+},#vec{Z}) in CM", 100,-1,1, weight,"Angles");
      hists->fill1DHist(co3, Form("co3_cut%i", num), ";cos(#Theta)",                 100,-1,1, weight,"Angles");
      hists->fill1DHist(phi, Form("phi_cut%i", num), ";#phi(l+)",50, -TMath::Pi(), TMath::Pi(),weight,"Angles");
    }

    //Phi* variable

    float phi_acop   = TMath::Pi() - fabs(l1.DeltaPhi(l2));
    float theta_star = TMath::ACos(TMath::TanH( (l2.Eta() - l1.Eta())/2));
    float phi_star   = sin(theta_star)*TMath::Tan(phi_acop/2);

    hists->fill1DHist(phi_acop, Form("star_phi_acop_cut%i", num), ";#phi_{acop}",  50, 0, TMath::Pi(), weight,"Angles");
    hists->fill1DHist(phi_star, Form("star_phi_star_cut%i", num), ";#phi*",        50, 0,         150, weight,"Angles");
    hists->fill1DHist(cos(theta_star), Form("star_cos_theta_cut%i", num), ";cos(#theta*)",  50, -1, 1, weight,"Angles");


    //Forward-Backward assymetry angles
    int sign = (diLep.Z() > 0) ? 1 : -1;
    float m = diLep.M();
    float cosCS = sign*2*(l2.Pz()*l1.E() - l1.Pz()*l2.E())/(m*sqrt(m + pow(diLep.Pt(),2)));

    hists->fill1DHist(cosCS, Form("FB_cosSC-ll_cut%i",  num), ";cosCS(l+,l-)", 100, -1, 1, weight,"Angles");
    hists->fill2DHist(cosCS, m, Form("FB_cosSC-ll_vs_Mll_cut%i",  num), ";cosCS(l+,l-);m_{ll}",  100, -1,1, 100,0,20, weight,"Angles");

    if (_isGammaSet){
      sign = (tri.Z() > 0) ? 1 : -1;
      m = tri.M();
      cosCS = sign*2*(diLep.Pz()*_gamma.E() - _gamma.Pz()*diLep.E())/(m*sqrt(m + pow(tri.Pt(),2)));

      hists->fill1DHist(cosCS, Form("FB_cosSC-diG_cut%i", num), ";cosCS(diLep,#gamma)",  100, -1, 1, weight,"Angles");
      hists->fill2DHist(cosCS, m, Form("FB_cosSC-diG_vs_Mllg_cut%i", num), ";cosCS(diLep,#gamma);m_{ll#gamma}",  100, -1,1, 100,50,140, weight,"Angles");

      float cosJ = _gamma.Dot(l1-l2)/_gamma.Dot(l1+l2);
      hists->fill1DHist(cosJ, Form("FB_cosJ_cut%i", num), ";cosJ",  100, -1, 1, weight,"Angles");
      hists->fill2DHist(cosJ, tri.M(), Form("FB_cosJ_vs_Mllg_cut%i", num), ";cosJ;m_{ll#gamma}",  100, -1,1, 100,50,140, weight,"Angles");
    }
  }
}


void HistMaker::MakePhotonEnergyCorrPlots(const TCPhoton& pho, Float_t corrPhoEnReco, Float_t corrPhoEnSC)
{
  const string dir = "Photon-EnergyCorr";
  hists->fill1DHist(pho.M(), "ph_recoMass",";M_{#gamma}", 100,-0.05,0.05, 1,dir);

  hists->fill1DHist((pho.E() - corrPhoEnReco)/pho.E(),"ph_energyCorrectionReco",
		    ";(E_{#gamma}^{reco} - E_{#gamma}^{reco corr})/E_{#gamma}^{reco}", 100,-0.1,0.1, 1,dir);

  hists->fill1DHist((pho.SCEnergy() - corrPhoEnSC)/pho.SCEnergy(),"ph_energyCorrectionSC",
		    ";(E_{#gamma}^{SC} - E_{#gamma}^{SC corr})/E_{#gamma}^{SC}", 100,-0.1,0.1, 1,dir);

  hists->fill1DHist((pho.E() - pho.SCEnergy())/pho.E(),"ph_energyRecoVsSC",
		    ";(E_{#gamma}^{reco} - E_{#gamma}^{SC})/E_{#gamma}^{reco}", 100,-0.1,0.1, 1,dir);

  if (pho.Pt()>40 && fabs(pho.SCEta())<1.444){
    hists->fill1DHist((pho.E() - pho.SCEnergy())/pho.E(),"ph_energyRecoVsSC_PtEtaCut",
		      ";(E_{#gamma}^{reco} - E_{#gamma}^{SC})/E_{#gamma}^{reco}", 100,-0.1,0.1, 1,dir);

    hists->fill1DHist((pho.E() - corrPhoEnReco)/pho.E(),"ph_energyCorrectionReco_PtEtaCut",
		      ";(E_{#gamma}^{reco} - E_{#gamma}^{reco corr})/E_{#gamma}^{reco}", 100,-0.1,0.1, 1,dir);

    hists->fill1DHist((pho.SCEnergy() - corrPhoEnSC)/pho.SCEnergy(),"ph_energyCorrectionSC_PtEtaCut",
		      ";(E_{#gamma}^{SC} - E_{#gamma}^{SC corr})/E_{#gamma}^{SC}", 100,-0.1,0.1, 1,dir);

    /* These are too much, no need anymore
    if (pho.R9() > 0.94)
      {
	hists->fill1DHist((pho.E() - corrPhoEnReco)/pho.E(),"ph_energyCorrectionReco_PtEtaCut_hightR9",
			  ";(E_{#gamma}^{reco} - E_{#gamma}^{reco corr})/E_{#gamma}^{reco}", 100,-0.1,0.1, 1,dir);

	hists->fill1DHist((pho.SCEnergy()-corrPhoEnSC)/pho.SCEnergy(),"ph_energyCorrectionSC_PtEtaCut_highR9",
			  ";(E_{#gamma}^{SC} - E_{#gamma}^{SC corr})/E_{#gamma}^{SC}", 100,-0.1,0.1, 1,dir);
      }
    else
      {
	hists->fill1DHist((pho.E() - corrPhoEnReco)/pho.E(),"ph_energyCorrectionReco_PtEtaCut_lowR9",
			  ";(E_{#gamma}^{reco} - E_{#gamma}^{reco corr})/E_{#gamma}^{reco}", 100,-0.1,0.1, 1,dir);

	hists->fill1DHist((pho.SCEnergy()-corrPhoEnSC)/pho.SCEnergy(),"ph_energyCorrectionSC_PtEtaCut_lowR9",
			  ";(E_{#gamma}^{SC} - E_{#gamma}^{SC corr})/E_{#gamma}^{SC}", 100,-0.1,0.1, 1,dir);

      }
    */
  }
}

void HistMaker::MakeZeePlots(const TCPhoton& p1, const TCPhoton& p2)
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


void HistMaker::MakeNPlots(Int_t num, Int_t nmu, Int_t nele1, Int_t nele2, Int_t nphoHZG, Int_t nphoTight, Int_t nphoMVA, Int_t nJets, Double_t w)
{
  hists->fill1DHist(nmu,       Form("size_mu_cut%i",  num), ";Number of muons",              5,0,5, w, "N");
  hists->fill1DHist(nele1,     Form("size_el_cut%i",  num), ";Number of electrons (HZZ)",    5,0,5, w, "N");
  hists->fill1DHist(nele2,     Form("size_el0_cut%i", num), ";Number of electrons (Loose)",  5,0,5, w, "N");
  hists->fill1DHist(nphoHZG,   Form("size_phHZG_cut%i",  num), ";Number of photons (HZG)",   5,0,5, w, "N");
  hists->fill1DHist(nphoTight, Form("size_phTight_cut%i",num), ";Number of photons (Tight)", 5,0,5, w, "N");
  hists->fill1DHist(nphoMVA,   Form("size_phMVA_cut%i",  num), ";Number of photons (MVA)",   5,0,5, w, "N");
  hists->fill1DHist(nJets,     Form("size_jets_cut%i",   num), ";Number of Jets",            5,0,5, w, "N");
}
