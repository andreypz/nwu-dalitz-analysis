TList *FileList;
TFile *Target;

#define nC 10
#define F0 5
#define F1 6
#define F2 9
#define F3 4

struct optimalCuts {
  Float_t SB;
  Float_t SrootB;
  Float_t SrootSB;
  Int_t bin1;
  Int_t bin2;
};

void makePlot(Int_t sel=1, TString hPath="00")
{
  gROOT->ProcessLine(".L ../data/tdrstyle.C");
  setTDRStyle();

  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLegendBorderSize(0);
  cout.precision(3); cout.setf(ios::fixed, ios::floatfield);
  TH1::SetDefaultSumw2(kTRUE);

  TString histoPath = hPath.Data();
  Float_t intLumi = 191.; //Note: in v11 and below the histograms are already normylized to 191, after that - to 1000
  //Float_t intLumi = 191.0 + 963.6; 
 
  cout<<"histoPath:  "<<histoPath.Data()<<"  int Lumi: "<<intLumi<<endl;

  //Types of met: met - pfMet, met1 - type1 corrected, met2 - pfMet passed Noise filters, 
  //met3 - projMet, met4 - puProj corrected met (those two are passed Noise filters) 
  TString metType("met3");   TString mtType("mt2"); 

  if (sel==1)  TString imgpath("~/afs/public_html/higgs/overview/muon/");  
  if (sel==2)  TString imgpath("~/afs/public_html/higgs/overview/electron/");  

  TString testpath("~/afs/public_html/test/");  
  Bool_t doPhotons = 1, doZjets=0, doEBEE=0, doOverview=0;
  Bool_t doSB =0, doTest=1;

  TFile* fda_2011A_DoubleMu  = new TFile(Form("./%s/dir_1_2011A_May10_DoubleMu_/hhhh.root",histoPath.Data()));
  TFile* fda_2011A_DoubleEl  = new TFile(Form("./%s/dir_2_2011A_May10_DoubleElectron_/hhhh.root",histoPath.Data()));
  if (sel==1) TFile  *fData = (TFile*)fda_2011A_DoubleMu;
  if (sel==2) TFile  *fData = (TFile*)fda_2011A_DoubleEl;
  TFile* fda_Data  = new TFile(Form("./m_Data_%i.root", sel));  //Merged Data
  // TFile  *fData = (TFile*)fda_Data;
 
  TFile *fA_Zj = new TFile(Form("./forAnton_Zjets_%i.root",sel), "RECREATE");
  TFile *fA_Da = new TFile(Form("./forAnton_Data_%i.root",sel), "RECREATE");
  TFile *fA_Ds = new TFile(Form("./forAnton_Data_sbtr_%i.root",sel), "RECREATE");

  TFile *fPhot = new TFile("./photonResult_070211_2.root", "OPEN");

  TFile* fmc_ZllG      = new TFile(Form("./%s/dir_%i_MC_ZllG_/hhhh.root", histoPath.Data(), sel));
  TFile* fmc_Wjets     = new TFile(Form("./%s/dir_%i_MC_Wjets_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_ZZtoAny   = new TFile(Form("./%s/dir_%i_MC_ZZtoAny_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_tt        = new TFile(Form("./%s/dir_%i_MC_tt2l2nu2b_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_WW        = new TFile(Form("./%s/dir_%i_MC_WW_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_WZ        = new TFile(Form("./%s/dir_%i_MC_WZ_/hhhh.root",histoPath.Data(), sel));

  TFile* fmc_ggH200    = new TFile(Form("./%s/dir_%i_MC_ggH200_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_ggH400    = new TFile(Form("./%s/dir_%i_MC_ggH400_/hhhh.root",histoPath.Data(), sel));
  TFile* fmc_ggH600    = new TFile(Form("./%s/dir_%i_MC_ggH600_/hhhh.root",histoPath.Data(), sel));
  //TFile* fmc_  = new TFile(Form("./%s/dir_%i__/hhhh.root",histoPath.Data(), sel));

  TFile* fmc_Zjets  = new TFile(Form("./m_Zjets_%i.root", sel));
  TFile* fmc_Top    = new TFile(Form("./m_Top_%i.root", sel));
  //TFile* fmc_ZQQ    = new TFile(Form("./m_ZQQ_%i.root", sel));

  //List of background samples to Stack
  list_bg = new TList();
  list_bg->Add(fmc_Top);
  list_bg->Add(fmc_tt);
  list_bg->Add(fmc_WZ);
  list_bg->Add(fmc_WW);
  list_bg->Add(fmc_ZZtoAny);
  list_bg->Add(fmc_Zjets);
  list_bg->Add(fmc_Wjets);

  list_bg2 = new TList();
  list_bg2->Add(fmc_Top);
  list_bg2->Add(fmc_tt);
  list_bg2->Add(fmc_WZ);
  list_bg2->Add(fmc_WW);
  list_bg2->Add(fmc_ZZtoAny);
  list_bg2->Add(fmc_Wjets);

  THStack *hs_met_et[nC], *hs_met2_et[nC], *hs_met_over_qt[nC], *hs_met2_over_qt[nC], *hs_di_qt[nC], *hs_met_et_ovQt[nC], *hs_mt[nC], *hs_mtZ[nC];
  THStack *hs_jet_N[nC], *hs_jet_b_N[nC], hs_jet_b_pt[nC];
  THStack *hs_di_mass[nC], *hs_di_mass_EB[nC], *hs_di_mass_EE[nC], *hs_di_mass_EX[nC];
  THStack *hs_met_dPhiLeadJet1[nC], *hs_met_dPhiLeadJet2[nC], *hs_met_dPhiClosJet1[nC], *hs_met_dPhiClosJet2[nC];
  for(Int_t n = 0; n<nC; n++)
    {
      hs_met_et[n]       = makeStack(list_bg, Form("%s_et_%i", metType.Data(), n), intLumi);
      hs_met2_et[n]      = makeStack(list_bg, Form("met2_et_%i", n), intLumi);
      hs_mt[n]           = makeStack(list_bg, Form("%s_%i", mtType.Data(), n), intLumi);
      hs_mtZ[n]          = makeStack(list_bg, Form("mtZ_%i", n), intLumi);

      hs_met_over_qt[n]  = makeStack(list_bg, Form("%s_over_qt_%i", metType.Data(), n), intLumi);
      hs_met2_over_qt[n] = makeStack(list_bg, Form("met2_over_qt_%i", n), intLumi);
      hs_met_et_ovQt[n]  = makeStack(list_bg, Form("%s_et_ovQt_%i", metType.Data(), n), intLumi);

      hs_di_qt[n]   = makeStack(list_bg, Form("di_qt_%i",n), intLumi);
      hs_di_mass[n] = makeStack(list_bg, Form("di_mass_%i",n), intLumi);
      hs_di_mass_EB[n] = makeStack(list_bg, Form("di_mass_EB_%i",n), intLumi);
      hs_di_mass_EE[n] = makeStack(list_bg, Form("di_mass_EE_%i",n), intLumi);
      hs_di_mass_EX[n] = makeStack(list_bg, Form("di_mass_EX_%i",n), intLumi);
      hs_jet_N[n]   = makeStack(list_bg, Form("jet_N_%i",n), intLumi);
      hs_jet_b_N[n] = makeStack(list_bg, Form("jet_b_N_%i",n), intLumi); 
    
      hs_met_dPhiLeadJet1[n] = makeStack(list_bg, Form("met2_dPhiLeadJet1_%i",n), intLumi);
      hs_met_dPhiLeadJet2[n] = makeStack(list_bg, Form("met2_dPhiLeadJet2_%i",n), intLumi);
      hs_met_dPhiClosJet1[n] = makeStack(list_bg, Form("met2_dPhiClosJet1_%i",n), intLumi);
      hs_met_dPhiClosJet2[n] = makeStack(list_bg, Form("met2_dPhiClosJet2_%i",n), intLumi);
    }

  //  hs_jet_b_pt[8]  = makeStack(list_bg, Form("jet_b_pt_%i",8)); //_ to be fixed


  hs_met_et[6] -> PrintYields(sel, 6);
  
  TIter next(hs_jet_N[F0] -> GetHists());
  TH1 * forLegend[21];
  Int_t aa = 0;
  TObject *obj;
  while ((obj = next()))
    {
      forLegend[aa] = (TH1F*)obj;
      //cout<<obj->GetName()<<endl;
      aa++;
    }

  fda_2011A_DoubleMu->cd();
  forLegend[20] = (TH1F*)met3_et_0;
  fmc_ZllG -> cd();
  forLegend[19] = (TH1F*)met3_et_0;
  fmc_ggH400 -> cd();
  forLegend[18] = (TH1F*)met3_et_0;
  fmc_ggH200 -> cd();
  met3_et_0 -> SetLineColor(kBlack); ///need a fix
  forLegend[17] = (TH1F*)met3_et_0;
  fmc_tt -> cd();
  forLegend[16] = (TH1F*)met3_et_0;

 
  leg01 = new TLegend(0.53,0.7,0.95,0.95);
  leg01 -> SetNColumns(2);
  leg01 -> SetTextSize(0.04);
  
  leg01->AddEntry(forLegend[20], "Data","epl");
  leg01->AddEntry(forLegend[5],  "Z + jets","f");
  leg01->AddEntry(forLegend[4],  "ZZ","f");
  leg01->AddEntry(forLegend[6],  "W + jets","f");
  leg01->AddEntry(forLegend[3],  "WW","f");
  leg01->AddEntry(forLegend[0],  "t#rightarrow l#nub","f");
  leg01->AddEntry(forLegend[2],  "WZ","f");
  leg01->AddEntry(forLegend[16], "tt#rightarrow 2l2#nu2b","f");  leg01->AddEntry(forLegend[19], "Z#gamma#rightarrow ll#gamma","f");
  leg01->AddEntry(forLegend[17], "10xH200","f");
  leg01->AddEntry(forLegend[17]," ","");  //empty slot
  leg01->AddEntry(forLegend[18], "10xH400","f");
  leg01->SetFillColor(kWhite);
  leg01->SetFillColor(kWhite);
   
  hs_met_et[F0] -> Draw("hist");
  c1 -> SaveAs(imgpath+"ov01.png");
  TCanvas *c2 = new TCanvas("c2","example",600,700);



  if(doSB){
    cout<<"doing S/B"<<endl;
    optimalCuts result;

    result =  calculateOptimalCuts(1, hs_met_over_qt[6]->Sum(), (TH1*)fmc_ggH400->Get("met3_over_qt_6"));
    cout<<"from "<<0.1*(result.bin1-1)<<" to "<< 0.1*result.bin2<<"\n S/B = "<<result.SB<<"  S/sqrt(B) = "<<result.SrootB<<"  S/sqrt(S+B) = "<<result.SrootSB<<endl;
    
    Float_t iSig = ((TH1*)fmc_ggH400->Get("met3_over_qt_6")) -> Integral(7,18); //from 0.6 to 1.8
    Float_t iBkg = (hs_met_over_qt[6]->Sum()) -> Integral(7,18);
    Float_t SB      = 0.1*iSig/iBkg;        //Higgs signal is 10*real in the histograms
    Float_t SrootB  = 0.1*iSig/sqrt(iBkg);
    Float_t SrootSB = 0.1*iSig/sqrt(iBkg + 0.1*iSig);
    cout<<"After cuts of projMET/qT from 0.6 to 1.8 "<<endl;
    cout<<iSig<<"  "<<iBkg<<" SB = "<<SB<<" SrootB = "<<SrootB<<" SrootSB = "<<SrootSB<<endl;
    cout<<endl;
    
    result =  calculateOptimalCuts(3, hs_mt[7]->Sum(), (TH1*)fmc_ggH400->Get("mt2_7"));
    cout<<"from "<<100+10*(result.bin1-1)<<" to "<< 100+10*result.bin2<<"\n S/B = "<<result.SB<<"  S/sqrt(B) = "<<result.SrootB<<"  S/sqrt(S+B) = "<<result.SrootSB<<endl;

    result =  calculateOptimalCuts(3, hs_mt[6]->Sum(), (TH1*)fmc_ggH400->Get("mt2_6"));
    cout<<"from "<<100+10*(result.bin1-1)<<" to "<< 100+10*result.bin2<<"\n S/B = "<<result.SB<<"  S/sqrt(B) = "<<result.SrootB<<"  S/sqrt(S+B) = "<<result.SrootSB<<endl;

    if(doTest){
      drawMuliPlot("projMET", 1, 0.001, 1000000, 0,5, hs_met_et[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(testpath+"p01.png");
      drawMuliPlot("projMET/q_{T}", 1, 0.001, 1000000, 0,3, hs_met_over_qt[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(testpath+"p02.png");
      drawMuliPlot("MT", 1, 0.001, 1000000, 0,3, hs_mt[7], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(testpath+"p03.png");
      drawMuliPlot("MT", 1, 0.001, 1000000, 0,3, hs_mt[6], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(testpath+"p05.png");
      
    }
  }
  //-----------------------------------------------//
  // Zjets  backgroung from Photon samples (Anton) //
  //-----------------------------------------------//
  if(doPhotons){
    if (sel==1)  TString imgpath("~/afs/public_html/nwu/110707_photon/muon/");  
    if (sel==2)  TString imgpath("~/afs/public_html/nwu/110707_photon/electron/");  
    
    THStack *ph_met_et[nC], *ph_met2_et[nC], *ph_met_over_qt[nC], *ph_met2_over_qt[nC], *ph_di_qt[nC], *ph_met_et_ovQt[nC], *ph_mt[nC], *ph_mtZ[nC];
    THStack *ph_jet_N[nC], *ph_jet_b_N[nC], ph_jet_b_pt[nC];
    for(Int_t n = 0; n<nC; n++)
      {
	ph_met2_et[n]      = makeStack(list_bg2, Form("met2_et_%i", n), intLumi);
	ph_mt[n]           = makeStack(list_bg2, Form("%s_%i", mtType.Data(), n), intLumi);
	
	ph_met_over_qt[n]  = makeStack(list_bg2, Form("met2_over_qt_%i", n), intLumi);
	ph_jet_N[n]        = makeStack(list_bg2, Form("jet_N_%i",n), intLumi);
	
	ph_di_qt[n]        = makeStack(list_bg2, Form("di_qt_%i",n), intLumi);
      }
    
    fPhot ->cd();
    
    if(sel==1){
      AddPhoton(photon2_weighted_pt, ph_di_qt[4], intLumi);
      AddPhoton(photon2_weighted_met_0, ph_met2_et[5], intLumi);
      AddPhoton(photon2_weighted_met_2, ph_met2_et[6], intLumi);
      
      AddPhoton(photon2_weighted_mt_0, ph_mt[5], intLumi);
      AddPhoton(photon2_weighted_mt_2, ph_mt[6], intLumi);
      
      AddPhoton(photon2_weighted_nJets_0, ph_jet_N[5], intLumi);
      AddPhoton(photon2_weighted_nJets_2, ph_jet_N[6], intLumi);
      
      AddPhoton(photon2_weighted_metOverQt_0, ph_met_over_qt[5], intLumi);
      AddPhoton(photon2_weighted_metOverQt_2, ph_met_over_qt[6], intLumi);
      
    }
    if(sel==2){
      AddPhoton(photon1_weighted_pt, ph_di_qt[4], intLumi);
      AddPhoton(photon1_weighted_met_0, ph_met2_et[5], intLumi);
      AddPhoton(photon1_weighted_met_2, ph_met2_et[6], intLumi);
      
      AddPhoton(photon1_weighted_mt_0, ph_mt[5], intLumi);
      AddPhoton(photon1_weighted_mt_2, ph_mt[6], intLumi);
      
      AddPhoton(photon1_weighted_nJets_0, ph_jet_N[5], intLumi);
      AddPhoton(photon1_weighted_nJets_2, ph_jet_N[6], intLumi);
      
      AddPhoton(photon1_weighted_metOverQt_0, ph_met_over_qt[5], intLumi);
      AddPhoton(photon1_weighted_metOverQt_2, ph_met_over_qt[6], intLumi);
    }
    
    
    drawMuliPlot("q_{T}", 1, 0.01, 100000, 0,3, ph_di_qt[4], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ph01.png");
    drawMuliPlot("q_{T}", 1, 0.01, 100000, 0,3, hs_di_qt[4], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ph02.png");
    
    drawMuliPlot("pfMET", 1, 0.01, 100000, 0,5, ph_met2_et[5], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ph03.png");
    drawMuliPlot("pfMET", 1, 0.01, 100000, 0,5, hs_met2_et[5], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ph04.png");
    
    drawMuliPlot("pfMET", 1, 0.001, 10000, 0,5, ph_met2_et[6], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ph05.png");
    drawMuliPlot("pfMET", 1, 0.001, 10000, 0,5, hs_met2_et[6], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ph06.png");
    
    drawMuliPlot("MT", 1, 0.001, 100000, 0,5, ph_mt[6], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ph07.png");
    drawMuliPlot("MT", 1, 0.001, 100000, 0,5, hs_mt[6], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ph08.png");
    
    drawMuliPlot("nJets", 1, 0.001, 100000, 0,5, ph_jet_N[6], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ph09.png");
    drawMuliPlot("nJets", 1, 0.001, 100000, 0,5, hs_jet_N[6], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ph10.png");
    
    drawMuliPlot("pfMet/q_{T}", 1, 0.001, 100000, 0,5, ph_met_over_qt[5], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ph11.png");
    drawMuliPlot("pfMet/q_{T}", 1, 0.001, 100000, 0,5, hs_met2_over_qt[5], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ph12.png");
    
    drawMuliPlot("pfMet/q_{T}", 1, 0.001, 100000, 0,5, ph_met_over_qt[6], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ph13.png");
    drawMuliPlot("pfMet/q_{T}", 1, 0.001, 100000, 0,5, hs_met2_over_qt[6], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ph14.png");
   
  }

  // Mll in electons, EE vs EB
  if(doEBEE){
    drawMuliPlot("Leptons in Barrel, |#eta|<1.444,  M(ll)", 1, 0.001, 1000000, 0,4, hs_di_mass_EB[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(testpath+"p01.png");
    drawMuliPlot("Leptons in Endcap, |#eta|>1.566, M(ll)", 1, 0.001, 1000000, 0,4, hs_di_mass_EE[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(testpath+"p02.png");
    drawMuliPlot("Leptons in EB/EE, mixed, M(ll)", 1, 0.001, 1000000, 0,4, hs_di_mass_EX[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(testpath+"p03.png");
  }

  if(doOverview){
    drawMuliPlot("projMET", 1, 0.001, 1000000, 0,5, hs_met_et[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ov01.png");
    drawMuliPlot("projMET/q_{T}", 1, 0.001, 1000000, 0,3, hs_met_over_qt[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ov02.png");
    
    drawMuliPlot("projMET", 1, 0.001, 1000000, 0,5, hs_met_et[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ov03.png");
    drawMuliPlot("projMET/q_{T}", 1, 0.001, 1000000, 0,5, hs_met_over_qt[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ov04.png");
    
    drawMuliPlot("projMET", 1, 0.001, 1000000, 0,3, hs_met_et[F2], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ov05.png");
    drawMuliPlot("projMET/q_{T}", 1, 0.001, 1000000, 0,3, hs_met_over_qt[F2], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ov06.png");
    
    drawMuliPlot("M(ll)", 1, 0.001, 1000000, 0,2, hs_di_mass[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ov07.png");
    drawMuliPlot("MT", 1, 0.001, 1000000, 0,2, hs_mt[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ov08.png");
    
    drawMuliPlot("M(ll)", 1, 0.001, 1000000, 0,5, hs_di_mass[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ov09.png");
    drawMuliPlot("MT", 1, 0.001, 1000000, 0,5, hs_mt[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
    c2 -> SaveAs(imgpath+"ov10.png");
    
    /*
      drawMuliPlot("M(ll)", 1, 0.001, 10000000, 0,2, hs_di_mass[F2], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov11.png");
      drawMuliPlot("MT", 1, 0.001, 1000000, 0,2, hs_mt[F2], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov12.png");
      
      drawMuliPlot("N jets", 1, 0.001, 1000000, 0,5, hs_jet_N[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov13.png");
      drawMuliPlot("N b-jets", 1, 0.001, 1000000, 0,5, hs_jet_b_N[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov14.png");
      
      drawMuliPlot("N jets", 1, 0.001, 1000000, 0,5, hs_jet_N[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov15.png");
      drawMuliPlot("N b-jets", 1, 0.001, 1000000, 0,5, hs_jet_b_N[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov16.png");
      
      
      drawMuliPlot("N jets", 1, 0.001, 1000000, 0,5, hs_jet_N[F3], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov17.png");
      drawMuliPlot("N b-jets", 1, 0.001, 1000000, 0,5, hs_jet_b_N[F3], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov18.png");
      
      
      // drawMuliPlot("N jets", 1, 0.001, 1000000, 0,5, hs_jet_N[F2], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      // c2 -> SaveAs(imgpath+"ov17.png");
      // drawMuliPlot("N b-jets", 1, 0.001, 1000000, 0,5, hs_jet_b_N[F2], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      //   c2 -> SaveAs(imgpath+"ov18.png");
      
      
      drawMuliPlot("N b-jets", 1, 0.001, 1000000, 0,5, hs_jet_b_N[7], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov19.png");
      drawMuliPlot("N b-jets", 1, 0.001, 1000000, 0,5, hs_jet_b_N[8], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov20.png");
      
      drawMuliPlot("#Delta#phi(MET, lead jet), p_{T}>20, |#eta|<2.4", 1, 0.001, 10000000, 0,2, hs_met_dPhiLeadJet1[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov21.png");
      drawMuliPlot("#Delta#phi(MET, closest jet), p_{T}>20, |#eta|<2.4", 1, 0.001, 10000000, 0,2, hs_met_dPhiClosJet1[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov22.png");
      
      drawMuliPlot("#Delta#phi(MET, lead jet), p_{T}>20, |#eta|<2.4", 1, 0.001, 1000000, 0,5, hs_met_dPhiLeadJet1[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov23.png");
      drawMuliPlot("#Delta#phi(MET, closest jet), p_{T}>20, |#eta|<2.4", 1, 0.001, 1000000, 0,5, hs_met_dPhiClosJet1[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov24.png");
      
      drawMuliPlot("#Delta#phi(MET, lead jet), p_{T}>20, |#eta|<2.4", 1, 0.001, 10000000, 0,2, hs_met_dPhiLeadJet1[F2], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov25.png");
      drawMuliPlot("#Delta#phi(MET, closest jet), p_{T}>20, |#eta|<2.4", 1, 0.001, 10000000, 0,2, hs_met_dPhiClosJet1[F2], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov26.png");
      
      
      drawMuliPlot("q_{T} (di-lepton p_{T})", 1, 0.001, 1000000, 0,5, hs_di_qt[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov27.png");
      drawMuliPlot("q_{T} (di_lepton p_{T})", 1, 0.001, 1000000, 0,5, hs_di_qt[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov28.png");
      
      drawMuliPlot("pfMET", 1, 0.001, 1000000, 0,5, hs_met2_et[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov29.png");
      drawMuliPlot("MT (using Z mass)", 1, 0.001, 1000000, 0,5, hs_mtZ[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov30.png");
      
      drawMuliPlot("pfMET", 1, 0.001, 1000000, 0,5, hs_met2_et[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov31.png");
      drawMuliPlot("MT (using Z mass)", 1, 0.001, 1000000, 0,5, hs_mtZ[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG, intLumi);
      c2 -> SaveAs(imgpath+"ov32.png");
    */
  
  }

  if(doZjets){
    fA_Zj -> cd();
    TH1F *Zjets_qt   = fmc_Zjets->Get("di_qt_4") ->Clone();
    TH1F *Zjets_mass = fmc_Zjets->Get("di_mass_4") ->Clone();
    TH1F *Zjets_met  = fmc_Zjets->Get("met3_et_5") ->Clone();
    TH1F *Zjets_met2 = fmc_Zjets->Get("met2_et_5") ->Clone();
    TH1F *Zjets_mt5  = fmc_Zjets->Get("mt2_5") ->Clone();
    TH1F *Zjets_mt6  = fmc_Zjets->Get("mt2_6") ->Clone();
    TH1F *Zjets_nj5  = fmc_Zjets->Get("jet_N_5") ->Clone();
    TH1F *Zjets_nj6  = fmc_Zjets->Get("jet_N_6") ->Clone();
    TH1F *Zjets_mOq5 = fmc_Zjets->Get("met2_over_qt_5") ->Clone();
    TH1F *Zjets_mOq6 = fmc_Zjets->Get("met2_over_qt_6") ->Clone();
    Zjets_qt   -> Write();
    Zjets_mass -> Write();
    Zjets_met  -> Write();
    Zjets_met2 -> Write();
    Zjets_mt5  -> Write();
    Zjets_mt6  -> Write();
    Zjets_nj5  -> Write();
    Zjets_nj6  -> Write();
    Zjets_mOq5 -> Write();
    Zjets_mOq6 -> Write();
    //fA_Zj -> Close();
    
    
    fA_Da -> cd();
    TH1F *da_qt   = fData->Get("di_qt_4") ->Clone();
    TH1F *da_mass = fData->Get("di_mass_4") ->Clone();
    TH1F *da_met  = fData->Get("met3_et_5") ->Clone();
    TH1F *da_met2 = fData->Get("met2_et_5") ->Clone();
    TH1F *da_mt5  = fData->Get("mt2_5") ->Clone();
    TH1F *da_mt6  = fData->Get("mt2_6") ->Clone();
    TH1F *da_nj5  = fData->Get("jet_N_5") ->Clone();
    TH1F *da_nj6  = fData->Get("jet_N_6") ->Clone();
    TH1F *da_mOq5 = fData->Get("met2_over_qt_5") ->Clone();
    TH1F *da_mOq6 = fData->Get("met2_over_qt_6") ->Clone();
    da_qt   -> Write();
    da_mass -> Write();
    da_met  -> Write();
    da_met2 -> Write();
    da_mt5  -> Write();
    da_mt6  -> Write();
    da_nj5  -> Write();
    da_nj6  -> Write();
    da_mOq5 -> Write();
    da_mOq6 -> Write();
    //fA_Da -> Close();
    
    
    fA_Ds -> cd();
    TH1F *ds_qt   = (TH1F*)hs_di_qt[4]->Sum()->Clone();
    TH1F *ds_mass = (TH1F*)hs_di_mass[4]->Sum()->Clone();
    TH1F *ds_met  = (TH1F*)hs_met_et[5]->Sum()->Clone();
    TH1F *ds_met2 = (TH1F*)hs_met2_et[5]->Sum()->Clone();
    TH1F *ds_mt5 = (TH1F*)hs_mt[5]->Sum()->Clone();
    TH1F *ds_mt6 = (TH1F*)hs_mt[6]->Sum()->Clone();
    TH1F *ds_nj5 = (TH1F*)hs_jet_N[5]->Sum()->Clone();
    TH1F *ds_nj6 = (TH1F*)hs_jet_N[6]->Sum()->Clone();
    TH1F *ds_mOq5 = (TH1F*)hs_met2_over_qt[5]->Sum()->Clone();
    TH1F *ds_mOq6 = (TH1F*)hs_met2_over_qt[6]->Sum()->Clone();
    
    da_qt -> Add(ds_qt,-1);
    da_mass -> Add(ds_mass,-1);
    da_met -> Add(ds_met,-1);
    da_met2 -> Add(ds_met2,-1);
    da_mt5 -> Add(ds_mt5,-1);
    da_mt6 -> Add(ds_mt6,-1);

    da_nj5 -> Add(ds_nj5,-1);
    da_nj6 -> Add(ds_nj6,-1);
    
    da_mOq5 -> Add(ds_mOq5,-1);
    da_mOq6 -> Add(ds_mOq6,-1);
    
    da_qt -> Write();
    da_mass -> Write();
    da_met  -> Write();
    da_met2 -> Write();
    da_mt5  -> Write();
    da_mt6  -> Write();
    da_nj5  -> Write();
    da_nj6  -> Write();
    da_mOq5 -> Write();
    da_mOq6 -> Write();
    fA_Ds -> Close();
    fA_Da -> Close();
    fA_Zj -> Close();
  }


  /*

  c1 -> SetLogy(0);
  fmc_ggH200 -> cd();
  met3_et_ovQt_9 -> Draw("colz");
  met3_et_ovQt_9 -> SetTitle(";projMET; projMET/q_{T}");
  c1 -> SaveAs(imgpath+"ov05.png");

  hs_met_et_ovQt[F1] -> Draw("colz");
  hs_met_et_ovQt[F1] -> SetTitle(";projMET; projMET/q_{T}");
    c1 -> SaveAs(imgpath+"ov06.png");
  */

  cout<<"\n\n   -- end of job  --"<<endl;
  
}

THStack* makeStack(TList *sourcelist, TString name, Float_t lumi)
{
  //cout<<"making stack!"<<endl;
  THStack * hs = new THStack(Form("hs_11_%i",1),"Stacked MT");

  TFile *first_source = (TFile*)sourcelist->First();
  first_source->cd();
  TDirectory *current_sourcedir = gDirectory;

  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key, *oldkey=0;
  while ( (key = (TKey*)nextkey())) 
    {
      if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;
      first_source->cd();
 
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1")) continue;
	
	if(strncmp(Form("%s",name.Data()),key->GetName(), name.Length()) ==0)
	  {
	    //cout<<"first "<<key->GetName()<<endl;
	    //hs -> Add((TH1*)key->ReadObj());
	    TH1 *hh1 = (TH1*)key->ReadObj()->Clone();
	    hh1 -> Scale(lumi/1000);
	    hs -> Add(hh1);
	    //hs -> Add((TH1*)key->ReadObj()->Clone());
	    //hs -> Add( ((TH1*)key->ReadObj()->Clone())->Scale(lumi/1000));
	  }
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      while ( nextsource )
	{
	  nextsource->cd();
	  TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(key->GetName());
	  
	  if (key2) 
	    {
	      //TH1 *h2 = (TH1*)key2->ReadObj();
	      if(strncmp(Form("%s",name.Data()),key2->GetName(), name.Length()) ==0)
		{
		  //cout<<"nextsource "<<key2->GetName()<<endl;
		  TH1 *hh2 = (TH1*)key2->ReadObj()->Clone();
		  hh2 -> Scale(lumi/1000);
		  hs -> Add(hh2);
		  //hs -> Add((TH1*)key2->ReadObj()->Clone());
		  //hs -> Add( ((TH1*)key2->ReadObj()->Clone())->Scale(lumi/1000));
		}
	      //delete h2;
	    }
	  
	  nextsource = (TFile*)sourcelist->After( nextsource );
	}
      
    }



  return hs;
}

void drawMuliPlot(TString xtitle, Int_t isLog, Float_t y1min, Float_t y1max, Float_t y2min, Float_t y2max, THStack *hs, TCanvas *cc, TLegend *leg, TFile * data, TFile* H200, TFile* H400, TFile* ZllG, Float_t lumi)
{

  //Find the samm objects as in hs stack
  TString name = Form("%s", hs->GetHists()->First()->GetName());
  //data -> ls();
  cout<<"drawMultiPlot:: " <<name<<endl;
  TH1 *h_data = (TH1*)data->Get(  name.Data()  )->Clone();
  TH1 *h_H200 = (TH1*)H200->Get(  name.Data()  )->Clone();
  TH1 *h_H400 = (TH1*)H400->Get(  name.Data()  )->Clone();
  TH1 *h_ZllG = (TH1*)ZllG->Get(  name.Data()  )->Clone();
  // cout<<name<<endl;
  h_H200 ->Scale(lumi/1000.);
  h_H400 ->Scale(lumi/1000.);
  h_ZllG ->Scale(lumi/1000.);

  h_data -> Print();

   cc ->cd();
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  cc -> cd();
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1 -> cd();
  pad1 -> SetLogy(isLog);

  hs -> Draw("hist");
  hs -> SetTitle(";; Events");
  hs -> SetMinimum(y1min);
  hs -> SetMaximum(y1max);
  h_data  -> Draw("same e1pl");
  h_H400 -> Draw("same hist"); 
  h_H200 -> Draw("same hist");
  h_ZllG   -> Draw("same hist");
  leg -> Draw();

  cc ->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetBottomMargin(0.25);
  pad2->SetTopMargin(0);
  pad2->Draw();
  pad2 -> cd();
  
  TH1F* h1 = (TH1F*) h_data->Clone("");
  TH1F* h2 = (TH1F*) hs->Sum();
  h1 -> Divide(h2);
  h1 -> SetTitle(Form(";%s; Data/MC",xtitle.Data()));
  // if(metType=="met3") h1 -> SetTitle("; projMET/q_{T}; Data/MC");
  //if(metType=="met4") h1 -> SetTitle("; PU corr projMET/q_{T}; Data/MC");
  h1 -> SetMaximum(y2max);
  h1 -> SetMinimum(y2min);
  h1->GetYaxis()->SetNdivisions(206);
  h1->GetYaxis()->SetTitleOffset(0.4);
  h1->SetTitleSize(0.1,"XYZ");
  h1 ->SetLabelSize(0.1,"XY");
  h1 -> Draw("ep");
 
  // cc -> SaveAs(Form("%s.png", pic.Data()));
}

TH1 *THStack::Sum()
{
  //cout<<"doing sum  " <<endl;
  TList * mylist = (TList*)this->GetHists();
  TIter next(mylist);
  TH1 *hh  = (TH1*) mylist -> First() ->Clone();
  hh -> SetLineColor(kBlack);
  hh -> SetFillStyle(0);
  TObject *obj; 
  while ((obj = next()))
    {
      // cout<<obj->GetName()<<endl;
      //skip first object since it's used by creating the histogram
      if(obj == mylist->First()) continue;
      hh -> Add((TH1*)obj);     
    }
  //cout<<"end of sum"<<endl;
  return hh;
}




void *THStack::PrintYields(Int_t sel, Int_t ver, string option="tex")
{
  Int_t bins = 0;
  ofstream oo;
  if (sel==1)  oo.open(Form("yields_muons_%i.txt",ver), ofstream::out);
  if (sel==2)  oo.open(Form("yields_electrons_%i.txt",ver), ofstream::out);
  oo.precision(3); oo.setf(ios::fixed, ios::floatfield);

  TList * mylist = (TList*)this->GetHists();
  TIter next(mylist);
  TH1 *hh  = (TH1*) mylist -> First() ->Clone();
  bins = hh -> GetNbinsX();
  if(option=="tex")    oo<<hh->Integral(0,bins+1)<<"\t& ";
  if(option=="twiki")  oo<<hh->Integral(0,bins+1)<<"\t& ";

  TObject *obj; 
  while ((obj = next()))
    {
      // cout<<obj->GetName()<<endl;
      if(obj == mylist->First()) continue;
      hh = (TH1*)obj;     
      if(option=="tex")    oo<<hh->Integral(0,bins+1)<<"\t& ";
      if(option=="twiki")  oo<<hh->Integral(0,bins+1)<<"\t& ";

    }
  if(option=="tex")    oo<<"\t \\\\ \n \\hline"<<endl;
  if(option=="twiki")  oo<<"\t \\\\ \n \\hline"<<endl;

  oo<<"end of printing yields"<<endl;
  oo.close();
}


 
void AddPhoton(TH1 *ph, THStack *st, Float_t lumi){
  ph -> SetLineColor(kRed);
  ph -> SetLineWidth(2);
  ph -> Scale(lumi/191.);
  st  -> Add(ph);
}

optimalCuts calculateOptimalCuts(Int_t opt = 1, TH1 *bkg, TH1 *sig) {
  //Option values opt can be 1,2,3
  // opt = 1 - optimize S/B
  // opt = 2 - optimize S/sqrt(B)
  // opt = 3 - optimize S/sqrt(S+B)

  cout<<  bkg -> GetName()<<endl;

  Int_t sig_nBins = sig->GetNbinsX();
  Int_t bkg_nBins = bkg->GetNbinsX();
  Int_t nBins=0;
  if (sig_nBins!=bkg_nBins) cout<<" IN calculateOptimalCuts \n WARNING:  different number of bis"<<endl;
  else nBins = sig_nBins;
  
  Float_t iSig = sig -> Integral(0,nBins+1);  //Count overflows
  Float_t iBkg = bkg -> Integral(0,nBins+1);
  Float_t SB      = 0.1*iSig/iBkg;        //Higgs signal is 10*real in the histograms
  Float_t SrootB  = 0.1*iSig/sqrt(iBkg);
  Float_t SrootSB = 0.1*iSig/sqrt(iBkg + 0.1*iSig);

  optimalCuts answer;
  answer.SB      = SB;
  answer.SrootB  = SrootB;
  answer.SrootSB = SrootSB;
  answer.bin1 = 0;
  answer.bin2 = nBins+1;
  cout<<"Initial values of S/B = "<<answer.SB<<"   S/sqrt(B) = "<<SrootB<<"   S/sqrt(B+S) = "<<SrootSB<<endl;

  for(Int_t a=0; a<nBins;a++)
    for(Int_t b=nBins+1; b>a+5; b--)
      {      
	iSig = sig -> Integral(a,b);
	iBkg = bkg -> Integral(a,b);
	SB      = 0.1*iSig/iBkg;        //Higgs signal is 10*real in the histograms
	SrootB  = 0.1*iSig/sqrt(iBkg);
	SrootSB = 0.1*iSig/sqrt(iBkg + 0.1*iSig);
	if(opt==1)
	  if (SB>answer.SB && SB<100){ //<100 to prevent infinity (zero BG)
	    answer.SB = SB;
	    answer.SrootB = SrootB;
	    answer.SrootSB = SrootSB;
	    answer.bin1 = a;
	    answer.bin2 = b;
	  }
	if(opt==2)
	  if (SrootB>answer.SrootB && SB<100){
	    answer.SB = SB;
	    answer.SrootB = SrootB;
	    answer.SrootSB = SrootSB;
	    answer.bin1 = a;
	    answer.bin2 = b;
	  }
	if(opt==3)
	  if (SrootSB>answer.SrootSB && SB<100){
	    answer.SB = SB;
	    answer.SrootB = SrootB;
	    answer.SrootSB = SrootSB;
	    answer.bin1 = a;
	    answer.bin2 = b;
	  }
	

      }
  
  return answer;
}
