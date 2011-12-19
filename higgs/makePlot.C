#define nC 20
#define F0 6
#define F1 16
#define F2 17
#define F3 4

void makePlot(Int_t sel=1, TString hPath="00")
{
  gROOT->ProcessLine(".L ./utils.C");
  gROOT->ProcessLine(".L ../data/tdrstyle.C");
  setTDRStyle();

  gStyle->SetPadGridX(1); gStyle->SetPadGridY(1);
  gStyle->SetHistLineWidth(2);  gStyle->SetLegendBorderSize(0);
  cout.precision(3); cout.setf(ios::fixed, ios::floatfield);
  TH1::SetDefaultSumw2(kTRUE);

  Float_t intLumi = 1;
  if(sel==1)  intLumi = 215.1 + 927.6 + 370.9 + 663.0 + 2511; //double mu	
  if(sel==2)  intLumi = 215.1 + 789.2 + 313.2 + 662.2 + 2511; //double ele	
  TString dir("test");

  Bool_t doTest    = 0;
  Bool_t doPhotons = 1;
  Bool_t doOverview= 1;

  TString ssel("none"), gsel("none");
  if (sel==1)  {
    ssel = "muon";  gsel ="muGamma";
    if(doPhotons)  TString imgpath(Form("~/afs/public_html/higgs/%s/%s/", dir.Data(), gsel.Data() ) ); 
    else           TString imgpath(Form("~/afs/public_html/higgs/%s/%s/", dir.Data(), ssel.Data() ) ); 
  }
  if (sel==2)  {
    ssel = "electron";  gsel ="eGamma";
    if(doPhotons)  TString imgpath(Form("~/afs/public_html/higgs/%s/%s/", dir.Data(), gsel.Data() ) ); 
    else           TString imgpath(Form("~/afs/public_html/higgs/%s/%s/", dir.Data(), ssel.Data() ) ); 
  }

  TString histoPath = Form("%s/%s", hPath.Data(), ssel.Data());
  cout<<"histoPath:  "<<histoPath.Data()<<"  int Lumi: "<<intLumi<<endl;

  TFile  *fda_Data  = new TFile(Form("./%s/m_Data_%i.root", hPath.Data(),  sel));  //Merged Data
  TFile  *fData = (TFile*)fda_Data;
  
  // TFile* fmc_Wjets     = new TFile(Form("./%s/hhhh_Wjets.root",histoPath.Data() ));

  TFile* fmc_ZZ        = new TFile(Form("./%s/hhhh_ZZ_1.root",histoPath.Data() ));
  TFile* fmc_WW        = new TFile(Form("./%s/hhhh_WW_1.root",histoPath.Data() ));
  TFile* fmc_WZ        = new TFile(Form("./%s/hhhh_WZ_1.root",histoPath.Data() ));

  TFile *fmc_ggH[10], *fmc_vbfH[10];
  for(Int_t i=1; i<9; i++){
    Int_t  m = 200+i*50;
    fmc_ggH[i]   = new TFile(Form("./%s/hhhh_ggHZZ%i_1.root",histoPath.Data(),m ));
    fmc_vbfH[i]  = new TFile(Form("./%s/hhhh_VBFHZZ%i_1.root",histoPath.Data(),m ));
  }

  TFile* fmc_tt     = new TFile(Form("./%s/m_ttbar_%i.root", hPath.Data(), sel ));
  TFile* fmc_Top    = new TFile(Form("./%s/m_Top_%i.root", hPath.Data(), sel ));

  TFile* fmc_Zjets  = new TFile(Form("./%s/m_Zjets_%i.root", hPath.Data(), sel ));
  //TFile* fmc_Zjets  = new TFile(Form("./%s/m_libZjets_%i.root", hPath.Data(), sel ));
  if(doPhotons)  fmc_Zjets  = new TFile(Form("./%s/m_DataPh_%i.root", hPath.Data(), sel ));

  /*
  RescaleToLumiAndColors(fmc_tt, 1000,1000, kMagenta+1, kBlue-3, 1001, 0);
  RescaleToLumiAndColors(fmc_Top,   1000,1000, kOrange+9, kOrange+6,1001, 0);
  RescaleToLumiAndColors(fmc_Zjets, 1000,1000, kRed+2, kRed+1,3004, 0);
  */    
 
  //List of background samples to Stack
  list_bg = new TList();
  list_bg->Add(fmc_Top);
  list_bg->Add(fmc_tt);
  list_bg->Add(fmc_WW);
  list_bg->Add(fmc_WZ);
  list_bg->Add(fmc_ZZ);
  list_bg->Add(fmc_Zjets);
  //list_bg->Add(fmc_Wjets);

  list_bg2 = new TList();
  list_bg2->Add(fmc_Top);
  list_bg2->Add(fmc_tt);
  list_bg2->Add(fmc_WW);
  list_bg2->Add(fmc_WZ);
  list_bg2->Add(fmc_ZZ);
  //list_bg2->Add(fmc_Wjets);

  list_overlay = new TList();
  list_overlay->Add(fData); 
  list_overlay->Add(fmc_ggH[1]);//250
  //list_overlay->Add(fmc_ggH[4]);//400

  list_ggH = new TList();
  for(Int_t i=1; i<9; i++)  
    list_ggH->Add(fmc_ggH[i]);

  list_vbfH = new TList();
  for(Int_t i=1; i<9; i++)  
    list_vbfH->Add(fmc_vbfH[i]);
  
  PrintYields(list_bg, list_ggH, list_vbfH, fData, intLumi, histoPath, "twiki");
  PrintYields(list_bg, list_ggH, list_vbfH, fData, intLumi, histoPath, "tex");
  
  
  TH1 * forLegend[21];
  list_legend = new TList();
  list_legend -> AddAll(list_bg);
  list_legend -> AddAll(list_overlay);

  TFile *ff[10]; 
  Int_t size = list_legend->GetSize();
  TH1 *hh[10]; 
  if(size>10) {cout<<"To many plots to overlay"<<endl; return 0;}
  for(Int_t n=0; n<size; n++)
    {
      ff[n] = (TFile*)list_legend->At(n);
      hh[n] = (TH1*) ff[n]->Get("Andrey/met0_et_0")->Clone();
      cout<<n<<"   "<<ff[n] -> GetName()<<endl;
      forLegend[n] = (TH1F*)hh[n];
    }
 
  leg01 = new TLegend(0.53,0.7,0.95,0.90);
  leg01 -> SetNColumns(2);
  leg01 -> SetTextSize(0.04);
  leg01->AddEntry(forLegend[6], "Data","epl");
  if(doPhotons) leg01->AddEntry(forLegend[5],  "Z+jets (#gamma data)","f");
  else          leg01->AddEntry(forLegend[5],  "Z + jets","f");
  leg01->AddEntry(forLegend[4],  "ZZ","f");
  //leg01->AddEntry(forLegend[10],  "W + jets","f");
  leg01->AddEntry(forLegend[3],  "WZ","f");
  leg01->AddEntry(forLegend[0],  "tW","f");
  leg01->AddEntry(forLegend[2],  "WW","f");
  leg01->AddEntry(forLegend[1], "ttbar","f");  
  leg01->AddEntry(forLegend[7], "H250","l");
  //leg01->AddEntry(forLegend[8], "H400","f");
  //leg01->AddEntry(forLegend[9], "H400","f");
  
  leg01->SetFillColor(kWhite);
  
  hh[0] -> Draw("hist");
  c1 -> SaveAs(imgpath+"ov01.png");
  TCanvas *c2 = new TCanvas("c2","example",600,700);
  
  
  if(doTest){
    TString testpath("~/afs/public_html/test/");  
    

    drawMuliPlot("","di qT", Form("di_qt_%i", F0), 1, 0.001, 1e6, 0,3, c2, leg01, list_overlay, list_bg, intLumi, sel); 
    c2 -> SaveAs(testpath+"p01.png");    

    c1->cd();

    TH1F * test1 = fmc_ggH[1]->Get("Andrey/jet_N_17")->Clone();
    TH1F * test2 = fmc_ggH[2]->Get("Andrey/jet_N_17")->Clone();

    test1->Scale(intLumi/1000./10);
    test2->Scale(intLumi/1000./10);
    test1->SetMaximum(25);
    test1->Draw("hist");
    test2->Draw("same hist");   
    c1 -> SaveAs(testpath+"p03.png");

    //TH1F * test1 = fmc_ggH[1]->Get("Andrey/mt2_17")->Clone();
    //TH1F * test2 = fmc_ggH[3]->Get("Andrey/mt2_17")->Clone();

    //test1->Draw();   
    //test2->Draw("same");   
    //c1 -> SaveAs(testpath+"p04.png");

  }
  
  if(doOverview){
      //-------------------------//
      //------- Met plots--------//
      //-------------------------//  
      drawMuliPlot("","pfMET", Form("met1_et_%i", F0), 1, 0.1, 1e6, 0,2.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Met/m01.png");

      drawMuliPlot("","pfMET Phi", Form("met1_phi_%i", F0), 1, 0.1, 1e6, 0,2.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Met/m02.png");

      drawMuliPlot("","MT", Form("mt2_%i", F0), 1, 0.1, 1e6, 0,1.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Met/m03.png");

      drawMuliPlot("","type 1 pfMET", Form("met2_et_%i", F0), 1, 0.1, 1e6, 0,2.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Met/m04.png");

      drawMuliPlot("","pfMET long to di-lepton", Form("met1_lg_%i", F0), 1, 0.1, 1e7, 0,2.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Met/m05.png");

      drawMuliPlot("","Long Recoil = -(pfMet+Z/G)", Form("met1_recoil_lg_%i", F0), 1, 0.1, 1e7, 0,2.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Met/m06.png");


      drawMuliPlot("","puCorrMET", Form("met3_et_%i", F0), 1, 0.1, 1e6, 0,2.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Met/m07.png");

      drawMuliPlot("","projMET", Form("met4_et_%i", F0), 1, 0.1, 1e6, 0,2.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Met/m08.png");
      drawMuliPlot("","ZprojMET", Form("met5_et_%i", F0), 1, 0.1, 1e6, 0,2.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Met/m09.png");

      drawMuliPlot("","redMET1", Form("met6_et_%i", F0), 1, 0.1, 1e6, 0,2.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Met/m10.png");
      drawMuliPlot("","redMET2", Form("met7_et_%i", F0), 1, 0.1, 1e6, 0,2.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Met/m11.png");

      //-------------------------//
      //------- Jet plots--------//
      //-------------------------//  

      drawMuliPlot("","N jets", Form("jet_N_%i", F0), 1, 0.1, 1e6, 0,2.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Jet/j01.png");
      drawMuliPlot("","#Delta#phi(MET, clos jet), p_{T}>30, |#eta|<4.8", Form("met1_dPhiClosJet1_%i", F0), 1, 0.1, 1e7, 0,2, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Jet/j02.png");

      drawMuliPlot("Post b-veto","N b-jets", Form("jet_b_N_%i", F0), 1, 0.1, 1e6, 0,1.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Jet/j03.png");
      drawMuliPlot("No b-veto, MET>70","N b-jets", Form("jet_b_N_%i", F1), 0, 0.1, 350, 0,1.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Jet/j04.png");

      drawMuliPlot("","pt of leading jet", Form("jet_pt_%i", F0), 1, 0.1, 1e6, 0,2.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Jet/j05.png");

      drawMuliPlot("","eta of leading jet", Form("jet_eta_%i", F0), 1, 0.1, 1e6, 0,2.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Jet/j06.png");

      drawMuliPlot("","dR(jet,lep1)", Form("jet_dRlep1_%i", F0), 1, 0.1, 1e6, 0,4.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Jet/j07.png");

      drawMuliPlot("","dR(jet,lep2)", Form("jet_dRlep2_%i", F0), 1, 0.1, 1e6, 0,4.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Jet/j08.png");

      drawMuliPlot("MET>40, b-veto", "N  jets", Form("jet_N_%i", F2), 0, 0.1, 3000, 0,1.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"Jet/j09.png");
     
      //-------------------------------//
      //------- Di-Lepton plots--------//
      //-------------------------------//  
      drawMuliPlot("","M(ll)", Form("di_mass_%i", F0), 1, 0.1, 1e6, 0,1.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"diLepton/di01.png");
      drawMuliPlot("","Leptons in Barrel, |#eta|<1.444,  M(ll)", Form("di_mass_EB_%i", F0), 1, 0.1, 1e6, 0,1.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"diLepton/di02.png");
      drawMuliPlot("","Leptons in Endcap, |#eta|>1.566, M(ll)", Form("di_mass_EE_%i", F0), 1, 0.1, 1e6, 0,1.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"diLepton/di03.png");
      drawMuliPlot("","Leptons in EB/EE, mixed, M(ll)", Form("di_mass_EX_%i", F0), 1, 0.1, 1e6, 0,1.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"diLepton/di04.png");

      drawMuliPlot("","q_{T} (di-lepton p_{T})", Form("di_qt_%i", F0), 1, 0.1, 1e6, 0,2.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"diLepton/di05.png");

      drawMuliPlot("","Di-lepton Eta", Form("di_eta_%i", F0), 1, 0.1, 1e6, 0.5,1.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"diLepton/di06.png");

      drawMuliPlot("","dPhi(Di-lep, Met)", Form("di_dPhiMet_%i", F0), 1, 0.1, 1e6, 0.5,1.9, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"diLepton/di07.png");

      //----------------------------//
      //------- Lepton plots--------//
      //----------------------------//  
      drawMuliPlot("","Leading Lepton eta", Form("l1_eta_%i", F0), 1, 0.1, 1e6, 0.5,1.4, c2, leg01, list_overlay, list_bg, intLumi, sel); 
      c2 -> SaveAs(imgpath+"Lepton/l01.png");
      drawMuliPlot("","Trailing Lepton eta", Form("l2_eta_%i", F0), 1, 0.1, 1e6, 0.5,1.4, c2, leg01, list_overlay, list_bg, intLumi, sel); 
      c2 -> SaveAs(imgpath+"Lepton/l02.png");

      drawMuliPlot("","Leading Lepton phi", Form("l1_phi_%i", F0), 1, 0.1, 1e6, 0.5,1.4, c2, leg01, list_overlay, list_bg, intLumi, sel); 
      c2 -> SaveAs(imgpath+"Lepton/l03.png");
      drawMuliPlot("","Trailing Lepton phi", Form("l2_phi_%i", F0), 1, 0.1, 1e6, 0.5,1.4, c2, leg01, list_overlay, list_bg, intLumi, sel); 
      c2 -> SaveAs(imgpath+"Lepton/l04.png");

      drawMuliPlot("","Leading Lepton pt", Form("l1_pt_%i", F0), 1, 0.1, 1e6, 0.5,1.4, c2, leg01, list_overlay, list_bg, intLumi, sel); 
      c2 -> SaveAs(imgpath+"Lepton/l05.png");
      drawMuliPlot("","Trailing Lepton pt", Form("l2_pt_%i", F0), 1, 0.1, 1e6, 0.5,1.4, c2, leg01, list_overlay, list_bg, intLumi, sel); 
      c2 -> SaveAs(imgpath+"Lepton/l06.png");

      drawMuliPlot("","dR(lep1,lep2)", Form("l0_dR_%i", F0), 1, 0.1, 1e6, 0.5,1.4, c2, leg01, list_overlay, list_bg, intLumi, sel); 
      c2 -> SaveAs(imgpath+"Lepton/l07.png");
      drawMuliPlot("","dPhi(lep1,lep2)", Form("l0_dPhi_%i", F0), 1, 0.1, 1e6, 0.5,1.4, c2, leg01, list_overlay, list_bg, intLumi, sel); 
      c2 -> SaveAs(imgpath+"Lepton/l08.png");
      drawMuliPlot("","dEta(lep1,lep2)", Form("l0_dEta_%i", F0), 1, 0.1, 1e6, 0.5,1.4, c2, leg01, list_overlay, list_bg, intLumi, sel); 
      c2 -> SaveAs(imgpath+"Lepton/l09.png");
      drawMuliPlot("","pt_{lep2}/pt_{lep1}", Form("l0_ptRatio_%i", F0), 1, 0.1, 1e6, 0.5,1.4, c2, leg01, list_overlay, list_bg, intLumi, sel); 
      c2 -> SaveAs(imgpath+"Lepton/l10.png");

      //----------------------------//
      //---Misc plots: vtx etc -----//
      //----------------------------//  
      drawMuliPlot("","evts cut by cut", "evt_byCut", 1, 0.1, 1e6, 0,2.9, c2, leg01, list_overlay, list_bg, intLumi, sel); 
      c2 -> SaveAs(imgpath+"Misc/mis01.png");
      drawMuliPlot("","nVtx total", Form("vtx_nPV_tot_%i", F0), 1, 0.1, 1e6, 0,1.9, c2, leg01, list_overlay, list_bg, intLumi, sel); 
      c2 -> SaveAs(imgpath+"Misc/mis02.png");

      drawMuliPlot("","nVtx raw", Form("vtx_nPV_raw_%i", F0), 1, 0.1, 1e6, 0,1.9, c2, leg01, list_overlay, list_bg, intLumi, sel); 
      c2 -> SaveAs(imgpath+"Misc/mis03.png");
      drawMuliPlot("","nVtx reweighted", Form("vtx_nPV_weight_%i", F0), 1, 0.1, 1e6, 0,1.9, c2, leg01, list_overlay, list_bg, intLumi, sel); 
      c2 -> SaveAs(imgpath+"Misc/mis04.png");


      drawMuliPlot("","vtx 1 nDof", Form("vtx_ndof_1_%i", F0), 1, 0.1, 1e6, 0.5,2.9, c2, leg01, list_overlay, list_bg, intLumi, sel); 
      c2 -> SaveAs(imgpath+"Misc/mis05.png");

      drawMuliPlot("","vtx 2 nDof", Form("vtx_ndof_2_%i", F0), 1, 0.1, 1e6, 0.5,2.9, c2, leg01, list_overlay, list_bg, intLumi, sel); 
      c2 -> SaveAs(imgpath+"Misc/mis06.png");

      //drawMuliPlot("","N photons", Form("ph_nGamma_%i", F0), 1, 0.1, 1e6, 0.5,2.9, c2, leg01, list_overlay, list_bg, intLumi, sel); 
      //c2 -> SaveAs(imgpath+"Misc/mis07.png");


      //--- Create hml page ---//
      //TString afs = "~/afs/public_html/higgs";
      //system (Form("python %s/writeIntexHTML.py %s/%s  > stdout.txt", afs.Data(), afs.Data(), dir.Data() ));
      

      //--------------------//
      //---Copy printouts---//
      //--------------------//

      //      system (Form("cp %s/events_printout_* %s/printouts/", histoPath.Data(), imgpath.Data()));
      //system (Form("cp -r %s/printout_Double* %s/printouts/", histoPath.Data(), imgpath.Data()));

    /*
    drawMuliPlot("","projMET/q_{T}", Form("_%i", F0), 1, 0.1, 1e6, 0,3, ____qt[F0], c2, leg01, list_overlay, list_bg, intLumi, sel);
    c2 -> SaveAs(imgpath+"ov02.png");
    
    drawMuliPlot("","projMET", Form("_%i", F0), 1, 0.1, 1e6, 0,5, ___[F1], c2, leg01, list_overlay, list_bg, intLumi, sel);
    c2 -> SaveAs(imgpath+"ov03.png");
    drawMuliPlot("","projMET/q_{T}", Form("_%i", F0), 1, 0.1, 1e6, 0,5, ____qt[F1], c2, leg01, list_overlay, list_bg, intLumi, sel);
    c2 -> SaveAs(imgpath+"ov04.png");
  

    drawMuliPlot("","projMET", Form("_%i", F0), 1, 0.1, 1e6, 0,3, ___[F2], c2, leg01, list_overlay, list_bg, intLumi, sel);
    c2 -> SaveAs(imgpath+"ov05.png");
    drawMuliPlot("","projMET/q_{T}", Form("_%i", F0), 1, 0.1, 1e6, 0,3, ____qt[F2], c2, leg01, list_overlay, list_bg, intLumi, sel);
    c2 -> SaveAs(imgpath+"ov06.png");
    */
    
      /*    
      
      drawMuliPlot("","#Delta#phi(MET, lead jet), p_{T}>20, |#eta|<2.4", Form("_%i", F0), 1, 0.1, 1e60, 0,2, c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"ov21.png");
      
      drawMuliPlot("","#Delta#phi(MET, lead jet), p_{T}>20, |#eta|<2.4", Form("_%i", F0), 1, 0.1, 1e6, 0,5, ___[F1], c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"ov23.png");
      drawMuliPlot("","#Delta#phi(MET, closest jet), p_{T}>20, |#eta|<2.4", Form("_%i", F0), 1, 0.1, 1e6, 0,5, ___[F1], c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"ov24.png");
      
      drawMuliPlot("","#Delta#phi(MET, lead jet), p_{T}>20, |#eta|<2.4", Form("_%i", F0), 1, 0.1, 1e60, 0,2, ___[F2], c2, leg01, list_overlay, list_bg, intLumi, sel);
      c2 -> SaveAs(imgpath+"ov25.png");
      drawMuliPlot("","#Delta#phi(MET, closest jet), p_{T}>20, |#eta|<2.4", Form("_%i", F0), 1, 0.1, 1e60, 0,2, ___[F2], c2, leg01, list_overlay, list_bg, intLumi, sel);       
      c2 -> SaveAs(imgpath+"ov26.png");
      
  */
  
  }

  cout<<"\n\n   -- end of job  --"<<endl;
  
}

/*struct optimalCuts {
  Float_t SB;
  Float_t SrootB;
  Float_t SrootSB;
  Int_t bin1;
  Int_t bin2;
  };*/

/*
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

*/


void twoScales(TGraph *g1, TGraph *g2, Float_t x1, Float_t x2, Float_t y1, Float_t y2, Float_t y3, TCanvas *c1) {
  c1 -> Clear();

    g1->SetFillColor(kBlue);
    g1->SetFillStyle(3010);
    g2->SetFillColor(kRed);
    g2->SetFillStyle(3010);


  //compute the pad range with suitable margins
  Double_t ymin = y1;
  Double_t ymax = y2;
  Double_t dy = (ymax-ymin)/0.76; //10 per cent margins top and bottom
  Double_t xmin = x1;
  Double_t xmax = x2;
  Double_t dx = (xmax-xmin)/0.76; //10 per cent margins left and right

  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  TPad *pad2 = new TPad("pad2","",0,0,1,1);
  pad2->SetFillStyle(4000); //will be transparent
  
  pad1->Draw();
  pad1->cd();
  // Margins:
  pad1->SetTopMargin(0.1);
  pad1->SetBottomMargin(0.14);
  pad1->SetLeftMargin(0.14);
  pad1->SetRightMargin(0.1);
  pad1->SetTicky(0);

  g1->Draw("AL3");
  g1-> SetMaximum(y3);
  g1-> SetMinimum(0);

  g1->GetHistogram()->SetAxisRange(xmin,xmax);
  g1 -> GetYaxis() -> SetLabelColor(kBlue);
  g1 -> GetYaxis() -> SetTitle("S/B");
  g1 -> GetYaxis() -> SetTitleOffset(1.0);
  g1 -> GetYaxis() -> SetTitleSize(0.05);
  pad1->Modified();
  c1->cd();

  pad2->Range(xmin-0.14*dx,ymin-0.14*dy,xmax+0.1*dx,ymax+0.1*dy);  //correspond to Margins in gstyle
  pad2->Draw();
  pad2->cd();
  
  g2->Draw("same L3");
   
  // draw axis on the right side of the pad
  TGaxis *axis1 = new TGaxis(x2,y1,x2,y2,y1,y2,510,"+L");
  axis1->SetLabelColor(kRed);
  axis1 -> SetTitle("S/#sqrt{S+B}");
  axis1->Draw();

  //TGaxis *axis2 = new TGaxis(x1,y1,x1,y2,y1,y2,510,"+L");
  //axis2->SetLabelColor(kRed);
  //axis2 -> SetTitle("S/#sqrt{S+B}");
  //axis2->Draw();

  pad2->Update();

  pad1 -> cd();
  leg03 = new TLegend(0.2,0.73,0.4,0.89);
  leg03 -> SetTextSize(0.04);
  leg03->AddEntry(g1, "S/B", "f");
  leg03->AddEntry(g2, "#frac{S}{#sqrt{S+B}}", "f");
  leg03->SetFillColor(kWhite);
  leg03 -> Draw();
  	

}

StoBPlot(string var= "pfMet", Float_t phiCut = 0, TFile *bgFile, Float_t sc1, TFile *sigFile, Float_t sc2,  Double_t StoB[], Double_t StoB_err[], Double_t sqrtStoB[],Double_t  sqrtStoB_err[], Double_t x_met[], Int_t startpoint){

  /*  bgFile -> cd("Andrey");
  TString ct_sample("a");

  TBranch *branch  = cutTree->GetBranch("ct_sample");
  branch->SetAddress(&ct_sample);
  cutTree -> GetEntry(0);
  cout<<"test:  "<<ct_sample<<endl;
  */

  Float_t dPhiCut = phiCut; //400 etc
  
  for(Int_t i=0; i<100; i++)
    {
      x_met[i]  = startpoint+ 1*i;

      //cout<<i<<"  "<<intLumi<<endl;
      bgFile -> cd("Andrey");
      //cutTree -> Draw(Form("ct_%s>>h1",var.c_str()),Form("ct_evtWeight*(ct_dPhiMetJet>%f && ct_MT>216 &&ct_MT<272 && ct_%s>%f)", dPhiCut, var.c_str(), x_met[i]),"hist");
      cutTree -> Draw(Form("ct_%s>>h1",var.c_str()),Form("ct_evtWeight*(ct_dPhiMetJet>%f && ct_%s>%f)", dPhiCut, var.c_str(), x_met[i]),"hist");
      //cutTree -> Draw(Form("ct_%s>>h1",var.c_str()),Form("ct_evtWeight*(ct_dPhiMetJet>%f && ct_%s>%f && ct_nJets>=2)", dPhiCut, var.c_str(), x_met[i]),"hist");
      
      TH1F * temp = (TH1F*)h1;
      temp->Scale(sc1);
      Double_t bg_err =0 ;
      //temp -> Print();
      Double_t bg =  temp->TH1::IntegralAndError(-1,2000, bg_err);
      
      sigFile -> cd("Andrey");
      //cutTree -> Draw(Form("ct_%s>>h2",var.c_str()),Form("ct_evtWeight*(ct_dPhiMetJet>%f && ct_MT>216 &&ct_MT<272 && ct_%s>%f)", dPhiCut, var.c_str(), x_met[i]),"hist");
      cutTree -> Draw(Form("ct_%s>>h2",var.c_str()),Form("ct_evtWeight*(ct_dPhiMetJet>%f && ct_%s>%f)", dPhiCut, var.c_str(), x_met[i]),"hist");
//cutTree -> Draw(Form("ct_%s>>h2",var.c_str()),Form("ct_evtWeight*(ct_dPhiMetJet>%f && ct_%s>%f && ct_nJets>=2)", dPhiCut, var.c_str(), x_met[i]),"hist");

            
      TH1F * temp = (TH1F*)h2;
      temp->Scale(sc2);
      Double_t sig_err = 0;
      //Double_t sig = 	temp -> Integral();
      Double_t sig = temp->TH1::IntegralAndError(-1,2000, sig_err);
            
      if(bg!=0) {
	StoB[i]    = sig/bg;
	StoB_err[i] = StoB[i]*sqrt( pow(sig_err/sig,2) + pow(bg_err/bg,2));
      }
      else {StoB[i]=0; StoB_err[i]=0;}
      //x_err[i] =0;
      
      if((sig+bg)  != 0) {
	sqrtStoB[i]     = sig/sqrt(sig+bg);
	sqrtStoB_err[i] =  sqrtStoB[i]*sqrt( pow(sig_err/sig,2) + 0.25* (pow(sig_err,2) + pow(bg_err,2))/pow(sig+bg,2)  );
      }
      else {sqrtStoB[i] = 0; sqrtStoB_err[i] =0;}
      
      //is this correct error propagation?	
      
      cout<<i<<"  "<<x_met[i]<<"    B="<<bg<<" S="<<sig<<"  S/B="<<StoB[i]<<"   S/sq(S+B)="<<sqrtStoB[i]<<endl;
      
    }
}

