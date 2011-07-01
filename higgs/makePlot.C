#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
//#include "Riostream.h"

TList *FileList;
TFile *Target;

#define nC 10
#define F0 5
#define F1 6
#define F2 9
#define F3 4


void makePlot(Int_t sel=1, TString hPath="00")
{
  gROOT->ProcessLine(".L ../data/tdrstyle.C");
  setTDRStyle();

  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLegendBorderSize(0);
  //gROOT->ForceStyle();
  cout.precision(7); cout.setf(ios::fixed, ios::floatfield);
  TH1::SetDefaultSumw2(kTRUE);

  TString histoPath = hPath.Data();
  cout<<"histoPath:  "<<histoPath.Data()<<endl;
  //  Int_t sel = 1; //1-muon, 2 -electron

  //Types of met: met - pfMet, met1 - type1 corrected, met2 - pfMet passed Noise filters, 
  //met3 - projMet, met4 - puProj corrected met (those two are passed Noise filters) 
  TString metType("met3");   TString mtType("mt2"); 


  if (sel==1)  TString imgpath("~/afs/public_html/higgs/overview/muon/");  
  if (sel==2)  TString imgpath("~/afs/public_html/higgs/overview/electron/");  
  //TString imgpath(Form("~/afs/public_html/higgs/%s/",sample.Data()));  

  TFile* fda_2011A_DoubleMu  = new TFile(Form("./%s/dir_1_2011A_May10_DoubleMu_/hhhh.root",histoPath.Data()));
  TFile* fda_2011A_DoubleEl  = new TFile(Form("./%s/dir_2_2011A_May10_DoubleElectron_/hhhh.root",histoPath.Data()));
  if (sel==1) TFile  *fData = (TFile*)fda_2011A_DoubleMu;
  if (sel==2) TFile  *fData = (TFile*)fda_2011A_DoubleEl;
 

  TFile *fA_Zj = new TFile(Form("./forAnton_Zjets_%i.root",sel), "RECREATE");
  TFile *fA_Da = new TFile(Form("./forAnton_Data_%i.root",sel), "RECREATE");
  TFile *fA_Ds = new TFile(Form("./forAnton_Data_sbtr_%i.root",sel), "RECREATE");

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

  THStack *hs_met_et[nC], *hs_met2_et[nC], *hs_met_over_qt[nC], *hs_di_qt[nC], *hs_met_et_ovQt[nC], *hs_mt[nC], *hs_mtZ[nC];
  THStack *hs_jet_N[nC], *hs_jet_b_N[nC], hs_jet_b_pt[nC];
  THStack *hs_di_mass[nC];
  THStack *hs_met_dPhiLeadJet1[nC], *hs_met_dPhiLeadJet2[nC], *hs_met_dPhiClosJet1[nC], *hs_met_dPhiClosJet2[nC];
  for(Int_t n = 0; n<nC; n++)
    {
      hs_met_et[n]       = makeStack(list_bg, Form("%s_et_%i", metType.Data(), n));
      hs_met2_et[n]      = makeStack(list_bg, Form("met2_et_%i", n));
      hs_mt[n]           = makeStack(list_bg, Form("%s_%i", mtType.Data(), n));
      hs_mtZ[n]          = makeStack(list_bg, Form("mtZ_%i", n));

      hs_met_over_qt[n]  = makeStack(list_bg, Form("%s_over_qt_%i", metType.Data(), n));
      hs_met_et_ovQt[n]  = makeStack(list_bg, Form("%s_et_ovQt_%i", metType.Data(), n));

      hs_di_qt[n]   = makeStack(list_bg, Form("di_qt_%i",n));
      hs_di_mass[n] = makeStack(list_bg, Form("di_mass_%i",n));
      hs_jet_N[n]   = makeStack(list_bg, Form("jet_N_%i",n));
      hs_jet_b_N[n] = makeStack(list_bg, Form("jet_b_N_%i",n)); 
    
      hs_met_dPhiLeadJet1[n] = makeStack(list_bg, Form("met2_dPhiLeadJet1_%i",n));
      hs_met_dPhiLeadJet2[n] = makeStack(list_bg, Form("met2_dPhiLeadJet2_%i",n));
      hs_met_dPhiClosJet1[n] = makeStack(list_bg, Form("met2_dPhiClosJet1_%i",n));
      hs_met_dPhiClosJet2[n] = makeStack(list_bg, Form("met2_dPhiClosJet2_%i",n));
    }
  //  hs_jet_b_pt[8]  = makeStack(list_bg, Form("jet_b_pt_%i",8)); //_ to be fixed

  /*
  fA_Ds -> cd();
  TH1F *Zjets_qt   = fData->Get("di_qt_4") ->Clone();
  TH1F *Zjets_mass = fData->Get("di_mass_4") ->Clone();
  TH1F *Zjets_met  = fData->Get("met3_et_5") ->Clone();
  TH1F *Zjets_met2 = fData->Get("met2_et_5") ->Clone();
  TH1F *Zjets_mt   = fData->Get("mt2_6") ->Clone();
  TH1F *Zjets_nj5  = fData->Get("jet_N_5") ->Clone();
  TH1F *Zjets_nj6  = fData->Get("jet_N_6") ->Clone();
  Zjets_qt   -> Write();
  Zjets_mass -> Write();
  Zjets_met  -> Write();
  Zjets_met2 -> Write();
  Zjets_mt   -> Write();
  Zjets_nj5  -> Write();
  Zjets_nj6  -> Write();
  fA_Ds -> Close();
  */

 
  fA_Zj -> cd();
  TH1F *Zjets_qt   = fmc_Zjets->Get("di_qt_4") ->Clone();
  TH1F *Zjets_mass = fmc_Zjets->Get("di_mass_4") ->Clone();
  TH1F *Zjets_met  = fmc_Zjets->Get("met3_et_5") ->Clone();
  TH1F *Zjets_met2 = fmc_Zjets->Get("met2_et_5") ->Clone();
  TH1F *Zjets_mt   = fmc_Zjets->Get("mt2_6") ->Clone();
  TH1F *Zjets_nj5  = fmc_Zjets->Get("jet_N_5") ->Clone();
  TH1F *Zjets_nj6  = fmc_Zjets->Get("jet_N_6") ->Clone();
  TH1F *Zjets_mOq5 = fmc_Zjets->Get("met2_over_qt_5") ->Clone();
  TH1F *Zjets_mOq6 = fmc_Zjets->Get("met2_over_qt_6") ->Clone();
  Zjets_qt   -> Write();
  Zjets_mass -> Write();
  Zjets_met  -> Write();
  Zjets_met2 -> Write();
  Zjets_mt   -> Write();
  Zjets_nj5  -> Write();
  Zjets_nj6  -> Write();
  Zjets_mOq5 -> Write();
  Zjets_mOq6 -> Write();
  fA_Zj -> Close();

  
  fA_Da -> cd();
  TH1F *da__qt   = fData->Get("di_qt_4") ->Clone();
  TH1F *da__mass = fData->Get("di_mass_4") ->Clone();
  TH1F *da__met  = fData->Get("met3_et_5") ->Clone();
  TH1F *da__met2 = fData->Get("met2_et_5") ->Clone();
  TH1F *da__mt   = fData->Get("mt2_6") ->Clone();
  TH1F *da__nj5  = fData->Get("jet_N_5") ->Clone();
  TH1F *da__nj6  = fData->Get("jet_N_6") ->Clone();
  TH1F *da__mOq5 = fData->Get("met2_over_qt_5") ->Clone();
  TH1F *da__mOq6 = fData->Get("met2_over_qt_6") ->Clone();
  da__qt   -> Write();
  da__mass -> Write();
  da__met  -> Write();
  da__met2 -> Write();
  da__mt   -> Write();
  da__nj5  -> Write();
  da__nj6  -> Write();
  da__mOq5 -> Write();
  da__mOq6 -> Write();
  fA_Da -> Close();
  



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
  leg01->AddEntry(forLegend[16], "tt#rightarrow 2l2#nu2b","f");
  leg01->AddEntry(forLegend[19], "Z#gamma#rightarrow ll#gamma","f");
  leg01->AddEntry(forLegend[17], "10xH200","f");
  leg01->AddEntry(forLegend[17]," ","");  //empty slot
  leg01->AddEntry(forLegend[18], "10xH400","f");
  leg01->SetFillColor(kWhite);
  leg01->SetFillColor(kWhite);
   
  hs_met_et[F0] -> Draw("hist");
  c1 -> SaveAs(imgpath+"ov01.png");

  TCanvas *c2 = new TCanvas("c2","example",600,700);
 
  drawMuliPlot("projMET", 1, 0.0001, 100000, 0,5, hs_met_et[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov01.png");
  drawMuliPlot("projMET/q_{T}", 1, 0.0001, 100000, 0,3, hs_met_over_qt[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov02.png");

  drawMuliPlot("projMET", 1, 0.0001, 100000, 0,5, hs_met_et[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov03.png");
  drawMuliPlot("projMET/q_{T}", 1, 0.0001, 100000, 0,5, hs_met_over_qt[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov04.png");

  drawMuliPlot("projMET", 1, 0.0001, 100000, 0,3, hs_met_et[F2], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov05.png");
  drawMuliPlot("projMET/q_{T}", 1, 0.0001, 100000, 0,3, hs_met_over_qt[F2], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov06.png");
  
  drawMuliPlot("M(ll)", 1, 0.0001, 100000, 0,2, hs_di_mass[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov07.png");
  drawMuliPlot("MT", 1, 0.0001, 100000, 0,2, hs_mt[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov08.png");

  drawMuliPlot("M(ll)", 1, 0.0001, 100000, 0,5, hs_di_mass[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov09.png");
  drawMuliPlot("MT", 1, 0.0001, 100000, 0,5, hs_mt[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov10.png");

  drawMuliPlot("M(ll)", 1, 0.0001, 1000000, 0,2, hs_di_mass[F2], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov11.png");
  drawMuliPlot("MT", 1, 0.0001, 100000, 0,2, hs_mt[F2], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov12.png");

  drawMuliPlot("N jets", 1, 0.0001, 100000, 0,5, hs_jet_N[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov13.png");
  drawMuliPlot("N b-jets", 1, 0.0001, 100000, 0,5, hs_jet_b_N[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov14.png");

  drawMuliPlot("N jets", 1, 0.0001, 100000, 0,5, hs_jet_N[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov15.png");
  drawMuliPlot("N b-jets", 1, 0.0001, 100000, 0,5, hs_jet_b_N[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov16.png");


  drawMuliPlot("N jets", 1, 0.0001, 100000, 0,5, hs_jet_N[F3], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov17.png");
  drawMuliPlot("N b-jets", 1, 0.0001, 100000, 0,5, hs_jet_b_N[F3], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov18.png");
  

 // drawMuliPlot("N jets", 1, 0.0001, 100000, 0,5, hs_jet_N[F2], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
 // c2 -> SaveAs(imgpath+"ov17.png");
 // drawMuliPlot("N b-jets", 1, 0.0001, 100000, 0,5, hs_jet_b_N[F2], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
 //   c2 -> SaveAs(imgpath+"ov18.png");
  

  drawMuliPlot("N b-jets", 1, 0.0001, 100000, 0,5, hs_jet_b_N[7], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov19.png");
  drawMuliPlot("N b-jets", 1, 0.0001, 100000, 0,5, hs_jet_b_N[8], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov20.png");

  drawMuliPlot("#Delta#phi(MET, lead jet), p_{T}>20, |#eta|<2.4", 1, 0.0001, 1000000, 0,2, hs_met_dPhiLeadJet1[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov21.png");
  drawMuliPlot("#Delta#phi(MET, closest jet), p_{T}>20, |#eta|<2.4", 1, 0.0001, 1000000, 0,2, hs_met_dPhiClosJet1[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov22.png");

  drawMuliPlot("#Delta#phi(MET, lead jet), p_{T}>20, |#eta|<2.4", 1, 0.0001, 100000, 0,5, hs_met_dPhiLeadJet1[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov23.png");
  drawMuliPlot("#Delta#phi(MET, closest jet), p_{T}>20, |#eta|<2.4", 1, 0.0001, 100000, 0,5, hs_met_dPhiClosJet1[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov24.png");

  drawMuliPlot("#Delta#phi(MET, lead jet), p_{T}>20, |#eta|<2.4", 1, 0.0001, 1000000, 0,2, hs_met_dPhiLeadJet1[F2], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov25.png");
  drawMuliPlot("#Delta#phi(MET, closest jet), p_{T}>20, |#eta|<2.4", 1, 0.0001, 1000000, 0,2, hs_met_dPhiClosJet1[F2], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov26.png");


  drawMuliPlot("q_{T} (di-lepton p_{T})", 1, 0.0001, 100000, 0,5, hs_di_qt[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov27.png");
  drawMuliPlot("q_{T} (di_lepton p_{T})", 1, 0.0001, 100000, 0,5, hs_di_qt[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov28.png");

  drawMuliPlot("pfMET", 1, 0.0001, 100000, 0,5, hs_met2_et[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov29.png");
  drawMuliPlot("MT (using Z mass)", 1, 0.0001, 100000, 0,5, hs_mtZ[F0], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov30.png");

  drawMuliPlot("pfMET", 1, 0.0001, 100000, 0,5, hs_met2_et[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov31.png");
  drawMuliPlot("MT (using Z mass)", 1, 0.0001, 100000, 0,5, hs_mtZ[F1], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  c2 -> SaveAs(imgpath+"ov32.png");
  

  // drawMuliPlot("pfMET", 1, 0.0001, 100000, 0,5, hs_jet_b_pt[8], c2, leg01, fData, fmc_ggH200, fmc_ggH400, fmc_ZllG);
  //c2 -> SaveAs(imgpath+"ov33.png");

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

THStack* makeStack(TList *sourcelist, TString name)
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
	    hs -> Add((TH1*)key->ReadObj()->Clone());
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
		   hs -> Add((TH1*)key2->ReadObj()->Clone());
		}
	      //delete h2;
	    }
	  
	  nextsource = (TFile*)sourcelist->After( nextsource );
	}
      
    }



  return hs;
}

void drawMuliPlot(TString xtitle, Int_t isLog, Float_t y1min, Float_t y1max, Float_t y2min, Float_t y2max, THStack *hs, TCanvas *cc, TLegend *leg, TFile * data, TFile* H200, TFile* H400, TFile* ZllG)
{

  //Find the samm objects as in hs stack
  TString name = Form("%s", hs->GetHists()->First()->GetName());
  //data -> ls();
  cout<<"drawMultiPlot:: " <<name<<endl;
  TH1 *h_data = (TH1F*)data->Get(  name.Data()  );
  TH1F *h_H200 = (TH1F*)H200->Get(  name.Data()  );
  TH1F *h_H400 = (TH1F*)H400->Get(  name.Data()  );
  TH1 *h_ZllG = (TH1F*)ZllG->Get(  name.Data()  );
  // cout<<name<<endl;

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
