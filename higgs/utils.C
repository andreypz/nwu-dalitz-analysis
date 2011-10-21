void RescaleToLumiAndColors( TFile *HistoFile, Float_t lumiInput=1, Float_t lumiToScale =1, Int_t lineColor, Int_t fillColor) {

  HistoFile -> cd("Andrey");
  TIter nextkey( gDirectory->GetListOfKeys() );
  //gDirectory->Print();

  TKey *key, *oldkey=0;
  while ( (key = (TKey*)nextkey()))
    {
      if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;
      TObject *obj = key->ReadObj();
      //cout<<obj->GetName()<<endl;
      
      if ( ! (obj->IsA()->InheritsFrom( TH1::Class() )) ) continue;
      TH1 *hist = (TH1*)key->ReadObj();
      hist -> Scale(lumiToScale/lumiInput);
      hist -> SetLineColor(lineColor);
      hist -> SetFillColor(fillColor);
      hist -> SetLineWidth(2);
      hist->Write("",TH1::kWriteDelete); //overwrite the object in the file
      delete obj;
	//cout<<sample1<<"  nEv: "<<nTotEv<<"   cs: "<<cs<<"   Rescaling all the histogram "<<hist->GetName()<<" by: "<<scaleFactor<<endl; 
       
  }
  
}
void PrintYields(TList *bgList, TList *ggHList, TList *vbfList, TFile *dataFile, Float_t lumi, TString path, string option="tex")
{
  Int_t bins = 0;
  ofstream oo;
  oo.open(Form("%s/yields_%s.txt",path.Data(), option.c_str()), ofstream::out);
  oo.precision(1); oo.setf(ios::fixed, ios::floatfield);

  TString beginLine("");
  TString endLine("");
  TString separator("");
  TString title1("");
  TString title2("");
  TString pmSign("");
  if(option=="tex"){
    title1 = " sel & WW & WZ   \\\\ \\hline";
    title2 = " sel & WW & WZ   \\\\ \\hline";
    beginLine = " ";
    endLine   = "\\\\ \\hline";
    separator = "\t &";
    pmSign = " $\pm$ ";
  }
  if(option=="twiki"){
    title1 = "| *cut*        | *top*  | *ttbar*  | *WZ*  | *WW* | *ZZ*  | *Zjets* | *Total bg* | *Data* |";
    title2 = "| *Higgs mass* | *top*  | *ttbar*  | *WZ*  | *WW* | *ZZ*  | *Zjets* | *Total bg* | *Data* | *ggH* | *vbf* | *S/&radic;S+B* | *S/B* |";
    beginLine = "| ";
    endLine   = "\t |";
    separator = "\t |";
    pmSign = " &plusmn; ";
  }

  string cutNames[18] = {
    "0. Total",
    "1. trigger",
    "2. Vtx,cosmic,3d lept veto",
    "3. Z mass",
    "4. soft 3d muon  veto",
    "5. qT > 25",
    "6. b-veto",
    "7.  H200",
    "8.  H250",
    "9.  H300",
    "10. H350",
    "11. H400",
    "12. H450",
    "13. H500",
    "14. H550",
    "15. H600",
    "16.",
    "17.",};

  oo<<title1<<endl;

  for(Int_t j = 0; j<=15; j++)
    {
      if(j<=6)  oo.precision(0);
      else  oo.precision(2);
      if(j==7) continue; //skip 200 mass
      if(j==8) oo<<"\n"<<title2<<endl;

      Int_t size = bgList->GetSize();
      TFile *ff[10];    TH1 *hh[10];
      if(size>10) {cout<<"Yields: to many MC samples to print-out"<<endl; return 0;}

      oo<<beginLine<<cutNames[j]<<separator;

      Float_t total_bg = 0, total_bgError =0;

      for(Int_t n=0; n<size; n++)
	{
	  ff[n] = (TFile*)bgList->At(n);
	  if(j<6) hh[n] = (TH1*)ff[n]->Get( Form("Andrey/met0_et_%i",j) )->Clone();
	  else hh[n] = (TH1*)ff[n]->Get( Form("Andrey/met0_et_%i",j) )->Clone();
	  //cout<<n<<"Yields:: file  "<<ff[n] -> GetName()<<"   histoname: "<<hh[n]->GetName()<<endl;
	  
	  hh[n] -> Scale(lumi/1000);
	  bins = hh[n] -> GetNbinsX();
	  //Float_t val = hh[n]->Integral(0,bins+1);
	  Double_t err;
	  Double_t val = hh[n]->TH1::IntegralAndError(0,bins+1, err);
	  if(j<=6) oo<<val<<separator;
	  else oo<<val<<pmSign<<err<<separator;
	  total_bg += val;
	  total_bgError += err;
	}      


      if(j<=7){
	TH1* data = (TH1*)dataFile->Get(  Form("Andrey/met0_et_%i",j) )->Clone();
	Int_t nBins = data->GetNbinsX();
	Float_t iBkg = total_bg;
	oo.precision(0);

	Float_t iData = data->Integral(0,nBins+1);
	oo<<iBkg<<separator<<iData<<endLine<<endl;

      }
      else{
	TH1* data = (TH1*)dataFile->Get(  Form("Andrey/met0_et_%i",j) )->Clone();

	TFile* ggHFile  = (TFile*)ggHList->At(j-8); //higgs 
	TFile* vbfFile  = (TFile*)vbfList->At(j-8); //higgs 

	//if(j==15)  
	//  TFile* ggHFile  = (TFile*)ggHList->At(1);
       
	TH1* ggH  = (TH1*)ggHFile ->Get(  Form("Andrey/met0_et_%i",j))->Clone();
	TH1* vbf  = (TH1*)vbfFile ->Get(  Form("Andrey/met0_et_%i",j))->Clone();
	ggH -> Scale(lumi/1000);
	vbf -> Scale(lumi/1000);
	Int_t sig_nBins = ggH->GetNbinsX();
	Int_t bkg_nBins = hh[0]->GetNbinsX();
	Int_t nBins=0;
	if (sig_nBins!=bkg_nBins) cout<<" IN yields \n WARNING:  different number of bis"<<endl;
	else nBins = sig_nBins;
	
	Float_t iggH = ggH->Integral(0,nBins+1); //Count overflows
	Float_t ivbf = vbf->Integral(0,nBins+1); //Count overflows
	Float_t iSig = iggH + ivbf;
	Float_t iBkg = total_bg;
	Float_t iBkg_err = total_bgError;
	Float_t SB   = 0.1*iSig/iBkg;        //Higgs signal is 10*real in the histograms
	//Float_t SrootB  = 0.1*iSig/sqrt(iBkg);
	Float_t SrootSB = 0.1*iSig/sqrt(iBkg + 0.1*iSig);

	Float_t iData = data->Integral(0,nBins+1);

	oo.precision(2);
	oo<<iBkg<<pmSign<<iBkg_err<<separator;
	oo.precision(0);
	oo<<iData<<separator;
	oo.precision(2);
	oo<<0.1*iggH<<separator<<0.1*ivbf<<separator;
	oo.precision(3);
	oo<<SrootSB<<separator<<SB<<endLine<<endl;  

      }
      
    }

  oo<<"end of printing yields"<<endl;
  oo.close();
  
}

THStack* makeStack(TList *sourcelist, TString name, Float_t lumi)
{
  //cout<<"making stack!"<<endl;
  THStack * hs = new THStack(Form("hs_11_%i",1),"Stacked hist");
  //sourcelist ->Print();

  TFile *first_source = (TFile*)sourcelist->First();
  first_source->cd("Andrey");
  TDirectory *current_sourcedir = gDirectory;

  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key, *oldkey=0;
  while ( (key = (TKey*)nextkey())) 
    {
      if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;
      //first_source->cd();
 
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1")) continue;
      
      if(strncmp(Form("%s",name.Data()),key->GetName(), name.Length()) ==0)
	{
	  //cout<<"first "<<key->GetName()<<endl;
	  //hs -> Add((TH1*)key->ReadObj());
	  TH1 *hh1 = (TH1*)key->ReadObj()->Clone();
	  hh1 -> Scale(lumi/1000);
	  hs  -> Add(hh1);
	  //hh1 -> Print();
	  //hs -> Add((TH1*)key->ReadObj()->Clone());
	  //hs -> Add( ((TH1*)key->ReadObj()->Clone())->Scale(lumi/1000));
	}
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      while ( nextsource )
	{
	  nextsource->cd("Andrey");
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

void drawMuliPlot(TString xtitle, Int_t isLog, Float_t y1min, Float_t y1max, Float_t y2min, Float_t y2max, THStack *hs, TCanvas *cc, TLegend *leg, TList *list, Float_t lumi)
{
  //Find the samm objects as in hs stack
  TString name = Form("Andrey/%s", hs->GetHists()->First()->GetName());
  //data -> ls();
  cout<<"drawMultiPlot:: " <<name<<endl;

  TFile *ff[5];   Int_t size = list->GetSize();
  TH1 *hh[5]; 
  if(size>5) {cout<<"To many plots to overlay"<<endl; return 0;}
  for(Int_t n=0; n<size; n++)
    {
      ff[n] = (TFile*)list->At(n);
      hh[n] = (TH1*)ff[n]->Get( name.Data() )->Clone();
      //cout<<n<<" file  "<<ff[n] -> GetName()<<"   histoname: "<<hh[n]->GetName()<<endl;
    }


  for(Int_t n=1; n<size; n++)
    hh[n] ->Scale(lumi/1000./10);

  TH1 *h_data = hh[0]->Clone();

  //h_data -> Print();

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
  h_data  -> Draw("same e1pl"); //DATA
  for(Int_t n=1; n<size; n++)
    hh[n] -> Draw("same hist");
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
  h1 -> SetMaximum(y2max);
  h1 -> SetMinimum(y2min);
  h1->GetYaxis()->SetNdivisions(206);
  h1->GetYaxis()->SetTitleOffset(0.4);
  h1->SetTitleSize(0.1,"XYZ");
  h1 ->SetLabelSize(0.1,"XY");
  h1 -> Draw("ep");

  pad1->cd();
  TLatex *prelim;
  prelim = new TLatex(0.30,0.95, Form("CMS Preliminary       #it{L_{int}} = %5.f pb^{-1}",lumi));
  prelim -> SetNDC(); prelim->SetTextSize(0.03); 
  prelim -> Draw();

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
