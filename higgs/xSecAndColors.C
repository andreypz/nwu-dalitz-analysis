#include <iostream>
#include <string>
using namespace std;


Float_t getCS(string dataset="MC_ggH400")
{
  Float_t cs =-1;

  cout<<"Getting cross section for  "<<dataset<<endl;

  if(dataset=="2011A_May10_DoubleMu") cs         = 0.;
 else if(dataset=="2011A_HZZSkim_DoubleMu") cs       = 0.;
 else if(dataset=="2011A_May10_DoubleElectron") cs   = 0.;
 else if(dataset=="2011A_HZZSkim_DoubleElectron") cs = 0.;

 else if(dataset=="MC_ggH200") cs   = 0.06277;
 else if(dataset=="MC_ggH400") cs   = 0.02229;
 else if(dataset=="MC_ggH600") cs   = 0.004232;

 else if(dataset=="MC_ggH200_WW2l2nu") cs   = 0.4721*4./9;
 else if(dataset=="MC_ggH400_WW2l2nu") cs   = 0.1253*4./9;
 else if(dataset=="MC_ggH600_WW2l2nu") cs   = 0.02256*4./9;

 else if(dataset=="MC_ggH200_WW2t2nu") cs   = 0.4721*1./9;
 else if(dataset=="MC_ggH400_WW2t2nu") cs   = 0.1253*1./9;
 else if(dataset=="MC_ggH600_WW2t2nu") cs   = 0.02256*1./9;

 else if(dataset=="MC_ZZtoAny") cs   = 4.37; //4.37 from Alexei's email. 5.3 LO? 5.9 is NLO
 else if(dataset=="MC_Wjets") cs     = 31314.0;
 else if(dataset=="MC_ZllG") cs      = 14.3;
 else if(dataset=="MC_tt2l2nu2b") cs = 16.5;
 else if(dataset=="MC_tSchannel") cs = 1.36;
 else if(dataset=="MC_tTchannel") cs = 20.9;
 else if(dataset=="MC_tWchannel") cs = 10.6;
 else if(dataset=="MC_DYmumu") cs    = 1666.;
 else if(dataset=="MC_DYee") cs      = 1666.;
 else if(dataset=="MC_DYtautau") cs  = 1666.;
 else if(dataset=="MC_WW") cs        =  4.51;
 else if(dataset=="MC_WZ") cs        =  0.596;

 else if(dataset=="MC_Zbb0") cs =  2.92;
 else if(dataset=="MC_Zbb1") cs =  1.65;
 else if(dataset=="MC_Zbb2") cs =  0.623;
 else if(dataset=="MC_Zbb3") cs =  0.274;
 else if(dataset=="MC_Zcc0") cs =  2.92;
 else if(dataset=="MC_Zcc1") cs =  1.63;
 else if(dataset=="MC_Zcc2") cs =  0.627;
 else if(dataset=="MC_Zcc3") cs =  0.281;
  
  return cs;
}


int getColors(string dataset="MC_ggH400", Int_t b =1)
{
  Int_t lCol, fCol;

  // cout<<"Getting Colors for  "<<dataset<<endl;

  if(dataset=="2011A_May10_DoubleMu")              {lCol = -1;  fCol = -1;}
  else if(dataset=="2011A_HZZSkim_DoubleMu")       {lCol = -1;  fCol = -1;}
  else if(dataset=="2011A_May10_DoubleElectron")   {lCol = -1;  fCol = -1;}
  else if(dataset=="2011A_HZZSkim_DoubleElectron") {lCol = -1;  fCol = -1;}

  else if(dataset=="MC_ggH200") {lCol = kWhite;  fCol  = -1;}
  else if(dataset=="MC_ggH400") {lCol = kCyan;    fCol  = -1;}
  else if(dataset=="MC_ggH600") {lCol = kCyan+1;  fCol  = -1;}

  else if(dataset=="MC_ggH200_WW2l2nu") {lCol = kWhite;  fCol  = -1;}
  else if(dataset=="MC_ggH400_WW2l2nu") {lCol = kCyan;    fCol  = -1;}
  else if(dataset=="MC_ggH600_WW2l2nu") {lCol = kCyan+1;  fCol  = -1;}
  else if(dataset=="MC_ggH200_WW2t2nu") {lCol = kWhite;  fCol  = -1;}
  else if(dataset=="MC_ggH400_WW2t2nu") {lCol = kCyan;    fCol  = -1;}
  else if(dataset=="MC_ggH600_WW2t2nu") {lCol = kCyan+1;  fCol  = -1;}

  else if(dataset=="MC_ZZtoAny")   {lCol = kBlue;     fCol = kMagenta;}
  else if(dataset=="MC_Wjets")     {lCol = kGreen-2;  fCol = kCyan+1;}
  else if(dataset=="MC_ZllG")      {lCol = kYellow;   fCol = -1;}
  else if(dataset=="MC_tt2l2nu2b") {lCol = kOrange+1; fCol = kGreen+2;}
  else if(dataset=="MC_WW")        {lCol = kBlue-2;  fCol = kSpring;}
  else if(dataset=="MC_WZ")        {lCol = kRed+2;   fCol = kBlue-1;}
  else if(dataset=="MC_tSchannel") {lCol = -1;  fCol = -1;}
  else if(dataset=="MC_tTchannel") {lCol = -1;  fCol = -1;}
  else if(dataset=="MC_tWchannel") {lCol = -1;  fCol = -1;}
  else if(dataset=="MC_DYmumu")    {lCol = -1;  fCol = -1;}
  else if(dataset=="MC_DYee")      {lCol = -1;  fCol = -1;}
  else if(dataset=="MC_DYtautau")  {lCol = -1;  fCol = -1;}
  else if(dataset=="MC_Zbb0") {lCol = -1;  fCol = -1;}
  else if(dataset=="MC_Zbb1") {lCol = -1;  fCol = -1;}
  else if(dataset=="MC_Zbb2") {lCol = -1;  fCol = -1;}
  else if(dataset=="MC_Zbb3") {lCol = -1;  fCol = -1;}
  else if(dataset=="MC_Zcc0") {lCol = -1;  fCol = -1;}
  else if(dataset=="MC_Zcc1") {lCol = -1;  fCol = -1;}
  else if(dataset=="MC_Zcc2") {lCol = -1;  fCol = -1;}
  else if(dataset=="MC_Zcc3") {lCol = -1;  fCol = -1;}

  //fCol = new TString(fColor.c_str());
  //lCol = new TString(lColor.c_str());
  int ddd;
  if(b==1) ddd = lCol;
  if(b==2) ddd = fCol;

  return ddd;
}
