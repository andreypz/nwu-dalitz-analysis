#include <iostream>
#include <string>
using namespace std;


Float_t getXsecOrColors(string dataset="ggHZZ300", Int_t b = 1)
{
  Int_t lCol, fCol;

  //cout<<"Getting cross section and colors for  "<<dataset<<endl;

  Float_t cs =-1, nevt = -1;

  if(dataset=="DATA")          {lCol = -1;        fCol = -1;       cs = -1;       nevt = -1; }
  else if(dataset=="ZZ")       {lCol = kBlue;     fCol = kMagenta; cs = 0.1787;   nevt = 1103468;}
  else if(dataset=="Wjets")    {lCol = kGreen-2;  fCol = kCyan+1;  cs = 31314.0;  nevt = 81.251e6;}
  else if(dataset=="ttbar")    {lCol = kOrange+1; fCol = kGreen+2; cs = 16.71;    nevt = 10.34e6;}
  else if(dataset=="ggWW")     {lCol = kBlue-2;   fCol = kSpring;  cs = 0.15048;  nevt = 109.99e3;}
  else if(dataset=="WW")       {lCol = kOrange+1; fCol = kSpring;  cs = 4.783;    nevt = 1197558;}
  else if(dataset=="WZ")       {lCol = kRed+2;    fCol = kBlue-1;  cs = 0.856;    nevt = 1221134;}
  else if(dataset=="tW")       {lCol = kRed;      fCol = kWhite;   cs = 7.87;     nevt = 814390 ;}
  else if(dataset=="tbarW")    {lCol = kRed;      fCol = kWhite;   cs = 7.87;     nevt = 809984;}
  else if(dataset=="DYmumu")   {lCol = kRed;      fCol = kWhite;   cs = 1666.0;   nevt = 29.7e6;}
  else if(dataset=="DYee")     {lCol = kRed;      fCol = kWhite;   cs = 1666.0;   nevt = 29.5e6;}
  else if(dataset=="DYtautau") {lCol = kRed;      fCol = kWhite;   cs = 1666.0;   nevt = 29.5e6;}
  else if(dataset=="DYjets")   {lCol = kRed;      fCol = kWhite;   cs = 3048.0;   nevt = 36.278e6;}

  else if(dataset=="ggHZZ200") {lCol = kPink;     fCol  = -1; cs = 0.05406;     nevt = 96.20e3;}
  else if(dataset=="ggHZZ250") {lCol = kCyan;     fCol  = -1; cs = 0.0398;      nevt = 99.993e3;}
  else if(dataset=="ggHZZ300") {lCol = kCyan-1;   fCol  = -1; cs = 0.0301;      nevt = 99.990e3;}
  else if(dataset=="ggHZZ350") {lCol = kCyan+2;   fCol  = -1; cs = 0.0286;      nevt = 99.998e3;}
  else if(dataset=="ggHZZ400") {lCol = kYellow;   fCol  = -1; cs = 0.0221;      nevt = 93.854e3;}
  else if(dataset=="ggHZZ450") {lCol = kYellow;   fCol  = -1; cs = 0.0143;      nevt = 96.642e3;}
  else if(dataset=="ggHZZ500") {lCol = kYellow;   fCol  = -1; cs = 0.00896;     nevt = 99.980e3;}
  else if(dataset=="ggHZZ550") {lCol = kYellow;   fCol  = -1; cs = 0.00568;     nevt = 94.585e3;}
  else if(dataset=="ggHZZ600") {lCol = kCyan+1;   fCol  = -1; cs = 0.00359;     nevt = 99.972e3;}

  else if(dataset=="ggHWW200") {lCol = kPink;     fCol  = -1; cs = 0.40840;    nevt = 109993;}
  else if(dataset=="ggHWW250") {lCol = kCyan;     fCol  = -1; cs = 0.24378;    nevt = 106038;}
  else if(dataset=="ggHWW300") {lCol = kCyan-1;   fCol  = -1; cs = 0.17598;    nevt = 109989;}
  else if(dataset=="ggHWW350") {lCol = kCyan+2;   fCol  = -1; cs = 0.16368;    nevt = 109991;}
  else if(dataset=="ggHWW400") {lCol = kYellow;   fCol  = -1; cs = 0.12418;    nevt = 109993;}
  else if(dataset=="ggHWW450") {lCol = kYellow;   fCol  = -1; cs = 0.07880;    nevt = 109977;}
  else if(dataset=="ggHWW500") {lCol = kYellow;   fCol  = -1; cs = 0.04868;    nevt = 109973;}
  else if(dataset=="ggHWW550") {lCol = kYellow;   fCol  = -1; cs = 0.03053;    nevt = 109973;}
  else if(dataset=="ggHWW600") {lCol = kCyan+1;   fCol  = -1; cs = 0.01914;    nevt = 109972;}

  else if(dataset=="VBFHZZ200"){lCol = kCyan+1;   fCol  = -1; cs = 6.562e-3;     nevt = 49.943e3;}
  else if(dataset=="VBFHZZ250"){lCol = kCyan+1;   fCol  = -1; cs = 5.163e-3;     nevt = 49.957e3;}
  else if(dataset=="VBFHZZ300"){lCol = kCyan+1;   fCol  = -1; cs = 3.734e-3;     nevt = 49.935e3;}
  else if(dataset=="VBFHZZ350"){lCol = kCyan+1;   fCol  = -1; cs = 2.644e-3;     nevt = 46.002e3;}
  else if(dataset=="VBFHZZ400"){lCol = kCyan+1;   fCol  = -1; cs = 1.760e-3;     nevt = 49.851e3;}
  else if(dataset=="VBFHZZ450"){lCol = kCyan+1;   fCol  = -1; cs = 1.299e-3;     nevt = 46.065e3;}
  else if(dataset=="VBFHZZ500"){lCol = kCyan+1;   fCol  = -1; cs = 1.001e-3;     nevt = 46.049e3;}
  else if(dataset=="VBFHZZ550"){lCol = kCyan+1;   fCol  = -1; cs = 0.876e-3;     nevt = 49.729e3;}
  else if(dataset=="VBFHZZ600"){lCol = kCyan+1;   fCol  = -1; cs = 0.634e-3;     nevt = 49.698e3;}


  Float_t ddd;
  if(b==1) ddd = (Float_t)lCol;
  if(b==2) ddd = (Float_t)fCol;
  if(b==3) ddd = cs; 
  if(b==4) ddd = nevt; 

  return ddd;
}
