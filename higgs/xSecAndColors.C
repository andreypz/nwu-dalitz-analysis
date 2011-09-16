#include <iostream>
#include <string>
using namespace std;


Float_t getXsecOrColors(string dataset="ggH200", Int_t b = 1)
{
  Int_t lCol, fCol;

  //cout<<"Getting cross section and colors for  "<<dataset<<endl;

  Float_t cs =-1;

  if(dataset=="DATA")        {lCol = -1;        fCol = -1;       cs = -1; }
  else if(dataset=="ZZ")     {lCol = kBlue;     fCol = kMagenta; cs = 0.214; }
  else if(dataset=="Wjets")  {lCol = kGreen-2;  fCol = kCyan+1;  cs = 31314.0; }
  else if(dataset=="ttbar")  {lCol = kOrange+1; fCol = kGreen+2; cs = 16.71; }
  else if(dataset=="WW")     {lCol = kBlue-2;   fCol = kSpring;  cs = 4.51; }
  else if(dataset=="WZ")     {lCol = kRed+2;    fCol = kBlue-1;  cs = 0.596; }
  else if(dataset=="tW")     {lCol = kRed;        fCol = kWhite;   cs = 7.87; }
  else if(dataset=="DYmumu") {lCol = kRed;        fCol = kWhite;   cs = 1666.0; }
  else if(dataset=="DYee")   {lCol = kRed;        fCol = kWhite;   cs = 1666.0; }
  else if(dataset=="DYtautau") {lCol = kRed;      fCol = kWhite;   cs = 1666.0; }
  else if(dataset=="ggHZZ140") {lCol = kBlue-1;  fCol  = -1; cs = 0.022; }
  else if(dataset=="ggHZZ200") {lCol = kPink;  fCol  = -1;   cs = 0.06277; }
  else if(dataset=="ggHZZ250") {lCol = kBlue-1;  fCol  = -1; cs = 0.0398; }
  else if(dataset=="ggHZZ300") {lCol = kCyan;    fCol  = -1; cs = 0.0301; }
  else if(dataset=="ggHZZ350") {lCol = kCyan;    fCol  = -1; cs = 0.0286; }
  else if(dataset=="ggHZZ400") {lCol = kYellow;  fCol  = -1; cs = 0.0221; }
  else if(dataset=="ggHZZ450") {lCol = kYellow;  fCol  = -1; cs = 0.0221; }
  else if(dataset=="ggHZZ600") {lCol = kCyan+1;  fCol  = -1; cs = 0.004232; }


  Float_t ddd;
  if(b==1) ddd = (Float_t)lCol;
  if(b==2) ddd = (Float_t)fCol;
  if(b==3) ddd = cs; 

  return ddd;
}
