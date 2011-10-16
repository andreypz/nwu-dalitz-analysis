{
  Float_t a, b, c, d, e;

  for(Int_t i = 0; i<8; i++)
    {
      a = 250 +50*i;
      b  = 4.121E-06*a*a*a - 5.364E-03*a*a + 2.860*a - 222.6;
      c  = 0.0008*a*a+0.7341*a +38.286;
      d = -2.606E-06*a*a*a + 2.828E-03*a*a - 6.259E-01*a + 79.82;
      e = 4.858E-06*a*a - 5.459E-03*a + 1.529;

      TString sep = ", ";
      
      cout.precision(0); cout.setf(ios::fixed, ios::floatfield);
      cout<<sep<<a<<sep<<b<<sep<<c<<sep<<d<<sep;
      cout.precision(2);
      cout<<e<<sep<<endl;
    }
}
