In order to use the weighting in root, one needs to create a shared library:
   
gfortran -shared -O2 *.f -o libShapeLine.so -fPIC

Then, load this library to root by gSystem->Load("libShapeLine.so");
 
Define an external function in .h file: 
 extern"C" {
   void pwhg_cphto_reweight_(double *mh, double *gh, double *mt, int *BWflag, double *m, double *w);
           }
