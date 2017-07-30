#include "2Drelated.h"

double integrate_numer(int r0ind, double bindrad, double rstar, gsl_matrix* pirr, double RstepSize, int Npt)
{
  
  int Ninterval=Npt-1;
  double sum;
  //cout <<"Npt: "<<Npt<<" Ninterval: "<<Ninterval<< " r0ind: "<<r0ind<<" r0val: "<<r0ind*RstepSize+bindrad<<" rval: "<<(Npt-1)*RstepSize+bindrad<<endl;
  // if(Ninterval%3==0)
//     sum=integrate_simp38(r0ind, bindrad, rstar, pirr, RstepSize, Npt);
//   if(Ninterval%2==0)
//     sum=integrate_simp(r0ind, bindrad, rstar, pirr, RstepSize, Npt);
  
  sum=integrate_trapz( r0ind,  bindrad,  rstar, pirr, RstepSize,  Npt);
  return sum;
}
