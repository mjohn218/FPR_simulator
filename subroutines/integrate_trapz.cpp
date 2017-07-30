#include "2Drelated.h"

double integrate_trapz(int r0ind, double bindrad, double rstar, gsl_matrix* pirr, double RstepSize, int N)
{
  int i;
  int indexr0=r0ind;
  double cumtrap=0;
  //int N=int( (rstar-bindrad)/RstepSize);
  double r1, pirrvalue;
  for(i=1;i<N-1;i++){
    r1=i*RstepSize+bindrad;
    pirrvalue=gsl_matrix_get(pirr, indexr0, i);
    cumtrap+=r1*pirrvalue*1.0*RstepSize;
  }
  i=0;
  r1=i*RstepSize+bindrad;
  pirrvalue=gsl_matrix_get(pirr, indexr0, i);
  cumtrap+=pirrvalue*r1*RstepSize*0.5;
  i=N-1;
  r1=i*RstepSize+bindrad;
  pirrvalue=gsl_matrix_get(pirr, indexr0, i);
  cumtrap+=pirrvalue*r1*RstepSize*0.5;
  cumtrap*=2.0*M_PI;
  //  cout <<"Integratl trap: "<<cumtrap<<endl;
  return cumtrap;

}
