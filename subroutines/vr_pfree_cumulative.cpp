#include "Vrnumer.h"
#include "GF_calls.h"
#include "numeric_GF.h"
#include "rand_gsl.h"

double vr_pfree_cumulative(double rcurr, double r0, double tcurr, double Dtot, double bindrad, double rc)
{
  /*Psurvive for the first time step is 1-passoc and should be exact. After
    more than one step the distribution will be off so passoc will be off as well*/
  
  double fDt=4.0*Dtot*tcurr;
  double sq_fDt=sqrt(fDt);
  double fterm=Dtot*tcurr*rc/(r0*r0);
  double f0=r0+fterm;
  double f1=1.0/(sqrt(4.0*M_PI*tcurr));
  double f2=1.0/(4.0*M_PI*f0*sqrt(Dtot));
  double cof=f1*f2;//equal to 1/( 8pir_0 sqrt(piDt))
  int i, j;
  double sep, dist;
  double sqrt_t=sqrt(tcurr);
  

  double r1, term1, term2, e1, ef1, sum;
  double adist;
  /*Normalization for no diffusion inside binding radius!
    with the addition of the potential, instead of integrating of exp(-(r-r0)^2/s) and exp(-(r+r0)^2/s) replace f0=r0+F(r0).
  */
  double c1=4.0*M_PI*cof;
  
  dist=bindrad-f0;
  adist=bindrad+f0;
  term1=-0.5*fDt*exp(-dist*dist/fDt)-0.5*sqrt(4.0*M_PI*Dtot*tcurr)*f0*erf(-dist/sqrt(4.0*Dtot*tcurr));
  term2=0.5*fDt*exp(-adist*adist/fDt)+0.5*sqrt(4.0*M_PI*Dtot*tcurr)*f0*erf(adist/sqrt(4.0*Dtot*tcurr));
  double pnorm=1.0-c1*(term1+term2);//this is then the normlization from sigma->inf 
  
  dist=rcurr-f0;
  adist=rcurr+f0;
  term1=-0.5*fDt*exp(-dist*dist/fDt)-0.5*sqrt(4.0*M_PI*Dtot*tcurr)*f0*erf(-dist/sqrt(4.0*Dtot*tcurr));
  term2=0.5*fDt*exp(-adist*adist/fDt)+0.5*sqrt(4.0*M_PI*Dtot*tcurr)*f0*erf(adist/sqrt(4.0*Dtot*tcurr));
  double pcume=c1*(term1+term2);//this is the value of integrated function at rcurr.
  return pcume/pnorm;

  
    

}
