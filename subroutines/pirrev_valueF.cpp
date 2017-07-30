#include "Vrnumer.h"
#include "GF_calls.h"
#include "numeric_GF.h"
#include "Faddeeva.hh"

double pirrev_valueF(double rcurr, double r0, double tcurr, double Dtot, double bindrad, double alpha)
{

  double fDt=4.0*Dtot*tcurr;
  double sq_fDt=sqrt(fDt);

  double f1=1.0/(sqrt(4.0*M_PI*tcurr));
  double f2=1.0/(4.0*M_PI*r0*sqrt(Dtot));

  double sep, dist;
  double sqrt_t=sqrt(tcurr);
  double a2=alpha*alpha;
  double r1, term1, term2, e1, ef1, sum;
  
  r1=rcurr;
  sep=r1+r0-2.0*bindrad;
  dist=r1-r0;
  term1=f1*( exp(-dist*dist/fDt)+exp(-sep*sep/fDt) );
  double a=sep/sq_fDt;
  double b=sqrt_t*alpha;
  e1=2.0*a*b+a2*tcurr;
  ef1=a+b;
  double ep1=exp(e1);
  if(isinf(ep1)){
    std::complex<double> z;
    z.real(0.0);//=0.0;
    z.imag(ef1);//=ef1;
    //cout <<"Complex number: "<<z<<endl;
    std::complex<double> value;
    double relerr=0;
    value=Faddeeva::w(z, relerr);
    double ea2=exp(-a*a);
    term2=ea2*real(value);
    
    
  }else{
    term2=exp(e1)*erfc(ef1);
  }
  
  sum=term1-alpha*term2;
  sum*=f2/r1;
  double pirr=sum;
  return pirr;
  
  
  
}
