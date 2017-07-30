#include "Vrnumer.h"
#include "GF_calls.h"
#include "numeric_GF.h"
#include "Faddeeva.hh"

using namespace std;

double peval_cumulativeF(double r1, double r0, double tcurr, double Dtot, double bindrad, double alpha, double kact)
{
  /*To sample from irreversible distribution, kact is finite, and alpha=(1+kact/kdiff)*sqrt(D)/sigma
    To sample from reflecting distribution, kact=0 and alpha=sqrt(D)/sigma
  */
  
  double fDt=4.0*Dtot*tcurr;
  double sqfDt=sqrt(fDt);
  double kdiff=4.0*M_PI*Dtot*bindrad;
  double cof=kact/(kact+kdiff);
  double spi=sqrt(M_PI);
  double c1=1.0/(2.0*r0*sqrt(M_PI*Dtot*tcurr));
  double sep, dist, sep2, d2;
  double sqrt_t=sqrt(tcurr);
  double a, b;
  double  term1, term2, term3, term4;
  double cumsum;
  double c3=alpha/(r0*sqrt(Dtot));
  double off=r0-2.0*bindrad;
  int j, jprev;
  jprev=0;
  b=alpha*sqrt_t;
  //or just skip all the probs less than pstart, so you don't need to know what p_survive is 
  
  
  sep=r1+off;//r0-2.0*bindrad;
  sep2=sep*sep;
  dist=r1-r0;
  d2=dist*dist;
  a=sep/sqfDt;
  term1=-0.5*erf(-dist/sqfDt)-0.5*fDt*c1*exp(-d2/fDt);
  term2=-0.5*(r0-2.0*bindrad)/r0*erf(sep/sqfDt)-0.5*fDt*c1*exp(-sep2/fDt);
  term3=c3*sqrt(4.0*Dtot)/(2.0*alpha*spi)*( -sqfDt*exp(-a*a)-spi*off*erf(a) );
  //term4=c3*-sqfDt/(4.0*b*b)*( sqfDt*erf(a)+(sqfDt-2.0*b*r1)*exp(2.0*a*b+b*b)*erfc(a+b) );
  double ef1=a+b;
  double ep1=2.0*a*b+b*b;
  double termF;
  if( isinf(exp(ep1) )){
  
    std::complex<double> z;
    z.real(0.0);//=0.0;
    z.imag(ef1);//=ef1;
    //cout <<"Complex number: "<<z<<endl;
    std::complex<double> value;
    double relerr=0;
    value=Faddeeva::w(z, relerr);
    double ea2=exp(-a*a);
    termF=ea2*real(value);
    
  }else{
    termF=exp(ep1)*erfc(ef1);
  }
  term4=-c3*sqfDt/(4.0*b*b)*( sqfDt*erf(a)+(sqfDt-2.0*b*r1)*termF );
  cumsum=term1+term2-term3-term4;
      
  /*this is the cumulative probability for a given r value, compare this sum to a random number to 
    generate the separation
	we want the probability to be evenly spaced
	It will sum up to the survival probability, since it survived, the random number will be between p_assoc=1-survival and 1. 
      */
  return cumsum;
  
}
