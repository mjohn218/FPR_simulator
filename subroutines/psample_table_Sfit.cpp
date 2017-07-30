#include "Vrnumer.h"
#include "GF_calls.h"
#include "numeric_GF.h"
#include "rand_gsl.h"

void psample_table_Sfit(double rmin, double r0, double tcurr, double Dtot, double bindrad, double alpha, double kact, double *table, int MAXP, double delp, double Sfit, double passoc, double Rlong)
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
  double r1, term1, term2, term3, term4;
  double cumsum;
  double c3=alpha/(r0*sqrt(Dtot));
  double off=r0-2.0*bindrad;
  //double passoc=;//survive_irr(r0, tcurr, Dtot, bindrad, alpha, cof);
  double probstart=0;
  int j, jprev;
  jprev=0;
  b=alpha*sqrt_t;
  //or just skip all the probs less than pstart, so you don't need to know what p_survive is 
  int i=0;
  double prob=0;
  table[0]=bindrad;
  
  int i1=i;
  probstart=0;
  //cout <<"first index: "<<i1<<" for reflecting should be zero "<<kact<<endl; 
  double psub;
  psub=peval_cumulativeF(bindrad,  r0,  tcurr,  Dtot,  bindrad,  alpha,  kact);
  psub*=Sfit;
  double psurvive=1-passoc;
  double plong=peval_cumulativeF(Rlong, r0,  tcurr,  Dtot,  bindrad,  alpha, kact);
  plong*=Sfit;
  double ps_numer=(plong-psub);
  delp=psurvive/(1.0*MAXP);
  cout <<"calculated psurvive from fit: "<<psurvive<<" integral over numer pirr, to Rlong: "<<ps_numer<<" Rlong: "<<Rlong<<" ratio: "<<psurvive/ps_numer<<endl;
  double eps=1E-8;
  double ps_ratio=(psurvive+eps)/ps_numer;

  // cout <<"value at sigma to be subtracted: "<<psub<<endl;
  for(i=i1;i<MAXP-1;i++){
    //now find which time this corresponds to
    prob=delp*(i+1);//at final value of i, prob will be equal=1
    j=jprev;
    //cout <<"time start: "<<tmin*1<<endl;
    cumsum=probstart;
    while(prob>cumsum){
      r1=rmin*j+bindrad;

      sep=r1+off;//r0-2.0*bindrad;
      sep2=sep*sep;
      dist=r1-r0;
      d2=dist*dist;
      a=sep/sqfDt;
      term1=-0.5*erf(-dist/sqfDt)-0.5*fDt*c1*exp(-d2/fDt);
      term2=-0.5*(r0-2.0*bindrad)/r0*erf(sep/sqfDt)-0.5*fDt*c1*exp(-sep2/fDt);
      term3=c3*sqrt(4.0*Dtot)/(2.0*alpha*spi)*( -sqfDt*exp(-a*a)-spi*off*erf(a) );
      term4=c3*-sqfDt/(4.0*b*b)*( sqfDt*erf(a)+(sqfDt-2.0*b*r1)*exp(2.0*a*b+b*b)*erfc(a+b) );
      cumsum=Sfit*(term1+term2-term3-term4)*ps_ratio+probstart-psub*ps_ratio;
      
      /*this is the cumulative probability for a given r value, compare this sum to a random number to 
	generate the separation
	we want the probability to be evenly spaced
	It will sum up to the survival probability, since it survived, the random number will be between p_assoc=1-survival and 1. 
      */
      j++;
    }
    if(j-jprev>10000){
      rmin*=10;
      j/=10;
      j-=10;
      
    }
    //    cout <<"iterations: "<<j-jprev<<" time: "<<time<<" tmin: "<<tmin<<endl;
    jprev=int(round(j-1));
    //now add this to the table, so it can be looked up.
    //except it is evenly space in time, not in prob, so interpolate?
    table[i+1]=r1;
    //cout <<r1<<'\t'<<prob-probstart<<endl;
  }
  cout <<r0<<'\t'<<r1<<'\t'<<prob-probstart<<" Final prob1? "<<prob<<endl;
  //last value
  prob=delp*(i+1)-1E-8;//at final value of i, prob will be equal=1
  j=jprev;
  double small=1E-16;
  double pirr1=1;
  cumsum=probstart;
  while(prob>cumsum && pirr1>small){
    r1=rmin*(j+1)+bindrad;
    
    sep=r1+off;//r0-2.0*bindrad;
    sep2=sep*sep;
    dist=r1-r0;
    d2=dist*dist;
    a=sep/sqfDt;
    term1=-0.5*erf(-dist/sqfDt)-0.5*fDt*c1*exp(-d2/fDt);
    term2=-0.5*(r0-2.0*bindrad)/r0*erf(sep/sqfDt)-0.5*fDt*c1*exp(-sep2/fDt);
    term3=c3*sqrt(4.0*Dtot)/(2.0*alpha*spi)*( -sqfDt*exp(-a*a)-spi*off*erf(a) );
    term4=c3*-sqfDt/(4.0*b*b)*( sqfDt*erf(a)+(sqfDt-2.0*b*r1)*exp(2.0*a*b+b*b)*erfc(a+b) );
    cumsum=Sfit*(term1+term2-term3-term4)*ps_ratio+probstart-psub*ps_ratio;
    pirr1= pirrev_valueF(r1, r0,tcurr,  Dtot,  bindrad,  alpha);
    pirr1*=Sfit*ps_ratio;
    /*this is the cumulative probability for a given r value, compare this sum to a random number to 
      generate the separation
      we want the probability to be evenly spaced
      It will sum up to the survival probability, since it survived, the random number will be between p_assoc=1-survival and 1. 
    */
    j++;
  }

  table[i+1]=r1;
  cout<<"Final r to get to prob1: "<<r1<<'\t'<<prob-probstart<<endl;

}
