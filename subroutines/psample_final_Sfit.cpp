#include "Vrnumer.h"
#include "GF_calls.h"
#include "numeric_GF.h"
#include "rand_gsl.h"

double psample_final_Sfit(double rmin, double r0, double tcurr, double Dtot, double bindrad, double kact, double Sfit, double passoc, double Rlong, double R_table_end, double rnum)
{
  /*To sample from irreversible distribution, kact is finite, and alpha=(1+kact/kdiff)*sqrt(D)/sigma
    To sample from reflecting distribution, kact=0 and alpha=sqrt(D)/sigma
  */
  
  double r1, term1, term2, term3, term4;
  double cumsum;
  
  
  double probstart=0;
  int j, jprev;
  jprev=0;

  //or just skip all the probs less than pstart, so you don't need to know what p_survive is 
  int i=0;
  double prob=0;

  int i1=i;
  //cout <<"first index: "<<i1<<" for reflecting should be zero "<<kact<<endl; 
  
  double kdiff=4.0*M_PI*Dtot*bindrad;
  double fact=1.0+kact/kdiff;
  double alpha=fact*sqrt(Dtot)/bindrad;
  
  double psub;
  psub=peval_cumulativeF(bindrad,  r0,  tcurr,  Dtot,  bindrad,  alpha,  kact);
  psub*=Sfit;//FitMat[2];

  double plong=peval_cumulativeF(Rlong,  r0,  tcurr,  Dtot,  bindrad,  alpha,  kact);
  plong*=Sfit;
  
  double ps_numer=(plong-psub);//first half

  double psurvive=1-passoc;
  //cout <<"calculated psurvive from fit: "<<psurvive<<" integral over numer pirr, to Rlong: "<<ps_numer<<" Rlong: "<<Rlong<<" ratio: "<<psurvive/ps_numer<<endl;
  double eps=1E-8;
  double ps_ratio=(psurvive+eps)/ps_numer;
  
  //cout <<"value at sigma to be subtracted: "<<psub<<endl;
  /*Since pirr may not integrate to psurvive, rescale it to. */
  probstart=0;
  prob=rnum;//rnum will be between delp*MAXP-passoc and psurvive.

  j=0;//jprev;
  //cout <<"time start: "<<tmin*1<<endl;
  cumsum=0;//probstart;
  double pcume;
  while(prob>cumsum){
    r1=rmin*j+R_table_end;
    /*need to use other fit parms, and subtract off sum at previous r1.*/
    
    pcume=peval_cumulativeF(r1, r0,  tcurr,  Dtot,  bindrad,  alpha,  kact);
    cumsum=Sfit*pcume*ps_ratio-psub*ps_ratio;
      
    /*this is the cumulative probability for a given r value, compare this sum to a random number to 
      generate the separation
      we want the probability to be evenly spaced
      It will sum up to the survival probability, since it survived, the random number will be between p_assoc=1-survival and 1. 
    */
    j++;
  }
  //  cout <<r0<<'\t'<<r1<<'\t'<<" Final prob1? "<<prob<<" Psurvive: "<<psurvive<<" Randonum Prob to match: "<<rnum<<endl;
  return r1;
}
