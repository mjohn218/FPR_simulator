#include "Vrnumer.h"
#include "GF_calls.h"
#include "numeric_GF.h"
#include "rand_gsl.h"

double psample_final_pfree_vr(double rmin, double r0, double tcurr, double Dtot, double bindrad, double rc, double passoc,double R_table_end, double rnum)
{
  double r1, term1;
  double cumsum;


  double psurvive=1.0-passoc;
  int j, jprev;
  jprev=0;

  //or just skip all the probs less than pstart, so you don't need to know what p_survive is 
  int i=0;
  double prob=0;
  int i1=i;
  //cout <<"first index: "<<i1<<" for reflecting should be zero "<<kact<<endl; 
  
  double psigma, pnorm;
  psigma=vr_pfree_cumulative(bindrad,  r0, tcurr, Dtot,  bindrad, rc);
  //  pnorm=1-psigma;//distibution at inf is 1. this is integral from sigma to inf.

  // cout <<"value at sigma to be subtracted: "<<psub<<endl;
  //now find which time this corresponds to
  prob=rnum;
  j=0;
  //cout <<"time start: "<<tmin*1<<endl;
  cumsum=0;
  double jmax=1E7;
  while(prob>cumsum && j<jmax){
    r1=rmin*j+R_table_end;

    term1=vr_pfree_cumulative(r1, r0, tcurr, Dtot, bindrad, rc);
    cumsum=(term1-psigma)*psurvive;
    
    /*this is the cumulative probability for a given r value, compare this sum to a random number to 
      generate the separation
      we want the probability to be evenly spaced
      It will sum up to the survival probability, since it survived, the random number will be between p_assoc=1-survival and 1. 
    */
    j++;
  }
  if(j>=jmax)
    cout <<"broke out of final sample loop, pfree_vr. final r:  "<<r1<<" prob: "<<prob<<" final cumulative sum: "<<cumsum<<endl;

  //cout <<r0<<'\t'<<r1<<'\t'<<" Final prob1? "<<prob<<" psurvive: "<<psurvive<<endl;
  return r1;
}
