#include "Vrnumer.h"
#include "GF_calls.h"
#include "numeric_GF.h"
#include "rand_gsl.h"

void psample_table_pfree_vr(double rmin, double r0, double tcurr, double Dtot, double bindrad, double rc,  double *table, int MAXP, double delp, double passoc)
{
  double r1, term1;
  double cumsum;

  
  double psurvive=1.0-passoc;
  int j, jprev;
  jprev=0;
  
  //or just skip all the probs less than pstart, so you don't need to know what p_survive is 
  int i=0;
  double prob=0;
  table[0]=bindrad;
  // while(prob<probstart){
//     i++;
//     prob=delp*(i+1);
//     table[i]=bindrad;//these probabilities will not be sample because they are <the prob of associating
//   }
  int i1=i;
  //cout <<"first index: "<<i1<<" for reflecting should be zero "<<kact<<endl; 
  
  double psigma;//, pnorm;
 
  double eps=1E-8;
  /*Inside vr_pfree_cumulative, value is divided already by pnorm*/
  psigma=vr_pfree_cumulative(bindrad,  r0, tcurr, Dtot,  bindrad, rc);
  //  pnorm=1-psigma;//distibution at inf is 1. this is integral from sigma to inf.
  delp=psurvive/(1.0*MAXP);
  psurvive+=eps;
  // cout <<"value at sigma to be subtracted: "<<psub<<endl;
  for(i=i1;i<MAXP-1;i++){
    //now find which time this corresponds to
    prob=delp*(i+1);//at final value of i, prob will be equal=1
    j=jprev;
    //cout <<"time start: "<<tmin*1<<endl;
    cumsum=0;
    while(prob>cumsum){
      r1=rmin*j+bindrad;

      term1=vr_pfree_cumulative(r1, r0, tcurr, Dtot, bindrad, rc);
      cumsum=(term1-psigma)*psurvive;
      
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
  cout <<r0<<'\t'<<r1<<'\t'<<" Final prob1? "<<prob<<endl;
  //last value
  // prob=delp*(i+1)-1E-8;//at final value of i, prob will be equal=1
//   j=jprev;
//   double pfree=1;
//   double small=1E-16;
//   cumsum=probstart;
//   while(prob>cumsum && pfree>small){
//     r1=rmin*(j+1)+bindrad;

//     term1=vr_pfree_cumulative(r1, r0, tcurr, Dtot, bindrad, rc);
//     cumsum=(term1-psigma)/pnorm*psurvive+probstart;
//     pfree=vr_pfree_value_norm(r1,  r0,  tcurr,  Dtot,  bindrad,  rc);
//     pfree*=psurvive;
//     /*this is the cumulative probability for a given r value, compare this sum to a random number to 
//       generate the separation
//       we want the probability to be evenly spaced
//       It will sum up to the survival probability, since it survived, the random number will be between p_assoc=1-survival and 1. 
//     */
//     j++;
//   }
//   //now add this to the table, so it can be looked up.
//   //except it is evenly space in time, not in prob, so interpolate?
//   table[i+1]=r1;
//   cout<<"Final r to get to prob1: "<<r1<<'\t'<<prob-probstart<<endl;

}
