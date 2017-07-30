#include "Vrnumer.h"
#include "GF_calls.h"
#include "numeric_GF.h"
#include "rand_gsl.h"

double psample_final_pfree_2Sfit(double rmin, double r0, double tcurr,  double bindrad, double passoc, Vrnumer &limit, double *FitMat, double Rlong, double R_table_end, double rnum)
{
  /*To sample from irreversible distribution, kact is finite, and alpha=(1+kact/kdiff)*sqrt(D)/sigma
    To sample from reflecting distribution, kact=0 and alpha=sqrt(D)/sigma
  */
  
  double r1, term1, term2, term3, term4;
  double cumsum;
  int j, jprev;
  jprev=0;

  //or just skip all the probs less than pstart, so you don't need to know what p_survive is 
  int i=0;
  double prob=0;
  int i1=i;
  //cout <<"first index: "<<i1<<" for reflecting should be zero "<<kact<<endl; 
  double psub;
  double psigma;
  double Dfit1=FitMat[0];
  double rcfit1=FitMat[1];
  double Sfree1=FitMat[2];
  double eps=1E-9;
  psigma=vr_pfree_cumulative(bindrad,  r0, tcurr, Dfit1,  bindrad, rcfit1);
  psigma*=Sfree1;
  double phalf=vr_pfree_cumulative(FitMat[6],  r0, tcurr, Dfit1,  bindrad, rcfit1);
  phalf*=Sfree1;
  double ps_numer1=(phalf-psigma);
  double Dfit2=FitMat[3];
  double rcfit2=FitMat[4];
  double Sfree2=FitMat[5];
  double plim2=vr_pfree_cumulative(FitMat[6],  r0, tcurr, Dfit2,  bindrad, rcfit2);
  plim2*=Sfree2;
  double tol=1E-8;
  double plong=vr_pfree_cumulative(Rlong,  r0, tcurr, Dfit2,  bindrad, rcfit2);
  //  if(abs(plong-1)>tol)cout <<"pfree long limit not 1: "<<plong<<endl;
  plong*=Sfree2;
  
  double ps_numer=ps_numer1+(plong-plim2);

  double psurvive=1-passoc;

  psurvive+=eps;
  double ps_ratio=psurvive/ps_numer;
  double prevsum;
  //  cout <<"calculated psurvive from fit: "<<psurvive<<" integral over numer pirr, to Rlong: "<<ps_numer<<" Rlong: "<<Rlong<<" ratio: "<<psurvive/ps_numer<<endl;
  //cout <<"value at sigma to be subtracted: "<<psub<<" Value of cumul at limit r: "<<limit.half<<" is; "<<phalf<<endl;
  /*Since pirr may not integrate to psurvive, rescale it to. */

  prob=rnum;//rnum will be between delp*MAXP-passoc and psurvive.

  j=0;//jprev;
  //cout <<"time start: "<<tmin*1<<endl;
  cumsum=0;
  if(R_table_end<FitMat[6]){
    //extending in first bin
    prevsum=-psigma*ps_ratio;
    Dfit1=FitMat[0];
    rcfit1=FitMat[1];
    Sfree1=FitMat[2];

  }else{
    prevsum=ps_numer1*ps_ratio-plim2*ps_ratio;
    Dfit1=FitMat[3];
    rcfit1=FitMat[4];
    Sfree1=FitMat[5];
    
  }
  double jmax=1E7;
  while(prob>cumsum &&j<jmax){
    r1=rmin*j+R_table_end;
    /*need to use other fit parms, and subtract off sum at previous r1.*/
    term1=vr_pfree_cumulative(r1, r0, tcurr, Dfit1, bindrad, rcfit1);
    cumsum=term1*Sfree1*ps_ratio+prevsum;
      
    /*this is the cumulative probability for a given r value, compare this sum to a random number to 
      generate the separation
      we want the probability to be evenly spaced
      It will sum up to the survival probability, since it survived, the random number will be between p_assoc=1-survival and 1. 
    */
    j++;
  }
  if(j>=jmax)
    cout <<"broke out of final sample loop, pfree_2sfiy. final r:  "<<r1<<" prob: "<<prob<<" final cumulative sum: "<<cumsum<<endl;

  //  cout <<r0<<'\t'<<r1<<'\t'<<" Final prob1? "<<prob<<" Psurvive: "<<psurvive<<" Randonum Prob to match: "<<rnum<<endl;
  return r1;
}