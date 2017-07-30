#include "Vrnumer.h"
#include "GF_calls.h"
#include "numeric_GF.h"
#include "rand_gsl.h"

double psample_final_hybrid(double rmin, double r0, double tcurr,  double bindrad, double passoc, Vrnumer &limit, double *FitMat, double Rlong, double R_table_end, double rnum)
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
  double Dtot=FitMat[0];
  double kact=FitMat[1];
  double Sfit=FitMat[2];
  
  double kdiff=4.0*M_PI*Dtot*bindrad;
  double fact=1.0+kact/kdiff;
  double alpha=fact*sqrt(Dtot)/bindrad;
  
  psub=peval_cumulativeF(bindrad,  r0,  tcurr,  Dtot,  bindrad,  alpha,  kact);
  psub*=Sfit;//FitMat[2];
  double plim1=peval_cumulativeF(FitMat[6],  r0,  tcurr,  Dtot,  bindrad,  alpha,  kact);
  plim1*=Sfit;//FitMat[2];
  double ps_numer1=(plim1-psub);//first half
  double first_half=ps_numer1;

  double Dfit2=FitMat[3];
  double rcfit2=FitMat[4];
  double Sfree2=FitMat[5];
  double plim2=vr_pfree_cumulative(FitMat[6],  r0, tcurr, Dfit2,  bindrad, rcfit2);
  plim2*=Sfree2;
  double plong=vr_pfree_cumulative(Rlong,  r0, tcurr, Dfit2,  bindrad, rcfit2);//plong will not be 1 because it's divided by pnorm. and pnorm is only for integral from sigma to inf.
  double tol=1E-8;
  //  if(abs(plong-1)>tol)cout <<"pfree long limit not 1: "<<plong<<endl;
  plong*=Sfree2;
  
  double ps_numer=ps_numer1+(plong-plim2);
  double psurvive=1-passoc;

  double eps=1E-9;
  psurvive+=eps;
  double ps_ratio=psurvive/ps_numer;
  double prevsum;
  //cout <<"calculated psurvive from fit: "<<psurvive<<" integral over numer pirr, to Rlong: "<<ps_numer<<" Rlong: "<<Rlong<<" ratio: "<<psurvive/ps_numer<<endl;
  //cout <<"value at sigma to be subtracted: "<<psub<<" Value of cumul at limit r: "<<limit.half<<" is; "<<phalf<<endl;
  /*Since pirr may not integrate to psurvive, rescale it to. */

  prob=rnum;//rnum will be between delp*MAXP-passoc and psurvive.
  double pcume;
  j=0;//jprev;
  //cout <<"time start: "<<tmin*1<<endl;
  cumsum=0;
  double jmax=1E7;
  if(R_table_end<FitMat[6]){
    //extending in first bin
    prevsum=-psub*ps_ratio;
    Dtot=FitMat[0];
    kact=FitMat[1];
    Sfit=FitMat[2];
    //pirr
    
    while(prob>cumsum&&j<jmax){
      r1=rmin*j+R_table_end;
      /*need to use other fit parms, and subtract off sum at previous r1.*/
      pcume=peval_cumulativeF(r1, r0,  tcurr,  Dtot,  bindrad,  alpha,  kact);
      cumsum=Sfit*pcume*ps_ratio+prevsum;
    
      
      /*this is the cumulative probability for a given r value, compare this sum to a random number to 
	generate the separation
	we want the probability to be evenly spaced
	It will sum up to the survival probability, since it survived, the random number will be between p_assoc=1-survival and 1. 
      */
      j++;
    }
    
  }else{
    prevsum=ps_numer1*ps_ratio-plim2*ps_ratio;
    Dfit2=FitMat[3];
    rcfit2=FitMat[4];
    Sfree2=FitMat[5];
    while(prob>cumsum&&j<jmax){
      r1=rmin*j+R_table_end;
      /*need to use other fit parms, and subtract off sum at previous r1.*/
      term1=vr_pfree_cumulative(r1, r0, tcurr, Dfit2, bindrad, rcfit2);
      cumsum=term1*Sfree2*ps_ratio+prevsum;
      
      /*this is the cumulative probability for a given r value, compare this sum to a random number to 
	generate the separation
	we want the probability to be evenly spaced
	It will sum up to the survival probability, since it survived, the random number will be between p_assoc=1-survival and 1. 
      */
      j++;
    }
    
  }
  if(j>=jmax)
    cout <<"broke out of final sample loop, phybrid. final r:  "<<r1<<" prob: "<<prob<<" final cumulative sum: "<<cumsum<<endl;

  //  cout <<r0<<'\t'<<r1<<'\t'<<" Final prob1? "<<prob<<" Psurvive: "<<psurvive<<" Randonum Prob to match: "<<rnum<<endl;
  return r1;
}
