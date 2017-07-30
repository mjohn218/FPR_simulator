#include "Vrnumer.h"
#include "GF_calls.h"
#include "numeric_GF.h"
#include "rand_gsl.h"

double psample_final_Sfit_2fit(double rmin, double r0, double tcurr,  double bindrad, double passoc, Vrnumer &limit, double *FitMat, double Rlong, double R_table_end, double rnum)
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
  double ps_numer=(plim1-psub);//first half
  double first_half=ps_numer;
  double pcume, prevsum;
  double Dfit, kfit;
  Dfit=FitMat[3];
  double kdiff_half=4*M_PI*Dfit*bindrad;
  kfit=FitMat[4];
  double alpha_half=(1+kfit/kdiff_half)*sqrt(Dfit)/bindrad;
	
  double phalf=peval_cumulativeF(FitMat[6], r0,  tcurr,  FitMat[3],  bindrad,  alpha_half,  FitMat[4]);
  phalf*=FitMat[5];//this is Sfit for second half.
  
  double plong=peval_cumulativeF(Rlong, r0,  tcurr,  FitMat[3],  bindrad,  alpha_half,  FitMat[4]);
  plong*=FitMat[5];//this is Sfit for second half.
  ps_numer+=(plong-phalf);
  double psurvive=1-passoc;
  //cout <<"calculated psurvive from fit: "<<psurvive<<" integral over numer pirr, to Rlong: "<<ps_numer<<" Rlong: "<<Rlong<<" ratio: "<<psurvive/ps_numer<<endl;
  double eps=1E-8;
  double ps_ratio=(psurvive+eps)/ps_numer;
  
  //cout <<"value at sigma to be subtracted: "<<psub<<" Value of cumul at limit r: "<<limit.half<<" is; "<<phalf<<endl;
  /*Since pirr may not integrate to psurvive, rescale it to. */

  prob=rnum;//rnum will be between delp*MAXP-passoc and psurvive.

  j=0;//jprev;
  //cout <<"time start: "<<tmin*1<<endl;
  cumsum=0;
  if(R_table_end<FitMat[6]){
    //extending in first bin
    prevsum=-psub*ps_ratio;
    Dfit=FitMat[0];
    kfit=FitMat[1];
    Sfit=FitMat[2];
    kdiff=4.0*M_PI*Dfit*bindrad;
    alpha=(1+kfit/kdiff)*sqrt(Dfit)/bindrad;
  }else{
    prevsum=first_half*ps_ratio-phalf*ps_ratio;
    Dfit=FitMat[3];
    kfit=FitMat[4];
    Sfit=FitMat[5];
    kdiff=4.0*M_PI*Dfit*bindrad;
    alpha=(1+kfit/kdiff)*sqrt(Dfit)/bindrad;
    
    // r1=rmin*j+R_table_end;
//     pcume=peval_cumulativeF(r1, r0,  tcurr,  Dfit,  bindrad,  alpha,  kfit);
//     cumsum=Sfit*pcume*ps_ratio+prevsum;
//     if(psurvive-cumsum>0.0001){
//       cout <<"at rend: "<<R_table_end <<" starting at low pcume! : "<<prevsum<<" first pcume: "<<cumsum<<" psurvive: "<<psurvive<<" r1: "<<r1<<" r0: "<<r0<<'\t';
//       cout <<"dfit: "<<Dfit<<" kfit: "<<kfit<<" sfit: "<<Sfit<<endl;
      
//     }
    
  }
  double jmax=1E7;
  while(prob>cumsum&&j<jmax){
    r1=rmin*j+R_table_end;
    /*need to use other fit parms, and subtract off sum at previous r1.*/
        
    pcume=peval_cumulativeF(r1, r0,  tcurr,  Dfit,  bindrad,  alpha,  kfit);
    cumsum=Sfit*pcume*ps_ratio+prevsum;
      
    /*this is the cumulative probability for a given r value, compare this sum to a random number to 
      generate the separation
      we want the probability to be evenly spaced
      It will sum up to the survival probability, since it survived, the random number will be between p_assoc=1-survival and 1. 
    */
    j++;
  }
  if(j>=jmax)
    cout <<"broke out of final sample loop, pirr_2sfi. final r:  "<<r1<<" prob: "<<prob<<" final cumulative sum: "<<cumsum<<endl;

  //  cout <<r0<<'\t'<<r1<<'\t'<<" Final prob1? "<<prob<<" Psurvive: "<<psurvive<<" Randonum Prob to match: "<<rnum<<endl;
  return r1;
}
