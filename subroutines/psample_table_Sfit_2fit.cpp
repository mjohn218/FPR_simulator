#include "Vrnumer.h"
#include "GF_calls.h"
#include "numeric_GF.h"
#include "rand_gsl.h"

void psample_table_Sfit_2fit(double rmin, double r0, double tcurr, double bindrad, double *table, int MAXP, double delp, double passoc, Vrnumer &limit, double *FitMat, double Rlong)
{
  /*To sample from irreversible distribution, kact is finite, and alpha=(1+kact/kdiff)*sqrt(D)/sigma
    To sample from reflecting distribution, kact=0 and alpha=sqrt(D)/sigma
  */
  double Dfit1=FitMat[0];
  double kfit1=FitMat[1];
  double Sfit1=FitMat[2];
  double kdiff=4.0*M_PI*Dfit1*bindrad;
  double fact=1.0+kfit1/kdiff;
  double alpha=fact*sqrt(Dfit1)/bindrad;
  double alpha2;
  
  double r1;
  double cumsum;
  

  int j, jprev;
  jprev=0;
  
  //or just skip all the probs less than pstart, so you don't need to know what p_survive is 
  int i=0;
  double prob=0;
  table[0]=bindrad;
  int i1=i;
  //cout <<"first index: "<<i1<<" for reflecting should be zero "<<kfit<<endl; 
  double psub;
  psub=peval_cumulativeF(bindrad,  r0,  tcurr,  Dfit1,  bindrad,  alpha,  kfit1);
  psub*=Sfit1;//FitMat[2];
  double plim1=peval_cumulativeF(FitMat[6],  r0,  tcurr,  Dfit1,  bindrad,  alpha,  kfit1);
  plim1*=Sfit1;//FitMat[2];
  double ps_numer=(plim1-psub);//first half
  double pcume, prevsum;
  double Dfit2, kfit2;
  Dfit2=FitMat[3];
  double kdiff_half=4*M_PI*Dfit2*bindrad;
  kfit2=FitMat[4];
  double alpha_half=(1+kfit2/kdiff_half)*sqrt(Dfit2)/bindrad;
	
  double phalf=peval_cumulativeF(FitMat[6], r0,  tcurr,  FitMat[3],  bindrad,  alpha_half,  FitMat[4]);
  phalf*=FitMat[5];//this is Sfit for second half.
  
  double plong=peval_cumulativeF(Rlong, r0,  tcurr,  FitMat[3],  bindrad,  alpha_half,  FitMat[4]);
  plong*=FitMat[5];//this is Sfit for second half.
  ps_numer+=(plong-phalf);
  double psurvive=1-passoc;
  double eps=1E-8;
  cout <<"calculated psurvive from fit: "<<psurvive<<" integral over numer pirr, to Rlong: "<<ps_numer<<" Rlong: "<<Rlong<<" ratio: "<<psurvive/ps_numer<<endl;
  double ps_ratio=(psurvive+eps)/ps_numer;
  
  cout <<"value at sigma to be subtracted: "<<psub<<" Value of cumul at limit r: "<<FitMat[6]<<" is; "<<phalf<<endl;
  /*Since pirr may not integrate to psurvive, rescale it to. */

  Dfit2=FitMat[3];
  double kdiff2=4*M_PI*Dfit2*bindrad;
  kfit2=FitMat[4];
  alpha2=(1+kfit2/kdiff2)*sqrt(Dfit2)/bindrad;
	
  delp=psurvive/(1.0*MAXP);
  for(i=i1;i<MAXP-1;i++){
    //now find which time this corresponds to
    prob=delp*(i+1);//at final value of i, prob will be equal=1
    j=jprev;
    //cout <<"time start: "<<tmin*1<<endl;
    cumsum=0;
    while(prob>cumsum){
      r1=rmin*j+bindrad;
      if(r1<FitMat[6]){
	
	pcume=peval_cumulativeF(r1,  r0,  tcurr,  Dfit1,  bindrad,  alpha,  kfit1);
	cumsum=Sfit1*pcume*ps_ratio-psub*ps_ratio;
	prevsum=cumsum;
      }else{
	/*need to use other fit parms, and subtract off sum at previous r1.*/
	
	pcume=peval_cumulativeF(r1, r0,  tcurr,  Dfit2,  bindrad,  alpha2,  kfit2);
	cumsum=FitMat[5]*pcume*ps_ratio-phalf*ps_ratio+prevsum;
      }
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
    //    cout <<r1<<'\t'<<prob-probstart<<endl;
  }
  cout <<r0<<'\t'<<r1<<'\t'<<" Final prob1? "<<prob<<" Psurvive: "<<psurvive<<" rmin at end: "<<rmin<<endl;
  
  /*We don't need the last bin anymore because after final r1 we sample directly from 
    cumulative, without interpolating.
  */
  //last value
  // prob=delp*(i+1)-1E-8;//at final value of i, prob will be equal=1
//   j=jprev;
  
//   double small=1E-16;
//   double pirr1=1;
//   cumsum=0;
//   while(prob>cumsum && pirr1>small){
//     r1=rmin*(j+1)+bindrad;
//     /*need to use other fit parms, and subtract off sum at previous r1.*/
  
    
//     pcume=peval_cumulativeF(r1, r0,  tcurr,  Dfit2,  bindrad,  alpha2,  kfit2);
//     cumsum=FitMat[5]*pcume*ps_ratio-phalf*ps_ratio+prevsum;

//     pirr1= pirrev_valueF(r1, r0,tcurr,  Dfit2,  bindrad,  alpha2);
//     pirr1*=FitMat[5]*ps_ratio;
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
//   cout<<"Final r to get to prob1: "<<r1<<'\t'<<prob<<endl;

}
