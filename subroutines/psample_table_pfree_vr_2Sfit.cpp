#include "Vrnumer.h"
#include "GF_calls.h"
#include "numeric_GF.h"
#include "rand_gsl.h"

void psample_table_pfree_vr_2Sfit(double rmin, double r0, double tcurr, double *FitMat, double bindrad,  double *table, int MAXP, double delp, double passoc, double Rlong)
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
  int i1=i;
  //cout <<"first index: "<<i1<<" for reflecting should be zero "<<kact<<endl; 
  
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
  double plong=vr_pfree_cumulative(Rlong,  r0, tcurr, Dfit2,  bindrad, rcfit2);//plong not to 1, but to pnorm
  //  if(abs(plong-1)>tol)cout <<"pfree long limit not 1: "<<plong<<endl;
  plong*=Sfree2;
  double ps_numer=ps_numer1+(plong-plim2);
  
  delp=psurvive/(1.0*MAXP);
  psurvive+=eps;
  double ps_ratio=psurvive/ps_numer;
  cout <<"calculated psurvive from fit: "<<psurvive<<" integral over numer pirr, to Rlong: "<<ps_numer<<" Rlong: "<<Rlong<<" ratio: "<<psurvive/ps_numer<<endl;
  double prevsum;
  // cout <<"value at sigma to be subtracted: "<<psub<<endl;
  for(i=i1;i<MAXP-1;i++){
    //now find which time this corresponds to
    prob=delp*(i+1);//at final value of i, prob will be equal=1
    j=jprev;
    //cout <<"time start: "<<tmin*1<<endl;
    cumsum=0;
    while(prob>cumsum){
      r1=rmin*j+bindrad;
      if(r1<FitMat[6]){
	term1=vr_pfree_cumulative(r1, r0, tcurr, Dfit1, bindrad, rcfit1);
	cumsum=(term1*Sfree1-psigma)*ps_ratio;
	prevsum=cumsum;
      }else{
	term1=vr_pfree_cumulative(r1, r0, tcurr, Dfit2, bindrad, rcfit2);
	cumsum=(term1*Sfree2-plim2)*ps_ratio+prevsum;
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
    //cout <<r1<<'\t'<<prob-probstart<<endl;
  }
  cout <<r0<<'\t'<<r1<<'\t'<<" Final prob1? "<<prob<<endl;
  //last value
//   prob=delp*(i+1)-1E-8;//at final value of i, prob will be equal=1
//   j=jprev;
//   double pfree=1;
//   double small=1E-16;
//   cumsum=0;
//   while(prob>cumsum && pfree>small){
//     r1=rmin*(j+1)+bindrad;
//     term1=vr_pfree_cumulative(r1, r0, tcurr, Dfit2, bindrad, rcfit2);
//     cumsum=(term1*Sfree2-plim2)*ps_ratio+prevsum;
    
//     pfree=vr_pfree_value_norm(r1,  r0,  tcurr,  Dfit2,  bindrad,  rcfit2);
//     pfree*=Sfree2;
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
