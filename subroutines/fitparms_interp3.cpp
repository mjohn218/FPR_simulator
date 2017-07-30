#include "Vrnumer.h"
#include "numeric_GF.h"
#include "GF_calls.h"

double fitparms_interp3(double **Dmat, double ct1, double ct2, double dih1, double delthet, double deldih, int Nthet, int dihbins )
{

  int cind1=int((ct1+1)/delthet);
  int cind2=int((ct2+1)/delthet);
  int dind=int(dih1/deldih);
  //if(cind1==Nthet-1)cind1=Nthet-2;//if costheta is exactly 1, this can happen
  //if(cind2==Nthet-1)cind2=Nthet-2;//if costheta is exactly 1, this can happen
  //if(dind==dihbins-1)dind=dihbins-2;//if dih is exactly pi, this can happen
  //  cout <<"c1ind: "<<cind1<<" c2ind: "<<cind2<<" dih ind: "<<dind<<endl;
  //interp ct1 for each value of ct2, and dih.
  double f1, f2;
  double slope;
  double ctlow;
  f1=Dmat[dind][cind1*Nthet+cind2];//at y1, z1
  f2=Dmat[dind][(cind1+1)*Nthet+cind2];
  slope=(f2-f1)/delthet;
  ctlow=cind1*delthet-1;
  double fint1=slope*(ct1-ctlow)+f1;
  //   cout <<"F2: "<<f2<<" f1: "<<f1<<" slope: "<<slope<<" cos_theta, low: "<<ctlow<<" f1interp: "<<fint1<<endl;
  f1=Dmat[dind][cind1*Nthet+cind2+1];//at y2, z1
  f2=Dmat[dind][(cind1+1)*Nthet+cind2+1];
  slope=(f2-f1)/delthet;
  double fint2=slope*(ct1-ctlow)+f1;
  // cout <<"F2: "<<f2<<" f1: "<<f1<<" slope: "<<slope<<" cos_theta, low: "<<ctlow<<" f2interp: "<<fint2<<endl;
  slope=(fint2-fint1)/delthet;
  ctlow=cind2*delthet-1;
  double fyint1=slope*(ct2-ctlow)+fint1;
  //cout <<"slope: "<<slope<<" cos_theta, low: "<<ctlow<<" ftotinterp, chi1: "<<fyint1<<endl;
  
  /*now repeat for the upper chi value*/
  f1=Dmat[dind+1][cind1*Nthet+cind2];//at y1, z1
  f2=Dmat[dind+1][(cind1+1)*Nthet+cind2];
  slope=(f2-f1)/delthet;
  ctlow=cind1*delthet-1;
  fint1=slope*(ct1-ctlow)+f1;
  f1=Dmat[dind+1][cind1*Nthet+cind2+1];//at y2, z1
  f2=Dmat[dind+1][(cind1+1)*Nthet+cind2+1];
  slope=(f2-f1)/delthet;
  fint2=slope*(ct1-ctlow)+f1;
  
  slope=(fint2-fint1)/delthet;
  ctlow=cind2*delthet-1;
  double fyint2=slope*(ct2-ctlow)+fint1;
  
  /*Now interp between these two values */
  slope=(fyint2-fyint1)/deldih;
  double dlow=dind*deldih;
  double Dinterp=slope*(dih1-dlow)+fyint1;

  return Dinterp;
}
