#include "Vrnumer.h"
#include "numeric_GF.h"

double pirr_simp(gsl_matrix* mpir, double RstepSize, double r, double r0, double a) 
{

    int indexr0 = int(floor((r0 - a) / RstepSize));
  int indexr = int(floor((r - a) / RstepSize));

  if (indexr < 0) {
    indexr = 0;
  }
  if (indexr0 < 0) {
    indexr0 = 0;
  }
  /*Interp for r0*/
	
  //double Ar0 = indexr0*RstepSize+a;//gsl_matrix_get(msur, 0, indexr0);
  //double Ar = indexr*RstepSize+a;//gsl_matrix_get(msur, 0, indexr);
  double valA = gsl_matrix_get(mpir, indexr0, indexr);
  //double distA = sqrt(pow(Ar0 - r0, 2) + pow(Ar - r, 2));

  //double Br0 = indexr0*RstepSize+a;//gsl_matrix_get(msur, 0, indexr0);
  //double Br = (indexr+1)*RstepSize+a;//gsl_matrix_get(msur, 0, indexr + 1);
  double valB = gsl_matrix_get(mpir, indexr0, indexr + 1);
  //double distB = sqrt(pow(Br0 - r0, 2) + pow(Br - r, 2));
  double slope=(valB-valA)/RstepSize;
  double delr=r-indexr*RstepSize-a;
  double fint1=valA+slope*delr;
  //  cout <<"r1: "<<r<<" rlow: "<<indexr*RstepSize+a<<" rhigh: "<<(indexr+1)*RstepSize+a<<" f at rlow: "<<valA<<" f at rhigh "<<valB<<" final f: "<<fint1<<" slope: "<<slope<<" delr: "<<delr<<endl;

  //double Cr0 = (indexr0+1)*RstepSize+a;//gsl_matrix_get(msur, 0, indexr0 + 1);
  //double Cr = (indexr+1)*RstepSize+a;//gsl_matrix_get(msur, 0, indexr + 1);
  double valC = gsl_matrix_get(mpir, indexr0 + 1, indexr + 1);
  //double distC = sqrt(pow(Cr0 - r0, 2) + pow(Cr - r, 2));

  //double Dr0 = (indexr0+1)*RstepSize+a;//gsl_matrix_get(msur, 0, indexr0 + 1);
  //double Dr = (indexr)*RstepSize+a;//gsl_matrix_get(msur, 0, indexr);
  double valD = gsl_matrix_get(mpir, indexr0 + 1, indexr);
  //double distD = sqrt(pow(Dr0 - r0, 2) + pow(Dr - r, 2));
  slope=(valC-valD)/RstepSize;
  delr=r-indexr*RstepSize-a;
  
  double fint2=valD+slope*delr;
  //cout <<"r1: "<<r<<" rlow: "<<indexr*RstepSize+a<<" rhigh: "<<(indexr+1)*RstepSize+a<<" f at rlow: "<<valD<<" f at rhigh "<<valC<<" final f: "<<fint2<<endl;
  //double prob = ((valA / distA) + (valB / distB) + (valC / distC) + (valD / distD)) / ((1 / distA) + (1 / distB) + (1 / distC) + (1 / distD));
  slope=(fint2-fint1)/RstepSize;
  delr=r0-indexr0*RstepSize-a;
  double prob=fint1+slope*delr;
  //  cout <<"r0: "<<r0<<" r0 low: "<<indexr0*RstepSize+a<<" r0high: "<<(indexr0+1)*RstepSize+a<<" f at r0: "<<fint1<<" f at r0up: "<<fint2<<" final f: "<<prob<<endl;
  return prob;
  
}
