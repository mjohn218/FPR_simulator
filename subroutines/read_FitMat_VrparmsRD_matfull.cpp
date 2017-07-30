#include "Vrnumer.h"
#include "numeric_GF.h"
#include "GF_calls.h"

void read_FitMat_VrparmsRD_matfull(double **FitMat, ifstream &infile, Vrnumer &limit, int ncol)
{
  int i, j;

  infile >>limit.rc;
  infile.ignore(400,'\n');
  infile >>limit.Dmat_ps[0] >>limit.Dmat_ps[1];
  infile.ignore(400,'\n');
  infile >>limit.kmat_ps[0]>>limit.kmat_ps[1];
  infile.ignore(400,'\n');
  infile >>limit.fit1;
  infile.ignore(400,'\n');

  
  infile >>limit.r0bins;
  infile.ignore(400,'\n');
  infile >>limit.delr0;
  
  infile.ignore(400,'\n');
  infile >>limit.r0_2hybridstart;
  infile.ignore(400,'\n');
  infile >>limit.r0_2freefitstart;
  infile.ignore(400,'\n');
  infile >>limit.r0_freeonlystart;
  infile.ignore(400,'\n');

  cout <<"Freefit limit: "<<limit.r0_2freefitstart<<" freefit end: "<<limit.r0_freeonlystart<<" r0bins: "<<limit.r0bins<<endl;
  double r0;
  for(i=0;i<limit.r0bins;i++){
    infile>>r0;
    for(j=0;j<ncol;j++)
      infile >>FitMat[i][j];
    cout <<"Fitmat[0]: "<<FitMat[i][0]<<endl;
  }
  

}
