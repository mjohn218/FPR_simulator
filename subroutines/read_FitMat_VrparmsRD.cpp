#include "Vrnumer.h"
#include "numeric_GF.h"
#include "GF_calls.h"

void read_FitMat_VrparmsRD(double **FitMat, ifstream &infile, Vrnumer &limit)
{
  int i, j;
  int ncol=6;
  infile >>limit.rc;
  infile.ignore(400,'\n');
  infile >>limit.Dfit_ps;
  infile.ignore(400,'\n');
  infile >>limit.kfit_ps;
  infile.ignore(400,'\n');
  
  infile >>limit.r0bins;
  infile.ignore(400,'\n');
  infile >>limit.delr0;
  infile.ignore(400,'\n');
  infile >>limit.r0_2fit;
  infile.ignore(400,'\n');
  //  infile >>limit.half;
  infile.ignore(400,'\n');
  infile >>limit.r0_onefit;
  infile.ignore(400,'\n');
  //  infile >>limit.persist;
  infile.ignore(400,'\n');
  infile >>limit.freefit;
  infile.ignore(400,'\n');
  cout <<"Freefit limit: "<<limit.freefit<<" r0bins: "<<limit.r0bins<<endl;
  double r0;
  for(i=0;i<limit.r0bins;i++){
    infile>>r0;
    for(j=0;j<ncol;j++)
      infile >>FitMat[i][j];
    cout <<"Fitmat[0]: "<<FitMat[i][0]<<endl;
  }
  

}
