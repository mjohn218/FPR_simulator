#include "Vrnumer.h"
#include "numeric_GF.h"
#include "GF_calls.h"

void load_fitparms(ifstream &infile, double **Dmat, int &Nthet, int &dihbins)
{
  int i, k, j;
  infile >>dihbins >>Nthet;
  for(i=0;i<dihbins;i++)
    Dmat[i]=(double *)realloc(Dmat[i], sizeof(double) *Nthet*Nthet);//new double[Nthet_fit*Nthet_fit];
  double Derr;
  for(i=0;i<dihbins;i++){
    cout <<"new dihedral "<<endl;
    for(k=0;k<Nthet;k++){
      for(j=0;j<Nthet;j++){
	infile >>Dmat[i][k*Nthet+j]>>Derr;
	cout <<Dmat[i][k*Nthet+j]<<'\t';
      }
    }
    cout <<endl;
  }


}
