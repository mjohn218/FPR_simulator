#include "Vrnumer.h"
#include "numeric_GF.h"
#include "GF_calls.h"

double lookup_pirr_pfree_ratio(double R1, double prevsep, double deltat, double ps_prev, double rtol, double rc, double Dtot, double bindrad, double **FitMat, Vrnumer &limit, double passoc)
{
  
  double pfree_vr=vr_pfree_value_norm(R1, prevsep, deltat,  Dtot,  bindrad,  rc);//this is exact free propagator with potential, renormalized to 1 (excludes r<sigma).
  double pirr=lookup_pirr( R1,  prevsep,  deltat, rc,  Dtot,  bindrad,  FitMat, limit, passoc );
  if(pirr<0){
    pirr=lookup_pirr( R1,  prevsep+1E-8,  deltat, rc,  Dtot,  bindrad,  FitMat, limit, passoc );
    cout <<"NEGATIVE PIRR, NEW VALUE: "<<pirr<<endl;
  }
  double ratio=pirr/(pfree_vr*ps_prev);
  //cout <<"Pfree: "<<pfree_vr<<" at R1: "<<R1<<" Rprev: "<<prevsep<<" ratio: "<<ratio<<" ps_prev: "<<ps_prev<<endl;
  return ratio;

}
