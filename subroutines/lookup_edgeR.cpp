#include "Vrnumer.h"
#include "GF_calls.h"
#include "numeric_GF.h"
#include "rand_gsl.h"

double lookup_edgeR(double rmin, double r0, int r0ind, double deltat, double Dtot, double bindrad,double Dfit_ps, double alpha_ps,double cof_ps, Vrnumer &limit, double **FitMat, double R_table_end, double rnum )
{

  double R;
  double Rlong=r0+5*sqrt(6*Dtot*deltat);
  double passoc=survive_irrF( r0, deltat,  Dfit_ps,  bindrad,  alpha_ps,  cof_ps);//survival prob is always from the Dfit_ps, kfit_ps parms 
  
  //cout <<"R0: "<<r0<<" r0ind: "<<r0ind<<endl;
  if(r0<limit.r0_2hybridstart){
    R=psample_final_Sfit_2fit( rmin,  r0, deltat,  bindrad,  passoc, limit,  FitMat[r0ind],  Rlong,  R_table_end, rnum);
  }else if(r0<limit.r0_2freefitstart){
    R=psample_final_hybrid( rmin,  r0, deltat,  bindrad,  passoc, limit,  FitMat[r0ind],  Rlong,  R_table_end, rnum);
  }else if(r0<limit.r0_freeonlystart){
    R=psample_final_pfree_2Sfit( rmin,  r0, deltat,  bindrad,  passoc, limit,  FitMat[r0ind],  Rlong,  R_table_end, rnum);
  }else{
    R=psample_final_pfree_vr( rmin,  r0, deltat, Dtot,  bindrad, limit.rc,  passoc, R_table_end, rnum);
  }
  return R;
}
