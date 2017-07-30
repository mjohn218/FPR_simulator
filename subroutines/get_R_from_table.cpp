#include "Vrnumer.h"
#include "GF_calls.h"
#include "numeric_GF.h"
#include "rand_gsl.h"

double get_R_from_table(int pind, int MAXP, double rcurr, double delp, int r0ind, double **table_irr, double Dtot,double deltat, double bindrad, Vrnumer limit, double **FitMat, double pvec, double alpha1, double alpha2, double cof1, double cof2, double rnum, int TABLE_LIM, double delr0)
{
  double Rlong;
  double rmin;
  double slope;
  double Rdown, Rup;
  double R, R_table_end;


  if(pind==TABLE_LIM || pind==0){
    //eed to use exact cumulative dist to get R.
    //cout <<"Rnum: "<<rnum<<" rnum-probvec: "<<rnum-probvec[1][0]<<" probvec: "<<probvec[1][0]<<" delp: "<<delp<<" pind:: "<<pind<<" r0: "<<R1<<" R0_ind: "<<r0ind<<" rep: "<<rep<<endl;
    rmin=1E-4;
    R_table_end=table_irr[r0ind][pind];
    Rlong=rcurr+5*sqrt(6.0*Dtot*deltat);
    if(rcurr<limit.fit1){
      Rdown=lookup_edgeR(rmin,  r0ind*delr0+bindrad, r0ind, deltat,  Dtot,  bindrad,  limit.Dmat_ps[0],  alpha1,  cof1, limit,  FitMat,  R_table_end, rnum-pvec);
      Rup=lookup_edgeR(rmin,  (r0ind+1)*delr0+bindrad, r0ind+1, deltat,  Dtot,  bindrad,  limit.Dmat_ps[0],  alpha1,  cof1, limit,  FitMat,  R_table_end, rnum-pvec);
      slope=(Rup-Rdown)/delr0;
      R=Rdown+slope*(rcurr-(r0ind*delr0+bindrad));
    }else{
      Rdown=lookup_edgeR(rmin,  r0ind*delr0+bindrad, r0ind, deltat,  Dtot,  bindrad,  limit.Dmat_ps[1],  alpha2,  cof2, limit,  FitMat,  R_table_end, rnum-pvec);
      Rup=lookup_edgeR(rmin,  (r0ind+1)*delr0+bindrad, r0ind+1, deltat,  Dtot,  bindrad,  limit.Dmat_ps[1],  alpha2,  cof2, limit,  FitMat,  R_table_end, rnum-pvec);
      slope=(Rup-Rdown)/delr0;
      R=Rdown+slope*(rcurr-(r0ind*delr0+bindrad));
    }
    if(R>Rlong){
      cout <<" Reject R as too long: "<<R<<" Rlong: "<<Rlong<<endl;
      rnum=rand_gsl()*1.0;
      while(rnum<pvec){
	rnum=rand_gsl()*1.0;
      }
      delp=(1.0-pvec)/(1.0*MAXP);
      pind=int((rnum-pvec)/delp);
      Rdown=lookup_psample(rnum-pvec, pind, delp, table_irr[r0ind]);//using current pos, you can calculate survival prob 
      Rup=lookup_psample(rnum-pvec, pind, delp, table_irr[r0ind+1]);//using current pos, you can calculate survival prob 
      slope=(Rup-Rdown)/delr0;
      R=Rdown+slope*(rcurr-(r0ind*delr0+bindrad));
    }
  }else{
    /*choose a new separation given your current separation from the irreversible propagator*/
    delp=(1.0-pvec)/(1.0*MAXP);
    Rdown=lookup_psample(rnum-pvec, pind, delp, table_irr[r0ind]);//using current pos, you can calculate survival prob 
    Rup=lookup_psample(rnum-pvec, pind, delp, table_irr[r0ind+1]);//using current pos, you can calculate survival prob 
    slope=(Rup-Rdown)/delr0;
    R=Rdown+slope*(rcurr-(r0ind*delr0+bindrad));
  }
  return R;
}
