#include "Vrnumer.h"
#include "numeric_GF.h"
#include "rand_gsl.h"

double get_R_from_table_2D_even(int MAXP, double rcurr, int r0ind, double **table_irr, double Rlong, double bindrad, double passoc, double rnum, int TABLE_LIM, double delr0, gsl_matrix *pcume, double RstepSize, double delp)
{

  double rmin;
  double slope;
  double Rdown, Rup;
  double R, R_table_end;

  double rnum_rescale=(rnum-passoc);///(1.0-passoc)*(1.0-passocvec[r0ind]);
  //double rnumup=(rnum-passoc);///(1.0-passoc)*(1.0-passocvec[r0ind+1]);
  //delpdown=delpvec[r0ind];//(1.0-passocvec[r0ind])/(1.0*MAXP);
  int pind=int(rnum_rescale/delp);
  //delpup=delpvec[r0ind+1];//(1.0-passocvec[r0ind+1])/(1.0*MAXP);
  //pindup=int(rnumup/delpup);
  if(pind==TABLE_LIM){
    rnum_rescale-=1E-9;
    //rnumup-=1E-9;
  }
  
  if(pind==TABLE_LIM || pind==0){
    //eed to use exact cumulative dist to get R.
    //cout <<"Rnum_Rescale: "<<rnum_rescale<<" rnum_rescale-probvec: "<<rnum_rescale-probvec[1][0]<<" probvec: "<<probvec[1][0]<<" delp: "<<delp<<" pind:: "<<pind<<" r0: "<<R1<<" R0_ind: "<<r0ind<<" rep: "<<rep<<endl;
    if(pind==0)rmin=1E-7;
    else rmin=1E-6;
    R_table_end=table_irr[r0ind][pind];
    //Rlong=rcurr+5*sqrt(6.0*Dtot*deltat);
    //    Rdown=psample_final_2D( rmin,  r0, double bindrad, double passoc, double Rlong, double R_table_end, double rnum_rescale, gsl_matrix *pcume, double RstepSize);
    Rdown=psample_final_2D( rmin,  r0ind*delr0+bindrad, bindrad,  passoc,   Rlong,  R_table_end, rnum_rescale, pcume, RstepSize);//lookup_edgeR(rmin,  r0ind*delr0+bindrad, r0ind, deltat,  Dtot,  bindrad,  limit.Dmat_ps[1],  alpha2,  cof2, limit,  FitMat,  R_table_end, rnum_rescale-pvec);
    Rup=psample_final_2D( rmin, (r0ind+1)*delr0+bindrad, bindrad, passoc, Rlong, R_table_end, rnum_rescale, pcume, RstepSize);//lookup_edgeR(rmin,  (r0ind+1)*delr0+bindrad, r0ind+1, deltat,  Dtot,  bindrad,  limit.Dmat_ps[1],  alpha2,  cof2, limit,  FitMat,  R_table_end, rnum_rescale-pvec);
    slope=(Rup-Rdown)/delr0;
    R=Rdown+slope*(rcurr-(r0ind*delr0+bindrad));
    
    if(R>Rlong){
      cout <<" Reject R as too long: "<<R<<" Rlong: "<<Rlong<<endl;
      rnum=rand_gsl()*1.0;
      while(rnum<passoc){
	rnum=rand_gsl()*1.0;
      }
      //delp=(1.0-pvec)/(1.0*MAXP);
      rnum_rescale=rnum-passoc;
      pind=int((rnum_rescale)/delp);
      
      Rdown=lookup_psample(rnum_rescale, pind, delp, table_irr[r0ind]);//using current pos, you can calculate survival prob 
      Rup=lookup_psample(rnum_rescale, pind, delp, table_irr[r0ind+1]);//using current pos, you can calculate survival prob 
      slope=(Rup-Rdown)/delr0;
      R=Rdown+slope*(rcurr-(r0ind*delr0+bindrad));
    }
  }else{
    /*choose a new separation given your current separation from the irreversible propagator*/
    //delp=(1.0-pvec)/(1.0*MAXP);
    Rdown=lookup_psample(rnum_rescale, pind, delp, table_irr[r0ind]);//using current pos, you can calculate survival prob 
    Rup=lookup_psample(rnum_rescale, pind, delp, table_irr[r0ind+1]);//using current pos, you can calculate survival prob 
    slope=(Rup-Rdown)/delr0;
    R=Rdown+slope*(rcurr-(r0ind*delr0+bindrad));
  }
  return R;
}
