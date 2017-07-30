#include "Vrnumer.h"
#include "numeric_GF.h"
#include "rand_gsl.h"

double get_R_from_table_2D(int MAXP, double rcurr, int r0ind, double **table_irr, double Rlong, double bindrad, double passoc_sample, double rnum, int TABLE_LIM, double delr0, gsl_matrix *pcume, double RstepSize, double *passocvec, double *delpvec)
{

  double rmin;
  double slope;
  double Rdown, Rup;
  double R, R_table_end;
  /*Because the value of passoc at the exact r0 value can vary from that at the actual bins, 
    need to rescale the random number to be over the range of that specific r0 bin value.
  */
  double rnumdown=(rnum-passoc_sample)/(1.0-passoc_sample)*(1.0-passocvec[r0ind]);
  double rnumup=(rnum-passoc_sample)/(1.0-passoc_sample)*(1.0-passocvec[r0ind+1]);
  double delpdown=delpvec[r0ind];//(1.0-passocvec[r0ind])/(1.0*MAXP);
  int pinddown=int(rnumdown/delpdown);
  double delpup=delpvec[r0ind+1];//(1.0-passocvec[r0ind+1])/(1.0*MAXP);
  int pindup=int(rnumup/delpup);
  if(pinddown==TABLE_LIM){
    rnumdown-=1E-9;
    rnumup-=1E-9;
  }
  if(pinddown>TABLE_LIM)cout <<"pind down too high! "<<pinddown<<endl;
  if(pindup>TABLE_LIM)cout <<"pind up too high! "<<pindup<<endl;
  
  if(pinddown==TABLE_LIM || pinddown==0){
    //eed to use exact cumulative dist to get R.
    //cout <<"Rnum: "<<rnum<<" rnum-probvec: "<<rnum-probvec[1][0]<<" probvec: "<<probvec[1][0]<<" delp: "<<delp<<" pind:: "<<pind<<" r0: "<<R1<<" R0_ind: "<<r0ind<<" rep: "<<rep<<endl;
    if(pinddown==0)rmin=1E-7;
    else rmin=1E-6;
    R_table_end=table_irr[r0ind][pinddown];
    //Rlong=rcurr+5*sqrt(6.0*Dtot*deltat);
    //    Rdown=psample_final_2D( rmin,  r0, double bindrad, double passoc, double Rlong, double R_table_end, double rnum, gsl_matrix *pcume, double RstepSize);
    Rdown=psample_final_2D( rmin,  r0ind*delr0+bindrad, bindrad,  passocvec[r0ind],   Rlong,  R_table_end, rnumdown, pcume, RstepSize);//lookup_edgeR(rmin,  r0ind*delr0+bindrad, r0ind, deltat,  Dtot,  bindrad,  limit.Dmat_ps[1],  alpha2,  cof2, limit,  FitMat,  R_table_end, rnum-pvec);
    R_table_end=table_irr[r0ind+1][pindup];
    Rup=psample_final_2D( rmin, (r0ind+1)*delr0+bindrad, bindrad, passocvec[r0ind+1], Rlong, R_table_end, rnumup, pcume, RstepSize);//lookup_edgeR(rmin,  (r0ind+1)*delr0+bindrad, r0ind+1, deltat,  Dtot,  bindrad,  limit.Dmat_ps[1],  alpha2,  cof2, limit,  FitMat,  R_table_end, rnum-pvec);
    slope=(Rup-Rdown)/delr0;
    R=Rdown+slope*(rcurr-(r0ind*delr0+bindrad));
    
    if(R>Rlong){
      cout <<" Reject R as too long: "<<R<<" Rlong: "<<Rlong<<endl;
      rnum=rand_gsl()*1.0;
      while(rnum<passoc_sample){
	rnum=rand_gsl()*1.0;
      }
      rnumdown=(rnum-passoc_sample)/(1.0-passoc_sample)*(1.0-passocvec[r0ind]);
      rnumup=(rnum-passoc_sample)/(1.0-passoc_sample)*(1.0-passocvec[r0ind+1]);
      delpdown=delpvec[r0ind];//(1.0-passocvec[r0ind])/(1.0*MAXP);
      pinddown=int(rnumdown/delpdown);
      delpup=delpvec[r0ind+1];//(1.0-passocvec[r0ind+1])/(1.0*MAXP);
      pindup=int(rnumup/delpup);
      Rdown=lookup_psample(rnumdown, pinddown, delpdown, table_irr[r0ind]);//using current pos, you can calculate survival prob 
      Rup=lookup_psample(rnumup, pindup, delpup, table_irr[r0ind+1]);//using current pos, you can calculate survival prob 
      slope=(Rup-Rdown)/delr0;
      R=Rdown+slope*(rcurr-(r0ind*delr0+bindrad));
    }
  }else{
    /*choose a new separation given your current separation from the irreversible propagator*/
    Rdown=lookup_psample(rnumdown, pinddown, delpdown, table_irr[r0ind]);//using current pos, you can calculate survival prob 
    Rup=lookup_psample(rnumup, pindup, delpup, table_irr[r0ind+1]);//using current pos, you can calculate survival prob 
    slope=(Rup-Rdown)/delr0;
    R=Rdown+slope*(rcurr-(r0ind*delr0+bindrad));
  }
  return R;
}
