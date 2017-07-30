#include <stdio.h>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

void TBLpirr(gsl_matrix*& mpir, int veclen, double bindrad, double Dtot, double kr, double deltat, double Rmax, double RstepSize, double Rlong);
void TBLsur(gsl_matrix*& msur, int veclen, double bindrad, double Dtot, double kr, double deltat, double Rmax, double RstepSize);
void TBLsur_hbin(gsl_matrix*& msur, int veclen, double bindrad, double Dtot, double kr, double deltat, double Rmax, double RstepSize);
void TBLnorm(gsl_matrix*& mnorm, int veclen, double bindrad, double Dtot, double kr, double deltat, double Rmax, double RstepSize);
void DDmatrixcreate(gsl_matrix*& msur, gsl_matrix*& mnorm, gsl_matrix*& mpir, double bindrad, double Dtot, double kr, double deltat, int veclen);
double DDpirr_pfree_ratio_ps(gsl_matrix* mpir, gsl_matrix* msur, gsl_matrix* mnorm, double r, double Dtot, double deltat, double r0, double ps_prev, double rtol, double bindrad, double RstepSize);
double pirr_pfree_ratio_ps(double rcurr, double r0, double tcurr, double Dtot,double bindrad, double alpha, double ps_prev, double rtol);
int sizelookup(double bindrad, double Dtot, double deltat);
double pirr(gsl_matrix* mpir, gsl_matrix* msur, double RstepSize, double r, double r0, double a);
double DDpsur(gsl_matrix* msur, double Dtot, double deltat, double r0, double a, double RstepSize);
double psur(gsl_matrix* msur, double RstepSize, double r0, double a);
double pnorm(gsl_matrix* mnorm, double RstepSize, double r0, double a);
double fsur(double x, void *p); //Gives us association probability not the survival
double fpir(double x, void *p);
double fnorm(double x, void *p);
double finfinity(double x, void *p);
double pirrfunc(double bindrad, double Dtot, double kr, double tval, double r0, double rcurr);
double ktfunc(double bindrad, double Dtot, double kr, double deltat, double r0);
int odefunc(double t, const double y[], double f[], void *params);
int jac(double t, const double y[], double *dfdy, double dfdt[], void *params);
double fkt(double x, void *p) ;


/*cumulative distribution of p(r,t|r0), over r.*/
double fpir_cum(double x, void *p) ;
void TBLpcume(gsl_matrix*& mpir, int r0bins, int r1bins, double bindrad, double Dtot, double kr, double deltat, double Rmax, double RstepSize, double Rlong);
void TBLpcume_allowerr(gsl_matrix*& mpir, int r0bins, int r1bins, double bindrad, double Dtot, double kr, double deltat, double Rmax, double RstepSize, double Rlong, gsl_matrix *pirr_mat);

/*Look up values in a table of r0, r values*/

double get_R_from_table_2D(int MAXP, double rcurr, int r0ind, double **table_irr, double Rlong, double bindrad, double passoc_sample, double rnum, int TABLE_LIM, double delr0, gsl_matrix *pcume, double RstepSize, double *passocvec, double *delpvec);
double get_R_from_table_2D_even(int MAXP, double rcurr, int r0ind, double **table_irr, double Rlong, double bindrad, double passoc, double rnum, int TABLE_LIM, double delr0, gsl_matrix *pcume, double RstepSize, double delp);
double lookup_psample_2D(double rnum, int pind, double delp, double *table);
void psample_table_2D_passoc(double rmin, double r0, double bindrad, double *table, int MAXP, double delp, double passoc, gsl_matrix *pcume, double RstepSize, double Rlongr );
double psample_final_2D(double rmin, double r0, double bindrad, double passoc, double Rlong, double R_table_end, double rnum, gsl_matrix *pcume, double RstepSize);
//void psample_table(double rmin, double r0, double tcurr, double Dtot, double bindrad, double kact, double *table, int MAXP, double delp, gsl_matrix *pcume, double passoc, double RstepSize);
double pirr_simp(gsl_matrix* mpir, double RstepSize, double r, double r0, double a) ;



/*numerical integration*/
double integrate_numer(int r0ind, double bindrad, double rstar, gsl_matrix* pirr, double RstepSize, int N);
double integrate_trapz(int r0ind, double bindrad, double rstar, gsl_matrix* pirr, double RstepSize, int N);
double integrate_simp(int r0ind, double bindrad, double rstar, gsl_matrix* pirr, double RstepSize, int N);
double integrate_simp38(int r0ind, double bindrad, double rstar, gsl_matrix* pirr, double RstepSize, int N);
