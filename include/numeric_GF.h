#include <fstream>
#include <iostream>
#include <cmath>
#include <gsl/gsl_matrix.h>

using namespace std;

void load_fitparms(ifstream &infile, double **Dmat, int &Nthet, int &dihbins);
double fitparms_interp3(double **Dmat, double ct1, double ct2, double dih1, double delthet, double deldih, int Nthet, int dihbins );
void read_FitMat_Vrparms(double **FitMat, ifstream &infile, Vrnumer &limit);
void read_FitMat_VrparmsRD(double **FitMat, ifstream &infile, Vrnumer &limit);
void read_FitMat_VrparmsRD_mat(double **FitMat, ifstream &infile, Vrnumer &limit);
void read_FitMat_VrparmsRD_matfull(double **FitMat, ifstream &infile, Vrnumer &limit, int ncol);


double lookup_pirr(double R1, double r0, double deltat, double rc, double Dtot, double bindrad, double **FitMat, Vrnumer &limit, double passoc );
double lookup_pirr_pfree_ratio(double R1, double prevsep, double deltat, double ps_prev, double rtol, double rc, double Dtot, double bindrad, double **FitMat, Vrnumer &limit, double passoc);
double vr_pfree_cumulative(double rcurr, double r0, double tcurr, double Dtot, double bindrad, double rc);
double vr_pfree_value_norm(double rcurr, double r0, double tcurr, double Dtot, double bindrad, double rc);



/*Need these for sampling from exact GF for trap problems*/
double get_R_from_table_2D(int MAXP, double rcurr, int r0ind, double **table_irr, double Rlong, double bindrad, double passoc_sample, double rnum, int TABLE_LIM, double delr0, gsl_matrix *pcume, double RstepSize, double *passocvec, double *delpvec);//2D
double get_R_from_table_2D_even(int MAXP, double rcurr, int r0ind, double **table_irr, double Rlong, double bindrad, double passoc, double rnum, int TABLE_LIM, double delr0, gsl_matrix *pcume, double RstepSize, double delp);//2D
double psample_final_2D(double rmin, double r0, double bindrad, double passoc, double Rlong, double R_table_end, double rnum, gsl_matrix *pcume, double RstepSize);
double pirr_simp(gsl_matrix* mpir, double RstepSize, double r, double r0, double a);
void psample_table_2D_passoc(double rmin, double r0, double bindrad, double *table, int MAXP, double delp, double passoc, gsl_matrix *pcume, double RstepSize, double Rlong );


double psample_final_Sfit_2fit(double rmin, double r0, double tcurr,  double bindrad, double passoc, Vrnumer &limit, double *FitMat, double Rlong, double R_table_end, double rnum);
double psample_final_hybrid(double rmin, double r0, double tcurr,  double bindrad, double passoc, Vrnumer &limit, double *FitMat, double Rlong, double R_table_end, double rnum);
double psample_final_pfree_2Sfit(double rmin, double r0, double tcurr,  double bindrad, double passoc, Vrnumer &limit, double *FitMat, double Rlong, double R_table_end, double rnum);
double psample_final_pfree_vr(double rmin, double r0, double tcurr, double Dtot, double bindrad, double rc, double passoc,double R_table_end, double rnum);

double lookup_psample(double rnum, int pind, double delp, double *table);
void psample_table(double rmin, double r0, double tcurr, double Dtot, double bindrad, double alpha, double kact, double *table, int MAXP, double delp, double Rlong);
double get_R_from_table(int MAXP, double rcurr, int r0ind, double **table_irr, double Rlong, double bindrad, double passoc_sample, double rnum, int TABLE_LIM, double delr0, double *passocvec, double *delpvec, double alpha, double kact, double Dtot, double deltat);//3D
double psample_final(double rmin, double r0, double bindrad, double passoc, double Rlong, double R_table_end, double rnum, double alpha, double kact, double Dtot, double deltat);
double lookup_edgeR(double rmin, double r0, int r0ind, double deltat, double Dtot, double bindrad,double Dfit_ps, double alpha_ps,double cof_ps, Vrnumer &limit, double **FitMat, double R_table_end, double rnum );

