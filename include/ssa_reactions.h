
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>



using namespace std;


void mem_to_sol(int c1, Fullmol *bases, Complex *ind_com, double *X, int *ihome, int **myrxn, int **Rlist, int *p_home, int *sol_prod);
void sol_to_mem(int c1, Fullmol *bases, Complex *ind_com, double *X, int *ihome, int **myrxn, int **Rlist, int *p_home, int *mem_prod);
void break_complex_nocrds(int p1, int mu, int kind, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, int p2, int i1, int i2, int *p_home, int **myrxn);
void associate_nocrds(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist);
void set_status_conc(ifstream &startfile, Protein *wholep, Fullmol *bases, Parms &plist, int *Ncopy, double *indivconc);
void sum_state(int Nprotypes, Protein *wholep, int *Nmyrxn, int **myrxn, int **Rlist, int **Del, double *X, int *rxntype);

