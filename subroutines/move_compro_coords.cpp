#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "reactions.h"
#include "vol_help.h"
#include "rand_gsl.h"
#include "Faddeeva.hh"
#include "utility_calls.h"
#include "vector_rot_calls.h"

void move_compro_coords(int c1, Complex *ind_com, Fullmol *bases, double *chg1) {
	/*Update the position of proteins in complex c1 by the displacement vector chg1
	 and correspondingly update the COM of the complex.
	 */
	double totalmassx = 0;
	double totalmassy = 0;
	double totalmassz = 0;
	int j, mp, k;

	ind_com[c1].xcom = 0;
	ind_com[c1].ycom = 0;
	ind_com[c1].zcom = 0;

	for (j = 0; j < ind_com[c1].mysize; j++) {
		mp = ind_com[c1].plist[j];
		//    bases[mp].mycomplexsize=ind_com[c1].mysize;
		totalmassx += bases[mp].massx;
		totalmassy += bases[mp].massy;
		totalmassz += bases[mp].massz;

		for (k = 0; k < bases[mp].ninterface; k++) {
			bases[mp].x[k] += chg1[0];
			bases[mp].y[k] += chg1[1];
			bases[mp].z[k] += chg1[2];
		}
		bases[mp].xcom += chg1[0];
		bases[mp].ycom += chg1[1];
		bases[mp].zcom += chg1[2];
		ind_com[c1].xcom += bases[mp].xcom * bases[mp].massx;
		ind_com[c1].ycom += bases[mp].ycom * bases[mp].massy;
		ind_com[c1].zcom += bases[mp].zcom * bases[mp].massz;
	}

	ind_com[c1].xcom /= totalmassx;
	ind_com[c1].ycom /= totalmassy;
	ind_com[c1].zcom /= totalmassz;

}
