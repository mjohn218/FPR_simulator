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

void rotate_onlyPBC(int p1, int c1, Complex *ind_com, Fullmol *bases, double *M, Parms plist) {

	/*First rotate all points in complex 1 around the com of protein p1, and translate them by vector d*/
	double *v = new double[3];
	double *v2 = new double[3];

	int mp;
	int i, k;
	int s1 = ind_com[c1].mysize;
	double pivx = bases[p1].xcom;
	double pivy = bases[p1].ycom;
	double pivz = bases[p1].zcom;

	for (i = 0; i < s1; i++) {
		mp = ind_com[c1].plist[i];
		for (k = 0; k < bases[mp].ninterface; k++) {
			v[0] = bases[mp].x[k] - pivx;
			v[1] = bases[mp].y[k] - pivy;
			v[2] = bases[mp].z[k] - pivz;
			v[0] -= plist.xboxl * round(v[0] / plist.xboxl);
			v[1] -= plist.yboxl * round(v[1] / plist.yboxl);
			v[2] -= plist.zboxl * round(v[2] / plist.zboxl);

			rotate(v, M, v2); //includes the interface that will align

			bases[mp].x[k] = pivx + v2[0];
			bases[mp].y[k] = pivy + v2[1];
			bases[mp].z[k] = pivz + v2[2];
			bases[mp].x[k] -= plist.xboxl * round(bases[mp].x[k] / plist.xboxl);
			bases[mp].y[k] -= plist.yboxl * round(bases[mp].y[k] / plist.yboxl);
			bases[mp].z[k] -= plist.zboxl * round(bases[mp].z[k] / plist.zboxl);

		}
		//rotate COM
		v[0] = bases[mp].xcom - pivx;
		v[1] = bases[mp].ycom - pivy;
		v[2] = bases[mp].zcom - pivz;
		v[0] -= plist.xboxl * round(v[0] / plist.xboxl);
		v[1] -= plist.yboxl * round(v[1] / plist.yboxl);
		v[2] -= plist.zboxl * round(v[2] / plist.zboxl);

		rotate(v, M, v2); //includes the interface that will align

		bases[mp].xcom = pivx + v2[0];
		bases[mp].ycom = pivy + v2[1];
		bases[mp].zcom = pivz + v2[2];
		bases[mp].xcom -= plist.xboxl * round(bases[mp].xcom / plist.xboxl);
		bases[mp].ycom -= plist.yboxl * round(bases[mp].ycom / plist.yboxl);
		bases[mp].zcom -= plist.zboxl * round(bases[mp].zcom / plist.zboxl);

	}
	delete[] v;
	delete[] v2;

}
