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

void rot_euler_1p(int p1, Fullmol *bases, double deltat, double *M) {
	/*rotate the proteins around the xyz axes*/
	double t1 = sqrt(bases[p1].Drx * 2.0 * deltat) * GaussV();
	double t2 = sqrt(bases[p1].Dry * 2.0 * deltat) * GaussV();
	double t3 = sqrt(bases[p1].Drz * 2.0 * deltat) * GaussV();
	int i;

	rotationEuler(t1, t2, t3, M);
	double *v = new double[3];
	double *v2 = new double[3];
	for (i = 0; i < bases[p1].ninterface; i++) {
		v[0] = bases[p1].x[i] - bases[p1].xcom;
		v[1] = bases[p1].y[i] - bases[p1].ycom;
		v[2] = bases[p1].z[i] - bases[p1].zcom;
		rotate(v, M, v2);
		bases[p1].x[i] = bases[p1].xcom + v2[0];
		bases[p1].y[i] = bases[p1].ycom + v2[1];
		bases[p1].z[i] = bases[p1].zcom + v2[2];
	}
	delete[] v;
	delete[] v2;
}
