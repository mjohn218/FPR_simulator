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

void reflect_traj_complex_rad(int p1, Fullmol *bases, Complex *ind_com, double xboxl, double yboxl, double zboxl, double **traj) {
	int i, j;
	int k = bases[p1].mycomplex;
	int s1 = ind_com[k].mysize;
	double xchg;
	double ychg;
	double zchg;
	double xtot = 0.0;
	double ytot = 0.0;
	double ztot = 0.0;
	int flag = 0;
	double rad = ind_com[k].radR;

	double currx = ind_com[k].xcom + traj[k][0];
	double curry = ind_com[k].ycom + traj[k][1];
	double currz = ind_com[k].zcom + traj[k][2];

	xchg = currx - xboxl / 2.0;

	if ((xchg + rad) > 0) {
		xtot = -(xchg + rad); //shift x coordinates back
		flag = 1;
	} else if ((xchg - rad) < -xboxl) {
		xtot = -(currx - rad + xboxl / 2.0);
		flag = 1;
	}
	ychg = curry - yboxl / 2.0;
	if ((ychg + rad) > 0) {
		ytot = -(ychg + rad); //shift x coordinates back
		flag = 1;
	} else if ((ychg - rad) < -yboxl) {
		ytot = -(curry - rad + yboxl / 2.0);
		flag = 1;
	}
	zchg = currz - zboxl / 2.0;
	if ((zchg + rad) > 0) {
		ztot = -(zchg + rad); //shift x coordinates back
		flag = 1;
	} else if ((zchg - rad) < -zboxl) {
		ztot = -(currz - rad + zboxl / 2.0);
		flag = 1;
	}

	int mp;
	if (flag == 1) {

		/*Put back inside the box*/
		traj[k][0] += 2.0 * xtot;
		traj[k][1] += 2.0 * ytot;
		traj[k][2] += 2.0 * ztot;

	}

}
