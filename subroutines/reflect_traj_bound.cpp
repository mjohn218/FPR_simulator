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

void reflect_traj_bound(int p1, Fullmol *bases, Complex *ind_com, double xboxl, double yboxl, double zboxl, double **traj) {
	/*correct if current displacement to particle (in traj) will move outside the box with reflection*/
	int k = bases[p1].mycomplex;
	double xchg;
	double ychg;
	double zchg;
	double xtot = 0.0;
	double ytot = 0.0;
	double ztot = 0.0;
	int flag = 0;
	double currx = bases[p1].xcom + traj[k][0];
	double curry = bases[p1].ycom + traj[k][1];
	double currz = bases[p1].zcom + traj[k][2];

	xchg = currx - xboxl / 2.0;
	if (xchg > 0) {
		xtot = -xchg; //shift x coordinates back
		flag = 1;
	} else if (xchg < -xboxl) {
		xtot = -(currx + xboxl / 2.0);
		flag = 1;
	}
	ychg = curry - yboxl / 2.0;
	if (ychg > 0) {
		ytot = -ychg; //shift x coordinates back
		flag = 1;
	} else if (ychg < -yboxl) {
		ytot = -(curry + yboxl / 2.0);
		flag = 1;
	}
	zchg = currz - zboxl / 2.0;
	if (zchg > 0) {
		ztot = -zchg; //shift x coordinates back
		flag = 1;
	} else if (zchg < -zboxl) {
		ztot = -(currz + zboxl / 2.0);
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
