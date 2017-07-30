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

void reflect_complex_rad(int p1, Fullmol *bases, Complex *ind_com, double xboxl, double yboxl, double zboxl) {
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
	int flagz = 0;
	xchg = ind_com[k].xcom - xboxl / 2.0;
	if ((xchg + rad) > 0) {
		xtot = -(xchg + rad); //shift x coordinates back
		flag = 1;
	} else if ((xchg - rad) < -xboxl) {
		xtot = -(ind_com[k].xcom - rad + xboxl / 2.0);
		flag = 1;
	}
	ychg = ind_com[k].ycom - yboxl / 2.0;
	if ((ychg + rad) > 0) {
		ytot = -(ychg + rad); //shift x coordinates back
		flag = 1;
	} else if ((ychg - rad) < -yboxl) {
		ytot = -(ind_com[k].ycom - rad + yboxl / 2.0);
		flag = 1;
	}
	zchg = ind_com[k].zcom - zboxl / 2.0;
	if ((zchg + rad) > 0) {

		flagz = 1;
	} else if ((zchg - rad) < -zboxl) {

		flagz = 1;
	}

	int mp;
	if (flag == 1) {

		/*Put back inside the box*/
		ind_com[k].xcom += 2.0 * xtot;
		ind_com[k].ycom += 2.0 * ytot;

		//update protein COM

		for (i = 0; i < s1; i++) {
			mp = ind_com[k].plist[i];
			bases[mp].xcom += 2.0 * xtot;
			bases[mp].ycom += 2.0 * ytot;

			//update interface coords
			for (j = 0; j < bases[mp].ninterface; j++) {
				bases[mp].x[j] += 2.0 * xtot;
				bases[mp].y[j] += 2.0 * ytot;

			}
		}

	}

	/*Z is separate to allow the interfaces to approach to the membrane
	 but don't need to test if the entire complex is far enough
	 away from the boundary.
	 */
	if (flagz == 1) {

		flag = 0;
		ztot = 0;
		for (i = 0; i < s1; i++) {
			mp = ind_com[k].plist[i];

			/*measure each interface to z plane*/
			for (j = 0; j < bases[mp].ninterface; j++) {
				zchg = bases[mp].z[j] - zboxl / 2.0;

				if (zchg > 0) {
					flag = 1;
					if (-zchg < ztot)
						ztot = -zchg;

				} else if (zchg < -zboxl) {
					flag = 1;
					if (-(zchg + zboxl) > ztot)
						ztot = -(zchg + zboxl);

				}
			}
		}
		if (flag == 1) {

			/*Put back inside the box*/
			ind_com[k].zcom += 2.0 * ztot;

			//update protein COM

			for (i = 0; i < s1; i++) {
				mp = ind_com[k].plist[i];
				bases[mp].zcom += 2.0 * ztot;
				//update interface coords
				for (j = 0; j < bases[mp].ninterface; j++) {
					bases[mp].z[j] += 2.0 * ztot;
				}
			}

		}
	}

}
