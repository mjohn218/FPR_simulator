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

void reflect_complex_rad_rot(int p1, Fullmol *bases, Complex *ind_com, double xboxl, double yboxl, double zboxl) {
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
	int flagy = 0;
	int flagx = 0;
	double xtot0, ytot0, ztot0;
	double tol = 1E-11;
	xchg = ind_com[k].xcom - xboxl / 2.0;
	if ((xchg + rad) > 0) {
		xtot0 = -(xchg + rad); //shift x coordinates back
		flagx = 1;
	} else if ((xchg - rad) < -xboxl) {
		xtot0 = -(ind_com[k].xcom - rad + xboxl / 2.0);
		flagx = 1;
	}
	ychg = ind_com[k].ycom - yboxl / 2.0;
	if ((ychg + rad) > 0) {
		ytot0 = -(ychg + rad); //shift x coordinates back
		flagy = 1;
	} else if ((ychg - rad) < -yboxl) {
		ytot0 = -(ind_com[k].ycom - rad + yboxl / 2.0);
		flagy = 1;
	}
	zchg = ind_com[k].zcom - zboxl / 2.0;
	if ((zchg + rad) > 0) {
		ztot0 = -(zchg + rad); //shift x coordinates back
		flagz = 1;
	} else if ((zchg - rad) < -zboxl) {
		ztot0 = -(ind_com[k].zcom - rad + zboxl / 2.0);
		flagz = 1;
	}

	int mp;

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
		//if(abs(ztot-ztot0)>tol)cout <<"Different calculations of z beyond box! "<<ztot0<<" final: "<<ztot<<endl;
	}
	if (flagy == 1) {

		flag = 0;
		ytot = 0;
		for (i = 0; i < s1; i++) {
			mp = ind_com[k].plist[i];

			/*measure each interface to y plane*/
			for (j = 0; j < bases[mp].ninterface; j++) {
				ychg = bases[mp].y[j] - yboxl / 2.0;

				if (ychg > 0) {
					flag = 1;
					if (-ychg < ytot)
						ytot = -ychg;

				} else if (ychg < -yboxl) {
					flag = 1;
					if (-(ychg + yboxl) > ytot)
						ytot = -(ychg + yboxl);

				}
			}
		}
		if (flag == 1) {

			/*Put back inside the box*/
			ind_com[k].ycom += 2.0 * ytot;

			//update protein COM

			for (i = 0; i < s1; i++) {
				mp = ind_com[k].plist[i];
				bases[mp].ycom += 2.0 * ytot;
				//update interface coords
				for (j = 0; j < bases[mp].ninterface; j++) {
					bases[mp].y[j] += 2.0 * ytot;
				}
			}

		}
		//if(abs(ytot-ytot0)>tol)cout <<"Different calculations of y beyond box! "<<ytot0<<" final: "<<ytot<<endl;
	}
	if (flagx == 1) {

		flag = 0;
		xtot = 0;
		for (i = 0; i < s1; i++) {
			mp = ind_com[k].plist[i];

			/*measure each interface to x plane*/
			for (j = 0; j < bases[mp].ninterface; j++) {
				xchg = bases[mp].x[j] - xboxl / 2.0;

				if (xchg > 0) {
					flag = 1;
					if (-xchg < xtot)
						xtot = -xchg;

				} else if (xchg < -xboxl) {
					flag = 1;
					if (-(xchg + xboxl) > xtot)
						xtot = -(xchg + xboxl);

				}
			}
		}
		if (flag == 1) {

			/*Put back inside the box*/
			ind_com[k].xcom += 2.0 * xtot;

			//update protein COM

			for (i = 0; i < s1; i++) {
				mp = ind_com[k].plist[i];
				bases[mp].xcom += 2.0 * xtot;
				//update interface coords
				for (j = 0; j < bases[mp].ninterface; j++) {
					bases[mp].x[j] += 2.0 * xtot;
				}
			}

		}
		//if(abs(xtot-xtot0)>tol)cout <<"Different calculations of x bexond box! "<<xtot0<<" final: "<<xtot<<endl;
	}

}
