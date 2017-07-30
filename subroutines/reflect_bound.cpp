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

void reflect_bound(int p1, Fullmol *bases, Complex *ind_com, double xboxl, double yboxl, double zboxl) {
	/*correct if protein complex COM position went outside box with reflection*/
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
	xchg = ind_com[k].xcom - xboxl / 2.0;
	if (xchg > 0) {
		xtot = -xchg; //shift x coordinates back
		flag = 1;
	} else if (xchg < -xboxl) {
		xtot = -(ind_com[k].xcom + xboxl / 2.0);
		flag = 1;
	}
	ychg = ind_com[k].ycom - yboxl / 2.0;
	if (ychg > 0) {
		ytot = -ychg; //shift x coordinates back
		flag = 1;
	} else if (ychg < -yboxl) {
		ytot = -(ind_com[k].ycom + yboxl / 2.0);
		flag = 1;
	}
	zchg = ind_com[k].zcom - zboxl / 2.0;
	if (zchg > 0) {
		ztot = -zchg; //shift x coordinates back
		flag = 1;
	} else if (zchg < -zboxl) {
		ztot = -(ind_com[k].zcom + zboxl / 2.0);
		flag = 1;
	}
	int mp;
	if (flag == 1) {
		//    cout <<"moved outside of the box after association/dissociation ! "<<p1<<endl;
//     write_crds(bases, p1);

		/*Put back inside the box*/
		ind_com[k].xcom += 2.0 * xtot;
		ind_com[k].ycom += 2.0 * ytot;
		ind_com[k].zcom += 2.0 * ztot;
		//update protein COM

		for (i = 0; i < s1; i++) {
			mp = ind_com[k].plist[i];
			bases[mp].xcom += 2.0 * xtot;
			bases[mp].ycom += 2.0 * ytot;
			bases[mp].zcom += 2.0 * ztot;

			//update interface coords
			for (j = 0; j < bases[mp].ninterface; j++) {
				bases[mp].x[j] += 2.0 * xtot;
				bases[mp].y[j] += 2.0 * ytot;
				bases[mp].z[j] += 2.0 * ztot;
			}
		}

	}
}
