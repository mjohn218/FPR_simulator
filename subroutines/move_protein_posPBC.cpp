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

void move_protein_posPBC(int p1, Fullmol *bases, Complex *ind_com, double **traj, int *movestat, Parms &plist) {
	int i, j;
	int k = bases[p1].mycomplex;
	int s1 = ind_com[k].mysize;
	double t1 = bases[p1].mytime;
	double dx = traj[k][0];
	double dy = traj[k][1];
	double dz = traj[k][2];

	traj[k][0] = 0;
	traj[k][1] = 0;
	traj[k][2] = 0;
	//cout <<"protein: "<<p1<<" complex: "<<k<<" size: "<<s1<<endl;
	//cout <<"traj "<<traj[p1][0]<<' '<<traj[p1][1]<<' '<<traj[p1][2]<<endl;

	ind_com[k].xcom += dx;
	ind_com[k].ycom += dy;
	ind_com[k].zcom += dz;
	ind_com[k].xcom -= plist.xboxl * round(ind_com[k].xcom / plist.xboxl);
	ind_com[k].ycom -= plist.yboxl * round(ind_com[k].ycom / plist.yboxl);
	ind_com[k].zcom -= plist.zboxl * round(ind_com[k].zcom / plist.zboxl);

	//update protein COM
	int mp;
	for (i = 0; i < s1; i++) {
		mp = ind_com[k].plist[i];
		/*We've moved them, don't move them again, by setting movestat=2*/
		movestat[mp] = 2;
		bases[mp].xcom += dx;
		bases[mp].ycom += dy;
		bases[mp].zcom += dz;
		bases[mp].xcom -= plist.xboxl * round(bases[mp].xcom / plist.xboxl);
		bases[mp].ycom -= plist.yboxl * round(bases[mp].ycom / plist.yboxl);
		bases[mp].zcom -= plist.zboxl * round(bases[mp].zcom / plist.zboxl);

		//update interface coords
		for (j = 0; j < bases[mp].ninterface; j++) {
			bases[mp].x[j] += dx;
			bases[mp].y[j] += dy;
			bases[mp].z[j] += dz;
			bases[mp].x[j] -= plist.xboxl * round(bases[mp].x[j] / plist.xboxl);
			bases[mp].y[j] -= plist.yboxl * round(bases[mp].y[j] / plist.yboxl);
			bases[mp].z[j] -= plist.zboxl * round(bases[mp].z[j] / plist.zboxl);

		}
	}
}
