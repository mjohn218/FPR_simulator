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

void move_rot_proteinsPBCCELL(int p1, Fullmol *bases, Complex *ind_com, double **traj, int *movestat, double **trajR, double *M, Parms plist) {
	int i, j;
	int k = bases[p1].mycomplex;
	int s1 = ind_com[k].mysize;
	double t1 = bases[p1].mytime;
	double *v = new double[3];
	double *v2 = new double[3];

	double dx = traj[k][0];
	double dy = traj[k][1];
	double dz = traj[k][2];
	rotationEuler(trajR[k][0], trajR[k][1], trajR[k][2], M);

	trajR[k][0] = 0;
	trajR[k][1] = 0;
	trajR[k][2] = 0;

	traj[k][0] = 0;
	traj[k][1] = 0;
	traj[k][2] = 0;

	double x0 = ind_com[k].xcom;
	double y0 = ind_com[k].ycom;
	double z0 = ind_com[k].zcom;

	ind_com[k].xcom += dx;
	ind_com[k].ycom += dy;
	ind_com[k].zcom += dz;
	ind_com[k].xcom -= plist.xboxl * round(ind_com[k].xcom / plist.xboxl);
	ind_com[k].ycom -= plist.yboxl * round(ind_com[k].ycom / plist.yboxl);
//	ind_com[k].zcom -= plist.zboxl * round(ind_com[k].zcom / plist.zboxl);

//update protein COM
	int mp;
	for (i = 0; i < s1; i++) {
		mp = ind_com[k].plist[i];
		/*We've moved them, don't move them again*/

		movestat[mp] = 2;
		v[0] = bases[mp].xcom - x0;
		v[1] = bases[mp].ycom - y0;
		v[2] = bases[mp].zcom - z0;
		v[0] -= plist.xboxl * round(v[0] / plist.xboxl);
		v[1] -= plist.yboxl * round(v[1] / plist.yboxl);
//		v[2] -= plist.zboxl * round(v[2] / plist.zboxl);

		rotate(v, M, v2);
		/*first would make xcom=x0+v2, then would also add dx */
		bases[mp].xcom = x0 + dx + v2[0];
		bases[mp].ycom = y0 + dy + v2[1];
		bases[mp].zcom = z0 + dz + v2[2];
		bases[mp].xcom -= plist.xboxl * round(bases[mp].xcom / plist.xboxl);
		bases[mp].ycom -= plist.yboxl * round(bases[mp].ycom / plist.yboxl);
//		bases[mp].zcom -= plist.zboxl * round(bases[mp].zcom / plist.zboxl);

//update interface coords
		for (j = 0; j < bases[mp].ninterface; j++) {
			v[0] = bases[mp].x[j] - x0;
			v[1] = bases[mp].y[j] - y0;
			v[2] = bases[mp].z[j] - z0;

			v[0] -= plist.xboxl * round(v[0] / plist.xboxl);
			v[1] -= plist.yboxl * round(v[1] / plist.yboxl);
//			v[2] -= plist.zboxl * round(v[2] / plist.zboxl);

			rotate(v, M, v2);
			/*first would make xcom=x0+v2, then would also add dx */
//			if((v2[0]+v2[1]+v2[2])!=0){
//				cout<<v[0]<<"\t"<<v[1]<<"\t"<<v[2]<<endl;
//				cout<<v2[0]<<"\t"<<v2[1]<<"\t"<<v2[2]<<endl;
//			}
			bases[mp].x[j] = x0 + dx + v2[0];
			bases[mp].y[j] = y0 + dy + v2[1];
			bases[mp].z[j] = z0 + dz + v2[2];
			bases[mp].x[j] -= plist.xboxl * round(bases[mp].x[j] / plist.xboxl);
			bases[mp].y[j] -= plist.yboxl * round(bases[mp].y[j] / plist.yboxl);
//			bases[mp].z[j] -= plist.zboxl * round(bases[mp].z[j] / plist.zboxl);

		}
	}
	delete[] v;
	delete[] v2;
}
