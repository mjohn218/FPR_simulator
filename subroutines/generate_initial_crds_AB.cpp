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

void generate_initial_crds_AB(Parms plist, Fullmol *bases, int *Ncopy, Complex *ind_com, double *bindrad) {

	int i, j, jj;
	for (i = 0; i < plist.Ntotalmol; i++) {
//		cout<<rand_gsl()<<endl;
		bases[i].xcom = plist.xboxl * rand_gsl() - plist.xboxl / 2.0;
		bases[i].ycom = plist.yboxl * rand_gsl() - plist.yboxl / 2.0;
		if (bases[i].Dz == 0) {
			bases[i].zcom = -plist.zboxl / 2.0;
		} else {
			bases[i].zcom = plist.zboxl * rand_gsl() - plist.zboxl / 2.0;
		}
	}
	int bflag = 1;
	int bit = 0;
	int noverlap = 0;

	/*This arbitrarily picks a binding radius*/
	double pbindr2 = bindrad[0] * bindrad[0];
	double stretch;
	int Maxit = 50;
	double d2;
	double dx, dy, dz;
	while (bflag == 1 && bit < Maxit) {
		bit++;
		bflag = 0;
		noverlap = 0;
		for (i = 0; i < plist.Ntotalmol; i++) {
			for (j = i + 1; j < plist.Ntotalmol; j++) { //for(j=Ncopy[0];j<plist.Ntotalmol;j++)
				dx = bases[j].xcom - bases[i].xcom;
				dy = bases[j].ycom - bases[i].ycom;
				dz = bases[j].zcom - bases[i].zcom;
				dx -= plist.xboxl * round(dx / plist.xboxl);
				dy -= plist.yboxl * round(dy / plist.yboxl);
//				dz -= plist.zboxl * round(dz / plist.zboxl);
				d2 = dx * dx + dy * dy + dz * dz;
				if (d2 < pbindr2) {

					bases[j].xcom = plist.xboxl * rand_gsl() - plist.xboxl / 2.0;
					bases[j].ycom = plist.yboxl * rand_gsl() - plist.yboxl / 2.0;
					if (bases[j].Dz == 0) {
						bases[j].zcom = -plist.zboxl / 2.0;
					} else {
						bases[j].zcom = plist.zboxl * rand_gsl() - plist.zboxl / 2.0;
					}
					noverlap++;
					bflag = 1;
				}
			}
		}
		cout << "Noverlap  " << noverlap << endl;
	}

	/*copy coords into interface pos and complex pos*/
	for (i = 0; i < plist.Ntotalmol; i++) {
		for (jj = 0; jj < bases[i].ninterface; jj++) {
			bases[i].x[jj] = bases[i].xcom;
			bases[i].y[jj] = bases[i].ycom;
			if (bases[i].Dz == 0) {
				bases[i].z[jj] = -plist.zboxl / 2.0;
			} else {
				bases[i].z[jj] = bases[i].zcom;
			}
		}
		ind_com[i].xcom = bases[i].xcom;
		ind_com[i].ycom = bases[i].ycom;
		if (bases[i].Dz == 0) {
			ind_com[i].zcom = -plist.zboxl / 2.0;
		} else {
			ind_com[i].zcom = bases[i].zcom;
		}
		ind_com[i].plist[0] = i;
		ind_com[i].mysize = 1;
		bases[i].mycomplex = i;
	}

}
