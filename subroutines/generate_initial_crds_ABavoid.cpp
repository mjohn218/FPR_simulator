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

void generate_initial_crds_ABavoid(Parms plist, Fullmol *bases, int *Ncopy, Complex *ind_com, double *bindrad, int zeroB) {

	int i, j;
	for (i = 0; i < plist.Ntotalmol; i++) {
		bases[i].xcom = plist.xboxl * rand_gsl() - plist.xboxl / 2.0;
		bases[i].ycom = plist.yboxl * rand_gsl() - plist.yboxl / 2.0;
		if (bases[i].Dz == 0)
			bases[i].zcom = -plist.zboxl / 2.0;
		else
			bases[i].zcom = plist.zboxl * rand_gsl() - plist.zboxl / 2.0;

	}
	if (zeroB == 1) {
		cout << "FORCING the initial particle to the origin !" << endl;
		bases[0].xcom = 0;
		bases[0].ycom = 0;
		bases[0].zcom = 0;
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
		for (i = 0; i < Ncopy[0]; i++) {
			/*THE DIFFERENCE TO THE ONE ABOVE< ONLY A AND B NEED AVOID, As can overlap each other    */
			for (j = Ncopy[0]; j < plist.Ntotalmol; j++) {
				dx = bases[j].xcom - bases[i].xcom;
				dy = bases[j].ycom - bases[i].ycom;
				dz = bases[j].zcom - bases[i].zcom;
				dx -= plist.xboxl * round(dx / plist.xboxl);
				dy -= plist.yboxl * round(dy / plist.yboxl);
				dz -= plist.zboxl * round(dz / plist.zboxl);
				d2 = dx * dx + dy * dy + dz * dz;
				if (d2 < pbindr2) {

					bases[j].xcom = plist.xboxl * rand_gsl() - plist.xboxl / 2.0;
					bases[j].ycom = plist.yboxl * rand_gsl() - plist.yboxl / 2.0;
					if (bases[j].Dz == 0)
						bases[j].zcom = -plist.zboxl / 2.0;
					else
						bases[j].zcom = plist.zboxl * rand_gsl() - plist.zboxl / 2.0;

					noverlap++;
					bflag = 1;
				}

			}
		}
		cout << "Noverlap  " << noverlap << endl;
	}

	/*copy coords into interface pos and complex pos*/
	for (i = 0; i < plist.Ntotalmol; i++) {
		bases[i].x[0] = bases[i].xcom;
		bases[i].y[0] = bases[i].ycom;
		bases[i].z[0] = bases[i].zcom;
		ind_com[i].xcom = bases[i].xcom;
		ind_com[i].ycom = bases[i].ycom;
		ind_com[i].zcom = bases[i].zcom;
		ind_com[i].plist[0] = i;
		ind_com[i].mysize = 1;
		bases[i].mycomplex = i;
	}

}
