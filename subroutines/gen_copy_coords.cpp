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

void gen_copy_coords(int Nprotypes, Protein *wholep, Fullmol* bases, Complex *ind_com, int *Ncopy) {
	int i, j, k;
	int t = 0;
	int p;
	int nint;
	double x, y, z;
	string name2;
	string name1;
	//at start, each complex is an unbound protein 
	for (p = 0; p < Nprotypes; p++) {
		nint = wholep[p].ninterface;
		cout << "protein: " << p << " ninterface: " << nint << endl;
		for (i = 0; i < Ncopy[p]; i++) {

			ind_com[t].xcom = bases[t].xcom;
			ind_com[t].ycom = bases[t].ycom;
			ind_com[t].zcom = bases[t].zcom;

			ind_com[t].mysize = 1;
			ind_com[t].plist[0] = t;
			ind_com[t].Dx = bases[t].Dx;
			ind_com[t].Dy = bases[t].Dy;
			ind_com[t].Dz = bases[t].Dz;
			ind_com[t].Drx = bases[t].Drx;
			ind_com[t].Dry = bases[t].Dry;
			ind_com[t].Drz = bases[t].Drz;

			bases[t].radR = wholep[p].radx;
			ind_com[t].radR = wholep[p].radx;

			for (k = 0; k < nint; k++) {

				bases[t].x[k] = bases[t].xcom;
				bases[t].y[k] = bases[t].ycom;
				bases[t].z[k] = bases[t].zcom;
			}
			//cout <<bases[t].x[0]<<' '<<bases[t].y[0]<<' '<<bases[t].z[0]<<endl;
			t++;
			//cout <<"t: "<<t<<" bases[t].nbnd: "<<bases[t].nbnd<<endl;
		}

	}

}
