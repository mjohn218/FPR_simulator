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

void get_bin2(Parms plist, Fullmol *bases, double cellx, double celly, double cellz, int Nx, int Ny, int Nz, int *binlist, int *npb, int MAXPERBIN, int Ncell, Complex *ind_com) {
	int i, j;
	int mybin;
	int xbin, ybin, zbin;
	double small = 1E-6;
	int k1, s1, f;
	double chgx, chgy, chgz;
	double pdx, pdy, pdz;
	int it = 0;
	int mp;
	for (j = 0; j < Ncell; j++)
		npb[j] = 0;

	for (i = 0; i < plist.Ntotalmol; i++) {
		xbin = int((bases[i].xcom + plist.xboxl / 2.0) / cellx);
		ybin = int((bases[i].ycom + plist.yboxl / 2.0) / celly);
		zbin = int(-(bases[i].zcom + small - plist.zboxl / 2.0) / cellz);
		mybin = xbin + ybin * Nx + zbin * Nx * Ny;

		if (mybin > (Ncell - 1)) {

			cout << "OUTSIDE ARRAY: " << i << " " << bases[i].xcom << ' ' << bases[i].ycom << ' ' << bases[i].zcom << endl;
			k1 = bases[i].mycomplex;
			cout << "Complex: " << k1 << ' ' << ind_com[k1].xcom << ' ' << ind_com[k1].ycom << ' ' << ind_com[k1].zcom << endl;
			cout << "Complex size and radius: " << ind_com[k1].mysize << ' ' << ind_com[k1].radR << endl;

			//So put back inside
			chgx = 0;
			chgy = 0;
			chgz = 0;
			mybin = 0;
			pdx = bases[i].xcom - plist.xboxl / 2.0;
			if (pdx > 0) {
				chgx = -pdx - small;
			} else if (pdx < -plist.xboxl) {
				chgx = -(pdx + plist.xboxl) + small;
			}
			pdy = bases[i].ycom - plist.yboxl / 2.0;
			if (pdy > 0) {
				chgy = -pdy - small;
			} else if (pdy < -plist.yboxl) {
				chgy = -(pdy + plist.yboxl) + small;
			}
			pdz = bases[i].zcom - plist.zboxl / 2.0;
			if (pdz > 0) {
				chgz = -pdz - small;
			} else if (pdz < -plist.zboxl) {
				chgz = -(pdz + plist.zboxl) + small;
			}
			s1 = ind_com[k1].mysize;
			for (j = 0; j < s1; j++) {
				mp = ind_com[k1].plist[j];

				bases[mp].xcom += chgx;
				bases[mp].ycom += chgy;
				bases[mp].zcom += chgz;

				//update interface coords
				for (f = 0; f < bases[mp].ninterface; f++) {
					bases[mp].x[f] += chgx;
					bases[mp].y[f] += chgy;
					bases[mp].z[f] += chgz;
				}
			}
			cout << "changed positions by: : " << chgx << ' ' << chgy << ' ' << chgz << endl;

			i = 0; //reassign everyone's box
			for (j = 0; j < Ncell; j++)
				npb[j] = 0;

			it++;
			if (it > 1000) {
				cout << "can't fit in box !" << endl;
				exit(1);
			}
		} else if (mybin < 0) {
			cout << "NEGATIVE OUTSIDE ARRAY: " << i << " " << bases[i].xcom << ' ' << bases[i].ycom << ' ' << bases[i].zcom << endl;
			k1 = bases[i].mycomplex;
			cout << "Complex: " << k1 << ' ' << ind_com[k1].xcom << ' ' << ind_com[k1].ycom << ' ' << ind_com[k1].zcom << endl;
			cout << "Complex size and radius: " << ind_com[k1].mysize << ' ' << ind_com[k1].radR << endl;

			//So put back inside
			chgx = 0;
			chgy = 0;
			chgz = 0;
			mybin = 0;
			pdx = bases[i].xcom - plist.xboxl / 2.0;
			if (pdx > 0) {
				chgx = -pdx - small;
			} else if (pdx < -plist.xboxl) {
				chgx = -(pdx + plist.xboxl) + small;
			}
			pdy = bases[i].ycom - plist.yboxl / 2.0;
			if (pdy > 0) {
				chgy = -pdy - small;
			} else if (pdy < -plist.yboxl) {
				chgy = -(pdy + plist.yboxl) + small;
			}
			pdz = bases[i].zcom - plist.zboxl / 2.0;
			if (pdz > 0) {
				chgz = -pdz - small;
			} else if (pdz < -plist.zboxl) {
				chgz = -(pdz + plist.zboxl) + small;
			}
			s1 = ind_com[k1].mysize;
			for (j = 0; j < s1; j++) {
				mp = ind_com[k1].plist[j];

				bases[mp].xcom += chgx;
				bases[mp].ycom += chgy;
				bases[mp].zcom += chgz;

				//update interface coords
				for (f = 0; f < bases[mp].ninterface; f++) {
					bases[mp].x[f] += chgx;
					bases[mp].y[f] += chgy;
					bases[mp].z[f] += chgz;
				}
			}
			cout << "changed positions by: : " << chgx << ' ' << chgy << ' ' << chgz << endl;
			i = 0;
			for (j = 0; j < Ncell; j++)
				npb[j] = 0;

			it++;
			if (it > 1000) {
				cout << "can't fit in box !" << endl;
				exit(1);
			}
		}
		//    cout <<" Protein : "<<i<<" mybin: "<<mybin<<endl;
		bases[i].mybin = mybin;
		bases[i].mybinind = npb[mybin];
		binlist[mybin * MAXPERBIN + npb[mybin]] = i;
		npb[mybin]++;

	}

}
