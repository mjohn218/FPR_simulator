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

void associate_translate_measure(int p1, int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist) {

	int prod = Rlist[mu][2];
	int iind = ihome[i1];
	int iind2 = ihome[i2];

	//get values of original associating complexes
	int c1 = bases[p1].mycomplex;
	int s1 = ind_com[c1].mysize;
	int c2 = bases[p2].mycomplex;
	int s2 = ind_com[c2].mysize;
	int newsize = s1 + s2;

//   double *v1=new double[3];
//   double *v2=new double[3];
//   double *M=new double[9];
//   double *Mneg=new double[9];//same axis, negative angle
//   double *u=new double[3];
	double *dtrans = new double[3];
	double *drev = new double[3];

	/*Now move protein*/
	double Dxtot = ind_com[c1].Dx + ind_com[c2].Dx;
	double Dytot = ind_com[c1].Dy + ind_com[c2].Dy;

	double Dztot = ind_com[c1].Dz + ind_com[c2].Dz;
	double tol = 1E-16;
	if (Dztot < tol)
		Dztot = 1; //otherwise you divide by zero 

	//distance between the associating proteins interfaces!
	/*The sign needs to be switched to ensure the COM's get rotated
	 to the outside of interfaces final spot*/
	double dx = -bases[p2].x[iind2] + bases[p1].x[iind];
	double dy = -bases[p2].y[iind2] + bases[p1].y[iind];
	double dz = -bases[p2].z[iind2] + bases[p1].z[iind];

	double theta;

	//distance to move to place interfaces, along vector v 
	dtrans[0] = -dx * ind_com[c1].Dx / Dxtot;
	dtrans[1] = -dy * ind_com[c1].Dy / Dytot;
	dtrans[2] = -dz * ind_com[c1].Dz / Dztot;

	drev[0] = +dx * ind_com[c2].Dx / Dxtot;
	drev[1] = +dy * ind_com[c2].Dy / Dytot;
	drev[2] = +dz * ind_com[c2].Dz / Dztot;

	translate_int(p1, c1, ind_com, bases, dtrans);

	translate_int(p2, c2, ind_com, bases, drev);

	/*Determine if you crashed two clathrins together
	 unless they are in the same complex, in which case you are closing a loop.
	 */
	int cancel;
	if (c1 != c2)
		cancel = measure_overlap(c1, c2, ind_com, bases, plist.pclath);
	else
		cancel = 0; //closing a loop, don't measure overlap becasue all proteins are in same complex 

	if (cancel == 0) {
		/*Update the status, free and bound lists of these two proteins*/
		//change status of the interface
		bases[p1].istatus[iind] = prod;
		bases[p2].istatus[iind2] = prod;
		//no longer free
		int i;
		for (i = 0; i < bases[p1].nfree; i++) {
			if (bases[p1].freelist[i] == i1) {
				bases[p1].freelist[i] = bases[p1].freelist[bases[p1].nfree - 1]; //put last reaction in place of this one
				i = bases[p1].nfree;
			}
		}
		bases[p1].nfree -= 1;
		for (i = 0; i < bases[p2].nfree; i++) {
			if (bases[p2].freelist[i] == i2) {
				bases[p2].freelist[i] = bases[p2].freelist[bases[p2].nfree - 1]; //put last reaction in place of this one
				i = bases[p2].nfree;
			}
		}
		bases[p2].nfree -= 1;

		/*add this as a possible dissociation reaction to both proteins
		 so we know who they are bound to.
		 */
		bases[p1].bndlist[bases[p1].nbnd] = prod; //put this reaction as last possible for uni 
		bases[p1].nbnd += 1; //now a dissociation reaction is possible

		bases[p2].bndlist[bases[p2].nbnd] = prod;
		bases[p2].nbnd += 1; //now a dissociation reaction is possible

		bases[p1].partner[iind] = p2;
		bases[p2].partner[iind2] = p1;

		/*add c2's proteins to c1's list*/

		/*Get new COM of complex*/

		double totalmassx = 0;
		double totalmassy = 0;
		double totalmassz = 0;
		ind_com[c1].xcom = 0;
		ind_com[c1].ycom = 0;
		ind_com[c1].zcom = 0;

		/*Copy final complex into the spot of now gone complex c2 */
		int mp;
		int j;
		int tar;
		if (c1 == c2) {
			cout << "CLOSING A LOOP! " << endl;
			plist.nloop++;
			ind_com[c1].mysize = s1;
			for (i = 0; i < s1; i++) {
				mp = ind_com[c1].plist[i];
				totalmassx += bases[mp].massx;
				totalmassy += bases[mp].massy;
				totalmassz += bases[mp].massz;
				ind_com[c1].xcom += bases[mp].xcom * bases[mp].massx;
				ind_com[c1].ycom += bases[mp].ycom * bases[mp].massy;
				ind_com[c1].zcom += bases[mp].zcom * bases[mp].massz;

			}
			ind_com[c1].xcom /= totalmassx;
			ind_com[c1].ycom /= totalmassy;
			ind_com[c1].zcom /= totalmassz;

			update_radius(c1, ind_com, bases);
			update_diffusion(c1, ind_com, bases);
			update_rot_diffusion(c1, ind_com, bases, plist.prerot);
		} else {
			ind_com[c1].mysize = newsize;
			int t = s1;
			for (i = 0; i < s2; i++) {
				mp = ind_com[c2].plist[i];
				ind_com[c1].plist[t] = mp; //add these proteins to the first complex
				t++;
				bases[mp].mycomplex = c1;
			}
			for (i = 0; i < newsize; i++) {
				mp = ind_com[c1].plist[i];
				totalmassx += bases[mp].massx;
				totalmassy += bases[mp].massy;
				totalmassz += bases[mp].massz;
				ind_com[c1].xcom += bases[mp].xcom * bases[mp].massx;
				ind_com[c1].ycom += bases[mp].ycom * bases[mp].massy;
				ind_com[c1].zcom += bases[mp].zcom * bases[mp].massz;

			}
			ind_com[c1].xcom /= totalmassx;
			ind_com[c1].ycom /= totalmassy;
			ind_com[c1].zcom /= totalmassz;

			update_radius(c1, ind_com, bases);
			update_diffusion(c1, ind_com, bases);
			update_rot_diffusion(c1, ind_com, bases, plist.prerot);

			plist.ntotalcomplex -= 1;
			tar = plist.ntotalcomplex;
			if (c2 != plist.ntotalcomplex) {

				/*otherwise, you are just deleting c2 entirely*/
				//copy element by element
				ind_com[c2].xcom = ind_com[tar].xcom;
				ind_com[c2].ycom = ind_com[tar].ycom;
				ind_com[c2].zcom = ind_com[tar].zcom;

				ind_com[c2].Dx = ind_com[tar].Dx;
				ind_com[c2].Dy = ind_com[tar].Dy;
				ind_com[c2].Dz = ind_com[tar].Dz;

				ind_com[c2].radR = ind_com[tar].radR;

				for (j = 0; j < ind_com[tar].mysize; j++) {
					ind_com[c2].plist[j] = ind_com[tar].plist[j];
				}
				ind_com[c2].mysize = ind_com[tar].mysize;

				for (j = 0; j < ind_com[c2].mysize; j++) {
					mp = ind_com[c2].plist[j];
					bases[mp].mycomplex = c2;
				}
			}
		}
	} else {
		/*IN THIS CASE< CLATHRINS CRASHED TOGETHER
		 un-bind them, allow to diffuse apart*/
		cout << "Unbind clathrins, crashed together ! " << endl;

		dtrans[0] *= -1;
		dtrans[1] *= -1;
		dtrans[2] *= -1;

		drev[0] *= -1;
		drev[1] *= -1;
		drev[2] *= -1;

		/*Update the position of proteins in complex one*/
		translate_int(p1, c1, ind_com, bases, dtrans);

		/*Update position of proteins in c2*/
		translate_int(p2, c2, ind_com, bases, drev);

	}

	cout << "final complex com: " << ind_com[c1].xcom << ' ' << ind_com[c1].ycom << ' ' << ind_com[c1].zcom << " radius: " << ind_com[c1].radR << " Dr: " << ind_com[c1].Drx << " Dtrans: " << ind_com[c1].Dx << endl;

//   delete[] u1;
//   delete[] u2;
//   delete[] v1;
//   delete[] v2;
//   delete[] u;
//   delete[] M;
//   delete[] Mneg;
	delete[] dtrans;
	delete[] drev;

}
