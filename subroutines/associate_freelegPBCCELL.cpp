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

void associate_freelegPBCCELL(int p1, int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist) {

	int prod = Rlist[mu][2];
	int iind = ihome[i1];
	int iind2 = ihome[i2];

	//get values of original associating complexes
	int c1 = bases[p1].mycomplex;
	int s1 = ind_com[c1].mysize;
	int c2 = bases[p2].mycomplex;
	int s2 = ind_com[c2].mysize;
	int newsize = s1 + s2;

	double *v1 = new double[3];
	double *v2 = new double[3];
	double *M = new double[9];
	double *Mneg = new double[9]; //same axis, negative angle
	double *u = new double[3];
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

	dx -= plist.xboxl * round(dx / plist.xboxl);
	dy -= plist.yboxl * round(dy / plist.yboxl);
//	dz -= plist.zboxl * round(dz / plist.zboxl);

	double theta;

	//distance to move to place interfaces, along vector v
	dtrans[0] = -dx * ind_com[c1].Dx / Dxtot;
	dtrans[1] = -dy * ind_com[c1].Dy / Dytot;
	dtrans[2] = -dz * ind_com[c1].Dz / Dztot;

	drev[0] = +dx * ind_com[c2].Dx / Dxtot;
	drev[1] = +dy * ind_com[c2].Dy / Dytot;
	drev[2] = +dz * ind_com[c2].Dz / Dztot;

	v1[0] = bases[p1].xcom - bases[p1].x[iind];
	v1[1] = bases[p1].ycom - bases[p1].y[iind];
	v1[2] = bases[p1].zcom - bases[p1].z[iind];
	v1[0] -= plist.xboxl * round(v1[0] / plist.xboxl);
	v1[1] -= plist.yboxl * round(v1[1] / plist.yboxl);
//	v1[2] -= plist.zboxl * round(v1[2] / plist.zboxl);

	v2[0] = bases[p2].xcom - bases[p2].x[iind2];
	v2[1] = bases[p2].ycom - bases[p2].y[iind2];
	v2[2] = bases[p2].zcom - bases[p2].z[iind2];
	v2[0] -= plist.xboxl * round(v2[0] / plist.xboxl);
	v2[1] -= plist.yboxl * round(v2[1] / plist.yboxl);
//	v2[2] -= plist.zboxl * round(v2[2] / plist.zboxl);

	/*Calculate rotation matrix*/

	//dotproduct(v2,v1, theta);
	double dp = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
	double l1 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
	double l2 = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2];
	double dpnorm = dp / sqrt(l1 * l2);
	if (dpnorm > 1) {
		cout << "dp out of range! " << dpnorm << endl;
		theta = 0;
	} else if (dpnorm < -1) {
		cout << "dp out of range! " << dpnorm << endl;
		theta = M_PI;
	} else {
		theta = acos(dpnorm);
	}

	//now calculate the axis of rotation
	crossproduct(v1, v2, u); //u is the UNIT vector of the rotation axis
	double trot = 0.5 * (M_PI - theta); //to push them 180 apart
	/*for rotating v1 around u the negative trot */

	//calc_Rmatrix(u, -trot,M);
	double ux = u[0];
	double uy = u[1];
	double uz = u[2];
	double cthet = cos(+trot);
	double sthet = sin(+trot);
	cout << "angle between: " << trot << " normal: " << ux << ' ' << uy << ' ' << uz << endl;
	M[0] = cthet + ux * ux * (1 - cthet);
	M[1] = ux * uy * (1 - cthet) - uz * sthet; //row 0, column 1, go across first!
	M[2] = ux * uz * (1 - cthet) + uy * sthet;
	M[3] = uy * ux * (1 - cthet) + uz * sthet;
	M[4] = cthet + uy * uy * (1 - cthet);
	M[5] = uy * uz * (1 - cthet) - ux * sthet;
	M[6] = uz * ux * (1 - cthet) - uy * sthet;
	M[7] = uz * uy * (1 - cthet) + ux * sthet;
	M[8] = cthet + uz * uz * (1 - cthet);

	Mneg[0] = M[0];
	Mneg[1] = ux * uy * (1 - cthet) + uz * sthet;
	Mneg[2] = ux * uz * (1 - cthet) - uy * sthet;
	Mneg[3] = uy * ux * (1 - cthet) - uz * sthet;
	Mneg[4] = M[4];
	Mneg[5] = uy * uz * (1 - cthet) + ux * sthet;
	Mneg[6] = uz * ux * (1 - cthet) + uy * sthet;
	Mneg[7] = uz * uy * (1 - cthet) - ux * sthet;
	Mneg[8] = M[8];

	/*Update the position of proteins in complex one*/
	rotate_and_translate_intPBCCELL(p1, c1, ind_com, bases, Mneg, dtrans, iind, plist);

	/*Update position of proteins in c2*/
	rotate_and_translate_intPBCCELL(p2, c2, ind_com, bases, M, drev, iind2, plist);

	/*Rotate to proper orientation*/
	int ind3;
	double *u1 = new double[3];
	double *u2 = new double[3];

	if (iind == 1)
		ind3 = 3;
	else
		ind3 = iind - 1;

	v2[0] = bases[p1].x[iind] - bases[p1].xcom;
	v2[1] = bases[p1].y[iind] - bases[p1].ycom;
	v2[2] = bases[p1].z[iind] - bases[p1].zcom;
	v2[0] -= plist.xboxl * round(v2[0] / plist.xboxl);
	v2[1] -= plist.yboxl * round(v2[1] / plist.yboxl);
//	v2[2] -= plist.zboxl * round(v2[2] / plist.zboxl);

	v1[0] = bases[p1].x[ind3] - bases[p1].xcom;
	v1[1] = bases[p1].y[ind3] - bases[p1].ycom;
	v1[2] = bases[p1].z[ind3] - bases[p1].zcom;
	v1[0] -= plist.xboxl * round(v1[0] / plist.xboxl);
	v1[1] -= plist.yboxl * round(v1[1] / plist.yboxl);
//	v1[2] -= plist.zboxl * round(v1[2] / plist.zboxl);

	/*2 cross 1 will be normal up*/
	crossproduct(v2, v1, u1);
	cout << "normal 1: " << u1[0] << ' ' << u1[1] << ' ' << u1[2] << endl;
	v2[0] = bases[p2].x[2] - bases[p2].xcom;
	v2[1] = bases[p2].y[2] - bases[p2].ycom;
	v2[2] = bases[p2].z[2] - bases[p2].zcom;
	v2[0] -= plist.xboxl * round(v2[0] / plist.xboxl);
	v2[1] -= plist.yboxl * round(v2[1] / plist.yboxl);
//	v2[2] -= plist.zboxl * round(v2[2] / plist.zboxl);

	v1[0] = bases[p2].x[1] - bases[p2].xcom;
	v1[1] = bases[p2].y[1] - bases[p2].ycom;
	v1[2] = bases[p2].z[1] - bases[p2].zcom;
	v1[0] -= plist.xboxl * round(v1[0] / plist.xboxl);
	v1[1] -= plist.yboxl * round(v1[1] / plist.yboxl);
//	v1[2] -= plist.zboxl * round(v1[2] / plist.zboxl);

	crossproduct(v2, v1, u2);
	cout << "normal 2: " << u2[0] << ' ' << u2[1] << ' ' << u2[2] << endl;
	//dotproduct(u2,u1, theta);

	dp = u1[0] * u2[0] + u1[1] * u2[1] + u1[2] * u2[2];
	l1 = u1[0] * u1[0] + u1[1] * u1[1] + u1[2] * u1[2];
	l2 = u2[0] * u2[0] + u2[1] * u2[1] + u2[2] * u2[2];
	dpnorm = dp / sqrt(l1 * l2);
	if (dpnorm > 1) {
		cout << "dp out of range! " << dpnorm << endl;
		theta = 0;
	} else if (dpnorm < -1) {
		cout << "dp out of range! " << dpnorm << endl;
		theta = M_PI;
	} else {
		theta = acos(dpnorm);
	}

	trot = theta / 2.0;
	double *M2 = new double[9];
	double *M2neg = new double[9];
	/*p1 rotate by positive theta/2, p2 rotate by negative theta/2
	 rotate around the bound leg segment
	 */
	crossproduct(u1, u2, u);
	ux = u[0];
	uy = u[1];
	uz = u[2];
	cthet = cos(+trot);
	sthet = sin(+trot);
	cout << "second rotation angle between: " << trot << " normal: " << ux << ' ' << uy << ' ' << uz << endl;
	M2[0] = cthet + ux * ux * (1 - cthet);
	M2[1] = ux * uy * (1 - cthet) - uz * sthet; //row 0, column 1, go across first!
	M2[2] = ux * uz * (1 - cthet) + uy * sthet;
	M2[3] = uy * ux * (1 - cthet) + uz * sthet;
	M2[4] = cthet + uy * uy * (1 - cthet);
	M2[5] = uy * uz * (1 - cthet) - ux * sthet;
	M2[6] = uz * ux * (1 - cthet) - uy * sthet;
	M2[7] = uz * uy * (1 - cthet) + ux * sthet;
	M2[8] = cthet + uz * uz * (1 - cthet);

	M2neg[0] = M2[0];
	M2neg[1] = ux * uy * (1 - cthet) + uz * sthet;
	M2neg[2] = ux * uz * (1 - cthet) - uy * sthet;
	M2neg[3] = uy * ux * (1 - cthet) - uz * sthet;
	M2neg[4] = M2[4];
	M2neg[5] = uy * uz * (1 - cthet) + ux * sthet;
	M2neg[6] = uz * ux * (1 - cthet) + uy * sthet;
	M2neg[7] = uz * uy * (1 - cthet) - ux * sthet;
	M2neg[8] = M2[8];

	rotate_onlyPBCCELL(p1, c1, ind_com, bases, M2, plist);
	rotate_onlyPBCCELL(p2, c2, ind_com, bases, M2neg, plist);

	/*Determine if you crashed two clathrins together
	 unless they are in the same complex, in which case you are closing a loop.
	 */
	int cancel;
	if (c1 != c2) {
		cancel = measure_overlap(c1, c2, ind_com, bases, plist.pclath);
	} else
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

		/*add c2's proteins to c1's list, unless closing a loop*/
		/*Get new COM of complex*/
		/*Copy final complex into the spot of now gone complex c2 */
		int mp;
		int j;
		int tar;
		int flagbndry = 0;
		if (c1 == c2) {
			cout << "CLOSING A LOOP! " << endl;
			plist.nloop++;
			if (ind_com[c1].radR * 2.0 > plist.xboxl / 2.0)
				flagbndry = 1;
			update_one_com_onlyPBCCELL(c1, ind_com, bases, plist, flagbndry);

			update_radiusPBCCELL(c1, ind_com, bases, plist); //requires correct ind_com
			update_trans_diffusion(c1, ind_com, bases, plist.pretrans);
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
			if ((ind_com[c1].radR + ind_com[c2].radR) * 2.0 > plist.xboxl / 2.0)
				flagbndry = 1;
			update_one_com_onlyPBCCELL(c1, ind_com, bases, plist, flagbndry);
			update_radiusPBCCELL(c1, ind_com, bases, plist); //requires correct ind_com
			update_trans_diffusion(c1, ind_com, bases, plist.pretrans);
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
		/*reverse previous rotation*/
		rotate_onlyPBCCELL(p1, c1, ind_com, bases, M2neg, plist);
		rotate_onlyPBCCELL(p2, c2, ind_com, bases, M2, plist);

		dtrans[0] *= -1;
		dtrans[1] *= -1;
		dtrans[2] *= -1;

		drev[0] *= -1;
		drev[1] *= -1;
		drev[2] *= -1;

		/*Update the position of proteins in complex one*/
		rotate_and_translate_intPBCCELL(p1, c1, ind_com, bases, M, dtrans, iind, plist);

		/*Update position of proteins in c2*/
		rotate_and_translate_intPBCCELL(p2, c2, ind_com, bases, Mneg, drev, iind2, plist);

	}

	cout << "final complex com: " << ind_com[c1].xcom << ' ' << ind_com[c1].ycom << ' ' << ind_com[c1].zcom << " radius: " << ind_com[c1].radR << " Dr: " << ind_com[c1].Drx << " Dtrans: " << ind_com[c1].Dx << endl;

	delete[] u1;
	delete[] u2;
	delete[] v1;
	delete[] v2;
	delete[] u;
	delete[] M;
	delete[] Mneg;
	delete[] dtrans;
	delete[] drev;

}
