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

int break_complexPBCCELLMEM(int p1, int mu, int kind, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double *bindrad, int p2, int i1, int i2, int *p_home, int **myrxn,Protein *wholep) {

	/*Dissociate two bound interfaces, make sure all the proteins each is bound to
	 get put in the correct final complex.
	 */

	int c1 = bases[p1].mycomplex;
	int s1 = ind_com[c1].mysize;
	int prod = Rlist[mu][0];

	//change status of the interface
	int iind = ihome[i1];
	int iind2 = ihome[i2];
	int i, k, j, membnd;
	int mp, s2;
	membnd = 0;

	/*assign each protein in original complex c1 to one of the two new complexes,
	 if the complex forms a loop, they will be put back together in c1, and the
	 individual interfaces that dissociated freed.
	 */
	int cancel = 0;
	/*continue on with the dissociation that creates two complexes*/

	bases[p1].istatus[iind] = i1;
	bases[p1].nbnd -= 1;
	bases[p1].bndlist[kind] = bases[p1].bndlist[bases[p1].nbnd];

	/*Add these protein into the bimolecular list*/
	bases[p1].nfree += 1;
	bases[p1].freelist[bases[p1].nfree - 1] = i1;

	/*create a new complex to hold partner and his partners*/
	//		plist.ntotalcomplex++;

	/*Now we have the correct lists of proteins in complex 1 and complex2,
	 */
//	update_diffusion(c1, ind_com, bases);

	//find the total number of lipids attached to c1
	for (i = 0; i < s1; i++) {
		p1 = ind_com[c1].plist[i];
		for (j=0; j<(bases[p1].nbnd); j++){
			if (bases[p1].partner[j]==-1){
				membnd += 1;
			}
		}
	}

	int size = ind_com[c1].mysize;
	double Dxinv = 0;
	double Dyinv = 0;
	double Dzinv = 0;
	double inf = 1E500;
	p1 = ind_com[c1].plist[0];
	ind_com[c1].Dx = bases[p1].Dx;
	ind_com[c1].Dy = bases[p1].Dy;
	ind_com[c1].Dz = bases[p1].Dz;
	for (i = 0; i < size; i++) {

		p1 = ind_com[c1].plist[i];
		//Dxinv+=1.0/bases[p1].Dx;
		//Dyinv+=1.0/bases[p1].Dy;

		if (bases[p1].Dx != 0) {
			Dxinv += 1.0 / bases[p1].Dx;
		} else {
			Dxinv = inf;
			//      ind_com[c1].Dx=0;
		}
		if (bases[p1].Dy != 0) {
			Dyinv += 1.0 / bases[p1].Dy;
		} else {
			Dyinv = inf;
			//ind_com[c1].Dy=0;
		}
		if (bases[p1].Dz != 0) {
			Dzinv += 1.0 / bases[p1].Dz;
		} else {
			Dzinv = inf;
			//ind_com[c1].Dz=0;
		}
	}
	for (j=0;j<membnd;j++){
		Dxinv += 1.0 / wholep[0].Dx;
		Dyinv += 1.0 / wholep[0].Dy;
		Dzinv = inf;
	}
	ind_com[c1].Dx = 1.0 / Dxinv;
	ind_com[c1].Dy = 1.0 / Dyinv;
	ind_com[c1].Dz = 1.0 / Dzinv;

	int flagbndry = 0;
	if (ind_com[c1].radR * 2.0 / plist.xboxl > 1)
		flagbndry = 1;
//	update_one_com_onlyPBCCELL(c1, ind_com, bases, plist, flagbndry);
	//update_com_only(c2, ind_com, bases);//will have to do this again after moving

	double stretch;
	double addx;
	double addy;
	double addz;
	double *chg1 = new double[3];
	double *chg2 = new double[3];

	double R = bindrad[mu];
	double Ry, R2, axe2;
	double sign = 1;
	double Dxtot = ind_com[c1].Dx + wholep[0].Dx;
	double Dytot = ind_com[c1].Dy + wholep[0].Dy;
	double Dztot = ind_com[c1].Dz + wholep[0].Dz;
	double tol = 1E-16;
	int dzflag = 0;
	if (Dztot < tol) {
		Dztot = 1;
		dzflag = 1;
	}
	if (Dxtot < tol)
		Dxtot = 1;
	if (Dytot < tol)
		Dytot = 1;

	/*Now, displace the complexes to the contact radius.
				 Move them apart based on the COM->interface vector,
				 otherwise if that vector is zero, sample positions uniformly on the sphere of diameter R.
	 */

	double dx = -bases[p1].x[iind] + ind_com[c1].xcom;
	double dy = -bases[p1].y[iind] + ind_com[c1].ycom;
	double dz = -bases[p1].z[iind] + ind_com[c1].zcom;
	dx -= plist.xboxl * round(dx / plist.xboxl);
	dy -= plist.yboxl * round(dy / plist.yboxl);
	//		dz -= plist.zboxl * round(dz / plist.zboxl);

	double vlen2, vlen;
	tol = 1E-13;
	if (dzflag == 1) {
		vlen2 = dx * dx + dy * dy;
		vlen = sqrt(vlen2);
	} else {
		vlen2 = dx * dx + dy * dy + dz * dz;
		vlen = sqrt(vlen2);
	}
	//cout <<"Breaking complex, current sep! : "<<vlen<<endl;
	double small = 1E-9;
	double rnum, phi, st, theta;
	if (vlen > tol) {

		addx = dx * R / vlen;
		if (addx > 0)
			addx += small;
		else
			addx -= small; //so final separation is not binding radius with precision issues
		addy = dy * R / vlen;
		addz = dz * R / vlen;
	} else {
		/*To sample uniformly on the sphere, choose
					 polar angle from -1:1, uniform in cos(theta),
					 and azimuthal uniform in 0:2pi.
		 */

		phi = rand_gsl() * 2.0 * M_PI; //between 0 and 2pi

		if (dzflag == 1) {
			rnum = 0; //polar angle
			theta = M_PI / 2.0;
			addx = R * cos(phi); //sin(theta)=1
			if (addx > 0)
				addx += small;
			else
				addx -= small;
			addy = R * sin(phi); //sin(theta)=1
			addz = 0; //cos(theta)=0
		} else {
			rnum = rand_gsl() * 2.0 - 1.0; //between -1 and 1 (pi:0)
			theta = acos(rnum);
			st = sin(theta);
			addx = R * st * cos(phi);
			if (addx > 0)
				addx += small;
			else
				addx -= small; //so final separation is not exactly binding radius.
			addy = R * st * sin(phi);
			addz = R * rnum; //rnum=cos(theta)
		}
	}

	/*The separation vector is of length R, move them in opposite directions based on their fraction
				 of the total diffusion, e.g. if D1=D2, they both move halfway in opposite directions
	 */
	chg1[0] = addx * ind_com[c1].Dx / Dxtot;
	chg1[1] = addy * ind_com[c1].Dy / Dytot;
	chg1[2] = addz * ind_com[c1].Dz / Dztot;

	move_pro_coordsPBCCELL(c1, ind_com, bases, chg1, plist);
	update_one_com_onlyPBCCELL(c1, ind_com, bases, plist, flagbndry);
	update_radiusPBCCELL(c1, ind_com, bases, plist);
	update_rot_diffusion(c1, ind_com, bases, plist.prerot);

	delete[] chg1;
	delete[] chg2;

	return cancel;

}
