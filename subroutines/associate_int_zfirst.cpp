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

void associate_int_zfirst(int p1, int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist) {

	/*Move the interfaces that are reacting with one another to 
	 the position geometrically between them, but move each
	 interface along this vector proportional to their diffusion constant
	 In this version, put the two interfaces into the same z-plane before calculating
	 the vector between them so they are not rotated around z.
	 
	 This version was implemented for the case that clathrin is fixed
	 in the plane, or when any molecule's z-orientation was fixed.
	 */
	int prod = Rlist[mu][2];
	int iind = ihome[i1];
	int iind2 = ihome[i2];

	/*change status of the interface from the free state to the bound state*/
	bases[p1].istatus[iind] = prod;
	bases[p2].istatus[iind2] = prod;

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

	/*Add the proteins as partners at their specific interfaces*/
	bases[p1].partner[iind] = p2;
	bases[p2].partner[iind2] = p1;

	/*store values for each protein's original pre-associating complex*/
	int c1 = bases[p1].mycomplex;
	int s1 = ind_com[c1].mysize;
	int c2 = bases[p2].mycomplex;
	int s2 = ind_com[c2].mysize;
	int newsize = s1 + s2;

	double *v0 = new double[3];
	double *v2 = new double[3];
	double *M = new double[9];
	double *u = new double[3];
	double *dtrans = new double[3];
	double *drev = new double[3];

	/*Now move proteins*/
	double Dxtot = ind_com[c1].Dx + ind_com[c2].Dx;
	double Dytot = ind_com[c1].Dy + ind_com[c2].Dy;
	double Dztot = ind_com[c1].Dz + ind_com[c2].Dz;

	/*For proteins bound to membrane, Dztot could be zero*/
	double tol = 1E-16;
	if (Dztot < tol)
		Dztot = 1; //otherwise you divide by zero 

	double dz = bases[p2].z[iind2] - bases[p1].z[iind];

	/*eliminate the differences in their z-position*/
	double delz1 = dz * ind_com[c1].Dz / Dztot;
	double delz2 = -dz * ind_com[c2].Dz / Dztot;

	move_zcrds(bases, p1, p2, delz1, delz2, ind_com);

	/*The sign ensures that the COMs of the protein get rotated 
	 to the outside of interfaces final spot*/
	v0[0] = -bases[p2].x[iind2] + bases[p1].x[iind];
	v0[1] = -bases[p2].y[iind2] + bases[p1].y[iind];
	v0[2] = -bases[p2].z[iind2] + bases[p1].z[iind]; //This should be zero now!
	v0[2] = 0; //it is possible this is not zero if both associating molecules were bound to the membrane already

	double theta;
	double pivx;
	double pivy;
	double pivz;

	//distance to move to place interfaces, along vector v 
	dtrans[0] = -v0[0] * ind_com[c1].Dx / Dxtot;
	dtrans[1] = -v0[1] * ind_com[c1].Dy / Dytot;
	dtrans[2] = -v0[2] * ind_com[c1].Dz / Dztot; //should be zero

	drev[0] = +v0[0] * ind_com[c2].Dx / Dxtot;
	drev[1] = +v0[1] * ind_com[c2].Dy / Dytot;
	drev[2] = +v0[2] * ind_com[c2].Dz / Dztot; //should be zero

	/*Below we have it set so that if the protein is single, it is rotated
	 and translated to the meeting point.
	 If the protein is part of a larger complex, then the complex is simply moved
	 to associate with the partner without rotating around the COM.
	 For single site proteins these are the same displacements
	 */
	if (s1 == 1) {
		pivx = bases[p1].x[iind];
		pivy = bases[p1].y[iind];
		pivz = bases[p1].z[iind];

		/*First calculate the angle between the central vector v and the vector from interface to center*/
		/*assume the center of mass is in the same z plane as the interface, because of rigid molecules*/
		v2[0] = -pivx + bases[p1].xcom;
		v2[1] = -pivy + bases[p1].ycom;
		v2[2] = 0; //-pivz+bases[p1].zcom;

		//  cout <<"starting vector com-p1: "<<v2[0]<<' '<<v2[1]<<' '<<v2[2]<<endl;

		/*Calculate rotation matrix*/

		dotproduct(v2, v0, theta);
		//now calculate the axis of rotation
		crossproduct(v2, v0, u); //u is the UNIT vector of the rotation axis
		//calculate the rotation matrix for rotating theta around u
		calc_Rmatrix(u, theta, M);

		/*Update the position of proteins in complex one*/
		rotate_and_translate_int(p1, c1, ind_com, bases, M, dtrans, iind);

	} else {
		//just translate
		translate_int(p1, c1, ind_com, bases, dtrans);
	}

	/*Same for complex2*********************************/

	v0[0] *= -1;
	v0[1] *= -1;
	v0[2] *= -1;

	if (s2 == 1) {
		/*First calculate the angle between the central vector v and the vector from interface to center*/
		pivx = bases[p2].x[iind2];
		pivy = bases[p2].y[iind2];
		pivz = bases[p2].z[iind2];

		/*again, ensure that vector from interface to com has zero z component */
		v2[0] = -pivx + bases[p2].xcom;
		v2[1] = -pivy + bases[p2].ycom;
		v2[2] = 0; //-pivz+bases[p2].zcom;

		/*Calculate Rotation matrix*/
		dotproduct(v2, v0, theta);
		crossproduct(v2, v0, u); //u is the UNIT vector of the rotation axis
		//calculate the rotation matrix for rotating theta around u
		calc_Rmatrix(u, theta, M);

		/*Update position of proteins in c2*/
		rotate_and_translate_int(p2, c2, ind_com, bases, M, drev, iind2);
	} else {
		//just translate
		translate_int(p2, c2, ind_com, bases, drev);
	}

	/*Combine the complexes into one complex by listing all c2 proteins into c1, 
	 unless you are closing a loop, which could happen for clathrin cages.
	 Update the complex COM and diffusion constants
	 */
	int mp;
	double totalmassx = 0;
	double totalmassy = 0;
	double totalmassz = 0;
	ind_com[c1].xcom = 0;
	ind_com[c1].ycom = 0;
	ind_com[c1].zcom = 0;

	int j, tar;
	if (c1 == c2) {
		cout << "CLOSING A LOOP! " << endl;
		ind_com[c1].mysize = s1;
		/*Size of complex does not change*/
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

		update_diffusion(c1, ind_com, bases);
		update_radius(c1, ind_com, bases);
	} else {
		int t = s1;
		for (i = 0; i < s2; i++) {
			mp = ind_com[c2].plist[i];
			ind_com[c1].plist[t] = mp; //add these proteins to the first complex
			t++;
			bases[mp].mycomplex = c1;
		}

		ind_com[c1].mysize = newsize;

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

		update_diffusion(c1, ind_com, bases);
		update_radius(c1, ind_com, bases);

		/*Complex c2 is merged with c1, so copy the last complex into the place
		 of c2.
		 */
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

	delete[] v0;
	delete[] v2;
	delete[] u;
	delete[] M;
	delete[] dtrans;
	delete[] drev;

}
