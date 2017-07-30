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

void associate_zsigmaPBC(int p1, int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom) {

	/*Move the interfaces that are reacting with one another to 
	 the position geometrically between them, but move each
	 interface along this vector proportional to their diffusion constant
	 In this version, put the two interfaces into the same z-plane before calculating
	 the vector between them so they are not rotated around z.
	 Enforce PBC.
	 */
	/*NOW WE ARE MOVING THE TWO PROTEINS TO A SEPARATION OF SIGMA, when they dissociate, they will not move*/

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

	/*Add each protein as the partner at the reacting interfaces*/
	bases[p1].partner[iind] = p2;
	bases[p2].partner[iind2] = p1;

	/*Store each protein's pre-association complex values*/
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

	/*Now move protein*/
	double Dxtot = ind_com[c1].Dx + ind_com[c2].Dx;
	double Dytot = ind_com[c1].Dy + ind_com[c2].Dy;
	double Dztot = ind_com[c1].Dz + ind_com[c2].Dz;

	/*For proteins bound to membrane, Dztot could be zero*/
	double tol = 1E-16;
	if (Dztot < tol)
		Dztot = 1; //otherwise you divide by zero 

	double dz = bases[p2].z[iind2] - bases[p1].z[iind];
	dz -= plist.zboxl * round(dz / plist.zboxl);

	/*eliminate the differences in their z-position*/
	double delz1 = dz * ind_com[c1].Dz / Dztot;
	double delz2 = -dz * ind_com[c2].Dz / Dztot;

	move_zcrdsPBC(bases, p1, p2, delz1, delz2, ind_com, plist);

	/*The sign ensures the COM's get rotated
	 to the outside of interfaces final spot*/
	v0[0] = -bases[p2].x[iind2] + bases[p1].x[iind];
	v0[0] -= plist.xboxl * round(v0[0] / plist.xboxl);
	v0[1] = -bases[p2].y[iind2] + bases[p1].y[iind];
	v0[1] -= plist.yboxl * round(v0[1] / plist.yboxl);
	v0[2] = -bases[p2].z[iind2] + bases[p1].z[iind]; //This should be zero now!
	v0[2] = 0; //it is possible this is not zero if both associating molecules were bound to the membrane already
	double R2 = v0[0] * v0[0] + v0[1] * v0[1];
	double R1 = sqrt(R2);

	double theta;
	double pivx;
	double pivy;
	double pivz;
	double small = 1E-9;
	//distance to move to place interfaces, along vector v 
	/*NOW WE ARE MOVING THE TWO PROTEINS TO A SEPARATION OF SIGMA*/
	dtrans[0] = v0[0] * ind_com[c1].Dx / Dxtot * (bindrad / R1 - 1.0);
	dtrans[1] = v0[1] * ind_com[c1].Dy / Dytot * (bindrad / R1 - 1.0);
	dtrans[2] = v0[2] * ind_com[c1].Dz / Dztot * (bindrad / R1 - 1.0); //should be zero

	drev[0] = -v0[0] * ind_com[c2].Dx / Dxtot * (bindrad / R1 - 1.0);
	drev[1] = -v0[1] * ind_com[c2].Dy / Dytot * (bindrad / R1 - 1.0);
	drev[2] = -v0[2] * ind_com[c2].Dz / Dztot * (bindrad / R1 - 1.0); //should be zero

	//just translate if larger complexes are coming together
	translate_intPBC(p1, c1, ind_com, bases, dtrans, plist);

	/*Same for complex2*********************************/
	//just translate
	translate_intPBC(p2, c2, ind_com, bases, drev, plist);

	//  cout <<"Initial separation: :"<<R1<<'\t';
	v0[0] = -bases[p2].x[iind2] + bases[p1].x[iind];
	v0[0] -= plist.xboxl * round(v0[0] / plist.xboxl);
	v0[1] = -bases[p2].y[iind2] + bases[p1].y[iind];
	v0[1] -= plist.yboxl * round(v0[1] / plist.yboxl);
	v0[2] = -bases[p2].z[iind2] + bases[p1].z[iind]; //This should be zero now!
	v0[2] = 0; //it is possible this is not zero if both associating molecules were bound to the membrane already
	R2 = v0[0] * v0[0] + v0[1] * v0[1];
	R1 = sqrt(R2);
	// cout <<"final sep: "<<R1<<endl;
//   cout <<"final crds: "<<bases[p1].xcom<<' '<<bases[p1].ycom<<' '<<bases[p1].zcom<<endl;
//   cout <<"final crds: "<<bases[p2].xcom<<' '<<bases[p2].ycom<<' '<<bases[p2].zcom<<endl;

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

		ind_com[c1].xcom -= plist.xboxl * round(ind_com[c1].xcom / plist.xboxl);
		ind_com[c1].ycom -= plist.yboxl * round(ind_com[c1].ycom / plist.yboxl);
		ind_com[c1].zcom -= plist.zboxl * round(ind_com[c1].zcom / plist.zboxl);

		update_diffusion(c1, ind_com, bases);
		update_radius(c1, ind_com, bases);
	} else {

		/*Move all c2 proteins into c1's list*/
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

		ind_com[c1].xcom -= plist.xboxl * round(ind_com[c1].xcom / plist.xboxl);
		ind_com[c1].ycom -= plist.yboxl * round(ind_com[c1].ycom / plist.yboxl);
		ind_com[c1].zcom -= plist.zboxl * round(ind_com[c1].zcom / plist.zboxl);

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
		  ncrosscom[c2] = ncrosscom[tar];
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
	//  cout <<"ind_com crds: "<<ind_com[c1].xcom<<' '<<ind_com[c1].ycom<<' '<<ind_com[c1].zcom<<endl;
	delete[] v0;
	delete[] v2;
	delete[] u;
	delete[] M;
	delete[] dtrans;
	delete[] drev;

}
