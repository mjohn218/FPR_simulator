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

int break_complex_trans(int p1, int mu, int kind, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double *bindrad, int p2, int i1, int i2, int *p_home, int **myrxn) {

	/*Dissociate two bound interfaces, make sure all the proteins each is bound to
	 get put in the correct final complex.
	 */

	int c1 = bases[p1].mycomplex;
	int s1 = ind_com[c1].mysize;
	int prod = Rlist[mu][0];

	/*change status of the interface*/
	int iind = ihome[i1];
	int iind2 = ihome[i2];
	int i, k;
	int mp, s2;
	int c2 = plist.ntotalcomplex;
	bases[p2].mycomplex = c2;
	ind_com[c2].mysize = 0; //haven't added anyone yet

	int cancel = 0;
	/*assign each protein in original complex c1 to one of the two new complexes,
	 if the complex forms a loop, they will be put back together in c1, and the 
	 individual interfaces that dissociated freed.
	 */
	//if (bases[p1].protype == bases[p2].protype && bases[p1].protype == plist.pclath)
	cancel = determine_which_complex_double(p1, p2, ind_com, c1, c2, bases, myrxn, Rlist, ihome, p_home, plist);
	//else
	//cancel = determine_which_complex_merge(p1, p2, ind_com, c1, c2, bases, myrxn, Rlist, ihome, p_home, plist);

	if (cancel == 0) {

		/*continue on with the dissociation that creates two complexes*/

		bases[p1].istatus[iind] = i1;
		bases[p1].nbnd -= 1;
		bases[p1].bndlist[kind] = bases[p1].bndlist[bases[p1].nbnd];

		bases[p2].istatus[iind2] = i2; //used to be the product state
		bases[p2].nbnd -= 1;

		for (k = 0; k < bases[p2].nbnd + 1; k++) {
			if (bases[p2].bndlist[k] == prod) {
				bases[p2].bndlist[k] = bases[p2].bndlist[bases[p2].nbnd]; //replace in list with last
				k = bases[p2].nbnd + 1; //only replace one of the products, could have multiple
			}
		}

		/*Add these protein into the bimolecular list*/
		bases[p1].nfree += 1;
		bases[p2].nfree += 1;
		bases[p1].freelist[bases[p1].nfree - 1] = i1;
		bases[p2].freelist[bases[p2].nfree - 1] = i2;

		/*create a new complex to hold partner and his partners*/
		plist.ntotalcomplex++;

		/*Now we have the correct lists of proteins in complex 1 and complex2 */
		update_diffusion(c1, ind_com, bases);
		update_diffusion(c2, ind_com, bases);

		update_one_com_only(c1, ind_com, bases); //get COM of the c1 complex so the two complexes are moved apart
		//update_com_only(c2, ind_com, bases);//will have to do this again after moving

		double stretch;
		double addx;
		double addy;
		double addz;
		double *chg1 = new double[3];
		double *chg2 = new double[3];

		double R = bindrad[mu]; //contact separation
		double Ry, R2, axe2;
		double sign = 1;
		double Dxtot = ind_com[c1].Dx + ind_com[c2].Dx;
		double Dytot = ind_com[c1].Dy + ind_com[c2].Dy;
		double Dztot = ind_com[c1].Dz + ind_com[c2].Dz;
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
		double vlen2, vlen;

		if (dzflag == 1) {
			vlen2 = dx * dx + dy * dy;
			vlen = sqrt(vlen2);
		} else {
			vlen2 = dx * dx + dy * dy + dz * dz;
			vlen = sqrt(vlen2);
		}
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
					addx -= small; //so final separation is not binding radius with precision issues
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
					addx -= small; //so final separation is not binding radius with precision issues
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
		chg2[0] = -addx * ind_com[c2].Dx / Dxtot;
		chg2[1] = -addy * ind_com[c2].Dy / Dytot;
		chg2[2] = -addz * ind_com[c2].Dz / Dztot;
		//    cout <<"displacement, p1: "<<addx1<<' '<<addy1<<' '<<addz1<<endl;
		//cout <<"displacement, p2: "<<addx2<<' '<<addy2<<' '<<addz2<<endl;

		/*Update the positions of the protein complexes and the COM of the complexes
		 then update the radius of each complex.
		 */
		move_compro_coords(c1, ind_com, bases, chg1);
		move_compro_coords(c2, ind_com, bases, chg2);
		update_radius(c1, ind_com, bases);
		update_radius(c2, ind_com, bases);

		update_rot_diffusion(c1, ind_com, bases, plist.prerot);
		update_rot_diffusion(c2, ind_com, bases, plist.prerot);

		delete[] chg1;
		delete[] chg2;

	} else {
		/*Reset all proteins back to complex c1, dissociation
		 will break the product state of the two proteins that dissociated but if they
		 are linked in a closed loop it will not create a new complex.
		 */

		s1 = ind_com[c1].mysize;
		s2 = ind_com[c2].mysize;
		for (i = 0; i < s2; i++) {
			mp = ind_com[c2].plist[i];
			bases[mp].mycomplex = c1;
			ind_com[c1].plist[s1] = mp;
			s1++;
		}
		ind_com[c1].mysize = s1;

		/*Put interfaces that dissociated into free state*/
		bases[p1].istatus[iind] = i1;
		bases[p1].nbnd -= 1;
		bases[p1].bndlist[kind] = bases[p1].bndlist[bases[p1].nbnd];

		bases[p2].istatus[iind2] = i2; //used to be the product state
		bases[p2].nbnd -= 1;

		for (k = 0; k < bases[p2].nbnd + 1; k++) {
			if (bases[p2].bndlist[k] == prod) {
				bases[p2].bndlist[k] = bases[p2].bndlist[bases[p2].nbnd]; //replace in list with last
				k = bases[p2].nbnd + 1; //only replace one of the products, could have multiple
			}
		}

		/*Add these protein into the bimolecular list*/
		bases[p1].nfree += 1;
		bases[p2].nfree += 1;
		bases[p1].freelist[bases[p1].nfree - 1] = i1;
		bases[p2].freelist[bases[p2].nfree - 1] = i2;

	}

	return cancel;

}
