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



int check_overlap(double &dr2, int p1, int p2, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, Parms plist, int *ihome, double *bindrad, double **trajR, double *M, int **Rlist) {
	/*
	  In this version, complex k1 is on the membrane. 
	  If both proteins are on the membrane (Dz==0), evaluate only xy displacement, not z.

	  In this _complex_ version, it tests overlap not just for each protein, but for each complex, so all the proteins in a complex, before performing
	  position updates.

	  This routine checks whether protein p1 is overlapping p2 in its reaction
	  zone. Enforce PBC in distance measurements.
	 */

	int i, j, p2test, flag=0, jj; 
	double *M2 = new double[9];
	double *v = new double[3];
	double *v2 = new double[3];
	double df1, df2, df3;
	int k1, k2, k, proi, proj;
	double dx1, dx2;
	double dy1, dy2;
	double dz1, dz2;
	//	int p2;
	int i1, i2;
	int rxn;
	int iind1, iind2;
//	double dr2;
	dr2 = 0;
	k1 = bases[p1].mycomplex;
	k2 = bases[p2].mycomplex;

	double x0 = ind_com[k1].xcom;
	double y0 = ind_com[k1].ycom;
	double z0 = ind_com[k1].zcom;
	rotationEuler(trajR[k1][0], trajR[k1][1], trajR[k1][2], M); //M is always for protein p1
	double x02;
	double y02;
	double z02;

	for (i = 0; i < ncross[p1]; i++) {
		p2test = cross_part[p1][i];
		if(p2test==p2){
			i1 = cross_int[p1][i];
			rxn = cross_rxn[p1][i];
			/// find what i2 is, p1 and p2 might be interacting through multiple interfaces
			if (Rlist[rxn][0] == i1)
				i2 = Rlist[rxn][1]; //this is the interface we're looking for on the other protein
			else
				i2 = Rlist[rxn][0];

			iind1 = ihome[i1];
			iind2 = ihome[i2];

			v[0] = bases[p1].x[iind1] - x0;
			v[1] = bases[p1].y[iind1] - y0;
			v[2] = bases[p1].z[iind1] - z0;
			v[0] -= plist.xboxl * round(v[0] / plist.xboxl);
			v[1] -= plist.yboxl * round(v[1] / plist.yboxl);
			//			v[2] -= plist.zboxl * round(v[2] / plist.zboxl);

			rotate(v, M, v2);
			dx1 = x0 + v2[0] + traj[k1][0];
			dy1 = y0 + v2[1] + traj[k1][1];
			dz1 = z0 + v2[2] + traj[k1][2];

			/*Now complex 2*/

			rotationEuler(trajR[k2][0], trajR[k2][1], trajR[k2][2], M2);
			reflect_traj_complex_rad_rotCELL(p2, bases, ind_com, plist.xboxl, plist.yboxl, plist.zboxl, traj, M2);
			x02 = ind_com[k2].xcom;
			y02 = ind_com[k2].ycom;
			z02 = ind_com[k2].zcom;

			v[0] = bases[p2].x[iind2] - x02;
			v[1] = bases[p2].y[iind2] - y02;
			v[2] = bases[p2].z[iind2] - z02;
			v[0] -= plist.xboxl * round(v[0] / plist.xboxl);
			v[1] -= plist.yboxl * round(v[1] / plist.yboxl);
			// v[2] -= plist.zboxl * round(v[2] / plist.zboxl);

			rotate(v, M2, v2);

			dx2 = x02 + v2[0] + traj[k2][0];
			dy2 = y02 + v2[1] + traj[k2][1];
			dz2 = z02 + v2[2] + traj[k2][2];

			/*separation*/
			df1 = dx1 - dx2;
			df2 = dy1 - dy2;
			df3 = dz1 - dz2;
			df1 -= plist.xboxl * round(df1 / plist.xboxl);
			df2 -= plist.yboxl * round(df2 / plist.yboxl);
			//df3 -= plist.zboxl * round(df3 / plist.zboxl);
			if(ind_com[k1].Dz==0 && ind_com[k2].Dz==0){
				dr2 = df1 * df1 + df2 * df2; //xy displacement only
			}else
				dr2 = df1 * df1 + df2 * df2 + df3 * df3;

			if (dr2 < bindrad[rxn] * bindrad[rxn]) {
				flag = 1;
			}
		}
	}
	delete[] M2;
	delete[] v;
	delete[] v2;

	return flag;

}
