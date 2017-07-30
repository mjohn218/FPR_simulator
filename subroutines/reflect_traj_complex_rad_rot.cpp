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

void reflect_traj_complex_rad_rot(int p1, Fullmol *bases, Complex *ind_com, double xboxl, double yboxl, double zboxl, double **traj, double *M) {
	int i, j;
	int k = bases[p1].mycomplex;
	int s1 = ind_com[k].mysize;
	double xchg;
	double ychg;
	double zchg;
	double xtot = 0.0;
	double ytot = 0.0;
	double ztot = 0.0;
	int flagx = 0;
	int flagy = 0;
	double rad = ind_com[k].radR;
	int flagz = 0;
	double currx = ind_com[k].xcom + traj[k][0];
	double curry = ind_com[k].ycom + traj[k][1];
	double currz = ind_com[k].zcom + traj[k][2];
	double xtot0, ytot0, ztot0;
	double tol = 1E-11;
	xchg = currx - xboxl / 2.0;

	if ((xchg + rad) > 0) {
		xtot0 = -(xchg + rad); //shift x coordinates back
		flagx = 1;
	} else if ((xchg - rad) < -xboxl) {
		xtot0 = -(currx - rad + xboxl / 2.0);
		flagx = 1;
	}
	ychg = curry - yboxl / 2.0;
	if ((ychg + rad) > 0) {
		ytot0 = -(ychg + rad); //shift x coordinates back
		flagy = 1;
	} else if ((ychg - rad) < -yboxl) {
		ytot0 = -(curry - rad + yboxl / 2.0);
		flagy = 1;
	}
	zchg = currz - zboxl / 2.0;
	if ((zchg + rad) > 0) {
		ztot0 = -(zchg + rad); //shift x coordinates back
		flagz = 1;
	} else if ((zchg - rad) < -zboxl) {
		ztot0 = -(currz - rad + zboxl / 2.0);
		flagz = 1;
	}

	int mp;

	/*Z is separate to allow the interfaces to approach to the membrane
	 but don't need to test if the entire complex is far enough
	 away from the boundary.
	 */

	double dx, dy, dz;
	double x0, y0, z0;
	double dzrot, dxrot, dyrot;
	double row[3];
	int flag;
	if (flagx == 1) {
		flag = 0;
		xtot = 0;
		x0 = ind_com[k].xcom;
		y0 = ind_com[k].ycom;
		z0 = ind_com[k].zcom;
		row[0] = M[0];
		row[1] = M[1];
		row[2] = M[2];

		/*these need to be what current positions
		 due to translation and rotation are*/
		for (i = 0; i < s1; i++) {
			mp = ind_com[k].plist[i];

			/*measure each interface to z plane*/
			for (j = 0; j < bases[mp].ninterface; j++) {
				dx = bases[mp].x[j] - x0;
				dy = bases[mp].y[j] - y0;
				dz = bases[mp].z[j] - z0;
				/*only need z component*/
				//rotate(v, M, v2);
				dxrot = row[0] * dx + row[1] * dy + row[2] * dz;
				currx = x0 + traj[k][0] + dxrot;

				xchg = currx - xboxl / 2.0;
				if (xchg > 0) {
					flag = 1;
					if (-xchg < xtot)
						xtot = -xchg;

				} else if (xchg < -xboxl) {
					flag = 1;
					if (-(xchg + xboxl) > xtot)
						xtot = -(xchg + xboxl);
				}
			}
		}
		if (flag == 1) {

			/*Put back inside the box*/
			traj[k][0] += 2.0 * xtot;

		}
		//if(abs(xtot-xtot0)>tol)cout <<"Different calculations of x beyond box! "<<xtot0<<" final: "<<xtot<<endl;

	}
	if (flagy == 1) {
		flag = 0;
		ytot = 0;
		x0 = ind_com[k].xcom;
		y0 = ind_com[k].ycom;
		z0 = ind_com[k].zcom;
		row[0] = M[3];
		row[1] = M[4];
		row[2] = M[5];

		/*these need to be what current positions
		 due to translation and rotation are*/
		for (i = 0; i < s1; i++) {
			mp = ind_com[k].plist[i];

			/*measure each interface to y plane*/
			for (j = 0; j < bases[mp].ninterface; j++) {
				dx = bases[mp].x[j] - x0;
				dy = bases[mp].y[j] - y0;
				dz = bases[mp].z[j] - z0;
				/*only need y component*/
				//rotate(v, M, v2);
				dyrot = row[0] * dx + row[1] * dy + row[2] * dz;
				curry = y0 + traj[k][1] + dyrot;

				ychg = curry - yboxl / 2.0;
				if (ychg > 0) {
					flag = 1;
					if (-ychg < ytot)
						ytot = -ychg;

				} else if (ychg < -yboxl) {
					flag = 1;
					if (-(ychg + yboxl) > ytot)
						ytot = -(ychg + yboxl);
				}
			}
		}
		if (flag == 1) {

			/*Put back inside the box*/
			traj[k][1] += 2.0 * ytot;

		}
		//if(abs(ytot-ytot0)>tol)cout <<"Different calculations of y beyond boy! "<<ytot0<<" final: "<<ytot<<endl;
	}
	if (flagz == 1) {
		flag = 0;
		ztot = 0;
		x0 = ind_com[k].xcom;
		y0 = ind_com[k].ycom;
		z0 = ind_com[k].zcom;
		row[0] = M[6];
		row[1] = M[7];
		row[2] = M[8];

		/*these need to be what current positions
		 due to translation and rotation are*/
		for (i = 0; i < s1; i++) {
			mp = ind_com[k].plist[i];

			/*measure each interface to z plane*/
			for (j = 0; j < bases[mp].ninterface; j++) {
				dx = bases[mp].x[j] - x0;
				dy = bases[mp].y[j] - y0;
				dz = bases[mp].z[j] - z0;
				/*only need z component*/
				//rotate(v, M, v2);
				dzrot = row[0] * dx + row[1] * dy + row[2] * dz;
				currz = z0 + traj[k][2] + dzrot;

				zchg = currz - zboxl / 2.0;
				if (zchg > 0) {
					flag = 1;
					if (-zchg < ztot)
						ztot = -zchg;

				} else if (zchg < -zboxl) {
					flag = 1;
					if (-(zchg + zboxl) > ztot)
						ztot = -(zchg + zboxl);
				}
			}
		}
		if (flag == 1) {

			/*Put back inside the box*/
			traj[k][2] += 2.0 * ztot;

		}
		//if(abs(ztot-ztot0)>tol)cout <<"Different calculations of z bezond boz! "<<ztot0<<" final: "<<ztot<<endl;
	}

}
