#include "reactions.h"
#include "vol_help.h"

void get_volumes(int nbins, Fullmol *bases, double **Vofr, Parms plist, double delr, int N) {
	//check distance to each wall when there are no PBC
	double asum;
	int j;
	double box_x = plist.xboxl;
	double box_y = plist.yboxl;
	double box_z = plist.zboxl;

	double dx;
	double dy;
	double dz;
	int flagx = 0;
	int flagy = 0;
	int flagz = 0;
	double Vnow;
	double Vcap;
	double Vprev = 0; //volume at r=0;
	double r1, r1sq;
	double dV;
	double xcom, ycom, zcom;
	int i;
	double Vfull;
	for (i = 0; i < N; i++) {
		Vofr[i][0] = 0;
		xcom = bases[i].xcom;
		ycom = bases[i].ycom;
		zcom = bases[i].zcom;
		dx = box_x / 2.0 - xcom;
		dy = box_y / 2.0 - ycom;
		dz = box_z / 2.0 - zcom;
		Vprev = 0;
		for (j = 1; j < nbins + 1; j++) {
			//calculate the volume around the particle at each radius
			r1 = j * delr;
			r1sq = r1 * r1;
			Vfull = 4.0 / 3.0 * M_PI * r1 * r1 * r1;
			asum = 0;
			flagx = 0;
			flagy = 0;
			flagz = 0;
			if (dx < r1)
				flagx = 1;
			else if (box_x / 2.0 + xcom < r1) {
				flagx = 1;
				dx = box_x / 2.0 + xcom;
			}
			if (dy < r1)
				flagy = 1;
			else if (box_y / 2.0 + ycom < r1) {
				flagy = 1;
				dy = box_y / 2.0 + ycom;
			}
			if (dz < r1)
				flagz = 1;
			else if (box_z / 2.0 + zcom < r1) {
				flagz = 1;
				dz = box_z / 2.0 + zcom;
			}
			asum = flagx + flagy + flagz;
			if (asum == 0)
				Vnow = Vfull;
			else if (asum == 1) {
				//only need to calculate one spherical cap
				Vcap = Vcap_one_dim(flagx, flagy, dx, dy, dz, r1);
				Vnow = Vfull - Vcap;
			} else if (asum == 2) {
				//there are two caps, also check for overlap
				Vcap = Vcap_two_dim(flagx, flagy, dx, dy, dz, r1);
				Vnow = Vfull - Vcap;
			} else {
				/*all three dimensions are overlapping*/
				Vcap = Vcap_three_dim(flagx, flagy, flagz, dx, dy, dz, r1);
				Vnow = Vfull - Vcap;
			}

			//now get volume fraction for the g(r) at this radius
			dV = Vnow - Vprev;
			Vofr[i][j - 1] = dV;
			Vprev = Vnow;
		}
	}
}
