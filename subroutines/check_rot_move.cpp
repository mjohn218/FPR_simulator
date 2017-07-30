#include "reactions.h"

void check_rot_move(int Nmol, Fullmol *bases, Complex *ind_com) {
	int i, k;
	double x0, y0, z0;
	double dx, dy, dz;
	double r2;
	double tol = 1E-12;
	for (i = 0; i < Nmol; i++) {
		k = bases[i].mycomplex;
		x0 = ind_com[k].xcom;
		y0 = ind_com[k].ycom;
		z0 = ind_com[k].zcom;
		dx = bases[i].xcom - x0;
		dy = bases[i].ycom - y0;
		dz = bases[i].zcom - z0;
		r2 = dx * dx + dy * dy + dz * dz;
		if (r2 > tol) {
			cout << " COM moved from complex, protein: " << i << " complex: " << k << " complex size: " << ind_com[k].mysize << " r2: " << r2 << endl;
		}
		dx = bases[i].x[0] - x0;
		dy = bases[i].y[0] - y0;
		dz = bases[i].z[0] - z0;
		r2 = dx * dx + dy * dy + dz * dz;
		if (r2 > tol) {
			cout << " Interface moved from complex, protein: " << i << " complex: " << k << " complex size: " << ind_com[k].mysize << " r2: " << r2 << endl;
		}

	}

}
