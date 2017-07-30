#include "reactions.h"

void check_movement(int Nmol, Fullmol *bases, double *crd) {
	int i;
	double dx, dy, dz;
	double r2;
	for (i = 0; i < Nmol; i++) {
		dx = bases[i].xcom - crd[i * 3];
		dy = bases[i].ycom - crd[i * 3 + 1];
		dz = bases[i].zcom - crd[i * 3 + 2];
		r2 = dx * dx + dy * dy + dz + dz;
		if (r2 == 0)
			cout << "MOLECULE " << i << " DID NOT MOVE! " << endl;
	}

}
