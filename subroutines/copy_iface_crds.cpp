#include "reactions.h"

void copy_iface_crds(int Nmol, Fullmol *bases, double *crd, Protein *wholep) {
	int j;
	int p;
	int t = 0;
	for (int i = 0; i < Nmol; i++) {
		p = bases[i].protype;
		for (j = 0; j < wholep[p].ninterface; j++) {
			crd[t * 3] = bases[i].x[j];
			crd[t * 3 + 1] = bases[i].y[j];
			crd[t * 3 + 2] = bases[i].z[j];
			t++;
		}
	}

}
