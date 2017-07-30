#include "reactions.h"

void check_com_to_pro(int Nmol, int Ncomplex, Fullmol *bases, Complex *ind_com) {
	int i, j;
	int size;
	int p1;
	int sum = 0;
	for (i = 0; i < Ncomplex; i++) {
		size = ind_com[i].mysize;
		sum += size;
		for (j = 0; j < size; j++) {
			p1 = ind_com[i].plist[j];

			if (bases[p1].mycomplex != i) {
				cout << "PROTEIN " << p1 << " not matched to complex " << i << " . protein complex is: " << bases[p1].mycomplex << endl;
			}
		}
	}
	if (sum != Nmol) {
		cout << "PROTEINS in complex does not sum to total number of proteins! " << endl;
	}

}
