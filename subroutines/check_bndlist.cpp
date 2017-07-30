#include "reactions.h"

void check_bndlist(int Ncomplex, Complex *ind_com, Fullmol *bases) {
	int i, j, k;
	int size, n;
	int p1;
	cout << "Checking bound list " << endl;
	for (i = 0; i < Ncomplex; i++) {
		size = ind_com[i].mysize;
		cout << "complex: " << i << " size: " << size << endl;
		for (j = 0; j < size; j++) {
			p1 = ind_com[i].plist[j];
			n = bases[p1].nbnd;
			cout << " protein: " << p1 << " nbnd: " << n << '\t';
			for (k = 0; k < n; k++) {
				cout << bases[p1].bndlist[k] << '\t';
			}
			cout << endl;
		}
	}

}
