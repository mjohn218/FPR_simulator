#include "reactions.h"

void check_pro_to_com(int Nmol, Fullmol *bases, Complex *ind_com) {
	int i, j;
	int c1, cs;
	int flag;
	for (i = 0; i < Nmol; i++) {
		c1 = bases[i].mycomplex;
		cs = ind_com[c1].mysize;
		flag = 0;
		for (j = 0; j < cs; j++) {
			if (ind_com[c1].plist[j] == i)
				flag = 1;
		}
		if (flag == 0) {
			cout << "COMPLEX " << c1 << " does not contain protein: " << i << " its size: " << cs << endl;
		}

	}

}
