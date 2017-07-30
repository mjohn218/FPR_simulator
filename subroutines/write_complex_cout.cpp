#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void write_complex_cout(int Nc, Complex *ind_com, Fullmol *bases) {
	int i, j;
	int size;
	int p1;
	for (i = 0; i < Nc; i++) {
		size = ind_com[i].mysize;
		cout << "complex: " << i << " size: " << size << " pros: " << '\t';
		for (j = 0; j < size; j++) {
			p1 = ind_com[i].plist[j];
			cout << p1 << '\t';
		}
		cout << " Dx " << ind_com[i].Dx << " Dy " << ind_com[i].Dy << " Dz " << ind_com[i].Dz << endl;
	}

}
