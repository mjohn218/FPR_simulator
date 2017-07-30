#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void write_complex_size(int Nc, Complex *ind_com, Fullmol *bases, ofstream &outfile, int it) {
	int i, j;
	int size;
	int p1;
	outfile << "iter: " << it << endl;
	for (i = 0; i < Nc; i++) {
		size = ind_com[i].mysize;
		outfile << "complex: " << i << " size: " << size << " Dx " << ind_com[i].Dx << " Dy " << ind_com[i].Dy << " Dz " << ind_com[i].Dz << " pros: " << '\t';
		for (j = 0; j < size; j++) {
			p1 = ind_com[i].plist[j];
			outfile << p1 << '\t';
		}
		outfile << endl;
	}

}
