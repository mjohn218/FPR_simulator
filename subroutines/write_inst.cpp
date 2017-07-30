#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void write_inst(ofstream &outfile, Parms &plist, int **inst_pro, Protein *wholep, int it) {
	int i, j;
	int nint;
	outfile << it << '\t';
	for (i = 0; i < plist.Nprotypes; i++) {
		nint = wholep[i].ninterface;
		outfile << "protein: " << i << '\t';
		for (j = 0; j < nint + 1; j++)
			outfile << inst_pro[i][j] << '\t';

	}
	outfile << endl;

}
