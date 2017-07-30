#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void write_complex(ofstream &outfile, Parms &plist, Complex *ind_com, double time) {
	int i;
	outfile << plist.ntotalcomplex << endl;
	outfile << "time: " << time << endl;
	for (i = 0; i < plist.ntotalcomplex; i++)
		outfile << "S " << ind_com[i].xcom << ' ' << ind_com[i].ycom << ' ' << ind_com[i].zcom << endl;

}
