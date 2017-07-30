#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void write_protein(ofstream &outfile, Parms &plist, Fullmol *bases, int *Ncopy, int it) {
	int i, j;
	int t = 0;
	outfile << plist.Ntotalmol << endl;
	outfile << "iter" << it << endl;
	char aname;
	char types[plist.Nprotypes];
	types[0] = 'H';
	types[1] = 'C';
	types[2] = 'O';
	for (i = 3; i < plist.Nprotypes; i++)
		types[i] = 'S';

	for (i = 0; i < plist.Nprotypes; i++) {
		aname = types[i];
		for (j = 0; j < Ncopy[i]; j++) {
			//outfile<<bases[t].xcom<<' '<<bases[t].ycom<<' '<<bases[t].zcom<<' '<<bases[t].mycomplex<<' '<<bases[t].mycomplexsize<<endl;
			outfile << aname << ' ' << bases[t].xcom << ' ' << bases[t].ycom << ' ' << bases[t].zcom << endl;
			t++;
		}
	}

}
