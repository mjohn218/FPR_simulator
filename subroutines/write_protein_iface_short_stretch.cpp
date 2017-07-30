#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>

void write_protein_iface_short_stretch(ofstream &outfile, Parms &plist, Fullmol *bases, int *Ncopy, int it, Protein *wholep, string *names, double stretch) {
	int i, j;
	int t = 0;
	outfile << plist.Natomwrite << endl;
	outfile << "iter:" << it << endl;
	char aname;
	int n, nint;
	int ind;
	for (i = 0; i < plist.Nprotypes; i++) {

		nint = wholep[i].nint_write;
		if (i == plist.pclath) {
			for (j = 0; j < Ncopy[i]; j++) {
				n = 0;
				outfile << names[i] << ' ' << bases[t].x[n] << ' ' << bases[t].y[n] << ' ' << bases[t].z[n] << endl;
				for (n = 1; n < 4; n++) {
					outfile << names[i] << ' ' << (bases[t].x[n] - bases[t].xcom) * stretch + bases[t].xcom << ' ' << (bases[t].y[n] - bases[t].ycom) * stretch + bases[t].ycom << ' ' << (bases[t].z[n] - bases[t].zcom) * stretch + bases[t].zcom << endl;
				}
				for (n = 4; n < nint; n++) {
					outfile << names[i] << ' ' << bases[t].x[n] << ' ' << bases[t].y[n] << ' ' << bases[t].z[n] << endl;
				}
				t++;
			}
		} else {
			for (j = 0; j < Ncopy[i]; j++) {

				for (n = 0; n < nint; n++) {
					outfile << names[i] << ' ' << bases[t].x[n] << ' ' << bases[t].y[n] << ' ' << bases[t].z[n] << endl;
				}
				t++;
			}
		}
	}

}
