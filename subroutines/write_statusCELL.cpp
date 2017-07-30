#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void write_statusCELL(ofstream &outfile, ofstream &outfile2, Protein *wholep, Fullmol *bases, Parms &plist, int *Ncopy, int iter, Complex *ind_com, double deltat) {
	int i, j, p;
	int ig;
	int f;
	int t = 0;
	int Ntotalmol = plist.Ntotalmol;
	outfile << "total complexes: " << plist.ntotalcomplex << " iter: " << iter << endl;
	for (i = 0; i < Ntotalmol; i++) {
		outfile << "protein: " << i << " mycomplex: " << bases[i].mycomplex << endl;
		outfile << "Nfree: " << bases[i].nfree << " Freelist: ";
		for (f = 0; f < bases[i].nfree; f++)
			outfile << bases[i].freelist[f] << '\t';
		outfile << endl;
		outfile << "Nbound: " << bases[i].nbnd << " Boundlist: ";
		for (f = 0; f < bases[i].nbnd; f++)
			outfile << bases[i].bndlist[f] << '\t';
		outfile << endl;
		outfile << "Niface: " << bases[i].ninterface << " Istatus: ";
		for (j = 0; j < bases[i].ninterface; j++) {
			outfile << bases[i].istatus[j] << '\t';
		}
		outfile << endl;
		outfile << "Partners: ";
		for (j = 0; j < bases[i].ninterface; j++)
			outfile << bases[i].partner[j] << '\t';
		outfile << endl;

	}
	/*From status of the proteins in the complex, should be able to reconstruct the
	 nature of each complex, its size, its prolist, and its diffusion and COM */

	int *nlistspecies = new int[8];
	int p1, p2, p3;
	for (i = 0; i < 8; i++) {
		nlistspecies[i] = 0;
	}
	int nmolecinside;
	for (i = 0; i < plist.ntotalcomplex; i++) {

		nmolecinside = ind_com[i].mysize;

		if (nmolecinside == 1) {
			p1 = ind_com[i].plist[nmolecinside - 1];
			if (p1 < Ncopy[0]) {
				nlistspecies[0] += 1; //PIP2
			} else if (p1 < (Ncopy[0] + Ncopy[1])) {
				nlistspecies[1] += 1; //AP2
			} else {
				nlistspecies[2] += 1; //EPN1
			}
		} else if (nmolecinside == 2) {
			p1 = ind_com[i].plist[nmolecinside - 2];
			p2 = ind_com[i].plist[nmolecinside - 1];

			if (p1 < Ncopy[0] || p2 < Ncopy[0]) {

				if (p1 >= Ncopy[0] + Ncopy[1] || p2 >= Ncopy[0] + Ncopy[1]) {
					nlistspecies[5] += 1; //EPN1-PIP2
				} else {
					nlistspecies[4] += 1; //AP2-PIP2
				}

			} else if (p1 >= Ncopy[0] || p2 >= Ncopy[0]) {
				nlistspecies[3] += 1; //AP2-EPS
			}

		} else if (nmolecinside == 3) {
			nlistspecies[6] += 1;

		} else if (nmolecinside == 4) {
			nlistspecies[7] += 1;
		}

	}

	outfile2 << iter * deltat * 1e-6 << "\t" << nlistspecies[0] << "\t" << nlistspecies[1] << "\t" << nlistspecies[2] << "\t" << nlistspecies[3] << "\t" << nlistspecies[4] << "\t" << nlistspecies[5] << "\t" << nlistspecies[6] << "\t" << nlistspecies[7] << endl;

}
