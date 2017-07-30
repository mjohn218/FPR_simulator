#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void write_timestat(ofstream &outfile, Protein *wholep, Fullmol *bases, Parms &plist, int *Ncopy, int iter, Complex *ind_com, double deltat, int Nprotypes) {
	int i, j, p, k;
	int ig;
	int f;
	int t = 0;
	int Ntotalmol = plist.Ntotalmol;
	int MaxComplexSize = 1000;
	int nmolecinside, typeofmol;
	int maxcompt = 0;

	int lims[Nprotypes];

	lims[0] = Ncopy[0];

	for (k = 1; k < Nprotypes; k++) {
		lims[k] = lims[k - 1] + Ncopy[k];
	}

	int sumdata[Nprotypes][MaxComplexSize];
	for (i = 0; i < Nprotypes; i++) {
		for (j = 0; j < MaxComplexSize; j++) {
			sumdata[i][j] = 0;
		}
	}

	for (i = 0; i < plist.ntotalcomplex; i++) {
		nmolecinside = ind_com[i].mysize;

		if (nmolecinside > maxcompt) {
			maxcompt = nmolecinside;
		}

		for (j = 0; j < nmolecinside; j++) {
			for (k = 0; k < Nprotypes; k++) {
//			  if(bases[ind_com[i].plist[j]].Dx==wholep[k].Dx){
//				  typeofmol = k;
//				  break;
//			  }
				if (ind_com[i].plist[j] < lims[k]) {
					typeofmol = k;
					break;
				}

			}
			sumdata[typeofmol][nmolecinside - 1] += 1;
		}

	}

	outfile << iter * deltat * 1e-6 << "\t";
	for (i = 0; i < maxcompt; i++) {
		for (k = 0; k < Nprotypes; k++) {
			outfile << sumdata[k][i] << "\t";
		}
	}
	outfile << endl;
}
