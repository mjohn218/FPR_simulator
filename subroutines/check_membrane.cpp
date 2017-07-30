#include "reactions.h"

void check_membrane(int *Ncopy, Fullmol *bases) {
	int i;
	int start = Ncopy[0] + Ncopy[1]; //start with pip2
	int Ntot = Ncopy[2] + Ncopy[3];
	double zlow = -50;
	for (i = start; i < Ntot; i++) {
		if (bases[i].zcom > zlow) {
			cout << "Membrane molecule: " << i << " zvalue: " << bases[i].zcom << endl;
		}
	}

}
