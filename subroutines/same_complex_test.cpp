#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "reactions.h"
#include "vol_help.h"
#include "rand_gsl.h"
#include "Faddeeva.hh"
#include "utility_calls.h"
#include "vector_rot_calls.h"

int same_complex_test(int p1, int p2, Fullmol *bases, int **crossint, int *i_home, double *bindrad, int rxn1, int ci1, int ci2) {

	int flag2 = 0;
	double dx, dy, dz, dist2;
	double multiple = 4.2; //don't let them travel too far to bind when within the same complex
	int iind, iind2;
	if (bases[p1].mycomplex == bases[p2].mycomplex) {

		/*these proteins are already connected, do not associate them if they are
		 too far apart*/
		iind = i_home[crossint[p1][ci1]];
		iind2 = i_home[crossint[p2][ci2]];

		dx = -bases[p2].x[iind2] + bases[p1].x[iind];
		dy = -bases[p2].y[iind2] + bases[p1].y[iind];
		dz = -bases[p2].z[iind2] + bases[p1].z[iind]; //This should be zero now!
		dist2 = dx * dx + dy * dy + dz * dz;

		if (dist2 > (bindrad[rxn1] * bindrad[rxn1] * multiple)) {
			flag2 = 1;
			cout << "BOTH PROTEINS IN SAME COMPLEX AND TOO FAR APART (SQRED): " << p1 << ' ' << p2 << ' ' << dist2 << endl;
		}
	}
	return flag2;
}
