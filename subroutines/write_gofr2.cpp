#include "reactions.h"
#include "vol_help.h"

void write_gofr2(int n1, int n2, Parms &plist, int nbins, double **gr, double delr, ofstream &outfile) {
	int i;
	double r;
	for (i = 0; i < nbins; i++) {
		r = (i + 0.5) * delr;
		outfile << r << ' ' << gr[n1 * plist.Nprotypes + n2][i] << endl;
	}

}
