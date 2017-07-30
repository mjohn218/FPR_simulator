#include "GF_calls.h"
#include "Faddeeva.hh"
#include "rand_gsl.h"

void pcorrection(int rbins, double delr, double deltat, double Dtot, double bindrad, double cof, double alpha_irr, double *prassoc, int r0bins, double delr0, double *pcorrect, double *pfreea) {
	double ifree = 0;
	double iexact = 0;
	double rval;
	double passoc;
	int i, n;
	double *pcomp = new double[rbins];
	double *pfree = new double[rbins];
	double r0;
	double pnorm;
	for (n = 0; n < r0bins; n++) {
		r0 = bindrad + (n + 0.5) * delr0;
		passoc = survive_irr(r0, deltat, Dtot, bindrad, alpha_irr, cof);
		pfree_dist(rbins, r0, deltat, Dtot, bindrad, pfree, alpha_irr, delr, passoc);
		pirrev_dist(rbins, r0, deltat, Dtot, bindrad, pcomp, alpha_irr, delr);
		ifree = 0;
		iexact = 0;
		pnorm = 0;
		for (i = 0; i < rbins; i++) {
			rval = bindrad + (i + 0.5) * delr;
			//freefile<<rval<<'\t'<<pfree[j]<<endl;
			ifree += rval * rval * pfree[i] * prassoc[i];
			pnorm += rval * rval * pfree[i];
			iexact += rval * rval * pcomp[i] * prassoc[i];
		} //done looping over all r, given r0

		pcorrect[n] = 4.0 * M_PI * (iexact - ifree) * delr;
		if (isnan(pcorrect[n])) {
			pcorrect[n] = 0;
			cout << "Replaced NAN wit zero " << endl;
		}
		pfreea[n] = 4.0 * M_PI * ifree * delr;
		//cout <<"r0: "<<r0<<" correction "<<pcorrect[n]<<" pfree: "<<pfreea[n]<<endl;
	}

	delete[] pfree;
	delete[] pcomp;
}
