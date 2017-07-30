#include "GF_calls.h"
#include "Faddeeva.hh"
#include "rand_gsl.h"

double pirr_pfree_diff(double rcurr, double r0, double deltat, double Dtot, double bindrad, double alpha_irr, double prevpassoc, double delr0) {

	/* version 1, very litte affect.
	 double pirr=pirrev_value(rcurr, r0, deltat, Dtot, bindrad, alpha_irr);
	 double pfree=pfree_value(rcurr, r0, deltat, Dtot, bindrad, alpha_irr, prevpassoc);
	 double diff=pirr-pfree;
	 */

	double pirr = pirrev_value(rcurr, r0, deltat, Dtot, bindrad, alpha_irr);
	double pfree = pfree_value(rcurr, r0, deltat, Dtot, bindrad, alpha_irr, prevpassoc);
	double diff = (pirr - pfree) * 4.0 * M_PI * rcurr * rcurr * delr0;
	return diff;
}
