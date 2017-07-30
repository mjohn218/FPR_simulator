#include "GF_calls.h"
#include "Faddeeva.hh"
#include "rand_gsl.h"

double survive_irr_abs(double r0, double tcurr, double Dtot, double bindrad) {
	double fDt = 4.0 * Dtot * tcurr;
	double sq_fDt = sqrt(fDt);

	int i, j;
	double sep, dist;

	double r1, term1, term2, e1, ef1, sum;
	double passoc;
	sep = (r0 - bindrad) / sq_fDt;

	term1 = bindrad / r0 * erfc(sep);

	//double psurvive=1.0-term1
	passoc = term1;
	return passoc;
	//  cout <<"s_irr: "<<sirr<<" time: "<<tcurr<<" unscaled: "<<term1-term2<<endl;
}
