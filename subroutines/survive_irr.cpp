#include "GF_calls.h"
#include "Faddeeva.hh"
#include "rand_gsl.h"

double survive_irr(double r0, double tcurr, double Dtot, double bindrad, double alpha, double cof) {
	double fDt = 4.0 * Dtot * tcurr;
	double sq_fDt = sqrt(fDt);

	double f1 = cof * bindrad / r0;

	int i, j;
	double sep, dist;
	double sqrt_t = sqrt(tcurr);
	double a2 = alpha * alpha;
	double r1, term1, term2, e1, ef1, sum;
	double onemsirr;
	sep = (r0 - bindrad) / sq_fDt;

	e1 = 2.0 * sep * sqrt_t * alpha + a2 * tcurr;
	ef1 = sep + alpha * sqrt_t;
	term1 = erfc(sep);
	term2 = exp(e1) * erfc(ef1);

	sum = term1 - term2;
	sum *= f1;
	onemsirr = sum;
	return onemsirr; //1-sirr=passoc
	//  cout <<"s_irr: "<<sirr<<" time: "<<tcurr<<" unscaled: "<<term1-term2<<endl;
}
