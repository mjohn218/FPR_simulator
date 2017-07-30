#include "GF_calls.h"
#include "Faddeeva.hh"
#include "rand_gsl.h"

void pirrev_dist(int rbins, double r0, double tcurr, double Dtot, double bindrad, double *pirrev, double alpha, double delr) {
	double fDt = 4.0 * Dtot * tcurr;
	double sq_fDt = sqrt(fDt);

	double f1 = 1.0 / (sqrt(4.0 * M_PI * tcurr));
	double f2 = 1.0 / (4.0 * M_PI * r0 * sqrt(Dtot));
	int i, j;
	double sep, dist;
	double sqrt_t = sqrt(tcurr);
	double a2 = alpha * alpha;
	double r1, term1, term2, e1, ef1, sum;

	for (i = 0; i < rbins; i++) {
		r1 = delr * (i + 0.5) + bindrad;
		sep = r1 + r0 - 2.0 * bindrad;
		dist = r1 - r0;

		term1 = f1 * (exp(-dist * dist / fDt) + exp(-sep * sep / fDt));

		e1 = 2.0 * sep / sq_fDt * sqrt_t * alpha + a2 * tcurr;
		ef1 = sep / sq_fDt + alpha * sqrt_t;
		term2 = alpha * exp(e1) * erfc(ef1);

		sum = term1 - term2;
		//    cout <<"irr 1: "<<term1<<" irr 2: "<<term2<<" sum: "<<sum<<" scaled: "<<sum*f2/r1<<endl;
		sum *= f2 / r1;
		pirrev[i] = sum;

	}

}
