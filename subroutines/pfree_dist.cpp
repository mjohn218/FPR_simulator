#include "GF_calls.h"
#include "Faddeeva.hh"
#include "rand_gsl.h"

void pfree_dist(int rbins, double r0, double tcurr, double Dtot, double bindrad, double *pfree, double alpha, double delr, double passoc) {
	/*Psurvive for the first time step is 1-passoc and should be exact. After
	 more than one step the distribution will be off so passoc will be off as well*/
	double psurvive = 1.0 - passoc;
	double fDt = 4.0 * Dtot * tcurr;
	double sq_fDt = sqrt(fDt);

	double f1 = 1.0 / (sqrt(4.0 * M_PI * tcurr));
	double f2 = 1.0 / (4.0 * M_PI * r0 * sqrt(Dtot));
	double cof = f1 * f2; //equal to 1/( 8pir_0 sqrt(piDt))
	int i, j;
	double sep, dist;
	double sqrt_t = sqrt(tcurr);
	double a2 = alpha * alpha;
	double r1, term1, term2, e1, ef1, sum;
	double adist;
	/*Normalization for no diffusion inside binding radius!*/
	double c1 = 4.0 * M_PI * cof;
	dist = bindrad - r0;
	adist = bindrad + r0;
	term1 = -0.5 * fDt * exp(-dist * dist / fDt) - 0.5 * sqrt(4.0 * M_PI * Dtot * tcurr) * r0 * erf(-dist / sqrt(4.0 * Dtot * tcurr));
	term2 = 0.5 * fDt * exp(-adist * adist / fDt) + 0.5 * sqrt(4.0 * M_PI * Dtot * tcurr) * r0 * erf(adist / sqrt(4.0 * Dtot * tcurr));
	double pnorm = 1.0 - c1 * (term1 + term2); //this is then the normlization from sigma->inf 

	//  cout <<"Normlization: integratl from sigma to infinity: "<<pnorm<<endl;
	for (i = 0; i < rbins; i++) {
		r1 = delr * (i + 0.5) + bindrad;
		//    sep=r1+r0-2.0*bindrad;
		dist = r1 - r0;
		adist = r1 + r0;
		term1 = exp(-dist * dist / fDt) - exp(-adist * adist / fDt);

		pfree[i] = cof / r1 * term1 * psurvive / pnorm;

	}

}
