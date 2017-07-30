#include "GF_calls.h"
#include "Faddeeva.hh"
#include "rand_gsl.h"

void psample_table(double rmin, double r0, double tcurr, double Dtot, double bindrad, double alpha, double kact, double *table, int MAXP, double delp, double Rlong) {
	/*To sample from irreversible distribution, kact is finite, and alpha=(1+kact/kdiff)*sqrt(D)/sigma
	 To sample from reflecting distribution, kact=0 and alpha=sqrt(D)/sigma
	 */

	double fDt = 4.0 * Dtot * tcurr;
	double sqfDt = sqrt(fDt);
	double kdiff = 4.0 * M_PI * Dtot * bindrad;
	double cof = kact / (kact + kdiff);
	double spi = sqrt(M_PI);
	double c1 = 1.0 / (2.0 * r0 * sqrt(M_PI * Dtot * tcurr));
	double sep, dist, sep2, d2;
	double sqrt_t = sqrt(tcurr);
	double a, b;
	double r1, term1, term2, term3, term4;
	double cumsum;
	double c3 = alpha / (r0 * sqrt(Dtot));
	double off = r0 - 2.0 * bindrad;
	double passoc = survive_irr(r0, tcurr, Dtot, bindrad, alpha, cof);

	int j, jprev;
	jprev = 0;
	b = alpha * sqrt_t;
	//or just skip all the probs less than pstart, so you don't need to know what p_survive is 
	int i = 0;
	double prob = 0;
	table[0] = bindrad;
	double psub = peval_cumulative(bindrad, r0, tcurr, Dtot, bindrad, alpha, kact);
	double plong = peval_cumulative(Rlong, r0, tcurr, Dtot, bindrad, alpha, kact);
	int i1 = i;
	//cout <<"first index: "<<i1<<" for reflecting should be zero "<<kact<<endl; 
	double psurvive = 1.0 - passoc;
	double ps_numer = plong - psub;
	double eps = 1E-9;
	cout << "calculated psurvive from fit: " << psurvive << " integral over numer pirr, to Rlong: " << ps_numer << " Rlong: " << Rlong << " ratio: " << psurvive / ps_numer << endl;
	double ps_ratio = (psurvive + eps) / ps_numer;

	double cumsum1;
	// cout <<"value at sigma to be subtracted: "<<psub<<endl;
	for (i = i1; i < MAXP - 1; i++) {
		//now find which time this corresponds to
		prob = delp * (i + 1); //at final value of i, prob will be equal=1
		j = jprev;
		//cout <<"time start: "<<tmin*1<<endl;
		cumsum = 0;
		while (prob > cumsum) {
			r1 = rmin * j + bindrad;
			cumsum1 = peval_cumulative(r1, r0, tcurr, Dtot, bindrad, alpha, kact);
			// sep=r1+off;//r0-2.0*bindrad;
//       sep2=sep*sep;
//       dist=r1-r0;
//       d2=dist*dist;
//       a=sep/sqfDt;
//       term1=-0.5*erf(-dist/sqfDt)-0.5*fDt*c1*exp(-d2/fDt);
//       term2=-0.5*(r0-2.0*bindrad)/r0*erf(sep/sqfDt)-0.5*fDt*c1*exp(-sep2/fDt);
//       term3=c3*sqrt(4.0*Dtot)/(2.0*alpha*spi)*( -sqfDt*exp(-a*a)-spi*off*erf(a) );
//       term4=c3*-sqfDt/(4.0*b*b)*( sqfDt*erf(a)+(sqfDt-2.0*b*r1)*exp(2.0*a*b+b*b)*erfc(a+b) );
//       cumsum=term1+term2-term3-term4-psub;
			cumsum = (cumsum1 - psub) * ps_ratio;
			//cumsum*=ps_ratio;
			/*this is the cumulative probability for a given r value, compare this sum to a random number to 
			 generate the separation
			 we want the probability to be evenly spaced
			 It will sum up to the survival probability, since it survived, the random number will be between p_assoc=1-survival and 1. 
			 */
			j++;
		}
		if (j - jprev > 10000) {
			rmin *= 10;
			j /= 10;
			j -= 10;

		}
		//    cout <<"iterations: "<<j-jprev<<" time: "<<time<<" tmin: "<<tmin<<endl;
		jprev = int(round(j - 1));
		//now add this to the table, so it can be looked up.
		//except it is evenly space in time, not in prob, so interpolate?
		table[i + 1] = r1;
		//cout <<r1<<'\t'<<prob-probstart<<endl;
	}
	cout << r0 << '\t' << r1 << '\t' << prob << " psurvive: " << psurvive << " rmin at end:" << rmin << endl;

}
