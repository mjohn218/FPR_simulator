#include "GF_calls.h"
#include "Faddeeva.hh"
#include "rand_gsl.h"

double psample_final(double rmin, double r0, double bindrad, double passoc, double Rlong, double R_table_end, double rnum, double alpha, double kact, double Dtot, double deltat) {
	/*To sample from irreversible distribution, kact is finite, and alpha=(1+kact/kdiff)*sqrt(D)/sigma
	 To sample from reflecting distribution, kact=0 and alpha=sqrt(D)/sigma
	 */

	double r1, term1, term2, term3, term4;
	double cumsum;

	int j, jprev;
	jprev = 0;

	//or just skip all the probs less than pstart, so you don't need to know what p_survive is 
	int i = 0;
	double prob = 0;
	int i1 = i;
	//cout <<"first index: "<<i1<<" for reflecting should be zero "<<kact<<endl; 
	double psub;
	psub = peval_cumulative(bindrad, r0, deltat, Dtot, bindrad, alpha, kact);
	double plong = peval_cumulative(Rlong, r0, deltat, Dtot, bindrad, alpha, kact);

	double pirr_cume, prevsum;
	double ps_numer = plong - psub;

	double psurvive = 1 - passoc;
	//cout <<"calculated psurvive from fit: "<<psurvive<<" integral over numer pirr, to Rlong: "<<ps_numer<<" Rlong: "<<Rlong<<" ratio: "<<psurvive/ps_numer<<endl;
	double eps = 1E-9;
	double ps_ratio = (psurvive + eps) / ps_numer;

	//cout <<"value at sigma to be subtracted: "<<psub<<" Value of cumul at limit r: "<<limit.half<<" is; "<<phalf<<endl;
	/*Since pirr may not integrate to psurvive, rescale it to. */

	prob = rnum; //rnum will be between delp*MAXP-passoc and psurvive.

	j = 0; //jprev;
	//cout <<"time start: "<<tmin*1<<endl;
	cumsum = 0;
	double cumsum1;
	double jmax = 1E7;
	while (prob > cumsum && j < jmax) {
		r1 = rmin * j + R_table_end;
		/*need to use other fit parms, and subtract off sum at previous r1.*/

		cumsum1 = peval_cumulative(r1, r0, deltat, Dtot, bindrad, alpha, kact);
		cumsum = (cumsum1 - psub) * ps_ratio;

		/*this is the cumulative probability for a given r value, compare this sum to a random number to 
		 generate the separation
		 we want the probability to be evenly spaced
		 It will sum up to the survival probability, since it survived, the random number will be between p_assoc=1-survival and 1. 
		 */
		j++;
	}
	if (j >= jmax)
		cout << "broke out of final sample loop, pirr_2sfi. final r:  " << r1 << " prob: " << prob << " final cumulative sum: " << cumsum << endl;

	cout << "Psample_final " << r0 << " r1-bindrad: " << r1 - bindrad << '\t' << " Final prob1? " << prob << " Psurvive: " << psurvive << " Randonum Prob to match: " << rnum << " Rtable_end: " << R_table_end << " cumsum: " << cumsum << " j: " << j << " pcume: " << cumsum1 << " psub: " << psub << " ps_ratio " << ps_ratio << endl;
	return r1;
}
