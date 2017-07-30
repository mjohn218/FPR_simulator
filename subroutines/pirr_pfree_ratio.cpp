#include "GF_calls.h"
#include "Faddeeva.hh"
#include "rand_gsl.h"

double pirr_pfree_ratio(double rcurr, double r0, double deltat, double Dtot, double bindrad, double alpha_irr, double prevpassoc) {

	double pirr = pirrev_value(rcurr, r0, deltat, Dtot, bindrad, alpha_irr);
	// double passoc=survive_irr( r0,  deltat,  Dtot,  bindrad, alpha_irr,  cof);
	double pfree = pfree_value(rcurr, r0, deltat, Dtot, bindrad, alpha_irr, prevpassoc);
	//  cout <<"pirr: " <<pirr<<" pfree: "<<pfree<<endl;
	double ratio = pirr / pfree;
	return ratio;
}
