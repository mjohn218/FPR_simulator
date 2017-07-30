#include "Vrnumer.h"
#include "numeric_GF.h"


double lookup_psample(double rnum, int pind, double delp, double *table) {
	double rlow = table[pind];
	double rhigh = table[pind + 1];
	double dp = rnum - pind * delp; //fraction along the interval
	double rmed = rlow + (rhigh - rlow) * dp / delp;
	return rmed;

}
