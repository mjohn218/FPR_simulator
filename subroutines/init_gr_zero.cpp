#include "reactions.h"
#include "vol_help.h"

void init_gr_zero(double *gr, int nbins) {
	int i;
	for (i = 0; i < nbins; i++)
		gr[i] = 0.0;
}
