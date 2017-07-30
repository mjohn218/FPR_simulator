#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "reactions.h"
#include "vol_help.h"
#include "rand_gsl.h"
#include "Faddeeva.hh"
#include "utility_calls.h"
#include "vector_rot_calls.h"

void update_trans_diffusion(int c1, Complex *ind_com, Fullmol *bases, double pretrans) {
	int i, p1, j;
	int size = ind_com[c1].mysize;
	double a = ind_com[c1].radR;
	if (a == 0)
		a = 1;
	ind_com[c1].Dx = pretrans / a;
	ind_com[c1].Dy = ind_com[c1].Dx;
	ind_com[c1].Dz = ind_com[c1].Dx;
	/*For Dz, need to check if it is bound to membrane*/
	for (i = 0; i < size; i++) {
		p1 = ind_com[c1].plist[i];
		if (bases[p1].Dz == 0) {
			ind_com[c1].Dz = 0;
			i = size;
		}
	}

}
