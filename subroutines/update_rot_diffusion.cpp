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

void update_rot_diffusion(int c1, Complex *ind_com, Fullmol *bases, double prerot) {
	int i, p1, j;
	int size = ind_com[c1].mysize;
	double a = ind_com[c1].radR;
	if (a == 0)
		a = 1;
	double a3 = a * a * a;
	ind_com[c1].Drx = prerot / a3;
	ind_com[c1].Dry = ind_com[c1].Drx;
	ind_com[c1].Drz = ind_com[c1].Drx;

}
