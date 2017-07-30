#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "rand_gsl.h"
#include "vol_help.h"

using namespace std;

double Vcap_one_dim(int flagx, int flagy, double dx, double dy, double dz, double R) {
	double h;
	if (flagx != 0)
		h = R - dx;
	else if (flagy != 0)
		h = R - dy;
	else
		h = R - dz;

	double Vcap = sphere_cap(R, h);
	return Vcap;
}
