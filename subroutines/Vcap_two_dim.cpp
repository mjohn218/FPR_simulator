#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "rand_gsl.h"
#include "vol_help.h"

using namespace std;

double Vcap_two_dim(int flagx, int flagy, double dx, double dy, double dz, double R) {
	double h;
	double Vcap1, Vcap2, Voverlap;
	if (flagx != 0) {
		h = R - dx;
		Vcap1 = sphere_cap(R, h);
		if (flagy != 0) {
			//x and y
			h = R - dy;
			Vcap2 = sphere_cap(R, h);
			//check for overlap
			Voverlap = vcap_overlap(dx, dy, R);
		} else {
			//x and z
			h = R - dz;
			Vcap2 = sphere_cap(R, h);
			Voverlap = vcap_overlap(dx, dz, R);
		}
	} else {
		//y and z are the caps

		h = R - dy;
		Vcap1 = sphere_cap(R, h);
		h = R - dz;
		Vcap2 = sphere_cap(R, h);
		Voverlap = vcap_overlap(dy, dz, R);
	}

	double Vcap = Vcap1 + Vcap2 - Voverlap;
	return Vcap;
}
