#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "rand_gsl.h"
#include "vol_help.h"

using namespace std;

double sphere_cap(double R, double h) {
	/*R is the sphere radius, H is the height of the cap*/
	double Vcap = 1.0 / 3.0 * M_PI * h * h * (3 * R - h);
	return Vcap;
}
