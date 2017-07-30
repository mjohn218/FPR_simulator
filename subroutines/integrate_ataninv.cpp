#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "rand_gsl.h"
#include "vol_help.h"

using namespace std;

double integrate_ataninv(double p, double R, double q1, double q2) {
	/*Integral of arctan[ (sqrt(R^2-p^2-y^2)/p ] dy
	 -p*arctan[x/sqrt(-p^2+R^2-x^2)]
	 + R*arctan(px/R/sqrt(-p^2+R^2-x^2)) +x*arctan((sqrt(-p^2+R^2-x^2)/p)]
	 

	 */
	double sroot1 = sqrt(R * R - p * p - q1 * q1);
	double sroot2 = sqrt(R * R - p * p - q2 * q2);
	double resq1 = -p * atan(q1 / sroot1) + R * atan(p * q1 / (R * sroot1)) + q1 * atan(sroot1 / p);
	double resq2 = -p * atan(q2 / sroot2) + R * atan(p * q2 / (R * sroot2)) + q2 * atan(sroot2 / p);
	double Voverlap = resq2 - resq1;
	return Voverlap;
}
