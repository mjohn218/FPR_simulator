#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "rand_gsl.h"
#include "vol_help.h"

using namespace std;

double integrate_atanx2(double p, double R, double q1, double q2) {
	/*Integral of y^2*arctan[ p/(sqrt(R^2-p^2-y^2)) ] dy
	 1/6*[ -px sqrt(-p^2+R^2-x^2) -p*(p^2-3*R^2)arctan[x/Bsqrt(-p^2+R^2-x^2)]
	 + 2*x^3 arctan(p/sqrt(-p^2+R^2-x^2)) -2*R^3arctan(p*x/R/(sqrt(-p^2+R^2-x^2)))]
	 
	 at q2 we have that R^2-p^2-x^2=0
	 */
	/*Result at q2, where arctan(Rlah/sqrt(R^2-p^2-x^2))=M_PI/2*/
	double t2 = p * (p * p - 3.0 * R * R);
	double t4 = 2.0 * R * R * R;
	//  double resq2=1.0/6.0*M_PI/2.0*( -t2+2*q2*q2*q2-t4 );
	double sroot1 = sqrt(R * R - p * p - q1 * q1);
	double sroot2 = sqrt(abs(R * R - p * p - q2 * q2));
	double resq1 = 1.0 / 6.0 * (-p * q1 * sroot1 - t2 * atan(q1 / sroot1) + 2.0 * q1 * q1 * q1 * atan(p / sroot1) - t4 * atan(p * q1 / (R * sroot1)));
	double resq2 = 1.0 / 6.0 * (-p * q2 * sroot2 - t2 * atan(q2 / sroot2) + 2.0 * q2 * q2 * q2 * atan(p / sroot2) - t4 * atan(p * q2 / (R * sroot2)));
	double Voverlap = resq2 - resq1;
	return Voverlap;
}