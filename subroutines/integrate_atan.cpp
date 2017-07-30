#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "rand_gsl.h"
#include "vol_help.h"

using namespace std;

double integrate_atan(double p, double R, double q1, double q2) {
	/*Integral of arctan[ p/(sqrt(R^2-p^2-y^2)) ] dy
	 x*arctan[p/sqrt(-p^2+R^2-x^2)]
	 + a*arctan(x/sqrt(-p^2+R^2-x^2)) -b*arctan(p*x/R/(sqrt(-p^2+R^2-x^2)))]
	 
	 at q2 we have that R^2-p^2-x^2=0
	 */
	/*Result at q2, where arctan(Rlah/sqrt(R^2-p^2-x^2))=M_PI/2*/

	//double resq2=M_PI/2.0*( q2+p-R);
	double sroot1 = sqrt(R * R - p * p - q1 * q1);
	double resq1 = q1 * atan(p / sroot1) + p * atan(q1 / sroot1) - R * atan(p * q1 / (R * sroot1));
	double sroot2 = sqrt(abs(R * R - p * p - q2 * q2));
	double resq2 = q2 * atan(p / sroot2) + p * atan(q2 / sroot2) - R * atan(p * q2 / (R * sroot2));
	double Voverlap = resq2 - resq1;
	return Voverlap;
}
