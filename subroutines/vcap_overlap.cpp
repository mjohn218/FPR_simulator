#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "rand_gsl.h"
#include "vol_help.h"

using namespace std;

double vcap_overlap(double dx, double dy, double R) {
	double diag = sqrt(dx * dx + dy * dy);
	double p;
	double Vtot = 0;
	double q1, q1sq;
	double q2, q2sq;

	if (diag < R) {
		//we have overlap, correct for it!
		p = dx; //lower bound
		//bounds integrate over y are q1 to q2
		q1 = dy;
		q1sq = q1 * q1;
		q2sq = (R * R - p * p);
		q2 = sqrt(q2sq);
		double A2 = q2sq; //because it is R^2-p^2
		double V1 = M_PI / 2.0 * (R * R * (q2 - q1) - 1.0 / 3.0 * (q2sq * q2 - q1sq * q1));
		double V2 = -p / 2.0 * (A2 * M_PI / 2.0) + p / 2.0 * (q1 * sqrt(A2 - q1sq) + A2 * atan(q1 / sqrt(A2 - q1sq)));

		double V3 = integrate_atanx2(p, R, q1, q2);
		double V4 = -R * R * integrate_atan(p, R, q1, q2);
		Vtot = V1 + V2 + V3 + V4;
	}
	return Vtot;
}
