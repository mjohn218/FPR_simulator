#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "rand_gsl.h"
#include "vol_help.h"

using namespace std;

double vcap_three(double dx, double dy, double dz, double R) {
	double diag = sqrt(dx * dx + dy * dy + dz * dz);
	double p;
	double Vtot = 0;
	double h;
	double q1, q1sq;
	double q2, q2sq;
	if (diag < R) {
		//we have overlap, correct for it!
		p = dx; //lower bound
		h = dz;
		//bounds integrate over y are q1 to q2
		q1 = dy;
		q1sq = q1 * q1;
		q2sq = (R * R - p * p - h * h);
		q2 = sqrt(q2sq);
		double A2 = R * R - p * p;
		double B2 = R * R - h * h;
		double R2 = R * R;
		double V1 = -1.0 / 4.0 * p * (q2 * sqrt(A2 - q2sq) - q1 * sqrt(A2 - q1sq) + A2 * atan(q2 / sqrt(A2 - q2sq)) - A2 * atan(q1 / sqrt(A2 - q1sq)));
		double V7 = -1.0 / 4.0 * h * (q2 * sqrt(B2 - q2sq) - q1 * sqrt(B2 - q1sq) + B2 * atan(q2 / sqrt(B2 - q2sq)) - B2 * atan(q1 / sqrt(B2 - q1sq)));
		double V2 = 0.5 * R2 * integrate_ataninv(h, R, q1, q2);
		double V3 = -0.5 * integrate_atanx2inv(h, R, q1, q2);
		double V4 = -0.5 * R2 * integrate_atan(p, R, q1, q2);
		double V5 = 0.5 * integrate_atanx2(p, R, q1, q2);
		double V6 = h * p * (q2 - q1);
		Vtot = V1 + V2 + V3 + V4 + V5 + V6 + V7;
	}
	return Vtot;
}
