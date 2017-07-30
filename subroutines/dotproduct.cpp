#include "vector_rot_calls.h"

void dotproduct(double *v1, double *v2, double &theta) {
	double dp = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
	double v1m = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
	double v2m = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2];
	double cthet = dp / (sqrt(v1m * v2m));
	double df = cthet - 1;
	double tol = 1E-12;
	if (abs(df) < tol)
		cthet = 1;
	theta = acos(cthet);
	if (dp == 0)
		theta = 0;

}
