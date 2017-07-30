/*
 2D related 
 */

#include <stdio.h>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_vector.h>
#include "2Drelated.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <cfloat>

using namespace std;

int jac(double t, const double y[], double *dfdy, double dfdt[], void *params) {
	double mu = *(double *) params;
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 1, 1);
	gsl_matrix * m = &dfdy_mat.matrix;
	gsl_matrix_set(m, 0, 0, -2 * mu * y[0]);
	dfdt[0] = 10.0;

	return GSL_SUCCESS;
}
