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

double finfinity(double x, void *p) {
	// This routine is the integrand for kt for the absorbing boundary conditions
	double f;
	f_params &params = *reinterpret_cast<f_params *>(p);
	f = 8 * params.D * exp(-params.D * params.t * x * x) / (M_PI * x * (y0(x * params.a) * y0(x * params.a) + j0(x * params.a) * j0(x * params.a)));
	return f;

}
