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

double ktintegrator(gsl_function F, f_params params, gsl_integration_workspace * w, double r0, double bindrad, double Dtot, double kr, double deltat, char* funcID, double (*f)(double, void*)) {

	const double xlow = 0, epsabs = 1e-8, epsrel = 1e-8;
	double publicationcrit = 1e-6;
	double result, error, umax, du, result2, umax2;
	int status, it, ctr = 1;

	status = gsl_integration_qagiu(&F, xlow, epsabs, epsrel, 1e6, w, &result, &error);
	it = 1.0;
	while (status != GSL_SUCCESS && epsrel * (it + 1.0) < publicationcrit) {
		status = gsl_integration_qagiu(&F, xlow, epsabs * (it + 1.0), epsrel * (it + 1.0), 1e6, w, &result, &error);
		it += 1.0;
	}
	umax = sqrt(-log(DBL_MIN) / (Dtot * deltat));
	gsl_integration_qags(&F, xlow, umax, epsabs, epsrel, 1e6, w, &result2, &error);

	if (status != GSL_SUCCESS || result2 > result) {
		for (double umax = DBL_MIN; umax <= DBL_MAX; umax *= 10) {
			ctr += 1;
		}

		gsl_vector * results = gsl_vector_alloc(ctr);
		gsl_vector * stats = gsl_vector_alloc(ctr);
		gsl_vector_set_zero(results);
		gsl_vector_set_zero(stats);
		ctr = 0;
		for (double umax = DBL_MIN; umax <= DBL_MAX; umax *= 10) {
			status = gsl_integration_qags(&F, xlow, umax, epsabs, epsrel, 1e6, w, &result, &error);
			if (isnormal(result)) {
				gsl_vector_set(results, ctr, result);
				gsl_vector_set(stats, ctr, umax);
			}
			ctr += 1;

			if (gsl_vector_max(results) > 0) {
				if (!isnormal(result) || result == 0) {
					break;
				}
			}
		}
		//refinement
		result = gsl_vector_max(results);
		umax2 = gsl_vector_get(stats, gsl_vector_max_index(results) - 1);
		du = umax2 / 10.0;
		it = 1.0;
		while (umax2 + du * it < 100.0 * umax2) {
			status = gsl_integration_qags(&F, xlow, umax2 + du * it, epsabs, epsrel, 1e6, w, &result2, &error);
			it += 1.0;
			if (status == GSL_SUCCESS && result2 > result) {
				result = result2;
			}
		}
	}

	return result;
}
