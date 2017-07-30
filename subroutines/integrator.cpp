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

double integrator(gsl_function F, f_params params, gsl_integration_workspace * w, double r0, double bindrad, double Dtot, double kr, double deltat, char* funcID, double (*f)(double, void*)) {

	int it, status;
	double result, error, umax;
	const double xlow = 0, epsabs = 1e-7, epsrel = 1e-7;
	double publicationcrit = 1e-6;

	status = gsl_integration_qagiu(&F, xlow, epsabs, epsrel, 1e6, w, &result, &error);

	it = 1.0;
	while (status != GSL_SUCCESS && epsrel * (it + 1.0) < publicationcrit) {
		status = gsl_integration_qagiu(&F, xlow, epsabs * (it + 1.0), epsrel * (it + 1.0), 1e6, w, &result, &error);
		it += 1.0;
	}

	if (status != GSL_SUCCESS) {
		umax = sqrt(-log(DBL_MIN) / (Dtot * deltat));
		while (abs((*f)(umax, &params)) > 1e-12) {
			umax = umax * 1.1;
		}
		cout << "______________" << endl;
		cout << funcID << " qagiu failed with status: " << status << endl;
		cout << "No solution found with rel/abs error smaller than: " << publicationcrit << endl;
		cout << "   @time: " << deltat << endl;
		cout << "   ka: " << kr << endl;
		cout << "   D: " << Dtot << endl;
		cout << "   r0: " << r0 << endl;
		cout << "Truncation will be performed on the semi-infinite domain..." << endl;

		it = 1.0;
		while (status != GSL_SUCCESS) {
			status = gsl_integration_qags(&F, xlow, umax, 0.0, publicationcrit, 1e6, w, &result, &error);
			umax = umax * 0.99;
		}

		cout << "New integration upper bound " << umax << " found after " << it << " iterations." << endl;
		cout << "______________" << endl;
	}

	return result;
}
