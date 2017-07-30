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

double callfassoc2dMF(double bindrad, double Dtot, double kr, double deltat, double rho) {
	/////////////////////ASSC(t)////////////////////////
	double result;
	char funcID[] = "callfassoc2dMF";

	gsl_function F;
	F.function = &fassoc2dMF;
	f_params params;
	params.a = bindrad;
	params.D = Dtot;
	params.k = kr;
	params.t = deltat;
	params.rho = rho;

	gsl_set_error_handler_off();
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1e6);
	F.params = reinterpret_cast<void *>(&params);

	if (kr < 1.0 / 0.0) {
		result = integrator(F, params, w, 0.0, bindrad, Dtot, kr, deltat, funcID, fassoc2dMF);

	} else {
//		result = integrator(F, params, w, r0, bindrad, Dtot, kr, deltat, funcID, fsur);
//		result = 1.0 - result;
		result = 0.0;//haven't tried this yet
	}

	gsl_integration_workspace_free(w);
	gsl_set_error_handler(NULL);
	return result;
}
