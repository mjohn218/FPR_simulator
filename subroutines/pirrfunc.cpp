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

double pirrfunc(double bindrad, double Dtot, double kr, double deltat, double r0, double rcurr) {
	/////////////////////MPIRR////////////////////////
	double result;
	char funcID[] = "pirrfunc";

	gsl_function F;
	F.function = &fpir;
	f_params params;
	params.a = bindrad;
	params.D = Dtot;
	params.k = kr;
	params.r0 = r0;
	params.t = deltat;
	params.r = rcurr;

	gsl_set_error_handler_off();
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1e6);
	F.params = reinterpret_cast<void *>(&params);
	result = integrator(F, params, w, rcurr, bindrad, Dtot, kr, deltat, funcID, fpir);
	gsl_integration_workspace_free(w);
	gsl_set_error_handler(NULL);

	return result;
}
