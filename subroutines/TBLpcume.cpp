#include <fstream>
#include <iostream>
#include "2Drelated.h"
using namespace std;

void TBLpcume(gsl_matrix*& mpir, double bindrad, double Dtot, double kr, double deltat, double Rmax, double RstepSize, double Rlong){
	/////////////////////MPIRR////////////////////////
	int ctr1 = 0;
	int ctr2 = 0;
	double result, error;
	const double xlow = 0;
	const double epsabs = 1e-6;
	const double epsrel = 1e-22;

	gsl_function F;
	F.function = &fpir_cum;
	f_params params;
	params.a = bindrad;
	params.D = Dtot;
	params.k = kr;
	params.t = deltat;
	
	int wspace=10000000;
	
	//cout <<"In table pcume, kr: "<<params.k<<endl;
	for (double Rindex = params.a; Rindex <= Rlong + RstepSize; Rindex += RstepSize) {
		params.r = Rindex;
		//cout <<"start with r1: "<<params.r<<endl;
		for (double R0index = params.a; R0index <= Rmax + RstepSize; R0index += RstepSize) {
			gsl_integration_workspace * w = gsl_integration_workspace_alloc(wspace);
			params.r0 = R0index;
			F.params = reinterpret_cast<void *>(&params);
			//int status=gsl_integration_qagiu(&F, xlow, epsabs, epsrel, wspace, w, &result, &error);
			gsl_integration_qagiu(&F, xlow, epsabs, epsrel, wspace, w, &result, &error);
			gsl_matrix_set(mpir, ctr1, ctr2, result*2.0*M_PI);
			gsl_integration_workspace_free(w);
			
			ctr1 = ctr1 + 1;
			
		}
		ctr1 = 0;
		ctr2 = ctr2 + 1;
		cout <<"pcume, gsl, finished r1: "<<params.r<<endl;
	}
	
}
