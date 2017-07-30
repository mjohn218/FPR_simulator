#include "2Drelated.h"

void TBLpcume_allowerr(gsl_matrix*& mpir, int r0bins, int r1bins, double bindrad, double Dtot, double kr, double deltat, double Rmax, double RstepSize, double Rlong, gsl_matrix *pirr_mat) {
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
	gsl_set_error_handler_off();
	int wspace=100000000;
	double value;
	double minerr=2E-5;
	double *brvalues=new double[r0bins];
	//cout <<"In table pcume, kr: "<<params.k<<endl;
	for (double Rindex = params.a; Rindex <= Rlong + RstepSize; Rindex += RstepSize) {
		params.r = Rindex;
		//cout <<"start with r1: "<<params.r<<endl;
		for (double R0index = params.a; R0index <= Rmax + RstepSize; R0index += RstepSize) {
		  // gsl_integration_workspace * w = gsl_integration_workspace_alloc(wspace);
// 			params.r0 = R0index;
// 			F.params = reinterpret_cast<void *>(&params);
// 			int status=gsl_integration_qagiu(&F, xlow, epsabs, epsrel, wspace, w, &result, &error);
// 			cout <<"status: "<<status<<" error: "<<error<<'\t';
// 			cout <<"r0: "<<params.r0<<" r: "<<params.r<<" pcume: "<<result*2.0*M_PI<<endl;
// 			value=result*2.0*M_PI;
// 			if(error>minerr &&params.r>bindrad){
		  if(ctr2>0)value=integrate_numer(ctr1, bindrad, params.r, pirr_mat, RstepSize, ctr2+1);
		  else value=0;
		  //cout <<"CORRECTED VALUE USING NUMERICAL INTEGRATION OF PIRR "<<value+brvalues[ctr1]<<endl;
		  //	}
		  
		  gsl_matrix_set(mpir, ctr1, ctr2, value);
		  // gsl_integration_workspace_free(w);
		  //if(Rindex==params.a)brvalues[ctr1]=value;
		  ctr1 = ctr1 + 1;
			
		}
		ctr1 = 0;
		ctr2 = ctr2 + 1;
		
	}
	delete[] brvalues;
}
