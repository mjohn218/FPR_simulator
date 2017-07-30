#include "GF_calls.h"
#include "Faddeeva.hh"
#include "rand_gsl.h"

void errormodel(double deltat, double Dtot, double bindrad, double kact, double &delr0, int r0bins, double *pcorrect, double &Rmax, double *pfreea) {
	char rname[300];
	sprintf(rname, "passoc_add_d%g_vs_r0.dat", deltat);
	ofstream pcfile(rname);

	int rbins; //=1000;
	double kdiff = 4.0 * M_PI * bindrad * Dtot;
	double fact = 1.0 + kact / kdiff;
	double alpha_irr = fact * sqrt(Dtot) / bindrad;
	double cof = kact / (kact + kdiff);

	double delr;
	delr = 0.01;
	double tol = 1E-17;
	double rval;
	int i;
	double ptemp;
	rval = bindrad + (0.5) * delr;
	ptemp = survive_irr(rval, deltat, Dtot, bindrad, alpha_irr, cof);
	i = 0;

	while (ptemp > tol) {
		i++;
		rval = bindrad + (i + 0.5) * delr;
		ptemp = survive_irr(rval, deltat, Dtot, bindrad, alpha_irr, cof);
		cout << "rval: " << rval << " prassoc: " << ptemp << endl;
	}
	int maxbins = i;
	cout << "MAXBINS: " << i << " rvalue: " << rval << " prvalue: " << ptemp << endl;
	/*Now use rmax to create smaller deltar*/
	double Rzero = rval;
	double range = Rzero - bindrad;
	delr = 0.01;
	rbins = int(range / delr);
	cout << "rbins: " << rbins << " delr: " << delr << endl;
	double *prassoc = new double[rbins];
	for (i = 0; i < rbins; i++) {
		rval = bindrad + (i + 0.5) * delr;
		prassoc[i] = survive_irr(rval, deltat, Dtot, bindrad, alpha_irr, cof);
	}

	i = 0;
	double tol2 = 1E-10;
	while (prassoc[i] > tol2)
		i++;
	double otherRmax = bindrad + (i + 0.5) * delr;
	cout << "Other Rmax, below tol of: " << tol2 << " r: " << otherRmax << endl;
	Rmax = otherRmax;
	//  double delr0=0.01;
	// int r0bins=5000;//int((Rmax-bindrad)/delr0);
	delr0 = (otherRmax - bindrad) / (1.0 * r0bins);
	cout << "delr0: " << delr0 << endl;
	double r0;

	pcorrection(rbins, delr, deltat, Dtot, bindrad, cof, alpha_irr, prassoc, r0bins, delr0, pcorrect, pfreea);

	for (i = 0; i < r0bins; i++) {
		r0 = bindrad + (i + 0.5) * delr0;
		pcfile << r0 << '\t' << pcorrect[i] << '\t' << pfreea[i] << endl;
	}

}
