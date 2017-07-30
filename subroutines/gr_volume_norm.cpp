#include "reactions.h"
#include "vol_help.h"

void gr_volume_norm(double *gr, int nbins, double delr, double rho0, int Nrep, double bindrad, int flagdim) {
	int i;
	double nfact;
	double vb, nideal;
	double r;
	double rup, rdown;
	cout << " gr[0]: " << gr[0] << endl;
	gr[0] = gr[0] / (Nrep * 1.0); //particles at the origin

	for (i = 1; i < nbins; i++) {
		r = delr * (i - 0.5) + bindrad;
		//2D
		if (flagdim == 2) {
			rup = (i) * delr + bindrad;
			rdown = (i - 1) * delr + bindrad;
			vb = rup * rup - rdown * rdown;
			nideal = M_PI * vb * rho0;
		} else {
			//3D
			rup = (i) * delr + bindrad;
			rdown = (i - 1) * delr + bindrad;

			vb = rup * rup * rup - rdown * rdown * rdown;
			nideal = 4.0 / 3.0 * M_PI * vb * rho0;
			//    cout <<"bin: "<<i<<" gr: "<<gr[i]<<" r: "<<r<<" vb: "<<vb<<" nideal: "<<nideal<<" nrep; "<<Nrep<<endl;

		}
		gr[i] = gr[i] / (nideal * Nrep);
	}
}
