#include "reactions.h"
#include "vol_help.h"

void calc_gr(Parms &plist, double delr, int nbins, Fullmol *bases, int *Ncopy, double **gr, double **Vofr) {
	/*Calculate the radial distribution function between species populations,
	 this version assumes the box boundary is reflective, so no PBC*/
	int i, j;
	int n, n2;
	double dx, dy, dz;
	double r, r2;
	int ind;
	double boxV = plist.xboxl * plist.yboxl * plist.zboxl;
	int Ntot = plist.Ntotalmol;
	int Np = plist.Nprotypes;
	double rlim = nbins * delr;
	for (i = 0; i < Np * Np; i++) {
		for (j = 0; j < nbins; j++) {
			gr[i][j] = 0;
		}
	}
	//calculate the distribution around each particle separately when no PBC
	get_volumes(nbins, bases, Vofr, plist, delr, Ntot);
	cout << "in g(r), finished volume calculation " << endl;
	int Nfree[Np];
	double vfrac1, vfrac2;
	for (i = 0; i < Np; i++)
		Nfree[i] = 0;
	int pind;
	for (i = 0; i < Ntot; i++) {
		if (bases[i].nbnd != 1) {
			//only look at free proteins
			n = bases[i].protype;
			Nfree[n]++;
			for (j = i + 1; j < Ntot; j++) {
				if (bases[j].nbnd != 1) {
					dx = bases[i].xcom - bases[j].xcom; //no PBC
					dy = bases[i].ycom - bases[j].ycom;
					dz = bases[i].zcom - bases[j].zcom;
					r2 = dx * dx + dy * dy + dz * dz;
					r = sqrt(r2);
					if (r < rlim) {
						ind = int(r / delr);
						vfrac1 = Vofr[i][ind];
						vfrac2 = Vofr[j][ind];

						n2 = bases[j].protype;
						if (n > n2)
							pind = n2 * Np + n;
						else
							pind = n * Np + n2;

						gr[pind][ind] += 1.0 / vfrac1;
						gr[pind][ind] += 1.0 / vfrac2;
					}
				}
			}
		}
	}
	//normalize by Nfree[n]*Nfree[n2]*2/V, factor of two because you add one for each protein
	//or by Nfree[n]*(Nfree-1)/V			
	double nfact;

	for (n = 0; n < Np; n++) {

		for (n2 = n; n2 < Np; n2++) {
			if (n == n2)
				nfact = Nfree[n] * (Nfree[n] - 1) / boxV;
			else
				nfact = Nfree[n] * Nfree[n2] * 2 / boxV;
			cout << "Nfree: " << Nfree[n] << " nfact: " << nfact << endl;
			for (i = 0; i < nbins; i++)
				gr[n * plist.Nprotypes + n2][i] /= nfact;
		}
	}
}
