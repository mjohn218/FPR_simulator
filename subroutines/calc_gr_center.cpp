#include "reactions.h"
#include "vol_help.h"

void calc_gr_center(Parms plist, double delr, int nbins, Fullmol *bases, double *gr, double bindrad) {
	/*Calculate the radial distribution function between species populations,
	 this version assumes the box boundary is reflective, so no PBC*/
	int i, j;
	int n, n2;
	double dx, dy, dz;
	double r, r2;
	int ind;

	double rlim = (nbins - 1) * delr + bindrad;
	//  for(i=0;i<Np*Np;i++){
	// for(j=0;j<nbins;j++){
//     gr[j]=0;
//   }
	//calculate the distribution around each particle separately when no PBC
	//  get_volumes(nbins, bases, Vofr, plist, delr, Ntot);

	int pind;
	j = 0;
	double tol = 1E-6;
	for (i = 1; i < plist.Ntotalmol; i++) {
		//if(bases[i].nbnd!=1){
		//only look at free proteins
		/*Distance is to the origin!*/

		dx = bases[i].xcom; //-bases[j].xcom;//no PBC
		dy = bases[i].ycom; //-bases[j].ycom;
		dz = bases[i].zcom; //-bases[j].zcom;
		r2 = dx * dx + dy * dy + dz * dz;
		r = sqrt(r2);
		if (r < tol) {
			gr[0] += 1.0;

		} else if (r < rlim) {
			ind = int((r - bindrad) / delr) + 1;

			gr[ind] += 1.0;
			//cout <<"final r: "<<r<<" index: "<<ind<<" delr: "<<delr<<" bindrad: "<<bindrad<<" particle: "<<i<<endl;
		}
	}
	//normalize by Nfree[n]*Nfree[n2]*2/V, factor of two because you add one for each protein
	//or by Nfree[n]*(Nfree-1)/V			

}
