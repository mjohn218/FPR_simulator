#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "reactions.h"
#include "vol_help.h"
#include "rand_gsl.h"
#include "Faddeeva.hh"
#include "utility_calls.h"
#include "vector_rot_calls.h"

void update_diffusionMF(int c1, Complex *ind_com, Fullmol *bases, Protein *wholep) {

	/*Update the diffusion of a complex bases on sum of their radii,
	 defined via the Einstein-stokes equation, such that diffusions sum via their
	 inverses.
	 */
	int i, p1, j, membnd=0;
	int size = ind_com[c1].mysize;
	double Dxinv = 0;
	double Dyinv = 0;
	double Dzinv = 0;
	double inf = 1E500;


	//find the total number of lipids attached to c1
	for (i = 0; i < size; i++) {
		p1 = ind_com[c1].plist[i];
		for (j=0; j<(bases[p1].nbnd); j++){
			if (bases[p1].partner[j]==-1){
				membnd += 1;
			}
		}
	}

	p1 = ind_com[c1].plist[0];
	ind_com[c1].Dx = bases[p1].Dx;
	ind_com[c1].Dy = bases[p1].Dy;
	ind_com[c1].Dz = bases[p1].Dz;
	for (i = 0; i < size; i++) {

		p1 = ind_com[c1].plist[i];
		//Dxinv+=1.0/bases[p1].Dx;
		//Dyinv+=1.0/bases[p1].Dy;

		if (bases[p1].Dx != 0) {
			Dxinv += 1.0 / bases[p1].Dx;
		} else {
			Dxinv = inf;
			//      ind_com[c1].Dx=0;
		}
		if (bases[p1].Dy != 0) {
			Dyinv += 1.0 / bases[p1].Dy;
		} else {
			Dyinv = inf;
			//ind_com[c1].Dy=0;
		}
		if (bases[p1].Dz != 0) {
			Dzinv += 1.0 / bases[p1].Dz;
		} else {
			Dzinv = inf;
			//ind_com[c1].Dz=0;
		}
	}
	for (j=0;j<membnd;j++){
		Dxinv += 1.0 / wholep[0].Dx;
		Dyinv += 1.0 / wholep[0].Dy;
		Dzinv = inf;
	}
	ind_com[c1].Dx = 1.0 / Dxinv;
	ind_com[c1].Dy = 1.0 / Dyinv;
	ind_com[c1].Dz = 1.0 / Dzinv;

}
