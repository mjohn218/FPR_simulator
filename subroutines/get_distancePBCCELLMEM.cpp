#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "reactions.h"
#include "vol_help.h"
#include "2Drelated.h"
#include "rand_gsl.h"
#include "Faddeeva.hh"
#include "utility_calls.h"
#include "vector_rot_calls.h"

void get_distancePBCCELLMEM(int *numpartners, int **Speclist, int **myrxn, double mydensity, Protein *wholep, Complex *ind_com, double *kr, Fullmol *bases, double deltat, double *bindrad, int *ncross, int **crosspart, double **prob, int **crossint, int **cross_rxn, Parms &plist) {

	int i, j, nc1, mu, myc, i1, i2, p, p1, n, np;

	double fact, alpha, aexp, bexp, kact, Dtot, kdiff, probvec1, sep, dz, Rmax;
	double fourpi = 4.0 * M_PI;

	//what is the index of lipid's interface

	for (p1=0;p1<plist.Ntotalmol;p1++){

		for (p = 0; p < bases[p1].nfree; p++) {
			i1 = bases[p1].freelist[p];
			/*test all of i1's binding partners to see whether they are on protein j */
			np = numpartners[i1];
			for (n = 0; n < np; n++) {
				i2 = Speclist[i1][n]; //binding interface partner
				if(i2==0){ //PIP2 interface
					myc = bases[p1].mycomplex;
					if (ind_com[myc].Dz>0){ // solution complex

						mu = myrxn[i1][n];
						probvec1 = 0.0;
						dz = ind_com[myc].zcom - (-plist.zboxl)/2.0;
						Dtot = ind_com[myc].Dx + wholep[0].Dx;

						Rmax = bindrad[mu]+3*sqrt(6*deltat*Dtot);

						/*Rmax should be the binding radius plus ~max diffusion distance, using 3*sqrt(6*Dtot*deltat)*/
						if (sqrt(dz * dz) < Rmax) {

							if (dz<=bindrad[mu]){
								dz=bindrad[mu];
							}

							/*in this case we evaluate the probability of this reaction*/
							nc1 = ncross[p1];
							crosspart[p1][nc1] = -1; //MEMBRANE
							crossint[p1][nc1] = i1;

							cross_rxn[p1][nc1] = mu;


							kdiff = fourpi * Dtot * bindrad[mu];
							kact = kr[mu];
							fact = 1.0 + kact / kdiff;
							alpha = fact * sqrt(Dtot) / bindrad[mu];
							//							aexp = sep / sqrt(4.0 * Dtot * deltat);
							//							bexp = alpha * sqrt(deltat);

							probvec1 = sqrt(Dtot)*exp(alpha*(alpha*sqrt(Dtot)*deltat-bindrad[mu]+dz)/sqrt(Dtot))*erfc((2*alpha*sqrt(Dtot)*deltat-bindrad[mu]+dz)/(2*sqrt(Dtot*deltat)))/alpha;
							probvec1 = probvec1 + (bindrad[mu]-dz-sqrt(Dtot)/alpha)*erfc((dz-bindrad[mu])/(2*sqrt(Dtot*deltat)));
							probvec1 = probvec1 + 2*sqrt(Dtot*deltat)*exp(-(bindrad[mu]-dz)*(bindrad[mu]-dz)/(4*Dtot*deltat))/sqrt(M_PI);
							probvec1 = 2*M_PI*bindrad[mu]*kact/(kact+kdiff)*probvec1*mydensity;

							prob[p1][nc1] = probvec1;
							ncross[p1]++;

						}
					}else {

						mu = myrxn[i1][n];
						Dtot = ind_com[myc].Dx + wholep[0].Dx;

						nc1 = ncross[p1];
						crosspart[p1][nc1] = -1; //DOUBLY MEMBRANE BOUND
						crossint[p1][nc1] = i1;
						cross_rxn[p1][nc1] = mu;

						prob[p1][nc1] = callfassoc2dMF(bindrad[mu], Dtot, kr[mu], deltat, mydensity); //callfassoc2dMF(bindrad[mu], Dtot, ktemp, deltat, mydensity);
						//THIS IS INEFFICIENT FIGURE OUT A WAY TO STORE THIS VALUE FOR FUTURE USE AND AVOID REDUNDANT RECALCULATIONS
						ncross[p1]++;

					}
				}
			}
		}
	}
}
