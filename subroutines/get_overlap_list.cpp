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

void get_overlap_list(int p1, Fullmol *bases, int *ncross, int **crosspart, double bindrad, int *binlist, int MAXPERBIN, int maxnbor, int *npb, int *nbor, int *nborrev, Parms plist, double maxsep2) {

	/*Look in all neighboring cells to find proteins that are within a distance
	 close to the reaction zone plus a little padding.
	 */
	int mybin = bases[p1].mybin;
	int i;
	ncross[p1] = 0;
	int pp, p2;

	double dx, dy, dz;
	int nc1;
	double rsq;

	//  cout <<"bindrad: "<<bindrad<<" maxsep2:  "<<maxsep<<endl;

	for (pp = 0; pp < npb[mybin]; pp++) {
		p2 = binlist[mybin * MAXPERBIN + pp];
		if (bases[p2].nfree > 0) {
			if (bases[p2].protype != bases[p1].protype) {
				//check for overlap
				dx = bases[p1].xcom - bases[p2].xcom;
				dy = bases[p1].ycom - bases[p2].ycom;
				dz = bases[p1].zcom - bases[p2].zcom;
				dx -= plist.xboxl * round(dx / plist.xboxl);
				dy -= plist.yboxl * round(dy / plist.yboxl);
				dz -= plist.zboxl * round(dz / plist.zboxl);

				rsq = dx * dx + dy * dy + dz * dz;
				if (rsq < maxsep2) {
					//if they are closer than twice the bind rad, test them for overlap
					nc1 = ncross[p1];
					crosspart[p1][nc1] = p2;
					ncross[p1]++;
					//cout <<"In Diss, overlap: "<<p1<<' '<<p2<<endl;
				}
			}
		}
	}
	//now look in neighboring bins
	int hh, qq, nb;
	for (hh = 0; hh < maxnbor; hh++) {
		nb = nbor[mybin * maxnbor + hh];
		for (qq = 0; qq < npb[nb]; qq++) {
			p2 = binlist[nb * MAXPERBIN + qq];
			if (bases[p2].nfree > 0) {
				if (bases[p2].protype != bases[p1].protype) {
					//check for overlap
					dx = bases[p1].xcom - bases[p2].xcom;
					dy = bases[p1].ycom - bases[p2].ycom;
					dz = bases[p1].zcom - bases[p2].zcom;
					dx -= plist.xboxl * round(dx / plist.xboxl);
					dy -= plist.yboxl * round(dy / plist.yboxl);
					dz -= plist.zboxl * round(dz / plist.zboxl);

					rsq = dx * dx + dy * dy + dz * dz;
					if (rsq < maxsep2) {
						//if they are closer than twice the bind rad, test them for overlap
						nc1 = ncross[p1];
						crosspart[p1][nc1] = p2;
						ncross[p1]++;
						//cout <<"In Diss, overlap: "<<p1<<' '<<p2<<endl;
					}
				}
			}
		}
	} //end looping over neighboring bins
	for (hh = 0; hh < maxnbor; hh++) {
		nb = nborrev[mybin * maxnbor + hh];
		for (qq = 0; qq < npb[nb]; qq++) {
			p2 = binlist[nb * MAXPERBIN + qq];
			if (bases[p2].nfree > 0) {
				if (bases[p2].protype == bases[p1].protype) {
					//check for overlap
					dx = bases[p1].xcom - bases[p2].xcom;
					dy = bases[p1].ycom - bases[p2].ycom;
					dz = bases[p1].zcom - bases[p2].zcom;
					dx -= plist.xboxl * round(dx / plist.xboxl);
					dy -= plist.yboxl * round(dy / plist.yboxl);
					dz -= plist.zboxl * round(dz / plist.zboxl);

					rsq = dx * dx + dy * dy + dz * dz;
					if (rsq < maxsep2) {
						//if they are closer than twice the bind rad, test them for overlap
						nc1 = ncross[p1];
						crosspart[p1][nc1] = p2;
						ncross[p1]++;
						//cout <<"In Diss, overlap: "<<p1<<' '<<p2<<endl;
					}
				}
			}
		}
	} //end looping over neighboring bins

}
