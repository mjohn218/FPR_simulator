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

int choose_one_reaction(double rnum, int p1, int *ncross, int **crosspart, double **prob, int &ci1, int &ci2, int **cross_rxn, int **cross_int, double irandmax) {

	/*Compare probability of reacting with all partners to random num, if 
	 reaction probability is larger than rnum, select the individual reaction
	 who's probability spans the rnum.
	 */
	int p2, pskip;
	int i, j;

	double pmatch;
	double pnow = 0;
	double pmax = 0;
	int flag = 0;
	double rnum2;
	for (i = 0; i < ncross[p1]; i++) {
		pmax += prob[p1][i];
	}

	if (rnum > pmax) {
		flag = 0;
		/*No reaction occurs, save current partners so can test not overlap*/

	} else {
		rnum2 = rnum + rand_gsl() * irandmax; //refine the random number so we don't get exactly zero!
		if (rnum2 > pmax) {
			flag = 0;
		} else {
			/*We will have an association reaction, establish which one*/
			i = 0;
			pnow = prob[p1][i];
			while (rnum2 > pnow) {
				i++;
				pnow += prob[p1][i];
			}
			flag = 1; //we chose an association reaction
			ci1 = i;
			// if(pmax>1)
			// 	cout <<"PROB TO ASSOCIATE>1: "<<p1<<" ncross: "<<ncross[p1]<<" prob individual: "<<prob[p1][i]<<" pmax: "<<pmax<<endl;
			//cout <<"Bind protein: "<<p1<<" ncross: "<<ncross[p1]<<" pmax: "<<pmax<<" prob select: "<<prob[p1][ci1]<<" index: "<<ci1<<" partner: "<<crosspart[p1][ci1]<<endl;
			p2 = crosspart[p1][ci1];
			pmatch = prob[p1][ci1];
			for (i = 0; i < ncross[p2]; i++) {
				if (prob[p2][i] == pmatch) {
					/*it is possible that you had prob 0 for multiple partners*/
					ci2 = i;
					if (cross_rxn[p1][ci1] == cross_rxn[p2][ci2])
						i = ncross[p2]; //break from loop, otherwise the next prob 0 point could be your reactant
				}
			}
		}
	}

	return flag;
}
