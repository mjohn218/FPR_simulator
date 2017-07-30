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

void remove_protein_reaction(int p1, int *ncross, int **crosspart, double **prob, int **cross_rxn, int **cross_int) {

	/*Remove just p1 from the list of all binding partners*/
	int pskip;
	int i, j;

	double pmatch;

	for (i = 0; i < ncross[p1]; i++) {

		pskip = crosspart[p1][i];

		/*Need to delete protein p1 from the pskip's list    */

		for (j = 0; j < ncross[pskip]; j++) {
			if (crosspart[pskip][j] == p1) {
				crosspart[pskip][j] = crosspart[pskip][ncross[pskip] - 1];
				prob[pskip][j] = prob[pskip][ncross[pskip] - 1];
				cross_rxn[pskip][j] = cross_rxn[pskip][ncross[pskip] - 1];
				cross_int[pskip][j] = cross_int[pskip][ncross[pskip] - 1];
				ncross[pskip] -= 1;
				j -= 1;
			}
		}
	}

	ncross[p1] = -1; //this protein is done

}
