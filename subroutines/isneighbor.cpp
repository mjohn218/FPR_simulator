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

int isneighbor(int p1, int p2, int **neighbors) {
	/*
	Find if p1 and p2 are neighbors
	 */
	
	int i, p2test;
	int decision=0;

	for(i = 1; i <= neighbors[p1][0]; i++){
		p2test = neighbors[p1][i];
		if(p2test==p2){
			decision = 1;
		}
	}

	return decision;
}
