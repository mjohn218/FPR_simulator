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



void islander_stats(int *islandsize, int *islandvec, int islandnum, int Ntotalmol, int &maxislandsize, double &meanislandsize) {
	/*
	Find cluster statistics: islandmembership, maxislandsize and meanislandsize
	 */
	
	int i, j;
	meanislandsize = 0.0;
	maxislandsize = 0;

	for(i = 1; i < islandnum+1; i++){
		islandsize[i] = 0;
		for(j = 0; j < Ntotalmol; j++){
			if(islandvec[j]==i){
				islandsize[i] += 1;
			}
		}
		meanislandsize += (1.0*islandsize[i])/islandnum;
		if(islandsize[i]>maxislandsize){
			maxislandsize = islandsize[i];
		}

	}

}
