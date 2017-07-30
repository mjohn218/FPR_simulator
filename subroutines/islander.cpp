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



void islander(int *movestat, int *islandsize, int *islandvec, int newsizeoverlap, int Ntotalmol, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **neighbors, int &maxislandsize, double &meanislandsize, int &islandnum) {
	/*
	Find contact islands and assign island idnumbers to each molecule
	 */
	
	int i, j, npar, partner, dumpstat = 0;
	int *visitedvec = new int[Ntotalmol];
	islandnum=0;

	for(i = 0; i < Ntotalmol; i++){
		islandvec[i]=0; // vector for islandid number assignments for each molecule, the zeroth island is for molecules with no partners at all and they have movestat=2
		visitedvec[i]=0; //vector to check if we visited this molecule along the way throughout the algorithm
	}

	for (i = 0; i < Ntotalmol; i++) { // assign cluster id's to individual molecules

		if(islandvec[i]==0){

	        npar = neighbors[i][0];
	        if(npar!=0 || movestat[i]!=2){

		        islandnum = islandnum + 1;
		        islandvec[i] = islandnum;

		        for(j=0;j<npar;j++){
		        	visitedvec[j] = neighbors[i][j+1];
		        }

		        while(npar>0){

		            partner = visitedvec[npar-1];
		            if(islandvec[partner]==0){
		                islandvec[partner] = islandnum;
		                npar = npar-1;
		                for(j=0;j<neighbors[partner][0];j++){
		                    if(islandvec[neighbors[partner][j+1]]!=islandnum){
		                    	npar = npar+1;
		                        visitedvec[npar-1] = neighbors[partner][j+1];
		                    }
		                }
		            }else{
		            	npar = npar-1;
		            }

		        }
	        }

		}

	}

	islander_stats(islandsize, islandvec, islandnum, Ntotalmol, maxislandsize, meanislandsize);

	if(dumpstat==1){//if you want to print out stats

		for (i = 0; i < Ntotalmol; i++) {
			for(j=0;j<neighbors[i][0]+1;j++){//PRINT OUT NEIGHBORLIST
				cout<<neighbors[i][j]<<'\t';
			}
			cout<<endl;
		}

		for(i = 0; i < Ntotalmol; i++){ //PRINT OUT crds
			cout<<bases[i].xcom<<' '<<bases[i].ycom<<' '<<bases[i].zcom<<endl;
		}

		for(i = 0; i < Ntotalmol; i++){ //PRINT OUT ISLANDS
			cout<<islandvec[i]<<endl;
		}

		cout<<maxislandsize<<endl;
		cout<<meanislandsize<<endl;
	}

	delete[] visitedvec;

}
