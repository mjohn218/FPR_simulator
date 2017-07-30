/*
 Generate coordinates for proteins.
 Allows many different types of (single-interface) proteins, many reactions.

 Reads in the same input files that are used to run the RD program last input is outfile name.

 compile using Makefile to include all subroutines.
 */

#include <fstream>
#include <iostream>
#include <ctime>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "rand_gsl.h"
#include "reactions.h"

using namespace std;

int main(int argc, char *argv[]) {

	int i, j, MAXPR = 100;
	double px, py, pz;
	timeval tim;
	gettimeofday(&tim, 0);
	double t1 = tim.tv_sec + tim.tv_usec;

	int seed = int(t1);
	cout << "seed: " << seed << endl;
	srand_gsl(seed);
	ifstream parmfile(argv[1]);
	Parms plist;
	plist.restart = 0; //in case you don't read it in.

	read_parms(parmfile, plist);
	int Nprotypes = plist.Nprotypes; //total distinct protein types, so 9
	int *Ncopy = new int[Nprotypes];
	/*Read in number of each molecules, and their coordinates*/
	ifstream numfile(argv[2]);
	int Ntotalmol = 0;
	plist.Natom = 0;
	int ntmp;
	for (i = 0; i < Nprotypes; i++) {
		numfile >> Ncopy[i];
		Ntotalmol += Ncopy[i];
	}
	cout << "Ntotal mols: " << Ntotalmol << endl;
	plist.Ntotalmol = Ntotalmol;
	Fullmol *bases = new Fullmol[Ntotalmol]; //contains information on each protein in the full system

	int *numpartners = new int[Nprotypes]; //assumes 1 interface per protein
	int **Speclist = new int*[Nprotypes];

	for (i = 0; i < Nprotypes; i++)
		Speclist[i] = new int[MAXPRTNER];

	ifstream netfile(argv[3]);

	double xboxl = plist.xboxl;
	double yboxl = plist.yboxl;
	double zboxl = plist.zboxl;
	double small = 0.01;

	int bflag = 1;
	int bit = 0;
	int noverlap = 0;

	int Nifaces = plist.Nifaces; //this is the number of interfaces

	int Nrxn = plist.Nrxn;
	int Nspecies = plist.Nspecies; //this will include product species

	int numcells = 1; //double sidelength=plist.boxl/1000;//micrometer: boxl is in nm
	double cellvol = plist.xboxl * plist.yboxl * plist.zboxl / (1.0 * 1E9); //sidelength*sidelength*sidelength;
	double V = cellvol * numcells; //micrometers^3
	double um_to_L = 1E15;
	double avagad = 6.022E23;
	V = V * avagad / um_to_L; //now V*Molar_Concentration gives a number of molecules
	plist.V = V; //kforward divided by this V gives per second, permolecule
	//double X0total=plist.X0total;//0.1 millimolar total concentration

	Complex *ind_com = new Complex[Ntotalmol]; //contains information on each complex

	/*Determine the constituents of the reactions*/
	ifstream protfile(argv[4]);
	char fname[100];

	Protein *wholep = new Protein[Nprotypes];
	int *p_home = new int[Nifaces]; //this reverses and tells you what protein a given interface belongs to
	int *i_home = new int[Nspecies]; //for both free and bound states, what index are you on the protein's list
	double *bindrad = new double[Nrxn]; //binding or unbinding radius for each reaction
	int *Ncoup = new int[Nrxn]; //list of reactions coupled to this one
	int **mycoupled = new int*[Nrxn];
	for (i = 0; i < Nrxn; i++)
		mycoupled[i] = new int[MAXRXN];
	/*The number of reactions is fixed and all the same reactions are possible in each spatial cell*/
	double *kr = new double[Nrxn]; //reaction rate (with dimensions)

	int **Rlist = new int*[Nrxn]; //the identity of the species in the reaction
	int *Npart = new int[Nrxn]; //The number of participant species in a reaction
	int **Del = new int*[Nrxn]; //The coeffiecients of the participants in the reaction

	/*Do not assume all reactions possible
	 *so instead, we have to figure out a way to read them in!
	 */
	int maxrctant = 5;
	for (i = 0; i < Nrxn; i++) {
		Rlist[i] = new int[maxrctant];
		Del[i] = new int[maxrctant];
	}
	int *Nmyrxn = new int[Nspecies];
	int **myrxn = new int*[Nspecies];
	for (i = 0; i < Nspecies; i++)
		myrxn[i] = new int[MAXRXN];
	int *cntrxn = new int[3]; //nfree, nbnd, nmut
	int *freelist = new int[Nspecies];
	int *bndlist = new int[Nspecies];
	int *zlist = new int[Nspecies];

	read_protlist(Nprotypes, wholep, Nifaces, p_home, protfile, i_home);

	ofstream cfile(argv[5]);
	int p;
	char *let = new char[MAXPR];
	let[0] = 'A';
	let[1] = 'B';
	let[2] = 'C';

	/*generate initial coordinates*/
	for (i = 0; i < Ntotalmol; i++) {
		bases[i].xcom = xboxl * rand_gsl() - xboxl / 2.0;
		bases[i].ycom = yboxl * rand_gsl() - yboxl / 2.0;
		if (bases[i].Dz == 0) {
			bases[i].zcom = -plist.zboxl / 2.0;
		} else {
			bases[i].zcom = zboxl * rand_gsl() - zboxl / 2.0;
		}
	}

	int t = 0;
	for (p = 0; p < Nprotypes; p++) {

		for (i = 0; i < Ncopy[p]; i++) {

			px = xboxl * rand_gsl() - xboxl / 2.0;
			py = yboxl * rand_gsl() - yboxl / 2.0;
			if (wholep[p].Dz == 0) {
				pz = -plist.zboxl / 2.0;
			} else {
				pz = zboxl * rand_gsl() - zboxl / 2.0;
			}

			if (p > 2) {
				cfile << "RC" << p - 3 << " " << px << ' ' << py << ' ' << pz << endl;

				for (j = 0; j < wholep[p].ninterface; j++) {
					cfile << "RC" << p - 3 << " " << px << ' ' << py << ' ' << pz << endl;
				}
			} else {
				cfile << let[p] << " " << px << ' ' << py << ' ' << pz << endl;

				for (j = 0; j < wholep[p].ninterface; j++) {
					cfile << let[p] << " " << px << ' ' << py << ' ' << pz << endl;
				}
			}
			t++;
		}
	}

}
