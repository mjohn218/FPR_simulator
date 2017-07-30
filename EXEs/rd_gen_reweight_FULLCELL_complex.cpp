/*
 Implement particle-based Reaction diffusion
 algorithm with trajectory reweighting.

 This general version allows any number of different
 particle types and reactions, with rigid body
 motion that includes translation and rotation.

 Periodic boundary conditions are enforced at boundaries.

 */
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "vol_help.h"
#include "reactions.h"
#include "rand_gsl.h"
#include "md_timer.h"
#include "GF_calls.h"
#include "utility_calls.h"
#include "vector_rot_calls.h"
#include "cell_neighbor_lists.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_bessel.h>
#include <2Drelated.h>
#include "Faddeeva.hh"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <iomanip>
#include <vector>
#include <string>
#include "evaluate_binding.h"


using namespace std;

struct MD_Timer totaltime;
struct MD_Timer bimoltime;

int main(int argc, char *argv[]) {

	/*Define the reactants, and all the species involved in the reactions
	 *this includes all the products, including the misbinding products
	 */
	int i, idd, j, n, k, recruitmentflag, crosscounter;
	recruitmentflag = 0;
	timeval tim;
	gettimeofday(&tim, 0);
	double t1 = tim.tv_sec + tim.tv_usec;

	vector<gsl_matrix *> contsur;
	vector<gsl_matrix *> contnorm;
	vector<gsl_matrix *> contpir;
	std::vector<std::string> infofilenames;

	int MAXALWTBL = 1000;
	cout << "Maximum number of unique 2D reactions allowed: " << MAXALWTBL << endl;
	double *TBLID=new double[MAXALWTBL*2]; //[type#][0]:ka, [tNx * Ny * (k + 1)ype#][1]:Dtot
	double Dtemp;
	double ktemp;
	int DDtableindex = 0;
	int tableexistflag = 0;
	int uniquetableindex = 0;

	int seed = int(t1);
	//seed=1445485654;
	cout << "seed: " << seed << endl;
	srand_gsl(seed);
	double randmax = pow(2.0, 32);
	double irandmax = 1.0 / randmax;

	ifstream parmfile(argv[1]);
	/*Read in number of each molecules, and their coordinates*/
//	ifstream numfile(argv[2]);
	/*Determine the constituents of the reactions*/
//	ifstream netfile(argv[3]);
//	ifstream protfile(argv[4]);
	ifstream rxnfile(argv[2]);
	char fname[100];
	char fnameComplexXYZ[100];
	char fnameProXYZ[100];

	Parms plist;
	plist.restart = 0; //in case you don't read it in.
	plist.pclath = -1; //in case you don't read it in, no clath
	plist.nloop=0;
//	read_parms(parmfile, plist); //  write_parms(plist);
	read_parmsN(parmfile, plist,infofilenames); //  write_parms(plist);

	initialize_timer(&totaltime);
	initialize_timer(&bimoltime);
	start_timer(&totaltime);

	int Nprotypes = plist.Nprotypes; //total distinct protein types, so 9
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

	int *Ncopy = new int[Nprotypes];

	int Ntotalmol = 0;
	plist.Natom = 0;
	int ntmp;
	for (i = 0; i < Nprotypes; i++) {
//		numfile >> Ncopy[i];
//		Ntotalmol += Ncopy[i];
		ifstream numfile(infofilenames[i].c_str());
		numfile.ignore(400, '\n');
		numfile >> Ncopy[i];
		Ntotalmol += Ncopy[i];
		numfile.close();
	}
	cout << "Ntotal mols: " << Ntotalmol << endl;
	plist.Ntotalmol = Ntotalmol;
	Fullmol *bases = new Fullmol[Ntotalmol]; //contains information on each protein in the full system
	Complex *ind_com = new Complex[Ntotalmol]; //contains information on each complex
	int *numpartners = new int[Nifaces]; //this should account for all free interfaces
	int **Speclist = new int*[Nifaces];

	for (i = 0; i < Nifaces; i++)
		Speclist[i] = new int[MAXPRTNER];

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

	int currentnumberofmolectypes = 0;
	int **molectypes; //define and allocate list of unique molecule types
	molectypes = new int*[MAXNUMCOMPLEXTYPES];
	for (i = 0; i < MAXNUMCOMPLEXTYPES; i++)
		molectypes[i] = new int[2+2*Nprotypes];
	for(i=0;i<MAXNUMCOMPLEXTYPES;i++){
		for(j=0;j<2+2*Nprotypes;j++){
			molectypes[i][j]=0;
		}
	}

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

//	read_protlist(Nprotypes, wholep, Nifaces, p_home, protfile, i_home);
//
//	cout << "read network: " << endl;
//	read_network(plist, numpartners, Speclist, netfile);

	int *rxtype = new int[plist.Nrxn];
//	double *Kd = new double[plist.Nrxn];
        double*** coordcont = new double**[Nprotypes];
        for(i=0; i < Nprotypes; i++){

    	    coordcont[i] = new double*[Nifaces+1];

	        for(j=0; j < Nifaces+1; j++){
	        	coordcont[i][j] = new double[3];

	            for(k=0; k < 3; k++){
	            	coordcont[i][j][k]= 0.0;
	            }
	        }
    	}

	int howmanylipids=0;
	read_inputs(rxnfile, bases, Ncopy, p_home, wholep, numpartners, Speclist,  plist, Rlist, Ncoup, mycoupled, Nmyrxn, myrxn, bindrad, kr, cntrxn, freelist, bndlist, zlist, rxtype, i_home, infofilenames, coordcont, howmanylipids);
	generate_initial_crds(plist, bases, Ncopy, ind_com, bindrad, wholep,Rlist,rxtype,p_home,howmanylipids,coordcont);
	gen_psf_system(plist, wholep, Ncopy, infofilenames);

	/*force it to have only a single partner, even though it has more than one binding interface partner*/
	int Ntotsite = 0;
	plist.Natomwrite = 0;
	for (i = 0; i < Nprotypes; i++) {
		ntmp = wholep[i].ninterface + 1;
		plist.Natom += Ncopy[i] * ntmp;
		Ntotsite += Ncopy[i] * wholep[i].ninterface;
		plist.Natomwrite += Ncopy[i] * wholep[i].nint_write;
	}
	cout << "N atoms: " << plist.Natom << " Nsites, not COM: " << Ntotsite << " Natoms to write out: " << plist.Natomwrite << endl;
	int t = 0;

//	cout << "read reactions " << endl;
//	int *rxtype = new int[plist.Nrxn];
//	double *Kd = new double[plist.Nrxn];
//	read_reactions(rxnfile, plist, Rlist, Ncoup, mycoupled, Nmyrxn, myrxn, bindrad, kr, cntrxn, freelist, bndlist, zlist, rxtype, Kd);
	/*Change rates to ka and kb, rather than kon and koff
	 unless you are reading in values of ka and kb
	 */
	//  update_rates(wholep, plist, Rlist,  bindrad, kr, rxtype, Kd, p_home);
	/*List reactant species, check reactants against the interface network*/
//	cout << "free to bind: " << endl;
//	renumber_list(cntrxn[0], freelist);
//	cout << "ready to unbind: " << endl;
//	renumber_list(cntrxn[1], bndlist);
//	cout << "free to mutate: " << endl;
//	if (cntrxn[2] > 0)
//		renumber_list(cntrxn[2], zlist);
//	cout << "Check reactions ! " << endl;
	check_reactions(plist, bases, numpartners, Speclist, Nmyrxn, Rlist, myrxn, Nspecies);

//	ifstream startfile(argv[6]);
//	cout << "now set protein status: " << endl;
//	set_status(startfile, wholep, bases, plist, Ncopy);

	/*Read in the coordinates*/
	cout << "now read in coordinates: " << endl;
	string *names = new string[Nprotypes];
//	ifstream crdfile(argv[7]);

	/*THIS CODE IS NOT USING READ IN COORDINATES, IT IS OVERWRITING, AND ASSUMES ALLL PROTEINS HAVE 1 SITE FOR BINDING SITES< 
	  AND ALL PROTEINS AVOID ONE ANOTHER.
	 */
//	read_coords(crdfile, Nprotypes, wholep, bases, ind_com, Ncopy, names);
//	generate_initial_crds_AB(plist, bases, Ncopy, ind_com, bindrad);

	double *savecrds = new double[Ntotalmol * 3]; //for x, y, z

	copy_crds(Ntotalmol, bases, savecrds);

	/*Print out specific reactions*/
	cout << "Print specific interaction network " << endl;
	int ncomplex = 0;
	for (i = 0; i < Nifaces; i++) {
		cout << i << '\t';
		for (j = 0; j < numpartners[i]; j++) {
			cout << Speclist[i][j] << '\t';
			ncomplex++;
		}
		cout << endl;
	}
	ncomplex /= 2;
	plist.nspec_complex = ncomplex;

	int ind, r1, m;
	int begin, end;

	double tau;

	double rnum;
	double rnum2, rnum3;
	int nfaces = 6; //for a cubic volume
	int direction, neighbor;

	int rxn;

	double curr_time = 0;

	int checkpoint = 10000000; /*How often to write out full current solution*/
	int stepwrite = 100; /*How often to write out species numbers*/
	char fnamemid[100];

	/*G(r) stuff*/
	double delr = 0.1;
	int nbins = 20;
	cout << "nbins: " << nbins << " Ntotal mol: " << Ntotalmol << " max g(r): " << nbins * delr << endl;
	double **gr = new double*[Nprotypes * Nprotypes];
	for (i = 0; i < Nprotypes * Nprotypes; i++)
		gr[i] = new double[nbins];

	int zeroctr=0;
	for(i=0;i<Nrxn;i++){
		if(kr[i]==0){
			zeroctr++;
		}
	}
//	MAXOVERLAP = int(ceil(MAXOVERLAP*(1+zeroctr/Nrxn)));

//	ofstream restart("restartsW.out",ios::trunc);

	char tname[100];
	sprintf(tname, "restartsW.out");
//	ofstream restart(tname);
	int newsizeoverlap = int(ceil(MAXOVERLAP*(1+(double)zeroctr/(double)Nrxn)));
	double **probvec = new double*[Ntotalmol];
	int **crosspart = new int*[Ntotalmol]; //Index of the protein partner in this reaction
	int **crossint = new int*[Ntotalmol]; //index of the interface species in this reaction
	int **cross_rxn = new int*[Ntotalmol]; //index of the reaction number of this reaction
	for (i = 0; i < Ntotalmol; i++) {
		probvec[i] = new double[newsizeoverlap];
		crosspart[i] = new int[newsizeoverlap];
		crossint[i] = new int[newsizeoverlap];
		cross_rxn[i] = new int[newsizeoverlap];
	}

	int *ncross = new int[Ntotalmol];
	int *ncrosscom=new int[Ntotalmol];
	int *movestat = new int[Ntotalmol];
	double h;

	/*The number of species does not include all separate interfaces, because it is for diffusible species*/

	/*Iterate over time steps until you hit a max time*/
	int mu;
	int flag2;
	double prob;
	double sum;
	double xchg, ychg, zchg;
	int icom;
	int pro_type, mp;
	int whichspecie;
	double rerand;
	double hfact;
	double dx, dy, dz;
	double dist2;
	double box_x = plist.xboxl;
	double box_y = plist.yboxl;
	double box_z = plist.zboxl; //nm
	double xtot, ytot, ztot;
	int nfree;
	int p, i1, i2;
	int np;
	double r2, r;
	double maxsep2 = bindrad[0] * bindrad[0]; //plist.maxsep2;
	cout << "squared distance cutoff: " << maxsep2 << endl;
	int iind, iind2, ppart;
	int twrite = plist.configwrite;
	int gwrite = plist.grwrite;
	// we have to do the following two lines so that timestat file and restart file are in time-sync
	int restartnum = (int) ceil((gwrite*1.0)/(twrite*1.0));
	gwrite = restartnum*gwrite;

	sprintf(fnameComplexXYZ, "ComplexW_COM_Np%d_Ni%d.xyz", plist.Nprotypes, Nifaces);
	ofstream compout(fnameComplexXYZ);
	sprintf(fnameProXYZ, "ProteinW_COM_Np%d_Ni%d.xyz", plist.Nprotypes, Nifaces);
	ofstream proout(fnameProXYZ);

	/*Find largest D (multiple by 2, even if larger than highest D), and largest sigma before establishing Rmax!!*/
	double Dlimit=wholep[0].Dx*2.0;
	for(i=0;i<Nprotypes;i++){
	  if(wholep[i].Dx*2.0>Dlimit)
	    Dlimit=wholep[i].Dx*2.0;
	}
	double siglimit=bindrad[0];
	for(i=0;i<Nrxn;i++){
	  if(bindrad[i]>siglimit)
	    siglimit=bindrad[i];
	}
	
	double deltat = plist.deltat;
	//double Rmax1 = 3.5 * sqrt(4.0 * D * deltat) + siglimit;
	double Rmaxlimit = 3.0 * sqrt(6.0 * Dlimit * deltat) + siglimit;
	double tx, ty, tz;

	/*cells*/
	double *M = new double[9];
	int myc1, myc2;
	/*Currently this is set to maximum number of cells, which is probably not optimal for speed.
	  Need to use largest Rmax to set cell size!!!!!!
	 */
	int Nx = int(floor(plist.xboxl / Rmaxlimit))/2; //10; //int(ceil(box_x/cellx));
	int Ny = int(floor(plist.yboxl / Rmaxlimit))/2; //10; //int(ceil(box_y/celly));
	int Nz = max(2,int(floor(plist.zboxl / Rmaxlimit))); //4; //int(ceil(box_z/cellz));
	int Ncell = Nx * Ny * Nz;
	double cellx = box_x / (Nx * 1.0);
	double celly = box_y / (Ny * 1.0);
	double cellz = box_z / (Nz * 1.0);
	cout << "Nx: " << Nx << " Ny " << Ny << " Nz " << Nz << " Ncell " << Ncell << endl;
	if (cellx < Rmaxlimit) {
	  cout << "CELL SIZE IS TOO SMALL " << " Rmax: "<<Rmaxlimit<<" cellsize: "<<cellx<<endl;
	  cout << "Nx: " << Nx << " Ny " << Ny << " Nz " << Nz << " Ncell " << Ncell << endl;
	  exit(1);
	}else if(Ncell > MAXNUMCELL){
		cout<< "Maximum number of cells " << MAXNUMCELL << " exceeded!"<< endl;
		cout<< "Scaling down number of cells" <<endl;

		Nx = max(2,int(floor(Nx * pow((MAXNUMCELL*1.0)/(Ncell*1.0),1.0/3.0))));
		Ny = max(2,int(floor(Ny * pow((MAXNUMCELL*1.0)/(Ncell*1.0),1.0/3.0))));
		Nz = max(2,int(floor(Nz * pow((MAXNUMCELL*1.0)/(Ncell*1.0),1.0/3.0))));

		Ncell = Nx * Ny * Nz;
		cout << "Nx: " << Nx << " Ny " << Ny << " Nz " << Nz << " New Ncell " << Ncell << endl;
		cellx = box_x / (Nx * 1.0);
		celly = box_y / (Ny * 1.0);
		cellz = box_z / (Nz * 1.0);
	}
	/*Each cell has 26 neighbors (unless at boundary and reflecting BC is used).
	 Since we will loop over all cell pairs, only keep track of half your neighbors.
	 */
	int maxnbor = 13; //geometry of cube
	//  int *Nnbor=new int[Ncell];
	int *nbor = new int[Ncell * maxnbor];
	int *nborrev = new int[Ncell * maxnbor];
	int *npb = new int[Ncell];
	int MAXPERBIN = 200; //int(Ntotalmol);
	int *binlist = new int[Ncell * MAXPERBIN];
	int MAXDISS = 100; //max number of proteins to dissociate in one time step
	int *disslist = new int[MAXDISS]; //holds index of dissociated proteins
	int ndiss;
	int mybin;
	int mybinind;
	int c;
	cell_neighbor_listPBCCELL(Nx, Ny, Nz, maxnbor, nbor, nborrev);
	cout << "N cell pairs (max is with PBC): " << Ncell * maxnbor << endl;

//	char tname3[100];
//	sprintf(tname3, "assembly.dat");
//	char tname2[100];
//	sprintf(tname2, "disassembly.dat");
//	ofstream disassemblyfile(tname2, std::ios::app); //disassemblyfile(tname2, std::ios::app);
//	ofstream assemblyfile(tname3, std::ios::app); // assemblyfile(tname3, std::ios::app);
	char tname4[100];
	sprintf(tname4, "time_molec_number.dat");
	ofstream timestatfile(tname4); //timestatfile(tname4, std::ios::app);
	sprintf(tname4, "time_molec_text.dat");
	ofstream timestatfiletext(tname4); //timestatfile(tname4, std::ios::app);
	sprintf(tname4, "molectypes.dat");
	ofstream molectypesfile(tname4); //timestatfile(tname4, std::ios::app);

//	assemblyfile.open(tname);
//	disassemblyfile.open(tname);

	/***************************/
	/*Begin RD simulation*/
	cout << "deltat: " << plist.deltat << endl;
	int it, iterass = 0;
	double pnormval, pirrval;
	ifstream restartf;
	plist.ntotalcomplex = Ntotalmol;
	if (plist.restart == 1) {
		/*update status of each protein and complex.*/

		restartf.open(argv[3]);
		iterass = read_restartCELL(restartf, Ntotalmol, bases, plist, ind_com, Ncopy, wholep);
		restartf.close();
		/*get complex com, complex radius, and complex diffusion*/
		update_complex_all(plist.ntotalcomplex, ind_com, bases);
	}
	int Nit = int(plist.Nit);
	if (Nit > 2.147E9) {
		cout << "ITERATIONS EXCEEDS INTEGER MAXIMUM! Exiting..." << endl;
		exit(1);
	}
	int s1;
	cout << "Ntotal complexes: " << plist.ntotalcomplex << endl;
//	write_complex(compout, plist, ind_com, 0);
//	write_protein_iface_short(proout, plist, bases, Ncopy, 0, wholep, names);
//	write_dcd(proout, plist, bases, Ncopy, 0, wholep, infofilenames);

	int amol, df;
	double us_to_s = 1E-6;
	int statwrite = plist.statwrite;

	double **traj = new double *[Ntotalmol];
	double **trajR = new double *[Ntotalmol];
	for (i = 0; i < Ntotalmol; i++) {
		traj[i] = new double[3];
		trajR[i] = new double[3];
	}

	double R1;

	double r0, passoc;



//	double Dtot = wholep[0].Dx + wholep[1].Dx;

	int veclen;

//	Dtot = 1/(1/wholep[0].Dx + 1/wholep[1].Dx + 1/wholep[2].Dx) + wholep[0].Dx;
//	veclen = sizelookup(bindrad[0], Dtot, deltat);
//	gsl_matrix * mpir0 = gsl_matrix_alloc(veclen, veclen);
//	gsl_matrix * msur0 = gsl_matrix_alloc(2, veclen);
//	gsl_matrix * mnorm0 = gsl_matrix_alloc(2, veclen);
//    DDmatrixcreate(msur0, mnorm0, mpir0, bindrad[0], Dtot, kr[0]/2/bindrad[0], deltat);
//
//	Dtot = 1/(1/wholep[0].Dx + 1/wholep[1].Dx + 1/wholep[2].Dx) + wholep[0].Dx;
//	veclen = sizelookup(bindrad[1], Dtot, deltat);
//	gsl_matrix * mpir1 = gsl_matrix_alloc(veclen, veclen);
//	gsl_matrix * msur1 = gsl_matrix_alloc(2, veclen);
//	gsl_matrix * mnorm1 = gsl_matrix_alloc(2, veclen);
//    DDmatrixcreate(msur1, mnorm1, mpir1, bindrad[1], Dtot, kr[1]/2/bindrad[1], deltat);
//
//	Dtot = 1/(1/wholep[0].Dx + 1/wholep[2].Dx) + 1/(1/wholep[0].Dx + 1/wholep[1].Dx);
//	veclen = sizelookup(bindrad[2], Dtot, deltat);
//	gsl_matrix * mpir2 = gsl_matrix_alloc(veclen, veclen);
//	gsl_matrix * msur2 = gsl_matrix_alloc(2, veclen);
//	gsl_matrix * mnorm2 = gsl_matrix_alloc(2, veclen);
//    DDmatrixcreate(msur2, mnorm2, mpir2, bindrad[2], Dtot, kr[2]/2/bindrad[2], deltat);

	int maxnbort;
	int BreakPBCduetomembraneonz = Nx * Ny * (Nz - 1);

	/*For rotating proteins*/
	double rleg2;
	double Dr1, Dr2;
	dx = bases[0].x[0] - bases[0].xcom;
	dy = bases[0].y[0] - bases[0].ycom;
	dz = bases[0].z[0] - bases[0].zcom;
	rleg2 = dx * dx + dy * dy + dz * dz;
	double rleg = sqrt(rleg2);
	double cf = cos(sqrt(6.0 * bases[0].Drx * deltat));
	Dr1 = 2.0 * rleg2 * (1.0 - cf);
	cout << "displacement2 from rotation on protein 1: " << Dr1 << " R: " << sqrt(rleg2) << " Dr " << bases[0].Drx << " diffusion from rot: " << Dr1 / (6.0 * deltat) << endl;

	/*Below use einstein-stokes to define D based on the particle
	 radius and the Temperature and viscosity.
	 To enforce the D read in from file, calculated via, e.g. bead models,
	 set scale>1 below.
	 */
	double Temp = 293; //K
	double nu = 0.001; //kg/(m*s)
	double scale = 3.0; //greater the one to correct for non-spherical
	double crad = wholep[0].radx;
	if (crad == 0)
		crad = 1;
	cout << "clathrin radius: " << crad << " nm. " << endl;
	plist.pretrans = trans_prefactor(Temp, nu, scale, crad, wholep[0].Dx);
	plist.prerot = rot_prefactor(Temp, nu, scale, crad, wholep[0].Drx);
	cout << "Diffusion prefactors: " << plist.pretrans << ' ' << plist.prerot << endl;
	cout << " Passoc_vs_separation " << endl;
	//passoc=survive_irr(bindrad[0], deltat, Dtot, bindrad[0], alpha, cof);
	//cout <<bindrad[0]<<'\t'<<passoc<<endl;
	int loop = 0;
	//i=0;
	//while(loop==0){
	// r0=(i+0.5)*0.001+bindrad[0];
	//  passoc=survive_irr(r0, deltat, Dtot, bindrad[1], alpha, cof);
	//  cout <<r0<<'\t'<<passoc<<endl;
	// if(passoc<1E-13){

	//Rmax=r0;
	//   loop=1;
	//   cout<<"passoc: "<<passoc<<" i: "<<i<<" r0: "<<r0<<endl;

	//  }
	//   i++;
	// }

	int *Nsum = new int[Nprotypes];
	Nsum[0] = 0;
	for (i = 1; i < Nprotypes; i++) {
		Nsum[i] = Nsum[i - 1] + Ncopy[i - 1];
	}
	int size = Ncopy[0] * Ncopy[1];
	cout << "Ncopy[0]: " << Ncopy[0] << endl;

	int *nprevpart = new int[Ntotalmol];
	int *ncurrpart = new int[Ntotalmol];
	int **prevlist = new int*[Ntotalmol];
	int **currlist = new int*[Ntotalmol];
	int **prevmyface = new int*[Ntotalmol];
	int **currmyface = new int*[Ntotalmol];
	int **prevpface = new int*[Ntotalmol];
	int **currpface = new int*[Ntotalmol];
	double **prevnorm = new double*[Ntotalmol];
	//int **previter=new int*[Ntotalmol];
	double **ps_prev = new double*[Ntotalmol];
	double **prevsep = new double*[Ntotalmol];
	double **currprevnorm = new double*[Ntotalmol];
	//int **previter=new int*[Ntotalmol];
	double **currps_prev = new double*[Ntotalmol];
	double **currprevsep = new double*[Ntotalmol];

	/*This size MAXNORM is keep track of how many particles are in each proteins reaction
	 zone at one time. Should be set large to buffer for fluctuations to large number
	 of partners, even though each interface (and protein) should ideally have only 1
	 partner in its reaction zone at each step.
	 */
	int MAXNORM = 200;
	int s;
	int ssave;
	for (i = 0; i < Ntotalmol; i++) {
		nprevpart[i] = 0;
		ncurrpart[i] = 0;

		currlist[i] = new int[MAXNORM];
		prevlist[i] = new int[MAXNORM];
		currmyface[i] = new int[MAXNORM];
		prevmyface[i] = new int[MAXNORM];
		currpface[i] = new int[MAXNORM];
		prevpface[i] = new int[MAXNORM];
		prevnorm[i] = new double[MAXNORM];
		ps_prev[i] = new double[MAXNORM];
		prevsep[i] = new double[MAXNORM];
		currprevnorm[i] = new double[MAXNORM];
		currps_prev[i] = new double[MAXNORM];
		currprevsep[i] = new double[MAXNORM];
		//previter[i]=new int[MAXNORM];
	}
	//cerr<<"allocated mem "<<endl;
	int myface, pface;
	int place2;
	double rtol = 1E-10;
	int flag;
	int p1, p2;
	double probvec1, p0_ratio, currnorm;

	cout << "Set prevnorms to one " << endl;
	for (i = 0; i < Ntotalmol; i++) { //previously this was Ncopy[0] by Dr. Johnson
		currlist[i][0] = 0;
		for (j = 0; j < MAXNORM; j++) {

			prevnorm[i][j] = 1.0;

			//previter[i][j]=-1;
			ps_prev[i][j] = 0;
			prevsep[i][j] = 0;
		}
	}
	/*For multi protein systems this will have to be calculated on-the-fly, because
	 Dtot will change.*/

	cout << "Rmaxlimit: " << Rmaxlimit << " Dlimit: " << Dlimit << " bindrad limit: " << siglimit << " cell size length: " << cellx << endl;
	cout << "IN GET DISTANCE, COMPARING RMAX TO THE ACTUAL DISTANCE R1, NOT JUST THE SEPARATION! " << endl;
	cout << "distance is separation beyond bind rad! " << endl;

	/*These only apply for two particle systems A and B*/
	double Vnm3 = plist.xboxl * plist.yboxl * plist.zboxl;
	double A0 = Ncopy[0];
	double B0 = Ncopy[1];
	double Keq = kr[0] * 1E6 / kr[1];
	double Bcof = (B0 / Vnm3 - A0 / Vnm3 + 1 / Keq);
	double Aeq = Vnm3 * (-Bcof / 2.0 + sqrt(Bcof * Bcof + 4.0 * A0 / Vnm3 / Keq) / 2.0);
	double A0mAeq = A0 - Aeq;
	double Nacurr;
	cout << "Aequil: " << Aeq << " A0; " << A0 << " B0: " << B0 << " V: " << Vnm3 << " Keq " << Keq << " nm^3 " << endl;

	double pact;
	double tmpx, tmpy, tmpz;

	int pnew;
	int c1, c2;
	int ci1, ci2;
	double tremain, tevent;
	int mu_ret;
	int rxn1;
	int go;
	double rate;
	int cancel;

	int nc1, nc2;
	double sep, ratio;
	double aexp;
	double bexp;

	int flagsep;
	
	int pp, qq, nb, hh;

	double maxrad = 150;
	int radbins = 3000;
	double delrad = maxrad / (1.0 * radbins);
	int ind_rad;
	double rad2, rad;

	int Ncsave = plist.ntotalcomplex;
	int Nccurr;

	double x0;
	int proa, pro2;
	double currx, curry, currz;
	/*initial separation between each pair*/

	int Nrep = 1;
	int rep = 0;
	int totpercell = 0;
	int dub;
	get_bin2(plist, bases, cellx, celly, cellz, Nx, Ny, Nz, binlist, npb, MAXPERBIN, Ncell, ind_com);
	for (i = 0; i < Ncell; i++) {
		cout << "cell: " << i << " npercell: " << npb[i] << endl;
		totpercell += npb[i];
	}
	cout << "Done defining bins. Total mols across cells " << totpercell << endl;
	if (iterass == 0) {
	  ofstream restart(tname);
      write_restart(restart, wholep, bases, plist, Ncopy, it + iterass, ind_com, deltat);
//	  write_timestat(timestatfile, wholep, bases, plist, Ncopy, iterass, ind_com, deltat, Nprotypes); ///was assemblyfile
	  write_timestat2(timestatfile, molectypesfile, timestatfiletext, wholep, bases, plist, Ncopy, it + iterass, ind_com, deltat, Nprotypes, molectypes, currentnumberofmolectypes);
	    //Have to uncomment these when you want to deal with restart files (coords data)
//		write_complex(compout, plist, ind_com, deltat*(it + iterass));
//		write_protein_iface_short(proout, plist, bases, Ncopy, (it + iterass), wholep, names);
	  write_dcd(proout, plist, bases, Ncopy, 0, wholep, infofilenames);
	  restart.close();
	}
	
	/*Begin iterations*/
	for (it = 1; it < Nit + 1; it++) {
	  
	  for (i = 0; i < Ntotalmol; i++) {
	    ncross[i] = 0;
	    ncrosscom[i]=0;
	    movestat[i] = 0;
	  }
	  
	  /*Test dissociation!! */
	  dissociate_sigmaPBC_com( Ntotalmol,  bases,  ind_com,  plist,  p_home,  i_home,  Rlist,  Nmyrxn, myrxn,  bindrad,  Ncoup,  mycoupled,  irandmax,  ndiss, disslist, ncross, kr, it, ncrosscom, movestat);
	 
   	  if (it % (gwrite * 10) == 0) { 
		  // precautionary: need to have this otherwise the restart file gets garbled when we're saving for a very large number of molecules and simulation abruptly were terminated
	    std::ifstream src0("restartsW.out", std::ios::binary);
	    std::ofstream dst0("SAVEDrestartsW.out", std::ios::binary);
	    dst0 << src0.rdbuf();
	    std::ifstream src1("timestat.dat", std::ios::binary);
	    std::ofstream dst1("SAVEDtimestat.dat", std::ios::binary);
	    dst1 << src1.rdbuf();

//	    Have to uncomment these when you want to deal with restart files (coords data)
//	    std::ifstream src2(fnameComplexXYZ, std::ios::binary);
//	    std::ofstream dst2("SAVEDfnameComplexXYZ.dat", std::ios::binary);
//	    dst2 << src2.rdbuf();
//	    std::ifstream src3(fnameProXYZ, std::ios::binary);
//	    std::ofstream dst3("SAVEDfnameComplexXYZ.dat", std::ios::binary);
//	    dst3 << src3.rdbuf();
	  }  
	  
	  get_bin2(plist, bases, cellx, celly, cellz, Nx, Ny, Nz, binlist, npb, MAXPERBIN, Ncell, ind_com);
	  
	  /*Either remove dissociated proteins from the binlist
	    so they do not try to react again this time step, or
	    ensure that they do not move again (D=0) and they
	    will not react, the other proteins will just avoid
	    overlapping with them.
	  */
	  for (i = 0; i < ndiss; i++) {
	    p1 = disslist[i];
	    k = bases[p1].mycomplex;
	    
	    traj[k][0]=0;
	    traj[k][1]=0;
	    traj[k][2]=0;
	    /*Comment out below to allow dissociated proteins to avoid overlap, but
	     need to update evaluate binding, to test if movestat!=2*/
	    /*
	    mybin = bases[p1].mybin;
	    mybinind = bases[p1].mybinind;
	    
	    pnew = binlist[mybin * MAXPERBIN + npb[mybin] - 1];
	    bases[pnew].mybinind = mybinind;
	    binlist[mybin * MAXPERBIN + mybinind] = pnew;
	    npb[mybin] -= 1;
	    */
	  }
	  
	  /*Keep track of partners in reaction zone for reweighting*/
	  for (i = 0; i < Ntotalmol; i++) {
	    ncurrpart[i] = 0;
	  }
	  
	  /*Measure separations between proteins in neighboring cells to identify
	    all possible reactions.
	  */
	  
	  for (c = 0; c < Ncell; c++) {
	    //			cout<<npb[c]<<endl;
	    for (pp = 0; pp < npb[c]; pp++) {
	      i = binlist[c * MAXPERBIN + pp];
	      
	      /*Test bimolecular reactions!*/
	      nfree = bases[i].nfree;
	      	      
	      if (nfree > 0) {
		/*first loop over proteins in your same cell.*/
		for (qq = pp + 1; qq < npb[c]; qq++) {
		  j = binlist[c * MAXPERBIN + qq];
		  //cout<<"protype"<<bases[j].protype<<endl;
		  
		  evaluate_binding_pair_com(i, j, bases,  ind_com,  plist,  wholep,  DDtableindex,  numpartners,  Speclist,  i_home,  ncross,  probvec,  TBLID, bindrad,  nprevpart,  prevlist,  prevmyface, prevpface,  prevnorm,  ncurrpart,  currlist,  currmyface,  currpface,  ps_prev, prevsep,  currprevnorm,  currps_prev,  currprevsep, myrxn, traj, crosspart, crossint, cross_rxn, kr, it, contsur, contnorm, contpir, MAXALWTBL, ncrosscom, movestat);
		  
		} //loop over protein partners in your same cell
		
		/*Now loop over all neighboring cells, and all proteins in those cells.
		  for PBC, all cells have maxnbor neighbor cells. For reflecting, edge have fewer.
		*/
		
		//if you are at z=Zmax layer you have only 4 neighbors
		if (c >= BreakPBCduetomembraneonz) {
		  maxnbort = 4;
		} else {
		  maxnbort = maxnbor;
		}
		
		for (hh = 0; hh < maxnbort; hh++) {
		  nb = nbor[c * maxnbor + hh];
		  for (qq = 0; qq < npb[nb]; qq++) {
		    j = binlist[nb * MAXPERBIN + qq];
		    
		    evaluate_binding_pair_com(i, j, bases,  ind_com,  plist,  wholep,  DDtableindex,  numpartners,  Speclist,  i_home,  ncross,  probvec,  TBLID, bindrad,  nprevpart,  prevlist,  prevmyface, prevpface,  prevnorm,  ncurrpart,  currlist,  currmyface,  currpface,  ps_prev, prevsep,  currprevnorm,  currps_prev,  currprevsep, myrxn, traj, crosspart, crossint, cross_rxn, kr, it, contsur, contnorm, contpir, MAXALWTBL, ncrosscom, movestat);
		    
		    
		  } //loop over all proteins in this neighbor cell
		} //loop over all neighbor cells
	      } //if protein i is free to bind
	    } //loop over all proteins in initial cell
	  } //End looping over all cells.
	  
	  /*Now that separations and reaction probabilities are calculated,
	    decide whether to perform reactions for each protein.
	  */
	  
	  for (i = 0; i < Ntotalmol; i++) {
	    
	    /*Skip any proteins that just dissociated during this time step*/
//		  disassemblyfile << i <<"\t"<<ncross[i]<<"\t";
//		  crosscounter = 0;
//		  for(int nc=0;nc<ncross[i];nc++){
//			  if (kr[cross_rxn[i][nc]]>0.0){
//				  crosscounter += 1;
//			  }
//		  }
//		  disassemblyfile <<crosscounter<<endl;

	    if (ncross[i] > 0) {

	      
	      	      	      
	      /*Evaluate whether to perform a reaction with protein i, and with whom. Flag=1 means
		reaction is performed. Returns correct ci1 and ci2 for this rxn.
	      */
	      /*Loop over all reactions individually, instead of summing probabilities*/
	      flag = choose_one_reaction_loop( i, ncross, crosspart, probvec, ci1, ci2, cross_rxn, crossint, irandmax);
	      
	      
	      if (flag == 1) {
		
		perform_association_sigma_com(i,  bases,  ind_com, plist,  crosspart,  cross_rxn, crossint,  Rlist,  i_home,  probvec,  Ncoup,  mycoupled,  p_home,  ci1,  ci2,  it, movestat, ncross, bindrad[cross_rxn[i][ci1]], ncrosscom, traj);
		
	      } else {
		
		/*This protein will not associate in this time step.
		  For this Sweeping version, just move the particle, don't check for overlap until everyone is done.
		  Store new move in the traj vector so you still know the original position store in bases[]
		*/
		
		/*movestat of zero means no traj value is selected.
		  movestat=1 means traj is selected, but particles have not moved
		  movestat=2 means particles have moved
		*/
		if (movestat[i] == 0) {
		  k = bases[i].mycomplex;
		  dx = sqrt(2.0 * deltat * ind_com[k].Dx) * GaussV();
		  dy = sqrt(2.0 * deltat * ind_com[k].Dy) * GaussV();
		  dz = sqrt(2.0 * deltat * ind_com[k].Dz) * GaussV();
		  
		  traj[k][0] = dx;
		  traj[k][1] = dy;
		  traj[k][2] = dz;
		  
		  tx = sqrt(2.0 * deltat * ind_com[k].Drx) * GaussV();
		  ty = sqrt(2.0 * deltat * ind_com[k].Dry) * GaussV();
		  tz = sqrt(2.0 * deltat * ind_com[k].Drz) * GaussV();
		  
		  trajR[k][0] = tx;
		  trajR[k][1] = ty;
		  trajR[k][2] = tz;
		  
		  rotationEuler(trajR[k][0], trajR[k][1], trajR[k][2], M);
		  reflect_traj_complex_rad_rotCELL(i, bases, ind_com, plist.xboxl, plist.yboxl, plist.zboxl, traj, M);
		  /*this displacement will apply to all the proteinss in this complex k.*/
		  for (j = 0; j < ind_com[k].mysize; j++) {
		    mp = ind_com[k].plist[j];
		    movestat[mp] = 1;
		  }
		  
		}
		
		/*Set probability of this protein to zero in all reactions so it doesn't try to
		  react again but the partners still will avoid overlapping.
		*/
		remove_one_prob_all(i, ncross, crosspart, probvec, cross_rxn, crossint);
		
	      }
	      
	    } else if (ncross[i] > -1) {
	      
	      /*this protein has ncross=0,
		meaning it neither dissociated nor tried to associate.
		however, it could have movestat=2 if it is part of a multi-protein
		complex that already displaced.
	      */
	      if (movestat[i] == 0) {
		/*don't move this protein if it already moved, or attached to someone who will
		  potentially have to resample their position (movestat=1)*/
		k = bases[i].mycomplex;
		dx = sqrt(2.0 * deltat * ind_com[k].Dx) * GaussV();
		dy = sqrt(2.0 * deltat * ind_com[k].Dy) * GaussV();
		dz = sqrt(2.0 * deltat * ind_com[k].Dz) * GaussV();
		
		traj[k][0] = dx;
		traj[k][1] = dy;
		traj[k][2] = dz;
		
		tx = sqrt(2.0 * deltat * ind_com[k].Drx) * GaussV();
		ty = sqrt(2.0 * deltat * ind_com[k].Dry) * GaussV();
		tz = sqrt(2.0 * deltat * ind_com[k].Drz) * GaussV();
		
		trajR[k][0] = tx;
		trajR[k][1] = ty;
		trajR[k][2] = tz;
		
		/*Don't automatically move the protein, in case another part of your protein complex has overlap
		 */
		rotationEuler(trajR[k][0], trajR[k][1], trajR[k][2], M);
		reflect_traj_complex_rad_rotCELL(i, bases, ind_com, plist.xboxl, plist.yboxl, plist.zboxl, traj, M);
		  
		//move_rot_proteinsPBCCELL(i, bases, ind_com, traj, movestat, trajR, M, plist);
		//reflect_complex_rad_rotCELL(i, bases, ind_com, plist.xboxl, plist.yboxl, plist.zboxl);
		//if other proteins in this complex have already moved, don't need to sample them as well;
		for(j=0;j<ind_com[k].mysize;j++){
		  mp=ind_com[k].plist[j];
		  movestat[mp]=1;
		}
	    
	      }
	      
	    }
	    
	  } //done testing all molecules for bimolecular reactions
	  
	  /*Now we have to check for overlap!!!*/
	  for (i = 0; i < Ntotalmol; i++) {
	    
	    /*
	      Now track each complex (ncrosscom), and test for overlap of all proteins in that complex before
	      performing final position updates. 
	    */
	    
	    if (ncrosscom[bases[i].mycomplex] > 0) {
	      if(movestat[i]!=2){
		/*For any protein that overlapped and did not react, check whether it overlaps with its partners,
		  do all proteins in the same complex at the same time.
		  Also, if both proteins are stuck to membrane, only do xy displacement, ignore z
		*/
		if(ind_com[bases[i].mycomplex].Dz==0){
		  //cout <<"Protein memtest: "<<i<<" ncross: "<<ncross[i]<<" ncrosscom: "<<ncrosscom[bases[i].mycomplex]<<endl;
		  sweep_separation_complex_rot_memtest_PBCCELL(deltat, i, bases, ind_com, ncross, crosspart, crossint, cross_rxn, traj, probvec, plist, movestat, i_home, bindrad, trajR, M, Rlist);
		}else
		  sweep_separation_complex_rot_PBCCELL(deltat, i, bases, ind_com, ncross, crosspart, crossint, cross_rxn, traj, probvec, plist, movestat, i_home, bindrad, trajR, M, Rlist);
		
		
	      }
	    } else if (ncrosscom[bases[i].mycomplex] == 0) {
	      
	      if (movestat[i] != 2) {
		
		/*For proteins with ncross=0, they either moved independently, or their displacements
		  were selected based on the complex they were part of, and they may not yet been moved.
		*/
		move_rot_proteinsPBCCELL(i, bases, ind_com, traj, movestat, trajR, M, plist);
		reflect_complex_rad_rotCELL(i, bases, ind_com, plist.xboxl, plist.yboxl, plist.zboxl);
		
	      }
	      
	    }
	  }
	  
	  if (it % statwrite == 0) {
	    //			Nacurr = Ncopy[0] - (Ncsave - plist.ntotalcomplex);
	    /*for clathrin. If clathrins close a loop ncomplex doesn't decreas but two more sites bind*/
	    Nccurr = Ncsave - plist.ntotalcomplex + plist.nloop; //total products formed
	    Nacurr = Ncsave * 3.0 - Nccurr * 2.0; //each product has 2 A sites in it, each clathrin has 3 As
	    //cout << "timestep: " << it * deltat << '\t' << plist.ntotalcomplex << " Nleg sites, if clathrin: " << Nacurr << " if single self: " << Ncsave - Nccurr * 2.0 << endl;
	    cout <<" timestep: "<<it*deltat<<" total complexes: "<<plist.ntotalcomplex<<endl;
	  }
	  
	  if (it % twrite == 0) {//timestat file written here
//	    write_timestat(timestatfile, wholep, bases, plist, Ncopy, it + iterass, ind_com, deltat, Nprotypes);
	    write_timestat2(timestatfile, molectypesfile, timestatfiletext, wholep, bases, plist, Ncopy, it + iterass, ind_com, deltat, Nprotypes, molectypes, currentnumberofmolectypes);
	  }
	  
	  if (it % gwrite == 0) {//restart file written here
	    ofstream restart(tname);
	    write_restart(restart, wholep, bases, plist, Ncopy, it + iterass, ind_com, deltat);
	    restart.close();
	    //Have to uncomment these when you want to deal with restart files (coords data)
//		write_complex(compout, plist, ind_com, deltat*(it + iterass));
//		write_protein_iface_short(proout, plist, bases, Ncopy, (it + iterass), wholep, names);
        write_dcd(proout, plist, bases, Ncopy, 0, wholep, infofilenames);
	  }

	  /*Now replace all currsep as prevsep, to keep track of
	    reweighting values for next time step*/
	  
	  for (i = 0; i < Ntotalmol; i++) {
	    nprevpart[i] = ncurrpart[i];
	    for (s = 0; s < ncurrpart[i]; s++) {
	      prevlist[i][s] = currlist[i][s];
	      prevmyface[i][s] = currmyface[i][s];
	      prevpface[i][s] = currpface[i][s];
	      prevnorm[i][s] = currprevnorm[i][s];
	      ps_prev[i][s] = currps_prev[i][s];
	      prevsep[i][s] = currprevsep[i][s];	      
	    }
	  }	  
  
	} //end iterating over time steps
	
	stop_timer(&totaltime);
	cout << timer_duration(totaltime) << " total time " << endl;
	
//	assemblyfile.close();
	//	disassemblyfile.close();
	
	/*Write out final result*/
	cout << "End Main, complete run " << endl;

} //end main
