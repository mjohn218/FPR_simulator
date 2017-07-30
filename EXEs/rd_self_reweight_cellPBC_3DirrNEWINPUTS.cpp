/*


Implement particle-based Reaction diffusion
algorithm with trajectory reweighting.

Allows only two types of particles A and B with
one reversible reaction.

Volume is split into cells to speed calculation.

In this SKIP version, A and B do not avoid positions with any bound complexes,
only with other free particles.

  
Periodic Boundary conditions at box edges.  
  
*/
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
#include "cell_neighbor_lists.h"
#include "Faddeeva.hh"

using namespace std;

struct MD_Timer totaltime;
struct MD_Timer bimoltime;

int main(int argc, char *argv[])
{
  
	/*Define the reactants, and all the species involved in the reactions
	 *this includes all the products, including the misbinding products
	 */
	int i, idd, j, n, k, recruitmentflag, crosscounter;
	recruitmentflag = 0;
	timeval tim;
	gettimeofday(&tim, 0);
	double t1 = tim.tv_sec + tim.tv_usec;

//	vector<gsl_matrix *> contsur;
//	vector<gsl_matrix *> contnorm;
//	vector<gsl_matrix *> contpir;
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
  
//  cout <<"read reactions "<<endl;
//  int *rxtype=new int[plist.Nrxn];
//  double *Kd=new double[plist.Nrxn];
//  read_reactions(rxnfile, plist, Rlist, Ncoup, mycoupled, Nmyrxn, myrxn, bindrad, kr, cntrxn, freelist, bndlist, zlist,  rxtype, Kd);
//  /*Change rates to ka and kb, rather than kon and koff
//    unless you are reading in values of ka and kb
//  */
//  //  update_rates(wholep, plist, Rlist,  bindrad, kr, kdiss2, kdiss3, rxtype, clrxn, Kd, Kd2, Kd3, p_home);
//
//  /*List reactant species, check reactants against the interface network*/
//  cout <<"free to bind: "<<endl;
//  renumber_list(cntrxn[0], freelist);
//  cout <<"ready to unbind: "<<endl;
//  renumber_list(cntrxn[1], bndlist);
//  cout <<"free to mutate: "<<endl;
//  if(cntrxn[2]>0)
//    renumber_list(cntrxn[2], zlist);
  cout <<"Check reactions ! "<<endl;
  check_reactions(plist, bases, numpartners, Speclist, Nmyrxn, Rlist, myrxn, Nspecies);
  
  
    
//  ifstream startfile(argv[6]);
//  cout <<"now set protein status: "<<endl;
//  set_status(startfile,wholep,bases, plist, Ncopy);
  /*Read in the coordinates*/
  cout <<"now read in coordinates: "<<endl;
  string *names=new string[Nprotypes];
//  ifstream crdfile(argv[7]);
//  read_coords(crdfile, Nprotypes, wholep, bases, ind_com, Ncopy, names);
//  cout <<"GENERATING INITIAL CONFIG ! "<<endl;
//  generate_initial_crds_AB_allavoid(plist, bases,Ncopy, ind_com, bindrad);
  
  double *savecrds=new double[Ntotalmol*3];//for x, y, z

  copy_crds(Ntotalmol, bases, savecrds);
  
  /*Print out specific reactions*/
  cout <<"Print specific interaction network "<<endl;
  int ncomplex=0;
  for(i=0;i<Nifaces;i++){
    cout <<i<<'\t';
    for(j=0;j<numpartners[i];j++){
      cout <<Speclist[i][j]<<'\t';
      ncomplex++;
    }
    cout <<endl;
  }
  ncomplex/=2;
  plist.nspec_complex=ncomplex;
  
      
  int ind, r1, m;
  int begin, end;
  
  double tau;

  double rnum;
  double rnum2, rnum3;
  int nfaces=6;//for a cubic volume
  int direction, neighbor;
  
  int rxn;
  
  double curr_time=0;

  int checkpoint=10000000; /*How often to write out full current solution*/
  int stepwrite=100; /*How often to write out species numbers*/
  char fnamemid[100];


  /*G(r) stuff*/
  int nbins=800;//int(plist.xboxl/2.0/delr);
  double grMax=30.0;
  double delr=(grMax-bindrad[0])/(1.0*nbins);//0.1;
  nbins+=1;//for zero spot
  cout <<"nbins gr: "<<nbins<<" Ntotal mol: "<<Ntotalmol<<" max g(r): "<<nbins*delr<<" delr "<<delr<<endl; 
  
  
  //int tpts=100;//Nit/grfreq;
  int tind;
  double tdub, targ;
  double Mintime=1E-3;//in microseconds
  double Maxtime=plist.Nit*plist.deltat;
  double deltat=plist.deltat;
  
  int Ndecades=int((log10(Maxtime)-log10(Mintime)))+1;
  int mindecade=int(-floor(log10(Mintime)));
  double *mydec=new double[Ndecades];
  int *nperdec=new int[Ndecades];
  double delbin;
  double dec_dub;
  int *binsum=new int[Ndecades];
  cout <<"Ndecades: "<<Ndecades<<" Maxtime: "<<Maxtime<<" us. "<<" mindecade: "<<mindecade<<endl;
  int tot=0;
  binsum[0]=0;
  int decroot;
  for(i=0;i<Ndecades;i++){
    mydec[i]=Mintime*pow(10.0, i);
    decroot=int(mydec[i]/deltat);
    nperdec[i]=9;
    // if(i<4)nperdec[i]=9;
//     else if(i>9) nperdec[i]=1000;
//     else nperdec[i]=100;
    tot+=nperdec[i];
    binsum[i+1]=binsum[i]+nperdec[i];
    cout <<"i: "<<i<<" mydec:" <<mydec[i]<<" nperdec:" <<nperdec[i]<<" sum to mystart bin: "<<binsum[i]<<" decroot: "<<decroot<<endl;
  }
  int dec;
  int dstep=1;
  int Nperdec=9*dstep;//1 through 9
  
  int tbins=tot+1;//Ndecades*Nperdec+1;//for the zero bin
  double *tvals=new double[tbins];
  //  int grfreq=plist.grfreq;
  int nn;
  i=0;
  int itcurr;
  double trev, trev2;
  tvals[i]=0;
  i++;
  for(n=0;n<Ndecades;n++){
    delbin=9.0/(1.0*nperdec[n]);
    for(nn=0;nn<nperdec[n];nn++){
      trev=i-1-binsum[n];
      trev2=(trev*delbin+1)*mydec[n];
      itcurr=int(round(trev2/deltat));
      tvals[i]=trev2;
      cout <<"Decade: "<<n<<" nperdec: "<<nperdec[n]<<" time: "<<tvals[i]<<" it:" <<itcurr<<" time at it: "<<itcurr*deltat<<endl;
      i++;
    }
  }
      
  // int tpts=Nit/grfreq+1;
  double rho0=1.0/(plist.xboxl*plist.yboxl*plist.zboxl);//don't normalize by particle numbers currently.
  cout <<"N gr time pts: "<<tbins<<" Gr time frequency, logsave 1/V: "<<rho0<<endl;
  //tpts+=1;
  double **gr=new double*[tbins];//
  double **grnn=new double*[tbins];//
  for(i=0;i<tbins;i++){
    gr[i]=new double[nbins];
    grnn[i]=new double[nbins];
  }
  for(i=0;i<tbins;i++){
    init_gr_zero(gr[i], nbins);
    init_gr_zero(grnn[i], nbins);
  }



  double **probvec=new double*[Ntotalmol];
  int **crosspart=new int*[Ntotalmol];//Index of the protein partner in this reaction
  int **crossint=new int*[Ntotalmol];//index of the interface species in this reaction
  int **cross_rxn=new int*[Ntotalmol];//index of the reaction number of this reaction
  for(i=0;i<Ntotalmol;i++){
    probvec[i]=new double[MAXOVERLAP];
    crosspart[i]=new int[MAXOVERLAP];
    crossint[i]=new int[MAXOVERLAP];
    cross_rxn[i]=new int[MAXOVERLAP];
  }
    
  int *ncross=new int[Ntotalmol];
  int *ncrosscom=new int[Ntotalmol];
  int *movestat=new int[Ntotalmol];
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
  double box_x=plist.xboxl;
  double box_y=plist.yboxl;
  double box_z=plist.zboxl;//nm
  double xtot, ytot, ztot;
  int nfree, wprot, wprot2;
  int p, i1, i2;
  int np;
  double r2, r;
  double maxsep2=bindrad[0]*bindrad[0];//plist.maxsep2;
  cout <<"squared distance cutoff: "<<maxsep2<<endl;
  int iind, iind2, ppart;
  

  sprintf(fname,"ComplexW_COM_Np%d_Ni%d.xyz",plist.Nprotypes, Nifaces); 
  ofstream compout(fname);
  sprintf(fname,"ProteinW_COM_Np%d_Ni%d.xyz",plist.Nprotypes, Nifaces); 
  ofstream proout(fname);
  /*cells*/
  double Dtot=2.0*wholep[0].Dx;

  double Rmax1=3.0*sqrt(6.0*Dtot*deltat)+bindrad[0];
  
  
  int Nx=20;//int(ceil(box_x/cellx));
  int Ny=20;//int(ceil(box_y/celly));
  int Nz=20;//int(ceil(box_z/cellz));
  int Ncell=Nx*Ny*Nz;
  double cellx=box_x/(Nx*1.0);
  double celly=box_y/(Ny*1.0);
  double cellz=box_z/(Nz*1.0);
  if(cellx<Rmax1){
    
    cout <<"CELL SIZE IS TOO SMALL "<<" resize"<<endl;
    Nx=int(box_x/(Rmax1+1.0));
    Ny=int(box_y/(Rmax1+1.0));
    Nz=int(box_z/(Rmax1+1.0));
    Ncell=Nx*Ny*Nz;
    cellx=box_x/(Nx*1.0);
    celly=box_y/(Ny*1.0);
    cellz=box_z/(Nz*1.0);
        //exit(1);
  }
  int c;
  cout <<"Nx: "<<Nx<<" Ny "<<Ny<<" Nz "<<Nz<<" Ncell "<<Ncell<<endl;

  /*Each cell has 26 neighbors (unless at boundary and reflecting BC is used). 
    Since we will loop over all cell pairs, only keep track of half your neighbors.
  */

  int maxnbor=13; //geometry of cube
  //  int *Nnbor=new int[Ncell];
  int *nbor=new int[Ncell*maxnbor];
  int *nborrev=new int[Ncell*maxnbor];

  int *npb=new int[Ncell];
  int MAXPERBIN=200;//int(Ntotalmol);
  int *binlist=new int[Ncell*MAXPERBIN];
  int MAXDISS=100;//max number of proteins to dissociate in one time step 
  int *disslist=new int[MAXDISS];//holds index of dissociated proteins
  int ndiss;
  int mybin;
  int mybinind;

  cell_neighbor_listPBC(Nx, Ny, Nz, maxnbor, nbor, nborrev);

  cout <<"N cell pairs (max is with PBC): "<<Ncell*maxnbor<<endl;
  
  

  /***************************/
  /*Begin RD simulation*/
  cout <<"deltat: "<<plist.deltat<<endl;

  int it;
  ifstream restartf;
  plist.ntotalcomplex=Ntotalmol;
  if(plist.restart==1){
    /*update status of each protein and complex.*/
    restartf.open(argv[8]);
    read_restart(restartf,  Ntotalmol,  bases, plist, ind_com);
    
    restartf.close();
    /*get complex com, complex radius, and complex diffusion*/
    update_complex_all(plist.ntotalcomplex, ind_com, bases);
  }
  
  int Nit=int(plist.Nit);
  if(Nit>2.147E9){
    cout <<"ITERATIONS EXCEEDS INTEGER MAXIMUM! Exiting..."<<endl;
    exit(1);
  }
  int s1;
  cout <<"Ntotal complexes: "<<plist.ntotalcomplex<<endl;
  
  write_complex(compout, plist, ind_com, 0);
  write_protein_iface_short(proout, plist, bases,Ncopy,0, wholep, names);
  
  int amol,df; 
  double us_to_s=1E-6;
  int statwrite=plist.statwrite;
  int instwrite=statwrite*100;
  //  ofstream statfile("Amols.out");
  double **traj=new double *[Ntotalmol];
  for(i=0;i<Ntotalmol;i++)
    traj[i]=new double[3];
  

  double kpi=4.0*M_PI*bindrad[0];
  double fourpi=4.0*M_PI;
  
  double kact=kr[0];
  double kdiff=4.0*M_PI*Dtot*bindrad[0];
  double fact=1.0+kr[0]/kdiff;
  
  double alpha=fact*sqrt(Dtot)/bindrad[0];
  double cof=kr[0]/(kr[0]+kdiff);
  /*Update all the rates to be ka and kb, rather than kon and koff*/
  cout <<" alpha: "<<alpha<<" cof: "<<cof<<" fact: "<<fact<<" "<<endl;
  cout <<"activation rate in nm^3/us: "<<kr[0]<<" kd: "<<kr[1]<<" Kd, uM: "<<kr[1]/kr[0]/1E6/6.022E23*1E24*1E6<<endl;
  cout <<"Dtot: "<<Dtot<<endl;
  //  check_bndlist(plist.ntotalcomplex, ind_com, bases);
  double R1;

  double Rmax=Rmax1;
  double r0, passoc;
  cout <<" Passoc_vs_separation "<<endl; 
  passoc=survive_irrF(bindrad[0], deltat, Dtot, bindrad[0], alpha, cof);
  cout <<bindrad[0]<<'\t'<<passoc<<endl;
  int loop=0;
  i=0;
  while(loop==0){
    r0=(i+0.5)*0.001+bindrad[0];
    passoc=survive_irrF(r0, deltat, Dtot, bindrad[0], alpha, cof);
    cout <<r0<<'\t'<<passoc<<endl;
    if(passoc<1E-13){
            
      loop=1;
      //cout<<"passoc: "<<passoc<<" i: "<<i<<" r0: "<<r0<<endl;
      
    }
    i++;
  }
  

  int *Nsum=new int[Nprotypes];
  Nsum[0]=0;
  for(i=1;i<Nprotypes;i++){
    Nsum[i]=Nsum[i-1]+Ncopy[i-1];
  }
  int size=Ncopy[0]*Ncopy[1];
  cout <<"Ncopy[0]: "<<Ncopy[0]<<endl;
  int *nprevpart=new int[Ncopy[0]];
  int *ncurrpart=new int[Ncopy[0]];
  int **prevlist=new int*[Ncopy[0]];
  int **currlist=new int*[Ncopy[0]];
  double **prevnorm=new double*[Ncopy[0]];
  //int **previter=new int*[Ncopy[0]];
  double **ps_prev=new double*[Ncopy[0]];
  double **prevsep=new double*[Ncopy[0]];
  double **currprevnorm=new double*[Ncopy[0]];
  //int **previter=new int*[Ncopy[0]];
  double **currps_prev=new double*[Ncopy[0]];
  double **currprevsep=new double*[Ncopy[0]];
  int MAXNORM=200;
  int s;
  int ssave;
  for(i=0;i<Ncopy[0];i++){
    nprevpart[i]=0;
    ncurrpart[i]=0;
    
    currlist[i]=new int[MAXNORM];
    prevlist[i]=new int[MAXNORM];
    prevnorm[i]=new double[MAXNORM];
    ps_prev[i]=new double[MAXNORM];
    prevsep[i]=new double[MAXNORM];
    currprevnorm[i]=new double[MAXNORM];
    currps_prev[i]=new double[MAXNORM];
    currprevsep[i]=new double[MAXNORM];
    //previter[i]=new int[MAXNORM];
  }
  //cerr<<"allocated mem "<<endl;
  int place2;
  double rtol=1E-10;
  int flag;
  int p1, p2;
  double probvec1, p0_ratio, currnorm;
  std::complex<double> z;
  std::complex<double> value;
  double term2, ep1, relerr;
  double ea2;
  
  cout <<"Set prevnorms to one "<<endl;
  for(i=0;i<Ncopy[0];i++){
    currlist[i][0]=0;
    for(j=0;j<MAXNORM;j++){
      
      prevnorm[i][j]=1.0;
      //previter[i][j]=-1;
      ps_prev[i][j]=0;
      prevsep[i][j]=0;
    }
  }
  /*For multi protein systems this will have to be calculated on-the-fly, because
    Dtot will change.*/

  cout <<"Rmax1: "<<Rmax1<<" Dtot: "<<Dtot<<" Rmaxfinal: "<<Rmax<<" cell size length: "<<cellx<<endl;
  cout <<"IN GET DISTANCE, COMPARING RMAX TO THE ACTUAL DISTANCE R1, NOT JUST THE SEPARATION! "<<endl;
  cout <<"error model "<<endl;
  cout <<"distance is separation beyond bind rad! "<<endl;
  
  /*These only apply for two particle systems A and B*/
  double Vnm3=plist.xboxl*plist.yboxl*plist.zboxl;
  double A0=Ncopy[0];
  double B0=Ncopy[1];
  double Keq=kr[0]*1E6/kr[1];
  double Bcof=(B0/Vnm3-A0/Vnm3+1/Keq);
  double Aeq=Vnm3*(-Bcof/2.0+sqrt(Bcof*Bcof+4.0*A0/Vnm3/Keq)/2.0 );
  double A0mAeq=1.0;//A0-Aeq;
  double Nacurr;
  cout <<"Aequil: "<<Aeq<<" A0; "<<A0<<" B0: "<<B0<<" V: "<<Vnm3<<" Keq "<<Keq<< " nm^3 "<<endl;
  

  int flagdim=3;
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
  double memfact=1.0;//If you are binding to PIP2, this should be 1/2.
  int pp, qq, nb, hh;
  
    
  double maxrad=150;
  int radbins=3000;
  double delrad=maxrad/(1.0*radbins);
  int ind_rad;
  double rad2, rad;
  
  int Ncsave=plist.ntotalcomplex;
  double x0;
  int proa, pro2;
  double currx, curry, currz;
  /*initial separation between each pair*/



  int totpercell=0;
  get_bin2( plist, bases, cellx,  celly,cellz,  Nx,  Ny,  Nz, binlist, npb, MAXPERBIN, Ncell, ind_com);  
  for(i=0;i<Ncell;i++){
    cout <<"cell: "<<i<<" npercell: "<<npb[i]<<endl;
    totpercell+=npb[i];
  }
  char tname[200];
  cout <<"Done defining bins. Total mols across cells "<<totpercell<<endl;
  sprintf(tname,"A_vs_time.dat");
  ofstream avfile;
  avfile.open(tname);
  
  int Nrep=10;//plist.Nrep;
  int rep=0;

  
  double *Aavg=new double[Nit+1];
  for(i=0;i<Nit+1;i++)
    Aavg[i]=0.0;
  
  for(rep=0;rep<Nrep;rep++){
    
    cout <<"rep: "<<rep<<endl;
    
    plist.ntotalcomplex=Ncsave;
    
    for(i=0;i<Ncopy[0];i++){
      currlist[i][0]=0;
      for(j=0;j<MAXNORM;j++){
	
	prevnorm[i][j]=1.0;
	//previter[i][j]=-1;
	ps_prev[i][j]=0;
	prevsep[i][j]=0;
      }
    }
    generate_initial_crds_AB_allavoid( plist, bases, Ncopy, ind_com, bindrad);
    gen_copy_coords( Nprotypes, wholep,  bases,ind_com, Ncopy);
    
    get_bin2( plist, bases, cellx,  celly,cellz,  Nx,  Ny,  Nz, binlist, npb, MAXPERBIN, Ncell, ind_com);  
    
    for(i=0;i<Ntotalmol;i++){
      bases[i].nfree=1;
      bases[i].nbnd=0;
      
    }
    calc_gr_self(plist, delr, nbins, bases,  gr[0], bindrad[0]);
    calc_gr_self_nnorm(plist, delr, nbins, bases,  grnn[0], bindrad[0]);
    
    it=0;
    
    Nacurr=Ncopy[0]-2.0*(Ncsave-plist.ntotalcomplex);//Nacurr=Ncopy[0]-(Ncsave-plist.ntotalcomplex);
//  
    cout<<"timestep us: "<<it*deltat<<'\t'<<plist.ntotalcomplex<<'\t'<<Nacurr<<'\t'<<(Nacurr-Aeq)/A0mAeq<<endl;
    Aavg[it]+=Nacurr;

  
    for(it=1;it<Nit+1;it++){
      
      
      for(i=0;i<Ntotalmol;i++){
	ncross[i]=0;
	movestat[i]=0;
	ncrosscom[i]=0;
      }
      ndiss=0;
      
      get_bin2( plist, bases, cellx,  celly,cellz,  Nx,  Ny,  Nz, binlist, npb, MAXPERBIN, Ncell, ind_com);
            
      for(i=0;i<Ntotalmol;i++){
	ncurrpart[i]=0;
      }
      for(c=0;c<Ncell;c++){
	for(pp=0;pp<npb[c];pp++){
	  i=binlist[c*MAXPERBIN+pp];
	  
	  /*MAKE SURE YOU LLOOP OVER PROTEINS IN YOUR OWN CELL AS WELL */
	  
	  /*Test bimolecular reactions!*/
	  nfree=bases[i].nfree;
	  wprot=bases[i].protype;
	  
	  if(nfree>0){
	    
	    /*If this protein just dissociated, don't check it for crossing*/
	    
	    for(qq=pp+1;qq<npb[c];qq++){
	      j=binlist[c*MAXPERBIN+qq];
	      /*test to see if they interact.*/
	      //if(bases[j].protype!=bases[i].protype){
	      
	      if(bases[j].nfree>0){
		
		i1=bases[i].protype;
		i2=bases[j].protype;
		mu=0;
		
		/*Here now we evaluate the probability of binding*/
		//integrate_path(bases, traj, i, j, deltat, bindrad[mu], ncross, treturn, crosspart, crossint, i1, i2, mu, cross_rxn, i_home); 
		//  Dtot=ind_com[bases[i].mycomplex].Dx+ind_com[bases[j].mycomplex].Dx;
		  /*CAME UP WITH A FIXED RMAX ALREADY*/
		  //Rmax=3.0*sqrt(6.0*Dtot*deltat);
		nc1=ncross[i];
		nc2=ncross[j];
		flagsep=get_distancePBC(bases, traj, i, j, deltat, bindrad[mu], ncross,  crosspart, crossint, i1, i2, mu, cross_rxn, i_home, it, Rmax, sep, R1, plist); 
		//cout <<"trying proteins i and j: "<<i<<' '<<j<<' '<<flagsep<<" Rmax: "<<Rmax<<" R1 "<< R1<<" rsep: "<<sep<<endl;
		if(flagsep==1){
		  //evaluate probability
		  ncrosscom[bases[i].mycomplex]++;
		  ncrosscom[bases[j].mycomplex]++;
		  
		  //sep=R1-bindrad[mu];
		  ratio=bindrad[mu]/R1;
		  if(sep<0){
		    cout <<"separation <0: "<<sep<<" r1 "<<R1<<" p1: "<<i<<" p2: "<<j<<" Rep: "<<rep<<" it "<<it<<endl;
		    cout <<bases[i].xcom<<' '<<bases[i].ycom<<' '<<bases[i].zcom<<" nfree: "<<bases[i].nfree<<endl;
		    cout <<bases[j].xcom<<' '<<bases[j].ycom<<' '<<bases[j].zcom<<" nfree: "<<bases[j].nfree<<endl;
		    
		    sep=0;
		    ratio=1;
		    R1=bindrad[mu];
		    
		  }
		  
		  //memfact=1.0;
		  //if(bases[i].Dz==0 || bases[j].Dz==0)memfact=0.5;
		  //kdiff=fourpi*Dtot*bindrad[mu]*memfact;
		  
		  //kact=kr[mu];
		  //fact=1.0+kact/kdiff;
		  //alpha=fact*sqrt(Dtot)/bindrad[mu];
		  aexp=sep/sqrt(4.0*Dtot*deltat);
		  bexp=alpha*sqrt(deltat);
		  /*add in error model*/
		  currnorm=1.0;
		  p0_ratio=1.0;
		  
		  //place=i1*Nifaces+i2;
		  //protein i is protype wprot and j
		  wprot2=bases[j].protype;
		  proa=i;
		  pro2=j;
		  if(i>j){
		    proa=j;
		    pro2=i;
		  }
		  /*Need other way to keep track of prevnorm*/
		  for(s=0;s<nprevpart[proa];s++){
		    if(prevlist[proa][s]==pro2){
		      p0_ratio=pirr_pfree_ratio_psF(R1, prevsep[proa][s], deltat, Dtot, bindrad[mu], alpha, ps_prev[proa][s], rtol);
		      currnorm=prevnorm[proa][s]*p0_ratio;
		      
		      ssave=s;
		      s=nprevpart[proa];
		    }
		  }
		  ep1=exp(2.0*aexp*bexp+bexp*bexp);
		  if(isinf(ep1)){
		    z.real(0.0);
		    z.imag(aexp+bexp);
		    relerr=0;
		    value=Faddeeva::w(z, relerr);
		    ea2=exp(-aexp*aexp);
		    term2=ea2*real(value);
		  }else{
		    term2=ep1*erfc(aexp+bexp);
		  }
		  probvec1=ratio*kact/(kact+kdiff)*(erfc(aexp)-term2);
		  //cout <<"Ratio: "<<ratio<<" kact: "<<kact<<" kdiff "<<kdiff<<" aexp: "<<aexp<<" bexp: "<<bexp<<" term2: "<<term2<<" ep1: "<<ep1<<" probvec: "<<probvec1<<endl;
		  probvec[i][nc1]=probvec1*currnorm;
		  probvec[j][nc2]=probvec[i][nc1];
		  //if(probvec[i][nc1]>1.0)cout <<"LARGE RATIO: "<<currnorm<<" probassoc: "<<probvec1<<" ratio for this step: "<<p0_ratio<<" protein a: " <<proa<<'\t'<<pro2<<" currdist: "<<R1<<" prev: "<<prevsep[proa][ssave]<<" previous survival prob: "<<ps_prev[proa][ssave]<<" previous norm factor: "<<prevnorm[proa][ssave]<<endl;
		  //cout <<" it: "<<it<<" "<<currnorm<<" probassoc: "<<probvec1<<" ratio for this step: "<<p0_ratio<<" protein a: " <<proa<<'\t'<<pro2<<" currdist: "<<R1<<" prev: "<<prevsep[proa][ssave]<<" previous survival prob: "<<ps_prev[proa][ssave]<<" previous norm factor: "<<prevnorm[proa][ssave]<<endl;
		  s=ncurrpart[proa];
		  currprevsep[proa][s]=R1;
		  //cout <<"proa: "<<proa<<" s: "<<s<<" pro2: "<<pro2<<" currlist: "<<currlist[proa][s]<<endl;
		  currlist[proa][s]=pro2;
		  currprevnorm[proa][s]=currnorm;
		  currps_prev[proa][s]=1.0-probvec1*currnorm;
		  ncurrpart[proa]++;
		}//close enough to evaluate
	      }
	      
	    }//done looping over other proteins in your cell
	  
	    
	    for(hh=0;hh<maxnbor;hh++){
	      nb=nbor[c*maxnbor+hh];
	      for(qq=0;qq<npb[nb];qq++){
		j=binlist[nb*MAXPERBIN+qq];
		/*test to see if these two proteins react*/
		
		//if(bases[j].protype!=bases[i].protype){
		if(bases[j].nfree>0){
		  
		  i1=bases[i].protype;
		  i2=bases[j].protype;
		  mu=0;
		  
		  /*Here now we evaluate the probability of binding*/
		  
		  //Dtot=ind_com[bases[i].mycomplex].Dx+ind_com[bases[j].mycomplex].Dx;
		  /*CAME UP WITH A FIXED RMAX ALREADY*/
		  //Rmax=3.0*sqrt(6.0*Dtot*deltat);
		  nc1=ncross[i];
		  nc2=ncross[j];
		  flagsep=get_distancePBC(bases, traj, i, j, deltat, bindrad[mu], ncross,  crosspart, crossint, i1, i2, mu, cross_rxn, i_home, it, Rmax, sep, R1, plist); 
		  //cout <<"trying proteins i and j: "<<i<<' '<<j<<' '<<flagsep<<" Rmax: "<<Rmax<<" R1 "<< R1<<" rsep: "<<sep<<endl;
		  if(flagsep==1){
		    ncrosscom[bases[i].mycomplex]++;
		    ncrosscom[bases[j].mycomplex]++;
		    
		    //evaluate probability
		    
		    ratio=bindrad[mu]/R1;
		    if(sep<0){
		      cout <<"separation <0: "<<sep<<" r1 "<<R1<<" p1: "<<i<<" p2: "<<j<<" Rep: "<<rep<<" it "<<it<<endl;
		      cout <<bases[i].xcom<<' '<<bases[i].ycom<<' '<<bases[i].zcom<<" nfree: "<<bases[i].nfree<<endl;
		      cout <<bases[j].xcom<<' '<<bases[j].ycom<<' '<<bases[j].zcom<<" nfree: "<<bases[j].nfree<<endl;
		      
		      sep=0;
		      ratio=1;
		      R1=bindrad[mu];
		      
		    }
		    
		    //memfact=1.0;
		    //if(bases[i].Dz==0 || bases[j].Dz==0)memfact=0.5;
		    //kdiff=fourpi*Dtot*bindrad[mu]*memfact;
		    
		    //kact=kr[mu];
		    //fact=1.0+kact/kdiff;
		    //alpha=fact*sqrt(Dtot)/bindrad[mu];
		    aexp=sep/sqrt(4.0*Dtot*deltat);
		    bexp=alpha*sqrt(deltat);
		    /*add in error model*/
		    currnorm=1.0;
		    p0_ratio=1.0;
		    
		    //place=i1*Nifaces+i2;
		    //protein i is protype wprot and j
		    wprot2=bases[j].protype;
		    proa=i;
		    pro2=j;
		    if(i>j){
		      proa=j;
		      pro2=i;
		    }
		    
		    /*Need other way to keep track of prevnorm*/
		    for(s=0;s<nprevpart[proa];s++){
		      if(prevlist[proa][s]==pro2){
			p0_ratio=pirr_pfree_ratio_psF(R1, prevsep[proa][s], deltat, Dtot, bindrad[mu], alpha, ps_prev[proa][s], rtol);
			currnorm=prevnorm[proa][s]*p0_ratio;
			ssave=s;
			s=nprevpart[proa];
		      }
		    }
		    ep1=exp(2.0*aexp*bexp+bexp*bexp);
		    if(isinf(ep1)){
		      z.real(0.0);
		      z.imag(aexp+bexp);
		      relerr=0;
		      value=Faddeeva::w(z, relerr);
		      ea2=exp(-aexp*aexp);
		      term2=ea2*real(value);
		    }else{
		      term2=ep1*erfc(aexp+bexp);
		    }
		    probvec1=ratio*kact/(kact+kdiff)*(erfc(aexp)-term2);
		    //cout <<"Ratio: "<<ratio<<" kact: "<<kact<<" kdiff "<<kdiff<<" aexp: "<<aexp<<" bexp: "<<bexp<<" term2: "<<term2<<" ep1: "<<ep1<<" probvec: "<<probvec1<<endl;
		    probvec[i][nc1]=probvec1*currnorm;
		    probvec[j][nc2]=probvec[i][nc1];
		    // cout <<" it: "<<it<<" "<<currnorm<<" probassoc: "<<probvec1<<" ratio for this step: "<<p0_ratio<<" protein a: " <<proa<<'\t'<<pro2<<" currdist: "<<R1<<" prev: "<<prevsep[proa][ssave]<<" previous survival prob: "<<ps_prev[proa][ssave]<<" previous norm factor: "<<prevnorm[proa][ssave]<<endl;
		    s=ncurrpart[proa];
		    currprevsep[proa][s]=R1;
		    currlist[proa][s]=pro2;
		    currprevnorm[proa][s]=currnorm;
		    currps_prev[proa][s]=1.0-probvec1*currnorm;
		    ncurrpart[proa]++;
		  }//close enough to evaluate
		  
		}//A and B bind
		
	      }//end looping over neighboring cells
	    }
	  }//end if i is free
	  
	}//end proteins in cell
      }//end cells
      
      
      for(i=0;i<Ntotalmol;i++){
	//cout <<"protein: "<<i<<" ncross; "<<ncross[i]<<endl;
	if(ncross[i]>0){
	  
	  /*might perform this reaction, depending on k_associate*/
	  
	  //if(i==0)
	  rnum=1.0*rand_gsl();
	  
	  p1=i;
	  ci1=0;//crossint[p1][0];
	  ci2=0;//crossint[p2][0];
	  p2=crosspart[p1][ci1]; //correct if this is who the reaction is for! 
	  
	  /*Now use this instead to decide which of the potential reactions to choose from*/

	  flag=choose_one_reaction(rnum, p1, ncross, crosspart, probvec, ci1, ci2, cross_rxn, crossint, irandmax);
	  /*if flag=0, no reaction, otherwise there was a reaction involving p1 and protein p2 */
	  
	  //cout <<"Try Associate: "<<p1<<' '<<p2<<" i1 "<<crossint[p1][ci1]<<" i2: "<<crossint[p2][ci2]<<" i: "<<i<<" iter: "<<it<<endl;
	  if(flag==1){
	    /*Make sure you choose the first thing that happens*/
	    p2=crosspart[p1][ci1];  
	    rxn1=cross_rxn[p1][ci1];
	    flag2=0;
	    //cout <<"Associate between proteins: p1 "<<p1<<' '<<p2<<" interfaces: "<<crossint[p1][ci1]<<' '<<crossint[p2][ci2]<<" reaction: "<<rxn1<< " pact; "<<probvec[p1][ci1]<<" rnum: "<<rnum<<" iter: "<<it<<endl;
	    
	    
	    /*Then perform association*/
	    /*Correct here for the different types of association*/
	    /*int this case both proteins p1 and p2 have position actually updated*/
	    associate_int_zfirstPBC(p1,  p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist);
	    
	    remove_reaction_all_ncom(p1, p2, ncross, crosspart, probvec, cross_rxn, crossint, ncrosscom, bases);
	    /*Set probability to zero */
	    //remove_one_prob_all(p1, ncross, crosspart, probvec, cross_rxn, crossint);
	    /*Set probability to zero */
	    //remove_one_prob_all(p2, ncross, crosspart, probvec, cross_rxn, crossint);
	    ncross[p1]=-1;
	    ncross[p2]=-1;
	    ncrosscom[bases[p1].mycomplex]=-1;//you won't avoid them, but they should avoid you. you won't move.
	    ncrosscom[bases[p2].mycomplex]=-1;
	    /*this will set ncross[p1] and [p2] to -1 so they stop moving*/
	    
	    
	    	    
	    /*STOPPING ALL THE PROTEINS YOU ARE BOUND TO>*/
	    set_movestat_zero(p1, bases, ind_com, movestat);
	    
	    //movestat[p1]=2;
	    //movestat[p2]=2;
	    
	  }else {
	    /*For this Sweeping version, just move the particle, don't check for overlap until everyone is done
	      store new move in the traj vector so you still know the original position
	    */
	    
	    /*movestat of zero means no traj value is selected.
	      movestat=1 means traj is selected, but particles have not moved
	      movestat=2 means particles have moved
	    */
	    if(movestat[i]==0){
	      k=bases[i].mycomplex;
	      dx=sqrt(2.0*deltat*ind_com[k].Dx)*GaussV();
	      dy=sqrt(2.0*deltat*ind_com[k].Dy)*GaussV();
	      dz=sqrt(2.0*deltat*ind_com[k].Dz)*GaussV();
	      //cout <<"dx: "<<k<<' '<<dx<<' '<<dy<<' '<<dz<<endl;
	      traj[k][0]=dx;
	      traj[k][1]=dy;
	      traj[k][2]=dz;
	    
	      //if other proteins in this complex have already moved, don't need to sample them as well;
	      for(j=0;j<ind_com[k].mysize;j++){
		mp=ind_com[k].plist[j];
		movestat[mp]=1;
	      }

	    }
	    
	    /*Set probability to zero */
	    remove_one_prob_all(p1, ncross, crosspart, probvec, cross_rxn, crossint);
	    
	  }
	  
	}else if(ncross[i]>-1){
	  
	  
	  /*this protein has ncross=0
	    meaning it neither dissociated or tried to associate
	    however, it could have movestat=2 if it was bound to a protein that already moved, multi-protein complex
	    
	    If it had a candidate partner but did not associate, it would also have movestat=2, to avoid overlap,
	    but it would also have ncross>0.
	  */
	  k=bases[i].mycomplex;
	  /*movestat of zero means no traj value is selected.
	    movestat=1 means traj is selected, but particles have not moved
	    movestat=2 means particles have moved
	  */
	  if(movestat[i]==0){
	    
	    dx=sqrt(2.0*deltat*ind_com[k].Dx)*GaussV();
	    dy=sqrt(2.0*deltat*ind_com[k].Dy)*GaussV();
	    dz=sqrt(2.0*deltat*ind_com[k].Dz)*GaussV();
	    
	    traj[k][0]=dx;
	    traj[k][1]=dy;
	    traj[k][2]=dz;
	    
	    //if other proteins in this complex have already moved, don't need to sample them as well;
	    for(j=0;j<ind_com[k].mysize;j++){
	      mp=ind_com[k].plist[j];
	      movestat[mp]=1;
	    }
	    
	  }
	  
	  // 	  if(movestat[i]==0){
// 	    //don't move if already moved, or attached to someone who will 
// 	    //potentially have to resample their position (movestat=1)
// 	    k=bases[i].mycomplex;
// 	    dx=sqrt(2.0*deltat*ind_com[k].Dx)*GaussV();
// 	    dy=sqrt(2.0*deltat*ind_com[k].Dy)*GaussV();
// 	    dz=sqrt(2.0*deltat*ind_com[k].Dz)*GaussV();
// 	    //cout <<"free dx: "<<k<<' '<<dx<<' '<<dy<<' '<<dz<<endl;
// 	    traj[k][0]=dx;
// 	    traj[k][1]=dy;
// 	    traj[k][2]=dz;
	    
// 	    move_protein_posPBC(i, bases, ind_com, traj, movestat, plist); 

	  //}
	  
	}
	
      }
      
      /*Now we have to check for overlap!!!*/
      
      for(i=0;i<Ntotalmol;i++){
	//in trapping case, only B can move, so skip As. 
	/*Any protein that reacted will have ncross=-1 and will have already moved. */
	
	if(ncrosscom[bases[i].mycomplex]>0){
	  if(movestat[i]!=2){
	    /*For any protein that overlapped and did not react, check whether it overlaps with its partners*/
	    //sweep_separation_allPBC(deltat, i, bases, ind_com, ncross, crosspart, crossint, cross_rxn, traj, probvec, plist, movestat, i_home, bindrad, Rlist); 
	    sweep_separation_complexPBC(deltat, i, bases, ind_com, ncross, crosspart, crossint, cross_rxn, traj, probvec, plist, movestat, i_home, bindrad, Rlist); 
	  }
	  
// 	if(ncross[i]>0){
	  
// 	  /*For any protein that overlapped and did not react, check whether it overlaps with its partners*/
// 	  sweep_separation_allPBC(deltat, i, bases, ind_com, ncross, crosspart, crossint, cross_rxn, traj, probvec, plist, movestat, i_home, bindrad, Rlist); 
// 	  /*This protein will not move again because removed from all its partners crossing lists*/
	  
	  
	  
// 	  /*else if you have one partner and it is lower than you, you already tested for overlap
// 	    so don't go through measuring the distance and moving the particles again
// 	  */
// 	}else if(ncross[i]==0){
	}else if(ncrosscom[bases[i].mycomplex]==0){
	  
	  if(movestat[i]!=2){
	    
	    
	    /*the movestat should be 1, if it were zero then
	      the traj would not have been selected. this could not be the case if it had no overlap, and if it did have overlap, it would've selected a new position 
	    */
	    move_protein_posPBC(i, bases, ind_com, traj, movestat, plist); 
	    //  check_complex_rad(i,  bases, ind_com,  plist.xboxl, plist.yboxl, plist.zboxl);
	  }
	  
	}
      }
      Nacurr=Ncopy[0]-2.0*(Ncsave-plist.ntotalcomplex);
      if(it%instwrite==0){
	
	cout<<"timestep: "<<it*deltat<<'\t'<<plist.ntotalcomplex<<'\t'<<Nacurr<<'\t'<<(Nacurr-Aeq)/A0mAeq<<endl;
      }
      curr_time=it*deltat;
      dec_dub=log10(curr_time)+mindecade;
      dec=int(dec_dub);
      decroot=int(mydec[dec]/deltat);//should be 1, 10, 100, etc, same as 10^dec
      if(it%decroot==0){
	delbin=9.0/nperdec[dec];
	targ=(curr_time/mydec[dec]-1.0)/delbin;
	if(int(targ)+1-targ<1E-10){
	  cout <<"PRECISION ISSUE: "<<targ<<" floored: "<<int(targ)<<endl;
	  tdub=int(targ)+2+binsum[dec];//curr_time/mydec is between 1 and 10, so interval sizes are all 9
	}else{
	  tdub=int(targ)+1+binsum[dec];//curr_time/mydec is between 1 and 10, so interval sizes are all 9
	}
	tind=int(tdub);
	//cout <<"it: "<<it<<" dec:" <<dec<<" tind: "<<tind<<" t_double: "<<tdub<<" inside of floor(): "<<targ<<" currtime: "<<curr_time<<" binsum: "<<binsum[dec]<<" floor: "<<floor(targ)<<endl;
	if(dec<0){
	  tind=0;//usually min bin is 1.
	  
	}
	
	//if(it%grfreq==0){
	calc_gr_self(plist, delr, nbins, bases, gr[tind], bindrad[0]);
	calc_gr_self_nnorm(plist, delr, nbins, bases, grnn[tind], bindrad[0]);
      }

      Aavg[it]+=Nacurr;
      // if(it%twrite==0){
// 	write_protein_iface_short(proout, plist, bases,Ncopy,it, wholep, names);
// 	write_status(restart, wholep, bases, plist, Ncopy, it);
// 	write_protein_iface(restart, plist, bases, Ncopy, it, wholep, names);

//       }
      /*Now replace all currsep as prevsep*/
      
      for(i=0;i<Ncopy[0];i++){
	nprevpart[i]=ncurrpart[i];
	for(s=0;s<ncurrpart[i];s++){
	  prevlist[i][s]=currlist[i][s];
	  prevnorm[i][s]=currprevnorm[i][s];
	  ps_prev[i][s]=currps_prev[i][s];
	  prevsep[i][s]=currprevsep[i][s];
	}
      }
      
    }//end iterating over time steps
    if(rep%10==0){
      avfile.open(tname);
      avfile <<0<<'\t'<<rep+1<<endl;
      for(i=0;i<Nit+1;i+=statwrite){
	avfile <<i*deltat<<'\t'<<Aavg[i]<<endl;
      }
      avfile.close();
    }
    
    
  }//end looping over reps  
  avfile.open(tname);
  avfile<<0<<'\t'<<rep<<endl;
  for(i=0;i<Nit+1;i+=statwrite){
    
    //if(i%statwrite==0)
    avfile <<i*deltat<<'\t'<<Aavg[i]<<endl;
  }
  avfile.close();
  
  ofstream grfile("gAself_vs_time.dat");
  grfile<<rho0<<'\t'<<0<<'\t';
  for(j=1;j<nbins;j++)
    grfile<<(j-0.5)*delr+bindrad[0]<<'\t';
  grfile<<endl;
  double Ntotpair=Ntotalmol*(Ntotalmol-1)/2.0;
    
  for(i=0;i<tbins;i++){
    if(tvals[i]<=Maxtime){
      itcurr=int(round(tvals[i]/deltat));
      double Natime=Aavg[itcurr]/Nrep;
      double Npair=Natime*(Natime-1)/2.0;
      gr_volume_norm(gr[i], nbins, delr, rho0, Nrep, bindrad[0], flagdim);
      grfile<<tvals[i]<<'\t';
      for(j=0;j<nbins;j++)
	grfile<<1.0-gr[i][j]/Ntotpair<<'\t';
      grfile<<endl;
    }
  }       
  ofstream grfile2("gAself_npairtime.dat");
  grfile2<<rho0<<'\t'<<0<<'\t';
  for(j=1;j<nbins;j++)
    grfile2<<(j-0.5)*delr+bindrad[0]<<'\t';
  grfile2<<endl;
  
  for(i=0;i<tbins;i++){
    if(tvals[i]<=Maxtime){
      itcurr=int(round(tvals[i]/deltat));
      double Natime=Aavg[itcurr]/Nrep;
      double Npair=Natime*(Natime-1)/2.0;
      
      // gr_volume_norm(gr[i], nbins, delr, rho0, Nrep, bindrad[0]);
      grfile2<<tvals[i]<<'\t' <<1.0-gr[i][0]/Ntotpair<<'\t';//divide survive prob by total pairs, not current density.
      for(j=1;j<nbins;j++)
	grfile2<<1.0-gr[i][j]/Npair<<'\t';
      grfile2<<endl;
    }
  }       
  ofstream grfile3("gAself_instdens.dat");
  grfile3<<rho0<<'\t'<<0<<'\t';
  for(j=1;j<nbins;j++)
    grfile3<<(j-0.5)*delr+bindrad[0]<<'\t';
  grfile3<<endl;

  for(i=0;i<tbins;i++){
    if(tvals[i]<=Maxtime){
      //double Natime=Aavg[i*plist.grfreq]/Nrep;
      //double Npair=Natime*(Natime-1)/2.0;
      gr_volume_norm(grnn[i], nbins, delr, rho0, Nrep, bindrad[0], flagdim);
      grfile3<<tvals[i]<<'\t';
      for(j=0;j<nbins;j++)
	grfile3<<1.0-grnn[i][j]<<'\t';
      grfile3<<endl;
    }
  }       

  
  
  
  
  stop_timer(&totaltime);
  cout <<timer_duration(totaltime)<<" total time "<<endl;  
  /*Write out final result*/
  cout <<"End Main, complete run "<<endl;
  
}//end main
