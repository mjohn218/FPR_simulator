/*


Implement particle-based Reaction diffusion
algorithm with trajectory reweighting.
Allows only two types of particles A and B with
one reversible reaction.

Volume is split into cells to speed calculation.

In this version, A and B avoid overlapping each other even
when they are in the bound state.
Also, the bound particle has the same Diffusion constant as the
A particle.
  
Periodic Boundary conditions at box edges.  

This version also runs multiple trajectories to collect statistics.  
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
#include "utility_calls.h"
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


using namespace std;

struct MD_Timer totaltime;
struct MD_Timer bimoltime;

int main(int argc, char *argv[])
{
  
  /*Define the reactants, and all the species involved in the reactions
   *this includes all the products, including the misbinding products
   */
  int i, idd, j, n, k;
  timeval tim;
  gettimeofday(&tim, 0);
  double t1=tim.tv_sec+tim.tv_usec;
  
	std::vector<gsl_matrix *> contsur;
	std::vector<gsl_matrix *> contnorm;
	std::vector<gsl_matrix *> contpir;
	int MAXALWTBL = 1000;
	cout<<"Maximum number of unique 2D reactions allowed: "<<MAXALWTBL<<endl;
	double TBLID[MAXALWTBL][2]; //[type#][0]:ka, [type#][1]:Dtot
	double Dtemp;
	double ktemp;
	int DDtableindex = 0;
	int tableexistflag = 0;
	int uniquetableindex = 0;

  ofstream AvstTHEO;
  char fname2[100];
  sprintf(fname2,"R1.dat");
  AvstTHEO.open(fname2);

//  struct timespec tim;
//  int retval=0;
//  retval=clock_gettime(CLOCK_REALTIME, &tim);
//  double t1=tim.tv_sec+tim.tv_nsec/10;
//
  int seed=int(t1);
  seed=1407581997;
//  while(seed<0){
//    retval=clock_gettime(CLOCK_REALTIME, &tim);
//    t1=tim.tv_sec+tim.tv_nsec/10;
//    seed=int(t1);
//  }
  
  //  seed=1364043256;
  cout <<"seed: "<<seed<<endl;
  srand_gsl(seed);
  double randmax=pow(2.0, 32);
  double irandmax=1.0/randmax;
  ifstream parmfile(argv[1]);
  Parms plist;
  plist.restart=0;//in case you don't read it in.
  read_parms(parmfile, plist);
  //  write_parms(plist);

  initialize_timer(&totaltime);
  initialize_timer(&bimoltime);
  start_timer(&totaltime);
  
  int Nprotypes=plist.Nprotypes;//total distinct protein types, so 9
  int Nifaces=plist.Nifaces;//this is the number of interfaces

  int Nrxn=plist.Nrxn;
  int Nspecies=plist.Nspecies;//this will include product species


  double cellvol=plist.xboxl*plist.yboxl*plist.zboxl/(1.0*1E9);//sidelength*sidelength*sidelength;
  double V=cellvol; //micrometers^3
  double um_to_L=1E15;
  double avagad=6.022E23;
  //  V=V*avagad/um_to_L;//now V*Molar_Concentration gives a number of molecules needed for Gillespie
  // plist.V=V;

  
  int *Ncopy=new int[Nprotypes];
  /*Read in number of each molecules, and their coordinates*/
  ifstream numfile(argv[2]);
  int Ntotalmol=0;
  plist.Natom=0;
  int ntmp;
  for(i=0;i<Nprotypes;i++){
    numfile >>Ncopy[i];
    Ntotalmol+=Ncopy[i];
  }
  cout <<"Ntotal mols: "<<Ntotalmol<<endl;
  plist.Ntotalmol=Ntotalmol;
  Fullmol *bases=new Fullmol[Ntotalmol];//contains information on each protein in the full system
  Complex *ind_com=new Complex[Ntotalmol];//contains information on each complex
  int *numpartners=new int[Nifaces];//this should account for all free interfaces
  int **Speclist=new int*[Nifaces];
  
  for(i=0;i<Nifaces;i++)
    Speclist[i]=new int[MAXPRTNER];
  

  /*Determine the constituents of the reactions*/
  
  ifstream netfile(argv[3]);
  ifstream protfile(argv[4]);
  ifstream rxnfile(argv[5]);
  char fname[100];
  
  /*ADD NEW REACTIONS SO COMPLEXES HAVE LArGEr BINDING RADII
    one for protein+complex
    one for complex+complex
  */
  
  Protein *wholep=new Protein[Nprotypes];
  int *p_home=new int[Nifaces];//this reverses and tells you what protein a given interface belongs to
  int *i_home=new int[Nspecies];//for both free and bound states, what index are you on the protein's list
  double *bindrad=new double[Nrxn+2];//binding or unbinding radius for each reaction
  int *Ncoup=new int[Nrxn+2];//list of reactions coupled to this one
  int **mycoupled=new int*[Nrxn+2];
  for(i=0;i<Nrxn+2;i++)
    mycoupled[i]=new int[MAXRXN];
  /*The number of reactions is fixed and all the same reactions are possible in each spatial cell*/
  double *kr=new double[Nrxn+2]; //reaction rate (with dimensions)
  
  int **Rlist=new int*[Nrxn+2]; //the identity of the species in the reaction
  int *Npart=new int[Nrxn+2]; //The number of participant species in a reaction
  int **Del=new int*[Nrxn+2]; //The coeffiecients of the participants in the reaction


  int maxrctant=5;
  for(i=0;i<Nrxn+2;i++){
    Rlist[i]=new int[maxrctant];
    Del[i]=new int[maxrctant];
  }
  int *Nmyrxn=new int[Nspecies];
  int **myrxn=new int*[Nspecies];
  for(i=0;i<Nspecies;i++)
    myrxn[i]=new int[MAXRXN];
  int *cntrxn=new int[3];//nfree, nbnd, nmut
  int *freelist=new int[Nspecies];
  int *bndlist=new int[Nspecies];
  int *zlist=new int[Nspecies];
  

  read_protlist(Nprotypes, wholep, Nifaces, p_home, protfile, i_home);
  cout <<"read network: "<<endl;

  
  read_network(plist, numpartners,Speclist, netfile);
  
  /*force it to have only a single partner, even though it has more than one binding interface partner*/
  int Ntotsite=0;
  for(i=0;i<Nprotypes;i++){  
    ntmp=wholep[i].ninterface+1;
    plist.Natom+=Ncopy[i]*ntmp;
    Ntotsite+=Ncopy[i]*wholep[i].ninterface;
  }
  cout <<"N atoms: "<<plist.Natom<<" Nsites, not COM: "<<Ntotsite<<endl;
  int t=0;
  
  cout <<"read reactions "<<endl;
  
  int *rxtype=new int[plist.Nrxn];
  double *Kd=new double[plist.Nrxn];
  read_reactions(rxnfile, plist, Rlist, Ncoup, mycoupled, Nmyrxn, myrxn, bindrad, kr, cntrxn, freelist, bndlist, zlist, rxtype, Kd);
  
  /*Change rates to ka and kb, rather than kon and koff
    unless you are reading in values of ka and kb
  */
  
  //  update_rates(wholep, plist, Rlist,  bindrad, kr, rxtype, Kd, p_home);

  /*List reactant species, check reactants against the interface network*/
  cout <<"free to bind: "<<endl;
  renumber_list(cntrxn[0], freelist);
  cout <<"ready to unbind: "<<endl;
  renumber_list(cntrxn[1], bndlist);
  cout <<"free to mutate: "<<endl;
  if(cntrxn[2]>0)
    renumber_list(cntrxn[2], zlist);
  cout <<"Check reactions ! "<<endl;
  check_reactions(plist, bases, numpartners, Speclist, Nmyrxn, Rlist, myrxn, Nspecies);
  
  /*ADD IN EXTRA REACTIONS*/
  Rlist[2][0]=0;
  Rlist[2][1]=0;
  Rlist[3][0]=0;
  Rlist[3][1]=0;
  double sigma=bindrad[0];
  bindrad[2]=sigma+sigma*0.5;
  bindrad[3]=2.0*sigma;
  
  ifstream startfile(argv[6]);
  cout <<"now set protein status: "<<endl;
  set_status(startfile,wholep,bases, plist, Ncopy);
  /*Read in the coordinates*/
  cout <<"now read in coordinates: "<<endl;
  string *names=new string[Nprotypes];
  ifstream crdfile(argv[7]);
  read_coords(crdfile, Nprotypes, wholep, bases, ind_com, Ncopy, names);
  
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
  double delr=0.1;
  int nbins=20;
  cout <<"nbins: "<<nbins<<" Ntotal mol: "<<Ntotalmol<<" max g(r): "<<nbins*delr<<endl; 
  double **gr=new double*[Nprotypes*Nprotypes];

  for(i=0;i<Nprotypes*Nprotypes;i++)
    gr[i]=new double[nbins];

  //  ofstream gaafile("gaaW.dat");
  ofstream restart("restartsW.out");
  //ofstream gabfile("gab.dat");
  //ofstream gbbfile("gbb.dat");
  


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
  int statwrite=plist.statwrite;
  int twrite=plist.configwrite;
  int outwrite=statwrite*100;
  sprintf(fname,"ComplexW_COM_Np%d_Ni%d.xyz",plist.Nprotypes, Nifaces); 
  ofstream compout(fname);
  sprintf(fname,"ProteinW_COM_Np%d_Ni%d.xyz",plist.Nprotypes, Nifaces); 
  ofstream proout(fname);
  
  
  double Dtot=wholep[0].Dx+wholep[1].Dx;//This is for Dtot between protein 0 and protein 1
  double deltat=plist.deltat;
  double Rmax1=3.5*sqrt(4.0*Dtot*deltat)+bindrad[0];//bindrad for reaction 0
  
  /*cells*/
	int Nx = int(floor(plist.xboxl/Rmax1));//10; //int(ceil(box_x/cellx));
	int Ny = int(floor(plist.yboxl/Rmax1));//10; //int(ceil(box_y/celly));
	int Nz = int(floor(plist.zboxl/Rmax1));//4; //int(ceil(box_z/cellz));
	int Ncell = Nx * Ny * Nz;
  double cellx=box_x/(Nx*1.0);
  double celly=box_y/(Ny*1.0);
  double cellz=box_z/(Nz*1.0);
  if(cellx<Rmax1){
    cout <<"CELL SIZE IS TOO SMALL "<<endl;
    exit(1);
  }
  int c;
  cout <<"Nx: "<<Nx<<" Ny "<<Ny<<" Nz "<<Nz<<" Ncell "<<Ncell<<endl;
  int grfreq=plist.grwrite;
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
  int *nassoclist=new int[MAXPERBIN];
  int nassoc;
  int ndiss;
  int mybin;
  int mybinind;

  cell_neighbor_listPBC(Nx, Ny, Nz, maxnbor, nbor, nborrev);

  cout <<"N cell pairs (max is with PBC): "<<Ncell*maxnbor<<endl;
  
  

  cout <<"deltat: "<<plist.deltat<<endl;

  int it;
  double currnorm, pnormval, pirrval;
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
  
//  write_complex(compout, plist, ind_com, 0);
//  write_protein_iface_short(proout, plist, bases,Ncopy,0, wholep, names);

  int amol,df; 
  double us_to_s=1E-6;

  //  ofstream statfile("Amols.out");
  double **traj=new double *[Ntotalmol];
  for(i=0;i<Ntotalmol;i++)
    traj[i]=new double[3];
  

  double kpi=4.0*M_PI*bindrad[0];
  double fourpi=4.0*M_PI;
  
  double kact;
  double kdiff=4.0*M_PI*Dtot*bindrad[0];
  double fact=1.0+kr[0]/kdiff;
  
  double alpha=fact*sqrt(Dtot)/bindrad[0];
  double cof=kr[0]/(kr[0]+kdiff);
  /*Update all the rates to be ka and kb, rather than kon and koff*/

  cout <<"activation rate in nm^3/us: "<<kr[0]<<" kd: "<<kr[1]<<" Kd, uM: "<<kr[1]/kr[0]/1E6/6.022E23*1E24*1E6<<endl;
  cout <<"Dtot: "<<Dtot<<endl;
  //  check_bndlist(plist.ntotalcomplex, ind_com, bases);
  double R1;

  
  double r0, passoc;
  double Dtot2D = Dtot;
  double Rmax2D=3.0*sqrt(4.0*Dtot2D*deltat)+bindrad[0];
  double Rmax=Rmax2D;
  
  int veclen;

//  double RstepSize = sqrt(Dtot2D*deltat)/50;//0.01;//(Rmax2D-bindrad[0]) / rnumsample;
//
//  veclen = sizelookup(bindrad[0], Rmax2D, RstepSize);
//  gsl_matrix * msur = gsl_matrix_alloc(2, veclen);
//  TBLsur(msur, veclen, bindrad[0], Dtot2D, kr[0], deltat, Rmax2D, RstepSize);
//  cout<<"TBLsur made successfully"<<endl;
//  gsl_matrix * mnorm = gsl_matrix_alloc(2, veclen);
//  TBLnorm(mnorm, veclen, bindrad[0], Dtot2D, kr[0], deltat, Rmax2D, RstepSize);
//  cout<<"TBLnorm made successfully"<<endl;
//  gsl_matrix * mpir = gsl_matrix_alloc(veclen, veclen);
//  TBLpirr(mpir, veclen, bindrad[0], Dtot2D, kr[0], deltat, Rmax2D, RstepSize, Rmax2D);
//  cout<<"TBLpirr made successfully"<<endl;
//
//  FILE * f1 = fopen ("testmsur.dat", "w");
//  gsl_matrix_fprintf (f1, msur,"%f");
//  fclose (f1);
//
//  FILE * f2 = fopen ("testmpir.dat", "w");
//  gsl_matrix_fprintf (f2, mpir,"%f");
//  fclose (f2);
//
//  FILE * f3 = fopen ("testmnorm.dat", "w");
//  gsl_matrix_fprintf (f3, mnorm,"%f");
//  fclose (f3);

//  cout <<" Passoc_vs_separation "<<endl;
//  passoc=survive_irr(bindrad[0], deltat, Dtot, bindrad[0], alpha, cof);
//  cout <<bindrad[0]<<'\t'<<passoc<<endl;
//  int loop=0;
//  i=0;
//  while(loop==0){
//    r0=(i+0.5)*0.001+bindrad[0];
//    passoc=survive_irr(r0, deltat, Dtot, bindrad[0], alpha, cof);
//    cout <<r0<<'\t'<<passoc<<endl;
//    if(passoc<1E-13){
//
//      loop=1;
//      cout<<"passoc: "<<passoc<<" i: "<<i<<" r0: "<<r0<<endl;
//
//    }
//    i++;
//  }
  

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
  double probvec1, p0_ratio;


  /*For multi protein systems this will have to be calculated on-the-fly, because
    Dtot will change.*/

  cout <<"Rmax1: "<<Rmax1<<" Dtot: "<<Dtot<<" Rmaxfinal: "<<Rmax<<" cell size length: "<<cellx<<endl;
  cout <<"IN GET DISTANCE, COMPARING RMAX TO THE ACTUAL DISTANCE R1, NOT JUST THE SEPARATION! "<<endl;
  cout <<"error model "<<endl;
  cout <<"distance is separation beyond bind rad! "<<endl;
  
  /*These only apply for two particle systems A and B*/
  double Vnm3=plist.xboxl*plist.yboxl;//*plist.zboxl;
  double A0=Ncopy[0];
  double B0=Ncopy[1];
  double Keq=kr[0]*1E6/kr[1];
  double Bcof=(B0/Vnm3-A0/Vnm3+1/Keq);
  double Aeq=Vnm3*(-Bcof/2.0+sqrt(Bcof*Bcof+4.0*A0/Vnm3/Keq)/2.0 );
  double A0mAeq=A0-Aeq;
  double Nacurr;
  cout <<"Aequil: "<<Aeq<<" A0; "<<A0<<" B0: "<<B0<<" V: "<<Vnm3<<" Keq "<<Keq<< " nm^3 "<<endl;
  

  
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
  
  
  
  double maxrad=150;
  int radbins=3000;
  double delrad=maxrad/(1.0*radbins);
  int ind_rad;
  double rad2, rad;
  
  int Ncsave=plist.ntotalcomplex;
  double x0;
  int proa, pro2;
  double currx, curry, currz;
  char tname[100];
  sprintf(tname,"A_vs_time.dat");
  ofstream avfile;
  avfile.open(tname);
  int xbin, ybin, zbin;
  int Nrep=plist.grwrite;
//   if (Ncopy[1]==1){
// 	Nrep = 1000;
//   }else if (Ncopy[1]==5){
// 	Nrep = 500;
//   }else if (Ncopy[1]==10){
// 	Nrep = 250;
//   }else if (Ncopy[1]==15){
// 	Nrep = 100;
//   }
  int zeroB=-1;//particles can go anywhere
  int rep=0;
  int totpercell=0;
  double *Aavg=new double[Nit+1];
  for(i=0;i<Nit+1;i++)
    Aavg[i]=0.0;
  int goflag=0;
  
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
    generate_initial_crds_ABavoid( plist, bases, Ncopy, ind_com, bindrad, zeroB);
    write_complex(compout, plist, ind_com, 0);
    write_protein_iface_short(proout, plist, bases,Ncopy,0, wholep, names);
    
    get_bin2( plist, bases, cellx,  celly,cellz,  Nx,  Ny,  Nz, binlist, npb, MAXPERBIN, Ncell, ind_com);  
    
    for(i=0;i<Ntotalmol;i++){
      bases[i].nfree=1;
      bases[i].nbnd=0;
      
    }
    it=0;
    Nacurr=Ncopy[0]-(Ncsave-plist.ntotalcomplex);
//    cout<<"timestep us: "<<it*deltat<<'\t'<<plist.ntotalcomplex<<'\t'<<Nacurr<<'\t'<<(Nacurr-Aeq)/A0mAeq<<endl;
    Aavg[it]+=Nacurr;
    
    
    /***************************/
    /*Begin RD simulation*/
    for(it=1;it<Nit+1;it++){
      
      for(i=0;i<Ntotalmol;i++){
	ncross[i]=0;
	ncrosscom[i]=0;
	movestat[i]=0;
      }
      ndiss=0;
      nassoc=0;
      /*Test dissociation!! */
      for(i=0;i<Ncopy[0];i++){
	/*Test this proteins unimolecular interactions*/
	
	for(j=0;j<bases[i].nbnd;j++){
	  icom=bases[i].bndlist[j];//index of the bound species
	  
	  
	  mu=1;//myrxn[icom][n];
	  
	  
	  i1=Rlist[mu][1];
	  i2=Rlist[mu][2];
	  iind=i_home[i1];
	  if(p_home[i1]!=bases[i].protype){
	    i1=Rlist[mu][2];
	    i2=Rlist[mu][1];
	    iind=i_home[i1];
	  }
	  rate=kr[mu];
	  ppart=bases[i].partner[iind];
	  
	  //then this is the first protein in the reaction, so do it
	  //ensures we do not try to dissociate this complex twice
	  
	  
	  prob=1-exp(-rate*deltat*us_to_s);
	  rnum=1.0*rand_gsl();
	  if(prob>rnum){
	    rnum2=rnum+rand_gsl()*irandmax;
	    if(prob>rnum2){
	      //cout <<"Dissociate at iter: "<<it<<" protein: "<<i<<" partner: "<<ppart<<" randomnum: "<<rnum2<<endl;
	      /*Perform this dissociation reaction*/
		      
	      cancel=break_sigmaPBC(i, mu, j, bases, Rlist, i_home, ind_com, plist, bindrad,  ppart, i1, i2, p_home, myrxn);
	      //cancel=break_complexPBC(i, mu, j, bases, Rlist, i_home, ind_com, plist, bindrad,  ppart, i1, i2, p_home, myrxn);
	      
	      j--;//replaced this reaction, so stay on this one
	      disslist[ndiss]=i;
	      disslist[ndiss+1]=ppart;
	      ndiss+=2;
	      
	      //ncross[i]=-1;
	      //ncross[ppart]=-1;
	      
	      /*Need to make sure that it hasn't been moved to overlap with another particle*/
	      movestat[i]=2;
	      movestat[ppart]=2;
	      
	      
	    }
	  }	      
	  
	}//finished trying out unimolecular reactions for this proteins
      }
      
      get_bin2( plist, bases, cellx,  celly,cellz,  Nx,  Ny,  Nz, binlist, npb, MAXPERBIN, Ncell, ind_com);

      
      /*remove dissociated proteins from the binlist*/
      for(i=0;i<ndiss;i++){
 	p1=disslist[i];
 	k=bases[p1].mycomplex;
	/*Keep these proteins in overlap lists, but don't let them move again!*/
 	ind_com[k].Dx=0;
 	ind_com[k].Dy=0;
 	ind_com[k].Dz=0;
	
	traj[k][0]=0;
	traj[k][1]=0;
	traj[k][2]=0;
	//cout <<"dissociated coords of p1: "<<p1<<" current rnum: "<<rnum<<endl;
	//write_crds(bases,p1);
	// 	// mybin=bases[p1].mybin;
	// // 	mybinind=bases[p1].mybinind;
// // 	//cout <<"ndiss: "<<ndiss<<" i: "<<i<<" p1: "<<p1<<" mybin: "<<mybin<<" nperbin: "<<npb[mybin]<<" mybinind: "<<mybinind<<endl;
// // 	//cout <<"replace: "<<binlist[mybin*MAXPERBIN+mybinind]<<" with: "<<endl;
// // 	cout <<binlist[mybin*MAXPERBIN+npb[mybin]-1]<<endl;
// // 	pnew=binlist[mybin*MAXPERBIN+npb[mybin]-1];
// // 	bases[pnew].mybinind=mybinind;
// // 	binlist[mybin*MAXPERBIN+mybinind]=pnew;
// 	//npb[mybin]-=1;
      }	    
      // cout<<"curr pos 516, 1245, iter: "<<it<<endl;
      //       write_crds(bases, 516);
       
      
      
      for(i=0;i<Ncopy[0];i++){
	ncurrpart[i]=0;
      }
      for(c=0;c<Ncell;c++){
	for(pp=0;pp<npb[c];pp++){
	  i=binlist[c*MAXPERBIN+pp];
	  
	  /*MAKE SURE YOU LLOOP OVER PROTEINS IN YOUR OWN CELL AS WELL */
	  
	  /*Test bimolecular reactions!*/
	  nfree=bases[i].nfree;
	  wprot=bases[i].protype;
	  
	  //	  if(nfree>0){
	  //if(ncross[i]>-1){
	  /*If this protein just dissociated, don't check it for crossing*/
	  /*Now test distance for bimolecular reactions*/
	  //	    for(j=i+1;j<Ntotalmol;j++){
	  //  for(j=Ncopy[0];j<Ntotalmol;j++){
	  for(qq=pp+1;qq<npb[c];qq++){
	    j=binlist[c*MAXPERBIN+qq];
	    /*test to see if they interact.*/
	    mu=0;
	    goflag=0;
	    //  if(bases[j].nfree>0){
	    // if(bases[i].nfree!=1 || bases[j].nfree!=1){
// 	      goflag=1;//test overlap of all particles to bound states.
// 	      if(bases[i].nfree!=1 &&bases[j].nfree!=1)
// 		mu=3;//both complex
// 	      else
// 		mu=2;//one complex
	    //}else 
	    if(bases[j].protype!=bases[i].protype){
	      goflag=1;//both are free, they have to be different to bind.
	    }


	    if(goflag==1){
	      // if(ncross[j]>-1){
	      i1=bases[i].protype;
	      i2=bases[j].protype;
	      
	      
	      if(bases[i].nbnd==1 && bases[i].partner[0]==j){
		/*don't test for overlap, they are bound!!*/
	      }else{
		
		/*Here now we evaluate the probability of binding*/
		
		Dtot=ind_com[bases[i].mycomplex].Dx+ind_com[bases[j].mycomplex].Dx;
		/*CAME UP WITH A FIXED RMAX ALREADY*/
		//Rmax=3.0*sqrt(6.0*Dtot*deltat)+bindrad[mu];
		nc1=ncross[i];
		nc2=ncross[j];
		//Rmax=3.0*sqrt(4.0*deltat*Dtot)+bindrad[mu];
		flagsep=get_distancePBC(bases, traj, i, j, deltat, bindrad[mu], ncross,  crosspart, crossint, i1, i2, mu, cross_rxn, i_home, it, Rmax, sep, R1, plist);
//		AvstTHEO<<R1<<endl;
		
		if(flagsep==1){
		  ncrosscom[bases[i].mycomplex]++;
		  ncrosscom[bases[j].mycomplex]++;
		  // if(sep<0){
		  //cout <<"separation <0: "<<sep<<" r1 "<<R1<<" p1: "<<i<<" p2: "<<j<<" Rep: "<<rep<<" it "<<it<<" nfree?: "<<bases[i].nfree<<'\t'<<bases[j].nfree<<'\t'<<movestat[i]<<'\t'<<movestat[j]<<" rmax: "<<Rmax<<" bindrad: "<<bindrad[mu]<<endl;
			
		  if(bases[i].nfree==1 && bases[j].nfree==1){
		    if(movestat[i]!=2 &&movestat[j]!=2){

		    	//Generate 2D tables unless they were not before
				Dtemp = ind_com[bases[i].mycomplex].Dx+ind_com[bases[j].mycomplex].Dx;
				ktemp = kr[mu]/2/bindrad[mu];

				for (idd = 0; idd < DDtableindex; idd++) {
					if(TBLID[idd][0]==ktemp && TBLID[idd][1]==Dtemp){
						tableexistflag=1;
						uniquetableindex = idd;
						break;
					}
				}

				if(tableexistflag==0){
					TBLID[DDtableindex][0]=ktemp;
					TBLID[DDtableindex][1]=Dtemp;
					veclen = sizelookup(bindrad[mu], Dtemp, deltat, Rmax);
					contsur.push_back(gsl_matrix_alloc(2, veclen));
					contnorm.push_back(gsl_matrix_alloc(2, veclen));
					contpir.push_back(gsl_matrix_alloc(veclen, veclen));
					DDmatrixcreate(contsur[DDtableindex], contnorm[DDtableindex], contpir[DDtableindex], bindrad[mu], Dtemp, ktemp, deltat, Rmax);
					uniquetableindex = DDtableindex;
					DDtableindex +=1;
					if(DDtableindex==MAXALWTBL){
						cout<<"You have hit the maximum number of unique 2D reactions allowed: "<<MAXALWTBL<<endl;
						cout<<"terminating...."<<endl;
						exit(1);
					}
				}
				tableexistflag=0; //reset

		      //evaluate probability
		      ratio=bindrad[mu]/R1;
		      if(sep<0){
			cout <<"separation <0: "<<sep<<" r1 "<<R1<<" p1: "<<i<<" p2: "<<j<<" Rep: "<<rep<<" it "<<it<<" nfree?: "<<bases[i].nfree<<'\t'<<bases[j].nfree<<'\t'<<movestat[i]<<'\t'<<movestat[j]<<" rmax: "<<Rmax<<" bindrad: "<<bindrad[mu]<<endl;
			cout <<"crds: "<<i<<' '<<bases[i].xcom<<' '<<bases[i].ycom<<' '<<bases[i].zcom<<endl;
			cout <<"crds: "<<j<<' '<<bases[j].xcom<<' '<<bases[j].ycom<<' '<<bases[j].zcom<<endl;
			sep=0;
			ratio=1;
			R1=bindrad[mu];
			
		      }
		  
		      
		      kdiff=fourpi*Dtot*bindrad[mu];
		      
		      kact=kr[mu];
		      fact=1.0+kact/kdiff;
		      alpha=fact*sqrt(Dtot)/bindrad[mu];
		      aexp=sep/sqrt(4.0*Dtot*deltat);
		      bexp=alpha*sqrt(deltat);
		      /*Check whether this pair was previously in reaction zone and 
			  therefore if reweighting should apply.
			  Only need to store reweighting values for one protein in the pair,
			  but for each pair need to know the protein partner and interface partner
			*/
		      currnorm=1.0;
		      p0_ratio=1.0;
		      
		      proa=i;
		      pro2=j;
		      if(i>j){
			proa=j;
			pro2=i;
		      }
		      /*Need other way to keep track of prevnorm*/
		      for(s=0;s<nprevpart[proa];s++){
			if(prevlist[proa][s]==pro2){
//			  p0_ratio=pirr_pfree_ratio_ps(R1, prevsep[proa][s], deltat, Dtot, bindrad[mu], alpha, ps_prev[proa][s], rtol);
			  p0_ratio = DDpirr_pfree_ratio_ps(contpir[uniquetableindex], contsur[uniquetableindex], contnorm[uniquetableindex], R1, Dtot, deltat, prevsep[proa][s], ps_prev[proa][s], rtol, bindrad[mu]);

//			  p0_ratio = DDpirr_pfree_ratio_ps(mpir, msur, mnorm, R1,  Dtot,  deltat, prevsep[proa][s],  ps_prev[proa][s],  rtol,  bindrad[mu], RstepSize); //pirr_pfree_ratio_ps(R1, Dtot, deltat, prevsep[proa][s], pirrval, pnormval, ps_prev[proa][s], rtol);pnormval = pnorm(mnorm, RstepSize, prevsep[proa][s], bindrad[mu]);
			  //pirrval = pirr(mpir, msur, RstepSize, R1, prevsep[proa][s], bindrad[mu]);
			  //p0_ratio = pirr_pfree_ratio_ps(R1, Dtot, deltat, prevsep[proa][s], pirrval, pnormval, ps_prev[proa][s], rtol);
			  currnorm=prevnorm[proa][s]*p0_ratio;
			  
			  ssave=s;
			  s=nprevpart[proa];
			}
		      }
//		      probvec1=ratio*kact/(kact+kdiff)*(erfc(aexp)-exp(2.0*aexp*bexp+bexp*bexp)*erfc(aexp+bexp));
//		      probvec1 = DDpsur(msur, Dtot, deltat, R1, bindrad[mu], RstepSize);//probvec1 = psur(msur, RstepSize, R1, bindrad[mu]);
		      probvec1 = DDpsur(contsur[uniquetableindex], Dtot, deltat, R1, bindrad[mu]);

		      probvec[i][nc1]=probvec1*currnorm;
		      probvec[j][nc2]=probvec[i][nc1];
		      
		      s=ncurrpart[proa];
		      currprevsep[proa][s]=R1;
		      
		      currlist[proa][s]=pro2;
		      currprevnorm[proa][s]=currnorm;
		      currps_prev[proa][s]=1.0-probvec1*currnorm;
		      ncurrpart[proa]++;
		    }else{
		      /*If one of the proteins just dissociated, you can still cross, but you can't reassociate*/
		      probvec[i][nc1]=0;
		      probvec[j][nc2]=0;
		    }
		  }else{
		    probvec[i][nc1]=0;
		    probvec[j][nc2]=0;
		  }
		}//close enough to evaluate
	      }
	    }
	  }
	  
	  for(hh=0;hh<maxnbor;hh++){
	    nb=nbor[c*maxnbor+hh];
	    for(qq=0;qq<npb[nb];qq++){
	      j=binlist[nb*MAXPERBIN+qq];
	      /*test to see if these two proteins react*/
	      //if(bases[j].nfree>0){
	      goflag=0;
	      mu=0;
	    //  if(bases[j].nfree>0){
	      // if(bases[i].nfree!=1 || bases[j].nfree!=1){
// 		goflag=1;//test overlap of all particles to bound states.
// 		if(bases[i].nfree!=1 &&bases[j].nfree!=1)
// 		  mu=3;//both complex
// 		else
// 		  mu=2;//one complex
// 	      }else 
	      if(bases[j].protype!=bases[i].protype){
		goflag=1;//both are free, they have to be different to bind.
	      }
	      
	      
	      if(goflag==1){
		
		
		i1=bases[i].protype;
		i2=bases[j].protype;
		
		if(bases[i].nbnd==1 && bases[i].partner[0]==j){
		  /*don't test for overlap, they are bound!!*/
		}else{
		  
		  /*Here now we evaluate the probability of binding*/
		
		  Dtot=ind_com[bases[i].mycomplex].Dx+ind_com[bases[j].mycomplex].Dx;
		  /*CAME UP WITH A FIXED RMAX ALREADY*/
		  //Rmax=3.0*sqrt(4.0*deltat*Dtot)+bindrad[mu];
		  nc1=ncross[i];
		  nc2=ncross[j];
		  flagsep=get_distancePBC(bases, traj, i, j, deltat, bindrad[mu], ncross,  crosspart, crossint, i1, i2, mu, cross_rxn, i_home, it, Rmax, sep, R1, plist);
//		  AvstTHEO<<R1<<endl;
		  if(flagsep==1){
		    ncrosscom[bases[i].mycomplex]++;
		    ncrosscom[bases[j].mycomplex]++;
		  
		    // if(sep<0){
// 		      cout <<"separation <0: "<<sep<<" r1 "<<R1<<" p1: "<<i<<" p2: "<<j<<" Rep: "<<rep<<" it "<<it<<" nfree?: "<<bases[i].nfree<<'\t'<<bases[j].nfree<<'\t'<<movestat[i]<<'\t'<<movestat[j]<<" Rmax: "<<Rmax<<" Bindrad: "<<bindrad[mu]<<endl;
			
		    if(bases[i].nfree==1 && bases[j].nfree==1){
		      if(movestat[i]!=2 &&movestat[j]!=2){

					//Generate 2D tables unless they were not before
					Dtemp = ind_com[bases[i].mycomplex].Dx+ind_com[bases[j].mycomplex].Dx;
					ktemp = kr[mu]/2/bindrad[mu];

					for (idd = 0; idd < DDtableindex; idd++) {
						if(TBLID[idd][0]==ktemp && TBLID[idd][1]==Dtemp){
							tableexistflag=1;
							uniquetableindex = idd;
							break;
						}
					}

					if(tableexistflag==0){
						TBLID[DDtableindex][0]=ktemp;
						TBLID[DDtableindex][1]=Dtemp;
						veclen = sizelookup(bindrad[mu], Dtemp, deltat, Rmax);
						contsur.push_back(gsl_matrix_alloc(2, veclen));
						contnorm.push_back(gsl_matrix_alloc(2, veclen));
						contpir.push_back(gsl_matrix_alloc(veclen, veclen));
						DDmatrixcreate(contsur[DDtableindex], contnorm[DDtableindex], contpir[DDtableindex], bindrad[mu], Dtemp, ktemp, deltat, Rmax);
						uniquetableindex = DDtableindex;
						DDtableindex +=1;
						if(DDtableindex==MAXALWTBL){
							cout<<"You have hit the maximum number of unique 2D reactions allowed: "<<MAXALWTBL<<endl;
							cout<<"terminating...."<<endl;
							exit(1);
						}
					}
					tableexistflag=0; //reset

			//evaluate probability

			ratio=bindrad[mu]/R1;
			if(sep<0){
			  cout <<"separation <0: "<<sep<<" r1 "<<R1<<" p1: "<<i<<" p2: "<<j<<" Rep: "<<rep<<" it "<<it<<" nfree?: "<<bases[i].nfree<<'\t'<<bases[j].nfree<<'\t'<<movestat[i]<<'\t'<<movestat[j]<<" Rmax: "<<Rmax<<" Bindrad: "<<bindrad[mu]<<endl;
			  cout <<"crds: "<<i<<' '<<bases[i].xcom<<' '<<bases[i].ycom<<' '<<bases[i].zcom<<endl;
			  cout <<"crds: "<<j<<' '<<bases[j].xcom<<' '<<bases[j].ycom<<' '<<bases[j].zcom<<endl;
			  
			  sep=0;
			  ratio=1;
			  R1=bindrad[mu];
			  
			}
		  
			
			
			
			kdiff=fourpi*Dtot*bindrad[mu];
			
			kact=kr[mu];
			fact=1.0+kact/kdiff;
			alpha=fact*sqrt(Dtot)/bindrad[mu];
			aexp=sep/sqrt(4.0*Dtot*deltat);
			bexp=alpha*sqrt(deltat);
			
			/*Check whether this pair was previously in reaction zone and 
			  therefore if reweighting should apply.
			  Only need to store reweighting values for one protein in the pair,
			  but for each pair need to know the protein partner and interface partner
			*/
			  
			currnorm=1.0;
			p0_ratio=1.0;
			
			proa=i;
			pro2=j;
			if(i>j){
			  proa=j;
			  pro2=i;
			}
			
			/*Find out what the previous reweighting was for this pair, if they 
			  were stored from previous step (inside reaction zone).
			*/
			for(s=0;s<nprevpart[proa];s++){
			  if(prevlist[proa][s]==pro2){
//			    p0_ratio=pirr_pfree_ratio_ps(R1, prevsep[proa][s], deltat, Dtot, bindrad[mu], alpha, ps_prev[proa][s], rtol);
					p0_ratio = DDpirr_pfree_ratio_ps(contpir[uniquetableindex], contsur[uniquetableindex], contnorm[uniquetableindex], R1, Dtot, deltat, prevsep[proa][s], ps_prev[proa][s], rtol, bindrad[mu]);
				//				pirrval = pirr(mpir, msur, RstepSize, R1, prevsep[proa][s], bindrad[mu]);
				//p0_ratio = pirr_pfree_ratio_ps(R1, Dtot, deltat, prevsep[proa][s], pirrval, pnormval, ps_prev[proa][s], rtol);
				currnorm=prevnorm[proa][s]*p0_ratio;
			    ssave=s;
			    s=nprevpart[proa];
			  }
			}
//			probvec1=ratio*kact/(kact+kdiff)*(erfc(aexp)-exp(2.0*aexp*bexp+bexp*bexp)*erfc(aexp+bexp));
			probvec1 = DDpsur(contsur[uniquetableindex], Dtot, deltat, R1, bindrad[mu]);
			//			probvec1 = psur(msur, RstepSize, R1, bindrad[mu]);
			probvec[i][nc1]=probvec1*currnorm;
			probvec[j][nc2]=probvec[i][nc1];
			s=ncurrpart[proa];
			currprevsep[proa][s]=R1;
			currlist[proa][s]=pro2;
			currprevnorm[proa][s]=currnorm;
			currps_prev[proa][s]=1.0-probvec1*currnorm;
			ncurrpart[proa]++;
		      }else{
			//if one of the proteins just dissociated, don't try to reassociate, but do avoid overlap
			probvec[i][nc1]=0;
			probvec[j][nc2]=0;
		      }
		    }else{
		      probvec[i][nc1]=0;
		      probvec[j][nc2]=0;
		    }
		  }//close enough to evaluate
		}
	      }
	    }
	    
	  }//end looping over neighboring cells
	  
	}//end proteins in cell
      }//end cells
      
      
      for(i=0;i<Ntotalmol;i++){
	
	if(ncross[i]>0){
	  
	  /*might perform this reaction, depending on k_associate*/
	  rnum=1.0*rand_gsl();
	  //cout <<"Molecule: "<<i<< " nross: "<<ncross[i]<<" rnum: "<<rnum<<endl;
	  p1=i;
	  ci1=0;
	  ci2=0;
	  p2=crosspart[p1][ci1]; //correct if this is who the reaction is for! 
	  
	  flag=choose_one_reaction(rnum, p1, ncross, crosspart, probvec, ci1, ci2, cross_rxn, crossint, irandmax);
	  
	  if(flag==1){
	    /*Make sure you choose the first thing that happens*/
	    p2=crosspart[p1][ci1];  
	    rxn1=cross_rxn[p1][ci1];
	    flag2=0;
	    //cout <<"Associate between proteins: p1 "<<p1<<' '<<p2<<" interfaces: "<<crossint[p1][ci1]<<' '<<crossint[p2][ci2]<<" reaction: "<<rxn1<< " pact; "<<probvec[p1][ci1]<<" rnum: "<<rnum<<" iter: "<<it<<endl;
	    
	    
	    /*Then perform association*/
	    associate_zsigmaPBC(p1,  p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist, bindrad[rxn1], ncrosscom);
	    //associate_int_zfirstPBC(p1,  p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist);
	    
	    /*Make product have same diffusion as proteins*/
	    k=bases[p1].mycomplex;
	    ind_com[k].Dx=0;//bases[p1].Dx;
	    ind_com[k].Dy=0;//bases[p1].Dy;
	    ind_com[k].Dz=0;//bases[p1].Dz;
	    nassoclist[nassoc]=p1;
	    traj[k][0]=0;
	    traj[k][1]=0;
	    traj[k][2]=0;
	    nassoc++;
	    
	    
	    
	    /*Set probability to zero */
	    remove_one_prob_all(p1, ncross, crosspart, probvec, cross_rxn, crossint);
	    /*Set probability to zero */
	    remove_one_prob_all(p2, ncross, crosspart, probvec, cross_rxn, crossint);
	    ncross[p1]=-1;
	    ncross[p2]=-1;
	    ncrosscom[bases[p1].mycomplex]=-1;//you won't avoid them, but they should avoid you. you won't move.
	    ncrosscom[bases[p2].mycomplex]=-1;
	    /*this will set ncross[p1] and [p2] to -1 so they stop moving*/
	    //remove_reaction_all(p1, p2, ncross, crosspart, probvec, ci1, ci2, cross_rxn, crossint);
	    
	    
	    /*STOPPING ALL THE PROTEINS YOU ARE BOUND TO>*/
	    set_movestat_zero(p1, bases, ind_com, movestat);
	    
	    
	  }else {
	    /*For this Sweeping version, just move the particle, don't check for overlap until everyone is done
	      store new move in the traj vector so you still know the original position
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
	    
	    	    
	    /*Set probability to zero */
	    remove_one_prob_all(p1, ncross, crosspart, probvec, cross_rxn, crossint);
	    
	  }
	  
	}else if(ncross[i]>-1){
	  
	  
	  /*this protein has ncross=0
	    meaning it neither dissociated or tried to associate
	    however, it could have movestat=2 if it was bound to a protein that already moved, multi-protein complex
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
	    
	    // if(movestat[i]==0){
	    
// 	    k=bases[i].mycomplex;
// 	    dx=sqrt(2.0*deltat*ind_com[k].Dx)*GaussV();
// 	    dy=sqrt(2.0*deltat*ind_com[k].Dy)*GaussV();
// 	    dz=sqrt(2.0*deltat*ind_com[k].Dz)*GaussV();
	    
// 	    traj[k][0]=dx;
// 	    traj[k][1]=dy;
// 	    traj[k][2]=dz;
	    
// 	    move_protein_posPBC(i, bases, ind_com, traj, movestat, plist); 
	    
// 	  }
	  
	}
	
      }
      // k=bases[469].mycomplex;
//        i=bases[1810].mycomplex;
// //       /*Now we have to check for overlap!!!*/
//        cout <<"ncross 469: "<<ncross[469]<<" 1810: "<<ncross[1810]<<" mycomplex:" <<k<<'\t'<<i<<endl;
//        cout <<"movestat 469: "<<movestat[469]<<" 1810: "<<movestat[1810]<<" traj: "<<traj[k][0]<<'\t'<<traj[k][1]<<" traj2: "<<traj[i][0]<<'\t'<<traj[i][1]<<endl;
//        cout<<"curr pos 469, 1810, iter: "<<it<<endl;
//        write_crds(bases, 469);
//        write_crds(bases,1810);

      
      for(i=0;i<Ntotalmol;i++){
	
	/*Any protein that reacted will have ncross=-1 and will have already moved. */
	
	if(ncrosscom[bases[i].mycomplex]>0){
	  if(movestat[i]!=2){
	    /*For any protein that overlapped and did not react, check whether it overlaps with its partners*/
	    //sweep_separation_allPBC(deltat, i, bases, ind_com, ncross, crosspart, crossint, cross_rxn, traj, probvec, plist, movestat, i_home, bindrad, Rlist); 
	    sweep_separation_complexPBC(deltat, i, bases, ind_com, ncross, crosspart, crossint, cross_rxn, traj, probvec, plist, movestat, i_home, bindrad, Rlist); 
	  }
	}else if(ncrosscom[bases[i].mycomplex]==0){
	  
	  if(movestat[i]!=2){
	    
	    /*For proteins with ncross=0, they either moved independently, or their displacements
	      were selected based on the complex they were part of, and they may not yet been moved.
	    */
	    move_protein_posPBC(i, bases, ind_com, traj, movestat, plist); 
	    
	  }
	  
	}
      }
      for(i=0;i<ndiss;i++){
	p1=disslist[i];
	k=bases[p1].mycomplex;
	/*Keep these proteins in overlap lists, but don't let them move again!*/
	ind_com[k].Dx=bases[p1].Dx;
	ind_com[k].Dy=bases[p1].Dy;
	ind_com[k].Dz=bases[p1].Dz;
	
	//write_crds(bases,p1);
	
      }
      //cout <<"End of iter: "<<it<<" dissociated coords of p1: "<<p1<<" current rnum: "<<rnum<<endl;
      for(i=0;i<nassoc;i++){
	p1=nassoclist[i];
	k=bases[p1].mycomplex;
	/*put this to zero temporarily to avoid overlap after association*/
	ind_com[k].Dx=bases[p1].Dx;
	ind_com[k].Dy=bases[p1].Dy;
	ind_com[k].Dz=bases[p1].Dz;
      }
      Nacurr=Ncopy[0]-(Ncsave-plist.ntotalcomplex);
      
      if(it%statwrite==0){
	cout<<"timestep: "<<it*deltat<<'\t'<<plist.ntotalcomplex<<'\t'<<Nacurr<<'\t'<<(Nacurr-Aeq)/A0mAeq<<" rnum: "<<rnum<<endl;
//	cout<<"timestep: "<<it*deltat<<'\t'<<plist.ntotalcomplex<<'\t'<<Nacurr<<'\t'<<(Nacurr-Aeq)/A0mAeq<<endl;
      }
      Aavg[it]+=Nacurr;
			if (it % twrite == 0) {
				write_protein_iface_short(proout, plist, bases, Ncopy, it, wholep, names);
//				write_status(restart, wholep, bases, plist, Ncopy, it);
//				write_protein_iface(restart, plist, bases, Ncopy, it, wholep, names);

			}
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
//      avfile <<0<<'\t'<<rep+1<<endl;
      for(i=0;i<Nit+1;i+=grfreq){
	avfile << std::setprecision(8) <<i*deltat<<'\t'<<Aavg[i]<<endl;
      }
      avfile.close();
    }
    
////    avfile.open(tname);
//	Nacurr=Ncopy[0]-(Ncsave-plist.ntotalcomplex);
//	avfile<<Nacurr<<endl;
////    avfile.close();

  }//end looping over reps  
  avfile.open(tname);
//  avfile<<0<<'\t'<<rep<<endl;
  for(i=0;i<Nit+1;i+=grfreq){

    //if(i%statwrite==0)
    avfile << std::setprecision(8) << i*deltat<<'\t'<<Aavg[i]<<endl;
  }
  avfile.close();
  stop_timer(&totaltime);

  cout <<timer_duration(totaltime)<<" total time "<<endl;  
  /*Write out final result*/
  cout <<"End Main, complete run "<<endl;
  
}//end main
