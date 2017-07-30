/*
single particle reaction diffusion
with trajectory reweighting.

This program specifically propagates an irreversible
reaction.  
Volume is split into cells to speed calculation,
assumes reacting species are all molecules, e.g.
A+A->0 only, no B particles. 

This version adds an interaction potential between the molecules.
  
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
#include "Vrnumer.h"
#include "GF_calls.h"
#include "Faddeeva.hh"
#include "utility_calls.h"
#include "cell_neighbor_lists.h"

#include "numeric_GF.h"

using namespace std;

struct MD_Timer totaltime;
struct MD_Timer bimoltime;

int main(int argc, char *argv[])
{
  
  /*Define the reactants, and all the species involved in the reactions
   *this includes all the products, including the misbinding products
   */
  int i, j, n, k;
  timeval tim;
  gettimeofday(&tim, 0);
  double t1=tim.tv_sec+tim.tv_usec;
  
  int seed=int(t1);
  
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

  int numcells=1;
  //  double sidelength=plist.boxl/1000;//micrometer: boxl is in nm
  double cellvol=plist.xboxl*plist.yboxl*plist.zboxl/(1.0*1E9);//sidelength*sidelength*sidelength;
  double V=cellvol*numcells; //micrometers^3
  double um_to_L=1E15;
  double avagad=6.022E23;
  V=V*avagad/um_to_L;//now V*Molar_Concentration gives a number of molecules  
  plist.V=V;
  //kforward divided by this V gives per second, permolecule
  
  //double X0total=plist.X0total;//0.1 millimolar total concentration

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
  
  Protein *wholep=new Protein[Nprotypes];
  int *p_home=new int[Nifaces];//this reverses and tells you what protein a given interface belongs to
  int *i_home=new int[Nspecies];//for both free and bound states, what index are you on the protein's list
  double *bindrad=new double[Nrxn];//binding or unbinding radius for each reaction
  int *Ncoup=new int[Nrxn];//list of reactions coupled to this one
  int **mycoupled=new int*[Nrxn];
  for(i=0;i<Nrxn;i++)
    mycoupled[i]=new int[MAXRXN];
  /*The number of reactions is fixed and all the same reactions are possible in each spatial cell*/
  double *kr=new double[Nrxn]; //reaction rate (with dimensions)
  
  int **Rlist=new int*[Nrxn]; //the identity of the species in the reaction
  int *Npart=new int[Nrxn]; //The number of participant species in a reaction
  int **Del=new int*[Nrxn]; //The coeffiecients of the participants in the reaction

  /*Do not assume all reactions possible
   *so instead, we have to figure out a way to read them in! 
   */
  int maxrctant=5;
  for(i=0;i<Nrxn;i++){
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
  /*Change rates to ka and kb, rather than kon and koff*/
  
  cout <<"free to bind: "<<endl;
  renumber_list(cntrxn[0], freelist);
  //  cout <<"ready to unbind: "<<endl;
  //renumber_list(cntrxn[1], bndlist);
  //cout <<"free to mutate: "<<endl;
  if(cntrxn[2]>0)
    renumber_list(cntrxn[2], zlist);
  cout <<"Check reactions ! "<<endl;
  check_reactions(plist, bases, numpartners, Speclist, Nmyrxn, Rlist, myrxn, Nspecies);

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
  

  /*get all the diffusion constants for the reactions
   *little d are D/l^2, where l is the length of the box: sidelength
   */
 
  
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
  double delr=0.6;
  int nbins=int(plist.zboxl/2.0/delr);
  cout <<"nbins: "<<nbins<<" Ntotal mol: "<<Ntotalmol<<" max g(r): "<<nbins*delr<<endl; 
  // double **gr=new double*[Nprotypes*Nprotypes];

//   for(i=0;i<Nprotypes*Nprotypes;i++)
//     gr[i]=new double[nbins];

//   double **Vofr=new double *[Ntotalmol];
//   for(i=0;i<Ntotalmol;i++)
//     Vofr[i]=new double[nbins];
//   for(i=0;i<Ntotalmol;i++){
//     for(j=0;j<nbins;j++)
//       Vofr[i][j]=0;
//   }
//   cout<<" vofr: "<<Vofr[0][0]<<endl;
  //  ofstream gaafile("gaaW.dat");
  // ofstream restart("restartsW.out");
  //  ofstream gabfile("gab.dat");
  //ofstream gbbfile("gbb.dat");
  ofstream cvstime("Ncomplex_vs_time.dat");
  
  /*Calculate some equilibrium averages*/

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

  int *movestat=new int[Ntotalmol];
  
  double h;
  //  sprintf(fname, "Correlation%d.out",Nifaces);
  //ofstream corrfile(fname);
  /*The number of species does not include all separate interfaces, because it is for diffusible species*/

  int **inst_pro=new int *[Nprotypes];
  for(i=0;i<Nprotypes;i++)
    inst_pro[i]=new int[MAXIFACE];
  for(i=0;i<Nprotypes;i++){
    inst_pro[i][0]=Ncopy[i];
    for(j=1;j<wholep[i].ninterface+1;j++)
      inst_pro[i][j]=0;
  }
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
  int twrite=plist.configwrite;

  sprintf(fname,"ComplexW_COM_Np%d_Ni%d.xyz",plist.Nprotypes, Nifaces); 
  ofstream compout(fname);
  sprintf(fname,"ProteinW_COM_Np%d_Ni%d.xyz",plist.Nprotypes, Nifaces); 
  ofstream proout(fname);


  /***************************/
  /*Begin RD simulation*/
  cout <<"deltat: "<<plist.deltat<<endl;
  double deltat=plist.deltat;
  int it;
  ifstream restartf;
  plist.ntotalcomplex=Ntotalmol;
  int Nit=int(plist.Nit);
  if(Nit>2.147E9){
    cout <<"ITERATIONS EXCEEDS INTEGER MAXIMUM! Exiting..."<<endl;
    exit(1);
  }
  int s1;
  cout <<"Ntotal complexes: "<<plist.ntotalcomplex<<endl;
  write_complex(compout, plist, ind_com, 0);
  write_protein_iface(proout, plist, bases,Ncopy,0, wholep, names);
  
  int amol,df; 
  double us_to_s=1E-6;
  int statwrite=plist.statwrite;
  //  ofstream statfile("Amols.out");
  double **traj=new double *[Ntotalmol];
  for(i=0;i<Ntotalmol;i++)
    traj[i]=new double[3];
  double **Ftraj=new double *[Ntotalmol];
  for(i=0;i<Ntotalmol;i++)
    Ftraj[i]=new double[3];
  
  double Dtot=wholep[0].Dx*2.0;//Da+Da
  double Rmax1=4.0*sqrt(6.0*Dtot*deltat)+bindrad[0];
  double Rmax=Rmax1;

  /*read in functional fit parameters to numerically calculated pirr(r,deltat|r0)
   */
  double delr0=0.05;
  int ncol=7;
  int r0bins_approx=int((Rmax-bindrad[0])/delr0)+10;
  double **FitMat=new double*[r0bins_approx];
  for(i=0;i<r0bins_approx;i++)
    FitMat[i]=new double[ncol];
  ifstream vrfile(argv[8]);
  Vrnumer limit;
  read_FitMat_VrparmsRD_matfull(FitMat, vrfile, limit, ncol);
  cout <<"Read in Vr parms. "<<endl;

  /*Establish parameters to define survival prob*/
  // double Dmat_ps=limit.Dmat_ps;
//   double limit.kmat_ps=limit.limit.kmat_ps;
  cout <<"Suvival Dfit and kfit parmts: "<<limit.Dmat_ps[0]<<' '<<limit.kmat_ps[0]<<endl;
  cout <<"Suvival Dfit and kfit parmts, second half: "<<limit.Dmat_ps[1]<<' '<<limit.kmat_ps[1]<<endl;
  
  double dudr;
  double Fx, Fy, Fz;
  double rc=limit.rc;
  cout <<"BetaU(r) rc factor: "<<rc<<endl;
  

  double kact;
  double fact;
  double kdiff;//will be kpi*D
  double aexp;
  double alpha;
  double fourpi=4.0*M_PI;
    
  mu=0;
  kdiff=4.0*M_PI*limit.Dmat_ps[0]*bindrad[mu];
  kact=limit.kmat_ps[0];
  fact=1.0+kact/kdiff;
  double alpha1=fact*sqrt(limit.Dmat_ps[0])/bindrad[mu];
  //double bexp=alpha*sqrt(deltat);
  
  double cof1=limit.kmat_ps[0]/(limit.kmat_ps[0]+kdiff);
  
  //region 2
  kdiff=4.0*M_PI*limit.Dmat_ps[1]*bindrad[mu];
  kact=limit.kmat_ps[1];
  fact=1.0+kact/kdiff;
  double alpha2=fact*sqrt(limit.Dmat_ps[1])/bindrad[mu];
  //double bexp=alpha*sqrt(deltat);
  
  double cof2=limit.kmat_ps[1]/(limit.kmat_ps[1]+kdiff);
  

  int myplace;
  /*Update all the rates to be ka and kb, rather than kon and koff*/

  cout <<"activation rate in nm^3/us: "<<kr[0]<<" kd: "<<kr[1]<<" Kd, uM: "<<kr[1]/kr[0]/1E6/6.022E23*1E24*1E6<<endl;
  cout <<"Dtot: "<<Dtot<<" kact: "<<kact<<" alpha1: "<<alpha1<<" cof 1 "<<cof1<<" alpha2: "<<alpha2<<" cof2: "<<cof2<<endl;
  //  check_bndlist(plist.ntotalcomplex, ind_com, bases);
  double R1, R2;
  double r0, passoc;
  r0=bindrad[0];
  passoc=survive_irrF(r0, deltat, limit.Dmat_ps[0], bindrad[0], alpha1, cof1);
    
  cout<<"passoc1: "<<passoc<<" bindrad: "<<bindrad[0]<<endl;
  for(i=0;i<1000;i++){
    r0=(i+0.5)*0.05+bindrad[0];
    passoc=survive_irrF(r0, deltat, Dtot, bindrad[0], alpha1, cof1);
    cout <<r0<<'\t'<<passoc<<endl;
    if(passoc<1E-13){
      i=1000;
      
    }
  }
  cout<<"passoc: "<<passoc<<" i: "<<i<<" r0: "<<r0<<endl;
  
  int Nx=10;//int(ceil(box_x/cellx));
  int Ny=10;//int(ceil(box_y/celly));
  int Nz=10;//int(ceil(box_z/cellz));
  int Ncell=Nx*Ny*Nz;
  double cellx=box_x/(Nx*1.0);
  double celly=box_y/(Ny*1.0);
  double cellz=box_z/(Nz*1.0);
  if(cellx<Rmax1){
    cout <<"CELL SIZE IS TOO SMALL "<<endl;
    exit(1);
  }
  
  cout <<"Nx: "<<Nx<<" Ny "<<Ny<<" Nz "<<Nz<<" Ncell "<<Ncell<<endl;

  int maxnbor=13; //geometry of cube
  //int *Nnbor=new int[Ncell];
  int *nbor=new int[Ncell*maxnbor];
  int *nborrev=new int[Ncell*maxnbor];
  
  int *npb=new int[Ncell];
  int MAXPERBIN=200;//int(Ntotalmol);
  int *binlist=new int[Ncell*MAXPERBIN];
  int MAXDISS=20;//max number of proteins to dissociate in one time step 
  int *disslist=new int[MAXDISS];//holds index of dissociated proteins
  int ndiss;
  int mybin;
  int mybinind;
  int c=0;
  
  cell_neighbor_listPBC(Nx, Ny, Nz, maxnbor, nbor, nborrev);
    
  cout <<"N cell pairs (max is with PBC): "<<Ncell*maxnbor<<endl;
    

  
  // double **prevnorm=new double*[Nifaces*Nifaces];
//   int **previter=new int*[Nifaces*Nifaces];
//   double **ps_prev=new double*[Nifaces*Nifaces];
//   double **prevsep=new double*[Nifaces*Nifaces];
//   int *arrsize=new int[Nifaces*Nifaces];


  int size=Ncopy[0]*Ncopy[0];//For this A+A-> reaction only!!
  double *prevnorm=new double[size];
  int *previter=new int[size];
  double *ps_prev=new double[size];
  double *prevsep=new double[size];
  //  int *arrsize=new int[Nifaces*Nifaces];

  //  int place, 
  int place2;
  double rtol=1E-10;
  int flag;
  int p1, p2;
  int tar, c1;
  double probvec1, p0_ratio, currnorm;
  // for(i=0;i<Nifaces;i++){
//     for(j=0;j<Nifaces;j++){
//       flag=0;
//       place=i*Nifaces+j;
//       for(k=0;k<numpartners[i];k++){
// 	if(Speclist[i][k]==j){
// 	  //then define an array based on the number of each protein copies
// 	  p1=p_home[i];
// 	  p2=p_home[j];
// 	  size=Ncopy[p1]*Ncopy[p2];
// 	  prevnorm[place]=new double[size];
// 	  previter[place]=new int[size];
// 	  ps_prev[place]=new double[size];
// 	  prevsep[place]=new double[size];
// 	  arrsize[place]=size;
// 	  flag=1;
// 	}
//       }
//       if(flag==0){
// 	prevnorm[place]=new double[0];
// 	previter[place]=new int[0];
// 	ps_prev[place]=new double[0];
// 	prevsep[place]=new double[0];
// 	arrsize[place]=0;
//       }
//     }
//   }


  int k1, k2;
//  for(i=0;i<Nifaces*Nifaces;i++){
  for(i=0;i<size;i++){
      prevnorm[i]=1.0;
      previter[i]=0;
      ps_prev[i]=0;
      prevsep[i]=0;
    
  }
  /*For multi protein systems this will have to be calculated on-the-fly, because
    Dtot will change.*/
  //  errormodel( deltat, Dtot,  bindrad[0],  kr[0],delr0,  r0bins, pcorrect, Rmax, pfreea);
  cout <<"Rmax1: "<<Rmax1<<" Dtot: "<<Dtot<<" Rmaxfinal: "<<Rmax<<endl;
  cout <<"IN GET DISTANCE, COMPARING RMAX TO THE ACTUAL DISTANCE R1, NOT JUST THE SEPARATION! "<<endl;
  cout <<"error model "<<endl;
  cout <<"distance is separation beyond bind rad! "<<endl;
  

  double pact;
  double tmpx, tmpy, tmpz;

  int pnew;
  int  c2;
  int ci1, ci2;
  double tremain, tevent;
  int mu_ret;
  int rxn1;
  int go;
  double rate;
  int cancel;

  int nc1, nc2;
  double sep, ratio;



  int flagsep;
  double memfact=1.0;//If you are binding to PIP2, this should be 1/2.
  int pp, qq, nb, hh;
  
  
  
  double maxrad=150;
  int radbins=3000;
  double delrad=maxrad/(1.0*radbins);
  int ind_rad;
  double rad2, rad;

  int Nrow=Ncopy[0]*Ncopy[1];
  cout <<"pairs of interacting particles: "<<Nrow<<endl;
  int row;


  int Nit1=Nit+1;
  int Rlen=radbins+1;

  int Ncsave=plist.ntotalcomplex;
  int Ntotsave=Ntotalmol;
  double x0;
  double *r0value=new double[Nrow];
  double currx, curry, currz;
  /*initial separation between each pair*/
  row=0;
  double NA=0;
  double *avgc=new double[Nit];
  int nc;
  for(i=0;i<Nit;i++)
    avgc[i]=0;

  double small=1E-3;
  std::complex<double> z;
  std::complex<double> value;
  double relerr;
  double term2, ea2, ep1;

  int rep;
  int Nrep=plist.grwrite;
  for(rep=0;rep<Nrep;rep++){
    /*start proteins separated by sigma */
    Ntotalmol=Ntotsave;
    cout <<"Repeat: "<<rep<<" molecules: "<<Ntotalmol<<endl;
    
    //for(i=0;i<Nifaces*Nifaces;i++){
    for(i=0;i<size;i++){
      prevnorm[i]=1.0;
      previter[i]=-1;
      ps_prev[i]=0;
      prevsep[i]=0;
    }
      
      /*generate initial coordinates*/
    for(i=0;i<Ntotalmol;i++){
      bases[i].xcom=plist.xboxl*rand_gsl()-plist.xboxl/2.0;
      bases[i].ycom=plist.yboxl*rand_gsl()-plist.yboxl/2.0;
      bases[i].zcom=plist.zboxl*rand_gsl()-plist.zboxl/2.0;
    }
    int bflag=1;
    int bit=0;
    int noverlap=0;
    double pbindr2=bindrad[0]*bindrad[0];
    double stretch;
    /*first just try resampling positions*/
    
    while(bflag==1 && bit<50 ){
      bit++;
      bflag=0;
      noverlap=0;
      for(i=0;i<Ntotalmol;i++){
	for(j=i+1;j<Ntotalmol;j++){
	  dx=bases[j].xcom-bases[i].xcom;
	  dy=bases[j].ycom-bases[i].ycom;
	  dz=bases[j].zcom-bases[i].zcom;
	  dx-=plist.xboxl*round(dx/plist.xboxl);
	  dy-=plist.yboxl*round(dy/plist.yboxl);
	  dz-=plist.zboxl*round(dz/plist.zboxl);
	  
	  double d2=dx*dx+dy*dy+dz*dz;
	  if(d2<pbindr2){
	    bases[j].xcom=plist.xboxl*rand_gsl()-plist.xboxl/2.0;
	    bases[j].ycom=plist.yboxl*rand_gsl()-plist.yboxl/2.0;
	    bases[j].zcom=plist.zboxl*rand_gsl()-plist.zboxl/2.0;
   
	    	    
	    noverlap++;
	    bflag=1;
	  }
	}
      }
      cout <<"Noverlap: "<<noverlap<<endl;
    }
    /*If that didn't work, push them apart*/
    noverlap=0;
    bflag=1;
    bit=0;
    while(bflag==1 && bit<50 ){
      bit++;
      bflag=0;
      noverlap=0;
      for(i=0;i<Ntotalmol;i++){
	for(j=i+1;j<Ntotalmol;j++){
	  dx=bases[j].xcom-bases[i].xcom;
	  dy=bases[j].ycom-bases[i].ycom;
	  dz=bases[j].zcom-bases[i].zcom;
	  dx-=plist.xboxl*round(dx/plist.xboxl);
	  dy-=plist.yboxl*round(dy/plist.yboxl);
	  dz-=plist.zboxl*round(dz/plist.zboxl);
	  
	  double d2=dx*dx+dy*dy+dz*dz;
	  
	  if(d2<pbindr2){
	    stretch=sqrt(pbindr2/d2)+small*4.0*rand_gsl();
	    bases[j].xcom+=dx*(stretch-1.0);
	    bases[j].ycom+=dy*(stretch-1.0);
	    bases[j].zcom+=dz*(stretch-1.0);
	    	    
	    noverlap++;
	    bflag=1;
	  }
	}
      }
      cout <<"Noverlap: "<<noverlap<<endl;
    }
    // cout <<"Config: "<<rep<<endl;
    // for(i=0;i<Ntotalmol;i++){
    //   cout <<bases[i].xcom<<' '<<bases[i].ycom<<' '<<bases[i].zcom<<endl;
    // }
    plist.ntotalcomplex=Ncsave;
    get_bin2( plist, bases, cellx,  celly,cellz,  Nx,  Ny,  Nz, binlist, npb, MAXPERBIN, Ncell, ind_com);  
    for(i=0;i<Ntotalmol;i++){
      /*Restart with a new configuration of A particles, and make sure they 
	  are a minimum of sigma apart from themselves*/
	
	// bases[i].xcom=savecrds[i*3+0];
// 	bases[i].ycom=savecrds[i*3+1];
// 	bases[i].zcom=savecrds[i*3+2];
	bases[i].x[0]=bases[i].xcom;
	bases[i].y[0]=bases[i].ycom;
	bases[i].z[0]=bases[i].zcom;
	bases[i].freelist[0]=bases[i].protype;//either 0 or 1 is also the interface number
	bases[i].istatus[0]=bases[i].freelist[0];//is this true?
	bases[i].mycomplex=i;
	bases[i].nfree=1;
	bases[i].nbnd=0;
	bases[i].npartner=0;
	ind_com[i].xcom=bases[i].xcom;
	ind_com[i].ycom=bases[i].ycom;
	ind_com[i].zcom=bases[i].zcom;
	ind_com[i].plist[0]=i;
	ind_com[i].mysize=1;
	ind_com[i].Dx=bases[i].Dx;
	ind_com[i].Dy=bases[i].Dy;
	ind_com[i].Dz=bases[i].Dz;
	
      }
      

    
      for(it=1;it<Nit+1;it++)
	{
	  
	  for(i=0;i<Ntotalmol;i++){
	    ncross[i]=0;
	    movestat[i]=0;
	    k=bases[i].mycomplex;
	    Ftraj[k][0]=0;
	    Ftraj[k][1]=0;
	    Ftraj[k][2]=0;
	    
	  }
	  
	  for(i=0;i<Ntotalmol;i++){
	    for(j=i+1;j<Ntotalmol;j++){
	      // get_bin2( plist, bases, cellx,  celly,cellz,  Nx,  Ny,  Nz, binlist, npb, MAXPERBIN, Ncell, ind_com);  
// 	  for(c=0;c<Ncell;c++){
// 	    for(pp=0;pp<npb[c];pp++){
// 	      i=binlist[c*MAXPERBIN+pp];
	      
	      
	      if(bases[i].nfree==1){
	    
	      /*MAKE SURE YOU LLOOP OVER PROTEINS IN YOUR OWN CELL AS WELL */
	      
	      /*Test bimolecular reactions!*/
	      
	      
	      //for(qq=pp+1;qq<npb[c];qq++){
		//j=binlist[c*MAXPERBIN+qq];
		if(bases[j].nfree==1){
		  
		  mu=0;
		  
		  /*Here now we evaluate the probability of binding*/
		  
		  nc1=ncross[i];
		  nc2=ncross[j];
		  i1=0;
		  i2=0;
		  flagsep=get_distancePBC(bases, traj, i, j, deltat, bindrad[mu], ncross,  crosspart, crossint, i1, i2, mu, cross_rxn, i_home, it, Rmax, sep, R1, plist); 
		  if(flagsep==1){
		    //evaluate probability
		    
		    ratio=bindrad[mu]/R1;
		    if(sep<0){
		      cout <<"separation <0: "<<sep<<" r1 "<<R1<<" p1: "<<i<<" p2: "<<j<<" Rep: "<<rep<<" it "<<it<<endl;
		      sep=0;
		      ratio=1;
		      R1=bindrad[mu];
		    }
		    //aexp=sep/sqrt(4.0*limit.Dmat_ps*deltat);
		    
		    /*add in error model*/
		    currnorm=1.0;
		    p0_ratio=1.0;
		    //place=0;
		    if(R1<limit.fit1){
		      probvec1=survive_irrF(R1, deltat, limit.Dmat_ps[0], bindrad[0], alpha1, cof1);
		    }else{
		      probvec1=survive_irrF(R1, deltat, limit.Dmat_ps[1], bindrad[0], alpha2, cof2);
		    }
		    
		    place2=i*Ncopy[0]+j;
		    if(previter[place2]==(it-1)){
		      p0_ratio=lookup_pirr_pfree_ratio( R1,  prevsep[place2],  deltat,  ps_prev[place2], rtol,  rc, Dtot,  bindrad[0],  FitMat, limit, probvec1);
		      if(isnan(p0_ratio)){
			cout <<"Ratio is NAN!: "<<p0_ratio<<" rep: "<<rep<<" iter: "<<it<<" R1: "<<R1<<" prevsep: "<<prevsep[place2]<<" ps_prev: "<<ps_prev[place2]<<endl;
		      }
		      currnorm=prevnorm[place2]*p0_ratio;
		    }
		    //cout <<"pair: "<<i<<' '<<j<<" prevnorm: "<<prevnorm[0][place2]<<" currnorm: "<<currnorm<<endl;
		    // ep1=exp(2.0*aexp*bexp+bexp*bexp);
		    // 		  if(isinf(ep1)){
		    // 		    real(z)=0.0;
		    // 		    imag(z)=aexp+bexp;
		    // 		    //cout <<"Complex number: "<<z<<endl;
		    // 		    relerr=0;
		    // 		    value=Faddeeva::w(z, relerr);
		    // 		    ea2=exp(-aexp*aexp);
		    // 		    term2=ea2*real(value);
		    // 		  }else{
		    // 		    term2=ep1*erfc(aexp+bexp);
		    // 		  }
		    //probvec1=ratio*kact/(kact+kdiff)*(erfc(aexp)-term2);
		    if(isnan(probvec1))cout <<"SURVIVE PROB IS NAN: "<<probvec1<<endl;
		    probvec[i][nc1]=probvec1*currnorm;
		    probvec[j][nc2]=probvec[i][nc1];
		    prevsep[place2]=R1;
		    previter[place2]=it;
		    prevnorm[place2]=currnorm;
		    ps_prev[place2]=1.0-probvec1*currnorm;
		    
		    //put in matrix element for i1<i2
		    place2=j*Ncopy[0]+i;
		    prevsep[place2]=R1;
		    previter[place2]=it;
		    prevnorm[place2]=currnorm;
		    ps_prev[place2]=1.0-probvec1*currnorm;
		    
		    R2=R1*R1;
		    dudr=-rc/(R2);
		    dx=bases[i].xcom-bases[j].xcom;
		    dy=bases[i].ycom-bases[j].ycom;
		    dz=bases[i].zcom-bases[j].zcom;
		    
		    dx-=plist.xboxl*round(dx/plist.xboxl);
		    dy-=plist.yboxl*round(dy/plist.yboxl);
		    dz-=plist.zboxl*round(dz/plist.zboxl);
		    
		    Fx=-dudr*dx/R1;
		    Fy=-dudr*dy/R1;
		    Fz=-dudr*dz/R1;
		    k1=bases[i].mycomplex;
		    k2=bases[j].mycomplex;
		    Ftraj[k1][0]+=ind_com[k1].Dx*deltat*Fx;
		    Ftraj[k1][1]+=ind_com[k1].Dy*deltat*Fy;
		    Ftraj[k1][2]+=ind_com[k1].Dz*deltat*Fz;
		    Ftraj[k2][0]-=ind_com[k2].Dx*deltat*Fx;
		    Ftraj[k2][1]-=ind_com[k2].Dy*deltat*Fy;
		    Ftraj[k2][2]-=ind_com[k2].Dz*deltat*Fz;
		  
		  }else{
		    /*Still need to calculate contribution to force!!!*/
		    //R1=prevsep[myplace];
		    R2=R1*R1;
		    dudr=-rc/(R2);
		    dx=bases[i].xcom-bases[j].xcom;
		    dy=bases[i].ycom-bases[j].ycom;
		    dz=bases[i].zcom-bases[j].zcom;
		    
		    dx-=plist.xboxl*round(dx/plist.xboxl);
		    dy-=plist.yboxl*round(dy/plist.yboxl);
		    dz-=plist.zboxl*round(dz/plist.zboxl);
		    
		    Fx=-dudr*dx/R1;
		    Fy=-dudr*dy/R1;
		    Fz=-dudr*dz/R1;
		    k1=bases[i].mycomplex;
		    k2=bases[j].mycomplex;
		    Ftraj[k1][0]+=ind_com[k1].Dx*deltat*Fx;
		    Ftraj[k1][1]+=ind_com[k1].Dy*deltat*Fy;
		    Ftraj[k1][2]+=ind_com[k1].Dz*deltat*Fz;
		    Ftraj[k2][0]-=ind_com[k2].Dx*deltat*Fx;
		    Ftraj[k2][1]-=ind_com[k2].Dy*deltat*Fy;
		    Ftraj[k2][2]-=ind_com[k2].Dz*deltat*Fz;
		  
		  }
		}//j is free
		
	      }//i is free
	    }//j proteins
	  }//all proteins
	  
	  //stop_timer(&bimoltime);
	  //cout <<"timer duration, bimolecular loop: "<<timer_duration(bimoltime)<<endl;
	  
	  for(i=0;i<Ntotalmol;i++){
	    //cout <<"protein: "<<i<<" ncross; "<<ncross[i]<<endl;
	    if(ncross[i]>0){
	      
	      /*might perform this reaction, depending on k_associate*/
	      if(ncross[i]>1){
		cout <<"Ncross: "<<ncross[i]<<" protein: "<<i<<" iter: "<<it<<" rep: "<<rep<<endl;
	      }
	      //if(i==0)
	      rnum=1.0*rand_gsl();
	      
	      p1=i;
	      ci1=0;//crossint[p1][0];
	      ci2=0;//crossint[p2][0];
	      p2=crosspart[p1][ci1]; //correct if this is who the reaction is for! 
	      
	      /*Now use this instead to decide which of the potential reactions to choose from*/
	      /*THIS IS ONLY HERE TO LIMIT THE NUMBER OF RNUMS TO SAME AS OTHER PROGRAM*/
	      flag=0;
	      // if(i==0)
	      flag=choose_one_reaction(rnum, p1, ncross, crosspart, probvec, ci1, ci2, cross_rxn, crossint, irandmax);
	      /*if flag=0, no reaction, otherwise there was a reaction involving p1 and protein p2 */
	      	      
	      if(flag==1){
		/*Make sure you choose the first thing that happens*/
		p2=crosspart[p1][ci1];  
		
		rxn1=cross_rxn[p1][ci1];
		cout <<"Associate between proteins: p1 "<<p1<<' '<<p2<<" interfaces: "<<crossint[p1][ci1]<<' '<<crossint[p2][ci2]<<" reaction: "<<rxn1<< " pact; "<<probvec[p1][ci1]<<" rnum: "<<rnum<<" iter: "<<it<<endl;
		associate_translate_measurePBC(p1,  p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist);
		// associate_int_zfirst(p1,  p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist, inst_pro);
// 		check_complex_rad(p1,  bases, ind_com,  plist.xboxl, plist.yboxl, plist.zboxl);
		// /*These two proteins will be deleted*/
// 		bases[p1].nfree=0;
// 		bases[p2].nfree=0;
// 		bases[p1].nbnd=1;
// 		bases[p2].nbnd=1;//
// 		bases[p2].xcom=bases[p1].xcom;
// 		bases[p2].ycom=bases[p1].ycom;
// 		bases[p2].zcom=bases[p1].zcom;
		//plist.ntotalcomplex-=1;
		/*this will set ncross[p1] and [p2] to -1 so they stop moving*/
		remove_reaction_all(p1, p2, ncross, crosspart, probvec, ci1, ci2, cross_rxn, crossint);
		
		movestat[p1]=2;
		movestat[p2]=2;

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
		  
		  // Ftraj[k][0]=ind_com[k].Dx*deltat*Fx;
// 		  Ftraj[k][1]=ind_com[k].Dy*deltat*Fy;
// 		  Ftraj[k][2]=ind_com[k].Dz*deltat*Fz;
		  //traj_complex_rad(i, bases, ind_com,  plist.xboxl, plist.yboxl, plist.zboxl, traj);
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
	      if(movestat[i]==0){
		//don't move if already moved, or attached to someone who will 
		//potentially have to resample their position (movestat=1)
		/*This protein is outside of the reaction zone, so it will move via free diffusion, without a potential.*/
		k=bases[i].mycomplex;
		dx=sqrt(2.0*deltat*ind_com[k].Dx)*GaussV();
		dy=sqrt(2.0*deltat*ind_com[k].Dy)*GaussV();
		dz=sqrt(2.0*deltat*ind_com[k].Dz)*GaussV();
		//cout <<"free dx: "<<k<<' '<<dx<<' '<<dy<<' '<<dz<<endl;
		traj[k][0]=dx+Ftraj[k][0];
		traj[k][1]=dy+Ftraj[k][1];
		traj[k][2]=dz+Ftraj[k][2];
		
		move_protein_posPBC(i, bases, ind_com, traj, movestat, plist); 
		//check_complex_rad(i,  bases, ind_com,  plist.xboxl, plist.yboxl, plist.zboxl);
		Ftraj[k][0]=0;//move_protein_pos sets traj to 0, so set Ftraj to zero as well.
		Ftraj[k][1]=0;
		Ftraj[k][2]=0;
		
	      }
	      
	    }
	    
	  }

	  /*Now we have to check for overlap!!!*/
	  
	  for(i=0;i<Ntotalmol;i++){
	    /*Any protein that reacted will have ncross=-1 and will have already moved. */
	    
	    
	    if(ncross[i]>0){
	     
		/*For any protein that overlapped and did not react, check whether it overlaps with its partners*/
	      sweep_separation_vr_allPBC(deltat, i, bases, ind_com, ncross, crosspart, crossint, cross_rxn, traj, probvec, plist, movestat, i_home, bindrad, Rlist, Ftraj); 
		/*This protein will not move again because removed from all its partners crossing lists*/
		//remove_protein_reaction(i, ncross, crosspart, probvec, cross_rxn, crossint);
	      
	      
	      /*else if you have one partner and it is lower than you, you already tested for overlap
		so don't go through measuring the distance and moving the particles again
	      */
	    }else if(ncross[i]==0){
	      
	      if(movestat[i]!=2){

		
		/*the movestat should be 1, if it were zero then
		  the traj would not have been selected. this could not be the case if it had no overlap, and if it did have overlap, it would've selected a new position 
		*/
		move_protein_posPBC(i, bases, ind_com, traj, movestat, plist); 
		//check_complex_rad(i,  bases, ind_com,  plist.xboxl, plist.yboxl, plist.zboxl);
		k=bases[i].mycomplex;
		Ftraj[k][0]=0;//move_protein_pos sets traj to 0, so set Ftraj to zero as well.
		Ftraj[k][1]=0;
		Ftraj[k][2]=0;
		
	      }
	      
	    }
	  }
	  NA=Ncsave-(Ncsave-plist.ntotalcomplex)*2.0;//NA=Ntotalmol;
	  
	  avgc[it]+=NA;
	  if(it%plist.statwrite==0){
	    cout <<"iteration: "<<it<<" time: "<<it*deltat<<" Na: "<<NA<<endl;
	    cvstime<<it*deltat<<'\t'<<NA<<endl;
	    
	    
	  } 
	  
	}//end iterating over time steps
      
      
      
  }//end looping over repetitions
  double tval;
  double Nhist;
  char tname[200];
  double space;
  ofstream probfile("avg_ncom_vs_time.dat");
  
  for(i=1;i<Nit;i+=plist.statwrite){
    probfile <<i*deltat<<'\t'<<avgc[i]/(1.0*Nrep)<<endl;
  }
  stop_timer(&totaltime);
  cout <<timer_duration(totaltime)<<" total time "<<endl;  
  /*Write out final result*/
  cout <<"End Main, complete run "<<endl;

}//end main
