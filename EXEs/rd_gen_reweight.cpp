/*


Implement particle-based Reaction diffusion
algorithm with trajectory reweighting.

This version is general, allows any number of 
particl types and reaction types, and rigid body
motion includes translation and rotation.

Particles reflect at box edges.
  
  
  
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
#include "vector_rot_calls.h"
#include "cell_neighbor_lists.h"

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
  
  //seed=1380925952;
  cout <<"seed: "<<seed<<endl;
  srand_gsl(seed);
  double randmax=pow(2.0, 32);
  double irandmax=1.0/randmax;
  ifstream parmfile(argv[1]);
  Parms plist;
  plist.restart=0;//in case you don't read it in.
  plist.pclath=-1;//in case you don't read it in, no clath.
  read_parms(parmfile, plist);
  //  write_parms(plist);
  cout <<"Pclath: "<<plist.pclath<<endl;
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
  plist.Natomwrite=0;
  for(i=0;i<Nprotypes;i++){  
    ntmp=wholep[i].ninterface+1;
    plist.Natom+=Ncopy[i]*ntmp;
    Ntotsite+=Ncopy[i]*wholep[i].ninterface;
    plist.Natomwrite+=Ncopy[i]*wholep[i].nint_write;
  }
  cout <<"N atoms: "<<plist.Natom<<" Nsites, not COM: "<<Ntotsite<<" Natoms to write out: "<<plist.Natomwrite<<endl;
  int t=0;
  
  cout <<"read reactions "<<endl;
  int *rxtype=new int[plist.Nrxn];
  double *Kd=new double[plist.Nrxn];
  read_reactions(rxnfile, plist, Rlist, Ncoup, mycoupled, Nmyrxn, myrxn, bindrad, kr, cntrxn, freelist, bndlist, zlist, rxtype, Kd);

  /*Change rates to ka and kb, rather than kon and koff
    unless you are reading in values of ka and kb
  */
  //  update_rates(wholep, plist, Rlist,  bindrad, kr, rxtype, Kd,p_home);

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
  int twrite=plist.configwrite;
  int gwrite=plist.grwrite;
  sprintf(fname,"ComplexW_COM_Np%d_Ni%d.xyz",plist.Nprotypes, Nifaces); 
  ofstream compout(fname);
  sprintf(fname,"ProteinW_COM_Np%d_Ni%d.xyz",plist.Nprotypes, Nifaces); 
  ofstream proout(fname);
  /*cells*/
  double Dtot=wholep[0].Dx+wholep[1].Dx;
  double Dtot1=Dtot;
  double deltat=plist.deltat;
  double Rmax1=3.0*sqrt(6.0*Dtot*deltat)+bindrad[0];
  double tx, ty, tz;
  double *M=new double[9];
  int myc1, myc2;
  
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

  /*Each cell has 26 neighbors (unless at boundary and reflecting BC is used). 
    Since we will loop over all cell pairs, only keep track of half your neighbors.
  */
  int maxnbor=13; //geometry of cube
  int *Nnbor=new int[Ncell];
  int *nbor=new int[Ncell*maxnbor];
  //int *nborrev=new int[Ncell*maxnbor];
  
  int *npb=new int[Ncell];
  int MAXPERBIN=200;//int(Ntotalmol);
  int *binlist=new int[Ncell*MAXPERBIN];
  int MAXDISS=100;//max number of proteins to dissociate in one time step 
  int *disslist=new int[MAXDISS];//holds index of dissociated proteins
  int ndiss;
  int mybin;
  int mybinind;
  int c;
  
  cell_neighbor_list(Nx, Ny, Nz, maxnbor, nbor, Nnbor);

  cout <<"N cell pairs (max is with PBC): "<<Ncell*maxnbor<<endl;

  char tname4[100];
  sprintf(tname4,"timestat.dat");
  ofstream timestatfile(tname4,std::ios::app);

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
  
  double **traj=new double *[Ntotalmol];
  double **trajR=new double *[Ntotalmol];
  for(i=0;i<Ntotalmol;i++){
    traj[i]=new double[3];
    trajR[i]=new double[3];
  }

  double kpi=4.0*M_PI*bindrad[0];
  double fourpi=4.0*M_PI;
  double kact;
  double kdiff=4.0*M_PI*Dtot*bindrad[0];
  double fact=1.0+kr[0]/kdiff;
  double alpha=fact*sqrt(Dtot)/bindrad[0];
  double cof=kr[0]/(kr[0]+kdiff);

  cout <<"activation rate in nm^3/us 1st reaction: "<<kr[0]<<" kd: "<<kr[1]<<" Kd, uM: "<<kr[1]/kr[0]/1E6/6.022E23*1E24*1E6<<endl;
  cout <<"Dtot: "<<Dtot<<" alpha: "<<alpha<<" cof: "<<cof<<endl;
  double R1;
  double Rmax=Rmax1;
  double r0, passoc;
  
  /*For rotating proteins*/
  double rleg2;
  double Dr1, Dr2;
  dx=bases[0].x[0]-bases[0].xcom;
  dy=bases[0].y[0]-bases[0].ycom;
  dz=bases[0].z[0]-bases[0].zcom;
  rleg2=dx*dx+dy*dy+dz*dz;
  double rleg=sqrt(rleg2);
  double cf=cos(sqrt(4.0*bases[0].Drx*deltat));
  Dr1=2.0*rleg2*(1.0-cf);
  cout <<"displacement2 from rotation on protein 1: "<<Dr1<<" R: "<<sqrt(rleg2)<<" Dr "<<bases[0].Drx<<" diffusion from rot: "<<Dr1/(6.0*deltat)<<endl;
  
  /*Below use einstein-stokes to define D based on the particle
    radius and the Temperature and viscosity. 
    To enforce the D read in from file, calculated via, e.g. bead models,
    set scale>1 below.
  */
  double Temp=293;//K
  double nu=0.001;//kg/(m*s)
  double scale=3.0;//greater the one to correct for non-spherical
  double crad=wholep[0].radx;
  if(crad==0)crad=1;
  cout <<"clathrin radius: "<<crad<<" nm. "<<endl;
  plist.pretrans=trans_prefactor(Temp, nu,  scale, crad, wholep[0].Dx);
  plist.prerot=rot_prefactor(Temp, nu,  scale, crad, wholep[0].Drx);
  cout <<"Diffusion prefactors: "<<plist.pretrans<<' '<<plist.prerot<<endl;
  cout <<" Passoc_vs_separation "<<endl; 
  cout <<"Dtot: "<<Dtot<<endl;
  passoc=survive_irr(bindrad[0], deltat, Dtot, bindrad[0], alpha, cof);
  cout <<bindrad[0]<<'\t'<<passoc<<endl;
  int loop=0;
  i=0;
  while(loop==0){
    r0=(i+0.5)*0.001+bindrad[0];
    passoc=survive_irr(r0, deltat, Dtot, bindrad[1], alpha, cof);
    cout <<r0<<'\t'<<passoc<<endl;
    if(passoc<1E-13){
      
      //Rmax=r0;
      loop=1;
      cout<<"passoc: "<<passoc<<" i: "<<i<<" r0: "<<r0<<endl;
      
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

  int *nprevpart=new int[Ntotalmol];
  int *ncurrpart=new int[Ntotalmol];
  int **prevlist=new int*[Ntotalmol];
  int **currlist=new int*[Ntotalmol];
  int **prevmyface=new int*[Ntotalmol];
  int **currmyface=new int*[Ntotalmol];
  int **prevpface=new int*[Ntotalmol];
  int **currpface=new int*[Ntotalmol];
  double **prevnorm=new double*[Ntotalmol];
  //int **previter=new int*[Ntotalmol];
  double **ps_prev=new double*[Ntotalmol];
  double **prevsep=new double*[Ntotalmol];
  double **currprevnorm=new double*[Ntotalmol];
  //int **previter=new int*[Ntotalmol];
  double **currps_prev=new double*[Ntotalmol];
  double **currprevsep=new double*[Ntotalmol];
  
  /*This size MAXNORM is keep track of how many particles are in each proteins reaction
    zone at one time. Should be set large to buffer for fluctuations to large number
    of partners, even though each interface (and protein) should ideally have only 1 
    partner in its reaction zone at each step.
  */
  int MAXNORM=100;
  int s;
  int ssave;
  for(i=0;i<Ntotalmol;i++){
    nprevpart[i]=0;
    ncurrpart[i]=0;
    
    currlist[i]=new int[MAXNORM];
    prevlist[i]=new int[MAXNORM];
    currmyface[i]=new int[MAXNORM];
    prevmyface[i]=new int[MAXNORM];
    currpface[i]=new int[MAXNORM];
    prevpface[i]=new int[MAXNORM];
    prevnorm[i]=new double[MAXNORM];
    ps_prev[i]=new double[MAXNORM];
    prevsep[i]=new double[MAXNORM];
    currprevnorm[i]=new double[MAXNORM];
    currps_prev[i]=new double[MAXNORM];
    currprevsep[i]=new double[MAXNORM];
    //previter[i]=new int[MAXNORM];
  }
  //cerr<<"allocated mem "<<endl;
  int myface, pface;
  int place2;
  double rtol=1E-10;
  
  int flag;
  int p1, p2;
  double probvec1, p0_ratio, currnorm;


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
  double memfact=1.0;//If you are binding to PIP2, this should be 1/2.
  int pp, qq, nb, hh;
  
  
  
  double maxrad=150;
  int radbins=3000;
  double delrad=maxrad/(1.0*radbins);
  int ind_rad;
  double rad2, rad;
  
  int Ncsave=plist.ntotalcomplex;
  int Nccurr;
  double x0;
  int proa, pro2;
  double currx, curry, currz;
  /*initial separation between each pair*/

  int Nrep=1;
  int rep=0;
  int totpercell=0;
  int dub;
  get_bin2( plist, bases, cellx,  celly,cellz,  Nx,  Ny,  Nz, binlist, npb, MAXPERBIN, Ncell, ind_com);  
  for(i=0;i<Ncell;i++){
    cout <<"cell: "<<i<<" npercell: "<<npb[i]<<endl;
    totpercell+=npb[i];
  }
  cout <<"Done defining bins. Total mols across cells "<<totpercell<<endl;
  
  write_timestat(timestatfile, wholep, bases, plist, Ncopy, it, ind_com, deltat,Nprotypes);///was assemblyfile

  for(it=1;it<Nit+1;it++)
    {
      
      for(i=0;i<Ntotalmol;i++){
	ncross[i]=0;
	movestat[i]=0;
      }
      
      ndiss=0;
      /*Test dissociation!! */
      for(i=0;i<Ntotalmol;i++){
	/*Test this proteins unimolecular interactions*/
	dub=0;
	for(j=0;j<bases[i].nbnd;j++){
	  icom=bases[i].bndlist[j];//index of the bound species
	  for(n=0;n<Nmyrxn[icom];n++){
	    mu=myrxn[icom][n];//reaction number
	    
	    i1=Rlist[mu][1];//determine interface for this reaction
	    i2=Rlist[mu][2];
	    iind=i_home[i1];
	    iind2=i_home[i2];
	    
	    if(p_home[i1]!=p_home[i2]){
	      if(p_home[i1]!=bases[i].protype){
		i1=Rlist[mu][2];
		i2=Rlist[mu][1];
		iind=i_home[i1];
		iind2=i_home[i2];
	      }
	    }else{
	      /*If two clathrins are dissociating, figure out which interface is on which protein.*/
	      if(i1!=i2){
		if(bases[i].istatus[iind]==icom){
		  if(bases[i].istatus[iind2]==icom){
		    if(dub==1){
		      /*check other interface as well */
		      i1=Rlist[mu][2];
		      i2=Rlist[mu][1];
		      iind=i_home[i1];
		    }
		    dub++;
		  }
		}else{
		  /*i1 is other interface */
		  i1=Rlist[mu][2];
		  i2=Rlist[mu][1];
		  iind=i_home[i1];
		  
		}
	      }//if they are the same, no need to switch anything
	    }
	    rate=kr[mu];
	    ppart=bases[i].partner[iind];
	    
	    if(i<ppart){
	      /*then this is the first protein in the reaction, so try it.
		This if statement ensures we do not try to dissociate the same complex twice
	      */
	      
	      
	      prob=1-exp(-rate*deltat*us_to_s);
	      rnum=1.0*rand_gsl();
	      if(prob>rnum){
		rnum2=rnum+rand_gsl()*irandmax;//to get higher resolution random numbers
		if(prob>rnum2){
		  cout <<"Dissociate at iter: "<<it<<" protein: "<<i<<" partner: "<<ppart<<" randomnum: "<<rnum2<<endl;
		  
		  /*Perform this dissociation reaction.
		    Sometimes it is a bond broken, not a full dissociation to two complexes if
		    the two interfaces are part of the same complex
		  */
		  
		  cancel=break_complex_trans(i, mu, j, bases, Rlist, i_home, ind_com, plist, bindrad,  ppart, i1, i2, p_home, myrxn);
		  ap2_coupledrxn_sub(mu, Ncoup, mycoupled, bases, Rlist , i_home, p_home, i,ppart);
		  
		  j--;//replaced this reaction, so stay on this one
		  disslist[ndiss]=i;
		  disslist[ndiss+1]=ppart;
		  ndiss+=2;
		  
		  /*If dissociated products are removed from overlap lists, use ncross=-1.
		    If they remain in list to avoid overlap, use movestat=2 and also
		    ensure that they are not allowed to diffuse again, by, for example,
		    temporarily setting D=0.
		  */
		  ncross[i]=-1;
		  ncross[ppart]=-1;
		  
		  reflect_complex_rad_rot(i,  bases, ind_com,  plist.xboxl, plist.yboxl, plist.zboxl);
		  reflect_complex_rad_rot(ppart,  bases, ind_com,  plist.xboxl, plist.yboxl, plist.zboxl);
		  
		  //movestat[i]=2;
		  //movestat[ppart]=2;
		  
		}
	      }
	      
	    }//only try each pair dissociating once
	  }
	}//finished trying out unimolecular reactions for this proteins
      }
      
      get_bin2( plist, bases, cellx,  celly,cellz,  Nx,  Ny,  Nz, binlist, npb, MAXPERBIN, Ncell, ind_com);

      
      /*Either remove dissociated proteins from the binlist
	so they do not try to react again this time step, or
	ensure that they do not move again (D=0) and they
	will not react, the other proteins will just avoid
	overlapping with them.
      */
      for(i=0;i<ndiss;i++){
 	p1=disslist[i];
 	k=bases[p1].mycomplex;

	/*ind_com[k].Dx=0;
 	ind_com[k].Dy=0;
 	ind_com[k].Dz=0;//these changes need to be reversed at the end of the timestep
	*/
	
	mybin=bases[p1].mybin;
	mybinind=bases[p1].mybinind;
	
	pnew=binlist[mybin*MAXPERBIN+npb[mybin]-1];
	bases[pnew].mybinind=mybinind;
	binlist[mybin*MAXPERBIN+mybinind]=pnew;
	npb[mybin]-=1;
      }	    
      
      /*Keep track of partners in reaction zone for reweighting*/
      for(i=0;i<Ncopy[0];i++){
	ncurrpart[i]=0;
      }
      
      /*Measure separations between proteins in neighboring cells to identify
	all possible reactions.
      */
      
      for(c=0;c<Ncell;c++){
	
	for(pp=0;pp<npb[c];pp++){
	  i=binlist[c*MAXPERBIN+pp];
	  
	  /*Test bimolecular reactions!*/
	  nfree=bases[i].nfree;
	  wprot=bases[i].protype;
	  
	  if(nfree>0){
	    	  
	    for(qq=pp+1;qq<npb[c];qq++){
	      j=binlist[c*MAXPERBIN+qq];
	      
	      /*test to see if proteins i and j interact.*/
	      
	      flag=0;
	      for(n=0;n<wholep[wprot].npropart;n++){
		if(wholep[wprot].propart[n]==bases[j].protype){
		  flag=1;//check this protein for binding
		  n=wholep[wprot].npropart;//break from loop
		}
	      }
	      
	      if(flag==1){
		/*These two proteins are capable of interacting,
		  test to see if those interfaces are free to bind
		*/
		
		for(p=0;p<bases[i].nfree;p++){
		  i1=bases[i].freelist[p];
		  /*test all of i1's binding partners to see whether they are on protein j */
		  np=numpartners[i1];
		  for(n=0;n<np;n++){
		    i2=Speclist[i1][n];//binding interface partner
		    for(k=0;k<bases[j].nfree;k++){
		      if(bases[j].freelist[k]==i2){
			//both binding interfaces are available!
			mu=myrxn[i1][n];
			
			/*Here now we evaluate the probability of binding*/
			/*Different interfaces can have different diffusion constants if they 
			  also rotate. <theta^2>=6Drdeltat.
			  In that case, <d>=sin(sqrt(6Drdeltat)/2)*2R
			  so add in <d>^2=4R^2sin^2(sqrt(6Drdeltat)/2) 
			  for a single clathrin, R=arm length,
			  otherwise R will be from pivot point of rotation, COM, 
			  calculate distance from interface to the complex COM
			*/
			myc1=bases[i].mycomplex; 
			myc2=bases[j].mycomplex;
			Dtot=ind_com[myc1].Dx+ind_com[myc2].Dx;
			iind=i_home[i1];
			dx=bases[i].x[iind]-ind_com[myc1].xcom;
			dy=bases[i].y[iind]-ind_com[myc1].ycom;
			dz=bases[i].z[iind]-ind_com[myc1].zcom;
			
			rleg2=dx*dx+dy*dy+dz*dz;
			cf=cos(sqrt(4.0*ind_com[myc1].Drx*deltat));
			Dr1=2.0*rleg2*(1.0-cf);
			iind2=i_home[i2];
			dx=bases[j].x[iind2]-ind_com[myc2].xcom;
			dy=bases[j].y[iind2]-ind_com[myc2].ycom;
			dz=bases[j].z[iind2]-ind_com[myc2].zcom;
			
			  
			rleg2=dx*dx+dy*dy+dz*dz;
			cf=cos(sqrt(4.0*ind_com[myc2].Drx*deltat));
			Dr2=2.0*rleg2*(1.0-cf);
			
			
			Dtot+=(Dr1+Dr2)/(6.0*deltat);//add in contributions from rotation
			//if(abs(Dtot-Dtot1)>rtol)cout <<"TOO HIGH DIFFUSION: "<<Dtot<<" p1: "<<i<<" p2: "<<j<<endl;
			/*Reaction zone radius between particle positions*/
			Rmax=3.0*sqrt(6.0*Dtot*deltat)+bindrad[mu];
			
			nc1=ncross[i];
			nc2=ncross[j];
			flagsep=get_distance(bases, traj, i, j, deltat, bindrad[mu], ncross,  crosspart, crossint, i1, i2, mu, cross_rxn, i_home, it, Rmax, sep, R1); 
			if(flagsep==1){
			  /*Evaluate reaction probability*/
			  ratio=bindrad[mu]/R1;
			  if(sep<0){
			    if(myc1!=myc2){
			       cout <<"separation <0: "<<sep<<" r1 "<<R1<<" p1: "<<i<<" p2: "<<j<<" Rep: "<<rep<<" it "<<it<<endl;
			       cout <<bases[i].xcom<<' '<<bases[i].ycom<<' '<<bases[i].zcom<<" nfree: "<<bases[i].nfree<<endl;
			       cout <<bases[j].xcom<<' '<<bases[j].ycom<<' '<<bases[j].zcom<<" nfree: "<<bases[j].nfree<<endl;
			    }
			    sep=0;
			    ratio=1;
			    R1=bindrad[mu];
			  }
			  /*If one particle (lipid) is bound in the membrane, reaction prob
			    is half due to flux across only top half of particle.
			  */
			  memfact=1.0;
			  if(bases[i].Dz==0 || bases[j].Dz==0)memfact=1.0;
			  kdiff=fourpi*Dtot*bindrad[mu]*memfact;
			  
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
			  
			  /*protein i is protype wprot and j is wprot2*/
			  wprot2=bases[j].protype;
			  proa=i;
			  myface=i1;
			  pro2=j;
			  pface=i2;
			  if(i>j){
			    proa=j;
			    myface=i2;
			    pro2=i;
			    pface=i1;
			  }
			  
			  /*Find out what the previous reweighting was for this pair, if they 
			    were stored from previous step (inside reaction zone).
			  */
			  if(myc1!=myc2){
			    for(s=0;s<nprevpart[proa];s++){
			      if(prevlist[proa][s]==pro2 && prevmyface[proa][s]==myface &&prevpface[proa][s]==pface){
				
				p0_ratio=pirr_pfree_ratio_psF(R1, prevsep[proa][s], deltat, Dtot, bindrad[mu], alpha, ps_prev[proa][s], rtol);
				currnorm=prevnorm[proa][s]*p0_ratio;
				
				ssave=s;
				s=nprevpart[proa];
			      }
			    }
			    /*don't renormalize if their positions are fixed by being in the same complex!*/
			  }
			  probvec1=ratio*kact/(kact+kdiff)*(erfc(aexp)-exp(2.0*aexp*bexp+bexp*bexp)*erfc(aexp+bexp));
			  probvec[i][nc1]=probvec1*currnorm;
			  probvec[j][nc2]=probvec[i][nc1];
			  
			  /*Store all the reweighting numbers for next step.*/
			  s=ncurrpart[proa];
			  currprevsep[proa][s]=R1;
			  //cout <<"proa: "<<proa<<" s: "<<s<<" pro2: "<<pro2<<" prevsep: "<<prevsep[proa][ssave]<<" probvec: "<<probvec1<<" normfact: "<<currnorm<<" R: "<<R1<<" myface: "<<myface <<" other: "<<pface<<" nprevpart: "<<nprevpart[proa]<<endl;
			  currlist[proa][s]=pro2;
			  currmyface[proa][s]=myface;
			  currpface[proa][s]=pface;
			  currprevnorm[proa][s]=currnorm;
			  currps_prev[proa][s]=1.0-probvec1*currnorm;
			  ncurrpart[proa]++;
			}//Within the reaction zone
		      }
		    }
		  }
		}
	      }//These protein partners interact
	    }//loop over protein partners in your same cell
	  
	    /*Now loop over all neighboring cells, and all proteins in those cells.
	      for PBC, all cells have maxnbor neighbor cells. For reflecting, edge have fewer. 
	    */
	    for(hh=0;hh<Nnbor[c];hh++){
	      nb=nbor[c*maxnbor+hh];
	      for(qq=0;qq<npb[nb];qq++){
		j=binlist[nb*MAXPERBIN+qq];
		
		/*test to see if  proteins  i and j react with one another*/
		
		flag=0;
		for(n=0;n<wholep[wprot].npropart;n++){
		  if(wholep[wprot].propart[n]==bases[j].protype){
		    flag=1;//check this protein for binding
		    n=wholep[wprot].npropart;//break from loop
		  }
		}
		
		if(flag==1){
		  /*These two proteins are capable of interacting,
		    test to see if those interfaces are free to bind
		  */
		  
		  for(p=0;p<bases[i].nfree;p++){
		    i1=bases[i].freelist[p];
		    /*test all of i1's binding partners to see whether they are on protein j */
		    np=numpartners[i1];
		    for(n=0;n<np;n++){
		      i2=Speclist[i1][n];//binding interface partner
		      for(k=0;k<bases[j].nfree;k++){
			if(bases[j].freelist[k]==i2){
			  //both binding interfaces are available!
			  mu=myrxn[i1][n];
			  
			  /*Here now we evaluate the probability of binding*/
			  /*Different interfaces can have different diffusion constants if they 
			    also rotate. <theta^2>=6Drdeltat.
			    In that case, <d>=sin(sqrt(6Drdeltat)/2)*2R
			    so add in <d>^2=4R^2sin^2(sqrt(6Drdeltat)/2) 
			    for a single clathrin, R=arm length,
			    otherwise R will be from pivot point of rotation, COM, 
			    calculate distance from interface to the complex COM
			  */
			  myc1=bases[i].mycomplex; 
			  myc2=bases[j].mycomplex;
			  Dtot=ind_com[myc1].Dx+ind_com[myc2].Dx;
			  iind=i_home[i1];
			  dx=bases[i].x[iind]-ind_com[myc1].xcom;
			  dy=bases[i].y[iind]-ind_com[myc1].ycom;
			  dz=bases[i].z[iind]-ind_com[myc1].zcom;
			  	
			  rleg2=dx*dx+dy*dy+dz*dz;
			  cf=cos(sqrt(4.0*ind_com[myc1].Drx*deltat));
			  Dr1=2.0*rleg2*(1.0-cf);
			  iind2=i_home[i2];
			  dx=bases[j].x[iind2]-ind_com[myc2].xcom;
			  dy=bases[j].y[iind2]-ind_com[myc2].ycom;
			  dz=bases[j].z[iind2]-ind_com[myc2].zcom;
			  
			  
			  rleg2=dx*dx+dy*dy+dz*dz;
			  cf=cos(sqrt(4.0*ind_com[myc2].Drx*deltat));
			  Dr2=2.0*rleg2*(1.0-cf);
			  
			  
			  Dtot+=(Dr1+Dr2)/(6.0*deltat);//add in contributions from rotation
			  
			  //if(abs(Dtot-Dtot1)>rtol)cout <<"TOO HIGH DIFFUSION: "<<Dtot<<" p1: "<<i<<" p2: "<<j<<endl;
			  
			  /*Reaction zone radius between particle positions*/
			  Rmax=3.0*sqrt(6.0*Dtot*deltat)+bindrad[mu];
			  
			  nc1=ncross[i];
			  nc2=ncross[j];
			  flagsep=get_distance(bases, traj, i, j, deltat, bindrad[mu], ncross,  crosspart, crossint, i1, i2, mu, cross_rxn, i_home, it, Rmax, sep, R1); 
			  if(flagsep==1){
			    
			    /*Evaluate probability of reaction, with reweighting*/
			    
			    ratio=bindrad[mu]/R1;
			    if(sep<0){
			      if(myc1!=myc2){
				 cout <<"separation <0: "<<sep<<" r1 "<<R1<<" p1: "<<i<<" p2: "<<j<<" Rep: "<<rep<<" it "<<it<<endl;
				 cout <<bases[i].xcom<<' '<<bases[i].ycom<<' '<<bases[i].zcom<<" nfree: "<<bases[i].nfree<<endl;
				 cout <<bases[j].xcom<<' '<<bases[j].ycom<<' '<<bases[j].zcom<<" nfree: "<<bases[j].nfree<<endl;
			      }
			      sep=0;
			      ratio=1;
			      R1=bindrad[mu];
			    }
			    
			    /*If one particle (lipid) is bound in the membrane, reaction prob
			      is half due to flux across only top half of particle.
			    */
			    memfact=1.0;
			    if(bases[i].Dz==0 || bases[j].Dz==0)memfact=1.0;
			    kdiff=fourpi*Dtot*bindrad[mu]*memfact;
			    
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
			    
			    /*protein i is protype wprot and j is wprot2*/
			    
			    proa=i;
			    myface=i1;
			    pro2=j;
			    pface=i2;
			    if(i>j){
			      proa=j;
			      myface=i2;
			      pro2=i;
			      pface=i1;
			    }
			    
			    /*Find out what the previous reweighting was for this pair, if they 
			      were stored from previous step (inside reaction zone).
			    */
			    if(myc1!=myc2){
			      for(s=0;s<nprevpart[proa];s++){
				if(prevlist[proa][s]==pro2 && prevmyface[proa][s]==myface &&prevpface[proa][s]==pface){
				  
				  p0_ratio=pirr_pfree_ratio_psF(R1, prevsep[proa][s], deltat, Dtot, bindrad[mu], alpha, ps_prev[proa][s], rtol);
				  currnorm=prevnorm[proa][s]*p0_ratio;
				  
				  ssave=s;
				  s=nprevpart[proa];
				}
			      }
			      /*don't renormalize if their positions are fixed by being in the same complex!*/
			    }
			    probvec1=ratio*kact/(kact+kdiff)*(erfc(aexp)-exp(2.0*aexp*bexp+bexp*bexp)*erfc(aexp+bexp));
			    probvec[i][nc1]=probvec1*currnorm;
			    probvec[j][nc2]=probvec[i][nc1];
			    
			    /*Store all the reweighting numbers for next step.*/
			    s=ncurrpart[proa];
			    currprevsep[proa][s]=R1;
			    //cout <<"proa: "<<proa<<" s: "<<s<<" pro2: "<<pro2<<" Prevsep: "<<prevsep[proa][ssave]<<" probvec: "<<probvec1<<" normfact: "<<currnorm<<" R: "<<R1<<" myface: "<<myface <<" other: "<<pface<<" nprevpart: "<<nprevpart[proa]<<endl;


			    currlist[proa][s]=pro2;
			    currmyface[proa][s]=myface;
			    currpface[proa][s]=pface;
			    currprevnorm[proa][s]=currnorm;
			    currps_prev[proa][s]=1.0-probvec1*currnorm;
			    ncurrpart[proa]++;
			  }//Within reaction zone
			}
		      }
		    }
		  }
		}//These proteins i and j react with one another
	      }//loop over all proteins in this neighbor cell
	    }//loop over all neighbor cells
	  }//if protein i is free to bind
	}//loop over all proteins in initial cell
      }//End looping over all cells.
      
      /*Now that separations and reaction probabilities are calculated,
	decide whether to perform reactions for each protein.
      */
      
      for(i=0;i<Ntotalmol;i++){
	
	/*Skip any proteins that just dissociated during this time step*/
	
	if(ncross[i]>0){
	  
	  rnum=1.0*rand_gsl();
	  p1=i;
	  ci1=0;//initialize to zero
	  ci2=0;//initialize to zero
	  p2=crosspart[p1][ci1]; //this will be update by choose_one_reaction below
	  
	  /*Evaluate whether to perform a reaction with protein i, and with whom. Flag=1 means
	    reaction is performed. Returns correct ci1 and ci2 for this rxn.
	  */
	  flag=choose_one_reaction(rnum, p1, ncross, crosspart, probvec, ci1, ci2, cross_rxn, crossint, irandmax);
	  
	  if(flag==1){
	    
	    p2=crosspart[p1][ci1];  
	    rxn1=cross_rxn[p1][ci1];
	    
	    /*test for being in same complex to avoid moving binding interfaces too far*/
	    //flag2=same_complex_test(p1, p2, bases, crossint, i_home, bindrad, rxn1, ci1, ci2);
	    
	    cout <<"Associate between proteins: p1 "<<p1<<' '<<p2<<" interfaces: "<<crossint[p1][ci1]<<' '<<crossint[p2][ci2]<<" reaction: "<<rxn1<< " pact; "<<probvec[p1][ci1]<<" rnum: "<<rnum<<" iter: "<<it<<endl;
	    
	    /*Associate proteins, move them to contact and update their free and bound lists*/
	    
	    if(bases[p1].protype==bases[p2].protype){
	      /*for self binding, use freeleg for clathrin type proteins.*/
	      if(bases[p1].protype==plist.pclath)
		associate_freeleg(p1,  p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist);
	      else
		associate_translate_measure(p1,  p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist);
	      
	      
	    }else{
	      /*two distinct proteins. This routine will roate them to align
		their COMs before translating, unless a complex is larger than a single protein. 
	      */
	      /*using associate_translate_measure here will also test for overlap with clathrins.*/
	      associate_translate_measure(p1,  p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist);

	    }
	    /*coupled reactions.*/
	    ap2_coupledrxn_add(rxn1, Ncoup, mycoupled, bases, Rlist , i_home, p_home, p1, p2);
	    
	    reflect_complex_rad_rot(p1,  bases, ind_com,  plist.xboxl, plist.yboxl, plist.zboxl);

	    /*Remove p1 and p2 from the list of potential reaction partners so they don't try to
	      associate again in this turn.
	    */
	    remove_reaction_all(p1, p2, ncross, crosspart, probvec, ci1, ci2, cross_rxn, crossint);
	    
	    
	    /*Since these proteins have moved to associate and taken their complex with them, 
	      don't allow any proteins in their complex to move again.
	    */
	    set_movestat_zero(p1, bases, ind_com, movestat);
	    
	  }else {
	    
	    /*This protein will not associate in this time step. 
	      For this Sweeping version, just move the particle, don't check for overlap until everyone is done.
	      Store new move in the traj vector so you still know the original position store in bases[]
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
	    
	      traj[k][0]=dx;
	      traj[k][1]=dy;
	      traj[k][2]=dz;
	      
	      tx=sqrt(2.0*deltat*ind_com[k].Drx)*GaussV();
	      ty=sqrt(2.0*deltat*ind_com[k].Dry)*GaussV();
	      tz=sqrt(2.0*deltat*ind_com[k].Drz)*GaussV();
	      
	      trajR[k][0]=tx;
	      trajR[k][1]=ty;
	      trajR[k][2]=tz;
	      
	      rotationEuler(trajR[k][0], trajR[k][1], trajR[k][2], M);
	      reflect_traj_complex_rad_rot(i, bases, ind_com,  plist.xboxl, plist.yboxl, plist.zboxl, traj, M);
	      /*this displacement will apply to all the proteinss in this complex k.*/
	      for(j=0;j<ind_com[k].mysize;j++){
		mp=ind_com[k].plist[j];
		movestat[mp]=1;
	      }
	    
	    }
	    	    
	    /*Set probability of this protein to zero in all reactions so it doesn't try to
	      react again but the partners still will avoid overlapping.
	    */
	    remove_one_prob_all(p1, ncross, crosspart, probvec, cross_rxn, crossint);
	    
	  }
	  
	}else if(ncross[i]>-1){
	  
	  
	  /*this protein has ncross=0,
	    meaning it neither dissociated nor tried to associate.
	    however, it could have movestat=2 if it is part of a multi-protein
	    complex that already displaced.
	  */
	  if(movestat[i]==0){
	    /*don't move this protein if it already moved, or attached to someone who will 
	      potentially have to resample their position (movestat=1)*/
	    k=bases[i].mycomplex;
	    dx=sqrt(2.0*deltat*ind_com[k].Dx)*GaussV();
	    dy=sqrt(2.0*deltat*ind_com[k].Dy)*GaussV();
	    dz=sqrt(2.0*deltat*ind_com[k].Dz)*GaussV();
	    
	    traj[k][0]=dx;
	    traj[k][1]=dy;
	    traj[k][2]=dz;
	    
	    tx=sqrt(2.0*deltat*ind_com[k].Drx)*GaussV();
	    ty=sqrt(2.0*deltat*ind_com[k].Dry)*GaussV();
	    tz=sqrt(2.0*deltat*ind_com[k].Drz)*GaussV();
	    
	    trajR[k][0]=tx;
	    trajR[k][1]=ty;
	    trajR[k][2]=tz;
	    
	    move_rot_proteins(i, bases, ind_com, traj, movestat, trajR, M); 
	    reflect_complex_rad_rot(i,  bases, ind_com,  plist.xboxl, plist.yboxl, plist.zboxl);
	  }
	  
	}
	
      }//done testing all molecules for bimolecular reactions
      
      
      /*Now we have to check for overlap!!!*/
      for(i=0;i<Ntotalmol;i++){
	
	/*Any protein that reacted will have ncross=-1 and will have already moved, skip
	  trying to avoid overlap with them.
	  However, this means that other interfaces on reacting proteins will not avoid overlap,
	  this could be changed? 
	*/
		
	if(ncross[i]>0){
	  
	  /*For any protein that overlapped and did not react, check whether it overlaps with its partners*/
	  sweep_separation_rot_all(deltat, i, bases, ind_com, ncross, crosspart, crossint, cross_rxn, traj, probvec, plist, movestat, i_home, bindrad, trajR, M, Rlist); 
	  //sweep_separation_all(deltat, i, bases, ind_com, ncross, crosspart, crossint, cross_rxn, traj, probvec, plist, movestat, i_home, bindrad); 
	  
	}else if(ncross[i]==0){
	  
	  if(movestat[i]!=2){
	    
	    /*For proteins with ncross=0, they either moved independently, or their displacements
	      were selected based on the complex they were part of, and they may not yet been moved.
	    */
	    move_rot_proteins(i, bases, ind_com, traj, movestat,trajR, M); 
	    reflect_complex_rad_rot(i,  bases, ind_com,  plist.xboxl, plist.yboxl, plist.zboxl);
	  }
	  
	}
      }
      //check_com_to_pro(Ntotalmol, plist.ntotalcomplex, bases, ind_com);
      //check_pro_to_com(Ntotalmol, bases, ind_com);
      //check_rot_move(Ntotalmol, bases, ind_com);
      
      if(it%statwrite==0){
	//Nacurr=Ncopy[0]-(Ncsave-plist.ntotalcomplex);
	/*for clathrin. If clathrins close a loop ncomplex doesn't decreas but two more sites bind*/
	Nccurr=Ncsave-plist.ntotalcomplex+plist.nloop;//total products formed
	Nacurr=Ncsave*3.0-Nccurr*2.0;//each product has 2 A sites in it
	cout<<"timestep: "<<it*deltat<<'\t'<<plist.ntotalcomplex<<" Nleg sites, if clathrin: "<<Nacurr<<endl;
      }
      if(it%twrite==0){
	write_protein_iface_short(proout, plist, bases,Ncopy, it, wholep, names);
	write_timestat(timestatfile, wholep, bases, plist, Ncopy, it, ind_com, deltat,Nprotypes);///was assemblyfile
	//write_status(restart, wholep, bases, plist, Ncopy, it);
	//write_protein_iface(restart, plist, bases, Ncopy, it, wholep, names);

      }
      
      /*Now replace all currsep as prevsep, to keep track of
	reweighting values for next time step*/
      
      for(i=0;i<Ntotalmol;i++){
	nprevpart[i]=ncurrpart[i];
	for(s=0;s<ncurrpart[i];s++){
	  prevlist[i][s]=currlist[i][s];
	  prevmyface[i][s]=currmyface[i][s];
	  prevpface[i][s]=currpface[i][s];
	  prevnorm[i][s]=currprevnorm[i][s];
	  ps_prev[i][s]=currps_prev[i][s];
	  prevsep[i][s]=currprevsep[i][s];
	}
      }
      
    }//end iterating over time steps
  
  
  
  
  
  
  stop_timer(&totaltime);
  cout <<timer_duration(totaltime)<<" total time "<<endl;  
  /*Write out final result*/
  cout <<"End Main, complete run "<<endl;
  
}//end main
