/*
  (formerly called rd_ptrajref_reweight2.cpp)
single particle reaction-diffusion for
Reversible Trap problem.
A is fixed at center and uniform concentration of non-interacting
B particles can bind to A, or they reflect off it. 
B particles avoid overlapping with A whether it is free or bound.

this program propagates using free diffusion, and reweights for particles
within the reaction zone.

This program runs multiple trajectories to collect statistics, generating
a new initial configuration for each new run.
  
For association, with all traj. reweighting programs, adds up probability from each B particle to compare
with one URN.

  
summer 2013  
maggie johnson
*/
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "vol_help.h"
#include "rand_gsl.h"
#include "md_timer.h"

#include "reactions.h"
#include "GF_calls.h"
#include "utility_calls.h"
#include "cell_neighbor_lists.h"
#include "2Drelated.h"
#include "Faddeeva.hh"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <iomanip>


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
  //  update_rates(wholep, plist, Rlist,  bindrad, kr, kdiss2, kdiss3, rxtype, clrxn, Kd, Kd2, Kd3, p_home);
  cout <<"free to bind: "<<endl;
  renumber_list(cntrxn[0], freelist);
  cout <<"ready to unbind: "<<endl;
  renumber_list(cntrxn[1], bndlist);
  cout <<"free to mutate: "<<endl;
  if(cntrxn[2]>0)
    renumber_list(cntrxn[2], zlist);
  cout <<"Check reactions ! "<<endl;
  check_reactions(plist, bases, numpartners, Speclist, Nmyrxn, Rlist, myrxn, Nspecies);
  int clathstop;//=wholep[0].ninterface+2;//index of nonbinding clathrin (3 partners already)


  ifstream startfile(argv[6]);
  cout <<"now set protein status: "<<endl;
  set_status(startfile,wholep,bases, plist, Ncopy);
  /*Read in the coordinates*/
  cout <<"now read in coordinates: "<<endl;
  string *names=new string[Nprotypes];
  ifstream crdfile(argv[7]);
  read_coords(crdfile, Nprotypes, wholep, bases, ind_com, Ncopy, names);
  
  
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

  int Nit=int(plist.Nit);
  /*G(r) stuff*/
  int nbins=1000;//int(plist.xboxl/2.0/delr);
  double delr=(plist.xboxl/2.0-bindrad[0])/(1.0*nbins);//0.1;
  nbins+=1;//for zero spot
  cout <<"nbins gr: "<<nbins<<" Ntotal mol: "<<Ntotalmol<<" max g(r): "<<nbins*delr<<endl; 
  int tpts=1;//plist.grwrite;
  //int tpts=100;//Nit/grfreq;
  int grfreq=Nit/tpts;
  cout <<"N gr time pts: "<<tpts<<" Gr time frequency: "<<grfreq<<endl;
  tpts+=1;
  double **gr=new double*[tpts];//
  for(i=0;i<tpts;i++)
    gr[i]=new double[nbins];
  for(i=0;i<tpts;i++)
    init_gr_zero(gr[i], nbins);
  

  //  double *gr=new double[nbins];

  // for(i=0;i<Nprotypes*Nprotypes;i++)
//     gr[i]=new double[nbins];

//   double **Vofr=new double *[Ntotalmol];
//   for(i=0;i<Ntotalmol;i++)
//     Vofr[i]=new double[nbins];
//   for(i=0;i<Ntotalmol;i++){
//     for(j=0;j<nbins;j++)
//       Vofr[i][j]=0;
//   }
//   cout<<" vofr: "<<Vofr[0][0]<<endl;
  ofstream gaafile("gA_B.dat");
  ofstream restart("restartsW.out");
  //  ofstream gabfile("gab.dat");
  //ofstream gbbfile("gbb.dat");

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
  double xboxl=plist.xboxl;
  double yboxl=plist.yboxl;
  double zboxl=plist.zboxl;//nm
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


  /***************************/
  /*Begin RD simulation*/
  cout <<"deltat: "<<plist.deltat<<endl;
  double deltat=plist.deltat;
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

   double rho0=(Ntotalmol-1)*1.0/(xboxl*yboxl-M_PI*bindrad[0]*bindrad[0]);//2D
//   calc_gr_center(plist, delr, nbins, bases,  Ncopy, gr, rho0);
   cout<< " Ntot complex "<<plist.ntotalcomplex<<" RHO 0: "<<rho0<<endl;
//   //  gabfile<<"iter: "<<0<<' '<<inst_pro[0][0]<<' '<<inst_pro[1][0]<< endl;
//   //gbbfile<<"iter: "<<0<<' '<<inst_pro[1][0]<<endl;
//   write_gofr( nbins, gr, delr, gaafile);
  // write_gofr(0, 1, plist, nbins, gr, delr, gabfile);
  //write_gofr(1, 1, plist, nbins, gr, delr, gbbfile);

  double kpi=4.0*M_PI*bindrad[0];
  double Dtot=wholep[0].Dx+wholep[1].Dx;
  double fourpi=4.0*M_PI;
  
  double kact=kr[0];
  double kdiff=4.0*M_PI*Dtot*bindrad[0];
  double fact=1.0+kr[0]/kdiff;
  cout <<"kact: "<<kact<<" kdiff: "<<kdiff<<" Dtot: "<<Dtot<<" bindrad "<<bindrad[0]<<endl;
  double alpha=fact*sqrt(Dtot)/bindrad[0];
  double alpharef=sqrt(Dtot)/bindrad[0];
  double cof=kr[0]/(kr[0]+kdiff);
  double bexp=alpha*sqrt(deltat);
  double sqfDt=sqrt(4.0*Dtot*deltat);
  /*Update all the rates to be ka and kb, rather than kon and koff*/

  cout <<"activation rate in nm^3/us: "<<kr[0]<<" kd: "<<kr[1]<<" Kd, uM: "<<kr[1]/kr[0]/1E6/6.022E23*1E24*1E6<<endl;
  cout <<"Dtot: "<<Dtot<<endl;
  //  check_bndlist(plist.ntotalcomplex, ind_com, bases);
  double R1;
  double Rmax1=3.5*sqrt(4.0*Dtot*deltat)+bindrad[0];
  double Rmax2D=Rmax1;

  double Dtot2D=Dtot;

  const double RstepSize=sqrt(Dtot*deltat)/50;
  int veclen=int((Rmax2D+RstepSize-bindrad[0])/RstepSize);
  veclen+=2;

  cout <<"Rmax: "<<Rmax2D<<" step size: "<<RstepSize<<" Veclength: "<<veclen<<endl;
  gsl_matrix * msur = gsl_matrix_alloc(2, veclen);
  TBLsur(msur, bindrad[0], Dtot2D, kr[0], deltat, Rmax2D, RstepSize);
  cout<<"TBLsur made successfully"<<endl;
  //FILE * f1 = fopen ("testmsur.dat", "w");
  //gsl_matrix_fprintf (f1, msur,"%f");
  //fclose (f1);
  char sname[300];
  ofstream s1file;
  s1file.precision(12);
  sprintf(sname, "passoc_vs_r_k%g_D%g_dt%g.dat", kr[0], Dtot2D, kr[0]);
  s1file.open(sname);
  for(i=0;i<veclen-2;i++){
    s1file<<bindrad[0]+i*RstepSize<<'\t'<<gsl_matrix_get(msur, 1, i)<<endl;
  }
  gsl_matrix * mnorm = gsl_matrix_alloc(2, veclen);
  TBLnorm(mnorm,  bindrad[0], Dtot2D, kr[0], deltat, Rmax2D, RstepSize);
  cout<<"TBLnorm made successfully"<<endl;
  gsl_matrix * mpir = gsl_matrix_alloc(veclen, veclen);
  TBLpirr_longR(mpir,  bindrad[0], Dtot2D, kr[0], deltat, Rmax2D, RstepSize, Rmax2D);
  cout<<"TBLpirr made successfully"<<endl;


  // FILE * f2 = fopen ("testmpir.dat", "w");
//   gsl_matrix_fprintf (f2, mpir,"%f");
//   fclose (f2);

  FILE * f3 = fopen ("testmnorm.dat", "w");
  gsl_matrix_fprintf (f3, mnorm,"%f");
  fclose (f3);

  i=0;
  
  double Rmaxsq=Rmax2D*Rmax2D;  
  cout <<"Rmaxsq: "<<Rmaxsq<<endl;
  int r0bins=1000;
  double delr0;
  double *pcorrect=new double[r0bins];
  double *pfreea=new double[r0bins];
  int r0ind;
  int *Nsum=new int[Nprotypes];
  Nsum[0]=0;
  for(i=1;i<Nprotypes;i++){
    Nsum[i]=Nsum[i-1]+Ncopy[i-1];
  }
  
  double *prevnorm=new double[Ntotalmol];
  int *previter=new int[Ntotalmol];
  double *ps_prev=new double[Ntotalmol];
  double *prevsep=new double[Ntotalmol];
  
  int size;
  int place, place2;
  double rtol=1E-10;
  int flag;
  int p1, p2;
  double probvec1, p0_ratio, currnorm;

  for(i=0;i<Ntotalmol;i++){
    prevnorm[i]=1.0;
    previter[i]=0;
    ps_prev[i]=0;
    prevsep[i]=0;
    
  }
  /*For multi protein systems this will have to be calculated on-the-fly, because
    Dtot will change.*/
  //  errormodel( deltat, Dtot,  bindrad[0],  kr[0],delr0,  r0bins, pcorrect, Rmax, pfreea);
  cout <<"Rmax1: "<<Rmax1<<" Dtot: "<<Dtot<<endl;
  cout <<"IN GET DISTANCE, COMPARING RMAX TO THE ACTUAL DISTANCE R1, NOT JUST THE SEPARATION! "<<endl;
  cout <<"error model "<<endl;
  cout <<"distance is separation beyond bind rad! "<<endl;
  

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


  int flagsep;
  double memfact=1.0;//If you are binding to PIP2, this should be 1/2.
  int pp, qq, nb, hh;
  
  
  
  double maxrad=150;
  int radbins=3000;
  double delrad=maxrad/(1.0*radbins);
  int ind_rad;
  double rad2, rad;

  
  int Nit1=Nit+1;
  int *phist=new int[Nit1];

  
  for(j=0;j<Nit1;j++){
    phist[j]=0;
  }

  int Ncsave=plist.ntotalcomplex;
  double x0;
  
  double currx, curry, currz;
  int rep;
  int Nrep=plist.Nrep;
  double d2;
  int ncurr=0;
  int bflag=1;
  int bit=0;
  int noverlap=0;
  double pbindr2=bindrad[0]*bindrad[0];
  double stretch;
  int *olist=new int[Ntotalmol];
  char tname[200];
  ofstream probfile;
  double phi;
  sprintf(tname, "prob_deltat%g.dat", deltat);
  probfile.open(tname);
  double addx, addy, addz;
  
  double axe2;
  double R, Ry;
  double small=1E-8;
  //  double alphacurr=alpha;
  double totprob=0;
  double pnow;
  for(rep=0;rep<Nrep;rep++){
    /*start proteins separated by sigma */
    cout <<"Repeat: "<<rep<<endl;
    //alphacurr=alpha;
    for(i=0;i<Ntotalmol;i++){
      prevnorm[i]=1.0;
      previter[i]=-1;
      ps_prev[i]=0;
      prevsep[i]=0;
    }
    
    plist.ntotalcomplex=Ncsave;
    int zeroB=1;//force B to zero
    generate_initial_crds_ABavoid( plist,  bases, Ncopy,  ind_com, bindrad,  zeroB);
    gen_copy_coords( Nprotypes,  wholep, bases,  ind_com, Ncopy);
    
    rho0=(plist.ntotalcomplex-1)*1.0/(plist.xboxl*plist.yboxl-M_PI*bindrad[0]*bindrad[0]);//2D
    calc_gr_center(plist, delr, nbins, bases,  gr[0], bindrad[0]);
    
    
    for(i=0;i<Ntotalmol;i++){
      bases[i].nfree=1;
      bases[i].nbnd=0;
      
    }
    for(it=1;it<Nit+1;it++){
      
      ncross[0]=0;
      
      
      /*test dissociation*/
      if(bases[0].nbnd==1){
	rate=kr[1];
	ppart=bases[0].partner[0];
	prob=1-exp(-rate*deltat*us_to_s);
	rnum=1.0*rand_gsl();
	if(prob>rnum){
	  rnum2=rnum+rand_gsl()*irandmax;
	  if(prob>rnum2){
	    //cout <<"Dissociate at iter: "<<it<<" protein: "<<i<<" partner: "<<ppart<<" randomnum: "<<rnum2<<endl;
	    
	    //move the two proteins to binding radius
	    //cancel=break_complex2(i, mu, j, bases, Rlist, i_home, ind_com, plist, bindrad, inst_pro, ppart, i1, i2, p_home, myrxn);
	    bases[0].nbnd=0;
	    bases[0].nfree=1;
	    
	    
	    ncross[0]=-1;
	    
	    plist.ntotalcomplex++;
	    R=bindrad[0];
	    /*To sample uniformly on the sphere, choose
	      polar angle from -1:1, uniform in cos(theta),
	      and azimuthal uniform in 0:2pi.
	    */
	    
	    phi=rand_gsl()*2.0*M_PI;//between 0 and 2pi
	    addx=R*cos(phi);//sin(theta)=1
	    if(addx>0)addx+=small;
	    else addx-=small;//so final separation is not binding radius with precision issues
	    addy=R*sin(phi);//sin(theta)=1
	    addz=0;//cos(theta)=0
	    bases[ppart].xcom=addx;
	    bases[ppart].ycom=addy;
	    bases[ppart].zcom=addz;
	    prevnorm[ppart]=1.0;
	    // bases[ppart].Dx=Dtot;
	    // 		bases[ppart].Dy=Dtot;
	    // 		bases[ppart].Dz=Dtot;
	    
	  }
	}
      }//Done testing Dissociation
      
      
      
      
      i=0;
      /*Measure overlap, even if A is bound*/
      totprob=0;
      /*Test bimolecular reactions!*/
      nfree=bases[0].nfree;
      mu=0;
      
      /*Now test distance for bimolecular reactions*/
      if(nfree==1){
	if(ncross[0]==0){
	  
	  /*Calc probability of associating*/
	  for(j=1;j<Ntotalmol;j++){
	    r2=bases[j].xcom*bases[j].xcom+bases[j].ycom*bases[j].ycom+bases[j].zcom*bases[j].zcom;
	    if(r2<Rmaxsq){
	      /*If A is free, calculate association probability, otherwise calculate trajectory probability*/
	      
	      nc1=ncross[0];
	      
	      R1=sqrt(r2);
	      sep=R1-bindrad[mu];
	      ratio=bindrad[mu]/R1;
	      if(sep<0){
		cout <<"separation <0: "<<sep<<" r1 "<<R1<<" p1: "<<i<<" p2: "<<j<<" Rep: "<<rep<<" it "<<it<<endl;
		cout <<"Current coordinates: "<<endl;
		cout <<bases[j].xcom<<' '<<bases[j].ycom<<' '<<bases[j].zcom<<endl;
		sep=0;
		ratio=1;
		R1=bindrad[mu];
		
	      }
	      	      
	      currnorm=1.0;
	      p0_ratio=1.0;
	      if(previter[j]==(it-1)){
		/*The next step after dissociation, this should be for a reflectin move, which you performed after dissociation*/
		p0_ratio = DDpirr_pfree_ratio_ps(mpir, msur, mnorm, R1,  Dtot,  deltat, prevsep[j],  ps_prev[j],  rtol,  bindrad[mu]); //pirr_pfree_ratio_ps(R1, Dtot, deltat, prevsep[proa][s], pirrval, pnormval, ps_prev[proa][s], rtol);
		//currnorm=prevnorm[proa][s]*p0_ratio;
		currnorm=prevnorm[j]*p0_ratio;
		// if(p0_ratio<1E-3)
// 		  cout <<"p0ratio: "<<p0_ratio<<" currnorm: "<<currnorm<<" particle: "<<j<<" iter: "<<it<< " rep: "<<rep<<endl;
// 		if(p0_ratio>30)
// 		  cout <<"p0ratio: "<<p0_ratio<<" currnorm: "<<currnorm<<" particle: "<<j<<" iter: "<<it<< " rep: "<<rep<<endl;
		
// 		if(currnorm<1E-3)
// 		  cout <<"p0ratio: "<<p0_ratio<<" currnorm: "<<currnorm<<" particle: "<<j<<" iter: "<<it<< " rep: "<<rep<<endl;
		
// 		if(currnorm>40)
// 		  cout <<"p0ratio: "<<p0_ratio<<" currnorm: "<<currnorm<<" particle: "<<j<<" iter: "<<it<< " rep: "<<rep<<endl;
	      }
	      
	      probvec1 = DDpsur(msur, Dtot, deltat, R1, bindrad[mu]);
	      probvec[0][nc1]=probvec1*currnorm;
	      crosspart[0][nc1]=j;
	      
	      totprob+=probvec[0][nc1];
	      //probvec[j][nc2]=probvec[i][nc1];
	      prevsep[j]=R1;
	      previter[j]=it;
	      prevnorm[j]=currnorm;
	      ps_prev[j]=1.0-probvec1*currnorm;
	      ncross[0]++;
	      
	    }else{
	      //just undergo free diffusion
	      dx=sqrt(2.0*deltat*bases[j].Dx)*GaussV();
	      dy=sqrt(2.0*deltat*bases[j].Dy)*GaussV();
	      dz=0.0;//sqrt(2.0*deltat*bases[j].Dz)*GaussV();
	      
	      bases[j].xcom+=dx;
	      bases[j].ycom+=dy;
	      bases[j].zcom+=dz;
	      /*put inside of box*/
	      flag=0;
	      xtot=0.0;
	      ytot=0.0;
	      ztot=0.0;
	      xchg=bases[j].xcom-xboxl/2.0;
	      if((xchg)>0){
		xtot=-(xchg);//shift x coordinates back
		flag=1;
	      }else if((xchg)<-xboxl){
		xtot=-(bases[j].xcom+xboxl/2.0);
		flag=1;
	      }
	      ychg=bases[j].ycom-yboxl/2.0;
	      if((ychg)>0){
		ytot=-(ychg);//shift x coordinates back
		flag=1;
	      }else if((ychg)<-yboxl){
		ytot=-(bases[j].ycom+yboxl/2.0);
		flag=1;
	      }
	      // zchg=bases[j].zcom-zboxl/2.0;
// 	      if((zchg)>0){
// 		ztot=-(zchg);//shift x coordinates back
// 		flag=1;
// 	      }else if((zchg)<-zboxl){
// 		ztot=-(bases[j].zcom+zboxl/2.0);
// 		flag=1;
// 	      }
	      if(flag==1){
		/*Put back inside the box*/
		bases[j].xcom+=2.0*xtot;
		bases[j].ycom+=2.0*ytot;
		//bases[j].zcom+=2.0*ztot;
	      }
	      
	      
	    }//close enough to evaluate
	  }//done looping over all other proteins
	  /*you just performed an association attempt*/
	  //	  alphacurr=alpha;
	}else{
	  /*You just dissociated, so avoid overlap, but don't try to associate
	    calculate ptraj of reflecting BCs
	  */
	  totprob=0;
	  for(j=1;j<Ntotalmol;j++){
	    r2=bases[j].xcom*bases[j].xcom+bases[j].ycom*bases[j].ycom+bases[j].zcom*bases[j].zcom;
	    if(r2<Rmaxsq){
	      /*If A is free, calculate association probability, otherwise calculate trajectory probability*/
	      if(j!=ppart){
		/*the protein that just dissociated should not be moved again*/
		/*Move particles and avoid overlap*/
		p1=j;
		dx=bases[p1].xcom+sqrt(2.0*deltat*bases[p1].Dx)*GaussV();
		dy=bases[p1].ycom+sqrt(2.0*deltat*bases[p1].Dy)*GaussV();
		dz=0;//bases[p1].zcom+sqrt(2.0*deltat*bases[p1].Dz)*GaussV();
		
		r2=dx*dx+dy*dy+dz*dz;
		while(r2<bindrad[0]*bindrad[0]){
		  dx=bases[p1].xcom+sqrt(2.0*deltat*bases[p1].Dx)*GaussV();
		  dy=bases[p1].ycom+sqrt(2.0*deltat*bases[p1].Dy)*GaussV();
		  dz=0;//bases[p1].zcom+sqrt(2.0*deltat*bases[p1].Dz)*GaussV();
		  r2=dx*dx+dy*dy+dz*dz;
		}
		/*move these proteins and avoid overlap*/
		bases[p1].xcom=dx;
		bases[p1].ycom=dy;
		bases[p1].zcom=dz;
		
	      }
	    }else{
	      //just undergo free diffusion
	      dx=sqrt(2.0*deltat*bases[j].Dx)*GaussV();
	      dy=sqrt(2.0*deltat*bases[j].Dy)*GaussV();
	      dz=0;//sqrt(2.0*deltat*bases[j].Dz)*GaussV();
	      
	      bases[j].xcom+=dx;
	      bases[j].ycom+=dy;
	      bases[j].zcom+=dz;
	      /*put inside of box*/
	      flag=0;
	      xtot=0.0;
	      ytot=0.0;
	      ztot=0.0;
	      
	      xchg=bases[j].xcom-xboxl/2.0;
	      if((xchg)>0){
		xtot=-(xchg);//shift x coordinates back
		flag=1;
	      }else if((xchg)<-xboxl){
		xtot=-(bases[j].xcom+xboxl/2.0);
		flag=1;
	      }
	      ychg=bases[j].ycom-yboxl/2.0;
	      if((ychg)>0){
		ytot=-(ychg);//shift x coordinates back
		flag=1;
	      }else if((ychg)<-yboxl){
		ytot=-(bases[j].ycom+yboxl/2.0);
		flag=1;
	      }
	      // zchg=bases[j].zcom-zboxl/2.0;
// 	      if((zchg)>0){
// 		ztot=-(zchg);//shift x coordinates back
// 		flag=1;
// 	      }else if((zchg)<-zboxl){
// 		ztot=-(bases[j].zcom+zboxl/2.0);
// 		flag=1;
// 	      }
	      if(flag==1){
		/*Put back inside the box*/
		bases[j].xcom+=2.0*xtot;
		bases[j].ycom+=2.0*ytot;
		//bases[j].zcom+=2.0*ztot;
	      }
	      
	    }//close enough to evaluate
	  }//done looping over all other proteins
	  //	  alphacurr=alpharef;
	}//whether A just dissociated or not
      }else{
	/*The A particle is bound, so update positions that reflect and correct ptraj
	  if they are within the binding radius.
	*/
	ncross[0]=0;
	totprob=0;
	for(j=1;j<Ntotalmol;j++){
	  r2=bases[j].xcom*bases[j].xcom+bases[j].ycom*bases[j].ycom+bases[j].zcom*bases[j].zcom;
	  if(r2<Rmaxsq){
	    /*If A is free, calculate association probability, otherwise calculate trajectory probability*/
	    if(j!=bases[0].partner[0]){
	      /*don't try to move the bound partner*/
	      /*Move particles and avoid overlap*/
	      p1=j;
	      dx=bases[p1].xcom+sqrt(2.0*deltat*bases[p1].Dx)*GaussV();
	      dy=bases[p1].ycom+sqrt(2.0*deltat*bases[p1].Dy)*GaussV();
	      dz=0;//bases[p1].zcom+sqrt(2.0*deltat*bases[p1].Dz)*GaussV();
		
	      r2=dx*dx+dy*dy+dz*dz;
	      while(r2<bindrad[0]*bindrad[0]){
		dx=bases[p1].xcom+sqrt(2.0*deltat*bases[p1].Dx)*GaussV();
		dy=bases[p1].ycom+sqrt(2.0*deltat*bases[p1].Dy)*GaussV();
		dz=0;//bases[p1].zcom+sqrt(2.0*deltat*bases[p1].Dz)*GaussV();
		r2=dx*dx+dy*dy+dz*dz;
	      }
	      /*move these proteins and avoid overlap*/
	      bases[p1].xcom=dx;
	      bases[p1].ycom=dy;
	      bases[p1].zcom=dz;
	      
	    }
	    
	  }else{
	    //just undergo free diffusion
	    dx=sqrt(2.0*deltat*bases[j].Dx)*GaussV();
	    dy=sqrt(2.0*deltat*bases[j].Dy)*GaussV();
	    dz=0;//sqrt(2.0*deltat*bases[j].Dz)*GaussV();
	    
	    bases[j].xcom+=dx;
	    bases[j].ycom+=dy;
	    bases[j].zcom+=dz;
	    
	    /*put inside of box*/
	    flag=0;
	    xtot=0.0;
	    ytot=0.0;
	    ztot=0.0;
	      
	    xchg=bases[j].xcom-xboxl/2.0;
	    if((xchg)>0){
	      xtot=-(xchg);//shift x coordinates back
	      flag=1;
	    }else if((xchg)<-xboxl){
	      xtot=-(bases[j].xcom+xboxl/2.0);
	      flag=1;
	    }
	    ychg=bases[j].ycom-yboxl/2.0;
	    if((ychg)>0){
	      ytot=-(ychg);//shift x coordinates back
	      flag=1;
	    }else if((ychg)<-yboxl){
	      ytot=-(bases[j].ycom+yboxl/2.0);
	      flag=1;
	    }
	    // zchg=bases[j].zcom-zboxl/2.0;
// 	    if((zchg)>0){
// 	      ztot=-(zchg);//shift x coordinates back
// 	      flag=1;
// 	    }else if((zchg)<-zboxl){
// 	      ztot=-(bases[j].zcom+zboxl/2.0);
// 	      flag=1;
// 	    }
	    if(flag==1){
	      /*Put back inside the box*/
	      bases[j].xcom+=2.0*xtot;
	      bases[j].ycom+=2.0*ytot;
	      //  bases[j].zcom+=2.0*ztot;
	    }
	    
	  }//close enough to evaluate
	}//done looping over all other proteins
	/*If A is bound, your move will be a reflecting move, even if it associated in this step, except
	  for the particle that associated*/
	//	alphacurr=alpharef;
      }//if free or bound	  
      
      
      if(ncross[0]>0){
	
	/*might perform this reaction, depending on k_associate*/
	//cout <<"ncross: " <<ncross[0]<<" iter: "<<it<<" rep: "<<rep<<endl;
	for(j=0;j<ncross[0];j++){
	  
	  rnum=1.0*rand_gsl();
	  if(probvec[0][j]>rnum){
	    rnum2=rnum+rand_gsl()*irandmax;//refine the random number so we don't get exactly zero!
	    if(probvec[0][j]>rnum2){
	      // i=0;
	      // 	    pnow=probvec[0][i];
	      // 	    while(rnum2>pnow){
	      // 	      i++;
	      // 	      pnow+=probvec[0][i];
	      // 	    }
	    
	      /*Make sure you choose the first thing that happens*/
	      p2=crosspart[0][j];  
	      
	      //cout <<"Associate between proteins: p1 "<<0<<' '<<p2<< " pact; "<<probvec[0][ci1]<<" rnum: "<<rnum2<<" iter: "<<it<<endl;
	      
	      /*Then perform association*/
	      bases[p2].xcom=0;
	      bases[p2].ycom=0;
	      bases[p2].zcom=0;
	      // bases[p2].Dx=0;
	      // 		bases[p2].Dy=0;
	      // 		bases[p2].Dz=0;
	      bases[0].nbnd=1;
	      bases[0].nfree=0;
	      bases[0].partner[0]=p2;
	      plist.ntotalcomplex-=1;
	      /*Move all the other particles that didn't associate*/
	      for(i=j+1;i<ncross[0];i++){
		p1=crosspart[0][i];
		probvec[0][i]=0;
		if(p1!=p2){
		  //don't try to move protein p2 that just associated
		  dx=bases[p1].xcom+sqrt(2.0*deltat*bases[p1].Dx)*GaussV();
		  dy=bases[p1].ycom+sqrt(2.0*deltat*bases[p1].Dy)*GaussV();
		  dz=0;//bases[p1].zcom+sqrt(2.0*deltat*bases[p1].Dz)*GaussV();
		  
		  r2=dx*dx+dy*dy+dz*dz;
		  while(r2<bindrad[0]*bindrad[0]){
		    dx=bases[p1].xcom+sqrt(2.0*deltat*bases[p1].Dx)*GaussV();
		    dy=bases[p1].ycom+sqrt(2.0*deltat*bases[p1].Dy)*GaussV();
		    dz=0;//bases[p1].zcom+sqrt(2.0*deltat*bases[p1].Dz)*GaussV();
		    r2=dx*dx+dy*dy+dz*dz;
		  }
		  /*move these proteins and avoid overlap*/
		  bases[p1].xcom=dx;
		  bases[p1].ycom=dy;
		  bases[p1].zcom=dz;
		}
	      }//moved the rest of the proteins.
	      j=ncross[0];
	    }else{
	      //don't associate
	      /*put all crossed particles to avoid the center*/
	      //for(i=0;i<ncross[0];i++){
	      p1=crosspart[0][j];
	      dx=bases[p1].xcom+sqrt(2.0*deltat*bases[p1].Dx)*GaussV();
	      dy=bases[p1].ycom+sqrt(2.0*deltat*bases[p1].Dy)*GaussV();
	      dz=0;//bases[p1].zcom+sqrt(2.0*deltat*bases[p1].Dz)*GaussV();
	      
	      r2=dx*dx+dy*dy+dz*dz;
	      while(r2<bindrad[0]*bindrad[0]){
		dx=bases[p1].xcom+sqrt(2.0*deltat*bases[p1].Dx)*GaussV();
		dy=bases[p1].ycom+sqrt(2.0*deltat*bases[p1].Dy)*GaussV();
		dz=0;//bases[p1].zcom+sqrt(2.0*deltat*bases[p1].Dz)*GaussV();
		r2=dx*dx+dy*dy+dz*dz;
	      }
	      /*move these proteins and avoid overlap*/
	      bases[p1].xcom=dx;
	      bases[p1].ycom=dy;
	      bases[p1].zcom=dz;
	      
	    }
	    
	    
	  }else {
	    //don't associate
	    /*put all crossed particles to avoid the center*/
	    //for(i=0;i<ncross[0];i++){
	    p1=crosspart[0][j];
	    dx=bases[p1].xcom+sqrt(2.0*deltat*bases[p1].Dx)*GaussV();
	    dy=bases[p1].ycom+sqrt(2.0*deltat*bases[p1].Dy)*GaussV();
	    dz=0;//bases[p1].zcom+sqrt(2.0*deltat*bases[p1].Dz)*GaussV();
	    
	    r2=dx*dx+dy*dy+dz*dz;
	    while(r2<bindrad[0]*bindrad[0]){
	      dx=bases[p1].xcom+sqrt(2.0*deltat*bases[p1].Dx)*GaussV();
	      dy=bases[p1].ycom+sqrt(2.0*deltat*bases[p1].Dy)*GaussV();
	      dz=0;//bases[p1].zcom+sqrt(2.0*deltat*bases[p1].Dz)*GaussV();
	      r2=dx*dx+dy*dy+dz*dz;
	    }
	    /*move these proteins and avoid overlap*/
	    bases[p1].xcom=dx;
	    bases[p1].ycom=dy;
	    bases[p1].zcom=dz;
	    
	  }
	  
	}//no dissociation
      }//no B particles were close enough to A
      
      /*If no one crossed A, then they were all freely diffused
	if they did cross A, they were either associated, or reflected.
      */
      /*For each time step, keep track of whether the protein associated or not */	  
      phist[it]+=bases[0].nfree;//either add 0 or 1, eventually will divide by the number of repetitions
      
    }//end iterating over time steps
    
    if(rep%10==0){
      probfile.open(tname);
      probfile<<0<<'\t'<<rep+1<<endl;//All trajectories are started out with A free
      for(i=1;i<Nit+1;i++){
	probfile <<i*deltat<<'\t'<<phist[i]<<endl;
      }
      probfile.close();
      
    }
  
    
  }//end looping over repetitions
  ofstream grfile("grA_center.dat");
  grfile<<"000"<<'\t'<<0<<'\t';
  for(j=1;j<nbins;j++)
    grfile<<(j-0.5)*delr<<'\t';
  grfile<<endl;
  int flagdim=2;
  for(i=0;i<tpts;i++){
    gr_volume_norm(gr[i], nbins, delr, rho0, Nrep, bindrad[0], flagdim);
    grfile<<i*grfreq*deltat<<'\t';
    for(j=0;j<nbins;j++)
      grfile<<gr[i][j]<<'\t';
    grfile<<endl;
  }       

  probfile.open(tname);
  probfile<<0<<'\t'<<rep<<endl;//All trajectories are started out with A free
  for(i=1;i<Nit+1;i++){
    probfile <<i*deltat<<'\t'<<phist[i]<<endl;
  }
  probfile.close();
  

  
  stop_timer(&totaltime);
  cout <<timer_duration(totaltime)<<" total time "<<endl;  
  /*Write out final result*/
  cout <<"End Main, complete run "<<endl;
  
}//end main
