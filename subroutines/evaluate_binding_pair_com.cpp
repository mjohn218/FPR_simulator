/*proteins in the same complex are not allowed to bind to each other
due to : line 230 if (myc1 != myc2)  
and same line in 3D section

Also, if proteins are already bound through a set of interfaces, it will not let them bind again.

Set Rmax in 2D to 3.5*sqrt(4Ddt)+sigma!
 */
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_bessel.h>
#include "GF_calls.h"
#include "2Drelated.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <iomanip>
#include <vector>
#include "reactions.h"
#include "evaluate_binding.h"

using namespace std;

void evaluate_binding_pair_com(int i, int j, Fullmol *bases, Complex *ind_com, Parms &plist, Protein *wholep, int &DDtableindex, int *numpartners, int **Speclist, int *i_home, int *ncross, double **probvec, double *TBLID, double *bindrad, int *nprevpart, int **prevlist, int **prevmyface, int **prevpface, double **prevnorm, int *ncurrpart, int **currlist, int **currmyface, int **currpface, double **ps_prev, double **prevsep, double **currprevnorm, double **currps_prev, double **currprevsep, int **myrxn, double **traj, int **crosspart, int **crossint, int **cross_rxn, double *kr, int it, vector<gsl_matrix *>& contsur, vector<gsl_matrix *>& contnorm, vector<gsl_matrix *>&contpir, int MAXALWTBL, int *ncrosscom, int *movestat)
{
  

  /*test to see if proteins i and j interact.
    Currently, if two proteins are already bound together, it will not test to see if they can bind in another set of interfaces.
    !!!
   */
  
  int n;
  int flag = 0;
  int wprot=bases[i].protype;
  for (n = 0; n < wholep[wprot].npropart; n++) {
    if(wholep[wprot].propart[n] == bases[j].protype){
      //&& bases[i].mycomplex!=bases[j].mycomplex){ //THIS WILL PREVENT PROTEINS IN THE SAME COMPLEX FROM TRYING TO BIND
      
      
      flag = 1; //check this protein for binding
      n = wholep[wprot].npropart; //break from loop
    }
  }

  /*If this pair of proteins are already bound together, don't test for binding OR overlap, set flag=0*/
  if(flag==1){
    /*If this pair of proteins are already bound together, don't test for binding OR overlap, set flag=0*/
    if(bases[i].nbnd>0 && bases[j].nbnd>0){
      for(n=0;n<bases[i].nbnd;n++){
	int prod=bases[i].bndlist[n];
	for(int s=0;s<bases[i].ninterface;s++){
	  if(bases[i].istatus[s]==prod && bases[i].partner[s]==j){
	    flag=0;
	    //cout <<"Proteins "<<i<<"  and "<<j <<" are already bound, do not test overlap between them ! Set Flag=0. "<<endl;
	  }
	}
      }
    }
    
  }
  
  int p;
  int i1, i2;
  int np, mu;
  int k;
  int myc1, myc2;
  double Dtot;
  double dx, dy, dz;
  int iind, iind2, ppart;
  double dist2;
  double rleg2, cf;
  double Dr1, Dr2;
  double Rmax;
  int nc1, nc2;
  int idd;
  int tableexistflag=0;
  double ratio, R1;
  int flagsep;
  double sep;
  int proa;
  int myface;
  int pro2;
  int pface;
  int recruitmentflag=0;
  double p0_ratio, currnorm;
  int s, ssave;
  int uniquetableindex=0;
  double deltat=plist.deltat;
  double ktemp;
  int veclen;
  double probvec1;
  double aexp, bexp, alpha, fact, kact, kdiff;
  double rtol=1E-10;
  double fourpi=4.0*M_PI;

  if (flag == 1) {
    /*These two proteins are capable of interacting,
      test to see if those interfaces are free to bind
    */
    
    for (p = 0; p < bases[i].nfree; p++) {
      i1 = bases[i].freelist[p];
      /*test all of i1's binding partners to see whether they are on protein j */
      np = numpartners[i1];
      for (n = 0; n < np; n++) {
	i2 = Speclist[i1][n]; //binding interface partner
	for (k = 0; k < bases[j].nfree; k++) {
	  if (bases[j].freelist[k] == i2) {
	    //both binding interfaces are available!
	    mu = myrxn[i1][n];
	    
	    /*Here now we evaluate the probability of binding*/
	    /*Different interfaces can have different diffusion constants if they
	      also rotate. <theta^2>=6Drdeltat.
	      In that case, <d>=sin(sqrt(6Drdeltat)/2)*2R
	      so add in <d>^2=4R^2sin^2(sqrt(6Drdeltat)/2)
	      for a single clathrin, R=arm length,
	      otherwise R will be from pivot point of rotation, COM,
	      calculate distance from interface to the complex COM
	    */
	    myc1 = bases[i].mycomplex;
	    myc2 = bases[j].mycomplex;
	    Dtot = ind_com[myc1].Dx + ind_com[myc2].Dx;
	    iind = i_home[i1];
	    dx = bases[i].x[iind] - ind_com[myc1].xcom;
	    dy = bases[i].y[iind] - ind_com[myc1].ycom;
	    dz = bases[i].z[iind] - ind_com[myc1].zcom;
	    dx -= plist.xboxl * round(dx / plist.xboxl);
	    dy -= plist.yboxl * round(dy / plist.yboxl);
	    //											dz -= plist.zboxl * round(dz / plist.zboxl);
	    // cout <<" interface : "<<i<<endl;
	    // 			write_crds(bases, i);
	    rleg2 = dx * dx + dy * dy + dz * dz;
	    cf = cos(sqrt(4.0 * ind_com[myc1].Drx * deltat));
	    Dr1 = 2.0 * rleg2 * (1.0 - cf);
	    iind2 = i_home[i2];
	    dx = bases[j].x[iind2] - ind_com[myc2].xcom;
	    dy = bases[j].y[iind2] - ind_com[myc2].ycom;
	    dz = bases[j].z[iind2] - ind_com[myc2].zcom;
	    dx -= plist.xboxl * round(dx / plist.xboxl);
	    dy -= plist.yboxl * round(dy / plist.yboxl);
	    //											dz -= plist.zboxl * round(dz / plist.zboxl);
	    // cout <<" interface j: "<<j<<endl;
	    // 			write_crds(bases, j);
	    
	    rleg2 = dx * dx + dy * dy + dz * dz;
	    cf = cos(sqrt(4.0 * ind_com[myc2].Drx * deltat));
	    Dr2 = 2.0 * rleg2 * (1.0 - cf);
	    
	    //cout <<"Dtrans: "<<Dtot<<" Dr1: "<<Dr1<<" Dr2: "<<Dr2<<" leg: "<<rleg2<<endl;
	    Dtot += (Dr1 + Dr2) / (6.0 * deltat); //add in contributions from rotation
	    /*Reaction zone radius between particle positions*/
	    
	    if (ind_com[myc1].Dz == 0 && ind_com[myc2].Dz == 0) { // This is a reaction on the membrane
	      
	      Rmax = 3.5 * sqrt(4.0 * Dtot * deltat) + bindrad[mu];
	     
	      nc1 = ncross[i];
	      nc2 = ncross[j];
	      //flagsep = get_distancePBCCELL(bases, traj, i, j, deltat, bindrad[mu], ncross, crosspart, crossint, i1, i2, mu, cross_rxn, i_home, it, Rmax, sep, R1, plist);
	      flagsep = get_distancePBCCELL_2D(bases, traj, i, j, deltat, bindrad[mu], ncross, crosspart, crossint, i1, i2, mu, cross_rxn, i_home, it, Rmax, sep, R1, plist, ncrosscom);
	      probvec[i][nc1]=0;//initialize to zero, in case kr[mu]=0
	      probvec[j][nc2]=0;
	      //cout <<" Evaluate 2D binding?: "<<i<<' '<<j<<' '<<R1<<" :R1 "<<Rmax<<" flagsep: "<<flagsep<<" Iter: "<<it<<" kr: "<<kr[mu]<<endl;
	      if(movestat[i]!=2 &&movestat[j]!=2){
		/*This movestat check is if you allow just dissociated proteins to avoid overlap*/
		if (flagsep == 1 && kr[mu]>0) {
		  /*Evaluate probability of reaction, with reweighting*/
		  
		  //Generate 2D tables unless they were not before
		  //Dtemp = ind_com[myc1].Dx + ind_com[myc2].Dx;
		  //Dtemp += (Dr1 + Dr2) / (6.0 * deltat); //add in contributions from rotation
		  ktemp = kr[mu] / 2 / bindrad[mu];
		  
		  for (idd = 0; idd < DDtableindex; idd++) {
		    if (TBLID[idd] == ktemp && TBLID[MAXALWTBL+idd] == Dtot) {
		      tableexistflag = 1;
		      uniquetableindex = idd;
		      break;
		    }
		  }
		  
		  if (tableexistflag == 0) {
		    cout <<"DDtableindex: "<<DDtableindex<<endl;
		    cout <<"Create new 2D table: "<<ktemp<<" D: "<<Dtot<<" Dtot: "<<Dtot<<endl;
		    TBLID[DDtableindex] = ktemp;//first dimension (+0*MAXALWTBL)
		    TBLID[DDtableindex+MAXALWTBL] = Dtot;//second dimension (+1*MAXALWTBL)
		    veclen = sizelookup(bindrad[mu], Dtot, deltat, Rmax);
		    contsur.push_back(gsl_matrix_alloc(2, veclen));
		    contnorm.push_back(gsl_matrix_alloc(2, veclen));
		    contpir.push_back(gsl_matrix_alloc(veclen, veclen));
		    DDmatrixcreate(contsur[DDtableindex], contnorm[DDtableindex], contpir[DDtableindex], bindrad[mu], Dtot, ktemp, deltat, Rmax);
		    uniquetableindex = DDtableindex;
		    DDtableindex += 1;
		    if (DDtableindex == MAXALWTBL) {
		      cout << "You have hit the maximum number of unique 2D reactions allowed: " << MAXALWTBL << endl;
		      cout << "terminating...." << endl;
		      exit(1);
		    }
		  }
		  tableexistflag = 0; //reset
		  
		  ratio = bindrad[mu] / R1;
		  if (sep < 0) {
		    if (myc1 != myc2) {
		      // && i1!=0 && i2!=0) {
		      cout << "separation <0: " << sep << " r1 " << R1 << " p1: " << i << " p2: " << j <<  " it " << it << " i1: " << i1 << " i2: " << i2 << endl;
		      cout << bases[i].xcom << ' ' << bases[i].ycom << ' ' << bases[i].zcom << " nfree: " << bases[i].nfree << endl;
		      cout << bases[j].xcom << ' ' << bases[j].ycom << ' ' << bases[j].zcom << " nfree: " << bases[j].nfree << endl;
		    }else{
		      cout <<"Protein interfaces within a complex, trying to bind. Separation <0: "<<sep<<" r1 " << R1 << " p1: " << i << " p2: " << j <<  " it " << it << " i1: " << i1 << " i2: " << i2 << endl;
		      cout << bases[i].xcom << ' ' << bases[i].ycom << ' ' << bases[i].zcom << " nfree: " << bases[i].nfree << endl;
		      cout << bases[j].xcom << ' ' << bases[j].ycom << ' ' << bases[j].zcom << " nfree: " << bases[j].nfree << endl;
		    }
		    sep = 0;
		    ratio = 1;
		    R1 = bindrad[mu];
		  }
		  
		  currnorm = 1.0;
		  p0_ratio = 1.0;
		  
		  /*protein i is protype wprot and j is wprot2*/
		  
		  proa = i;
		  myface = i1;
		  pro2 = j;
		  pface = i2;
		  if (i > j) {
		    proa = j;
		    myface = i2;
		    pro2 = i;
		    pface = i1;
		  }
		  
		  if (myc1 != myc2) {
		    //recruitmentflag = 0;//SET TO ZERO HERE
		    for (s = 0; s < nprevpart[proa]; s++) {
		      if (prevlist[proa][s] == pro2 && prevmyface[proa][s] == myface && prevpface[proa][s] == pface) {
			
			//recruitmentflag = 0;
			if (prevsep[proa][s] >= Rmax) {
			  p0_ratio=1.0;//BEcause previous reweighting was for 3D, now restart reweighting in 2D.
			  currnorm=1.0;
			  //recruitmentflag = 1; //previous step was a recruitment, new complex formed with PIP2 so we will skip this step
			} else {
			  p0_ratio = DDpirr_pfree_ratio_ps(contpir[uniquetableindex], contsur[uniquetableindex], contnorm[uniquetableindex], R1, Dtot, deltat, prevsep[proa][s], ps_prev[proa][s], rtol, bindrad[mu]);
			  currnorm = prevnorm[proa][s] * p0_ratio;
			  //																	if(isnan(currnorm)){
			}
			ssave = s;
			s = nprevpart[proa];
		      }
		    }
		    /*don't renormalize if their positions are fixed by being in the same complex!*/
		    
		    //if (recruitmentflag != 1) {
		    
		    probvec1 = DDpsur(contsur[uniquetableindex], Dtot, deltat, R1, bindrad[mu]);
		    
		    probvec[i][nc1] = probvec1 * currnorm;
		    probvec[j][nc2] = probvec[i][nc1];
		    
		    /*Store all the reweighting numbers for next step.*/
		    s = ncurrpart[proa];
		    currprevsep[proa][s] = R1;
		    //cout <<"proa: "<<proa<<" s: "<<s<<" pro2: "<<pro2<<" currlist: "<<currlist[proa][s]<<" probvec: "<<probvec1<<" normfact: "<<currnorm<<" R: "<<R1<<endl;
		    
		    currlist[proa][s] = pro2;
		    currmyface[proa][s] = myface;
		    currpface[proa][s] = pface;
		    currprevnorm[proa][s] = currnorm;
		    currps_prev[proa][s] = 1.0 - probvec1 * currnorm;
		    ncurrpart[proa]++;
		  } //Within reaction zone
		  //}
		}
	      }
	    } else {
	      
	      /*3D reaction*/
	      
	      Rmax = 3.0 * sqrt(6.0 * Dtot * deltat) + bindrad[mu];
	      //				cout<<Rmax<<endl;
	      
	      nc1 = ncross[i];
	      nc2 = ncross[j];
	      flagsep = get_distancePBCCELL(bases, traj, i, j, deltat, bindrad[mu], ncross, crosspart, crossint, i1, i2, mu, cross_rxn, i_home, it, Rmax, sep, R1, plist, ncrosscom);
	      probvec[i][nc1]=0;//initialize to zero, in case kr[mu]=0
	      probvec[j][nc2]=0;
	      //cout <<" Evaluate 3D binding?: "<<i<<' '<<j<<' '<<R1<<" :R1 "<<Rmax<<" flagsep: "<<flagsep<<" Iter: "<<it<<" kr[mu]: "<<endl;
	      if(movestat[i]!=2 &&movestat[j]!=2){
		/*This movestat check is if you allow just dissociated proteins to avoid overlap*/
		
		if (flagsep == 1 && kr[mu]>0) {
		  
		  /*Evaluate probability of reaction, with reweighting*/
		  
		  ratio = bindrad[mu] / R1;
		  if (sep < 0) {
		    if (myc1 != myc2 ){
		      //&& i1!=0 && i2!=0) {//only writes out for nonlipids?
		      cout << "separation <0: " << sep << " r1 " << R1 << " p1: " << i << " p2: " << j <<" it " << it << " i1: " << i1 << " i2: " << i2 << endl;
		      cout << bases[i].xcom << ' ' << bases[i].ycom << ' ' << bases[i].zcom << " nfree: " << bases[i].nfree << endl;
		      cout << bases[j].xcom << ' ' << bases[j].ycom << ' ' << bases[j].zcom << " nfree: " << bases[j].nfree << endl;
		    }else{
		      cout <<"Protein interfaces within a complex trying to bind: "<< "separation <0: " << sep << " r1 " << R1 << " p1: " << i << " p2: " << j <<" it " << it << " i1: " << i1 << " i2: " << i2 << endl;
		      cout << bases[i].xcom << ' ' << bases[i].ycom << ' ' << bases[i].zcom << " nfree: " << bases[i].nfree << endl;
		      cout << bases[j].xcom << ' ' << bases[j].ycom << ' ' << bases[j].zcom << " nfree: " << bases[j].nfree << endl;
		      
		    }
		    sep = 0;
		    ratio = 1;
		    R1 = bindrad[mu];
		  }
		  
		  /*If one particle (lipid) is bound in the membrane, reaction prob
		    is half due to flux across only top half of particle.
		  */
		  //memfact = 1.0;
		  //if (bases[i].Dz == 0 || bases[j].Dz == 0)
		  //memfact = 1.0; //0.5;
		  kdiff = fourpi * Dtot * bindrad[mu];// * memfact;
		  
		  kact = kr[mu];
		  
		  /*WHAT IS THE PURPOSE OF THIS BELOW, DOESN'T CHANGE THE FACTOR OF KACT?*/
		  // if ((ind_com[myc1].Dz == 0 && ind_com[myc1].mysize == 2) || (ind_com[myc2].Dz == 0 && ind_com[myc2].mysize == 2)) {
		  //   kact = kact / 1.0; //0.33;
		  // }
		  
		  fact = 1.0 + kact / kdiff;
		  alpha = fact * sqrt(Dtot) / bindrad[mu];
		  aexp = sep / sqrt(4.0 * Dtot * deltat);
		  bexp = alpha * sqrt(deltat);
		  
		  /*Check whether this pair was previously in reaction zone and
		    therefore if reweighting should apply.
		    Only need to store reweighting values for one protein in the pair,
		    but for each pair need to know the protein partner and interface partner
		  */
		  currnorm = 1.0;
		  p0_ratio = 1.0;
		  
		  /*protein i is protype wprot and j is wprot2*/
		  
		  proa = i;
		  myface = i1;
		  pro2 = j;
		  pface = i2;
		  if (i > j) {
		    proa = j;
		    myface = i2;
		    pro2 = i;
		    pface = i1;
		  }
		  /*Find out what the previous reweighting was for this pair, if they
		    were stored from previous step (inside reaction zone).
		  */
		  if (myc1 != myc2) {
		    for (s = 0; s < nprevpart[proa]; s++) {
		      if (prevlist[proa][s] == pro2 && prevmyface[proa][s] == myface && prevpface[proa][s] == pface) {
			p0_ratio = pirr_pfree_ratio_psF(R1, prevsep[proa][s], deltat, Dtot, bindrad[mu], alpha, ps_prev[proa][s], rtol);
			currnorm = prevnorm[proa][s] * p0_ratio;
			
			ssave = s;
			s = nprevpart[proa];
		      }
		    }
		    /*don't renormalize if their positions are fixed by being in the same complex!*/
		  }
		  probvec1 = ratio * kact / (kact + kdiff) * (erfc(aexp) - exp(2.0 * aexp * bexp + bexp * bexp) * erfc(aexp + bexp));
		  //													probvec1 = psur(msur, RstepSize, R1, bindrad[mu]);
		  //cout<<" probvec and randomu num: "<<probvec1<<endl;
		  probvec[i][nc1] = probvec1 * currnorm;
		  probvec[j][nc2] = probvec[i][nc1];
		  
		  /*Store all the reweighting numbers for next step.*/
		  s = ncurrpart[proa];
		  currprevsep[proa][s] = R1;
		  //cout <<"proa: "<<proa<<" s: "<<s<<" pro2: "<<pro2<<" currlist: "<<currlist[proa][s]<<" probvec: "<<probvec1<<" normfact: "<<currnorm<<" R: "<<R1<<endl;
		  
		  currlist[proa][s] = pro2;
		  currmyface[proa][s] = myface;
		  currpface[proa][s] = pface;
		  currprevnorm[proa][s] = currnorm;
		  currps_prev[proa][s] = 1.0 - probvec1 * currnorm;
		  ncurrpart[proa]++;
		} //Within reaction zone
	      }//did not just dissociate
	    }
	  }
	}
      }
    }
  } //These protein partners interact
  


}
