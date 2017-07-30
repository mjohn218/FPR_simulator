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

void associate_lipid_sigmaPBCCELL(int p1, int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom) 
{

	/*Move the interfaces that are reacting with one another to
	 the position geometrically between them, but move each
	 interface along this vector proportional to their diffusion constant
	 
	 In this version, 
	 Move the two interfaces to be vertically aligned, so the protein is directly above the lipid.
	 Enforce PBC.
	 */

	int prod = Rlist[mu][2];
	int iind = ihome[i1];
	int iind2 = ihome[i2];
	
	/*Store each protein's pre-association complex values*/
	int c1 = bases[p1].mycomplex;
	int s1 = ind_com[c1].mysize;
	int c2 = bases[p2].mycomplex;
	int s2 = ind_com[c2].mysize;
	int newsize = s1 + s2;

	
	double dx, dy, dz;

	double *dtrans = new double[3];
	double *drev = new double[3];

	/*Now move protein*/
	double Dxtot = ind_com[c1].Dx + ind_com[c2].Dx;
	double Dytot = ind_com[c1].Dy + ind_com[c2].Dy;
	double Dztot = ind_com[c1].Dz + ind_com[c2].Dz;

	/*For proteins bound to membrane, Dztot could be zero*/
	double tol = 1E-16;
	if (Dztot < tol)
		Dztot = 1; //otherwise you divide by zero

	
	
	dx = bases[p2].x[iind2] - bases[p1].x[iind];//Opposite sign from zsigmaPBCCELL.cpp (now its p2-p1)
	dx -= plist.xboxl * round(dx / plist.xboxl);
	dy = bases[p2].y[iind2] - bases[p1].y[iind];
	dy -= plist.yboxl * round(dy / plist.yboxl);
	dz = -bases[p2].z[iind2] + bases[p1].z[iind];//This one is still (p1-p2) 
	
	
	/*Eliminate difference in x and y positions, then move to contact in z!!
	 */
	double R1 = abs(dz);//So R1 is separation in z dimension
	double theta;
	double pivx;
	double pivy;
	double pivz;
	double small = 1E-9;
	//distance to move to place interfaces, along vector v
	dtrans[0] = dx * ind_com[c1].Dx / Dxtot;
	dtrans[1] = dy * ind_com[c1].Dy / Dytot;
	dtrans[2] = dz * ind_com[c1].Dz / Dztot * (bindrad / R1 - 1.0);
   

	drev[0] = -dx * ind_com[c2].Dx / Dxtot;
	drev[1] = -dy * ind_com[c2].Dy / Dytot;
	drev[2] = -dz * ind_com[c2].Dz / Dztot * (bindrad / R1 - 1.0); 
	
	if(dz==0){//this shouldn't happen 
	  dtrans[2] = bindrad * ind_com[c1].Dz / Dztot;
	  drev[2] = -bindrad * ind_com[c2].Dz / Dztot;
	}
	/*NOW WE ARE MOVING THE TWO PROTEINS TO A SEPARATION OF SIGMA*/
	//just translate if larger complexes are coming together
	translate_intPBCCELL(p1, c1, ind_com, bases, dtrans, plist);

	/*Same for complex2*********************************/
	//just translate
	translate_intPBCCELL(p2, c2, ind_com, bases, drev, plist);

	cout <<"Initial separation in z (Lipids) :"<<R1<<'\t';
	dx = -bases[p2].x[iind2] + bases[p1].x[iind];
 	dx -= plist.xboxl * round(dx / plist.xboxl);
 	dy = -bases[p2].y[iind2] + bases[p1].y[iind];
 	dy -= plist.yboxl * round(dy / plist.yboxl);
 	dz = -bases[p2].z[iind2] + bases[p1].z[iind]; //This should be zero now!
	
 	double R2 = dx * dx + dy * dy + dz*dz;
	R1 = sqrt(R2);
	cout <<"final sep: "<<R1<<endl;
	cout <<"final crds: "<<bases[p1].xcom<<' '<<bases[p1].ycom<<' '<<bases[p1].zcom<<endl;
	cout <<"final crds: "<<bases[p2].xcom<<' '<<bases[p2].ycom<<' '<<bases[p2].zcom<<endl;

	/*Combine the complexes into one complex by listing all c2 proteins into c1,
	 unless you are closing a loop, which could happen for clathrin cages.
	 Update the complex COM and diffusion constants
	 */

	int mp;
	double totalmassx = 0;
	double totalmassy = 0;
	double totalmassz = 0;
	ind_com[c1].xcom = 0;
	ind_com[c1].ycom = 0;
	ind_com[c1].zcom = 0;
	
	/*Determine if you crashed two clathrins together
	 unless they are in the same complex, in which case you are closing a loop.
	 */
	int cancel;
	if (c1 != c2)
	  cancel = measure_overlap(c1, c2, ind_com, bases, plist.pclath);
	else
	  cancel = 0; //closing a loop, don't measure overlap becasue all proteins are in same complex
	
	if (cancel == 0) {
	  /*Update the status, free and bound lists of these two proteins*/
	  //change status of the interface
	  bases[p1].istatus[iind] = prod;
	  bases[p2].istatus[iind2] = prod;
	  //no longer free
	  int i;
	  for (i = 0; i < bases[p1].nfree; i++) {
	    if (bases[p1].freelist[i] == i1) {
	      bases[p1].freelist[i] = bases[p1].freelist[bases[p1].nfree - 1]; //put last reaction in place of this one
	      i = bases[p1].nfree;
	    }
	  }
	  bases[p1].nfree -= 1;
	  for (i = 0; i < bases[p2].nfree; i++) {
	    if (bases[p2].freelist[i] == i2) {
	      bases[p2].freelist[i] = bases[p2].freelist[bases[p2].nfree - 1]; //put last reaction in place of this one
	      i = bases[p2].nfree;
	    }
	  }
	  bases[p2].nfree -= 1;
	  
	  /*add this as a possible dissociation reaction to both proteins
	    so we know who they are bound to.
	  */
	  bases[p1].bndlist[bases[p1].nbnd] = prod; //put this reaction as last possible for uni
	  bases[p1].nbnd += 1; //now a dissociation reaction is possible
	  
	  bases[p2].bndlist[bases[p2].nbnd] = prod;
	  bases[p2].nbnd += 1; //now a dissociation reaction is possible
	  
	  bases[p1].partner[iind] = p2;
	  bases[p2].partner[iind2] = p1;
	  
	  /*add c2's proteins to c1's list*/
	  /*Get new COM of complex*/
	  /*Copy final complex into the spot of now gone complex c2 */
	  int mp;
	  int j;
	  int tar;
	  int flagbndry = 0;
	  if (c1 == c2) {
	    cout << "CLOSING A LOOP! " << endl;
	    plist.nloop++;
	    if (ind_com[c1].radR * 2.0 > plist.xboxl / 2.0)
	      flagbndry = 1;
	    update_one_com_onlyPBCCELL(c1, ind_com, bases, plist, flagbndry);
	    update_radiusPBCCELL(c1, ind_com, bases, plist);
	    update_diffusion(c1, ind_com, bases);
	    update_rot_diffusion(c1, ind_com, bases, plist.prerot);
	  } else {
	    ind_com[c1].mysize = newsize;
	    int t = s1;
	    for (i = 0; i < s2; i++) {
	      mp = ind_com[c2].plist[i];
	      ind_com[c1].plist[t] = mp; //add these proteins to the first complex
	      t++;
	      bases[mp].mycomplex = c1;
	    }
	    if ((ind_com[c1].radR + ind_com[c2].radR) * 2.0 > plist.xboxl / 2.0)
	      flagbndry = 1;
	    update_one_com_onlyPBCCELL(c1, ind_com, bases, plist, flagbndry);
	    update_radiusPBCCELL(c1, ind_com, bases, plist);
	    update_diffusion(c1, ind_com, bases);
	    update_rot_diffusion(c1, ind_com, bases, plist.prerot);
	    
	    plist.ntotalcomplex -= 1;
	    tar = plist.ntotalcomplex;
	    if (c2 != plist.ntotalcomplex) {
	      
	      /*otherwise, you are just deleting c2 entirely*/
	      //copy element by element
	      ncrosscom[c2] = ncrosscom[tar];
	      ind_com[c2].xcom = ind_com[tar].xcom;
	      ind_com[c2].ycom = ind_com[tar].ycom;
	      ind_com[c2].zcom = ind_com[tar].zcom;
	      
	      ind_com[c2].Dx = ind_com[tar].Dx;
	      ind_com[c2].Dy = ind_com[tar].Dy;
	      ind_com[c2].Dz = ind_com[tar].Dz;
	      
	      ind_com[c2].radR = ind_com[tar].radR;
	      
	      for (j = 0; j < ind_com[tar].mysize; j++) {
		ind_com[c2].plist[j] = ind_com[tar].plist[j];
	      }
	      ind_com[c2].mysize = ind_com[tar].mysize;
	      
	      for (j = 0; j < ind_com[c2].mysize; j++) {
		mp = ind_com[c2].plist[j];
		bases[mp].mycomplex = c2;
	      }
	    }
	  }
	} else {
	  /*IN THIS CASE< CLATHRINS CRASHED TOGETHER
	    un-bind them, allow to diffuse apart*/
	  cout << "Unbind clathrins, crashed together ! " << endl;
	  
	  dtrans[0] *= -1;
	  dtrans[1] *= -1;
	  dtrans[2] *= -1;
	  
	  drev[0] *= -1;
	  drev[1] *= -1;
	  drev[2] *= -1;
	  
	  /*Update the position of proteins in complex one*/
	  translate_intPBCCELL(p1, c1, ind_com, bases, dtrans, plist);
	  
	  /*Update position of proteins in c2*/
	  translate_intPBCCELL(p2, c2, ind_com, bases, drev, plist);
	  
	}
	
	
	//  cout <<"ind_com crds: "<<ind_com[c1].xcom<<' '<<ind_com[c1].ycom<<' '<<ind_com[c1].zcom<<endl;

	

	delete[] dtrans;
	delete[] drev;

}
