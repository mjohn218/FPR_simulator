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

void sweep_separation_allPBC(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, int **Rlist) {
	/*This routine checks whether protein p1 is overlapping any partners in its reaction 
	 zone at its new position that is given by its current position +traj. If it does
	 overlap, the displacement by traj is rejected and a new position for itself and any overlapping partners are selected. Once it
	 no longer overlaps anyone, this protein and its complex are moved and the partners
	 retain their stored new displacements. If a protein has already updated its position (done in sequential order)
	 then it cannot resample a new position, the current protein must still continue to avoid overlapping, however.
	 Enforce PBC in distance measurements.
	 */

	int i, j;

	int k1, k2;
	double dx1, dx2;
	double dy1, dy2;
	double dz1, dz2;
	int p2;
	int i1, i2;
	int rxn;
	int iind1, iind2;
	double dr2;
	int ilist[ncross[p1]];
	int olaplist[ncross[p1]];
	k1 = bases[p1].mycomplex;

	/*The sampled displacement for p1 is stored in traj. the position from the 
	 previous step is still stored in bases[p1].xcom, etc, and will be updated
	 at the end of this routine*/

	double df1, df2, df3;
	/*figure out i2*/
	for (i = 0; i < ncross[p1]; i++) {
		p2 = cross_part[p1][i];
		k2 = bases[p2].mycomplex;
		i1 = cross_int[p1][i];
		rxn = cross_rxn[p1][i];
		if (Rlist[rxn][0] == i1)
			i2 = Rlist[rxn][1]; //this is the interface we're looking for on the other protein
		else
			i2 = Rlist[rxn][0];

		ilist[i] = i2;
	}
	int it = 0;
	int maxit = 50;
	int saveit = 0;
	int t = 0;
	int flag = 0;
	int tsave = 0;
	int flags = 0;
	while (it < maxit) {
		t = 0;
		flag = 0;
		for (i = 0; i < ncross[p1]; i++) {
			p2 = cross_part[p1][i];
			k2 = bases[p2].mycomplex;
			i1 = cross_int[p1][i];
			rxn = cross_rxn[p1][i];
			i2 = ilist[i];
			iind1 = ihome[i1];
			iind2 = ihome[i2];
			dx1 = traj[k1][0] + bases[p1].x[iind1];
			dy1 = traj[k1][1] + bases[p1].y[iind1];
			dz1 = traj[k1][2] + bases[p1].z[iind1];

			dx2 = traj[k2][0] + bases[p2].x[iind2];
			dy2 = traj[k2][1] + bases[p2].y[iind2];
			dz2 = traj[k2][2] + bases[p2].z[iind2];

			df1 = dx1 - dx2;
			df2 = dy1 - dy2;
			df3 = dz1 - dz2;
			df1 -= plist.xboxl * round(df1 / plist.xboxl);
			df2 -= plist.yboxl * round(df2 / plist.yboxl);
			df3 -= plist.zboxl * round(df3 / plist.zboxl);
			dr2 = df1 * df1 + df2 * df2 + df3 * df3;

			if (dr2 < bindrad[rxn] * bindrad[rxn]) {
				/*reselect positions for protein p2*/
//    cout<<it<<' '<<sqrt(dr2)<<endl;
				olaplist[t] = p2;
				t++;
				flag = 1;
			}
		}
		/*Now resample positions of p1 and olaplist, if t>0, otherwise no overlap, so
		 break from loop*/
		if (flag == 1) {
			it++;
			traj[k1][0] = sqrt(2.0 * deltat * ind_com[k1].Dx) * GaussV();
			traj[k1][1] = sqrt(2.0 * deltat * ind_com[k1].Dy) * GaussV();
			traj[k1][2] = sqrt(2.0 * deltat * ind_com[k1].Dz) * GaussV();
//	  cout<<'p1'<<' '<<traj[k1][0]<<' '<<traj[k1][1]<<' '<<ind_com[k1].Dx<<' '<<ind_com[k1].Dy<<endl;
			for (j = 0; j < t; j++) {
				p2 = olaplist[j];
				k2 = bases[p2].mycomplex;
				if (p2 > p1 && movestat[p2] != 2) {
					/*
					 We loop over proteins sequentially, so earlier proteins have already moved and avoided
					 their neighbors and should not be moved again.
					 These new positions selected for proteins not yet moved will be stored and
					 then used when they test for overlap themselves.
					 */

					/*If p2 just dissociated, also don't try to move again*/
					traj[k2][0] = sqrt(2.0 * deltat * ind_com[k2].Dx) * GaussV();
					traj[k2][1] = sqrt(2.0 * deltat * ind_com[k2].Dy) * GaussV();
					traj[k2][2] = sqrt(2.0 * deltat * ind_com[k2].Dz) * GaussV();
//	  cout<<'p2'<<' '<<traj[k2][0]<<' '<<traj[k2][1]<<' '<<ind_com[k2].Dx<<' '<<ind_com[k2].Dy<<endl;
				}
			}
			tsave = t;
		} else {
			saveit = it;
			it = maxit; //break from loop
		}
		if (it == maxit - 1) {
			if (k1 != k2) {
				cout << "can't solve overlap: " << p1 << " max cross: " << ncross[p1] << " n overlap: " << tsave << " pro1: " << olaplist[0] << " D: " << ind_com[k1].Dx << " " << ind_com[k2].Dx << endl;
			}
			//cout <<bases[p1].xcom<<' '<<bases[p1].ycom<<' '<<bases[p1].zcom<<endl;
			// cout <<bases[olaplist[0]].xcom<<' '<<bases[olaplist[0]].ycom<<' '<<bases[olaplist[0]].zcom<<endl;

			flags = 1;
			//exit(1);
		}

	} //end maximum iterations
	/*
	 if(flags==1){
	 //didn't solve, put molecules at contact 
	 dx1=traj[k1][0]+bases[p1].x[iind1];
	 dy1=traj[k1][1]+bases[p1].y[iind1];
	 dz1=traj[k1][2]+bases[p1].z[iind1];
	 
	 dx2=traj[k2][0]+bases[p2].x[iind2];
	 dy2=traj[k2][1]+bases[p2].y[iind2];
	 dz2=traj[k2][2]+bases[p2].z[iind2];
	 
	 df1=dx1-dx2;
	 df2=dy1-dy2;
	 df3=dz1-dz2;
	 df1-=plist.xboxl*round(df1/plist.xboxl);
	 df2-=plist.yboxl*round(df2/plist.yboxl);
	 df3-=plist.zboxl*round(df3/plist.zboxl);
	 
	 dr2=df1*df1+df2*df2+df3*df3;
	 double stretch=bindrad[rxn]/sqrt(dr2)+1E-4*rand_gsl();
	 traj[k1][0]+=df1*(stretch-1.0);
	 traj[k1][1]+=df2*(stretch-1.0);
	 traj[k1][2]+=df3*(stretch-1.0);
	 
	 }
	 */

	int s1 = ind_com[k1].mysize;

	dx1 = traj[k1][0];
	dy1 = traj[k1][1];
	dz1 = traj[k1][2];
	int mp;
	ind_com[k1].xcom += dx1;
	ind_com[k1].ycom += dy1;
	ind_com[k1].zcom += dz1;

	ind_com[k1].xcom -= plist.xboxl * round(ind_com[k1].xcom / plist.xboxl);
	ind_com[k1].ycom -= plist.yboxl * round(ind_com[k1].ycom / plist.yboxl);
	ind_com[k1].zcom -= plist.zboxl * round(ind_com[k1].zcom / plist.zboxl);

	for (i = 0; i < s1; i++) {
		mp = ind_com[k1].plist[i];
		movestat[mp] = 2; //physically updated position of these proteins
		bases[mp].xcom += dx1;
		bases[mp].ycom += dy1;
		bases[mp].zcom += dz1;
		bases[mp].xcom -= plist.xboxl * round(bases[mp].xcom / plist.xboxl);
		bases[mp].ycom -= plist.yboxl * round(bases[mp].ycom / plist.yboxl);
		bases[mp].zcom -= plist.zboxl * round(bases[mp].zcom / plist.zboxl);

		//update interface coords
		for (j = 0; j < bases[mp].ninterface; j++) {
			bases[mp].x[j] += dx1;
			bases[mp].y[j] += dy1;
			bases[mp].z[j] += dz1;
			bases[mp].x[j] -= plist.xboxl * round(bases[mp].x[j] / plist.xboxl);
			bases[mp].y[j] -= plist.yboxl * round(bases[mp].y[j] / plist.yboxl);
			bases[mp].z[j] -= plist.zboxl * round(bases[mp].z[j] / plist.zboxl);

		}
	}
	/*Reset displacements to zero so distance is measured to your current
	 updated position that won't change again this turn 
	 */

	traj[k1][0] = 0;
	traj[k1][1] = 0;
	traj[k1][2] = 0;

}
