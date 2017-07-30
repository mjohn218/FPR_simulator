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



void all_sweep_separation_complex_rot_memtest_PBCCELL(int &overlaps, int myislandsize, int **neighbors, int islandid, int Ntotalmol, int *islandvec, double deltat, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, double **trajR, double *M, int **Rlist) {
	/*
	  In this version, complex k1 is on the membrane. 
	  If both proteins are on the membrane (Dz==0), evaluate only xy displacement, not z.
		  
	  In this _complex_ version, it tests overlap not just for each protein, but for each complex, so all the proteins in a complex, before performing
	  position updates.

	  This routine checks whether protein p1 is overlapping any partners in its reaction
	  zone at its new position that is given by its current position +traj. If it does
	 overlap, the displacement by traj is rejected and a new position for itself and any overlapping partners are selected. Once it
	 no longer overlaps anyone, this protein and its complex are moved and the partners
	 retain their stored new displacements. If a protein has already updated its position (done in sequential order)
	 then it cannot resample a new position, the current protein must still continue to avoid overlapping, however.
	 Enforce PBC in distance measurements.
	 */
  
	int i, j, jj, existflag=0;
	double *M2 = new double[9];
	double *v = new double[3];
	double *v2 = new double[3];
	int *olaplist = new int[myislandsize];
	int mp, s1;
	int p1, p2, k1, k2, k, proi, proj, ind=0, neighborflag=0, overlapcount = 0, t=0; 
	double x0, dx1, dx2, dr2=0;
	double y0, dy1, dy2;
	double z0, dz1, dz2;

	int *myisland = new int[myislandsize]; //find the molecules in islandid
	for(j = 0; j < Ntotalmol; j++){
		if(islandvec[j]==islandid){
			myisland[ind] = j;
			ind++;
		}
	}
	if(ind!=myislandsize){
		cout<<"islandsize mismatch"<<endl;
	}

	/*The sampled displacement for p1 is stored in traj. the position from the
	 previous step is still stored in bases[p1].xcom, etc, and will be updated
	 at the end of this routine*/

	double df1, df2, df3;
	/*figure out i2*/
	
	int it = 0;
	int maxit = 100;
	int saveit = 0;
	int flag = 0;
	int tsave = 0;
	int flags = 0;

	while (it < maxit) {

		overlapcount=0;

		for(i=0;i<myislandsize-1;i++){
			p1 = myisland[i];

			for(j=i+1;j<myislandsize;j++){
				p2 = myisland[j];
				neighborflag = isneighbor(p1, p2, neighbors);
				if(neighborflag==1){ // p1&p2 are neighbors, check if they overlap
					t = check_overlap(dr2, p1, p2, bases, ind_com, ncross, cross_part, cross_int, cross_rxn, traj, plist, ihome, bindrad, trajR, M, Rlist);

					if(t==1){
//						cout<<it<<" "<<p1<<" "<<p2<<" "<<dr2<<endl;
						existflag = 0;
						for(jj=0;jj<overlapcount;jj++){//if this p1 is not in the existing overlaping protein list, add it
							if(olaplist[jj]==p1){
								existflag = 1;
								jj=overlapcount;//break out the loop
							}
						}
						if(existflag==0){
							olaplist[overlapcount] = p1;
							overlapcount += t;
						}
						existflag = 0;
						for(jj=0;jj<overlapcount;jj++){//if this p2 is not in the existing overlaping protein list, add it
							if(olaplist[jj]==p2){
								existflag = 1;
								jj=overlapcount;
							}
						}
						if(existflag==0){
							olaplist[overlapcount] = p2;
							overlapcount += t;
						}

					}
				}
			}
		}

		/*Now resample positions of p1 and olaplist, if t>0, otherwise no overlap, so
		 break from loop*/
		if (overlapcount > 0) {
			//cout<<"it: "<<it<<" #olap_pro: "<<overlapcount<<" in island "<<islandid<<" with size: "<<myislandsize<<endl;
			it++;
			for(i=0;i<overlapcount;i++){
				p1 = olaplist[i];
				if(movestat[p1] != 2) {
					k1 = bases[p1].mycomplex;
					traj[k1][0] = sqrt(2.0 * deltat * ind_com[k1].Dx) * GaussV();
					traj[k1][1] = sqrt(2.0 * deltat * ind_com[k1].Dy) * GaussV();
					traj[k1][2] = sqrt(2.0 * deltat * ind_com[k1].Dz) * GaussV();
					trajR[k1][0] = sqrt(2.0 * deltat * ind_com[k1].Drx) * GaussV();
					trajR[k1][1] = sqrt(2.0 * deltat * ind_com[k1].Dry) * GaussV();
					trajR[k1][2] = sqrt(2.0 * deltat * ind_com[k1].Drz) * GaussV();

					rotationEuler(trajR[k1][0], trajR[k1][1], trajR[k1][2], M);
					reflect_traj_complex_rad_rotCELL(p1, bases, ind_com, plist.xboxl, plist.yboxl, plist.zboxl, traj, M);

				}
			}
		} else {
			saveit = it;
			it = maxit; //break from loop
		}

		if (it == maxit - 1) {
			cout<<"can't solve overlap islandnum: "<<islandid<<" islandsize: "<<myislandsize<<" number of overlaps: "<<overlapcount<<endl;
			overlaps = overlaps + 1;
		}

	} //end maximum iterations

	for(jj=0;jj<myislandsize;jj++){
		p1 = myisland[jj];
		k1 = bases[p1].mycomplex;
		s1 = ind_com[k1].mysize;

		if(movestat[p1] != 2){ //don't move, don't let your complex partners move

			x0 = ind_com[k1].xcom;
			y0 = ind_com[k1].ycom;
			z0 = ind_com[k1].zcom;

			dx1 = traj[k1][0];
			dy1 = traj[k1][1];
			dz1 = traj[k1][2];

			ind_com[k1].xcom += dx1;
			ind_com[k1].ycom += dy1;
			ind_com[k1].zcom += dz1;

			ind_com[k1].xcom -= plist.xboxl * round(ind_com[k1].xcom / plist.xboxl);
			ind_com[k1].ycom -= plist.yboxl * round(ind_com[k1].ycom / plist.yboxl);
			//ind_com[k1].zcom -= plist.zboxl * round(ind_com[k1].zcom / plist.zboxl);

			for (i = 0; i < s1; i++) {
				mp = ind_com[k1].plist[i];
				movestat[mp] = 2; //physically updated position of these proteins
//				cout<<0<<"\t"<<mp<<endl;
				v[0] = bases[mp].xcom - x0;
				v[1] = bases[mp].ycom - y0;
				v[2] = bases[mp].zcom - z0;
				v[0] -= plist.xboxl * round(v[0] / plist.xboxl);
				v[1] -= plist.yboxl * round(v[1] / plist.yboxl);
				//		v[2] -= plist.zboxl * round(v[2] / plist.zboxl);

				rotate(v, M, v2); //M is one used in last step with no overlap.
				bases[mp].xcom = x0 + v2[0] + traj[k1][0];
				bases[mp].ycom = y0 + v2[1] + traj[k1][1];
				bases[mp].zcom = z0 + v2[2] + traj[k1][2];

				bases[mp].xcom -= plist.xboxl * round(bases[mp].xcom / plist.xboxl);
				bases[mp].ycom -= plist.yboxl * round(bases[mp].ycom / plist.yboxl);
		//		bases[mp].zcom -= plist.zboxl * round(bases[mp].zcom / plist.zboxl);

				//update interface coords
				for (j = 0; j < bases[mp].ninterface; j++) {
					v[0] = bases[mp].x[j] - x0;
					v[1] = bases[mp].y[j] - y0;
					v[2] = bases[mp].z[j] - z0;
					v[0] -= plist.xboxl * round(v[0] / plist.xboxl);
					v[1] -= plist.yboxl * round(v[1] / plist.yboxl);
		//			v[2] -= plist.zboxl * round(v[2] / plist.zboxl);

					rotate(v, M, v2);
					bases[mp].x[j] = x0 + v2[0] + traj[k1][0];
					bases[mp].y[j] = y0 + v2[1] + traj[k1][1];
					bases[mp].z[j] = z0 + v2[2] + traj[k1][2];
					bases[mp].x[j] -= plist.xboxl * round(bases[mp].x[j] / plist.xboxl);
					bases[mp].y[j] -= plist.yboxl * round(bases[mp].y[j] / plist.yboxl);
					//bases[mp].z[j] -= plist.zboxl * round(bases[mp].z[j] / plist.zboxl);

				}
			}
			/*Reset displacements to zero so distance is measured to your current
			  updated position that won't change again this turn
			*/
			traj[k1][0] = 0;
			traj[k1][1] = 0;
			traj[k1][2] = 0;
			trajR[k1][0] = 0;
			trajR[k1][1] = 0;
			trajR[k1][2] = 0;
		}
	}
	
	delete[] M2;
	delete[] v;
	delete[] v2;
	delete[] myisland;
	delete[] olaplist;
}
