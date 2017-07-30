#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstring>
#include <sys/time.h>
#include "reactions.h"
#include "vol_help.h"
#include "rand_gsl.h"
#include "Faddeeva.hh"
#include "utility_calls.h"
#include "vector_rot_calls.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

void generate_initial_crdsN(Parms plist, Fullmol *bases, int *Ncopy, Complex *ind_com, double *bindrad, int Nprotypes, int Nifaces, Protein *wholep, int **Rlist,int *rxtype, int *p_home, int *ihome, std::vector<std::string> &infofilenames) {

	int i1,i2,i,j,jj,ii,k,howmanylipids=0, MAXPR = 10,p,t = 0,mu=0, ct = 0, iface;
	int bflag = 1, Maxit = 50, bit = 0, noverlap = 0,rind, nint=0;
	int theyinteract = 0;
	double coordcont[Nprotypes][Nifaces+1][3];
	double pbindr2, stretch, d2, dx, dy, dz,px,py,pz;
	double dx0,dx1,dy0,dy1,dz0,dz1,d0arm,d1arm;
	double d20,d21,d20arm,d21arm,dx0arm,dx1arm,dy0arm,dy1arm,dz0arm,dz1arm;
	char ctype[3];
	char buffer[32], buffer2[32];
	char *let = new char[MAXPR];
	let[0] = 'A';
	let[1] = 'B';
	let[2] = 'C';
	let[3] = 'D';
	let[4] = 'E';
	let[5] = 'F';
	let[6] = 'G';
	let[7] = 'H';
	let[8] = 'I';
	let[9] = 'J';
	ofstream cfile("coordsALL.out");
	ofstream comfile("coordsCOM.out");

	//extract the structure from .info files
	for(i=0;i<Nprotypes;i++){
		ifstream sfile(infofilenames[i].c_str()); // Use the constructor rather than `open`
		cout<<"Reading molecule info file "<<infofilenames[i]<<endl;
		if (sfile) // Verify that the file was opened successfully
		{
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');

			sfile >> wholep[i].Dx >> wholep[i].Dy >> wholep[i].Dz;
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');
			sfile >> wholep[i].Drx >> wholep[i].Dry >> wholep[i].Drz;
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');
			sfile >> wholep[i].radx >> wholep[i].rady >> wholep[i].radz;
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');
			sfile >> nint; //number of protein i's interfaces
			wholep[i].ninterface = nint;
			wholep[i].nint_write = nint; //this can change below
			cout << "Protein: " << i << " Numinterfaces: " << wholep[i].ninterface << endl;
			cout << "Protein: " << i << " Dx: " << wholep[i].Dx << " Dz: " << wholep[i].Dz << endl;
			cout << "Protein: " << i << " radx: " << wholep[i].radx << " radz: " << wholep[i].radz << endl;
			for (j = 0; j < nint; j++) {
				sfile >> iface;
				cout << " index: " << iface << endl;
				wholep[i].valiface[j] = iface;
				p_home[iface] = i;
				ihome[iface] = j;
			}
			ct += nint;
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');
			if(wholep[i].Dz == 0){
				howmanylipids += Ncopy[i];
			}
			//	initialize structure container to zero
			for(j=0;j<Nifaces+1;j++){
				for (k=0;k<3;k++){
					coordcont[i][j][k] = 0.0;
				}
			}

			for(j=0;j<wholep[i].ninterface+1;j++){
				sfile >> ctype;
				sfile >> coordcont[i][j][0] >> coordcont[i][j][1] >> coordcont[i][j][2];
				cout<< coordcont[i][j][0] << coordcont[i][j][1] << coordcont[i][j][2]<<endl;
				if(j>0){
					for(k=0;k<3;k++){
						coordcont[i][j][k] -= coordcont[i][0][k]; //make all interfaces with respect to the center of mass
					}
				}
			}

			for (j = 0; j < Ncopy[i]; j++) {
				bases[t].ninterface = nint;
				bases[t].protype = i;
				bases[t].nfree = nint;

				bases[t].Dx = wholep[i].Dx;
				bases[t].Dy = wholep[i].Dy;
				bases[t].Dz = wholep[i].Dz;
				bases[t].Drx = wholep[i].Drx;
				bases[t].Dry = wholep[i].Dry;
				bases[t].Drz = wholep[i].Drz;

				bases[t].massx = wholep[i].radx; //plist.mass;
				if (wholep[i].radx == 0)
					bases[t].massx = 1;
				bases[t].massy = wholep[i].rady; //plist.mass;
				if (wholep[i].rady == 0)
					bases[t].massy = 1;
				bases[t].massz = wholep[i].radz; //plist.mass;
				if (wholep[i].radz == 0)
					bases[t].massz = 1;

				bases[t].nbnd = 0;
				bases[t].npartner = 0;
				bases[t].mycomplex = t;

				t++;
			}

		}
		else
		{
			cerr << "File could not be opened!\n"; // Report error
			cerr << "Error code: " << strerror(errno); // Get some info as to why
			exit(1);
		}

	}

	// create random center of masses
	for (i = 0; i < plist.Ntotalmol; i++) {
		bases[i].xcom = plist.xboxl * rand_gsl() - plist.xboxl / 2.0;
		bases[i].ycom = plist.yboxl * rand_gsl() - plist.yboxl / 2.0;
		if (bases[i].Dz == 0) {
			bases[i].zcom = -plist.zboxl / 2.0;
		} else {
			bases[i].zcom = plist.zboxl * rand_gsl() - plist.zboxl / 2.0;
		}
	}

	while (bflag == 1 && bit < Maxit) {

		bit++;
		bflag = 0;
		noverlap = 0;

		for (i = 0; i < plist.Ntotalmol; i++) {
			for (j = howmanylipids; j < plist.Ntotalmol; j++) {
				if(i!=j){
					for(ii=0;ii<bases[i].ninterface;ii++){
						for(jj=0;jj<bases[j].ninterface;jj++){

							dx = bases[i].xcom + coordcont[bases[i].protype][ii+1][0] - (bases[j].xcom + coordcont[bases[j].protype][jj+1][0]);
							dy = bases[i].ycom + coordcont[bases[i].protype][ii+1][1] - (bases[j].ycom + coordcont[bases[j].protype][jj+1][1]);
							dz = bases[i].zcom + coordcont[bases[i].protype][ii+1][2] - (bases[j].zcom + coordcont[bases[j].protype][jj+1][2]);
							dx -= plist.xboxl * round(dx / plist.xboxl);
							dy -= plist.yboxl * round(dy / plist.yboxl);
							d2 = dx * dx + dy * dy + dz * dz;

							//do ii and jj interact? if so find mu
							theyinteract = 0;
							for (rind=0; rind<plist.Nrxn; rind++) {
								if(rxtype[rind] == 0){

									i1 = Rlist[rind][0];//interface1
									i2 = Rlist[rind][1];//interface2

									if ((bases[i].protype==p_home[i1] && bases[j].protype==p_home[i2]) || (bases[i].protype==p_home[i2] && bases[j].protype==p_home[i1])){
										theyinteract = 1;
										mu = rind;
										rind = plist.Nrxn;
									}

								}
							}

							if(theyinteract==1){
								pbindr2 = bindrad[mu] * bindrad[mu];
							}
							else{
								pbindr2 = 0.0;
							}

							//make necessary measurements to check for overlap
							dx0 = bases[i].xcom + coordcont[bases[i].protype][ii+1][0] - (bases[j].xcom);
							dy0 = bases[i].ycom + coordcont[bases[i].protype][ii+1][1] - (bases[j].ycom);
							dz0 = bases[i].zcom + coordcont[bases[i].protype][ii+1][2] - (bases[j].zcom);
							dx0 -= plist.xboxl * round(dx / plist.xboxl);
							dy0 -= plist.yboxl * round(dy / plist.yboxl);
							d20 = dx0 * dx0 + dy0 * dy0 + dz0 * dz0;

							dx0arm = coordcont[bases[j].protype][ii+1][0];
							dy0arm = coordcont[bases[j].protype][ii+1][1];
							dz0arm = coordcont[bases[j].protype][ii+1][2];
							d20arm = dx0arm * dx0arm + dy0arm * dy0arm + dz0arm * dz0arm;

							dx1 = bases[i].xcom - (bases[j].xcom + coordcont[bases[j].protype][jj+1][0]);
							dy1 = bases[i].ycom - (bases[j].ycom + coordcont[bases[j].protype][jj+1][1]);
							dz1 = bases[i].zcom - (bases[j].zcom + coordcont[bases[j].protype][jj+1][2]);
							dx1 -= plist.xboxl * round(dx / plist.xboxl);
							dy1 -= plist.yboxl * round(dy / plist.yboxl);
							d21 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;

							dx1arm = coordcont[bases[i].protype][ii+1][0];
							dy1arm = coordcont[bases[i].protype][ii+1][1];
							dz1arm = coordcont[bases[i].protype][ii+1][2];
							d21arm = dx1arm * dx1arm + dy1arm * dy1arm + dz1arm * dz1arm;

							if(d21<d21arm || d20<d20arm){

								noverlap++;
								bflag = 1;
								ii=bases[i].ninterface;
								jj=bases[j].ninterface;
//								cout<<i<<'\t'<<j<<endl;

							}else if (d2 < pbindr2) {

								bases[j].xcom = plist.xboxl * rand_gsl() - plist.xboxl / 2.0;
								bases[j].ycom = plist.yboxl * rand_gsl() - plist.yboxl / 2.0;
								if (bases[j].Dz == 0) {
									bases[j].zcom = -plist.zboxl / 2.0;
								} else {
									bases[j].zcom = plist.zboxl * rand_gsl() - plist.zboxl / 2.0;
								}
								noverlap++;
								bflag = 1;
								ii=bases[i].ninterface;
								jj=bases[j].ninterface;
//								cout<<i<<'\t'<<j<<endl;

							}
						}
					}
				}

			}
		}
		cout << "Noverlap  " << noverlap << endl;
	}

	/*copy coords into interface pos and complex pos*/
	for (i = 0; i < plist.Ntotalmol; i++) {

		for (jj = 0; jj < bases[i].ninterface; jj++) {
			bases[i].x[jj] = bases[i].xcom+coordcont[bases[i].protype][jj+1][0];
			bases[i].y[jj] = bases[i].ycom+coordcont[bases[i].protype][jj+1][1];
			if (bases[i].Dz == 0) {
				bases[i].z[jj] = -plist.zboxl / 2.0 + coordcont[bases[i].protype][jj+1][2];
			} else {
				bases[i].z[jj] = bases[i].zcom+coordcont[bases[i].protype][jj+1][2];
			}
		}
		ind_com[i].xcom = bases[i].xcom;
		ind_com[i].ycom = bases[i].ycom;
		if (bases[i].Dz == 0) {
			ind_com[i].zcom = -plist.zboxl / 2.0;
		} else {
			ind_com[i].zcom = bases[i].zcom;
		}
		ind_com[i].plist[0] = i;
		ind_com[i].mysize = 1;
		bases[i].mycomplex = i;

		ind_com[i].Dx = bases[i].Dx;
		ind_com[i].Dy = bases[i].Dy;
		ind_com[i].Dz = bases[i].Dz;
		ind_com[i].Drx = bases[i].Drx;
		ind_com[i].Dry = bases[i].Dry;
		ind_com[i].Drz = bases[i].Drz;

	}

	t=0;
	for (p = 0; p < Nprotypes; p++) {
		for (i = 0; i < Ncopy[p]; i++) {

			bases[i].radR = wholep[p].radx;
			ind_com[i].radR = wholep[p].radx;
			cfile << let[p] << " " << bases[t].xcom << ' ' << bases[t].ycom << ' ' << bases[t].zcom << endl;
			comfile << let[p] << " " << bases[t].xcom << ' ' << bases[t].ycom << ' ' << bases[t].zcom << endl;

			for (j = 0; j < wholep[p].ninterface; j++) {
				px = bases[t].xcom+coordcont[bases[t].protype][j+1][0];
				py = bases[t].ycom+coordcont[bases[t].protype][j+1][1];
				pz = bases[t].zcom+coordcont[bases[t].protype][j+1][2];
				sprintf(buffer2,"%d",j);
				cfile << let[p] <<buffer2<< " " << px << ' ' << py << ' ' << pz << endl;
			}
			t++;
		}
	}

	if (ct != Nifaces) {
		cerr << "Number of interfaces assigned to proteins does not match network interface numbers!" << endl;
		exit(1);
	}
}
