#include "reactions.h"
#include "utility_calls.h"

void read_coords(ifstream &crdfile, int Nprotypes, Protein *wholep, Fullmol* bases, Complex *ind_com, int *Ncopy, string *names) {
	int i, j, k;
	int t = 0;
	int p;
	int nint;
	double x, y, z;
	string name2;
	string name1;
	//at start, each complex is an unbound protein 
	for (p = 0; p < Nprotypes; p++) {
		nint = wholep[p].ninterface;
		cout << "protein: " << p << " ninterface: " << nint << " Ncopy " << Ncopy[p] << endl;
		for (i = 0; i < Ncopy[p]; i++) {
			crdfile >> name1 >> bases[t].xcom >> bases[t].ycom >> bases[t].zcom;
			ind_com[t].xcom = bases[t].xcom;
			ind_com[t].ycom = bases[t].ycom;
			ind_com[t].zcom = bases[t].zcom;

			ind_com[t].mysize = 1;
			ind_com[t].plist[0] = t;
			ind_com[t].Dx = bases[t].Dx;
			ind_com[t].Dy = bases[t].Dy;
			ind_com[t].Dz = bases[t].Dz;
			ind_com[t].Drx = bases[t].Drx;
			ind_com[t].Dry = bases[t].Dry;
			ind_com[t].Drz = bases[t].Drz;

			bases[t].radR = wholep[p].radx;
			ind_com[t].radR = wholep[p].radx;
			cout << bases[t].xcom << ' ' << bases[t].ycom << ' ' << bases[t].zcom << endl;
			for (k = 0; k < nint; k++) {
				crdfile >> name2 >> x >> y >> z;
				bases[t].x[k] = x;
				bases[t].y[k] = y;
				bases[t].z[k] = z;
			}
			cout << bases[t].x[0] << ' ' << bases[t].y[0] << ' ' << bases[t].z[0] << endl;
			t++;
			//cout <<"t: "<<t<<" bases[t].nbnd: "<<bases[t].nbnd<<endl;
		}
		names[p] = name1;
		cout << name1 << endl;
	}

}
