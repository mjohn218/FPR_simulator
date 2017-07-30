#include "reactions.h"
#include "utility_calls.h"

void set_status(ifstream &startfile, Protein *wholep, Fullmol *bases, Parms &plist, int *Ncopy) {
	int i, j, p;
	int ig, nint, nfree;
	int f;
	int t = 0;
	int flist[MAXPRTNER];
	int start[MAXIFACE];
	cout << "in set status " << endl;
	startfile.ignore(600, '\n');
	cout << "Nprotypes: " << plist.Nprotypes << endl;
	for (p = 0; p < plist.Nprotypes; p++) {
		nint = wholep[p].ninterface;
		startfile >> ig >> nfree;
		f = 0;
		cout << "p: " << p << " ninterface: " << nint << " nfree interfaces: " << nfree << endl;
		for (j = 0; j < nint; j++) {
			startfile >> start[j];
			//check to see if this interface is free for binding
			if (start[j] < plist.Nifaces) {
				flist[f] = start[j];
				f++;
			}
		}
		int npp = wholep[p].npropart;
		cout << "Pro: " << p << " Ncopy: " << Ncopy[p] << " Dx: " << wholep[p].Dx << " npropart: " << npp << endl;
		for (i = 0; i < Ncopy[p]; i++) {
			bases[t].ninterface = nint;
			bases[t].protype = p;
			bases[t].nfree = nfree;
			bases[t].npropart = npp;
			bases[t].Dx = wholep[p].Dx;
			bases[t].Dy = wholep[p].Dy;
			bases[t].Dz = wholep[p].Dz;
			bases[t].Drx = wholep[p].Drx;
			bases[t].Dry = wholep[p].Dry;
			bases[t].Drz = wholep[p].Drz;

			bases[t].massx = wholep[p].radx; //plist.mass;
			if (wholep[p].radx == 0)
				bases[t].massx = 1;
			bases[t].massy = wholep[p].rady; //plist.mass;
			if (wholep[p].rady == 0)
				bases[t].massy = 1;
			bases[t].massz = wholep[p].radz; //plist.mass;
			if (wholep[p].radz == 0)
				bases[t].massz = 1;

			// if(bases[t].Dz==0)
// 	bases[t].massz=0;

			for (j = 0; j < npp; j++)
				bases[t].propart[j] = wholep[p].propart[j];
			for (f = 0; f < nfree; f++)
				bases[t].freelist[f] = flist[f];
			for (j = 0; j < nint; j++) {
				bases[t].istatus[j] = start[j];
			}
			bases[t].nbnd = 0;
			bases[t].npartner = 0;
			bases[t].mycomplex = t;
			//bases[t].mycomplexsize=1;
			//cout <<"bases[: "<<t<<" ninterface: "<<bases[t].ninterface<<" nbnd: "<<bases[t].nbnd<<endl;
			t++;
			//cout <<"t: "<<t<<endl;
		}
	}

}
