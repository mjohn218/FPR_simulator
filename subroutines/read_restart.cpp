#include "reactions.h"
#include "utility_calls.h"

void read_restart(ifstream &restart, int Ntotalmol, Fullmol *bases, Parms &plist, Complex *ind_com) {
	int i, j, p;
	int ig;
	int f;
	char str[300];
	cout << "RESTART FROM PRIOR CONFIG " << endl;
	//startfile.ignore(600,'\n');

	restart >> str >> str >> plist.ntotalcomplex >> str >> ig;
	for (i = 0; i < plist.ntotalcomplex; i++) {
		ind_com[i].mysize = 0;
	}
	int c1, n1, nfree, nbnd, niface;
	for (i = 0; i < Ntotalmol; i++) {
		restart >> str >> ig >> str >> c1;
		bases[i].mycomplex = c1;
		n1 = ind_com[c1].mysize;
		ind_com[c1].mysize++;
		ind_com[c1].plist[n1] = i;
		restart >> str >> nfree >> str;
		bases[i].nfree = nfree;
		for (j = 0; j < nfree; j++)
			restart >> bases[i].freelist[j];
		restart >> str >> nbnd >> str;
		bases[i].nbnd = nbnd;
		bases[i].npartner = bases[i].nbnd;
		for (j = 0; j < nbnd; j++)
			restart >> bases[i].bndlist[j];
		restart >> str >> niface >> str;

		for (j = 0; j < niface; j++)
			restart >> bases[i].istatus[j];
		restart >> str;
		for (j = 0; j < niface; j++)
			restart >> bases[i].partner[j];

	}

	cout << "total complexes: " << plist.ntotalcomplex << endl;
	for (i = 0; i < plist.ntotalcomplex; i++) {
		cout << "complex: " << i << " size: " << ind_com[i].mysize << endl;
	}

}
