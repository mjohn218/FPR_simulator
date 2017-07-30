#include "reactions.h"
#include "utility_calls.h"

void read_protlist_twostate(int nwhole, Protein *wholep, int Nifaces, int *p_home, ifstream &protfile, int *ihome) {
	//read in each proteins number of interfaces and their identities
	int i, j;
	int nint;
	int ig;
	int ct = 0;
	int iface;
	protfile.ignore(500, '\n');
	for (i = 0; i < nwhole; i++) {
		protfile >> ig; //protein number, =i
		protfile >> wholep[i].Dx >> wholep[i].Dy >> wholep[i].Dz;
		protfile >> wholep[i].radx >> wholep[i].rady >> wholep[i].radz;
		protfile >> wholep[i].Drx >> wholep[i].Dry >> wholep[i].Drz;
		protfile >> nint; //number of protein i's interfaces
		wholep[i].ninterface = nint;
		wholep[i].nint_write = nint; //this can change below
		cout << "Protein: " << i << " Numinterfaces: " << wholep[i].ninterface << endl;
		cout << "Protein: " << i << " Dx: " << wholep[i].Dx << " Dz: " << wholep[i].Dz << endl;
		cout << "Protein: " << i << " radx: " << wholep[i].radx << " radz: " << wholep[i].radz << endl;
		for (j = 0; j < nint; j++) {
			protfile >> iface;
			cout << " index: " << iface << endl;
			wholep[i].valiface[j] = iface;
			p_home[iface] = i;
			ihome[iface] = j;
			/*Now membrane bound state!*/
			p_home[iface + 1] = i;
			ihome[iface + 1] = j;
		}
		protfile.ignore(500, '\n');
		ct += nint * 2;
	}
	//  protfile.ignore(600,'\n');
	protfile.ignore(600, '\n');
	int Nmore;
	int i1, p1, imatch;
	protfile >> Nmore;
	cout << "Nextra interface: " << Nmore << endl;
	ct += Nmore;
	for (i = 0; i < Nmore; i++) {
		//read in the identity of the interface, the protein it's on,
		//and the other interface that it replaces
		protfile >> i1 >> p1 >> imatch;
		p_home[i1] = p1;
		ihome[i1] = ihome[imatch];
		cout << "extra interface: " << i1 << " p_home: " << p1 << " ihome: " << ihome[i1] << endl;
	}
	protfile.ignore(600, '\n');
	protfile.ignore(600, '\n');
	for (i = 0; i < nwhole; i++) {
		protfile >> ig >> wholep[i].npropart;
		for (j = 0; j < wholep[i].npropart; j++) {
			protfile >> wholep[i].propart[j];
		}
		cout << "protein: " << i << " Npropart: " << wholep[i].npropart << endl;
		protfile.ignore(500, '\n');
	}
	protfile.ignore(600, '\n');
	protfile.ignore(600, '\n');
	for (i = 0; i < nwhole; i++) {
		protfile >> ig >> wholep[i].nint_write;
//     for(j=0;j<wholep[i].nint_write;j++){
//       protfile >>wholep[i].wrlist[j];
//     }
		cout << "Protein: " << i << " nint_write: " << wholep[i].nint_write << endl;
	}
	if (ct != Nifaces) {
		cerr << "Number of interfaces assigned to proteins does not match network interface numbers!" << endl;
		// exit(1);
	}

}
