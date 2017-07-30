#include "reactions.h"
#include "utility_calls.h"

void read_network_comment(Parms &plist, int *numpartner, int **Speclist, ifstream &infile) {

	int i, j;
	int ig;
	int prot, prot2;
	int iface2;
	infile.ignore(600, '\n');
	cout << "Reading in network topology from file..." << endl;
	int ind;
	cout << "Ninterfaces: " << plist.Nifaces << endl;
	for (i = 0; i < plist.Nifaces; i++) {
		infile >> ig >> numpartner[i];
		cout << "interface: " << i << " numpartners: " << numpartner[i] << endl;
		for (j = 0; j < numpartner[i]; j++) {
			infile >> iface2;
			Speclist[i][j] = iface2;
		}
		infile.ignore(600, '\n');
	}

}
