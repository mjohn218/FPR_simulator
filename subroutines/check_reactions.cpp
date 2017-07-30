#include "reactions.h"
#include "utility_calls.h"

void check_reactions(Parms plist, Fullmol *bases, int *numpartners, int **Speclist, int *Nmyrxn, int **Rlist, int **myrxn, int Nspecies) {
	/*This routine makes sure that the network that defines which protein interfaces bind is
	 consistent with the definition between reactions between protein interfaces*/
	int i, j;
	int mu;
	int partner;
	for (i = 0; i < plist.Nifaces; i++) {
		cout << "Interface: " << i << " npartners: " << numpartners[i] << " Nrxns: " << Nmyrxn[i] << endl;
		for (j = 0; j < numpartners[i]; j++) {
			mu = myrxn[i][j];
			cout << " partner1: " << Speclist[i][j] << " rxn1: " << mu << " rxnpartners: " << Rlist[mu][0] << '\t' << Rlist[mu][1] << endl;
			if (Rlist[mu][0] == i)
				partner = Rlist[mu][1];
			else
				partner = Rlist[mu][0];
			if (Speclist[i][j] != partner)
				cout << "REACTION MISMATCH! interface: " << i << " reaction: " << mu << endl;
		}
		cout << " additional partners " << endl;
		for (j = numpartners[i]; j < Nmyrxn[i]; j++) {
			mu = myrxn[i][j];
			cout << " rxn1: " << mu << " rxnpartners: " << Rlist[mu][0] << '\t' << Rlist[mu][1] << endl;
		}
	}
	for (i = plist.Nifaces; i < Nspecies; i++) {
		cout << " Specie: " << i << " Nrxns: " << Nmyrxn[i] << endl;
		for (j = 0; j < Nmyrxn[i]; j++) {
			mu = myrxn[i][j];
			cout << "rxnnum: " << mu << " rxnpartners: " << Rlist[mu][0] << '\t' << Rlist[mu][1] << endl;
		}
	}

}
