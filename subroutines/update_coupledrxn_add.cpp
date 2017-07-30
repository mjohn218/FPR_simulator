#include "reactions.h"

void update_coupledrxn_add(int mu, int *Nmycoup, int **mycoupled, Fullmol *bases, int **Rlist, int *i_home, int *p_home, int p1, int p2) {
	/*If there are any reactions coupled to the ongoing reaction, loop over them and change
	 their reactants into products.
	 Here it is assumed that the coupled reaction is zeroth order A->B and creates a free-to-bind interface
	 */
	int i;
	int r1, i1, prod, iind;
	int ptype;
	int wprot;
	for (i = 0; i < Nmycoup[mu]; i++) {
		r1 = mycoupled[mu][i];
		//this reaction will also occur, probably a conformational change
		//order of reaction, zero

		i1 = Rlist[r1][0];
		prod = Rlist[r1][1]; //free interface
		iind = i_home[prod];
		ptype = p_home[prod];

		if (ptype == bases[p1].protype)
			wprot = p1;
		else
			wprot = p2;

		bases[wprot].istatus[iind] = prod;
		bases[wprot].freelist[bases[wprot].nfree] = prod;
		bases[wprot].nfree += 1;
		//cout <<"Added free interface on protein: "<<wprot<<" interface: "<<prod<<" index: "<<iind<<endl;

	}

}
