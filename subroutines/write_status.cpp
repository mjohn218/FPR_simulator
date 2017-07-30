#include "reactions.h"

void write_status(ofstream &outfile, Protein *wholep, Fullmol *bases, Parms &plist, int *Ncopy, int iter) {
	int i, j, p;
	int ig;
	int f;
	int t = 0;
	int Ntotalmol = plist.Ntotalmol;
	outfile << "total complexes: " << plist.ntotalcomplex << " iter: " << iter << endl;
	for (i = 0; i < Ntotalmol; i++) {
		outfile << "protein: " << i << " mycomplex: " << bases[i].mycomplex << endl;
		outfile << "Nfree: " << bases[i].nfree << " Freelist: ";
		for (f = 0; f < bases[i].nfree; f++)
			outfile << bases[i].freelist[f] << '\t';
		outfile << endl;
		outfile << "Nbound: " << bases[i].nbnd << " Boundlist: ";
		for (f = 0; f < bases[i].nbnd; f++)
			outfile << bases[i].bndlist[f] << '\t';
		outfile << endl;
		outfile << "Niface: " << bases[i].ninterface << " Istatus: ";
		for (j = 0; j < bases[i].ninterface; j++) {
			outfile << bases[i].istatus[j] << '\t';
		}
		outfile << endl;
		outfile << "Partners: ";
		for (j = 0; j < bases[i].ninterface; j++)
			outfile << bases[i].partner[j] << '\t';
		outfile << endl;

	}
	/*From status of the proteins in the complex, should be able to reconstruct the
	 nature of each complex, its size, its prolist, and its diffusion and COM */

}
