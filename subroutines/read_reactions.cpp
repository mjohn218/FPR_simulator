#include "reactions.h"
#include "utility_calls.h"

void read_reactions(ifstream &rxnfile, Parms &plist, int **Rlist, int *Ncoup, int **mycoupled, int *Nmyrxn, int **myrxn, double *bindrad, double *kr, int *cntrxn, int *freelist, int *bndlist, int *zlist, int *rxtype, double *Kd) {

	/*Read in formatted list of reactions*/
	char rxntype;
	int i;
	int rxnnum;

	int p1, p2, p3;
	int s;
	/*Only count reactions where you are a reactant!!!
	 as  a product, you just get created, you don't initiate any reactions 
	 */
	int j;
	int nfree = 0;
	int nbnd = 0;
	int nmut = 0;
	for (j = 0; j < plist.Nspecies; j++) {
		Nmyrxn[j] = 0;
	}
	rxnfile.ignore(600, '\n');
	cout << "N reactions: " << plist.Nrxn << endl;
	for (j = 0; j < plist.Nrxn; j++) {
		rxnfile >> rxnnum >> rxntype;
		//    cout <<"rxnchar: "<<endl;
		//cout <<rxnnum<<' '<<rxntype<<endl;
		if (rxntype == 'B') {
			rxtype[rxnnum] = 0;
			//binary reaction
			rxnfile >> Rlist[rxnnum][0] >> Rlist[rxnnum][1] >> Rlist[rxnnum][2];
			rxnfile >> Ncoup[rxnnum];
			for (i = 0; i < Ncoup[rxnnum]; i++) {
				rxnfile >> mycoupled[rxnnum][i];
			}

			p1 = Rlist[rxnnum][0];
			p2 = Rlist[rxnnum][1];
			//      p3=Rlist[rxnnum][2];
			s = Nmyrxn[p1];
			myrxn[p1][s] = rxnnum;
			Nmyrxn[p1]++;
			freelist[nfree] = p1;
			nfree++;
			if (p1 != p2) {
				s = Nmyrxn[p2];
				myrxn[p2][s] = rxnnum;
				Nmyrxn[p2]++;
				freelist[nfree] = p2;
				nfree++;
			}
			// s=Nmyrxn[p3];
			//myrxn[p3][s]=rxnnum;
			//Nmyrxn[p3]++;

		} else if (rxntype == 'U') {
			//unimolecular reaction
			rxtype[rxnnum] = 1;
			rxnfile >> Rlist[rxnnum][0] >> Rlist[rxnnum][1] >> Rlist[rxnnum][2];
			rxnfile >> Ncoup[rxnnum];
			for (i = 0; i < Ncoup[rxnnum]; i++) {
				rxnfile >> mycoupled[rxnnum][i];
			}
			/*For unbinding, also need to read in Kd in uM!!!!!!*/
			rxnfile >> Kd[rxnnum];
			p1 = Rlist[rxnnum][0];
			//p2=Rlist[rxnnum][1];
			//p3=Rlist[rxnnum][2];

			s = Nmyrxn[p1];
			myrxn[p1][s] = rxnnum;
			Nmyrxn[p1]++;
			bndlist[nbnd] = p1;
			nbnd++;
			//  s=Nmyrxn[p2];
//       myrxn[p2][s]=rxnnum;
//       Nmyrxn[p2]++;
//       s=Nmyrxn[p3];
//       myrxn[p3][s]=rxnnum;
//       Nmyrxn[p3]++;

		} else if (rxntype == 'Z') {
			//zeroth order
			rxtype[rxnnum] = 2;
			rxnfile >> Rlist[rxnnum][0] >> Rlist[rxnnum][1];
			rxnfile >> Ncoup[rxnnum];
			for (i = 0; i < Ncoup[rxnnum]; i++) {
				rxnfile >> mycoupled[rxnnum][i];
			}
			p1 = Rlist[rxnnum][0];
			// p2=Rlist[rxnnum][1];
			//      s=Nmyrxn[p1];
			//myrxn[p1][s]=rxnnum;
			//Nmyrxn[p1]++;
			// s=Nmyrxn[p2];
//       myrxn[p2][s]=rxnnum;
//       Nmyrxn[p2]++;
			zlist[nmut] = p1;
			nmut++;

		}
		// cout <<rxnnum<<' '<<rxntype<<endl;
		//    cout <<"Nmy reactions: "<<endl;
		//cout <<"p1: "<<p1<<' '<<Nmyrxn[p1]<<" p2: "<<p2<<' '<<Nmyrxn[p2]<<endl;
		/*Need reaction rates for the dissociation reactions*/
		rxnfile >> bindrad[rxnnum] >> kr[rxnnum];

	}
	cntrxn[0] = nfree;
	cntrxn[1] = nbnd;
	cntrxn[2] = nmut;

}
