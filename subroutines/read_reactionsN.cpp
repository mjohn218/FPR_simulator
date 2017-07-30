#include "reactions.h"
#include "utility_calls.h"

void read_reactionsN(ifstream &rxnfile, Fullmol *bases, int *Ncopy, int *p_home, Protein *wholep, int *numpartner, int **Speclist, Parms &plist, int **Rlist, int *Ncoup, int **mycoupled, int *Nmyrxn, int **myrxn, double *bindrad, double *kr, int *cntrxn, int *freelist, int *bndlist, int *zlist, int *rxtype, double *Kd) {

	/*Read in formatted list of reactions*/
	char rxntype;
	int i, pro1, pro2, npp, t=0;
	int rxnnum;
	int p1, p2, p3;
	int s;
	/*Only count reactions where you are a reactant!!!
	 as  a product, you just get created, you don't initiate any reactions 
	 */
	int j, flagexist=0;
	int nfree = 0;
	int nbnd = 0;
	int nmut = 0;
	for (j = 0; j < plist.Nspecies; j++) {
		Nmyrxn[j] = 0;
	}
	for (j = 0; j < plist.Nifaces; j++) {
		numpartner[j] = 0;
	}
	for (i = 0; i < plist.Nprotypes; i++) {
		wholep[i].npropart = 0;
	}
	rxnfile.ignore(600, '\n');
	cout << "N reactions: " << plist.Nrxn << endl;
	for (j = 0; j < plist.Nrxn; j++) {
		rxnfile >> rxnnum >> rxntype;
		if (rxntype == 'B') {
			rxtype[rxnnum] = 0;
			//binary reaction
			rxnfile >> Rlist[rxnnum][0] >> Rlist[rxnnum][1] >> Rlist[rxnnum][2];

			pro1 = p_home[Rlist[rxnnum][0]];
			pro2 = p_home[Rlist[rxnnum][1]];

			if (wholep[pro1].npropart<1){
				wholep[pro1].propart[wholep[pro1].npropart] = pro2;
				wholep[pro1].npropart++;
			}else{
				flagexist = 0;
				for(i=0;i<wholep[pro1].npropart;i++){
					if (wholep[pro1].propart[i]==pro2){
						flagexist = 1;
					}
				}
				if(flagexist==0){
					wholep[pro1].propart[wholep[pro1].npropart] = pro2;
					wholep[pro1].npropart++;
				}
			}
			if (wholep[pro2].npropart<1){
				wholep[pro2].propart[wholep[pro2].npropart] = pro1;
				wholep[pro2].npropart++;
			}else{
				flagexist = 0;
				for(i=0;i<wholep[pro2].npropart;i++){
					if (wholep[pro2].propart[i]==pro1){
						flagexist = 1;
					}
				}
				if(flagexist==0){
					wholep[pro2].propart[wholep[pro2].npropart] = pro1;
					wholep[pro2].npropart++;
				}
			}

			Speclist[Rlist[rxnnum][0]][numpartner[Rlist[rxnnum][0]]] = Rlist[rxnnum][1];
			Speclist[Rlist[rxnnum][1]][numpartner[Rlist[rxnnum][1]]] = Rlist[rxnnum][0];
			numpartner[Rlist[rxnnum][0]]++;
			numpartner[Rlist[rxnnum][1]]++;

			rxnfile >> Ncoup[rxnnum];
			for (i = 0; i < Ncoup[rxnnum]; i++) {
				rxnfile >> mycoupled[rxnnum][i];
			}

			p1 = Rlist[rxnnum][0];
			p2 = Rlist[rxnnum][1];
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

		} else if (rxntype == 'U') {
			//unimolecular reaction
			rxtype[rxnnum] = 1;
			rxnfile >> Rlist[rxnnum][0] >> Rlist[rxnnum][1] >> Rlist[rxnnum][2];
			rxnfile >> Ncoup[rxnnum];
			for (i = 0; i < Ncoup[rxnnum]; i++) {
				rxnfile >> mycoupled[rxnnum][i];
			}
//			/*For unbinding, also need to read in Kd in uM!!!!!!*/
//			rxnfile >> Kd[rxnnum];
			p1 = Rlist[rxnnum][0];
			s = Nmyrxn[p1];
			myrxn[p1][s] = rxnnum;
			Nmyrxn[p1]++;
			bndlist[nbnd] = p1;
			nbnd++;

		} else if (rxntype == 'Z') {
			//zeroth order
			rxtype[rxnnum] = 2;
			rxnfile >> Rlist[rxnnum][0] >> Rlist[rxnnum][1];
			rxnfile >> Ncoup[rxnnum];
			for (i = 0; i < Ncoup[rxnnum]; i++) {
				rxnfile >> mycoupled[rxnnum][i];
			}
			p1 = Rlist[rxnnum][0];
			zlist[nmut] = p1;
			nmut++;

		}
		/*Need reaction rates for the dissociation reactions*/
		rxnfile >> bindrad[rxnnum] >> kr[rxnnum];

	}
	cntrxn[0] = nfree;
	cntrxn[1] = nbnd;
	cntrxn[2] = nmut;

	for (int p = 0; p < plist.Nprotypes; p++) {
		npp = wholep[p].npropart;

		for (i = 0; i < Ncopy[p]; i++) {
			bases[t].npropart = npp;
			for (j = 0; j < npp; j++)
				bases[t].propart[j] = wholep[p].propart[j];

			for (j = 0; j < wholep[p].ninterface; j++) {
				bases[t].istatus[j] = wholep[p].valiface[j];
				bases[t].freelist[j] = wholep[p].valiface[j];
			}
			t++;
		}
	}

	for (i = 0; i < plist.Nprotypes; i++) {
		cout << "Protein: " << i << " Npropart: " << wholep[i].npropart << " Partners:" << endl;
		for (j = 0; j < wholep[i].npropart; j++) {
			cout << wholep[i].propart[j] << endl;
		}
	}

	for (i = 0; i < plist.Nifaces; i++) {
		cout << "Interface: " << i << " Nintpartner: " << numpartner[i] << " Partners:" << endl;
		for (j = 0; j < numpartner[i]; j++) {
			cout << Speclist[i][j] << endl;
		}

	}
	rxnfile.clear();
	rxnfile.seekg(0,rxnfile.beg);
}
