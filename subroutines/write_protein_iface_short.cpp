#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void write_protein_iface_short(ofstream &outfile, Parms &plist, Fullmol *bases, int *Ncopy, int it, Protein *wholep, string *names) {
	int i, j;
	int t = 0;
//  outfile<<2*plist.Natomwrite-10<<endl;
	outfile << plist.Natomwrite << endl;
	outfile << "iter:" << it << endl;
	char aname;
	int n, nint;
	int ind;
	for (i = 0; i < plist.Nprotypes; i++) {

		nint = wholep[i].nint_write;
		for (j = 0; j < Ncopy[i]; j++) {
			//outfile<<bases[t].xcom<<' '<<bases[t].ycom<<' '<<bases[t].zcom<<' '<<bases[t].mycomplex<<' '<<bases[t].mycomplexsize<<endl;

			//outfile<<names[i]<<' '<<bases[t].xcom<<' '<<bases[t].ycom<<' '<<bases[t].zcom<<endl;
			for (n = 0; n < nint; n++) {
				//ind=wholep[i].wrlist[n];
				outfile << names[i] << ' ' << bases[t].x[n] << ' ' << bases[t].y[n] << ' ' << bases[t].z[n] << endl;
			}
			t++;
		}
	}

//  t = 0;
//  for(i=0;i<1;i++){
//
//    nint=wholep[i].nint_write;
//    for(j=0;j<Ncopy[i];j++){
//
//      for(n=0;n<nint;n++){
//    	  outfile<<'D'<<' '<<bases[t].x[n]<<' '<<bases[t].y[n]<<' '<<-200.0<<endl;
//      }
//      t++;
//    }
//  }

//  for(i=1;i<2;i++){
//
//    nint=wholep[i].nint_write;
//    for(j=0;j<Ncopy[i];j++){
//
//      for(n=0;n<nint;n++){
//    	  outfile<<'E'<<' '<<bases[t].x[n]<<' '<<bases[t].y[n]<<' '<<-200.0<<endl;
//      }
//      t++;
//    }
//  }

}
