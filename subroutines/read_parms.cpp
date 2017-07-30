#include "reactions.h"
#include "utility_calls.h"

void read_parms(ifstream &parmfile, Parms &plist) {
	parmfile >> plist.Nprotypes;
	parmfile.ignore(400, '\n');
	parmfile >> plist.Nifaces;
	parmfile.ignore(400, '\n');
	parmfile >> plist.Nit;
	parmfile.ignore(400, '\n');
	parmfile >> plist.deltat;
	parmfile.ignore(400, '\n');
	parmfile >> plist.Nrxn;
	parmfile.ignore(400, '\n');
	parmfile >> plist.Nspecies;
	parmfile.ignore(400, '\n');
	parmfile >> plist.mass;
	parmfile.ignore(400, '\n');
	parmfile >> plist.xboxl;
	parmfile.ignore(400, '\n');
	parmfile >> plist.yboxl;
	parmfile.ignore(400, '\n');
	parmfile >> plist.zboxl;
	parmfile.ignore(400, '\n');
	parmfile >> plist.statwrite;
	parmfile.ignore(400, '\n');
	parmfile >> plist.configwrite;
	parmfile.ignore(400, '\n');
	parmfile >> plist.grwrite;
	parmfile.ignore(400, '\n');
	parmfile >> plist.restart;
	parmfile.ignore(400, '\n');
	parmfile >> plist.pclath;
	parmfile.ignore(400, '\n');

}
