#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "rand_gsl.h"
#include "vol_help.h"

using namespace std;

void runmctest(double r1, double xcom, double ycom, double zcom, double Vcalc, double box_x) {
	/*Run an MC test*/
	int Nit = 1000000;
	double delx = box_x / 2.0 - xcom;
	double dely = box_x / 2.0 - ycom;
	double delz = box_x / 2.0 - zcom;
	double xr, yr, zr, dr2, dr;
	double Vsphere = 0;
	double Vint = 0;
	double R2 = r1 * r1;
	int i;
	for (i = 0; i < Nit; i++) {
		xr = 2.0 * r1 * rand_gsl() - r1;
		yr = 2.0 * r1 * rand_gsl() - r1;
		zr = 2.0 * r1 * rand_gsl() - r1;
		dr2 = xr * xr + yr * yr + zr * zr;
		if (dr2 < R2) {
			Vsphere++;
			if (xr < delx && yr < dely && zr < delz)
				Vint++;
		}
	}
	double Vcube = 8.0 * r1 * r1 * r1;
	double Vfull = 4.0 / 3.0 * M_PI * r1 * r1 * r1;
	cout << "MC Ratio Sphere:Square: " << Vsphere * 1.0 / (Nit * 1.0) << " Actual: " << Vfull / Vcube << endl;
	cout << "MC Ratio Truncated:Sphere: " << Vint * 1.0 / (Vsphere * 1.0) << " Calculated: " << Vcalc / Vfull << endl;

}
