#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "rand_gsl.h"
#include "vol_help.h"

using namespace std;

void numerical_test(double dx, double dy, double dz, double r1) {

	double R = sqrt(r1 * r1 - dy * dy);
	double xr, yr, zr, dr2, dr;
	double Vsphere = 0;
	double Vint = 0;
	double R2 = r1 * r1;

	double area = area_slice(dx, dy, dz, r1);
	double Vcir = M_PI * R * R;
	double Rad2 = R * R;
	Vsphere = 0;
	Vint = 0;
	int i;
	int Nit = 10000000;
	for (i = 0; i < Nit; i++) {
		xr = 2.0 * R * rand_gsl() - R;
		yr = 2.0 * R * rand_gsl() - R;
		dr2 = xr * xr + yr * yr;
		if (dr2 < Rad2) {
			Vsphere++;
			if (xr > dx && yr > dz)
				Vint++;
		}
	}
	cout << "MC Circle Area, Ratio Cutout:Sphere: " << Vint * 1.0 / (Vsphere * 1.0) << " Calculated: " << area / Vcir << endl;
	/*Try to numerically integrate*/
	int Nslice = 100;
	double yend = sqrt(r1 * r1 - dz * dz - dx * dx);
	double range = yend - dy;
	double step = range / Nslice;
	double total = 0;
	double myy;
	for (i = 0; i < Nslice; i++) {
		myy = i * step + dy;
		area = area_slice(dx, myy, dz, r1);
		total += area;
	}
	cout << "numerical integral: " << total * step << endl;

}
