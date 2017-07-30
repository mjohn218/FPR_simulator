#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "rand_gsl.h"
#include "vol_help.h"

using namespace std;

double Vcap_three_dim(int flagx, int flagy, int flagz, double dx, double dy, double dz, double R) {
	double h;
	double Vcap1, Vcap2, Vcap3;
	double Voverlap1, Voverlap2, Voverlap3;
	double Vremain;

	h = R - dx;
	Vcap1 = sphere_cap(R, h); //cap in x
	h = R - dy;
	Vcap2 = sphere_cap(R, h); //cap in y
	h = R - dz;
	Vcap3 = sphere_cap(R, h); //cap in z

	int check = 0;
	//check for overlap
	Voverlap1 = vcap_overlap(dy, dx, R);
	Voverlap2 = vcap_overlap(dz, dx, R);
	Voverlap3 = vcap_overlap(dy, dz, R);
	Vremain = vcap_three(dx, dy, dz, R); //this one is if yousubtract off same part twice  
	double Vcap = Vcap1 + Vcap2 + Vcap3 + Vremain - Voverlap1 - Voverlap2 - Voverlap3;
	return Vcap;
}
