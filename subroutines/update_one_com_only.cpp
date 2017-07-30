#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "reactions.h"
#include "vol_help.h"
#include "rand_gsl.h"
#include "Faddeeva.hh"
#include "utility_calls.h"
#include "vector_rot_calls.h"

void update_one_com_only(int c1, Complex *ind_com, Fullmol *bases) {
	/*Update Complex COM of a single proteins*/
	double totalmassx = 0;
	double totalmassy = 0;
	double totalmassz = 0;
	int s1;
	int mp;
	int i, j;

	ind_com[c1].xcom = 0;
	ind_com[c1].ycom = 0;
	ind_com[c1].zcom = 0;

	s1 = ind_com[c1].mysize;
	/*Size of complex does not change*/
	for (i = 0; i < s1; i++) {
		mp = ind_com[c1].plist[i];
		totalmassx += bases[mp].massx;
		totalmassy += bases[mp].massy;
		totalmassz += bases[mp].massz;
		ind_com[c1].xcom += bases[mp].xcom * bases[mp].massx;
		ind_com[c1].ycom += bases[mp].ycom * bases[mp].massy;
		ind_com[c1].zcom += bases[mp].zcom * bases[mp].massz;

	}
	ind_com[c1].xcom /= totalmassx;
	ind_com[c1].ycom /= totalmassy;
	ind_com[c1].zcom /= totalmassz;

}
