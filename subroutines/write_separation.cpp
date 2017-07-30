#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void write_separation(Fullmol *bases, int p1, int p2) {
	double dx, dy, dz;
	dx = bases[p1].xcom - bases[p2].xcom;
	dy = bases[p1].ycom - bases[p2].ycom;
	dz = bases[p1].zcom - bases[p2].zcom;
	double d2 = dx * dx + dy * dy + dz * dz;
	cout << "Calulated Separation: " << sqrt(d2) << endl;

}
