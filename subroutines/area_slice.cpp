#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "rand_gsl.h"
#include "vol_help.h"

using namespace std;

double area_slice(double dx, double dy, double dz, double R) {
	double p = dx; //lower bound
	double h = dz;
	//bounds integrate over y are q1 to q2
	double g2 = (R * R - dy * dy);

	double g = sqrt(g2);
	//cout <<"g2: "<<g2<<" g: "<<g<<endl;
	double area = -0.5 * p * sqrt(g2 - p * p) - 0.5 * h * sqrt(g2 - h * h) + 0.5 * g2 * atan(sqrt(g2 - h * h) / h) - 0.5 * g2 * atan(p / sqrt(g2 - p * p)) + h * p;
	double f = sqrt(g2 - h * h);
	//cout <<"f: "<<f<<endl;
	double at = atan(f / sqrt(g2 - f * f));
	double at2 = atan(p / sqrt(g2 - p * p));

	double a2 = 0.5 * (f * sqrt(g2 - f * f) + g2 * at) - 0.5 * (p * sqrt(g2 - p * p) + g2 * at2) - h * (f - p);
	//  cout <<"atanf: "<<at<<" atanp: "<<at2<<" area2: "<<a2<<endl;
	return area;

}
