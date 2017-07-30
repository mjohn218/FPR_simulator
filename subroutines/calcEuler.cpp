#include "vector_rot_calls.h"

void calcEuler(double &alpha, double &beta, double &gamma, double z1, double z2, double z3, double x3, double y3) {
	//z-axis is [0]. and x is [1] axis
	//cout <<"In calc euler: "<<z1<<" -z2: "<<-z2<<" x3, y3: "<<x3<<' '<<y3<<endl;
	alpha = atan2(z1, -z2);
	if (alpha < 0)
		alpha += 2 * M_PI;
	beta = acos(z3);
	gamma = atan2(x3, y3);
	if (gamma < 0)
		gamma += 2 * M_PI;

}
