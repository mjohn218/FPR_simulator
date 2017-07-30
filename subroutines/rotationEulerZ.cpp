#include "vector_rot_calls.h"

void rotationEulerZ(double tz, double *M) {
	/*Rx=[1 0 0; 0 costhet -sinthet ; 0 sinthet costhet]
	 Ry=[costhet 0 sinthet; 0 1 0 ; -sinthet 0 costhet]
	 Rz=[costhet -sinthet 0 ; sinthet costhet 0 ; 0 0 1]*/

	double sz = sin(tz);
	double cz = cos(tz);
	M[0] = cz;
	M[1] = -sz;
	M[2] = 0.0;
	M[3] = sz;
	M[4] = cz;
	M[5] = 0.0;
	M[6] = 0.0;
	M[7] = 0.0;
	M[8] = 1.0;

}
