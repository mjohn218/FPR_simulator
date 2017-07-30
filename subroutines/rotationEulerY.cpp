#include "vector_rot_calls.h"

void rotationEulerY(double ty, double *M) {
	/*Rx=[1 0 0; 0 costhet -sinthet ; 0 sinthet costhet]
	 Ry=[costhet 0 sinthet; 0 1 0 ; -sinthet 0 costhet]
	 Rz=[costhet -sinthet 0 ; sinthet costhet 0 ; 0 0 1]*/

	double sy = sin(ty);
	double cy = cos(ty);
	M[0] = cy;
	M[1] = 0.0;
	M[2] = sy;
	M[3] = 0.0;
	M[4] = 1.0;
	M[5] = 0.0;
	M[6] = -sy;
	M[7] = 0.0;
	M[8] = cy;

}
